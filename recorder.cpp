#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

/**
 * Length of one packet from `esp32_firmware.cpp`, in bytes.
 *
 * Two sync bytes, a 32-bit microsecond stamp, twelve int16 channels and one
 * checksum. The firmware writes the struct field by field rather than
 * memcpy-ing a packed type, so this is the contract between the two
 * programs and neither can change it alone.
 */
static constexpr int PKT_LEN = 31;

/**
 * Accelerometer counts per g at the MPU6050's default +/-2 g range.
 */
static constexpr double ACC_LSB_PER_G = 16384.0;

/**
 * Gyro counts per degree per second at the default +/-250 dps range.
 */
static constexpr double GYR_LSB_PER_DPS = 131.0;

/**
 * Set once by SIGINT or SIGTERM so the read loop can leave cleanly.
 *
 * A recording ends when the operator stops it, so the interesting exit path
 * is the interrupt rather than end of input. Flushing the CSV on the way out
 * is what makes a Ctrl-C'd take usable instead of truncated mid-row.
 */
static volatile sig_atomic_t g_stop = 0;

/**
 * Ask the read loop to finish at the next iteration.
 *
 * @param[in] sig  signal number, ignored
 * @exceptsafe no-throw
 */
static void request_stop(int sig) {
  (void)sig;
  g_stop = 1;
}

/**
 * Map a baud rate to its termios constant.
 *
 * Only the rates an ESP32 will hold on a USB-serial bridge are listed. The
 * firmware runs at 921600, which carries a 31-byte packet every millisecond
 * with room to spare; anything slower drops samples rather than slowing the
 * sender, since the firmware never blocks on the write.
 *
 * @param[in] baud  rate in bits per second
 * @returns the termios speed, or `B0` if the rate is not supported
 * @exceptsafe no-throw
 */
static speed_t baud_to_speed(int baud) {
  switch (baud) {
    case 115200:
      return B115200;
    case 230400:
      return B230400;
    case 460800:
      return B460800;
    case 921600:
      return B921600;
    case 1000000:
      return B1000000;
    case 2000000:
      return B2000000;
    default:
      return B0;
  }
}

/**
 * Open a serial port in raw mode at the given rate.
 *
 * Raw mode matters: the default line discipline would eat 0x0D and expand
 * 0x0A, and the packet is binary, so a canonical terminal corrupts roughly
 * one sample in a hundred with no error anywhere. `VMIN` 0 with `VTIME` 1
 * makes `read` return after a decisecond of silence instead of blocking, so
 * an interrupt is noticed even when the board has stopped sending.
 *
 * @param[in] port  device path, normally `/dev/ttyUSB0` or `/dev/ttyACM0`
 * @param[in] baud  rate in bits per second
 * @returns the file descriptor, or -1 with the reason on stderr
 * @exceptsafe no-throw
 */
static int open_serial(
    const std::string& port,
    int baud
) {
  const speed_t sp = baud_to_speed(baud);
  if (sp == B0) {
    fprintf(stderr, "unsupported baud %d\n", baud);
    return -1;
  }
  const int fd = open(port.c_str(), O_RDONLY | O_NOCTTY);
  if (fd < 0) {
    fprintf(stderr, "cannot open %s: %s\n", port.c_str(), strerror(errno));
    return -1;
  }
  termios tio{};
  if (tcgetattr(fd, &tio) != 0) {
    fprintf(stderr, "tcgetattr %s: %s\n", port.c_str(), strerror(errno));
    close(fd);
    return -1;
  }
  cfmakeraw(&tio);
  cfsetispeed(&tio, sp);
  cfsetospeed(&tio, sp);
  tio.c_cflag |= (CLOCAL | CREAD);
  tio.c_cflag &= ~CRTSCTS;
  tio.c_cc[VMIN] = 0;
  tio.c_cc[VTIME] = 1;
  if (tcsetattr(fd, TCSANOW, &tio) != 0) {
    fprintf(stderr, "tcsetattr %s: %s\n", port.c_str(), strerror(errno));
    close(fd);
    return -1;
  }
  tcflush(fd, TCIFLUSH);
  return fd;
}

/**
 * Find the next 0xAA 0x55 sync pair at or after an offset.
 *
 * The stream is resynchronised on every packet rather than trusted to stay
 * aligned, because a USB hiccup drops bytes rather than whole packets and a
 * reader that assumed alignment would emit garbage for the rest of the run.
 *
 * @param[in] buf   bytes received so far
 * @param[in] from  index to start scanning at
 * @returns index of the 0xAA, or -1 if no pair is present
 * @exceptsafe no-throw
 */
static long find_sync(
    const std::vector<uint8_t>& buf,
    size_t from
) {
  for (size_t i = from; i + 1 < buf.size(); i++)
    if (buf[i] == 0xAA && buf[i + 1] == 0x55) return (long)i;
  return -1;
}

/**
 * Print the command line synopsis.
 *
 * @param[in] prog  program name as invoked, normally `argv[0]`
 * @exceptsafe no-throw
 */
static void usage(const char* prog) {
  fprintf(
      stderr,
      "usage: %s [--port /dev/ttyUSB0] [--baud 921600] [--out mpu_log.csv]\n",
      prog
  );
}

/**
 * Record one session from the ESP32 to a CSV the analysis can read.
 *
 * Every accepted packet becomes one row in the thirteen column layout
 * `main.cpp` parses, with the counts converted to g and degrees per second
 * here so the CSV is readable on its own. A packet whose XOR checksum fails
 * is counted and skipped rather than written, and the scan resumes two bytes
 * past the bad sync so a corrupt payload cannot hide the next good packet.
 * The running count on stderr is what tells the operator the rig is alive
 * before a take rather than after it.
 *
 * @param[in] argc  argument count
 * @param[in] argv  arguments, `argv[0]` skipped
 * @returns 0 on a clean stop, 1 if the port or the output could not be opened
 * @exceptsafe no-throw
 */
int main(
    int argc,
    char** argv
) {
  std::string port = "/dev/ttyUSB0";
  std::string out = "mpu_log.csv";
  int baud = 921600;

  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--port" && i + 1 < argc) {
      port = argv[++i];
    } else if (a == "--baud" && i + 1 < argc) {
      baud = atoi(argv[++i]);
    } else if (a == "--out" && i + 1 < argc) {
      out = argv[++i];
    } else {
      usage(argv[0]);
      return 1;
    }
  }

  signal(SIGINT, request_stop);
  signal(SIGTERM, request_stop);

  const int fd = open_serial(port, baud);
  if (fd < 0) return 1;

  FILE* csv = fopen(out.c_str(), "w");
  if (!csv) {
    fprintf(stderr, "cannot open %s: %s\n", out.c_str(), strerror(errno));
    close(fd);
    return 1;
  }
  fprintf(
      csv,
      "t_us,"
      "m1_ax_g,m1_ay_g,m1_az_g,m1_gx_dps,m1_gy_dps,m1_gz_dps,"
      "m2_ax_g,m2_ay_g,m2_az_g,m2_gx_dps,m2_gy_dps,m2_gz_dps\n"
  );

  fprintf(
      stderr,
      "reading %s @ %d  ->  %s   (Ctrl-C to stop)\n",
      port.c_str(),
      baud,
      out.c_str()
  );

  std::vector<uint8_t> buf;
  buf.reserve(4096);
  uint8_t rd[4096];
  uint64_t n = 0;
  uint64_t bad = 0;

  const auto t0 = std::chrono::steady_clock::now();
  auto last_report = t0;

  while (!g_stop) {
    const ssize_t r = read(fd, rd, sizeof(rd));
    if (r > 0) buf.insert(buf.end(), rd, rd + r);

    size_t consumed = 0;
    for (;;) {
      const long i = find_sync(buf, consumed);
      if (i < 0) {
        consumed = buf.empty() ? 0 : buf.size() - 1;
        break;
      }
      if (buf.size() - (size_t)i < PKT_LEN) {
        consumed = (size_t)i;
        break;
      }
      const uint8_t* p = &buf[i];
      uint8_t cs = 0;
      for (int k = 2; k < 30; k++) cs ^= p[k];
      if (cs != p[30]) {
        bad++;
        consumed = (size_t)i + 2;
        continue;
      }
      uint32_t t_us;
      int16_t v[12];
      memcpy(&t_us, p + 2, 4);
      memcpy(v, p + 6, 24);
      fprintf(
          csv,
          "%u,"
          "%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,"
          "%.5f,%.5f,%.5f,%.4f,%.4f,%.4f\n",
          t_us,
          v[0] / ACC_LSB_PER_G,
          v[1] / ACC_LSB_PER_G,
          v[2] / ACC_LSB_PER_G,
          v[3] / GYR_LSB_PER_DPS,
          v[4] / GYR_LSB_PER_DPS,
          v[5] / GYR_LSB_PER_DPS,
          v[6] / ACC_LSB_PER_G,
          v[7] / ACC_LSB_PER_G,
          v[8] / ACC_LSB_PER_G,
          v[9] / GYR_LSB_PER_DPS,
          v[10] / GYR_LSB_PER_DPS,
          v[11] / GYR_LSB_PER_DPS
      );
      n++;
      consumed = (size_t)i + PKT_LEN;
    }
    if (consumed) buf.erase(buf.begin(), buf.begin() + consumed);

    const auto now = std::chrono::steady_clock::now();
    if (now - last_report >= std::chrono::seconds(1)) {
      const double secs = std::chrono::duration<double>(now - t0).count();
      fprintf(
          stderr,
          "\r%llu samples  |  %6.1f /s  |  %llu bad",
          (unsigned long long)n,
          n / secs,
          (unsigned long long)bad
      );
      fflush(stderr);
      last_report = now;
    }
  }

  fflush(csv);
  fclose(csv);
  close(fd);

  const double secs =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
          .count();
  fprintf(
      stderr,
      "\nstopped. %llu samples in %.1fs (%.1f/s avg), %llu bad. saved %s\n",
      (unsigned long long)n,
      secs,
      secs > 0 ? n / secs : 0.0,
      (unsigned long long)bad,
      out.c_str()
  );
  return 0;
}
