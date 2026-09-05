#include <Wire.h>

TwoWire I2C_A = TwoWire(0);
TwoWire I2C_B = TwoWire(1);
const uint8_t MPU = 0x68;
const uint32_t PERIOD_US = 1000;

static volatile int16_t g_b_acc[3], g_b_gyro[3];
static TaskHandle_t s_reader = nullptr;
static TaskHandle_t s_main = nullptr;

void initMPU(TwoWire& b) {
  b.beginTransmission(MPU);
  b.write(0x6B);
  b.write(0x00);
  b.endTransmission();
  b.beginTransmission(MPU);
  b.write(0x1A);
  b.write(0x03);
  b.endTransmission();
  b.beginTransmission(MPU);
  b.write(0x19);
  b.write(0x00);
  b.endTransmission();
}

bool readRaw(
    TwoWire& b,
    int16_t* a,
    int16_t* g
) {
  b.beginTransmission(MPU);
  b.write(0x3B);
  if (b.endTransmission(false) != 0) return false;
  if (b.requestFrom((uint8_t)MPU, (uint8_t)14) != 14) return false;
  uint8_t d[14];
  for (int i = 0; i < 14; i++) d[i] = b.read();
  a[0] = (int16_t)((d[0] << 8) | d[1]);
  a[1] = (int16_t)((d[2] << 8) | d[3]);
  a[2] = (int16_t)((d[4] << 8) | d[5]);
  g[0] = (int16_t)((d[8] << 8) | d[9]);
  g[1] = (int16_t)((d[10] << 8) | d[11]);
  g[2] = (int16_t)((d[12] << 8) | d[13]);
  return true;
}

void readerTask(void*) {
  I2C_B.begin(25, 26, 400000);
  initMPU(I2C_B);
  for (;;) {
    ulTaskNotifyTake(pdTRUE, portMAX_DELAY);
    int16_t a[3], g[3];
    if (!readRaw(I2C_B, a, g)) {
      a[0] = a[1] = a[2] = g[0] = g[1] = g[2] = 0;
    }
    for (int i = 0; i < 3; i++) {
      g_b_acc[i] = a[i];
      g_b_gyro[i] = g[i];
    }
    xTaskNotifyGive(s_main);
  }
}

void setup() {
  Serial.setTxBufferSize(2048);
  Serial.begin(921600);
  I2C_A.begin(32, 33, 400000);
  initMPU(I2C_A);
  s_main = xTaskGetCurrentTaskHandle();
  xTaskCreatePinnedToCore(readerTask, "reader", 4096, nullptr, 1, &s_reader, 0);
  delay(100);
}

void loop() {
  static uint32_t next = micros();
  uint32_t now = micros();
  if ((int32_t)(now - next) < 0) return;
  next += PERIOD_US;

  xTaskNotifyGive(s_reader);
  int16_t a1[3], g1[3];
  readRaw(I2C_A, a1, g1);
  ulTaskNotifyTake(pdTRUE, portMAX_DELAY);

  uint8_t pkt[31];
  pkt[0] = 0xAA;
  pkt[1] = 0x55;
  memcpy(&pkt[2], &now, 4);
  int16_t vals[12] = {
      a1[0],
      a1[1],
      a1[2],
      g1[0],
      g1[1],
      g1[2],
      g_b_acc[0],
      g_b_acc[1],
      g_b_acc[2],
      g_b_gyro[0],
      g_b_gyro[1],
      g_b_gyro[2]
  };
  memcpy(&pkt[6], vals, 24);
  uint8_t cs = 0;
  for (int i = 2; i < 30; i++) cs ^= pkt[i];
  pkt[30] = cs;

  Serial.write(pkt, sizeof(pkt));
}
