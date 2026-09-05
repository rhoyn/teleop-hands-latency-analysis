CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

STB_PREFIX ?= /opt/stb

CXXFLAGS += -I/usr/include/eigen3 -I$(STB_PREFIX)/include

LDLIBS := -pthread

BIN := build/teleop-hands-latency-analysis
REC := build/teleop-hands-recorder
SRC := main.cpp $(wildcard algos/*.cpp)
OBJ := $(patsubst %.cpp,build/obj/%.o,$(SRC))
DEP := $(OBJ:.o=.d)

.PHONY: all clean

all: $(BIN) $(REC)

$(BIN): $(OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

$(REC): recorder.cpp | build
	$(CXX) -std=c++17 -O2 -Wall -Wextra -o $@ recorder.cpp

build/obj/%.o: %.cpp | build
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -Ialgos -MMD -MP -c -o $@ $<

build:
	@mkdir -p $@

-include $(DEP)

clean:
	rm -rf build
