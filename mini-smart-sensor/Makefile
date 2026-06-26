# Makefile for mini-smart-sensor
# Smart Sensor Measurement and Fusion Library
#
# Reference: IEEE 1451.0-2007, JCGM 100:2008 (GUM)
#            Frank (2013) "Understanding Smart Sensors", 3rd ed.
#
# Targets:
#   make          - build library and tests
#   make test     - build and run tests
#   make examples - build all examples
#   make clean    - remove build artifacts

CC      := gcc
CFLAGS  := -Wall -Wextra -std=c11 -O2 -I./include
LDFLAGS := -lm

# Source files
SRC_DEFS   := src/sensor_definitions.c
SRC_COND   := src/sensor_conditioning.c
SRC_CAL    := src/sensor_calibration.c
SRC_FUSION := src/sensor_fusion.c
SRC_APP    := src/sensor_applications.c
SRC_ADV    := src/sensor_advanced.c

ALL_SRC := $(SRC_DEFS) $(SRC_COND) $(SRC_CAL) $(SRC_FUSION) $(SRC_APP) $(SRC_ADV)

# Objects
ALL_OBJ := $(ALL_SRC:.c=.o)

# Test
TEST_SRC := tests/test_smart_sensor.c
TEST_BIN := test_smart_sensor

# Examples
EX1_SRC := examples/example_thermistor_cal.c
EX2_SRC := examples/example_imu_fusion.c
EX3_SRC := examples/example_smart_node.c
EX1_BIN := example_thermistor_cal
EX2_BIN := example_imu_fusion
EX3_BIN := example_smart_node

.PHONY: all test examples clean

all: $(ALL_OBJ) $(TEST_BIN)

# Compile C source to object
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link test binary
$(TEST_BIN): $(ALL_OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(ALL_OBJ) -o $(TEST_BIN) $(LDFLAGS)

# Build examples
examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(EX1_BIN): $(ALL_OBJ) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) $(ALL_OBJ) -o $(EX1_BIN) $(LDFLAGS)

$(EX2_BIN): $(ALL_OBJ) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) $(ALL_OBJ) -o $(EX2_BIN) $(LDFLAGS)

$(EX3_BIN): $(ALL_OBJ) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) $(ALL_OBJ) -o $(EX3_BIN) $(LDFLAGS)

# Run tests
test: $(TEST_BIN)
	./$(TEST_BIN)

# Clean
clean:
	rm -f $(ALL_OBJ) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	rm -f src/*.o
