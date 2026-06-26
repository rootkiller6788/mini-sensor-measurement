# Makefile for mini-sensor-principle
# Sensor measurement principles — transfer functions, noise, bridge, calibration

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude -MMD -MP
LDFLAGS = -lm

SRC_DIR = src
TEST_DIR = tests
EXAMPLE_DIR = examples
BUILD_DIR = build

# Source files
SRCS = $(SRC_DIR)/sensor_defs.c \
       $(SRC_DIR)/sensor_transfer.c \
       $(SRC_DIR)/sensor_noise.c \
       $(SRC_DIR)/sensor_bridge.c \
       $(SRC_DIR)/sensor_calibration.c \
       $(SRC_DIR)/sensor_signal_cond.c

# Test files
TESTS = $(TEST_DIR)/test_sensor_core.c \
        $(TEST_DIR)/test_sensor_bridge_cal.c

# Example files
EXAMPLES = $(EXAMPLE_DIR)/ex_thermocouple_cjc.c \
           $(EXAMPLE_DIR)/ex_strain_gauge_bridge.c \
           $(EXAMPLE_DIR)/ex_accelerometer_sim.c

# Objects
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
TEST_BINS = $(TESTS:$(TEST_DIR)/%.c=$(BUILD_DIR)/%)
EXAMPLE_BINS = $(EXAMPLES:$(EXAMPLE_DIR)/%.c=$(BUILD_DIR)/%)

.PHONY: all clean test examples

all: test examples

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Build test binaries
$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

# Build example binaries
$(BUILD_DIR)/ex_%: $(EXAMPLE_DIR)/ex_%.c $(OBJS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $< $(OBJS) -o $@ $(LDFLAGS)

# Run all tests
test: $(TEST_BINS)
	@echo "=== Running Tests ==="
	@for t in $(TEST_BINS); do \
		echo; \
		$$t || exit 1; \
	done
	@echo "\n=== All tests passed ==="

# Build examples
examples: $(EXAMPLE_BINS)

# Run all examples
run-examples: examples
	@for e in $(EXAMPLE_BINS); do \
		echo; \
		echo "--- Running $$e ---"; \
		$$e; \
	done

clean:
	rm -rf $(BUILD_DIR)

-include $(BUILD_DIR)/*.d
