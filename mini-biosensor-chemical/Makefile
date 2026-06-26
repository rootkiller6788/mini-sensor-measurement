# Makefile for mini-biosensor-chemical
# Builds the biosensor C library, test suite, and examples.
#
# Targets:
#   make          — build library and test
#   make test     — build and run all tests
#   make examples — build all example programs
#   make all      — build everything
#   make clean    — remove build artifacts

CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -O2 -g
INCLUDES = -Iinclude
LDFLAGS  = -lm

SRC_DIR  = src
INC_DIR  = include
TST_DIR  = tests
EX_DIR   = examples
BLD_DIR  = build

# Source files
SRCS = $(SRC_DIR)/biosensor_types.c \
       $(SRC_DIR)/biosensor_electrochemical.c \
       $(SRC_DIR)/biosensor_optical.c \
       $(SRC_DIR)/biosensor_kinetics.c \
       $(SRC_DIR)/biosensor_signal.c \
       $(SRC_DIR)/biosensor_calibration.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)

# Test and examples
TEST_SRC = $(TST_DIR)/test_biosensor.c
TEST_BIN = $(BLD_DIR)/test_biosensor

EXAMPLE_SRCS = $(EX_DIR)/example_glucose_biosensor.c \
               $(EX_DIR)/example_elisa_calibration.c \
               $(EX_DIR)/example_sensor_array.c

EXAMPLE_BINS = $(EXAMPLE_SRCS:$(EX_DIR)/%.c=$(BLD_DIR)/%)

# Default target
.PHONY: all
all: test examples

# Create build directory
$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# Compile library objects
$(BLD_DIR)/%.o: $(SRC_DIR)/%.c | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Link test binary
$(TEST_BIN): $(OBJS) $(TEST_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJS) $(TEST_SRC) -o $@ $(LDFLAGS)

# Build and run tests
.PHONY: test
test: $(TEST_BIN)
	@echo "=== Running biosensor test suite ==="
	$(TEST_BIN)

# Build examples
$(BLD_DIR)/%: $(EX_DIR)/%.c $(OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(OBJS) $< -o $@ $(LDFLAGS)

.PHONY: examples
examples: $(EXAMPLE_BINS)

# Run individual examples
run-glucose: $(BLD_DIR)/example_glucose_biosensor
	$(BLD_DIR)/example_glucose_biosensor

run-elisa: $(BLD_DIR)/example_elisa_calibration
	$(BLD_DIR)/example_elisa_calibration

run-sensor-array: $(BLD_DIR)/example_sensor_array
	$(BLD_DIR)/example_sensor_array

# Run all examples
run-examples: examples
	@echo "=== Glucose Biosensor ==="
	$(BLD_DIR)/example_glucose_biosensor
	@echo ""
	@echo "=== ELISA Calibration ==="
	$(BLD_DIR)/example_elisa_calibration
	@echo ""
	@echo "=== Sensor Array ==="
	$(BLD_DIR)/example_sensor_array

# Count lines in include/ and src/
.PHONY: count
count:
	@echo "include/ headers:"
	@find include -name "*.h" -exec wc -l {} +
	@echo "src/ implementations:"
	@find src -name "*.c" -o -name "*.lean" | xargs wc -l
	@echo "Total include/ + src/:"
	@find . \( -path "*/include/*.h" -o -path "*/src/*.c" -o -path "*/src/*.lean" \) -exec cat {} + | wc -l

# Safety scan (per SKILL.md §10)
.PHONY: safety-scan
safety-scan:
	@echo "=== Filler detection scan ==="
	@grep -rn "_fn[0-9]" --include="*.c" --include="*.h" . || echo "  _fnN: CLEAN"
	@grep -rn "_aux[0-9]" --include="*.c" --include="*.h" . || echo "  _auxN: CLEAN"
	@grep -rn "_ext[0-9][^a-z]" --include="*.c" --include="*.h" . || echo "  _extN: CLEAN"
	@grep -rn "algorithm variant" --include="*.c" . || echo "  variant: CLEAN"
	@grep -rn "auxiliary function" --include="*.c" . || echo "  auxiliary: CLEAN"
	@grep -rn "extension point" --include="*.c" . || echo "  extension: CLEAN"
	@grep -rn "supplemental assert" --include="*.c" . || echo "  supplemental: CLEAN"
	@echo ""
	@echo "=== Lean filler detection ==="
	@grep -rn "SystemMetric" --include="*.lean" . || echo "  SystemMetric: CLEAN"
	@grep -rn "traceability_matrix.*:=\s*\[\]" --include="*.lean" . || echo "  traceability: CLEAN"
	@grep -rn ":= 0\.0" --include="*.lean" . || echo "  const return: CLEAN"
	@grep -rn "by trivial" --include="*.lean" . || echo "  by trivial: CLEAN"
	@grep -rn "\bsorry\b" --include="*.lean" . || echo "  sorry: CLEAN"
	@echo ""
	@echo "=== Stub detection ==="
	@find . -name "*.c" -size -200c -exec echo "  SUSPECT (small): {}" \;
	@echo "=== Safety scan complete ==="

# Clean
.PHONY: clean
clean:
	rm -rf $(BLD_DIR)

.PHONY: help
help:
	@echo "mini-biosensor-chemical Makefile"
	@echo "  make              - build library + test"
	@echo "  make test         - build and run test suite"
	@echo "  make examples     - build all examples"
	@echo "  make all          - build test + examples"
	@echo "  make run-glucose  - run glucose biosensor example"
	@echo "  make run-elisa    - run ELISA calibration example"
	@echo "  make run-sensor-array - run sensor array example"
	@echo "  make run-examples - run all examples"
	@echo "  make count        - count lines in include/ + src/"
	@echo "  make safety-scan  - filler/stub detection (SKILL §10)"
	@echo "  make clean        - remove build artifacts"
