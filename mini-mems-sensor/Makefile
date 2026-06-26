# Makefile for mini-mems-sensor
# MEMS Sensors — L1-L8 knowledge coverage
# IEEE 1293/2700/1139/952 MEMS standards

CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -O2 -g
LDFLAGS  = -lm
INCLUDES = -Iinclude

SRCDIR   = src
TESTDIR  = tests
EXAMPLEDIR = examples
BENCHDIR = benches

SRCS     = src/mems_basics.c src/mems_dynamics.c src/mems_noise.c \
           src/mems_calibration.c src/mems_fusion.c src/mems_applications.c
OBJS     = $(SRCS:.c=.o)

TESTS    = tests/test_basics.c tests/test_dynamics.c \
           tests/test_noise.c tests/test_calibration.c
TESTBINS = $(TESTS:.c=.test)

EXAMPLES = examples/example_imu_calib.c examples/example_accel_tilt.c \
           examples/example_sensor_fusion.c
EXAMPLEBINS = $(EXAMPLES:.c=.exe)

.PHONY: all clean test examples bench check

all: $(OBJS)
	@echo "mini-mems-sensor: all objects built"

src/%.o: src/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test: $(TESTBINS)
	@echo "=== Running MEMS sensor tests ==="
	@for t in $(TESTBINS); do echo "--- $$(basename $$t) ---" && ./$$t && echo "  PASS" || echo "  FAIL"; done

tests/%.test: tests/%.c $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) $(LDFLAGS) -o $@

examples: $(EXAMPLEBINS)
	@for e in $(EXAMPLEBINS); do echo "--- $$(basename $$e) ---" && ./$$e && echo ""; done

examples/%.exe: examples/%.c $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) $(LDFLAGS) -o $@

check: $(OBJS) $(TESTBINS)

clean:
	rm -f src/*.o tests/*.test examples/*.exe
