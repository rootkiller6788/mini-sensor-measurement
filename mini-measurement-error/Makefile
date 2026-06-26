# Makefile for mini-measurement-error
# Measurement Error Analysis - L1-L8 knowledge coverage
# ISO 5725 / GUM JCGM 100:2008 / IEC 60051

CC       = gcc
CFLAGS   = -Wall -Wextra -std=c11 -O2 -g
LDFLAGS  = -lm
INCLUDES = -Iinclude

SRCDIR   = src
TESTDIR  = tests
EXAMPLEDIR = examples
BENCHDIR = benches

SRCS     = src/measurement_error.c src/measurement_uncertainty.c src/measurement_calibration.c src/measurement_statistics.c src/measurement_filtering.c src/measurement_compensation.c src/measurement_applications.c
OBJS     = $(SRCS:.c=.o)

TESTS    = tests/test_error.c tests/test_uncertainty.c tests/test_calibration.c
TESTBINS = $(TESTS:.c=.test)

EXAMPLES = examples/example_thermocouple.c examples/example_strain_gauge.c examples/example_adc.c examples/example_sensor_fusion.c
EXAMPLEBINS = $(EXAMPLES:.c=.exe)

.PHONY: all clean test examples bench check

all: $(OBJS)

src/%.o: src/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

test: $(TESTBINS)
	@echo "=== Running measurement error tests ==="
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
