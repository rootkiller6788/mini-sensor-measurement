# Makefile for mini-adc-dac-converter
# ADC/DAC Converter Theory & Implementation
#
# Usage:
#   make          - Build library and tests
#   make test     - Build and run tests
#   make examples - Build all examples
#   make clean    - Remove build artifacts

CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -O2 -g -Iinclude -Isrc
LDFLAGS  = -lm

# ── Source files ──────────────────────────────────────────────
SRC = src/adc_core.c        \
      src/dac_core.c        \
      src/adc_metrics.c     \
      src/quantization.c    \
      src/sampling.c        \
      src/sigma_delta.c     \
      src/sar_adc.c

# ── Object files ──────────────────────────────────────────────
OBJ = $(SRC:.c=.o)

# ── Library ───────────────────────────────────────────────────
LIB = build/libadc_dac.a

# ── Test files ────────────────────────────────────────────────
TEST_SRC = tests/test_adc_core.c       \
           tests/test_dac_core.c       \
           tests/test_quantization.c   \
           tests/test_sampling.c       \
           tests/test_sigma_delta.c    \
           tests/test_sar_adc.c        \
           tests/test_adc_metrics.c

TEST_BIN = $(TEST_SRC:.c=.test)

# ── Example files ─────────────────────────────────────────────
EXAMPLE_SRC = examples/example_adc_sweep.c    \
              examples/example_sar_sim.c      \
              examples/example_sigma_delta.c  \
              examples/example_sampling.c

EXAMPLE_BIN = $(EXAMPLE_SRC:.c=.exe)

# ── Targets ───────────────────────────────────────────────────
.PHONY: all test examples clean

all: $(LIB) test

$(LIB): $(OBJ)
	@mkdir -p build
	ar rcs $@ $(OBJ)

# Pattern rule for C compilation
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# ── Tests ─────────────────────────────────────────────────────
test: $(TEST_BIN)
	@echo "=== Running ADC/DAC Converter Tests ==="
	@passed=0; failed=0; \
	for t in $(TEST_BIN); do \
		if [ -f $$t ]; then \
			echo "--- $$t ---"; \
			if ./$$t; then \
				passed=$$((passed+1)); \
			else \
				failed=$$((failed+1)); \
				echo "FAIL: $$t"; \
			fi; \
		fi; \
	done; \
	echo "=== Results: $$passed passed, $$failed failed ==="; \
	test $$failed -eq 0

$(TEST_BIN): %.test: %.c $(LIB)
	$(CC) $(CFLAGS) $< -Lbuild -ladc_dac $(LDFLAGS) -o $@

# ── Examples ──────────────────────────────────────────────────
examples: $(EXAMPLE_BIN)

$(EXAMPLE_BIN): %.exe: %.c $(LIB)
	$(CC) $(CFLAGS) $< -Lbuild -ladc_dac $(LDFLAGS) -o $@

# ── Clean ─────────────────────────────────────────────────────
clean:
	rm -rf build
	rm -f src/*.o
	rm -f tests/*.test
	rm -f examples/*.exe
