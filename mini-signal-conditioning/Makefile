# Makefile - mini-signal-conditioning
# Signal Conditioning: filters, linearization, excitation, isolation, CJC, bridges

CC       = gcc
CFLAGS   = -std=gnu11 -Wall -Wextra -Wno-unused-parameter -O2 -g -Iinclude
LDFLAGS  = -lm

SRC_DIR  = src
TST_DIR  = tests
EX_DIR   = examples
DM_DIR   = demos
BN_DIR   = benches
BLD_DIR  = build

SRCS = $(SRC_DIR)/sigcond_filter.c \
       $(SRC_DIR)/sigcond_linearize.c \
       $(SRC_DIR)/sigcond_excitation.c \
       $(SRC_DIR)/sigcond_isolation.c \
       $(SRC_DIR)/sigcond_cjc.c \
       $(SRC_DIR)/sigcond_bridge.c \
       $(SRC_DIR)/sigcond_selfcal.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)
LIB  = $(BLD_DIR)/libsigcond.a

TEST_SRC = $(TST_DIR)/test_sigcond.c
TEST_BIN = $(BLD_DIR)/test_sigcond

EX1_SRC = $(EX_DIR)/ex1_thermocouple_chain.c
EX2_SRC = $(EX_DIR)/ex2_strain_gage_system.c
EX3_SRC = $(EX_DIR)/ex3_rtd_measurement.c
EX1_BIN = $(BLD_DIR)/ex1_thermocouple_chain
EX2_BIN = $(BLD_DIR)/ex2_strain_gage_system
EX3_BIN = $(BLD_DIR)/ex3_rtd_measurement

DEMO_SRC = $(DM_DIR)/demo_filter_response.c
DEMO_BIN = $(BLD_DIR)/demo_filter_response
BENCH_SRC = $(BN_DIR)/bench_sigcond.c
BENCH_BIN = $(BLD_DIR)/bench_sigcond

.PHONY: all test examples demo bench clean

all: $(LIB) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(LIB): $(BLD_DIR) $(OBJS)
	ar rcs $@ $(OBJS)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(LIB) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) -L$(BLD_DIR) -lsigcond $(LDFLAGS) -o $@

$(EX1_BIN): $(LIB) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) -L$(BLD_DIR) -lsigcond $(LDFLAGS) -o $@

$(EX2_BIN): $(LIB) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) -L$(BLD_DIR) -lsigcond $(LDFLAGS) -o $@

$(EX3_BIN): $(LIB) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) -L$(BLD_DIR) -lsigcond $(LDFLAGS) -o $@

$(DEMO_BIN): $(LIB) $(DEMO_SRC)
	$(CC) $(CFLAGS) $(DEMO_SRC) -L$(BLD_DIR) -lsigcond $(LDFLAGS) -o $@

$(BENCH_BIN): $(LIB) $(BENCH_SRC)
	$(CC) $(CFLAGS) $(BENCH_SRC) -L$(BLD_DIR) -lsigcond $(LDFLAGS) -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

demo: $(DEMO_BIN)

bench: $(BENCH_BIN)

clean:
	rm -rf $(BLD_DIR)
