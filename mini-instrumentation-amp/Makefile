# Makefile - mini-instrumentation-amp
# Instrumentation Amplifier: topologies, analysis, calibration, sensor interfaces

CC       = gcc
CFLAGS   = -std=gnu11 -Wall -Wextra -O2 -g -Iinclude
LDFLAGS  = -lm

SRC_DIR  = src
TST_DIR  = tests
EX_DIR   = examples
DM_DIR   = demos
BN_DIR   = benches
BLD_DIR  = build

SRCS = $(SRC_DIR)/inst_amp_core.c \
       $(SRC_DIR)/inst_amp_topology.c \
       $(SRC_DIR)/inst_amp_analysis.c \
       $(SRC_DIR)/inst_amp_calibration.c \
       $(SRC_DIR)/inst_amp_sensorif.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)
LIB  = $(BLD_DIR)/libinstamp.a

TEST_SRC = $(TST_DIR)/test_inst_amp.c
TEST_BIN = $(BLD_DIR)/test_inst_amp

EX1_SRC = $(EX_DIR)/ex1_bridge_sensor.c
EX2_SRC = $(EX_DIR)/ex2_thermocouple.c
EX3_SRC = $(EX_DIR)/ex3_ecg_frontend.c
EX1_BIN = $(BLD_DIR)/ex1_bridge_sensor
EX2_BIN = $(BLD_DIR)/ex2_thermocouple
EX3_BIN = $(BLD_DIR)/ex3_ecg_frontend

DEMO_SRC = $(DM_DIR)/demo_ia_response.c
DEMO_BIN = $(BLD_DIR)/demo_ia_response
BENCH_SRC = $(BN_DIR)/bench_inst_amp.c
BENCH_BIN = $(BLD_DIR)/bench_inst_amp

.PHONY: all test examples demo bench clean

all: $(LIB) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(LIB): $(BLD_DIR) $(OBJS)
	ar rcs $@ $(OBJS)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(LIB) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) -L$(BLD_DIR) -linstamp $(LDFLAGS) -o $@

$(EX1_BIN): $(LIB) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) -L$(BLD_DIR) -linstamp $(LDFLAGS) -o $@

$(EX2_BIN): $(LIB) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) -L$(BLD_DIR) -linstamp $(LDFLAGS) -o $@

$(EX3_BIN): $(LIB) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) -L$(BLD_DIR) -linstamp $(LDFLAGS) -o $@

$(DEMO_BIN): $(LIB) $(DEMO_SRC)
	$(CC) $(CFLAGS) $(DEMO_SRC) -L$(BLD_DIR) -linstamp $(LDFLAGS) -o $@

$(BENCH_BIN): $(LIB) $(BENCH_SRC)
	$(CC) $(CFLAGS) $(BENCH_SRC) -L$(BLD_DIR) -linstamp $(LDFLAGS) -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

demo: $(DEMO_BIN)

bench: $(BENCH_BIN)

clean:
	rm -rf $(BLD_DIR)
