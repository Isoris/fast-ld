CC ?= gcc
CFLAGS ?= -O3 -march=native -fopenmp -Wall -Wextra \
         -Wno-unused-parameter -Wno-unused-result \
         -Wno-format-truncation -Wno-missing-field-initializers
LDFLAGS ?= -fopenmp
LIBS ?= -lz -lm

SRC = src/fast_ld.c
BIN = src/fast_ld

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

# Static build for distribution. Requires static zlib + libgomp on the
# build host; on Ubuntu install zlib1g-dev libgomp1.
static:
	$(CC) $(CFLAGS) -static -o $(BIN).static $(SRC) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(BIN) $(BIN).static

test: $(BIN)
	cd tests && python3 test_fast_ld.py && python3 test_wrapper.py

bench: $(BIN)
	cd tests && python3 bench.py

.PHONY: all static clean test bench
