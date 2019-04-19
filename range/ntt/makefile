GMP_DIR = $(HOME)/local/gmp-install
MPIR_DIR = $(HOME)/local/mpir-install

FLAGS = -fopenmp -m64 -fgnu89-inline -std=c99 -O2 -lgmp -lm

CC = gcc

HEADERS = ntt.h ntt-internal.h test.h
SOURCES = profile.c misc.c modarith.c memory.c fft_main.c fft_base.c fft_array.c intmult.c tunetab.c

all: tune-gmp tune-mpir test-gmp test-mpir

tune-gmp: $(SOURCES) $(HEADERS) tune.c
	$(CC) $(FLAGS) -o tune-gmp tune.c $(SOURCES) -I$(GMP_DIR)/include -L$(GMP_DIR)/lib -Wl,-rpath -Wl,$(GMP_DIR)/lib

tune-mpir: $(SOURCES) $(HEADERS) tune.c
	$(CC) $(FLAGS) -o tune-mpir tune.c $(SOURCES) -I$(MPIR_DIR)/include -L$(MPIR_DIR)/lib -Wl,-rpath -Wl,$(MPIR_DIR)/lib

test-gmp: $(SOURCES) $(HEADERS) test.c
	$(CC) $(FLAGS) -o test-gmp test.c $(SOURCES) -I$(GMP_DIR)/include -L$(GMP_DIR)/lib -Wl,-rpath -Wl,$(GMP_DIR)/lib

test-mpir: $(SOURCES) $(HEADERS) test.c
	$(CC) $(FLAGS) -o test-mpir test.c $(SOURCES) -I$(MPIR_DIR)/include -L$(MPIR_DIR)/lib -Wl,-rpath -Wl,$(MPIR_DIR)/lib

clean:
	rm -rf test-gmp test-mpir
