CC = gcc

CFLAGS = -fopenmp -O3 -Wall
NTT_CFLAGS = -fopenmp -O3
LIBS = -lm -lgmp -lflint

NTT_OBJS = ntt/profile.o ntt/misc.o ntt/modarith.o ntt/memory.o ntt/fft_main.o \
           ntt/fft_base.o ntt/fft_array.o ntt/intmult.o ntt/tunetab.o
TARGET = lfrs

all: $(NTT_OBJS) $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) -o $(TARGET) *.c gmp-ntt/*.c ntt/*.o $(LIBS)

$(NTT_OBJS): %.o: %.c
	$(CC) $(NTT_CFLAGS) -c $< -o $@

clean:
	rm -rf $(NTT_OBJS) $(TARGET)
