CC    = icx
COPTS = -Wall -O3 -qopenmp

ifdef PAD_LD
COPTS += -DPAD_LD
endif

ifndef MKLROOT
$(error "MKLROOT is undefined")
endif

MKL_LIB_DIR = $(MKLROOT)/lib/intel64
MKL_LIBS    = -Wl,--start-group $(MKL_LIB_DIR)/libmkl_intel_lp64.a $(MKL_LIB_DIR)/libmkl_intel_thread.a $(MKL_LIB_DIR)/libmkl_core.a -Wl,--end-group
MKL_INC     = $(MKLROOT)/include
AUX_LIBS    = -lpthread -lm

ALL_LIBS    = $(MKL_LIBS) $(AUX_LIBS)

SRC = gemm_bench.c
OBJ = sgemm_bench.o dgemm_bench.o cgemm_bench.o zgemm_bench.o

sgemm.bin: sgemm_bench.o
	$(CC) $(COPTS) -o $@ $< $(ALL_LIBS)
sgemm_bench.o: gemm_bench.c
	$(CC) -c -DSINGLE $(COPTS) -I$(MKL_INC) -o $@ $<

zgemm.bin: zgemm_bench.o
	$(CC) $(COPTS) -o $@ $< $(ALL_LIBS)
zgemm_bench.o: gemm_bench.c
	$(CC) -c -DCOMPLEX16 $(COPTS) -I$(MKL_INC) -o $@ $<

dgemm.bin: dgemm_bench.o
	$(CC) $(COPTS) -o $@ $< $(ALL_LIBS)
dgemm_bench.o: gemm_bench.c
	$(CC) -c -DDOUBLE $(COPTS) -I$(MKL_INC) -o $@ $<

cgemm.bin: cgemm_bench.o
	$(CC) $(COPTS) -o $@ $< $(ALL_LIBS)
cgemm_bench.o: gemm_bench.c
	$(CC) -c -DCOMPLEX8 $(COPTS) -I$(MKL_INC) -o $@ $<

all: sgemm.bin dgemm.bin cgemm.bin zgemm.bin

clean:
	rm -rf $(OBJ) *.bin

.PHONY: all clean run
