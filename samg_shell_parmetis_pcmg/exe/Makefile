# makefile real petsc without openmp
# module load generators/make/4.3/gcc-10.3.0 mpi/mpich/3.2.1/gcc-10.3.0 compilers/gcc/10.3.0
# module load tools/cJSON/1.7.18 tpl/gklib/gklib tpl/metis/metis tpl/parmetis/parmetis
include config.mk

CC = ${MPICH_DIR}/bin/mpicc
SRC = ../src

CFLAGS = -c
CMD = -fPIC -Wall -Wwrite-strings \
      -Wno-unknown-pragmas -Wno-lto-type-mismatch \
      -fstack-protector -fvisibility=hidden -g -O3 \
      -march=native -mtune=native

ifneq ($(GCC_DIR),)
	GCC_ARCH    := $(shell basename $(wildcard $(GCC_DIR)/lib/gcc/*))
	GCC_VERSION := $(shell basename $(wildcard $(GCC_DIR)/lib/gcc/$(GCC_ARCH)/*))

	GCC_LIBPATH = \
        -L$(GCC_DIR)/lib \
        -L$(GCC_DIR)/lib64 \
        -L$(GCC_DIR)/lib/gcc/$(GCC_ARCH)/$(GCC_VERSION) \
        -Wl,-rpath,$(GCC_DIR)/lib \
        -Wl,-rpath,$(GCC_DIR)/lib64 \
        -Wl,-rpath,$(GCC_DIR)/lib/gcc/$(GCC_ARCH)/$(GCC_VERSION)
else
	GCC_LIBPATH =
endif

INCL = -I${PETSC_DIR}/include \
       -I${MPICH_DIR}/include \
	   -I${CJSON_DIR}/include \
	   -I${METIS_DIR}/include \
	   -I${PARMETIS_DIR}/include \
	   -I${GKLIB_DIR}/include

LIB_PATH = \
    -L$(MPICH_DIR)/lib \
    -Wl,-rpath,$(MPICH_DIR)/lib \
    -L$(PETSC_DIR)/lib \
    -Wl,-rpath,$(PETSC_DIR)/lib \
    -L$(METIS_DIR)/lib \
    -L$(PARMETIS_DIR)/lib \
    -L$(CJSON_DIR)/lib64 \
    -L$(GKLIB_DIR)/lib \
    -Wl,-rpath,$(METIS_DIR)/lib \
    -Wl,-rpath,$(PARMETIS_DIR)/lib \
    -Wl,-rpath,$(CJSON_DIR)/lib64 \
    $(GCC_LIBPATH)

LIB_FLAG = -lpetsc -lHYPRE -ldmumps -lsmumps -lzmumps -lcmumps -lmumps_common \
	   -lpord -lpthread -lscalapack -lopenblas \
	   -lpthread -lparmetis -lmetis -lGKlib -lmetis -lGKlib \
	   -lm -lX11 -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran \
	   -lm -lgcc_s -lquadmath -lstdc++ -lquadmath -ldl -lcjson

LIB = $(LIB_PATH) $(LIB_FLAG)
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_petsc_exe

app_petsc_exe: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_petsc_exe $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm -f app_petsc_exe $(OBJECT)
