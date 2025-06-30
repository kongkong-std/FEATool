# makefile real petsc without openmp
# module load generators/make/4.3/gcc-10.3.0 mpi/mpich/3.2.1/gcc-10.3.0 compilers/gcc/10.3.0
# module load tools/cJSON/1.7.18 tpl/gklib/gklib tpl/metis/metis tpl/parmetis/parmetis
CC = mpicc
SRC = ../src

CFLAGS = -c
CMD = -Wall -Wwrite-strings \
      -Wno-unknown-pragmas -Wno-lto-type-mismatch \
      -fstack-protector -fvisibility=hidden -g -O3 \
      -march=native -mtune=native -fPIC

INCL = -I/home/kongkong/software/petsc/3.21.4/real/gklib-metis/include \
       -I${MPICH_HOME}/include \
	   -I${CJSON_INCLUDE_PATH} \
	   -I${METIS_INCLUDE} \
	   -I${PARMETIS_INCLUDE}

OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: libMLA.a libMLA.so

# Static library
libMLA.a: $(OBJECT)
	ar rcs libMLA.a $(OBJECT)

# Shared library
libMLA.so: $(OBJECT)
	$(CC) -shared -fPIC -o libMLA.so $(OBJECT)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm -f libMLA.a libMLA.so $(OBJECT)
