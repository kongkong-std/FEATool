# makefile real petsc without openmp
# module load make/4.3/gcc-10.3.0 mpi/mpich/3.2.1/gcc-10.3.0 compilers/gcc/10.3.0
# module load cjson/1.7.18
CC = mpicc
SRC = ../src

CFLAGS = -c
CMD = -fPIC -Wall -Wwrite-strings \
      -Wno-unknown-pragmas -Wno-lto-type-mismatch \
      -fstack-protector -fvisibility=hidden -g -O3 \
      -march=native -mtune=native

INCL = -I/home/kongkong/software/petsc/3.21.4/real/0openmp/include \
       -I/opt/apps/mpich/3.2.1/gcc-10.3.0/include \
	   -I${CJSON_INCLUDE_PATH}

LIB_PATH = -Wl,-rpath,/home/kongkong/software/petsc/3.21.4/real/0openmp/lib \
	   -L/home/kongkong/software/petsc/3.21.4/real/0openmp/lib \
	   -Wl,-rpath,/opt/apps/mpich/3.2.1/gcc-10.3.0/lib \
	   -L/opt/apps/mpich/3.2.1/gcc-10.3.0/lib \
	   -Wl,-rpath,/opt/apps/compilers/gcc-10.3.0/lib/gcc/x86_64-pc-linux-gnu/10.3.0 \
	   -L/opt/apps/compilers/gcc-10.3.0/lib/gcc/x86_64-pc-linux-gnu/10.3.0 \
	   -Wl,-rpath,/opt/apps/compilers/gcc-10.3.0/lib64 \
	   -L/opt/apps/compilers/gcc-10.3.0/lib64 \
	   -Wl,-rpath,/opt/apps/compilers/gcc-10.3.0/lib \
	   -L/opt/apps/compilers/gcc-10.3.0/lib \
	   -L${CJSON_LIB_PATH} \
	   -Wl,-rpath,${CJSON_LIB_PATH}

LIB_FLAG = -lpetsc -lHYPRE -ldmumps -lsmumps -lzmumps -lcmumps -lmumps_common \
	   -lpord -lpthread -lscalapack -lopenblas \
	   -lptesmumps -lptscotchparmetisv3 -lptscotch -lptscotcherr \
	   -lesmumps -lscotch -lscotcherr -lpthread \
	   -lparmetis -lmetis -lm -lX11 -ldl \
	   -lmpifort -lmpi -lgfortran -lm -lgfortran -lm \
	   -lgcc_s -lquadmath -lstdc++ -lrt -lquadmath -ldl -lcjson

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
