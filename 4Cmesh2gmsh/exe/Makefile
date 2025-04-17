# module load generators/make/4.3/gcc-10.3.0
CC = gcc
SRC = ../src

CFLAGS = -c
CMD = -fPIC -Wall -Wwrite-strings \
      -Wno-unknown-pragmas -Wno-lto-type-mismatch \
      -fstack-protector -fvisibility=hidden -g -O3 \
      -march=native -mtune=native

INCL = 

LIB_PATH = 

LIB_FLAG = 

LIB = $(LIB_PATH) $(LIB_FLAG)
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_msh_exe

app_msh_exe: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_msh_exe $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm -f app_msh_exe $(OBJECT)
