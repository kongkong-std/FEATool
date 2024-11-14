# makefile for FEATool
# module load generators/make/4.3/gcc-10.3.0
CC = gcc
CFLAGS = -Wall -g -O3
SRC = ../src
OBJECT = $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))
# LIB = -lm
LIB = 

.PHONY: all clean
all: app_exe

app_exe: $(OBJECT)
	$(CC) $(CFLAGS) -o app_exe $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) -c "$<" -o "$@"

clean:
	rm -f app_exe $(OBJECT)