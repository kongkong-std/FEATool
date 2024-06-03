# makefile for FEATool
CC = gcc
CFLAGS = -Wall -g
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
	rm app_exe $(OBJECT)