CC = gcc
CFLAGS = -Wall -Wextra -O2

SRC = main.c vec.c objfile.c floatn.c
OBJ = $(SRC:.c=.o)
DEPS = vec.h objfile.h floatn.h

TARGET = svd_solver

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
