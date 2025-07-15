CC = gcc
CFLAGS = -I./easel -I./src -Wall -O2
LDFLAGS = -lm

SRC_DIR = src
EASEL_DIR = easel

SRC_FILES = $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES = $(SRC_FILES:.c=.o)

BIN = dmmbuild

all: $(BIN)

# -> build easel first via its Makefile
$(OBJ_FILES): | easel

easel:
	$(MAKE) -C $(EASEL_DIR)

# final binary
$(BIN): $(OBJ_FILES)
	$(CC) $(OBJ_FILES) $(EASEL_DIR)/*.o -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(MAKE) -C $(EASEL_DIR) clean
	rm -f $(OBJ_FILES) $(BIN)

.PHONY: all clean easel
