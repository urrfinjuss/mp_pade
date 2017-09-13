CC=gcc
RM=rm -f
LIB_PATH=$(HOME)/usr/lib
INC_PATH=$(HOME)/usr/include
FLAGS=-std=gnu99 -march=native -Ofast -flto -I/home/orange/code/mppade/src/ -I$(INC_PATH)
CFLAGS=-std=gnu99 -fPIC -Wall -Wextra -Wno-unused-result -Ofast -flto -I/home/orange/code/mppade/src/ -I$(INC_PATH)
LDFLAGS=-shared

TARGET=libmppade.so

SOURCES=$(wildcard src/*.c)

HEADERS=$(wildcard src/*.h)

OBJECTS=$(SOURCES:.c=.o)

.PHONY: all
all: ${TARGET} ${TARGET_THREAD}

$(TARGET): $(OBJECTS)
	$(CC) ${FLAGS} ${LDFLAGS} -o $@ $^
	cp $(TARGET) $(LIB_PATH)
	cp $(HEADERS) $(INC_PATH)

$(SOURCES:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@
include $(SOURCES:.c=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET} ${OBJECTS} $(SOURCES:.c=.d)


