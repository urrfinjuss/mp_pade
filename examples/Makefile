CC=gcc
FLAGS=-std=gnu99 -I$(HOME)/usr/include -L$(HOME)/usr/lib -Ofast -flto
LFLAGS= -lmppade -lmpfr

SRCS1=main.c
OBJS1=$($SRCS1:.c=.o)
TARGET1=test_pade.x

.PHONY: all
all: $(TARGET1) $(TARGET2)

test_pade.x: main.c
	$(CC) $(FLAGS) -o $(TARGET1) $(SRCS1) $(LFLAGS)
	cp $(TARGET1) ./KdV-soliton/ 
	cp $(TARGET1) ./branch-cut/ 
	cp $(TARGET1) ./three-poles/ 

clean:
	rm -f $(TARGET1) ./KdV-soliton/$(TARGET1) ./branch-cut/$(TARGET1) ./three-poles/$(TARGET1)
