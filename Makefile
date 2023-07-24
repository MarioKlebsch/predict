.PHONY: all clean

OBJS=SDP4.o SGP4.o deep.o more_math.o predict.o
LDLIBS=-lncurses -ltinfo -lm

all: predict

predict: $(OBJS)


clean:
	-rm -rf *.o predict
