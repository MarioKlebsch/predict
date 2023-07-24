.PHONY: all clean

LDLIBS=-lncurses -ltinfo -lm
OBJS=SDP4.o SGP4.o deep.o more_math.o predict.o

all: predict

predict: $(OBJS)


clean:
	-rm -rf *.o predict
