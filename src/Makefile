CC = g++-4.8
CPPFLAGS = -std=c++11 -fopenmp -O3
OBJECTS = proRec.o proNet.o DeepWalk.o LINE.o HPE.o

all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproRec.a $(OBJECTS)

clean:
	rm $(OBJECTS)
	rm ../bin/libproRec.a
