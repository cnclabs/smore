CC = g++-4.9
CPPFLAGS = -std=c++11 -fPIC -fopenmp -O3
OBJECTS = proRec.o proNet.o DeepWalk.o Walklets.o LINE.o HPE.o FINE.o

all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproRec.a $(OBJECTS)

clean:
	rm $(OBJECTS)
	rm ../bin/libproRec.a
