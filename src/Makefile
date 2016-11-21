CC = g++-4.9
CPPFLAGS = -std=c++11 -fPIC -fopenmp -O3
OBJECTS = proNet.o DeepWalk.o Walklets.o LINE.o HPE.o FINE.o MF.o

all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproNet.a $(OBJECTS)

clean:
	rm $(OBJECTS)
	rm ../bin/libproNet.a
