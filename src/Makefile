CC = g++-5
CPPFLAGS = -std=c++11 -fPIC -fopenmp -lm -O3
OBJECTS = proNet.o ./model/DeepWalk.o ./model/Walklets.o ./model/LINE.o ./model/HPE.o ./model/FINE.o ./model/msFINE.o ./model/MF.o ./model/ProximityEmbedding.o

all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproNet.a $(OBJECTS)

clean:
	rm $(OBJECTS)
	rm ../bin/libproNet.a
