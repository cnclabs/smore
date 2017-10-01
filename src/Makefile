CC = g++
CPPFLAGS = -std=c++11 -fPIC -fopenmp -lm -O3
OBJECTS = util.o random.o proNet.o ./model/DeepWalk.o ./model/Walklets.o ./model/LINE.o ./model/HPE.o ./model/APP.o ./model/MF.o ./model/ProximityEmbedding.o ./model/msFINE.o ./model/BPR.o

all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproNet.a $(OBJECTS)

clean:
	rm $(OBJECTS)
	rm ../bin/libproNet.a
