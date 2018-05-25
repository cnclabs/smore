CC = g++
CPPFLAGS = -std=c++11 -fPIC -fopenmp -lm -Ofast -w
OBJECTS = util.o random.o proNet.o ./model/DeepWalk.o ./model/Walklets.o ./model/LINE.o ./model/HPE.o ./model/APP.o ./model/MF.o ./model/msFINE.o ./model/BPR.o ./model/BPRNE.o ./model/GMF.o ./model/HBPR.o ./model/NEMF.o ./model/WARP.o ./model/fHPE.o
#OBJECTS = util.o random.o proNet.o ./model/MF.o
all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproNet.a $(OBJECTS)

clean:
	rm $(OBJECTS)
	rm ../bin/libproNet.a
