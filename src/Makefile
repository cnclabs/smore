CC = g++
CPPFLAGS = -std=c++11 -fPIC -fopenmp -lm -Ofast
OBJECTS = util.o random.o proNet.o ./model/DeepWalk.o ./model/Walklets.o ./model/LINE.o ./model/HPE.o ./model/APP.o ./model/MF.o ./model/BPR.o ./model/HBPR.o ./model/NEMF.o ./model/WARP.o ./model/NERANK.o ./model/ECO.o ./model/GCN.o ./model/TEXTGCN.o ./model/TEXTGCNdev.o
all: $(OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libproNet.a $(OBJECTS)

%.o:%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

clean:
	rm $(OBJECTS)
	rm ../bin/libproNet.a
