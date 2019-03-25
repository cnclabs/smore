CC= g++
CPPFLAGS= -std=c++11 -fopenmp -lm -Ofast
OBJS= deepwalk line walklets hpe app mf bpr hoprec warp nemf nerank eco gcn textgcn
LIBS= -L ../bin -lproNet

all: $(OBJS)

$(OBJS):
	$(CC) $(CPPFLAGS) $@.cpp $(LIBS) -o $@
clean:
	rm $(OBJS)
