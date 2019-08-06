CC= g++
CPPFLAGS= -std=c++11 -fopenmp -lm -Ofast
OBJS= tester
LIBS= -L ../bin -lsmore

all: $(OBJS)

$(OBJS):
	$(CC) $(CPPFLAGS) $@.cpp $(LIBS) -o $@
clean:
	rm $(OBJS)
