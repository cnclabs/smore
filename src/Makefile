CC = g++
CPPFLAGS = -std=c++11 -fPIC -fopenmp -lm -Ofast
UTIL_OBJECTS = util/util.o util/smore_graph.o util/hash.o util/random.o
SAMPLER_OBJECTS = sampler/alias_sampler.o sampler/vc_sampler.o

all: $(UTIL_OBJECTS) $(SAMPLER_OBJECTS)
	mkdir -p ../bin
	ar rcs ../bin/libsmore.a $(UTIL_OBJECTS)
	ar rcs ../bin/libsmore.a $(SAMPLER_OBJECTS)

%.o: %.cpp %.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

clean:
	rm $(UTIL_OBJECTS)
	rm ../bin/libsmore.a
