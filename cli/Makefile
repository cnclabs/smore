CC= g++-5
CPPCFLAGS= -std=c++11 -fopenmp -lm -O3

cli:
	$(CC) cli.cpp $(CPPCFLAGS) -L../bin -lproNet -lboost_random -o cli
clean:
	rm cli
