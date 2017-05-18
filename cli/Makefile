CC= g++-5
CPPFLAGS= -std=c++11 -fopenmp -lm -O3

cli:
	$(CC) $(CPPFLAGS) cli.cpp -L../bin -lproNet -lboost_random -o cli
clean:
	rm cli
