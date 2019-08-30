#ifndef UTIL_H
#define UTIL_H
#include <sys/stat.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <vector>

int ArgPos(char *str, int argc, char **argv);
bool isDirectory(std::string path);
double dot_similarity(std::vector<double>& embeddingA, std::vector<double>& embeddingB, int dimension);

class Monitor {
    public:
        unsigned long long total_step;

        // constructor
        Monitor(unsigned long long total_step);

        // count
        void progress(unsigned long long current_step);
        void end();
};

#endif
