#include "util.h"

int ArgPos(char *str, int argc, char **argv) {
    int a;
    for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
        if (a == argc - 1) {
            printf("Argument missing for %s\n", str);
            exit(1);
        }
        return a;
    }
    return -1;
}

bool isDirectory(std::string path) {
    struct stat statbuf;
    if (stat(path.c_str(), &statbuf) != 0)
        return 0;
    return S_ISDIR(statbuf.st_mode);
}

double dot_similarity(std::vector<double>& embeddingA, std::vector<double>& embeddingB, int dimension) {
    double prediction=0;
    for (int d=0; d<dimension; d++)
    {
        prediction += embeddingA[d]*embeddingB[d];
    }
    return prediction;
}

Monitor::Monitor(unsigned long long total_step) {
    this->total_step = total_step;
}

void Monitor::progress(unsigned long long current_step) {
    printf("\tProgress: %.3f %%%c", (double)current_step/this->total_step*100.0, 13);
}

void Monitor::end() {
    printf("\tProgress: %.3f", 100.0);
}

