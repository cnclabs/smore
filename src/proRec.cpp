#ifndef PROREC_H
#include "proRec.h"
#endif

double random_gen(const int & min, const int & max) {
    thread_local mt19937 generator(clock());
    uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}

