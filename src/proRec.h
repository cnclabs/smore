#ifndef PROREC_H
#define PROREC_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <thread>
#include <string>
#include <string.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <utility>
#include <omp.h>

using namespace std;

#define MONITOR 10000
#define HASH_TABLE_SIZE 30000000

double random_gen(const int&, const int&);

struct cmp_char
{
    bool operator()(char const *a, char const *b)
    {
        return strcmp(a, b) < 0;
    }
};

class Vertex {
    public:
        int offset, branch;
        double out_degree=0.0;
        double in_degree=0.0;
        int alias=-1;
        double prob=0.0;
        int neg_alias=-1;
        double neg_prob=0.0;
        double subsample=0.0;
};

class Field {
    public:
        int field=-1;
        vector<int> vids;
};

class Context {
    public:
        int vid=-1;
        double in_degree=0.0;
        int alias=-1;
        double prob=0.0;
};

#endif
