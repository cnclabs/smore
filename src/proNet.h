#ifndef PRONET_H
#define PRONET_H

using namespace std;
#if defined (_MSC_VER)  // Visual studio
    #define thread_local __declspec( thread )
#elif defined (__GCC__) // GCC
    #define thread_local __thread
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <random>
#include <string>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <map>
#include <utility>
#include <stdio.h>
#include <omp.h>
#include <dirent.h>
#include "random.h"
#include "util.h"

#define MONITOR 10000
#define POWER_SAMPLE 0.75
#define HASH_TABLE_SIZE 30000000
#define SIGMOID_TABLE_SIZE 1000
#define MAX_SIGMOID 8.0
#define MAX_NEG 1e8

struct cmp_char
{
    bool operator()(char const *a, char const *b)
    {
        return strcmp(a, b) < 0;
    }
};

class Vertex {
    public:
        long offset, branch;
        double out_degree;
        double in_degree;
        Vertex() { out_degree=0.0; in_degree=0.0; }
};

class Field {
    public:
        vector<int> fields = {0};
        vector<int> vids;
};

class Context {
    public:
        long vid;
        double in_degree;
        Context() { vid=-1; in_degree=0.0; }
};

class AliasTable {
    public:
        long alias;
        double prob;
        AliasTable() { alias=-1; prob=0.0;}
};

class proNet {

    public:
        
        proNet();
        ~proNet();
        
        // MAX index number
        unsigned long long MAX_line;
        long MAX_vid;
        long MAX_fvid;
        long MAX_field;
        
        // graph basics
        vector< long > hash_table;
        vector< char* > keys;
        map< char*, long, cmp_char > kmap;     // vertex map to index number

        // Alias Graph
        vector< Vertex > vertex;
        vector< Context > context;
        vector< AliasTable > vertex_AT;
        vector< AliasTable > context_AT;
//        vector< vector<AliasTable> > context_AT;
        vector< AliasTable > negative_AT;
        vector< Field > field;
        vector< long > neg_table;

        // cahce
        vector< double > cached_sigmoid;
        void InitSigmoid();
        void InitNegTable();
        void BuildAliasMethod(unordered_map<long, vector<long>>&, unordered_map<long, vector<double>>&);
        vector<AliasTable> AliasMethod(vector<double>& distribution, double power);

        // Key Process
        unsigned int BKDRHash(char*);
        int InsertHashTable(char*);
        int SearchHashTable(char*);
        
        // Math-related Process
        double fastSigmoid(double);

        // Data Process
        void LoadEdgeList(string, bool);
        void LoadFieldMeta(string);
        
        vector< int > dynamic_walk;
        void LoadWalkMeta(string);

        // Network Process
        long SourceSample();
        long TargetSample();
        long TargetSample(long);
        long NegativeSample();
        long NegativeFieldSample(long);
        //long AliasSample();
        //long AliasSample(vector<AliasTable>&);
        vector< long > RandomWalk(long, int);
        vector< vector< long > > SkipGrams(vector<long>&, int, int);
        vector< vector< long > > ScaleSkipGrams(vector<long>&, int, int, int);

        // Optimizer
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdatePair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateDirectedPair(vector< vector<double> >&, vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
       
        // vertex vector, context vector, vertex series, context series, dimension, negative samples, alpha
        void UpdatePairs(vector< vector<double> >&, vector< vector<double> >&, vector<long>&, vector<long>&, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, regularization, negative samples, community walk steps, alpha
        void UpdateCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, community walk steps, bfs, alpha
        void UpdateDCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateFieldCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, walk_steps, negative samples, alpha
        void UpdateMSFieldCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);
        
        // vertex vector, context vector, vertex, context, dimension_f1, dimension_f2, walk_steps, negative samples, alpha
        //void UpdateMSDiffFieldCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateFieldsCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);

};

#endif
