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
#include <array>
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

struct cmp_char {
    bool operator()(char const *a, char const *b) {
        return strcmp(a, b) < 0;
    }
};

struct char_cmp {
    bool operator()(char const *a, char const *b) const {
        return strcmp(a, b) == 0;
    }
};

struct bkdr_hash{
    unsigned int operator()(char *key) const {
        unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
        unsigned int hash = 0;
        while(*key)
            hash = hash * seed + (*key++);
        return (hash % HASH_TABLE_SIZE);
    }
};


class Vertex {
    public:
        long offset, branch;
        double out_degree, in_degree;
        Vertex() { out_degree=0.0; in_degree=0.0; }
};

class devField {
    public:
        int field;
        vector<long> offset, brnach;
};

class Context {
    public:
        long vid;
        double in_degree;
        Context() { vid=-1; in_degree=0.0; }
};

class Field {
    public:
        vector<int> fields = {0};
        vector<int> vids;
};

class AliasTable {
    public:
        long alias;
        double prob;
        AliasTable() { alias=-1; prob=0.0;}
};

class HashTable {
    public:
        vector< long > table;
        vector< char* > keys;
};

class HashTable2 {
    public:
        unordered_map<char*, long, bkdr_hash, char_cmp> table;
        unordered_map<char*, long, bkdr_hash, char_cmp>::iterator it;
        long Find(char*);
        long Insert(char*);
};

class proNet {

    private:

        void InitSigmoid();
        //void InitNegTable();
        void BuildAliasMethod(unordered_map<long, vector<long>>&, unordered_map<long, vector<double>>&);
        
    public:
        
        proNet();
        ~proNet();
        
        void SetNegativeMethod(char*);
        void SetVertexMethod(char*);

        // MAX index number
        char vertex_method[20];
        char context_method[20];
        char negative_method[20];
        unsigned long long MAX_line;
        long MAX_vid;
        long MAX_fvid;
        long MAX_field;

        // Alias Graph
        vector< Vertex > vertex;
        vector< Context > context;
        vector< AliasTable > vertex_AT;
        vector< AliasTable > context_AT;
        vector< AliasTable > negative_AT;
        vector< Field > field;
        //vector< long > neg_table;

        // cahce
        vector< double > cached_sigmoid;
        vector<AliasTable> AliasMethod(vector<double>& distribution, double power);

        // Key Process
        HashTable vertex_hash;

		    // Pretrain
		    unordered_map< long, vector< double >> pretrain;
        //HashTable2 vertex_hash;
        unsigned int BKDRHash(char*);
        void InsertHashTable(HashTable&, char*);
        long SearchHashTable(HashTable&, char*);
        
        // Math-related Process
        double fastSigmoid(double);

        // Data Process
        void LoadEdgeList(string, bool);
        void LoadFieldMeta(string);
        
        vector< int > dynamic_walk;
        void LoadWalkMeta(string);
        void LoadPreTrain(string,int);

        // Network Process
        long SourceSample();
        long TargetSample();
        long TargetSample(long);
        long NegativeSample();
        long NegativeFieldSample(long);
        vector< long > RandomWalk(long, int);
        vector< long > JumpingRandomWalk(long, double);
        vector< vector< long > > CBOWs(vector<long>&, int, int);
        vector< vector< long > > SkipGrams(vector<long>&, int, int);
        vector< vector< long > > OrdinalSkipGrams(vector<long>&);
        vector< vector< long > > ScaleSkipGrams(vector<long>&, int, int, int);

        // Optimizer

        // vertex representation, context representation, label, alpha, regularization, vertex loss, context loss, alpha
        void Opt_SGD(vector<double>&, vector<double>&, double, double, double, vector<double>&, vector<double>&);

        // vertex representation, context representation, alpha, vertex loss, context loss, alpha
        void Opt_BPRSGD(vector<double>&, vector<double>&, double, vector<double>&, vector<double>&);
        int Opt_FBPRSGD(vector<double>&, vector<double>&, double, vector<double>&, vector<double>&, double);

        // vertex representation, context representation, label, alpha, vertex loss, context loss, alpha
        void Opt_SigmoidSGD(vector<double>&, vector<double>&, double, int, double, vector<double>&, vector<double>&);
        void Opt_CosineSGD(vector<double>&, vector<double>&, double, int, double, vector<double>&, vector<double>&);
        void Opt_LengthSGD(vector<double>&, vector<double>&, double, int, double, vector<double>&, vector<double>&);
        void Opt_SigmoidSGD1(double*, double*, double, int, double, double*, double*);

        // vertex representation, context representation, label, alpha, regularization, vertex loss, context loss, alpha
        void Opt_SigmoidRegSGD(vector<double>&, vector<double>&, double, double, double, vector<double>&, vector<double>&);

        // Drawing Pairs
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdatePair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
        void UpdateCosinePair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
        void UpdateLengthPair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
        void UpdatePair1(double*, double*, long, long, int, int, double);
        void UpdateFreezePair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateBPRPair(vector< vector<double> >&, vector< vector<double> >&, long, long, long, int, double, double);
        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateWARPPair(vector< vector<double> >&, vector< vector<double> >&, long, long, long, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateCBPRPair(vector< vector<double> >&, vector< vector<double> >&, long, long, long, int, double, double);
        // vertex vector, context vector, vertex, context, dimension, negative samples, margin
        void UpdateFBPRPair(vector< vector<double> >&, vector< vector<double> >&, long, long, long, int, double, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateBPRPairs(vector< vector<double> >&, vector< vector<double> >&, vector<long>&, vector<long>&, vector<long>&, int, double, double);

        // vertex vector, context vector, vertex, context, dimension, regularization, negative samples, alpha
        void UpdateFactorizedPair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, double);
        void UpdateChoice(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, double);
        void UpdateRAWChoice(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, double);
        void UpdateGroupingPair(vector< vector<double> >&, vector< vector<double> >&, long, long, double, int, double, int, double);

        // vertex vector, context vector, vertex, context, dimension, regularization, walk steps, negative samples, alpha
        void UpdateCBOW(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double);
        void UpdateCBOWs(vector< vector<double> >&, vector< vector<double> >&, vector<long>&, vector<long>&, int, double, int, int, double);

        // vertex vector, context vector, vertex series, context series, dimension, negative samples, alpha
        void UpdatePairs(vector< vector<double> >&, vector< vector<double> >&, vector<long>&, vector<long>&, int, int, double);

        // user vertex vector, item vertex vector, context vector, vertex, context, dimension, regularization, negative samples, community walk steps, alpha
        void UpdateUIPair(vector< vector<double> >&, vector< vector<double> >&, vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double);
      
        // vertex vector, context vector, vertex, context, dimension, regularization, negative samples, community walk steps, alpha
        void UpdateCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double);
        void UpdateFCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, regularization, negative samples, community walk steps, alpha
        void UpdateBatchCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double);
        // vertex vector, context vector, vertex, context, dimension, regularization, negative samples, community walk steps, alpha, u or i
        void UpdateUICommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, double, int, int, double, int);

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
