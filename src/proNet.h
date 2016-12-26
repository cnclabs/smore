#ifndef PRONET_H
#define PRONET_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <thread>
#include <string>
#include <string.h>
#include <vector>
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
        //int alias=-1;
        //double prob=0.0;
        //int neg_alias=-1;
        //double neg_prob=0.0;
};

class Field {
    public:
        vector<int> fields;
        vector<int> vids;
};

class Context {
    public:
        int vid=-1;
        double in_degree=0.0;
        //int alias=-1;
        //double prob=0.0;
};

class AliasTable {
    public:
        int alias=-1;
        double prob=0.0;
};

class proNet {

    public:
        
        proNet();
        ~proNet();
        
        // MAX index number
        long long int MAX_line=0;
        int MAX_vid=0;
        int MAX_fvid=0;
        int MAX_field=0;
       
        // graph basics
        vector< int > hash_table;
        vector< char* > keys;
        map< char*, int, cmp_char > kmap;     // vertex map to index number

        // Alias Graph
        vector< Vertex > vertex;
        vector< Context > context;
        vector< AliasTable > vertex_AT;
        vector< AliasTable > context_AT;
        vector< AliasTable > negative_AT;
        vector< Field > field;
        
        // cahce
        void BuildAliasMethod(unordered_map<long, vector<long>>&, unordered_map<long, vector<double>>&);
        void BuildNegativeAliasTable();
        void BuildSourceAliasTable();
        void BuildTargetAliasTable();
        
        // Key Process
        unsigned int BKDRHash(char*);
        int InsertHashTable(char*);
        int SearchHashTable(char*);

        // Data Process
        void LoadEdgeList(string, bool);
        void LoadFieldMeta(string);

        // Network Process
        long SourceSample();
        long TargetSample();
        long TargetSample(long);
        long NegativeSample();
        long NegativeFieldSample(long);
        vector< long > RandomWalk(long, int);
        vector< vector< long > > SkipGrams(vector<long>&, int, int);
        vector< vector< long > > ScaleSkipGrams(vector<long>&, int, int, int);

        // Optimizer
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdatePair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateDirectedPair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
       
        // vertex vector, context vector, vertex series, context series, dimension, negative samples, alpha
        void UpdatePairs(vector< vector<double> >&, vector< vector<double> >&, vector<long>&, vector<long>&, int, int, double);
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, community walk steps, alpha
        void UpdateCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateFieldCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);

        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdateFieldsCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);

};

#endif
