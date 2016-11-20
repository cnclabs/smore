#ifndef PRONET_H
#define PRONET_H

#include "proRec.h"

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
        vector< Field > field;
        
        // cahce
        void BuildAliasMethod(unordered_map<long, vector<long>>, unordered_map<long, vector<double>>);
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
        long TargetSample(long);
        long NegativeSample();
        vector< long > RandomWalk(long, int);
        vector< vector< long > > SkipGrams(vector<long>&, int, int);
        vector< vector< long > > ScaleSkipGrams(vector<long>&, int, int, int);

        // Optimizer
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, alpha
        void UpdatePair(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, double);
        
        // vertex vector, context vector, vertex series, context series, dimension, negative samples, alpha
        void UpdatePairs(vector< vector<double> >&, vector< vector<double> >&, vector<long>&, vector<long>&, int, int, double);
        
        // vertex vector, context vector, vertex, context, dimension, negative samples, community walk steps, alpha
        void UpdateCommunity(vector< vector<double> >&, vector< vector<double> >&, long, long, int, int, int, double);

};

#endif
