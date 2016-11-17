#ifndef DEEPWALK_H
#define DEEPWALK_H

#include "proRec.h"
#include "proNet.h"

/*****
 * DeepWalk
 * **************************************************************/

class DeepWalk {

    public:
        
        DeepWalk();
        ~DeepWalk();

        proNet rgraph;

        // parameters
        int dim;                // representation dimensions
        //unordered_map<long, double*> w_vertex;
        //unordered_map<long, double*> w_context;
        //vector< double* > w_vertex;
        //vector< double* > w_context;
        //double** w_vertex;
        //double** w_context;
        vector< vector<double> > w_vertex;
        vector< vector<double> > w_context;
     
        // data function
        void LoadEdgeList(string, bool);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Update(vector<long>&, vector<long>&, int, double);
        void Train(int, int, int, int, double, int);

};

#endif
