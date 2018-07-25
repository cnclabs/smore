#ifndef NERANK_H
#define NERANK_H

#include "../proNet.h"

/*****
 * NERANK
 * **************************************************************/

class NERANK {

    public:
        
        NERANK();
        ~NERANK();
        
        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        vector< vector<double> > w_vertexU;
        vector< vector<double> > w_vertexI;
        vector< vector<double> > w_contextU;
        vector< vector<double> > w_contextI;
        vector< vector<double> > w_contextUU;
        vector< vector<double> > w_contextII;

        // data function
        void LoadEdgeList(string, bool);
        void LoadFieldMeta(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, double, int);

};


#endif
