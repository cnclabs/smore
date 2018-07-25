#ifndef HBPR_H
#define HBPR_H

#include "../proNet.h"

/*****
 * HBPR
 * **************************************************************/

class HBPR {

    public:
        
        HBPR();
        ~HBPR();
        
        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        vector< vector<double> > w_vertex;
        vector< vector<double> > w_context;

        // data function
        void LoadEdgeList(string, bool);
        void LoadFieldMeta(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, double, int);

};


#endif
