#ifndef MF_H
#define MF_H

#include "../proNet.h"

/*****
 * MF
 * **************************************************************/

class MF {

    public:
        
        MF();
        ~MF();
        
        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        vector< vector<double> > w_vertex;

        // data function
        void LoadEdgeList(string, bool);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, double, double, int);

};


#endif
