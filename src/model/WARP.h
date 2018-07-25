#ifndef WARP_H
#define WARP_H

#include "../proNet.h"

/*****
 * WARP
 * **************************************************************/

class WARP {

    public:
        
        WARP();
        ~WARP();
        
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
