#ifndef SPR_H
#define SPR_H

#include "../proNet.h"

/*****
 * SkewOPT
 * **************************************************************/

class SPR {

    public:

        SPR();
        ~SPR();

        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        vector< vector<double> > w_vertex;

        // data function
        void LoadEdgeList(string, bool);
        void SaveWeights(string);

        // model function
        void Init(int);
        void Train(int, int, double, double, double, double, int, int);

};


#endif
