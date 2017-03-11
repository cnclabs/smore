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
        vector< vector<double> > w_context;

        // data function
        void LoadFieldMeta(string);
        void LoadEdgeList(string, bool);
        void SaveFirstWeights(string);
        void SaveSecondWeights(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, double, int);

};


#endif
