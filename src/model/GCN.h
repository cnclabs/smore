#ifndef GCN_H
#define GCN_H

#include "LINE.h"

/*****
 * GCN
 * **************************************************************/

class GCN: public LINE {

    public:
        
        GCN();
        ~GCN();

        void SaveWeights(string);
        
        // additional model parameters
        int num_aggregation = 5;

        // model function
        void Init(int);
        void LoadFieldMeta(string);
        void Train(int, int, int, double, double, int);

};


#endif
