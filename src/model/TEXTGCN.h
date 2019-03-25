#ifndef TEXTGCN_H
#define TEXTGCN_H

#include "LINE.h"

/*****
 * TEXTGCN
 * **************************************************************/

class TEXTGCN: public LINE {

    public:
        
        TEXTGCN();
        ~TEXTGCN();

        void SaveWeights(string);
        
        // additional model parameters
        int num_aggregation = 5;

        // model function
        void Init(int);
        void LoadFieldMeta(string);
        void Train(int, int, int, double, double, int);

};


#endif
