#ifndef FINE_H
#define FINE_H

#include "LINE.h"

/*****
 * FINE
 * **************************************************************/

class FINE: public LINE {

    public:
        
        FINE();
        ~FINE();

        int dim_1;
        int dim_2;

        // data function
        void LoadFieldMeta(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, int, double, int);

};


#endif
