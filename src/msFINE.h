#ifndef msFINE_H
#define msFINE_H

#include "LINE.h"

/*****
 * msFINE
 * **************************************************************/

class msFINE: public LINE {

    public:
        
        msFINE();
        ~msFINE();

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
