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

        // data function
        void LoadFieldMeta(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, int, double, int);

};


#endif
