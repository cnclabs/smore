#ifndef PE_H
#define PE_H

#include "LINE.h"

/*****
 * ProximityEmbedding
 * **************************************************************/

class PE: public LINE {

    public:
        
        PE();
        ~PE();

        void LoadWalkMeta(string);
        
        void SaveWeights(string);

        // model function
        void Init(int);
        void Train(int, int, int, double, double, int);

};


#endif
