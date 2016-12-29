#ifndef ProximityEmbedding_H
#define ProximityEmbedding_H

#include "LINE.h"

/*****
 * ProximityEmbedding
 * **************************************************************/

class ProximityEmbedding: public LINE {

    public:
        
        ProximityEmbedding();
        ~ProximityEmbedding();

        void SaveWeights(string);

        // model function
        void Init(int);
        void Train(int, int, int, double, int);

};


#endif
