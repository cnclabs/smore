#ifndef WALKLETS_H
#define WALKLETS_H

#include "DeepWalk.h"

/*****
 * Walklets
 * **************************************************************/

class Walklets: public DeepWalk {

    public:
        
        Walklets();
        ~Walklets();

        // model function
        void Train(int, int, int, int, int, double, int);

};


#endif
