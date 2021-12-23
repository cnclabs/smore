#ifndef HPE_H
#define HPE_H

#include "LINE.h"

/*****
 * HPE
 * **************************************************************/

class HPE: public LINE {

    public:

        HPE();
        ~HPE();

        void SaveWeights(string);

        // model function
        void Init(int);
        void Init(int, char*, char*);
        void Train(int, int, int, double, double, int);

};


#endif
