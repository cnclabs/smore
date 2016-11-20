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

        // model function
        void Train(int, int, int, double, int);

};


#endif
