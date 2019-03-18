#ifndef ECO_H
#define ECO_H

#include "../proNet.h"

/*****
 * ECO
 * **************************************************************/

class ECO {

    public:
        
        ECO();
        ~ECO();
        
        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        int dimf;                // representation dimensions
        vector< vector<double> > w_vertex;
        vector< vector<double> > w_ignore;

        // data function
        void LoadEdgeList(string, bool);
        void LoadFieldMeta(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, double, double, int);

};


#endif
