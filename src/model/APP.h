#ifndef APP_H
#define APP_H

#include "../proNet.h"

/*****
 * APP
 * **************************************************************/

class APP {

    public:
        
        APP();
        ~APP();

        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        vector< vector<double> > w_vertex;
        vector< vector<double> > w_context;
     
        // data function
        void LoadEdgeList(string, bool);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void Train(int, int, double, int, double, int);

};

#endif
