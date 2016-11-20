#ifndef LINE_H
#define LINE_H

#include "proNet.h"

/*****
 * LINE
 * **************************************************************/

class LINE {

    public:
        
        LINE();
        ~LINE();
        
        proNet pnet;

        // parameters
        int dim;                // representation dimensions
        vector< vector<double> > w_vertex;
        vector< vector<double> > w_context;

        // data function
        void LoadEdgeList(string, bool);
        void SaveFirstWeights(string);
        void SaveSecondWeights(string);
        void SaveWeights(string);
        
        // model function
        void Init(int);
        void UpdateFirstOrder(long, long, int, double);
        //void UpdateSecondOrder(long, long, int, double);
        void Train(int, int, double, int);

};


#endif
