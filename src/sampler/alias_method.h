#ifndef ALIAS_METHOD_H
#define ALIAS_METHOD_H
#include <vector>
#include <cmath>
#include "../util/random.h"

class AliasMethod {
    /* AliasMethod is an efficient implementation of weighted sampling
     * Reference: https://en.wikipedia.org/wiki/Alias_method
     */
    public:
        // constuctor
        AliasMethod();
        AliasMethod(std::vector<double>& distribution, const double power);

        // variables
        long distribution_size = 0.0;
        std::vector<long> alias_position;
        std::vector<double> alias_probability;

        // initialization
        void build(std::vector<double>& distribution, const double power);
        
        // functions
        long draw(); // return a position
};
#endif
