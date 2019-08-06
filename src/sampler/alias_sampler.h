#ifndef ALIAS_SAMPLER_H
#define ALIAS_SAMPLER_H
#include <vector>
#include <cmath>
#include "../util/random.h"

class AliasSampler {
    /* AliasSampler is an implementation of Alias Method for weighted sampling purpose
     * It draws a sample from a weighted distribution
     * Reference: https://en.wikipedia.org/wiki/Alias_method
     */
    public:
        // constuctor
        AliasSampler();
        AliasSampler(std::vector<double>& distribution, const double power);

        // variables
        long distribution_size = 0.0;
        std::vector<long> alias_position;
        std::vector<long> alias_probability;

        // initialization
        void build(std::vector<double>& distribution, const double power);
        
        // functions
        long draw(); // return a position
};
#endif
