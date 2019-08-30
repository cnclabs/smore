#ifndef PAIRWISE_OPTIMIZER_H
#define PAIRWISE_OPTIMIZER_H
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>

#define SIGMOID_TABLE_SIZE 1000
#define MAX_SIGMOID 8.0

class PairwiseOptimizer {
    private:
        void init_sigmoid();

    public:
        // constructor
        PairwiseOptimizer();
        
        // variables
        std::vector<double> cached_sigmoid;

        // functions
        double fast_sigmoid(double value);

        // loss
        void feed_l2_loss(std::vector<double>& embedding, int dimension, std::vector<double>& loss);
        void feed_loglikelihood_loss(std::vector<double>& from_embedding,
                                     std::vector<double>& to_embedding,
                                     double label,
                                     int dimension,
                                     std::vector<double>& from_loss,
                                     std::vector<double>& to_loss);
};
#endif
