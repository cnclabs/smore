#include "pairwise_optimizer.h"

PairwiseOptimizer::PairwiseOptimizer() {
    this->init_sigmoid();
}

void PairwiseOptimizer::init_sigmoid() {
    cached_sigmoid.resize(SIGMOID_TABLE_SIZE);
    for (int i = 0; i != SIGMOID_TABLE_SIZE + 1; i++)
    {
        double x = i * 2.0 * MAX_SIGMOID / SIGMOID_TABLE_SIZE - MAX_SIGMOID;
        cached_sigmoid[i] = 1.0 / (1.0 + exp(-x));
    }
}

double PairwiseOptimizer::fast_sigmoid(double value) {
    if (value < -MAX_SIGMOID)
    {
        return 0.0;
    }
    else if (value > MAX_SIGMOID)
    {
        return 1.0;
    }
    else
    {
        return this->cached_sigmoid[ int((value + MAX_SIGMOID) * SIGMOID_TABLE_SIZE / MAX_SIGMOID / 2) ];
    }
}

void PairwiseOptimizer::feed_loglikelihood_loss(std::vector<double>& from_embedding, std::vector<double>& to_embedding, double label, int dimension, std::vector<double>& from_loss, std::vector<double>& to_loss) {
    double gradient, prediction=0;
    for (int d=0; d<dimension;d++)
    {
        prediction += from_embedding[d] * to_embedding[d];
    }
    gradient = label - this->fast_sigmoid(prediction);
    for (int d=0; d<dimension; ++d)
    {
        from_loss[d] += gradient * to_embedding[d];
        to_loss[d] += gradient * from_embedding[d];
    }
}
