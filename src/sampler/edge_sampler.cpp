#include <vector>
#include <cmath>
#include "../util/random.h"

class EdgeSampler {
    /* Edge Sampler
     * EdgeSampler performs edge-based sampling
     */
    public:
        EdgeSampler();
        EdgeSampler(std::vector<double>& distribution, const double power);

        // variables
        std::vector<long> alias;
        std::vector<long> probability;

        // operations
        std::vector<long> sample_an_edge();
}
