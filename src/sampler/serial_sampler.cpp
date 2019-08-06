#include <vector>
#include <cmath>
#include "../util/random.h"

class VCASampler {
    /* Vertex-Context Alias (VCA) Sampler
     * VCASampler performs vertex-context-style sampling
     */
    public:
        VCASampler();
        VCASampler(std::vector<double>& distribution, const double power);

        // variables
        std::vector<long> alias;
        std::vector<long> probability;

        // operations
        long sample_a_vertex();
        long sample_a_context(long vertex);
        std::vector<long> sample_a_path(long vertex, const int walk_steps);
        std::vector<std::vector<long>> sample_by_skipgram(long vertex, const int walk_steps, const int window_size);
}
