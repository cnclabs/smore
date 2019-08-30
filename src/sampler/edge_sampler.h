#ifndef EDGE_SAMPLER_H
#define EDGE_SAMPLER_H

#include <cmath>
#include <unordered_map>
#include <vector>
#include "../util/file_graph.h"
#include "../util/random.h"
#include "alias_method.h"

class EdgeSampler {
    /* EdgeSampler performs edge-style sampling
     */
    private:
        void build_from_file_graph(FileGraph* file_graph);

    public:
        EdgeSampler(FileGraph* file_graph);

        // variables
        long vertex_size=0, edge_size=0;
        AliasMethod edge_sampler, negative_sampler;
        std::vector<long> from_edge;
        std::vector<long> to_edge;

        // functions
        long draw_an_edge();
        long draw_a_negative();
};
#endif
