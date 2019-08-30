#ifndef VC_SAMPLER_H
#define VC_SAMPLER_H

#include <cmath>
#include <unordered_map>
#include <vector>
#include "../util/file_graph.h"
#include "../util/random.h"
#include "alias_method.h"

class VCSampler {
    /* VCSampler performs vertex-context-style sampling
     */
    private:
        void build_from_file_graph(FileGraph* file_graph);

    public:
        VCSampler(FileGraph* file_graph);

        // variables
        long vertex_size=0, context_size=0;
        AliasMethod vertex_sampler;
        AliasMethod negative_sampler;
        std::unordered_map<long, AliasMethod> context_sampler;
        std::unordered_map<long, std::vector<long>> adjacency; // context ref.

        // functions
        long draw_a_vertex();
        long draw_a_context(long vertex);
        long draw_a_negative();
        std::vector<long> draw_a_vertex_sequence();
        std::vector<long> draw_a_walk(long vertex, const int walk_steps);
        std::vector<std::vector<long>> draw_skipgram_pairs(long vertex, const int walk_steps, const int window_size);
};
#endif
