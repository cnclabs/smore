#ifndef VC_SAMPLER_H
#define VC_SAMPLER_H

#include <cmath>
#include <unordered_map>
#include <vector>
#include "../util/hash.h"
#include "../util/random.h"
#include "../util/smore_graph.h"
#include "alias_sampler.h"

class VCSampler {
    /* VCSampler performs vertex-context-style sampling
     */
    public:
        VCSampler();

        // variables
        long vertex_size=0, context_size=0;
        std::unordered_map<std::string, long index> 
        AliasSampler vertex_sampler;
        std::vector<AliasSampler> context_sampler;

        // initialization
        void build_from_edge_list(std::string path, bool undirected);
        void build_from_adjacency_list(std::string path, bool undirected);
        void build_from_stream_list(std::string path, bool undirected);

        // functions
        long draw_a_vertex();
        long draw_a_context(long vertex);
        long draw_a_negative();
        std::vector<long> draw_a_walk(long vertex, const int walk_steps);
        std::vector<std::vector<long>> draw_skipgram_pairs(long vertex, const int walk_steps, const int window_size);
};
#endif
