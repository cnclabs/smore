#include "edge_sampler.h"

EdgeSampler::EdgeSampler(FileGraph* file_graph) {
    this->build_from_file_graph(file_graph);
}

void EdgeSampler::build_from_file_graph(FileGraph* file_graph) {
    std::cout << "Build Sampler:" << std::endl;
    long from_index, to_index, to_offset;
    int branch_size, current_index=0;
    double weight;
    std::vector<double> edge_distribution, neg_distribution;
    this->edge_size = file_graph->num_edge;
    this->vertex_size = file_graph->index2vertex.size();

    // build vertex alias
    std::cout << "\tbuild alias method" << std::endl;
    edge_distribution.resize(this->edge_size, 0.0);
    neg_distribution.resize(this->vertex_size, 0.0);
    for (from_index=0; from_index<this->vertex_size; from_index++)
    {
        branch_size = file_graph->adjacency[from_index].size();
        for (to_offset=0; to_offset<branch_size; to_offset++)
        {
            to_index = file_graph->adjacency[from_index][to_offset];
            weight = file_graph->weight[from_index][to_offset];
            edge_distribution[current_index] = weight;
            neg_distribution[to_index] += weight;
        }
    }
    this->edge_sampler.build(edge_distribution, 1.0);
    this->negative_sampler.build(neg_distribution, 0.75);
    std::cout << "\tdone" << std::endl;
}

long EdgeSampler::draw_an_edge() {
    double rand_prob = random_range(0, 1);
    long rand_index = random_range(0, this->edge_size);

    if (rand_prob < this->edge_sampler.alias_probability[rand_index])
        return rand_index;
    else
        return this->edge_sampler.alias_position[rand_index];
}

long EdgeSampler::draw_a_negative() {
    double rand_prob = random_range(0, 1);
    long rand_index = random_range(0, this->vertex_size);

    if (rand_prob < this->negative_sampler.alias_probability[rand_index])
        return rand_index;
    else
        return this->negative_sampler.alias_position[rand_index];
}

