#include "vc_sampler.h"

VCSampler::VCSampler(FileGraph* file_graph) {
    this->build_from_file_graph(file_graph);
    this->adjacency = file_graph->adjacency;
}

void VCSampler::build_from_file_graph(FileGraph* file_graph) {
    std::cout << "Build Sampler:" << std::endl;
    long from_index, to_index, to_offset, num_edges=0;
    int branch_size;
    double weight;
    std::vector<double> vertex_distribution, context_distribution, neg_distribution;
    this->vertex_size = file_graph->index2vertex.size();
    this->context_size = 0;

    // build vertex alias
    std::cout << "\tbuild alias method" << std::endl;
    vertex_distribution.resize(this->vertex_size, 0.0);
    neg_distribution.resize(this->vertex_size, 0.0);
    for (from_index=0; from_index<this->vertex_size; from_index++)
    {
        branch_size = file_graph->adjacency[from_index].size();
        context_distribution.resize(branch_size, 0);
        for (to_offset=0; to_offset<branch_size; to_offset++)
        {
            to_index = file_graph->adjacency[from_index][to_offset];
            weight = file_graph->weight[from_index][to_offset];
            vertex_distribution[from_index] += weight;
            context_distribution[to_offset] = weight;
            neg_distribution[to_index] += weight;
            this->context_size++;
        }
        this->context_sampler[from_index] = AliasMethod(context_distribution, 1.0);
    }
    this->vertex_sampler.build(vertex_distribution, 1.0);
    this->negative_sampler.build(neg_distribution, 0.75);
    std::cout << "\tdone" << std::endl;
}

long VCSampler::draw_a_vertex() {
    double rand_prob = random_range(0, 1);
    long rand_index = random_range(0, this->vertex_size);

    if (rand_prob < vertex_sampler.alias_probability[rand_index])
        return rand_index;
    else
        return vertex_sampler.alias_position[rand_index];
}

long VCSampler::draw_a_context(long vertex_index) {
    double rand_prob = random_range(0, 1);
    long rand_index = random_range(0, this->adjacency[vertex_index].size());

    if (rand_prob < this->context_sampler[vertex_index].alias_probability[rand_index])
        return this->adjacency[vertex_index][rand_index];
    else
        return this->adjacency[vertex_index][this->context_sampler[vertex_index].alias_position[rand_index]];
}

long VCSampler::draw_a_negative() {
    double rand_prob = random_range(0, 1);
    long rand_index = random_range(0, this->vertex_size);

    if (rand_prob < this->negative_sampler.alias_probability[rand_index])
        return rand_index;
    else
        return this->negative_sampler.alias_position[rand_index];
}

