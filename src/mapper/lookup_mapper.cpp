#include "lookup_mapper.h"

LookupMapper::LookupMapper(int size, int dimension) {
    this->embedding.reserve(size);
    this->size = size;
    this->dimension = dimension;
    for (long index=0; index<size; ++index)
    {
        this->embedding[index].reserve(dimension);
        for (int d=0; d<dimension; ++d)
        {
            this->embedding[index][d] = (rand()/(double)RAND_MAX - 0.5);
        }
    }
}

void LookupMapper::update(long index, std::vector<double>& loss_vector, double alpha) {
    for (int d=0; d<this->dimension; d++)
        this->embedding[index][d] += alpha*loss_vector[d];
}

void LookupMapper::update_with_l2(long index, std::vector<double>& loss_vector, double alpha, double lambda) {
    for (int d=0; d<this->dimension; d++)
        this->embedding[index][d] += alpha*(loss_vector[d] - lambda*embedding[index][d]);
}

void LookupMapper::save_to_file(std::unordered_map<long, char*>& index2vertex, std::string file_name) {
    std::cout << "Save Model:" << std::endl;
    std::ofstream embedding_file(file_name);
    if (embedding_file)
    {
        embedding_file << this->size << " " << this->dimension << std::endl;
        for (long index=0; index!=this->size; index++)
        {
            embedding_file << index2vertex[index];
            for (int dim=0; dim!=this->dimension; dim++)
            {
                embedding_file << " " << this->embedding[index][dim];
            }
            embedding_file << std::endl;
        }
        std::cout << "\tSave to <" << file_name << ">" << std::endl;
    }
    else
    {
        std::cout << "\tfail to open file" << std::endl;
    }
}

std::vector<double>& LookupMapper::operator[](long index) {
    return this->embedding[index];
}
