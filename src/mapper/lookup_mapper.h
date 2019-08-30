#ifndef LOOKUP_MAPPER_H
#define LOOKUP_MAPPER_H
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class LookupMapper {
    public:
        //variable
        int size, dimension;
        std::vector<std::vector<double>> embedding;

        // constructor
        LookupMapper(int size, int dimension);
        
        // update function
        void update(long index, std::vector<double>& loss_vector, double alpha);
        void update_with_l2(long index, std::vector<double>& loss_vector, double alpha, double lambda);
        
        // overload operator
        void save_to_file(std::unordered_map<long, char*>& index2vertex, std::string file_name);

        // overload operator
        std::vector<double>& operator[](long index);
};
#endif
