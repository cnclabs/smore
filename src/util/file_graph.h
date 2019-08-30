#ifndef BASE_SAMPLER_H
#define BASE_SAMPLER_H
#include <string>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include "hash.h"

#define SAMPLER_MONITOR 10000

class FileGraph {
    /* FileGraph loads file-based data as a graph.
     */
    private:
        int is_directory(std::string path);
        void load_file_stat(std::string path);
        void load_from_edge_list(std::string path, int undirected);
        // TODO: implement other ways to read from grpah files
        //void load_from_adjacency_list(std::string path);

    public:
        // file-related variables
        std::vector<std::string> file_names;
        std::vector<unsigned long long> file_lines;

        // graph-related variables
        long max_index, num_edge;
        Hash vertex2index; // index map
        std::unordered_map<long, char*> index2vertex; // index map
        //std::unordered_map<long, std::unordered_map<long, double>> graph;
        std::unordered_map<long, std::vector<long>> adjacency;
        std::unordered_map<long, std::vector<double>> weight;

        // constuctor
        FileGraph(std::string path, int undirected);
};
#endif
