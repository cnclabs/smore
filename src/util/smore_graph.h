#ifndef SMORE_GRAPH_H
#define SMORE_GRAPH_H

#include <vector>
#include <string>
#include <string.h>
#include <unordered_map>
#include <dirent.h>
#include <iostream>
#include "../util/util.h"

#define MONITOR 10000

typedef std::unordered_map<std::string, std::unordered_map<std::string, double>> string_graph;
typedef std::unordered_map<long, std::unordered_map<long, double>> index_graph;

class GraphReader {
    /* Custom file reader
     * Read from a path and transform received data into a 2d-map graph
     */
    private:
        // functions
        void get_file_names(std::string);

    public:
        // constructor
        Reader();
        
        // variables
        std::vector<std::string> file_names;
        std::vector<unsigned long long> file_lines;

        // functions
        smore_graph* smore_graph_from_edge_list(std::string, bool);
        smore_graph* smore_graph_from_adjacency_list(std::string);
};
#endif
