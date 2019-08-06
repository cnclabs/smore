#include "smore_graph.h"
#include <string>

Reader::Reader() {
}

void Reader::get_file_names(std::string path) {
    // get file_names
    if (isDirectory(path.c_str()))
    {
        DIR *dir;
        struct dirent *ent;
        dir = opendir(path.c_str());
        while ((ent = readdir (dir)) != NULL) {
            std::string fname = path + "/" + ent->d_name;
            this->file_names.push_back(fname);
        }
        closedir(dir);
    }
    else
    {
        this->file_names.push_back(path.c_str());
    }
    // get total lines
    FILE *fin;
    char c_line[1000];
    unsigned long long num_lines = 0;
    std::cout << "Lines Preview:" << std::endl;
    for (auto fname: this->file_names)
    {
        fin = fopen(fname.c_str(), "rb");
        while (fgets(c_line, sizeof(c_line), fin))
        {
            if (num_lines % MONITOR == 0)
            {
                printf("\t# of lines:\t%lld%c", num_lines, 13);
            }
            ++num_lines;
        }
        fclose(fin);
        this->file_lines.push_back(num_lines);
        printf("\t# of lines:\t%lld%c\n", num_lines, 13);
    }
}

smore_graph* Reader::smore_graph_from_edge_list(std::string path, bool undirected) {
    this->get_file_names(path);
    
    smore_graph* graph = new smore_graph();
    
    std::cout << "Loading Lines:" << std::endl;
    FILE *fin;
    char v1[160], v2[160];
    double weight;
    unsigned long long line = 0;
    for (int i=0; i<this->file_names.size(); i++)
    {
        fin = fopen(this->file_names[i].c_str(), "rb");
        for (; line != this->file_lines[i]; line++)
        {
            if ( fscanf(fin, "%s %s %lf", v1, v2, &weight) != 3 )
            {
                std::cout << "\t[WARNING] skip line " << line << std::endl; 
                continue;
            }
            if (line % MONITOR == 0)
            {
                printf("\tProgress:\t%.2f %%%c", (double)(line)/(this->file_lines[this->file_names.size()-1]) * 100, 13);
                fflush(stdout);
            }
            (*graph)[v1][v2] = weight;
            if (undirected)
            {
                (*graph)[v2][v1] = weight;
            }
        }
    }
    std::cout << "\tProgress:\t100.00 %\r" << std::endl;
    return graph;
}

