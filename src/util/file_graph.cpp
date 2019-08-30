#include "file_graph.h"

FileGraph::FileGraph(std::string path, int undirected) {
    this->load_from_edge_list(path, undirected);
}

int FileGraph::is_directory(std::string path) {
    struct stat info;
    if( stat( path.c_str(), &info ) != 0 )
        return -1;
    else if( info.st_mode & S_IFDIR )
        return 1;
    return 0;
}

void FileGraph::load_file_stat(std::string path) {
    /* Load file names and lines.
     */
    // clear variables
    this->file_names.clear();
    this->file_lines.clear();
    int is_directory = this->is_directory(path);

    if (is_directory==-1) // fail
    {
        std::cout << "cannot access " << path << std::endl;
        exit(0);
    }
    else if (is_directory==1) // folder with multiple files
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
    else // single file
    {
        this->file_names.push_back(path.c_str());
    }
    // get lines
    FILE *fin;
    char c_line[1000];
    unsigned long long num_lines = 0;
    std::cout << "Lines Preview:" << std::endl;
    for (auto fname: this->file_names)
    {
        fin = fopen(fname.c_str(), "rb");
        while (fgets(c_line, sizeof(c_line), fin))
        {
            if (num_lines % SAMPLER_MONITOR == 0)
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

void FileGraph::load_from_edge_list(std::string path, int undirected) {
    this->load_file_stat(path);
    std::cout << "Loading Lines:" << std::endl;
    FILE *fin;
    char from_vertex[256], to_vertex[256];
    long from_index, to_index, max_index=0;
    double weight;
    unsigned long long line = 0;
    // read from file list
    this->num_edge = 0;
    for (int i=0; i<this->file_names.size(); i++)
    {
        fin = fopen(this->file_names[i].c_str(), "rb");
        for (; line != this->file_lines[i]; line++)
        {
            // read from file
            if ( fscanf(fin, "%s %s %lf", from_vertex, to_vertex, &weight) != 3 )
            {
                std::cout << "\t[WARNING] skip line " << line << std::endl; 
                continue;
            }
            if (line % SAMPLER_MONITOR == 0)
            {
                printf("\tProgress:\t%.2f %%%c", (double)(line)/(this->file_lines[this->file_names.size()-1]) * 100, 13);
                fflush(stdout);
            }
            
            // generate index map
            from_index = this->vertex2index.search_key(from_vertex);
            if (from_index==-1)
            {
                this->vertex2index.insert_key(from_vertex);
                from_index = this->vertex2index.search_key(from_vertex);
                this->index2vertex[from_index] = strdup(from_vertex);
            }
            to_index = this->vertex2index.search_key(to_vertex);
            if (to_index==-1)
            {
                this->vertex2index.insert_key(to_vertex);
                to_index = this->vertex2index.search_key(to_vertex);
                this->index2vertex[to_index] = strdup(to_vertex);
            }

            // store in graph
            //this->graph[from_index][to_index] = weight;
            this->adjacency[from_index].push_back(to_index);
            this->weight[from_index].push_back(weight);
            this->num_edge++;
            if (undirected)
            {
                //this->graph[to_index][from_index] = weight;
                this->adjacency[to_index].push_back(from_index);
                this->weight[to_index].push_back(weight);
                this->num_edge++;
            }
        }
    }
    std::cout << "\tProgress:\t100.00 %\r" << std::endl;
    std::cout << "\t# of vertex:\t" << this->index2vertex.size() << std::endl;
}

