#include "proNet.h"

double random_gen(const int & min, const int & max) {
    thread_local mt19937 generator(clock());
    uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}

bool isDirectory(const char *path) {
    struct stat statbuf;
    if (stat(path, &statbuf) != 0)
        return 0;
    return S_ISDIR(statbuf.st_mode);
}

proNet::proNet() {
    hash_table.resize(HASH_TABLE_SIZE, -1);
    InitSigmoid();
}

proNet::~proNet() {
}

unsigned int proNet::BKDRHash(char *key) {
    
    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int hash = 0;
    while (*key)
    {
        hash = hash * seed + (*key++);
    }
    return (hash % HASH_TABLE_SIZE);
}

void proNet::InitSigmoid() {

    cached_sigmoid.resize(SIGMOID_TABLE_SIZE);
    for (int i = 0; i < SIGMOID_TABLE_SIZE + 1; i++) {
        double x = double(i * 2 * MAX_SIGMOID) / SIGMOID_TABLE_SIZE - MAX_SIGMOID;
        cached_sigmoid[i] = 1.0 / (1.0 + exp(-x));
    }
}

double proNet::fastSigmoid(double x) {
    if (x < -MAX_SIGMOID) {
        return 0.0;
    } else if (x > MAX_SIGMOID) {
        return 1.0;
    } else {
        return cached_sigmoid[ int((x + MAX_SIGMOID) * SIGMOID_TABLE_SIZE / MAX_SIGMOID / 2) ];
    }
}

int proNet::InsertHashTable(char *key)
{
    unsigned int pos = BKDRHash(key);
    while (hash_table[pos] != -1)
        pos = (pos + 1) % HASH_TABLE_SIZE;
    hash_table[pos] = MAX_vid;
    kmap[ strdup(key) ] = MAX_vid;
    keys.push_back( strdup(key) );
    MAX_vid++;

    return MAX_vid-1;
}

int proNet::SearchHashTable(char *key)
{
    unsigned int pos = BKDRHash(key);
    while (1)
    {
        if (hash_table[pos] == -1)
            return -1;
        if ( !strcmp(key, keys[ hash_table[pos] ]) )
            return hash_table[pos];
        pos = (pos + 1) % HASH_TABLE_SIZE;
    }
}


void proNet::LoadEdgeList(string filename, bool undirect) {

    // calculate the total connections
    FILE *fin;
    char c_line[1000];
    vector< string > filenames;
    vector< int > filelines;
    
    // load from a folder or from a file
    if (isDirectory(filename.c_str()))
    {
        DIR *dir;
        struct dirent *ent;
        dir = opendir(filename.c_str());
        while ((ent = readdir (dir)) != NULL) {
            string fname = filename + "/" + ent->d_name;
            filenames.push_back(fname);
        }
        closedir(dir);
    }
    else
    {
        filenames.push_back(filename.c_str());
    }
    
    cout << "Connections Preview:" << endl; 
    for (auto fname: filenames)
    {
        fin = fopen(fname.c_str(), "rb");
        while (fgets(c_line, sizeof(c_line), fin))
        {
            if (MAX_line % MONITOR == 0)
            {
                printf("\t# of connection:\t%lld%c", MAX_line, 13);
            }
            ++MAX_line;
        }
        fclose(fin);
        filelines.push_back(MAX_line);
    }
    cout << "\t# of connection:\t" << MAX_line << endl;

    // load the connections
    char v1[160], v2[160];
    int vid1, vid2;
    double w;
    vector< int > v_in, v_out;
    vector< double > e_w;
    if (undirect)
    {
        v_in.resize(MAX_line*2);
        v_out.resize(MAX_line*2);
        e_w.resize(MAX_line*2);
    }
    else
    {
        v_in.resize(MAX_line);
        v_out.resize(MAX_line);
        e_w.resize(MAX_line);
    }

    cout << "Connections Loading:" << endl;
    long long int line = 0;
    for (int i=0; i<filenames.size();i++)
    {
        fin = fopen(filenames[i].c_str(), "rb");
        for (; line != filelines[i]; line++)
        {
            if ( fscanf(fin, "%s %s %lf", v1, v2, &w) != 3 )
            {
                cout << "\t[ERROR] line " << line << " contains wrong number of column data" << endl; 
                continue;
            }

            // generate keys lookup table (kmap)
            vid1 = SearchHashTable(v1);
            if (vid1 == -1)
            {
                vid1 = InsertHashTable(v1);
            }
            vid2 = SearchHashTable(v2);
            if (vid2 == -1)
            {
                vid2 = InsertHashTable(v2);
            }

            v_in[line] = vid1;
            v_out[line] = vid2;
            e_w[line] = w;

            if (undirect)
            {
                v_in[line] = vid2;
                v_out[line] = vid1;
                e_w[line] = w;
            }

            if (line % MONITOR == 0)
            {
                printf("\tProgress:\t\t%.2f %%%c", (double)(line)/(MAX_line+1) * 100, 13);
                fflush(stdout);
            }
        }
        fclose(fin);
    }
    cout << "\tProgress:\t\t100.00 %\r" << endl;
    cout << "\t# of vertex:\t\t" << MAX_vid << endl;

    unordered_map< long, vector< long > > graph;
    unordered_map< long, vector< double > > edge;

    // convert the connections to a graph structure
    cout << "Graph Re-constructing:" << endl;
    for (long long int line=0; line!=MAX_line; line++)
    {
        vid1 = v_in[line];
        vid2 = v_out[line];
        w = e_w[line];

        graph[ vid1 ].push_back( vid2 );
        edge[ vid1 ].push_back( w );

        if (undirect)
        {
            graph[ vid2 ].push_back( vid1 );
            edge[ vid2 ].push_back( w );
        }

        if (line % MONITOR == 0)
        {
            printf("\tProgress:\t\t%.2f %%%c", (double)(line)/(MAX_line+1) * 100, 13);
            fflush(stdout);
        }
    }
    cout << "\tProgress:\t\t100.00 %\r" << endl;

    // release the occupied memory
    v_in = vector<int>();
    v_out = vector<int>();
    e_w = vector<double>();

    cout << "Build the Alias Method:" << endl;
    if (undirect)
        MAX_line *= 2;
    BuildAliasMethod(graph, edge);
    cout << "\tFinished." << endl;

}


void proNet::LoadFieldMeta(string filename) {

    // calculate the # of meta data
    FILE *fin;
    char c_line[1000];
    int max_line=0;
    
    field.resize(MAX_vid);
    MAX_fvid = MAX_vid;

    cout << "Meta Data Preview:" << endl;
    fin = fopen(filename.c_str(), "rb");
    while (fgets(c_line, sizeof(c_line), fin))
    {
        if (max_line % MONITOR == 0)
        {
            printf("\t# of meta data:\t\t%d%c", max_line, 13);
        }
        ++max_line;
    }
    fclose(fin);
    cout << "\t# of meta data:\t\t" << max_line << endl;
    
    char v[160], meta[160];
    int vid;
    map< char*, int, cmp_char > meta_idx;

    cout << "Meta Data Loading:" << endl;
    fin = fopen(filename.c_str(), "rb");
    for (int line = 0; line != max_line; line++)
    {
        if ( fscanf(fin, "%s %s", v, meta)!=2 )
        {
            cout << "line " << line << " contains wrong number of data" << endl; 
            continue;
        }
        
        // generate keys lookup table (meta_idx)
        if (meta_idx.find(meta) == meta_idx.end())
        {
            meta_idx[ strdup(meta) ] = MAX_field;
            MAX_field++;
        }
        vid = kmap[strdup(v)];
        if (vid != -1)
            field[ vid ].fields.push_back(meta_idx[meta]);
        else
            cout << "vertex " << v << " does not appear in given network" << endl;
       
        if (line % MONITOR == 0)
        {
            printf("\tProgress:\t\t%.2f %%%c", (double)(line)/(max_line+1) * 100, 13);
            fflush(stdout);
        }
    }
    fclose(fin);
    cout << "\tProgress:\t\t100.00 %\r" << endl;
    cout << "\t# of field:\t\t" << MAX_field << endl;

    cout << "Init Field Index:" << endl;
    for (int vid=0; vid<MAX_vid; vid++)
    {
        field[ vid ].vids.resize(MAX_field);
        for (int i=0; i<MAX_field; i++)
        {
            if (i == field[vid].fields[0])
            {
                field[vid].vids[i] = vid;
            }
            else
            {
                field[vid].vids[i] = MAX_fvid;
                MAX_fvid++;
            }
        }
    }
    cout << "\tFinished." << endl;

}

void proNet::BuildAliasMethod(unordered_map< long, vector< long > > &graph, unordered_map< long, vector< double > > &edge) {

    // re-construct the graph
    // source === (weight) === target
    
    cout << "\tInitializing ..." << endl;
    
    vertex.resize(MAX_vid);
    context.resize(MAX_line);
    
    int vid;
    long long int line_g=0, line_e=0;
    long offset = 0;
    for (int v1=0; v1!=MAX_vid; v1++)
    {
        vertex[v1].offset = offset;
        vertex[v1].branch = graph[v1].size();
        offset += graph[v1].size();

        for (auto v2: graph[v1])
        {
            context[line_g].vid = v2;
            line_g++;
        }
        for (auto w: edge[v1])
        {
            vertex[v1].out_degree += w;
            context[line_e].in_degree = w;
            line_e++;
        }
    }

    for (long long int line=0; line!=MAX_line; line++)
    {
        vid = context[line].vid;
        vertex[vid].in_degree += context[line].in_degree;
    }
    
    graph.clear();
    edge.clear();

    // compute alias table
    cout << "\tAlias Table Constructing ..." << endl;
    BuildNegativeAliasTable();
    BuildSourceAliasTable();
    BuildTargetAliasTable();

}

void proNet::BuildNegativeAliasTable() {

    // normalization of vertices weights
    double sum, norm;
    vector <double> norm_prob;
    negative_AT.resize(MAX_vid);
    
    sum = 0;
    for (long v1=0; v1!=MAX_vid; v1++)
    {
        sum += pow((vertex[v1].in_degree+vertex[v1].out_degree), 0.75);
        //sum += pow((vertex[v1].in_degree), 0.75);
        //sum += vertex[v1].in_degree;
        //sum += vertex[v1].in_degree+vertex[v1].out_degree;
    }
    norm = MAX_vid/sum;

    for (long v1=0; v1!=MAX_vid; v1++)
    {
        norm_prob.push_back( pow((vertex[v1].in_degree+vertex[v1].out_degree), 0.75)*norm );
        //norm_prob.push_back( pow((vertex[v1].in_degree), 0.75)*norm );
        //norm_prob.push_back( vertex[v1].in_degree*norm );
        //norm_prob.push_back( (vertex[v1].in_degree+vertex[v1].out_degree)*norm );
    }
 
    // block divison
    vector <long> small_block, large_block;
    long num_small_block = 0, num_large_block = 0;
    
    for (long v1=0; v1!=MAX_vid; v1++)
    {
        if ( norm_prob[v1]<1 )
        {
            small_block.push_back( v1 );
            num_small_block++;
        }
        else
        {
            large_block.push_back( v1 );
            num_large_block++;
        }
    }

    // assign alias table
    long small_v, large_v;

    while (small_block.size() && large_block.size())
    {
        small_v = small_block.back();
        small_block.pop_back();
        large_v = large_block.back();
        large_block.pop_back();

        negative_AT[small_v].alias = large_v;
        negative_AT[small_v].prob = norm_prob[small_v];
        norm_prob[large_v] = norm_prob[large_v] + norm_prob[small_v] - 1;
        if (norm_prob[large_v] < 1)
        {
            small_block.push_back( large_v );
        }
        else
        {
            large_block.push_back( large_v );
        }
    }

    while (large_block.size())
    {
        large_v = large_block.back();
        large_block.pop_back();
        negative_AT[large_v].prob = 1.0;
    }

    while (small_block.size())
    {
        small_v = small_block.back();
        small_block.pop_back();
        negative_AT[small_v].prob = 1.0;
    }

}

void proNet::BuildSourceAliasTable() {

    // normalization of vertices weights
    double sum, norm;
    vector <double> norm_prob;
    vertex_AT.resize(MAX_vid);
    
    sum = 0;
    for (long v1=0; v1<MAX_vid; v1++)
    {
        sum += vertex[v1].out_degree;
    }
    norm = MAX_vid/sum;

    for (long v1=0; v1<MAX_vid; v1++)
    {
        norm_prob.push_back( vertex[v1].out_degree*norm );
    }
 
    // block divison
    vector <long> small_block, large_block;
    long num_small_block = 0, num_large_block = 0;
    
    for (long v1=0; v1<MAX_vid; v1++)
    {
        if ( norm_prob[v1]<1 )
        {
            small_block.push_back( v1 );
            num_small_block++;
        }
        else
        {
            large_block.push_back( v1 );
            num_large_block++;
        }
    }
    
    // assign alias table
    long small_v, large_v;

    while (small_block.size() && large_block.size())
    {
        small_v = small_block.back();
        small_block.pop_back();
        large_v = large_block.back();
        large_block.pop_back();

        vertex_AT[small_v].alias = large_v;
        vertex_AT[small_v].prob = norm_prob[small_v];
        norm_prob[large_v] = norm_prob[large_v] + norm_prob[small_v] - 1;
        if (norm_prob[large_v] < 1)
        {
            small_block.push_back( large_v );
        }
        else
        {
            large_block.push_back( large_v );
        }
    }

    while (large_block.size())
    {
        large_v = large_block.back();
        large_block.pop_back();
        vertex_AT[large_v].prob = 1.0;
    }

    while (small_block.size())
    {
        small_v = small_block.back();
        small_block.pop_back();
        vertex_AT[small_v].prob = 1.0;
    }

}

void proNet::BuildTargetAliasTable() {
    
    context_AT.resize(MAX_line);

    for (int vid=0; vid<MAX_vid;vid++)
    {
        // normalization of vertices weights
        long offset, branch;
        offset = vertex[vid].offset;
        branch = vertex[vid].branch;

        double sum, norm;
        vector <double> norm_prob;
        sum = 0;
        for (int i=0; i<branch; i++)
        {
            sum += context[i+offset].in_degree;
        }
        norm = branch/sum;
        for (int i=0; i<branch; i++)
        {
            norm_prob.push_back( context[i+offset].in_degree*norm );
        }

        // block divison
        vector <long> small_block, large_block;
        long num_small_block = 0, num_large_block = 0;
        for (int i=0; i<branch; i++)
        {
            if ( norm_prob[i]<1 )
            {
                small_block.push_back( i );
                num_small_block++;
            }
            else
            {
                large_block.push_back( i );
                num_large_block++;
            }
        }

        // assign alias table
        long small_i, large_i;
        while (small_block.size() && large_block.size())
        {
            small_i = small_block.back();
            small_block.pop_back();
            large_i = large_block.back();
            large_block.pop_back();

            context_AT[small_i+offset].alias = context[large_i+offset].vid;
            context_AT[small_i+offset].prob = norm_prob[small_i];
            norm_prob[large_i] = norm_prob[large_i] + norm_prob[small_i] - 1;
            if (norm_prob[large_i] < 1)
            {
                small_block.push_back( large_i );
            }
            else
            {
                large_block.push_back( large_i );
            }
        }
        
        while (large_block.size())
        {
            large_i = large_block.back();
            large_block.pop_back();
            context_AT[large_i+offset].prob = 1.0;
        }

        while (small_block.size())
        {
            small_i = small_block.back();
            small_block.pop_back();
            context_AT[small_i+offset].prob = 1.0;
        }

    }
}

long proNet::NegativeSample() {
    
    long rand_v = random_gen(0, MAX_vid);
    double rand_p = random_gen(0, 1);

    //double rand_v = random_gen(0, MAX_vid);
    //double rand_p = rand_v - (int)rand_v;
   
    if (rand_p < negative_AT[rand_v].prob)
        return rand_v;
    else
        return negative_AT[rand_v].alias;

}

long proNet::NegativeFieldSample(long fid) {
    
    long rand_v = random_gen(0, MAX_vid);
    double rand_p = random_gen(0, 1);
    
    if (rand_p < negative_AT[rand_v].prob)
        return field[rand_v].vids[fid];
    else
        return field[negative_AT[rand_v].alias].vids[fid];

}

long proNet::SourceSample() {
    
    long rand_v = random_gen(0, MAX_vid);
    double rand_p = random_gen(0, 1);
    
    if (rand_p < vertex_AT[rand_v].prob)
        return rand_v;
    else
        return vertex_AT[rand_v].alias;

}

long proNet::TargetSample() {
    
    long rand_v = random_gen(0, MAX_line);
    double rand_p = random_gen(0, 1);

    if (rand_p < context_AT[rand_v].prob)
        return context[rand_v].vid;
    else
        return context_AT[rand_v].alias;

}

long proNet::TargetSample(long vid) {
    
    if (vertex[vid].branch==0) return -1;

    long rand_v = random_gen(0, vertex[vid].branch) + vertex[vid].offset;
    double rand_p = random_gen(0, 1);

    if (rand_p < context_AT[rand_v].prob)
        return context[rand_v].vid;
    else
        return context_AT[rand_v].alias;

}

vector< long > proNet::RandomWalk(long start, int steps) {

    long next = start;
    vector< long > walk;

    walk.push_back(next);
    for (int s=0; s<steps; ++s)
    {   
        if (vertex[next].branch == 0)
        {
            if (next==start)
                return walk;
            else
                next = start;
        }
        next = TargetSample(next);
        walk.push_back(next);
    }

    return walk;
}

vector< vector< long > > proNet::SkipGrams(vector< long > &walk, int window_size, int negative_samples){

    vector< vector< long > > wraper;
    vector< long > vertices;
    vector< long > contexts;
    vector< long > labels;

    int length = walk.size();
    int left;
    int right;
    vector< long > couple;
    for (int i=0; i<length; ++i)
    {
        left = i-(window_size-int(random_gen(0, window_size)));
        if (left < 0) left = 0;
        right = i+(window_size-int(random_gen(0, window_size)));
        if (right >= length) right = length-1;

        for (int j=left; j<=right; j++)
        {
            if (i==j) continue;
            vertices.push_back(walk[i]);
            contexts.push_back(walk[j]);
            labels.push_back(1);
            
            for (int n=0; n<negative_samples; ++n){
                vertices.push_back(walk[i]);
                contexts.push_back(NegativeSample());
                labels.push_back(0);
            }
        }
    }
    wraper.push_back(vertices);
    wraper.push_back(contexts);
    wraper.push_back(labels);
    
    return wraper; // [ [vertices], [contexts], [labels] ]

}


vector< vector< long > > proNet::ScaleSkipGrams(vector< long > &walk, int window_min, int window_max, int negative_samples){

    vector< vector< long > > wraper;
    vector< long > vertices;
    vector< long > contexts;
    vector< long > labels;

    int length = walk.size();
    int left;
    int right;
    vector< long > couple;
    for (int i=0; i<length; ++i)
    {
        left = i-window_max;
        if (left < 0) left = 0;
        right = i-window_min;
        if (right < 0) right = 0;

        for (int j=left; j<=right; j++)
        {
            if (i==j) continue;
            vertices.push_back(walk[i]);
            contexts.push_back(walk[j]);
            labels.push_back(1);
            
            for (int n=0; n<negative_samples; ++n){
                vertices.push_back(walk[i]);
                contexts.push_back( int(random_gen(0, MAX_vid)) );
                labels.push_back(0);
            }
        }

        left = i+window_min;
        if (left >= length) left = length-1;
        right = i+window_max;
        if (right >= length) right = length-1;

        for (int j=left; j<=right; j++)
        {
            if (i==j) continue;
            vertices.push_back(walk[i]);
            contexts.push_back(walk[j]);
            labels.push_back(1);
            
            for (int n=0; n<negative_samples; ++n){
                vertices.push_back(walk[i]);
                contexts.push_back( int(random_gen(0, MAX_vid)) );
            }
        }

    }
    wraper.push_back(vertices);
    wraper.push_back(contexts);
    wraper.push_back(labels);
    
    return wraper; // [ [vertices], [contexts], [labels] ]

}

// Optimizer
void proNet::UpdatePair(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int negative_samples, double alpha){
    
    vector< double >* w_vertex_ptr;
    vector< double >* w_context_ptr;
    vector< double > back_err;
    back_err.resize(dimension, 0.0);

    int d;
    long rand_v;
    double label, g, f, rand_p;
    
    label = 1.0;
    w_vertex_ptr = &w_vertex[vertex];
    w_context_ptr = &w_context[context];

    // 0 for postive sample, others for negative sample
    for (int neg=0; neg<=negative_samples; ++neg)
    {
        // negative sampling
        if (neg!=0){
            label = 0.0;
            w_context_ptr = &w_context[ NegativeSample() ]; // Negative Sample
        }

        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
        //f = 1.0/(1.0+exp(-f)); // sigmoid(prediction)
        //f = f/(1.0 + fabs(f)); // fast sigmoid(prediction)
        //f = tanh(f); // fast sigmoid(prediction)
        f = fastSigmoid(f); // fast sigmoid(prediction)
        //f = min(1.0, max(-1.0, f)); // relu(prediction)
        g = (label - f) * alpha; // gradient
        for (d=0; d<dimension; ++d) // store the back propagation error
            back_err[d] += g * (*w_context_ptr)[d];
        for (d=0; d<dimension; ++d) // update context
            (*w_context_ptr)[d] += g * (*w_vertex_ptr)[d];
    }
    for (d=0; d<dimension; ++d)
        (*w_vertex_ptr)[d] += back_err[d];

}

void proNet::UpdateDirectedPair(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int negative_samples, double alpha){
    
    vector< double >* w_vertex_ptr;
    vector< double >* w_context_ptr;
    vector< double > back_err;
    back_err.resize(dimension, 0.0);

    int d;
    long rand_v;
    double label, g, f, rand_p, reg;
    
    label = 1.0;
    w_vertex_ptr = &w_vertex[vertex];
    w_context_ptr = &w_context[context];

    // 0 for postive sample, others for negative sample
    for (int neg=0; neg<=negative_samples; ++neg)
    {
        // negative sampling
        if (neg!=0){
            label = 0.0;
            w_context_ptr = &w_context[ TargetSample() ]; // Negative Target Sample
        }

        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
        //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
        //f = tanh(f); // fast sigmoid(prediction)
        f = fastSigmoid(f); // fast sigmoid(prediction)
        g = (label - f) * alpha; // gradient
        for (d=0; d<dimension; ++d) // store the back propagation error
        {
            reg = 1.0*alpha*(*w_vertex_ptr)[d];
            back_err[d] += (g-reg) * (*w_context_ptr)[d];
        }
        for (d=0; d<dimension; ++d) // update context
        {
            reg = 1.0*alpha*(*w_vertex_ptr)[d];
            (*w_context_ptr)[d] += (g+reg) * (*w_vertex_ptr)[d];
        }
    }
    for (d=0; d<dimension; ++d)
        (*w_vertex_ptr)[d] += back_err[d];
    

    // opposite opt
    label = 1.0;
    w_vertex_ptr = &w_vertex[context];
    w_context_ptr = &w_context[vertex];

    for (d=0; d<dimension; ++d)
        back_err[d] = 0.0;
    // 0 for postive sample, others for negative sample
    for (int neg=0; neg<=negative_samples; ++neg)
    {
        // negative sampling
        if (neg!=0){
            label = 0.0;
            w_context_ptr = &w_context[ SourceSample() ]; // Negative Source Sample
        }

        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
        //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
        //f = tanh(f); // fast sigmoid(prediction)
        f = fastSigmoid(f); // fast sigmoid(prediction)
        g = (label - f) * alpha; // gradient
        for (d=0; d<dimension; ++d) // store the back propagation error
            back_err[d] += g * (*w_context_ptr)[d];
        for (d=0; d<dimension; ++d) // update context
            (*w_context_ptr)[d] += g * (*w_vertex_ptr)[d];
    }
    for (d=0; d<dimension; ++d)
        (*w_vertex_ptr)[d] += back_err[d];
}

void proNet::UpdatePairs(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, vector<long>& vertex, vector<long>& context, int dimension, int negative_samples, double alpha){

    vector<long>::iterator it_v = vertex.begin();
    vector<long>::iterator it_c = context.begin();
    
    while( it_v != vertex.end() )
    {
        UpdatePair(w_vertex, w_context, (*it_v), (*it_c), dimension, negative_samples, alpha);
        ++it_v;
        ++it_c;
    }
    
}


void proNet::UpdateCommunity(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int walk_steps, int negative_samples, double alpha){

    vector<double>* w_vertex_ptr;
    vector<double>* w_context_ptr;
    vector<double> back_err;
    back_err.resize(dimension, 0.0);

    int d;
    long rand_v;
    double label, g, f, rand_p;
    
    w_vertex_ptr = &w_vertex[vertex];
    w_context_ptr = &w_context[context];

    // 0 for postive sample, others for negative sample
    for (int s = -1; s < walk_steps; s++)
    {
        label = 1.0;
        if (s != -1)
        {
            context = TargetSample(context);
            if (context==-1) break;
            w_context_ptr = &w_context[ context ];
        }

        for (d=0; d<dimension; ++d)
            back_err[d] = 0.0;
        for (int neg=0; neg<=negative_samples; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0.0;
                w_context_ptr = &w_context[ NegativeSample() ];
            }

            f = 0;
            for (d=0; d<dimension; ++d) // prediciton
                f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
            //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
            //f = tanh(f); // fast sigmoid(prediction)
            f = fastSigmoid(f); // fast sigmoid(prediction)
            g = (label - f) * alpha; // gradient
            for (d=0; d<dimension; ++d) // store the back propagation error
                back_err[d] += g * (*w_context_ptr)[d];
            for (d=0; d<dimension; ++d) // update context
                (*w_context_ptr)[d] += g * (*w_vertex_ptr)[d];
        }
        for (d=0; d<dimension; ++d)
            (*w_vertex_ptr)[d] += back_err[d];

    }

}


void proNet::UpdateFieldCommunity(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int walk_steps, int negative_samples, double alpha){

    vector<double>* w_vertex_ptr;
    vector<double>* w_context_ptr;
    vector<double> back_err;
    back_err.resize(dimension, 0.0);

    int d, vid, fid;
    long rand_v;
    double label, g, f, rand_p;
    
    // vertex
    fid = field[context].fields[0];
    vid = field[vertex].vids[fid];
    w_vertex_ptr = &w_vertex[vid];

    // context
    fid = field[vertex].fields[0];
    vid = field[context].vids[fid];
    w_context_ptr = &w_context[vid];
        
    // 0 for postive sample, others for negative sample
    for (int s = 0; s <= walk_steps; s++) {
        label = 1.0;
        if (s != 0)
        {
            context = TargetSample(context);
            if (context==-1) break; 
            vid = field[context].vids[fid];
            w_context_ptr = &w_context[vid];
            fid = field[context].fields[0];
            vid = field[vertex].vids[fid];
            w_vertex_ptr = &w_vertex[vid];
            fid = field[vertex].fields[0];
        }

        for (d=0; d<dimension; ++d)
            back_err[d] = 0.0;
        for (int neg=0; neg<=negative_samples; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0.0;
                w_context_ptr = &w_context[ NegativeFieldSample(fid) ];
            }

            f = 0;
            for (d=0; d<dimension; ++d) // prediciton
                f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
            //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
            //f = tanh(f); // fast sigmoid(prediction)
            f = fastSigmoid(f); // fast sigmoid(prediction)
            g = (label - f) * alpha; // gradient
            for (d=0; d<dimension; ++d) // store the back propagation error
                back_err[d] += g * (*w_context_ptr)[d];
            for (d=0; d<dimension; ++d) // update context
                (*w_context_ptr)[d] += g * (*w_vertex_ptr)[d];
        }
        for (d=0; d<dimension; ++d)
            (*w_vertex_ptr)[d] += back_err[d];

        //if (random_gen(0, 1) < 0.2)
        //    break;
    }

}


void proNet::UpdateFieldsCommunity(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int walk_steps, int negative_samples, double alpha){

    vector<double>* w_vertex_ptr;
    vector<double>* w_context_ptr;
    vector<double> back_err;
    back_err.resize(dimension, 0.0);

    int d, vid, v_fid, c_fid;
    long rand_v;
    double label, g, f, rand_p;
    
    // 0 for postive sample, others for negative sample
    for (int s = 0; s <= walk_steps; s++) {
        label = 1.0;
        if (s != 0)
        {
            context = TargetSample(context);
            if (context==-1) break; 
        }

        for (auto c_fid: field[context].fields)
        {
            // vertex
            vid = field[vertex].vids[c_fid];
            w_vertex_ptr = &w_vertex[vid];

            for (auto v_fid: field[vertex].fields)
            {
                // context
                vid = field[context].vids[v_fid];
                w_context_ptr = &w_context[vid];

                // optimization
                for (d=0; d<dimension; ++d)
                    back_err[d] = 0.0;
                for (int neg=0; neg<=negative_samples; ++neg)
                {
                    // negative sampling
                    if (neg!=0){
                        label = 0.0;
                        w_context_ptr = &w_context[ NegativeFieldSample(v_fid) ];
                    }

                    f = 0;
                    for (d=0; d<dimension; ++d) // prediciton
                        f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
                    //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
                    //f = tanh(f); // fast sigmoid(prediction)
                    f = fastSigmoid(f); // fast sigmoid(prediction)
                    g = (label - f) * alpha; // gradient
                    for (d=0; d<dimension; ++d) // store the back propagation error
                        back_err[d] += g * (*w_context_ptr)[d];
                    for (d=0; d<dimension; ++d) // update context
                        (*w_context_ptr)[d] += g * (*w_vertex_ptr)[d];
                }
                for (d=0; d<dimension; ++d)
                    (*w_vertex_ptr)[d] += back_err[d];

            }
        }
    }
}
