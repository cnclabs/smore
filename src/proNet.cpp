#include "proNet.h"

proNet::proNet() {

    MAX_line=0;
    MAX_vid=0;
    MAX_fvid=0;
    MAX_field=0;

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
    for (int i = 0; i != SIGMOID_TABLE_SIZE + 1; i++) {
        double x = i * 2.0 * MAX_SIGMOID / SIGMOID_TABLE_SIZE - MAX_SIGMOID;
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

void proNet::InitNegTable() {

    double sum = 0, cur_sum = 0, por = 0;
    long vid = 0;
    neg_table.resize(MAX_NEG);
    for (long k = 0; k != MAX_vid; k++) sum += pow(vertex[k].in_degree+vertex[k].out_degree, POWER_SAMPLE);
    for (long k = 0; k != MAX_NEG; k++)
    {
        if ((double)(k + 1) / MAX_NEG > por)
        {
            cur_sum += pow(vertex[vid].in_degree+vertex[vid].out_degree, POWER_SAMPLE);
            por = cur_sum / sum;
            vid++;
        }
        neg_table[k] = vid - 1;
    }
}

int proNet::InsertHashTable(char *key) {

    unsigned int pos = BKDRHash(key);
    while (hash_table[pos] != -1)
        pos = (pos + 1) % HASH_TABLE_SIZE;
    hash_table[pos] = MAX_vid;
    kmap[ strdup(key) ] = MAX_vid;
    keys.push_back( strdup(key) );
    MAX_vid++;

    return MAX_vid-1;
}

int proNet::SearchHashTable(char *key) {

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
    long vid1, vid2;
    double w;
    unordered_map< long, vector< long > > graph;
    unordered_map< long, vector< double > > edge;

    cout << "Connections Loading:" << endl;
    unsigned long long line = 0;
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
            
            graph[vid1].push_back(vid2);
            edge[vid1].push_back(w);

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
        fclose(fin);
    }
    cout << "\tProgress:\t\t100.00 %\r" << endl;
    cout << "\t# of vertex:\t\t" << MAX_vid << endl;

    cout << "Build the Alias Method:" << endl;
    if (undirect)
        MAX_line *= 2;
    BuildAliasMethod(graph, edge);
    InitNegTable();
    cout << "\tFinished." << endl;

}

void proNet::LoadWalkMeta(string filename) {

    FILE *fin;
    char c_line[1000];
    unsigned long long max_line=0;

    cout << "Walk Data Loading:" << endl;
    fin = fopen(filename.c_str(), "rb");
    while (fgets(c_line, sizeof(c_line), fin))
    {
        if (max_line % MONITOR == 0)
        {
            printf("\t# of walking data:\t%llu%c", max_line, 13);
        }
        ++max_line;
    }
    fclose(fin);
    cout << "\t# of walking data:\t" << max_line << endl;

    char v[160];
    int w;
    long vid;
    dynamic_walk.resize(MAX_vid, 3);

    fin = fopen(filename.c_str(), "rb");
    for (unsigned long long line = 0; line != max_line; line++)
    {
        if ( fscanf(fin, "%s %d", v, &w)!=2 )
        {
            cout << "line " << line << " contains wrong number of data" << endl; 
            continue;
        }

        vid = SearchHashTable(v);
        if (vid != -1)
            dynamic_walk[vid] = w;
        else
            cout << "vertex " << v << " is not in given network" << endl;
    }

}

void proNet::LoadFieldMeta(string filename) {

    // calculate the # of meta data
    FILE *fin;
    char c_line[1000];
    unsigned long long max_line=0;
    
    field.resize(MAX_vid);
    MAX_fvid = MAX_vid;

    cout << "Meta Data Preview:" << endl;
    fin = fopen(filename.c_str(), "rb");
    while (fgets(c_line, sizeof(c_line), fin))
    {
        if (max_line % MONITOR == 0)
        {
            printf("\t# of meta data:\t\t%llu%c", max_line, 13);
        }
        ++max_line;
    }
    fclose(fin);
    cout << "\t# of meta data:\t\t" << max_line << endl;
    
    char v[160], meta[160];
    long vid;
    map< char*, long, cmp_char > meta_idx;

    cout << "Meta Data Loading:" << endl;
    fin = fopen(filename.c_str(), "rb");
    for (unsigned long long line = 0; line != max_line; line++)
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
        vid = SearchHashTable(v);
        if (vid != -1)
            field[ vid ].fields[0] = meta_idx[ strdup(meta) ];
        else
            cout << "vertex " << v << " is not in given network" << endl;
       
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
    for (long vid=0; vid<MAX_vid; vid++)
    {
        field[ vid ].vids.resize(MAX_field);
        for (long i=0; i<MAX_field; i++)
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
    
    long vid;
    unsigned long long line_g=0, line_e=0;
    long offset = 0;
    for (long v1=0; v1!=MAX_vid; v1++)
    {
        vertex[v1].offset = offset;
        vertex[v1].branch = graph[v1].size();
        offset += graph[v1].size();

        //for (auto v2: graph[v1])
        for (int i=0; i<graph[v1].size(); i++)
        {
            context[line_g].vid = graph[v1][i];
            line_g++;
        }
        //for (auto w: edge[v1])
        for (int i=0; i<edge[v1].size(); i++)
        {
            vertex[v1].out_degree += edge[v1][i];
            context[line_e].in_degree = edge[v1][i];
            line_e++;
        }
    }

    for (unsigned long long line=0; line!=MAX_line; line++)
    {
        vid = context[line].vid;
        vertex[vid].in_degree += context[line].in_degree;
    }
    
    graph.clear();
    edge.clear();

    // compute alias table
    cout << "\tAlias Table Constructing ..." << endl;
    vector<double> distribution;
    
    // Alias table for source vertices
    distribution.resize(MAX_vid);
    for (long v=0; v<MAX_vid; v++)
    {
        distribution[v] = vertex[v].out_degree;
    }
    vertex_AT = AliasMethod(distribution, 1.0);
    
    // Alias table for negative sampling
    distribution.resize(MAX_vid);
    for (long v=0; v<MAX_vid; v++)
    {
        distribution[v] = vertex[v].in_degree + vertex[v].out_degree;
    }
    negative_AT = AliasMethod(distribution, POWER_SAMPLE);

    // Alias table for context vertices
    long branch;
    for (long vid=0; vid<MAX_vid;vid++)
    {
        offset = vertex[vid].offset;
        branch = vertex[vid].branch;
        
        distribution.resize(branch);
        for (long i=0; i<branch; i++)
        {
            distribution[i] = context[i+offset].in_degree;
        }
        vector<AliasTable> sub_at = AliasMethod(distribution, 1.0);
        for (long i=0; i<branch; i++)
        {
            sub_at[i].alias += offset;
        }
        context_AT.insert(context_AT.end(), sub_at.begin(), sub_at.end());
    }

}

vector<AliasTable> proNet::AliasMethod(vector<double>& distribution, double power) {
    
    vector<AliasTable> alias_table;

    // normalization of vertices weights
    double sum, norm;
    vector<double> norm_prob;
    alias_table.resize(distribution.size());
    
    sum = 0;
    vector<double>::iterator distribution_i;

    for (distribution_i=distribution.begin(); distribution_i!=distribution.end(); ++distribution_i)
    {
        sum += pow(*distribution_i, POWER_SAMPLE);
    }
    norm = distribution.size()/sum;
    
    for (distribution_i=distribution.begin(); distribution_i!=distribution.end(); ++distribution_i)
    {
        norm_prob.push_back( pow(*distribution_i, POWER_SAMPLE)*norm );
    }

    // block divison
    vector<long> small_block, large_block;
    
    for (long pos=0; pos!=norm_prob.size(); ++pos)
    {
        if ( norm_prob[pos]<1 )
        {
            small_block.push_back( pos );
        }
        else
        {
            large_block.push_back( pos );
        }
    }

    // assign alias table
    long small_pos, large_pos;

    while (small_block.size() && large_block.size())
    {
        small_pos = small_block.back();
        small_block.pop_back();
        large_pos = large_block.back();
        large_block.pop_back();

        alias_table[small_pos].alias = large_pos;
        alias_table[small_pos].prob = norm_prob[small_pos];
        norm_prob[large_pos] = norm_prob[large_pos] + norm_prob[small_pos] - 1;
        if (norm_prob[large_pos] < 1)
        {
            small_block.push_back( large_pos );
        }
        else
        {
            large_block.push_back( large_pos );
        }
    }

    while (large_block.size())
    {
        large_pos = large_block.back();
        large_block.pop_back();
        alias_table[large_pos].prob = 1.0;
    }

    while (small_block.size())
    {
        small_pos = small_block.back();
        small_block.pop_back();
        alias_table[small_pos].prob = 1.0;
    }
    
    return alias_table;
}


long proNet::NegativeSample() {
    
    long rand_v = random_gen(0, MAX_vid);
    double rand_p = random_gen(0, 1);
      
    if (rand_p < negative_AT[rand_v].prob)
        return rand_v;
    else
        return negative_AT[rand_v].alias;

}

long proNet::NegativeFieldSample(long fid) {
    
    double rand_p = random_gen(0, 1);
    long rand_v = random_gen(0, MAX_vid);
   
    if (rand_p < negative_AT[rand_v].prob)
        return field[rand_v].vids[fid];
    else
        return field[negative_AT[rand_v].alias].vids[fid];

}

long proNet::SourceSample() {
    
    double rand_p = random_gen(0, 1);
    long rand_v = random_gen(0, MAX_vid);
        
    if (rand_p < vertex_AT[rand_v].prob)
        return rand_v;
    else
        return vertex_AT[rand_v].alias;

}

long proNet::TargetSample() {
    
    double rand_p = random_gen(0, 1);
    long rand_v = random_gen(0, MAX_line);

    if (rand_p < context_AT[rand_v].prob)
        return context[rand_v].vid;
    else
        return context_AT[rand_v].alias;

}

long proNet::TargetSample(long vid) {
    
    if (vertex[vid].branch==0) return -1;

    double rand_p = random_gen(0, 1);
    long rand_v = random_gen(0, vertex[vid].branch) + vertex[vid].offset;
    
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
    int reduce;
    vector< long > couple;
    for (int i=0; i<length; ++i)
    {
        reduce = random_gen(0, window_size) + 1;
        //reduce = gsl_rng_uniform(gsl_r)*window_size + 1;
        left = i-reduce;
        if (left < 0) left = 0;
        right = i+reduce;
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
                //contexts.push_back( int(gsl_rng_uniform(gsl_r)* MAX_vid) );
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
                //contexts.push_back( int(gsl_rng_uniform(gsl_r)* MAX_vid) );
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

    negative_samples += 1;
    // 0 for postive sample, others for negative sample
    for (int neg=0; neg!=negative_samples; ++neg)
    {
        // negative sampling
        if (neg!=0){
            label = 0.0;
            w_context_ptr = &w_context[ NegativeSample() ]; // Negative Sample
            //w_context_ptr = &w_context[ neg_table[random_gen(0, MAX_NEG)] ]; // Negative Sample
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

void proNet::UpdateDirectedPair(vector< vector<double> >& w_vertexA, vector< vector<double> >& w_vertexB, vector< vector<double> >& w_context, long vertex, long context, int dimension, int negative_samples, double alpha){
    
    vector< double >* w_vertexA_ptr;
    vector< double >* w_vertexB_ptr;
    vector< double >* w_context_ptr;
    vector< double > back_err;

    int d;
    long rand_v;
    double label, g, f, rand_p, reg;
    
    label = 1.0;
    w_vertexA_ptr = &w_vertexA[vertex];
    w_vertexB_ptr = &w_vertexB[context];
    w_context_ptr = &w_context[context];

/*
    f = 0;
    for (d=0; d<dimension; ++d) // prediciton
        f += (*w_vertexA_ptr)[d] * (*w_vertexB_ptr)[d];
    g = (label - f); // gradient
    for (d=0; d<dimension; ++d) // store the back propagation error
    {
        back_err[d] += alpha*( g*(*w_vertexB_ptr)[d] - 0.01*(*w_vertexA_ptr)[d] );
    }
    for (d=0; d<dimension; ++d) // update context
    {
        (*w_vertexB_ptr)[d] += alpha*( g*(*w_vertexA_ptr)[d] - 0.01*(*w_vertexB_ptr)[d] );
    }
    for (d=0; d<dimension; ++d)
        (*w_vertexA_ptr)[d] += back_err[d];

    label = 0.0;
    for (int neg=0; neg<negative_samples; ++neg)
    {
        f = 0.0;
        back_err.resize(dimension, 0.0);
        long vid = (long)random_gen(0, MAX_vid);
        while(field[vid].fields[0]==field[vertex].fields[0])
        {
            vid = (long)random_gen(0, MAX_vid);
        }
        w_vertexB_ptr = &w_vertexB[ vid ]; // Negative Target Sample


        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertexA_ptr)[d] * (*w_vertexB_ptr)[d];
        g = (label - f); // gradient
        for (d=0; d<dimension; ++d) // store the back propagation error
        {
            back_err[d] += alpha*( g*(*w_vertexB_ptr)[d] - 0.01*(*w_vertexA_ptr)[d] );
        }
        for (d=0; d<dimension; ++d) // update context
        {
            (*w_vertexB_ptr)[d] += alpha*( g*(*w_vertexA_ptr)[d] - 0.01*(*w_vertexB_ptr)[d] );
        }
        for (d=0; d<dimension; ++d)
            (*w_vertexA_ptr)[d] += back_err[d];
    }
*/

    back_err.resize(dimension, 0.0);
    // 0 for postive sample, others for negative sample
    for (int neg=0; neg<=negative_samples; ++neg)
    {
        // negative sampling
        if (neg!=0){
            label = 0.0;
            long vid = (long)random_gen(0, MAX_vid);
            while(field[vid].fields[0]==field[vertex].fields[0])
            {
                vid = (long)random_gen(0, MAX_vid);
            }
            w_vertexB_ptr = &w_vertexB[ vid ]; // Negative Target Sample
        }

        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertexA_ptr)[d] * (*w_vertexB_ptr)[d];
        //f = fastSigmoid(f); // fast sigmoid(prediction)
        g = (label - f); // gradient
        for (d=0; d<dimension; ++d) // store the back propagation error
        {
            back_err[d] += alpha*( g*(*w_vertexB_ptr)[d] - 0.0*(*w_vertexA_ptr)[d] );
        }
        for (d=0; d<dimension; ++d) // update context
        {
            (*w_vertexB_ptr)[d] += alpha*( g*(*w_vertexA_ptr)[d] - 0.0*(*w_vertexB_ptr)[d] );
        }
    }
    //for (d=0; d<dimension; ++d)
    //    (*w_vertexA_ptr)[d] += back_err[d];
    //back_err.resize(dimension, 0.0);
    
    w_vertexB_ptr = &w_vertexB[context];
    //w_vertexB_ptr = &w_context[context];
    for (int s = 0; s < 5; s++)
    {
        //label = 1.0;
        //if (s!=0)
        //{
            context = TargetSample(context);
            if (context==-1) break;
            context = TargetSample(context);
            if (context==-1) break;
            w_vertexB_ptr = &w_context[ context ];
        //}
        //for (int neg=0; neg<=negative_samples; ++neg)
        for (int neg=0; neg<=0; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0.0;
                long vid = (long)random_gen(0, MAX_vid);
                while(field[vid].fields[0]==field[vertex].fields[0])
                {
                    vid = (long)random_gen(0, MAX_vid);
                }
                //long vid = NegativeSample();
                w_vertexB_ptr = &w_context[ vid ]; // Negative Target Sample
            }

            f = 0;
            for (d=0; d<dimension; ++d) // prediciton
                f += (*w_vertexA_ptr)[d] * (*w_vertexB_ptr)[d];
            g = (label - f) * alpha * 0.2; // gradient

            for (d=0; d<dimension; ++d) // store the back propagation error
            {
                back_err[d] += alpha*( g*(*w_vertexB_ptr)[d] - 0.01*(*w_vertexA_ptr)[d] );
            }
            for (d=0; d<dimension; ++d) // update context
            {
                (*w_vertexB_ptr)[d] += alpha*( g*(*w_vertexA_ptr)[d] - 0.01*(*w_vertexB_ptr)[d] );
                //(*w_vertexB_ptr)[d] += alpha*( g*(*w_vertexA_ptr)[d] );
            }
        }
    }

    /*
    label = 0.0;
    for (int s = 0; s < 3; s++)
    {
        w_vertexB_ptr = &w_vertexB[ random_gen(0, MAX_vid) ];
        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertexA_ptr)[d] * (*w_vertexB_ptr)[d];
        g = (label - f) * alpha * 0.1; // gradient

        for (d=0; d<dimension; ++d) // store the back propagation error
        {
            back_err[d] += alpha*( g*(*w_vertexB_ptr)[d] - 0.01*(*w_vertexA_ptr)[d] );
        }
        for (d=0; d<dimension; ++d) // update context
        {
            (*w_vertexB_ptr)[d] += alpha*( g*(*w_vertexA_ptr)[d] - 0.01*(*w_vertexB_ptr)[d] );
        }
    }
    */
    for (d=0; d<dimension; ++d)
        (*w_vertexA_ptr)[d] += back_err[d];
    

/*
    // 0 for postive sample, others for negative sample
    for (int s = 0; s <= 3; s++)
    {
        label = 1.0;
        //if (s != -1)
        //{
            context = TargetSample(context);
            if (context==-1) break;
            w_context_ptr = &w_context[ context ];
        //}

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
                f += (*w_vertexA_ptr)[d] * (*w_context_ptr)[d];
            //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
            //f = tanh(f); // fast sigmoid(prediction)
            //f = fastSigmoid(f); // fast sigmoid(prediction)
            g = (label - f) * alpha; // gradient
            for (d=0; d<dimension; ++d) // store the back propagation error
                back_err[d] += g * (*w_context_ptr)[d];
            for (d=0; d<dimension; ++d) // update context
                (*w_context_ptr)[d] += g * (*w_vertexA_ptr)[d];
        }
        for (d=0; d<dimension; ++d)
            (*w_vertexA_ptr)[d] += back_err[d];
    }
*/

/*
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
            long vid = (long)random_gen(0, MAX_vid);
            while(field[vid].fields[0]==field[context].fields[0])
            {
                vid = (long)random_gen(0, MAX_vid);
            }
            w_context_ptr = &w_context[ vid ];
            //w_context_ptr = &w_context[ NegativeSample() ]; // Negative Source Sample
            //w_context_ptr = &w_context[ SourceSample() ]; // Negative Source Sample
        }

        f = 0;
        for (d=0; d<dimension; ++d) // prediciton
            f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
        //f = f/(1.0 + fabs(f)); // sigmoid(prediction)
        //f = tanh(f); // fast sigmoid(prediction)
        //f = fastSigmoid(f); // fast sigmoid(prediction)
        f = min(1.0, f);
        f = max(0.0, f);
        g = (label - f); // gradient
        for (d=0; d<dimension; ++d) // store the back propagation error
        {
            back_err[d] += alpha*( g*(*w_context_ptr)[d] - 0.01*(*w_vertex_ptr)[d] );
//            reg = 1.0*alpha*(*w_vertex_ptr)[d];
//            back_err[d] += (g-reg) * (*w_context_ptr)[d];
        }
        for (d=0; d<dimension; ++d) // update context
        {
            (*w_context_ptr)[d] += alpha*( g*(*w_vertex_ptr)[d] - 0.01*(*w_context_ptr)[d] );
//            reg = 1.0*alpha*(*w_context_ptr)[d];
//            (*w_context_ptr)[d] += (g-reg) * (*w_vertex_ptr)[d];
        }
    }
    for (d=0; d<dimension; ++d)
        (*w_vertex_ptr)[d] += back_err[d];
*/

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


void proNet::UpdateDCommunity(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int negative_samples, double bfs, double alpha){

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
    for (int s = -1; s < dynamic_walk[context]; s++)
    {
        label = 1.0;
        if (s != -1)
        {
            context = TargetSample(vertex);
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

    int d, vid, v_fid, c_fid;
    long rand_v;
    double label, g, f, rand_p;
    
    v_fid = field[vertex].fields[0];
    c_fid = field[context].fields[0];

    // vertex
    vid = field[vertex].vids[c_fid];
    w_vertex_ptr = &w_vertex[vid];

    // context
    vid = field[context].vids[v_fid];
    w_context_ptr = &w_context[vid];
    
    // 0 for postive sample, others for negative sample
    for (int s = 0; s <= walk_steps; s++) {
        label = 1.0;

        if (s != 0)
        {
            context = TargetSample(context);
            if (context==-1) break;
            c_fid = field[context].fields[0];
            
            vid = field[context].vids[v_fid];
            w_context_ptr = &w_context[vid];
            
            vid = field[vertex].vids[c_fid];
            w_vertex_ptr = &w_vertex[vid];
        }

        for (d=0; d<dimension; ++d)
            back_err[d] = 0.0;
        for (int neg=0; neg<=negative_samples; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0.0;

                vid = NegativeSample();
                while(field[vid].fields[0]!=c_fid)
                {
                    vid = NegativeSample();
                }
                vid = field[vid].vids[v_fid];
                w_context_ptr = &w_context[vid];
                //w_context_ptr = &w_context[ vid ];
                //w_context_ptr = &w_context[ NegativeSample() ];
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


void proNet::UpdateMSFieldCommunity(vector< vector<double> >& w_vertex, vector< vector<double> >& w_context, long vertex, long context, int dimension, int walk_steps, int negative_samples, double alpha){

    vector<double>* w_vertex_ptr;
    vector<double>* w_context_ptr;
    vector<double> back_err;
    back_err.resize(dimension, 0.0);

    int d, vid, v_fid, c_fid;
    long rand_v;
    double label, g, f, rand_p;
    
    v_fid = field[vertex].fields[0];
    c_fid = field[context].fields[0];

    // vertex
    vid = field[vertex].vids[c_fid];
    w_vertex_ptr = &w_vertex[vid];

    // context
    w_context_ptr = &w_context[context];
    
    // 0 for postive sample, others for negative sample
    for (int s = 0; s <= walk_steps; s++) {
        label = 1.0;

        if (s != 0)
        {
            context = TargetSample(context);
            if (context==-1) break;
            w_context_ptr = &w_context[context];
            c_fid = field[context].fields[0];
            vid = field[vertex].vids[c_fid];
            w_vertex_ptr = &w_vertex[vid];
        }

        for (d=0; d<dimension; ++d)
            back_err[d] = 0.0;
        for (int neg=0; neg<=negative_samples; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0.0;
                //w_context_ptr = &w_context[ NegativeFieldSample(0) ];
                vid = NegativeSample();
                while(field[vid].fields[0]!=c_fid)
                {
                    vid = NegativeSample();
                }
                w_context_ptr = &w_context[ vid ];
            }

            f = 0;
            for (d=0; d<dimension; ++d) // prediciton
                f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
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
