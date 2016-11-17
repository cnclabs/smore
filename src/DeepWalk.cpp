#include "DeepWalk.h"

DeepWalk::DeepWalk () {
}

DeepWalk::~DeepWalk () {
}

void DeepWalk::LoadEdgeList(string filename, bool undirect) {
    rgraph.LoadEdgeList(filename, undirect);
}

void DeepWalk::SaveWeights(string model_name){
    
    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        for (auto k: rgraph.keys)
        {
            model << k;
            for (int d=0; d<dim; ++d)
                model << " " << w_vertex[rgraph.kmap[k]][d];
            model << endl;
        }
        cout << "\tSave to <" << model_name << ">" << endl;
    }
    else
    {
        cout << "\tfail to open file" << endl;
    }
}

void DeepWalk::Init(int dim) {
    
    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;
    
    w_vertex.resize(rgraph.MAX_vid);
    w_context.resize(rgraph.MAX_vid);
    //w_vertex = new double* [rgraph.MAX_vid];
    //w_context = new double* [rgraph.MAX_vid];

    for (long vid=0; vid<rgraph.MAX_vid; ++vid)
    {
        w_vertex[vid].resize(dim);
        //w_vertex[vid] = new double [dim];
        for (int d=0; d<dim;++d)
            w_vertex[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }

    for (long vid=0; vid<rgraph.MAX_vid; ++vid)
    {
        w_context[vid].resize(dim);
        //w_context[vid] = new double [dim];
        for (int d=0; d<dim;++d)
            w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }
}

void DeepWalk::Update(vector<long>& vertex, vector<long>& context, int negative_samples, double alpha){
    
    vector<long>::iterator it_v = vertex.begin();
    vector<long>::iterator it_c = context.begin();

    vector<double>* w_vertex_ptr;
    vector<double>* w_context_ptr;
    double* back_err = new double[dim];

    int d, label;
    double g, f;
    while( it_v != vertex.end() )
    {
        label = 1;
        w_vertex_ptr = &w_vertex[(*it_v)];
        w_context_ptr = &w_context[(*it_c)];
        for (d=0; d<dim; ++d)
            back_err[d] = 0.0;
        
        // 0 for postive sample, others for negative sample
        for (int neg=0; neg<=negative_samples; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0;
                w_context_ptr = &w_context[ rgraph.NegativeSample() ];
            }

            f = 0;
            for (d=0; d<dim; ++d) // prediciton
                f += (*w_vertex_ptr)[d] * (*w_context_ptr)[d];
            f = f/(1.0 + fabs(f)); // sigmoid(prediction)
            g = (label - f) * alpha; // gradient
            for (d=0; d<dim; ++d) // store the back propagation error
                back_err[d] += g * (*w_context_ptr)[d];
            for (d=0; d<dim; ++d) // update context
                (*w_context_ptr)[d] += g * (*w_vertex_ptr)[d];
        }
        for (d=0; d<dim; ++d)
            (*w_vertex_ptr)[d] += back_err[d];

        ++it_v;
        ++it_c;
    }
    
    delete [] back_err;

}

void DeepWalk::Train(int walk_times, int walk_steps, int window_size, int negative_samples, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\tDeepWalk" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\twalk_times:\t\t" << walk_times << endl;
    cout << "\twalk_steps:\t\t" << walk_steps << endl;
    cout << "\twindow_size:\t\t" << window_size << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;


    long total = walk_times*rgraph.MAX_vid;
    double alpha_min = alpha*0.0001;
    double _alpha = alpha;
    long count = 0;

    for (int t=0; t<walk_times; ++t)
    {
        // shuffle the order for random keys access        
        std::vector<long> random_keys(rgraph.MAX_vid);
        for (long vid = 0; vid < rgraph.MAX_vid; vid++) {
            random_keys[vid] = vid;
        }
        for (long vid = 0; vid < rgraph.MAX_vid; vid++) {
            int rdx = vid + rand() % (rgraph.MAX_vid - vid); // careful here!
            swap(random_keys[vid], random_keys[rdx]);
        }

        #pragma omp parallel for
        for (long vid=0; vid<rgraph.MAX_vid; ++vid)
        {
            if (_alpha < alpha_min) _alpha = alpha_min;

            vector<long> walks = rgraph.RandomWalk(random_keys[vid], walk_steps);
            vector<vector<long>> train_data = rgraph.SkipGrams(walks, window_size, 0);
            Update(train_data[0], train_data[1], negative_samples, _alpha);
            
            count++;
            if (count % MONITOR == 0)
            {
                _alpha = alpha* ( 1.0 - (double)(count)/total );
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(count)/total * 100, 13);
                fflush(stdout);
            }

        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", _alpha);

}

/*
int main(int argc, char **argv){
    
    DeepWalk *dw;
    dw = new DeepWalk();
    dw->LoadEdgeList("../../../text8.pair", 0);
    dw->Init(64);
    
    if (argc == 2)
        dw->Train(100, 40, 1, 5, 0.025, atoi(argv[1]));
    else
        dw->Train(100, 40, 1, 5, 0.025, 4);

    dw->SaveWeights("GraphRecDW.model");

    return 0;
}
*/
