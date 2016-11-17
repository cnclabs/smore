#include "LINE.h"
#include <omp.h>

LINE::LINE() {
}
LINE::~LINE() {
}

void LINE::LoadEdgeList(string filename, bool undirect) {
    rgraph.LoadEdgeList(filename, undirect);
}

void LINE::SaveWeights(string model_name){
    
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

void LINE::Init(int dim) {
   
    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;

    w_vertex.resize(rgraph.MAX_vid);
    w_context.resize(rgraph.MAX_vid);

    for (long vid=0; vid<rgraph.MAX_vid; ++vid)
    {
        w_vertex[vid].resize(dim);
        for (int d=0; d<dim;++d)
            w_vertex[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }

    for (long vid=0; vid<rgraph.MAX_vid; ++vid)
    {
        w_context[vid].resize(dim);
        for (int d=0; d<dim;++d)
            w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }
}

void LINE::UpdateSecondOrder(long vertex, long context, int negative_samples, double alpha){
    
    vector< double >* w_vertex_ptr;
    vector< double >* w_context_ptr;
    vector< double > back_err;
    back_err.resize(dim, 0.0);

    int d;
    long rand_v;
    double label, g, f, rand_p;
    
    label = 1;
    w_vertex_ptr = &w_vertex[vertex];
    w_context_ptr = &w_context[context];

    // 0 for postive sample, others for negative sample
    for (int neg=0; neg<=negative_samples; ++neg)
    {
        // negative sampling
        if (neg!=0){
            label = 0;
            w_context_ptr = &w_context[ rgraph.NegativeSample() ]; // Negative Sample
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

}

void LINE::Train(int sample_times, int negative_samples, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[LINE]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;

    sample_times *= 1000000;
    double alpha_min = alpha * 0.0001;
    double alpha_last;
    
    int current_sample = 0;
    int jobs = sample_times/workers;

    //for (int samples=0; samples<sample_times; ++samples)
    #pragma omp parallel for
    for (int worker=0; worker<workers; ++worker)
    {
        
        int count = 0;
        double _alpha = alpha;
        
        while (count<jobs)
        {
            count++;
            if (count % MONITOR == 0)
            {
                current_sample += MONITOR;
                _alpha = alpha* ( 1.0 - (double)(count)/jobs );
                alpha_last = _alpha;
                if (_alpha < alpha_min) _alpha = alpha_min;
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(current_sample)/sample_times * 100, 13);
                fflush(stdout);
            }
            
            long v1 = rgraph.SourceSample();
            long v2 = rgraph.TargetSample(v1);
            UpdateSecondOrder(v1, v2, negative_samples, _alpha);
        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", alpha_last);

}

/*
int main(int argc, char **argv){
    
    LINE *line;
    line = new LINE();
    //line->LoadEdgeList("../ml-1m/result/ml-1m.graph.train", 0);
    line->LoadEdgeList("../../../text8.pair", 0);
    //line->LoadEdgeList("in", 0);
    line->Init(200);
    
    if (argc == 2)
        line->Train(1000, 5, 0.025, atoi(argv[1]));
    else
        line->Train(1000, 5, 0.025, 4);

    line->SaveWeights("GraphRecLINE.model");

    return 0;
}
*/
