#include "FINE.h"
#include <omp.h>

FINE::FINE() {}
FINE::~FINE() {}

void FINE::LoadFieldMeta(string filename) {
    rgraph.LoadFieldMeta(filename);
}

void FINE::SaveWeights(string model_name){
    
    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        for (auto k: rgraph.keys)
        {
            model << k;
            int o_vid = rgraph.kmap[k];
            for (auto vid: rgraph.field[o_vid].vids)
            {
                for (int d=0; d<dim; ++d)
                    model << " " << w_vertex[rgraph.kmap[k]][d];
            }
            model << endl;
        }
        cout << "\tSave to <" << model_name << ">" << endl;
    }
    else
    {
        cout << "\tfail to open file" << endl;
    }
}

void FINE::Init(int dim) {
   
    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;
    
    w_vertex.resize(rgraph.MAX_fvid);
    w_context.resize(rgraph.MAX_fvid);

    for (long vid=0; vid<rgraph.MAX_vid; ++vid)
    {
        for (auto fvid: rgraph.field[vid].vids)
        {
            w_vertex[fvid].resize(dim);
            for (int d=0; d<dim;++d)
                w_vertex[fvid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
        }
    }

    for (long vid=0; vid<rgraph.MAX_vid; ++vid)
    {
        for (auto fvid: rgraph.field[vid].vids)
        {
            w_context[fvid].resize(dim);
            for (int d=0; d<dim;++d)
                w_context[fvid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
        }
    }
}

void FINE::Train(int sample_times, int negative_samples, int walk_steps, double alpha, int workers) {
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[FINE]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\twalk_steps:\t\t" << walk_steps << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;

    sample_times *= 1000000;
    double alpha_min = alpha * 0.0001;
    double _alpha;
    
    int count = 0;

    #pragma omp parallel for
    for (int samples=0; samples<sample_times; ++samples)
    {
        
        count++;
        if (count % MONITOR == 0)
        {
            _alpha = alpha* ( 1.0 - (double)(count)/sample_times );
            if (_alpha < alpha_min) _alpha = alpha_min;
            printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(count)/sample_times * 100, 13);
            fflush(stdout);
        }
            
        long v1 = rgraph.SourceSample();
        long v2 = rgraph.TargetSample(v1);
        //rgraph.UpdateFieldCommunity(w_vertex, w_context, v1, v2, dim, negative_samples, walk_steps, _alpha);

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", _alpha);

}

