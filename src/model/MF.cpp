#include "MF.h"
#include <omp.h>

MF::MF() {
}
MF::~MF() {
}

void MF::LoadFieldMeta(string filename) {
    pnet.LoadFieldMeta(filename);
}

void MF::LoadEdgeList(string filename, bool undirect) {
    pnet.LoadEdgeList(filename, undirect);
}

void MF::SaveWeights(string model_name){
    
    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        model << pnet.MAX_vid << " " << dim << endl;
        for (auto k: pnet.keys)
        {
            model << k;
            for (int d=0; d<dim; ++d)
                model << " " << w_vertex[pnet.kmap[k]][d];
            model << endl;
        }
        cout << "\tSave to <" << model_name << ">" << endl;
    }
    else
    {
        cout << "\tfail to open file" << endl;
    }
}

void MF::Init(int dim) {
   
    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;

    w_vertex.resize(pnet.MAX_vid);
    w_context.resize(pnet.MAX_vid);

    //random_device rd;
    //mt19937 e2(rd());
    //normal_distribution<double> dist(0.0, 0.01);

    for (long vid=0; vid<pnet.MAX_vid; ++vid)
    {
        w_vertex[vid].resize(dim);
        w_context[vid].resize(dim);
        for (int d=0; d<dim;++d)
            //w_vertex[vid][d] = ran_gaussian(0.0, 0.1);
            w_vertex[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
        for (int d=0; d<dim;++d)
            //w_context[vid][d] = ran_gaussian(0.0, 0.1);
            w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }

}


void MF::Train(int sample_times, int negative_samples, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[MF]" << endl;

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
                if (_alpha < alpha_min) _alpha = alpha_min;
                alpha_last = _alpha;
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(current_sample)/sample_times * 100, 13);
                fflush(stdout);
            }
            
            long v1 = pnet.SourceSample();
            long v2 = pnet.TargetSample(v1);
            //pnet.UpdatePair(w_vertex, w_vertex, v1, v2, dim, negative_samples, _alpha);
            pnet.UpdateDirectedPair(w_vertex, w_vertex, v1, v2, dim, negative_samples, _alpha);
            //pnet.UpdateCommunity(w_vertex, w_context, v1, v2, dim, 3, negative_samples, _alpha);
        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", alpha_last);

}

