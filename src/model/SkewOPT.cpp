#include "SkewOPT.h"
#include <omp.h>

SPR::SPR() {
    char negative_method[15] = "no_degrees";
    pnet.SetNegativeMethod(negative_method);
}
SPR::~SPR() {
}

void SPR::LoadEdgeList(string filename, bool undirect) {
    pnet.LoadEdgeList(filename, undirect);
}

void SPR::SaveWeights(string model_name){

    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        model << pnet.MAX_vid << " " << dim << endl;
        for (long vid=0; vid!=pnet.MAX_vid; vid++)
        {
            model << pnet.vertex_hash.keys[vid];
            for (int d=0; d<dim; ++d)
                model << " " << w_vertex[vid][d];
            model << endl;
        }
        cout << "\tSave to <" << model_name << ">" << endl;
    }
    else
    {
        cout << "\tfail to open file" << endl;
    }
}

void SPR::Init(int dim) {

    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;

    w_vertex.resize(pnet.MAX_vid);

    for (long vid=0; vid<pnet.MAX_vid; ++vid)
    {
        w_vertex[vid].resize(dim);
        for (int d=0; d<dim;++d)
            w_vertex[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim + 0.01;
    }

}


void SPR::Train(int sample_times, int negative_samples, double alpha, double reg, double xi, double omega, int eta, int workers){

    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[Skew-OPT]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tregularization:\t\t" << reg << endl;
    cout << "\txi:\t\t\t" << xi << endl;
    cout << "\tomega:\t\t\t" << omega << endl;
    cout << "\teta:\t\t\t" << eta << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;

    unsigned long long total_sample_times = (unsigned long long)sample_times*1000000;
    double alpha_min = alpha * 0.0001;
    double alpha_last;

    unsigned long long current_sample = 0;
    unsigned long long jobs = total_sample_times/workers;

    #pragma omp parallel for
    for (int worker=0; worker<workers; ++worker)
    {

        long v1, v2, v3;
        unsigned long long count = 0;
        double _alpha = alpha;

        while (count<jobs)
        {
            v1 = pnet.SourceSample();
            v2 = pnet.TargetSample(v1);

            pnet.UpdateSBPRPair(w_vertex, w_vertex, v1, v2, dim, reg, xi, omega, eta, _alpha);

            count ++;
            if (count % MONITOR == 0)
            {
                _alpha = alpha* ( 1.0 - (double)(current_sample)/total_sample_times );
                current_sample += MONITOR;
                if (_alpha < alpha_min) _alpha = alpha_min;
                alpha_last = _alpha;
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(current_sample)/total_sample_times * 100, 13);
                fflush(stdout);
            }
        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", alpha_last);

}

