#include "HPE.h"

HPE::HPE() {}
HPE::~HPE() {}

void HPE::SaveWeights(string model_name){
    
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

void HPE::Init(int dim) {
   
    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;

    w_vertex.resize(pnet.MAX_vid);
    w_context.resize(pnet.MAX_vid);

    for (long vid=0; vid<pnet.MAX_vid; ++vid)
    {
        w_vertex[vid].resize(dim);
        for (int d=0; d<dim;++d)
        {
            w_vertex[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
        }
    }

    for (long vid=0; vid<pnet.MAX_vid; ++vid)
    {
        w_context[vid].resize(dim);
        for (int d=0; d<dim;++d)
            w_context[vid][d] = 0.0;
            //w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }
}


void HPE::Train(int sample_times, int walk_steps, int negative_samples, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[HPE]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\twalk_steps:\t\t" << walk_steps << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;

    unsigned int total_sample_times = sample_times*1000000;
    double alpha_min = alpha * 0.0001;
    double alpha_last = alpha;
    
    unsigned int current_sample = 0;
    unsigned int jobs = total_sample_times/workers;

    #pragma omp parallel for
    for (int worker=0; worker<workers; ++worker)
    {
        unsigned int count = 0;
        double _alpha = alpha;
        long v1, v2;
        
        while (count<jobs)
        {
            count ++;
            if (count % MONITOR == 0)
            {
                current_sample += MONITOR;
                _alpha = alpha* ( 1.0 - (double)(count)/jobs );
                if (_alpha < alpha_min)
                    _alpha = alpha_min;
                if (_alpha < alpha_last)
                    alpha_last = _alpha;
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", alpha_last, (double)(current_sample)/total_sample_times * 100, 13);
                fflush(stdout);
            }
            
            v1 = pnet.SourceSample();
            v2 = pnet.TargetSample(v1);
            pnet.UpdateCommunity(w_vertex, w_context, v1, v2, dim, walk_steps, negative_samples, _alpha);
            pnet.UpdatePair(w_vertex, w_context, v2, v1, dim, negative_samples, _alpha);
        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", alpha_last);

}

