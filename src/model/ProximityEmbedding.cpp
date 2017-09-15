#include "ProximityEmbedding.h"

PE::PE() {}
PE::~PE() {}

void PE::LoadWalkMeta(string filename) {
    pnet.LoadWalkMeta(filename);
}

void PE::SaveWeights(string model_name){
    
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

void PE::Init(int dim) {
   
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
            w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }
}


void PE::Train(int sample_times, int walk_steps, int negative_samples, double bfs, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[PE]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\twalk_steps:\t\t" << "given" << endl;
//    cout << "\tbfs:\t\t\t" << bfs << endl;
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
            pnet.UpdateDCommunity(w_vertex, w_context, v1, v2, dim, negative_samples, bfs, _alpha);
            //pnet.UpdatePair(w_vertex, w_context, v2, v1, dim, negative_samples, _alpha);
        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", alpha_last);

}

