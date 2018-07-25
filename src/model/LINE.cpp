#include "LINE.h"
#include <omp.h>

LINE::LINE() {
}
LINE::~LINE() {
}

void LINE::LoadEdgeList(string filename, bool undirect) {
    pnet.LoadEdgeList(filename, undirect);
}

void LINE::SaveWeights(string model_name){
    
    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        model << pnet.MAX_vid << " " << dim << endl;

        if (order==1)
        {
            for (long vid=0; vid!=pnet.MAX_vid; vid++)
            {
                model << pnet.vertex_hash.keys[vid];
                for (int d=0; d<dim; ++d)
                    model << " " << w_vertex_o1[vid][d];
                model << endl;
            }
        }
        else
        {
            for (long vid=0; vid!=pnet.MAX_vid; vid++)
            {
                model << pnet.vertex_hash.keys[vid];
                for (int d=0; d<dim; ++d)
                    model << " " << w_vertex[vid][d];
                model << endl;
            }
        }
        cout << "\tSave to <" << model_name << ">" << endl;
    }
    else
    {
        cout << "\tfail to open file" << endl;
    }
}

void LINE::Init(int dimension, int order) {
   
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dimension << endl;
    this->dim = (int)(dimension);
    if (order == 1)
        this->order = order;

    if (order==1)
    {
        w_vertex_o1.resize(pnet.MAX_vid);
        for (long vid=0; vid<pnet.MAX_vid; ++vid)
        {
            w_vertex_o1[vid].resize(dim);
            for (int d=0; d<dim;++d)
            {
                w_vertex_o1[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
            }
        }
    }
    else
    {
        w_vertex.resize(pnet.MAX_vid);
        w_context.resize(pnet.MAX_vid);
        int alignment;
        //alignment = posix_memalign((void **)&w_vertex1, 128, pnet.MAX_vid * dim * sizeof(double));
        //alignment = posix_memalign((void **)&w_context1, 128, pnet.MAX_vid * dim * sizeof(double));
        w_vertex1 = (double *)malloc(pnet.MAX_vid * dim * sizeof(double));
        w_context1 = (double *)malloc(pnet.MAX_vid * dim * sizeof(double));
        for (long vid=0; vid<pnet.MAX_vid; ++vid)
        {
            w_vertex[vid].resize(dim);
            for (int d=0; d<dim;++d)
            {
                w_vertex[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
                w_vertex1[vid*dim+d] = (rand()/(double)RAND_MAX - 0.5) / dim;
            }
        }
        for (long vid=0; vid<pnet.MAX_vid; ++vid)
        {
            w_context[vid].resize(dim);
            for (int d=0; d<dim;++d)
            {
                w_context[vid][d] = 0.0;
                w_context1[vid*dim+d] = 0.0;
            }
        }
    }
}


void LINE::Train(int sample_times, int negative_samples, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[LINE]" << endl;

    cout << "Learning Parameters:" << endl;
    if (order==1)
        cout << "\torder:\t\t\t" << order << "st" << endl;
    else
        cout << "\torder:\t\t\t" << order << "nd" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;

    unsigned long long total_sample_times = (unsigned long long)sample_times*1000000;
    double alpha_min = alpha * 0.0001;
    double alpha_last;
    
    unsigned long long current_sample = 0;
    unsigned long long jobs = total_sample_times/workers;
    
    if (order==1)
    {
        #pragma omp parallel for
        for (int worker=0; worker<workers; ++worker)
        {

            unsigned long long count = 1;
            double _alpha = alpha;
            long v1, v2;

            while (count<jobs)
            {            
                v1 = pnet.SourceSample();
                v2 = pnet.TargetSample(v1);
                pnet.UpdatePair(w_vertex_o1, w_vertex_o1, v1, v2, dim, negative_samples, _alpha);
                //pnet.UpdateCosinePair(w_vertex_o1, w_vertex_o1, v1, v2, dim, negative_samples, _alpha);
                //pnet.UpdateLengthPair(w_vertex_o1, w_vertex_o1, v1, v2, dim, negative_samples, _alpha);

                count++;
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
    else
    {
        #pragma omp parallel for
        for (int worker=0; worker<workers; ++worker)
        {

            unsigned long long count = 1;
            double _alpha = alpha;
            long v1, v2;

            while (count<jobs)
            {            
                v1 = pnet.SourceSample();
                v2 = pnet.TargetSample(v1);
                pnet.UpdatePair(w_vertex, w_context, v1, v2, dim, negative_samples, _alpha);
                //pnet.UpdateCosinePair(w_vertex, w_context, v1, v2, dim, negative_samples, _alpha);
                //pnet.UpdatePair1(w_vertex1, w_context1, v1, v2, dim, negative_samples, _alpha);

                count++;
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

}

