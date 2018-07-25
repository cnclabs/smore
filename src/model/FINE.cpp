#include "FINE.h"
#include <omp.h>

FINE::FINE() {}
FINE::~FINE() {}

void FINE::LoadFieldMeta(string filename) {
    pnet.LoadFieldMeta(filename);
}

void FINE::SaveWeights(string model_name){
    
    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        model << pnet.MAX_vid << " " << dim_2 << endl;
        for (long vid=0; vid!=pnet.MAX_vid; vid++)
        {
            model << pnet.vertex_hash.keys[vid];
            for (int f=0; f<2; f++)
            {
                long fvid = pnet.field[vid].vids[f];
                for (int d=0; d<dim_2; ++d)
                    model << " " << w_vertex[ fvid ][d];
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

void FINE::Init(int dimension) {
   
    cout << "Model Setting:" << endl;
    //cout << "\tdimension:\t\t" << dimension*(1+pnet.MAX_field) << endl;
    cout << "\tdimension:\t\t" << dimension << endl;
    dim_1 = int(dimension);
    dim_2 = int(dimension);
    
    //w_vertex_o1.resize(pnet.MAX_vid);
    //w_vertex_o1.resize(pnet.MAX_fvid);
    //w_vertex.resize(pnet.MAX_vid);
    //w_context.resize(pnet.MAX_vid);
    w_vertex.resize(pnet.MAX_fvid);
    w_context.resize(pnet.MAX_fvid);

    for (long vid=0; vid<pnet.MAX_vid; ++vid)
    {
        //w_vertex_o1[vid].resize(dim_1);
        //w_vertex[vid].resize(dim_2);
        //w_context[vid].resize(dim_2);
        
        //for (int d=0; d<dim_1;++d)
        //    w_vertex_o1[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim_1;
        
        //for (int d=0; d<dim_2;++d)
        //    w_vertext[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim_2;
        //for (int d=0; d<dim_2;++d)
        //    w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim_2;
   
        for (auto fvid: pnet.field[vid].vids)
        {
            //w_vertex_o1[fvid].resize(dim_1);
            w_vertex[fvid].resize(dim_2);
            w_context[fvid].resize(dim_2);
            //for (int d=0; d<dim_1;++d)
            //    w_vertex_o1[fvid][d] = (rand()/(double)RAND_MAX - 0.5) / dim_1;
            for (int d=0; d<dim_2;++d)
                w_vertex[fvid][d] = (rand()/(double)RAND_MAX - 0.5) / dim_2;
            for (int d=0; d<dim_2;++d)
                w_context[fvid][d] = (rand()/(double)RAND_MAX - 0.5) / dim_2;
        }
    }
    
}

void FINE::Train(int sample_times, int walk_steps, int negative_samples, double alpha, int workers) {
    
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

    unsigned long long total_sample_times = (unsigned long long)sample_times*1000000;
    double alpha_min = alpha * 0.0001;
    double alpha_last;
    
    unsigned long long current_sample = 0;
    unsigned long long jobs = total_sample_times/workers;

    #pragma omp parallel for
    for (int worker=0; worker<workers; ++worker)
    {
        unsigned long long count = 0;
        double _alpha = alpha;
        long v1, v2;

        while (count<jobs)
        {
            count++;
            if (count % MONITOR == 0)
            {
                current_sample += MONITOR;
                _alpha = alpha* ( 1.0 - (double)(count)/jobs );
                if (_alpha < alpha_min) _alpha = alpha_min;
                alpha_last = _alpha;
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(current_sample)/total_sample_times * 100, 13);
                fflush(stdout);
            }
            
            v1 = pnet.SourceSample();
            v2 = pnet.TargetSample(v1);
            pnet.UpdateFieldCommunity(w_vertex, w_context, v1, v2, dim_2, walk_steps, negative_samples, _alpha);

            //v1 = pnet.SourceSample();
            //v2 = pnet.TargetSample(v1);
            //pnet.UpdateCommunity(w_vertex_o1, w_vertex_o1, v1, v2, dim_1, walk_steps, negative_samples, _alpha);
            //pnet.UpdatePair(w_vertex_o1, w_vertex_o1, v1, v2, dim_1, negative_samples, _alpha);
            //pnet.UpdateFieldCommunity(w_vertex_o1, w_vertex_o1, v1, v2, dim_1, 0, negative_samples, _alpha);
            //pnet.UpdateFieldCommunity(w_vertex_o1, w_vertex_o1, v1, v2, dim, walk_steps, negative_samples, _alpha);
        }

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", alpha_last);

}

