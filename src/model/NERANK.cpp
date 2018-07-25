#include "NERANK.h"
#include <omp.h>

NERANK::NERANK() {
    char vertex_method[15] = "out_degrees";
    pnet.SetVertexMethod(vertex_method);
    char negative_method[15] = "degrees";
    pnet.SetNegativeMethod(negative_method);
}
NERANK::~NERANK() {
}

void NERANK::LoadEdgeList(string filename, bool undirect) {
    pnet.LoadEdgeList(filename, undirect);
}

void NERANK::LoadFieldMeta(string filename) {
    pnet.LoadFieldMeta(filename);
}

void NERANK::SaveWeights(string model_name){
    
    cout << "Save Model:" << endl;
    ofstream model(model_name);
    if (model)
    {
        model << pnet.MAX_vid << " " << dim << endl;
        for (long vid=0; vid!=pnet.MAX_vid; vid++)
        {
            model << pnet.vertex_hash.keys[vid];
            
            if (pnet.field[vid].fields[0]==0)
                for (int d=0; d<dim; ++d)
                    model << " " << w_vertexU[vid][d];
            if (pnet.field[vid].fields[0]==1)
                for (int d=0; d<dim; ++d)
                    model << " " << w_vertexI[vid][d];

            model << endl;
        }
        cout << "\tSave to <" << model_name << ">" << endl;
    }
    else
    {
        cout << "\tfail to open file" << endl;
    }
}

void NERANK::Init(int dim) {
   
    this->dim = dim;
    cout << "Model Setting:" << endl;
    cout << "\tdimension:\t\t" << dim << endl;

    w_vertexU.resize(pnet.MAX_vid);
    w_vertexI.resize(pnet.MAX_vid);
    w_contextU.resize(pnet.MAX_vid);
    w_contextI.resize(pnet.MAX_vid);
    w_contextUU.resize(pnet.MAX_vid);
    w_contextII.resize(pnet.MAX_vid);

    for (long vid=0; vid<pnet.MAX_vid; ++vid)
    {
        w_vertexU[vid].resize(dim);
        w_vertexI[vid].resize(dim);
        for (int d=0; d<dim;++d)
        {
            w_vertexU[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
            w_vertexI[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
        }
        w_contextU[vid].resize(dim);
        w_contextI[vid].resize(dim);
        w_contextUU[vid].resize(dim);
        w_contextII[vid].resize(dim);
        for (int d=0; d<dim;++d)
        {
            w_contextU[vid][d] = 0.0;
            w_contextI[vid][d] = 0.0;
            w_contextUU[vid][d] = 0.0;
            w_contextII[vid][d] = 0.0;
        }
    }

}


void NERANK::Train(int sample_times, int walk_steps, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[NERANK]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\twalk steps:\t\t" << walk_steps << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training[1]:" << endl;

    unsigned long long total_sample_times = (unsigned long long)sample_times*1000000;
    double alpha_min = alpha * 0.0001;
    double alpha_last;
    
    unsigned long long current_sample = 0;
    unsigned long long jobs = total_sample_times/workers;
    int negative_samples = 5;
    
    #pragma omp parallel for
    for (int worker=0; worker<workers; ++worker)
    {
        
        long vid, cid, cid2, nid;
        unsigned long long count = 0;
        double _alpha = alpha;
        double margin;

        while (count<jobs)
        {
            vid = pnet.SourceSample();
            while (pnet.field[vid].fields[0]!=0)
                vid = pnet.SourceSample();
            cid = pnet.TargetSample(vid);
            
            {
                pnet.UpdateBatchCommunity(w_vertexI, w_contextI, cid, vid, dim, 0.0, walk_steps, negative_samples, _alpha*0.05);
                pnet.UpdateBatchCommunity(w_vertexU, w_contextU, vid, cid, dim, 0.0, walk_steps, negative_samples, _alpha*0.05);
                pnet.UpdateUIPair(w_vertexU, w_vertexI, w_contextU, w_contextI, vid, cid, dim, 0.025, walk_steps, negative_samples, _alpha);
            }

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

