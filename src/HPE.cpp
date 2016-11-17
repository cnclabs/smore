#include "HPE.h"

HPE::HPE() {}
HPE::~HPE() {}

void HPE::UpdateSecondOrder(long vertex, long context, int negative_samples, int walk_steps, double alpha){
    
    vector< double >* w_vertex_ptr;
    vector< double >* w_context_ptr;
    double* back_err = new double[dim];

    int d;
    long rand_v;
    double label, g, f, rand_p;
    
    w_vertex_ptr = &w_vertex[context];
    w_context_ptr = &w_context[vertex];

    // 0 for postive sample, others for negative sample
    for (int s = -1; s <= walk_steps; s++) {
        label = 1;
        if (s == 0)
        {
            w_vertex_ptr = &w_vertex[vertex];
            w_context_ptr = &w_context[context];
        }
        if (s > 0)
        {
            context = rgraph.TargetSample(context);
            if (context==-1) break;
            w_context_ptr = &w_context[ context ];
        }

        for (d=0; d<dim; ++d)
            back_err[d] = 0.0;
        for (int neg=0; neg<=negative_samples; ++neg)
        {
            // negative sampling
            if (neg!=0){
                label = 0;
                rand_v = rgraph.NegativeSample(); // Negative Sample
                w_context_ptr = &w_context[ rand_v ];
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

    delete [] back_err;

}

void HPE::Train(int sample_times, int negative_samples, int walk_steps, double alpha, int workers){
    
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
        UpdateSecondOrder(v1, v2, negative_samples, walk_steps, _alpha);

    }
    printf("\tAlpha: %.6f\tProgress: 100.00 %%\n", _alpha);

}

/*
int main(int argc, char **argv){
    
    HPE *hpe;
    hpe = new HPE();
    //hpe->LoadEdgeList("../ml-1m/result/ml-1m.graph.train", 0);
    hpe->LoadEdgeList("../../../text8.pair", 0);
    //hpe->LoadEdgeList("in", 0);
    hpe->Init(200);
    
    if (argc == 2)
        hpe->Train(1000, 5, 1, 0.025, atoi(argv[1]));
    else
        hpe->Train(1000, 5, 1, 0.025, 4);

    hpe->SaveWeights("GraphRecHPE.model");

    return 0;
}
*/
