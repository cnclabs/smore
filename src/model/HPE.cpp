#include "HPE.h"
#include <msgpack.hpp>

HPE::HPE() {}
HPE::~HPE() {}

void HPE::SaveWeights(string model_name, int save_binary){

    cout << "Save Model:" << endl;
    if (save_binary)
    {
        std::string key_file_name = model_name + ".key";
        std::string vec_file_name = model_name + ".vec";
        FILE* key_file = fopen(key_file_name.c_str(), "wb");
        FILE* vec_file = fopen(vec_file_name.c_str(), "wb");
        if (key_file && vec_file)
        {
            msgpack::sbuffer key_buffer;
            msgpack::pack(key_buffer, pnet.vertex_hash.keys);
            std::fwrite(key_buffer.data(), key_buffer.size(), 1, key_file);
            fclose(key_file);
            cout << "\tSave keys to <" << key_file_name << ">" << endl;

            msgpack::sbuffer vec_buffer;
            msgpack::pack(vec_buffer, w_vertex);
            std::fwrite(vec_buffer.data(), vec_buffer.size(), 1, vec_file);
            fclose(vec_file);
            cout << "\tSave vectors to <" << vec_file_name << ">" << endl;
        }
        else
        {
            cout << "\tfail to open file" << endl;
        }
    }
    else
    {
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
        model.close();
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
            w_context[vid][d] = (rand()/(double)RAND_MAX - 0.5) / dim;
    }
}


void HPE::Train(int sample_times, int walk_steps, int negative_samples, double reg, double alpha, int workers){

    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[HPE]" << endl;

    cout << "Learning Parameters:" << endl;
    cout << "\tsample_times:\t\t" << sample_times << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\twalk_steps:\t\t" << walk_steps << endl;
    cout << "\tregularization:\t\t" << reg << endl;
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
            v1 = pnet.SourceSample();
            v2 = pnet.TargetSample(v1);
            pnet.UpdateCommunity(w_vertex, w_context, v1, v2, dim, reg, walk_steps, negative_samples, _alpha);
            pnet.UpdatePair(w_vertex, w_context, v2, v1, dim, negative_samples, _alpha);
            //pnet.UpdateCBOW(w_vertex, w_context, v2, v1, dim, reg, walk_steps, negative_samples, _alpha);

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

