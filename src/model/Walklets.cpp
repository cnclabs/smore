#include "Walklets.h"

Walklets::Walklets() {}
Walklets::~Walklets() {}

void Walklets::Train(int walk_times, int walk_steps, int window_min, int window_max, int negative_samples, double alpha, int workers){
    
    omp_set_num_threads(workers);

    cout << "Model:" << endl;
    cout << "\t[Walklets]" << endl;

    cout << "Parameters:" << endl;
    cout << "\twalk_times:\t\t" << walk_times << endl;
    cout << "\twalk_steps:\t\t" << walk_steps << endl;
    cout << "\twindow_min:\t\t" << window_min << endl;
    cout << "\twindow_max:\t\t" << window_max << endl;
    cout << "\tnegative_samples:\t" << negative_samples << endl;
    cout << "\talpha:\t\t\t" << alpha << endl;
    cout << "\tworkers:\t\t" << workers << endl;

    cout << "Start Training:" << endl;

    long total = walk_times*pnet.MAX_vid;
    double alpha_min = alpha*0.0001;
    double _alpha = alpha;
    unsigned long long count = 0;

    for (int t=0; t<walk_times; ++t)
    {

        // for random keys access
        std::vector<long> random_keys(pnet.MAX_vid);
        for (long vid = 0; vid < pnet.MAX_vid; vid++) {
             random_keys[vid] = vid;
        }
        for (long vid = 0; vid < pnet.MAX_vid; vid++) {
            int rdx = vid + rand() % (pnet.MAX_vid - vid); // careful here!
            swap(random_keys[vid], random_keys[rdx]);
        }

        #pragma omp parallel for
        for (long vid=0; vid<pnet.MAX_vid; ++vid)
        {

            vector<long> walks = pnet.RandomWalk(vid, walk_steps);
            vector<vector<long>> train_data = pnet.ScaleSkipGrams(walks, window_min, window_max, 0);
            pnet.UpdatePairs(w_vertex, w_context, train_data[0], train_data[1], dim, negative_samples, _alpha);
            
            count++;
            if (count % MONITOR == 0)
            {
                _alpha = alpha* ( 1.0 - (double)(count)/total );
                if (_alpha < alpha_min) _alpha = alpha_min;
                printf("\tAlpha: %.6f\tProgress: %.3f %%%c", _alpha, (double)(count)/total * 100, 13);
                fflush(stdout);
            }

        }

    }
    cout << "\tProgress:\t\t100.00 %\r" << endl;

}

