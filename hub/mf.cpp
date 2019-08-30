#define _GLIBCXX_USE_CXX11_ABI 1
#include <omp.h>
#include "../src/util/util.h"                       // arguments
#include "../src/util/file_graph.h"                 // file-based graph reader
#include "../src/sampler/vc_sampler.h"              // sampler
#include "../src/mapper/lookup_mapper.h"            // mapper
#include "../src/optimizer/pairwise_optimizer.h"    // optimizer

int main(int argc, char **argv){
 
    // arguments
    if (argc == 1) {
        printf("[SMORe-MF]\n");
        printf("Options Description:\n");
        printf("\t-train <string>\n");
        printf("\t\tTrain the Network data\n");
        printf("\t-dimension <double>\n");
        printf("\t\tDimension of Embedding\n");
        printf("\t-negative <int>\n");
        printf("\t\tNumber of Negative Samples\n");
        printf("\t-worker <int>\n");
        printf("\t\tNumber of Workers\n");
        printf("\t-update_times <int>\n");
        printf("\t\tNumber of Updates (*million)\n");

        printf("Usage:\n");

        return 0;
    }

    int i, dimension=64, worker=1, negative=5, update_times=10;
    std::string path, model_name="model.rep";
    if ((i = ArgPos((char *)"-train", argc, argv)) > 0) path.assign(argv[i + 1]);
    if ((i = ArgPos((char *)"-save", argc, argv)) > 0) model_name.assign(argv[i + 1]);
    if ((i = ArgPos((char *)"-dimension", argc, argv)) > 0) dimension = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-negative", argc, argv)) > 0) negative = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-update_times", argc, argv)) > 0) update_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-worker", argc, argv)) > 0) worker = atoi(argv[i + 1]);
    
    unsigned long long total_update_times = (unsigned long long)update_times*1000000;
    unsigned long long worker_update_times = total_update_times/worker;
    unsigned long long finished_update_times = 0;
    unsigned long long report_period = 10000;
    Monitor monitor(total_update_times);

    // main
    // 0. [Graph] read from file-based graph
    FileGraph *file_graph = new FileGraph(path, 0);

    // 1. [Sampler] determine what sampler to be used
    VCSampler sampler(file_graph);

    // 2. [Mapper] define what embedding mapper to be used
    LookupMapper mapper(sampler.vertex_size, dimension);
    
    // 3. [Optimizer] claim the optimizer
    PairwiseOptimizer optimizer;

    // 4. building the blocks [MF]
    omp_set_num_threads(worker);
    #pragma omp parallel for
    for (int w=0; w<worker; w++)
    {
        static thread_local double alpha=0.25, decay=alpha*0.9999/total_update_times;
        static thread_local long user, item;
        static thread_local std::vector<double> user_batch_loss(dimension, 0.0);
        static thread_local std::vector<double> item_loss(dimension, 0.0);
        static thread_local unsigned long long update=0;

        while (update < worker_update_times)
        {
            // 4.0 reset user batch loss
            user_batch_loss.assign(dimension, 0.0);
            item_loss.assign(dimension, 0.0);
    
            // 4.1 sample positive (user, item) pair, feed the loss, update
            user = sampler.draw_a_vertex();
            item = sampler.draw_a_context(user);
            optimizer.feed_loglikelihood_loss(mapper[user], mapper[item], 1.0, dimension, user_batch_loss, item_loss);
            mapper.update_with_l2(item, item_loss, alpha, 0.001);
            item_loss.assign(dimension, 0.0);
    
            // 4.2 sampler negative (user, item) pair, feed the loss, update
            for (int n=0; n<negative; n++)
            {
                item = sampler.draw_a_negative();
                optimizer.feed_loglikelihood_loss(mapper[user], mapper[item], 0.0, dimension, user_batch_loss, item_loss);
                mapper.update_with_l2(item, item_loss, alpha, 0.001);
                item_loss.assign(dimension, 0.0);
            }
            mapper.update_with_l2(user, user_batch_loss, alpha, 0.01);
            
            // 5. print progress
            alpha -= decay;
            update++;
            if (update % report_period == 0) {
                finished_update_times += report_period;
                monitor.progress(finished_update_times);
            }
        }
        monitor.end();
    }
    mapper.save_to_file(file_graph->index2vertex, model_name);
    
    return 0;
}
