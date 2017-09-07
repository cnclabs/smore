#define _GLIBCXX_USE_CXX11_ABI 1
#include "../src/model/APP.h"

int ArgPos(char *str, int argc, char **argv) {
    int a;
    for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
        if (a == argc - 1) {
            printf("Argument missing for %s\n", str);
            exit(1);
        }
        return a;
    }
    return -1;
}

int main(int argc, char **argv){
    
    int i;

    if (argc == 1) {
        printf("[proNet-core]\n");
        printf("\tcommand line interface for proNet-core\n\n");
        printf("Options Description:\n");
        printf("\t-train <string>\n");
        printf("\t\tTrain the Network data\n");
        printf("\t-save <string>\n");
        printf("\t\tSave the representation data\n");
        printf("\t-dimensions <int>\n");
        printf("\t\tDimension of vertex representation; default is 64\n");
        printf("\t-undirected <int>\n");
        printf("\t\tWhether the edge is undirected; default is 1\n");
        printf("\t-walk_times <int>\n");
        printf("\t\tTimes of being staring vertex; default is 100\n");
        printf("\t-sample_times <int>\n");
        printf("\t\tNumber of paths generated in a walk period; default is 20\n");
        printf("\t-jump <double>\n");
        printf("\t\tjumping facor, the probability to stop a walk; default is 0.15\n");
        printf("\t-negative_samples <int>\n");
        printf("\t\tNumber of negative examples; default is 5\n");
        printf("\t-alpha <float>\n");
        printf("\t\tInit learning rate; default is 0.025\n");
        printf("\t-threads <int>\n");
        printf("\t\tNumber of training threads; default is 1\n");

        printf("Usage:\n");
        printf("./app -train net.txt -save rep.txt -undirected 1 -dimensions 64 -walk_times 100 -sample_times 20 -jump 0.5 -negative_samples 5 -alpha 0.025 -threads 1\n");

        return 0;
    }
    
    char network_file[100], rep_file[100];
    int dimensions=64, undirected=1, window_size=5, negative_samples=5, walk_times=10, sample_times=10, threads=1;
    double init_alpha=0.025, jump=0.15;

    if ((i = ArgPos((char *)"-train", argc, argv)) > 0) strcpy(network_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-save", argc, argv)) > 0) strcpy(rep_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-undirected", argc, argv)) > 0) undirected = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-dimensions", argc, argv)) > 0) dimensions = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-negative_samples", argc, argv)) > 0) negative_samples = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-walk_times", argc, argv)) > 0) walk_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-sample_times", argc, argv)) > 0) sample_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-jump", argc, argv)) > 0) jump = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-alpha", argc, argv)) > 0) init_alpha = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) threads = atoi(argv[i + 1]);

    APP *app;
    app = new APP();
    app->LoadEdgeList(network_file, undirected);
    app->Init(dimensions);
    app->Train(walk_times, sample_times, jump, negative_samples, init_alpha, threads);
    app->SaveWeights(rep_file);

    return 0;


}
