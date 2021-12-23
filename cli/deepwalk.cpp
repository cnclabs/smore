#define _GLIBCXX_USE_CXX11_ABI 1
#include "../src/model/DeepWalk.h"

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
        printf("\t-load_v <string>\n");
        printf("\t\tPath to vertex pretrain\n");
        printf("\t-load_c <string>\n");
        printf("\t\tPath to context pretrain\n");
        printf("\t-dimensions <int>\n");
        printf("\t\tDimension of vertex representation; default is 64\n");
        printf("\t-undirected <int>\n");
        printf("\t\tWhether the edge is undirected; default is 1\n");
        printf("\t-negative_samples <int>\n");
        printf("\t\tNumber of negative examples; default is 5\n");
        printf("\t-window_size <int>\n");
        printf("\t\tSize of skip-gram window; default is 5\n");
        printf("\t-walk_times <int>\n");
        printf("\t\tTimes of being staring vertex; default is 10\n");
        printf("\t-walk_steps <int>\n");
        printf("\t\tStep of random walk; default is 40\n");
        printf("\t-threads <int>\n");
        printf("\t\tNumber of training threads; default is 1\n");
        printf("\t-alpha <float>\n");
        printf("\t\tInit learning rate; default is 0.025\n");

        printf("Usage:\n");
        printf("./deepwalk -train net.txt -save rep.txt -undirected 1 -dimensions 64 -walk_times 10 -walk_steps 40 -window_size 5 -negative_samples 5 -alpha 0.025 -threads 1\n");

        return 0;
    }

    char network_file[100], rep_file[100], load_v[100]="", load_c[100]="";
    int dimensions=64, undirected=1, window_size=5, negative_samples=5, walk_times=10, walk_steps=5, sample_times=10, threads=1;
    double init_alpha=0.025, bfs=0.0;

    if ((i = ArgPos((char *)"-train", argc, argv)) > 0) strcpy(network_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-save", argc, argv)) > 0) strcpy(rep_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-load_v", argc, argv)) > 0) strcpy(load_v, argv[i + 1]);
    if ((i = ArgPos((char *)"-load_c", argc, argv)) > 0) strcpy(load_c, argv[i + 1]);
    if ((i = ArgPos((char *)"-undirected", argc, argv)) > 0) undirected = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-dimensions", argc, argv)) > 0) dimensions = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-window_size", argc, argv)) > 0) window_size = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-negative_samples", argc, argv)) > 0) negative_samples = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-walk_times", argc, argv)) > 0) walk_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-walk_steps", argc, argv)) > 0) walk_steps = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-alpha", argc, argv)) > 0) init_alpha = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) threads = atoi(argv[i + 1]);

    DeepWalk *dw;
    dw = new DeepWalk();
    dw->LoadEdgeList(network_file, undirected);
    dw->Init(dimensions, load_v, load_c);
    dw->Train(walk_times, walk_steps, window_size, negative_samples, init_alpha, threads);
    dw->SaveWeights(rep_file);

    return 0;


}
