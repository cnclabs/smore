#define _GLIBCXX_USE_CXX11_ABI 1
#include "../src/model/Walklets.h"

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
        printf("\t-negative_samples <int>\n");
        printf("\t\tNumber of negative examples; default is 5\n");
        printf("\t-window_min <int>\n");
        printf("\t\tBegin window of skip-gram; default is 2\n");
        printf("\t-window_max <int>\n");
        printf("\t\tEnd window of skip-gram; default is 5\n");
        printf("\t-walk_times <int>\n");
        printf("\t\tTimes of being staring vertex; default is 10\n");
        printf("\t-walk_steps <int>\n");
        printf("\t\tStep of random walk; default is 40\n");
        printf("\t-threads <int>\n");
        printf("\t\tNumber of training threads; default is 1\n");
        printf("\t-alpha <float>\n");
        printf("\t\tInit learning rate; default is 0.025\n");

        printf("Usage:\n");
        printf("./walklets -train net.txt -save rep.txt -undirected 1 -dimensions 64 -walk_times 10 -walk_steps 40 -window_min 2 -window_max 5 -negative_samples 5 -alpha 0.025 -threads 1\n");

        return 0;
    }
    
    char network_file[100], rep_file[100];
    int dimensions=64, undirected=1, window_min=2, window_max=5, negative_samples=5, walk_times=10, walk_steps=40, threads=1;
    double init_alpha=0.025;

    if ((i = ArgPos((char *)"-train", argc, argv)) > 0) strcpy(network_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-save", argc, argv)) > 0) strcpy(rep_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-undirected", argc, argv)) > 0) undirected = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-dimensions", argc, argv)) > 0) dimensions = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-window_min", argc, argv)) > 0) window_min = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-window_max", argc, argv)) > 0) window_max = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-negative_samples", argc, argv)) > 0) negative_samples = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-walk_times", argc, argv)) > 0) walk_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-walk_steps", argc, argv)) > 0) walk_steps = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-alpha", argc, argv)) > 0) init_alpha = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) threads = atoi(argv[i + 1]);

    int window_left=1, window_right=2;
    Walklets *wl;
    wl = new Walklets();
    wl->LoadEdgeList(network_file, undirected);
    wl->Init(dimensions);
    wl->Train(walk_times, walk_steps, window_min, window_max, negative_samples, init_alpha, threads);
    wl->SaveWeights(rep_file);

   return 0;


}
