#define _GLIBCXX_USE_CXX11_ABI 1
#include "../src/model/TEXTGCNdev.h"

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
        printf("\t-field <string>\n");
        printf("\t\tField data\n");
        printf("\t-dimensions <int>\n");
        printf("\t\tDimension of vertex representation; default is 64\n");
        printf("\t-undirected <int>\n");
        printf("\t\tWhether the edge is undirected; default is 1\n");
        printf("\t-num_events <int>\n");
        printf("\t\tNumber of events; default is 1\n");
        printf("\t-num_words <int>\n");
        printf("\t\tNumber of words; default is 5\n");
        printf("\t-sample_times <int>\n");
        printf("\t\tNumber of training samples *Million; default is 5\n");
        printf("\t-threads <int>\n");
        printf("\t\tNumber of training threads; default is 1\n");
        printf("\t-reg <float>\n");
        printf("\t\tRegularization term; default is 0.01\n");
        printf("\t-alpha <float>\n");
        printf("\t\tInit learning rate; default is 0.025\n");

        printf("Usage:\n");
        printf("./textgcn -train net.txt -field field.txt -save rep.txt -undirected 0 -dimensions 64 -reg 0.01 -sample_times 5 -walk_steps 5 -negative_samples 5 -alpha 0.025 -threads 1\n");

        return 0;
    }
    
    char network_file[100], rep_file[100], field_file[100];
    int dimensions=64, undirected=1, num_events=1, num_words=5, sample_times=10, threads=1;
    double init_alpha=0.025, reg=0.01;

    if ((i = ArgPos((char *)"-train", argc, argv)) > 0) strcpy(network_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-save", argc, argv)) > 0) strcpy(rep_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-field", argc, argv)) > 0) strcpy(field_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-undirected", argc, argv)) > 0) undirected = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-dimensions", argc, argv)) > 0) dimensions = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-num_events", argc, argv)) > 0) num_events = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-num_words", argc, argv)) > 0) num_words = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-sample_times", argc, argv)) > 0) sample_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-reg", argc, argv)) > 0) reg = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-alpha", argc, argv)) > 0) init_alpha = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) threads = atoi(argv[i + 1]);

    TEXTGCNdev *textgcn;
    textgcn = new TEXTGCNdev();
    textgcn->LoadEdgeList(network_file, undirected);
    textgcn->LoadFieldMeta(field_file);
    textgcn->Init(dimensions);
    textgcn->Train(sample_times, num_events, num_words, reg, init_alpha, threads);
    textgcn->SaveWeights(rep_file);

   return 0;


}
