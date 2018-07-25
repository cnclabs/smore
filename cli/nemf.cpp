#define _GLIBCXX_USE_CXX11_ABI 1
#include "../src/model/NEMF.h"

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
        printf("\tcommand nemf interface for proNet-core\n\n");
        printf("Options Description:\n");
        printf("\t-train <string>\n");
        printf("\t\tTrain the Network data\n");
        printf("\t-save <string>\n");
        printf("\t\tSave the representation data\n");
        printf("\t-field <string>\n");
        printf("\t\tField data\n");
        printf("\t-dimensions <int>\n");
        printf("\t\tDimension of vertex representation; default is 64\n");
        printf("\t-sample_times <int>\n");
        printf("\t\tNumber of training samples *Million; default is 10\n");
        printf("\t-walk_steps <int>\n");
        printf("\t\tWalk steps; default is 5\n");
        printf("\t-threads <int>\n");
        printf("\t\tNumber of training threads; default is 1\n");
        printf("\t-alpha <float>\n");
        printf("\t\tInit learning rate; default is 0.025\n");

        printf("Usage:\n");
        printf("\n[NEMF]\n");
        printf("./nemf -train net.txt -save rep.txt -dimensions 64 -sample_times 10 -alpha 0.025 -threads 1\n");

        return 0;
    }
    
    char network_file[100], rep_file[100], field_file[100];
    int dimensions=64, negative_samples=5, sample_times=10, threads=1, walk_steps=5;
    double init_alpha=0.025, reg=0.01;

    if ((i = ArgPos((char *)"-train", argc, argv)) > 0) strcpy(network_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-save", argc, argv)) > 0) strcpy(rep_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-field", argc, argv)) > 0) strcpy(field_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-dimensions", argc, argv)) > 0) dimensions = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-sample_times", argc, argv)) > 0) sample_times = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-walk_steps", argc, argv)) > 0) walk_steps = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-alpha", argc, argv)) > 0) init_alpha = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) threads = atoi(argv[i + 1]);
    
    NEMF *nemf;
    nemf = new NEMF();
    nemf->LoadEdgeList(network_file, 1);
    nemf->LoadFieldMeta(field_file);
    nemf->Init(dimensions);
    nemf->Train(sample_times, walk_steps, init_alpha, threads);
    nemf->SaveWeights(rep_file);

    return 0;


}
