#include "../src/DeepWalk.h"
#include "../src/Walklets.h"
#include "../src/LINE.h"
#include "../src/HPE.h"
#include "../src/FINE.h"

int main(int argc, char **argv){
    
    /*
    DeepWalk *dw;
    dw = new DeepWalk();
    dw->LoadEdgeList(argv[1], 0);
    dw->Init(64);
    
    if (argc == 3)
    {
        dw->Train(100, 40, 1, 5, 0.025, atoi(argv[2]));
    }
    else
    {
        cout << "./cli [FilePath] [# of Threads]";
        return 1;
    }

    dw->SaveWeights("GraphRecDW.model");

    return 0;
    */
    
    /*
    LINE *line;
    line = new LINE();
    line->LoadEdgeList(argv[1], 0);
    line->Init(200);
    
    if (argc == 3)
    {
        line->Train(1000, 5, 0.025, atoi(argv[2]));
        line->SaveWeights("GraphRecLINE.model");
    }

    return 0;

    HPE *hpe;
    hpe = new HPE();
    hpe->LoadEdgeList(argv[1], 0);
    hpe->Init(200);
    
    if (argc == 3)
    {
        hpe->Train(1000, 5, 1, 0.025, atoi(argv[2]));
        hpe->SaveWeights("GraphRecHPE.model");
    }

    return 0;

    Walklets *wl;
    wl = new Walklets();
    wl->LoadEdgeList(argv[1], 0);
    wl->Init(200);

    if (argc == 3)
        wl->Train(10, 40, 2, 3, 5, 0.025, atoi(argv[2]));

    wl->SaveWeights("GraphRecWL.model");

    return 0;
    */

    FINE *fine;
    fine = new FINE();
    fine->LoadEdgeList(argv[1], 1);
    fine->LoadFieldMeta(argv[2]);
    fine->Init(10);
    
    if (argc == 4)
        fine->Train(1000, 5, 1, 0.025, atoi(argv[3]));

    fine->SaveWeights("GraphRecFINE.model");

   return 0;


}
