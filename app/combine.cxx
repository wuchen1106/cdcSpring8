#include <stdlib.h>

#include "TChain.h"
#include "TString.h"

int main(int argc, char ** argv){
    if (argc<3) return -1;
    int runNo = atoi(argv[1]);
    TString runname = argv[2];
    int theLayer = 0;
    if (argc==4) theLayer = atoi(argv[3]);
    int N = 0;
    if (runNo == 1012) N = 493190;
    else if (runNo == 117) N = 605460;
    else if (runNo == 118) N = 769003;

    int iLayerStart = 1;
    int iLayerStop = 8;
    if (theLayer){
        iLayerStart = theLayer;
        iLayerStop = theLayer;
    }
    for (int iLayer = iLayerStart; iLayer<=iLayerStop; iLayer++){
        TChain * c = new TChain("t");
        for (int i = 0; i<N; i+=5000){
            c->Add(Form("t_%d.%s.%d-%d.layer%d.root",runNo,runname.Data(),i,i+4999>=N?N-1:i+4999,iLayer));
        }
        if (c->GetEntries())
            c->Merge(Form("t_%d.%s.layer%d.root",runNo,runname.Data(),iLayer));
    }
    return 0;
}
