#include <stdlib.h>

#include "TChain.h"
#include "TString.h"

int main(int argc, char ** argv){
    if (argc<3) return -1;
    int runNo = atoi(argv[1]);
    TString runname = argv[2];
    int theLayer = 0;
    if (argc==4) theLayer = atoi(argv[3]);
    TChain * iChain_h = new TChain("t","t");
    iChain_h->Add(Form("h_%d.root",runNo));
    int N = iChain_h->GetEntries();

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
        c->GetEntries();
        if (c->GetEntries())
            c->Merge(Form("t_%d.%s.layer%d.root",runNo,runname.Data(),iLayer));
        printf("finished combining t_%d.%s.layer%d.root\n",runNo,runname.Data(),iLayer);
    }
    return 0;
}
