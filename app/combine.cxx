#include <stdlib.h>

#include "TChain.h"
#include "TString.h"

#include "header.h"

int main(int argc, char ** argv){
    if (argc<3) return -1;
    int runNo = atoi(argv[1]);
    TString runname = argv[2];
    int nPerRun = 0;
    if (argc>=4) nPerRun = atoi(argv[3]);
    int theLayer = 0;
    if (argc>=5) theLayer = atoi(argv[4]);
    TChain * iChain_h = new TChain("t","t");
    iChain_h->Add(Form("h_%d.root",runNo));
    int N = iChain_h->GetEntries();

    int iLayerStart = 0;
    int iLayerStop = NLAY-1;
    if (theLayer){
        iLayerStart = theLayer;
        iLayerStop = theLayer;
    }
    for (int iLayer = iLayerStart; iLayer<=iLayerStop; iLayer++){
        TChain * c = new TChain("t");
        if (nPerRun){
            for (int i = 0; i<N; i+=nPerRun){
                int j = i+nPerRun-1;
                if (j>=N) j = N-1;
                c->Add(Form("t_%d.%s.%d-%d.layer%d.root",runNo,runname.Data(),i,j,iLayer));
            }
        }
        else{
            c->Add(Form("t_%d.%s.*-*.layer%d.root",runNo,runname.Data(),iLayer));
        }
        c->GetEntries();
        if (c->GetEntries())
            c->Merge(Form("t_%d.%s.layer%d.root",runNo,runname.Data(),iLayer));
        printf("finished combining t_%d.%s.layer%d.root\n",runNo,runname.Data(),iLayer);
    }
    return 0;
}
