#include "TChain.h"

int main(int argc, char ** argv){
    if (argc<=1) return 1;
    TChain * ichain = new TChain("t","t");
    ichain->Add(argv[1]);
    int N = ichain->GetEntries();
    printf("%d\n",N);
    return 0;
}
