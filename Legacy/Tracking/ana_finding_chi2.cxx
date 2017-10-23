#include "TChain.h"

int main(int argc, char ** argv){
    TChain * ichain  = new TChain("t","t");
    ichain->Add(argv[1]);
    int iev;
    int icombi;
    int ifit;
    int nGood = 5;
    double chi2z;
    double chi2x;
    ichain->SetBranchAddress("i",&iev);
    ichain->SetBranchAddress("c",&icombi);
    ichain->SetBranchAddress("f",&ifit);
    ichain->SetBranchAddress("n",&nGood);
    ichain->SetBranchAddress("cz",&chi2z);
    ichain->SetBranchAddress("cx",&chi2x);
    int pre_iev = -1;
    int min_icombi[3];
    double min_chi2z[3];
    double min2_chi2z[3];
    double min_chi2x[3];
    double min2_chi2x[3];
    for (int i = 0; i<3; i++){
        min_icombi[i] = -1;
        min_chi2z[i] = 1e9;
        min2_chi2z[i] = 1e9;
        min_chi2x[i] = 1e9;
        min2_chi2x[i] = 1e9;
    }
    printf("i/I ng/I ic1/I ic2/I ic3/I cz_1st_1 cz_1st_2 cz_1st_3 cz_2nd_1 cz_2nd_2 cz_2nd_3 cx_1st_1 cx_1st_2 cx_1st_3 cx_2nd_1 cx_2nd_2 cx_2nd_3\n");
    int max_nGood[3] = {3,3,3};
    for (int i = 0; i<ichain->GetEntries(); i++){
        ichain->GetEntry(i);
        // FIXME: only valid for new runs. In order to compare with old runs
        chi2x*=nGood;
        chi2z*=nGood;
        if (iev==pre_iev){
            if (nGood>max_nGood[ifit]){
                max_nGood[ifit] = nGood;
                min_chi2z[ifit] = chi2z;
                min_chi2x[ifit] = chi2x;
                min_icombi[ifit] = icombi;
                min2_chi2z[ifit] = 1e9;
                min2_chi2x[ifit] = 1e9;
            }
            else if (nGood<max_nGood[ifit]) continue;
            else{
                if (min_chi2z[ifit] > chi2z){
                    min_chi2z[ifit] = chi2z;
                    min_chi2x[ifit] = chi2x;
                    min_icombi[ifit] = icombi;
                }
                else if (min2_chi2z[ifit] > chi2z){
                    min2_chi2z[ifit] = chi2z;
                    min2_chi2x[ifit] = chi2x;
                }
            }
        }
        else{
            if (i!=0) printf("%d %d %d %d %d %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",pre_iev,max_nGood[2],min_icombi[0],min_icombi[1],min_icombi[2],min_chi2z[0],min_chi2z[1],min_chi2z[2],min2_chi2z[0],min2_chi2z[1],min2_chi2z[2],min_chi2x[0],min_chi2x[1],min_chi2x[2],min2_chi2x[0],min2_chi2x[1],min2_chi2x[2]);
            for (int ifit = 0; ifit<3; ifit++){
                max_nGood[ifit] = 3;
                min_chi2z[ifit] = 1e9;
                min2_chi2z[ifit] = 1e9;
                min_chi2x[ifit] = 1e9;
                min2_chi2x[ifit] = 1e9;
            }
            max_nGood[ifit] = nGood;
            min_chi2z[ifit] = chi2z;
            min_chi2x[ifit] = chi2x;
            min_icombi[ifit] = icombi;
            pre_iev = iev;
        }
    }
    return 0;
}
