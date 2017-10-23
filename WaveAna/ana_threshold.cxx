#include <iostream>
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"

int main(int argc, char ** argv){
    TFile * ifile = new TFile(argv[1]);
    TTree * itree = (TTree*) ifile->Get("t");
    int nh[50];
    int Nh[50][12]={0};
    itree->SetBranchAddress("nh",nh);

    Long64_t nEntries = itree->GetEntries();
    for (Long64_t iEntry = 0; iEntry<nEntries; iEntry++){
        itree->GetEntry(iEntry);
        if (iEntry%10000==0) std::cout<<iEntry/((double)nEntries)*100<<" %"<<std::endl;
        for (int i = 0; i<50; i++){
            if (nh[i]<=10)
                Nh[i][nh[i]]++;
            else
                Nh[i][11]++;
        }
    }

    std::cout<<"i/I n0/I n1/I n2/I n3/I n4/I n5/I n6/I n7/I n8/I n9/I n10 n11/I"<<std::endl;
    for (int i = 0; i<50; i++){
        std::cout<<i<<" "<<Nh[i][0]<<" "<<Nh[i][1]<<" "<<Nh[i][2]<<" "<<Nh[i][3]<<" "<<Nh[i][4]<<" "<<Nh[i][5]<<" "<<Nh[i][6]<<" "<<Nh[i][7]<<" "<<Nh[i][8]<<" "<<Nh[i][9]<<" "<<Nh[i][10]<<" "<<Nh[i][11]<<std::endl;
    }

    int y[8][50];
    int x[50];
    int ymax = 0;
    for (int i = 0; i<50; i++){
        x[i] = i;
        y[0][i] = Nh[8][i]+Nh[9][i];
        y[1][i] = Nh[0][i]+Nh[1][i]+Nh[2][i]+Nh[3][i]+Nh[4][i]+Nh[5][i];
        y[2][i] = Nh[6][i];
        y[3][i] = Nh[7][i];
        y[4][i] = Nh[8][i];
        y[5][i] = Nh[9][i];
        y[6][i] = Nh[10][i];
        y[7][i] = Nh[11][i];
        if (y[0][i]>ymax) ymax = y[0][i];
    }
    TGraph * g[8];
    int color[8];
    int icolor = 0;
    color[icolor++] = kBlack;
    color[icolor++] = kGreen;
    color[icolor++] = kGreen+2;
    color[icolor++] = kCyan;
    color[icolor++] = kBlue;
    color[icolor++] = kMagenta;
    color[icolor++] = kRed-2;
    color[icolor++] = kRed;
    TCanvas * c = new TCanvas("c","c",1024,768);
    for (int i = 0; i<8; i++){
        g[i] = new TGraph(50,x,y[i]);
        g[i]->SetMarkerStyle(20);
        g[i]->SetMarkerColor(color[i]);
        if (i == 0){
            g[i]->GetYaxis()->SetRangeUser(0,ymax*1.1);
            g[i]->Draw("ALP");
        }
        else{
            g[i]->Draw("LPSAME");
        }
    }
    c->SaveAs("threshold.pdf");
    c->SaveAs("threshold.png");

    return 0;
}
