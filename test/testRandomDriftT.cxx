#include <iostream>

#include <TRandom.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TFile.h>

#include "XTManager.hxx"

int main(int argc, char ** argv){
    if (argc<2){
        std::cout<<"Should provide a file containing the 2D dt relation h2_dt_*"<<std::endl;
        return 1;
    }
    XTManager::Get().SetInputFileXT(argv[1]);
    XTManager::Get().Initialize();

    TFile * ofile = new TFile("testRandomDriftT.root","recreate");
    TH2D * histXT = new TH2D("hdt","",792,-23.4375,801.5625,1001,-10.01,10.01);
    for (int iEvent = 0; iEvent<1000000; iEvent++){
        if (iEvent%1000==0) std::cout<<iEvent<<std::endl;
        double doca = gRandom->Uniform(-8,8);
        double driftT = XTManager::Get().RandomDriftT(doca,4,4);
        histXT->Fill(driftT,doca);
    }

    TCanvas * canv = new TCanvas("canv","",1024,768);
    canv->Divide(2,1);
    canv->cd(1);
    XTManager::Get().GetXTHistDefault()->Draw("COLZ");
    canv->cd(2);
    histXT->Draw("COLZ");
    canv->SaveAs("testRandomDriftT.png");

    return 0;
}
