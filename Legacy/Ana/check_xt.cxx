#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TChain.h"
#include "TF1.h"
#include <sstream>
#include <stdlib.h>
#include <iostream>

int main(int argc, char** argv){

	if (argc<2){
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	//if (argc>=4) nEventMax = (int)strtol(argv[3],NULL,10);
	TString suffix = "";
	if (argc>=3){
		suffix  = argv[2];
		suffix="."+suffix;
	}

	TChain * ichain = new TChain("t","t");
	ichain->AddFile(Form("../info/xt.%d",runNo)+suffix+".root");
	ichain->SetMarkerStyle(7);
	ichain->SetMarkerColor(kBlack);
	ichain->SetLineColor(kBlack);

	TCanvas * canvas = new TCanvas("c","c",1024,768);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	gStyle->SetStatX(0.6);

	for (int theLayer = 1; theLayer<=8; theLayer++){
		TFile * ifile = new TFile(Form("../root/t_%d.layer%d",runNo,theLayer)+suffix+".root");
		if (!ifile) continue;
		TTree * it = (TTree*) ifile->Get("t");
		if (!it) continue;
		// FIXME
		for (int theWire = 2; theWire<=6; theWire++){
//		for (int theWire = 3+(theLayer>3); theWire<=3+(theLayer>4); theWire++){
			std::cout<<"################################"<<std::endl;
			std::cout<<"       "<<theLayer<<", "<<theWire<<std::endl;
			it->Draw("driftT:fitD>>h(500,-9,9,260,-10,250)",Form("chi2<10&&layerID==%d&&wireID==%d&&abs(slZ)<0.145",theLayer,theWire),"COLZ");
			it->SetMarkerStyle(7);
			it->SetMarkerColor(kRed);
			it->Draw("driftT:driftD",Form("chi2<3&&layerID==%d&&wireID==%d&&abs(slZ)<0.17",theLayer,theWire),"SAME");
			ichain->Draw("t:x",Form("lid==%d&&wid==%d",theLayer,theWire),"LPSAME");
			canvas->SaveAs(Form("xt.%d.%d.%d.png",runNo,theLayer,theWire));
		}
	}

	return 0;
}
