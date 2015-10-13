#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TChain.h"
#include "TF1.h"
#include <sstream>
#include <stdlib.h>
#include <math.h>
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
	int theLayer = 7;

	std::vector<int> colors;
	colors.push_back(kBlack);
	colors.push_back(kGreen+1);
	colors.push_back(kCyan+1);
	colors.push_back(kBlue);
	colors.push_back(kRed);
	colors.push_back(kMagenta);
	colors.push_back(kYellow+2);
	colors.push_back(kOrange-3);

	TH1D * h_deltaslz = new TH1D("deltaslz","deltaslz",128,0,0.2);
	TH1D * h_slz[8];
	for (int ilayer = 1; ilayer<=7; ilayer++){
		if (ilayer == theLayer) continue;
		h_slz[ilayer] = new TH1D(Form("slz%d",ilayer),Form("slz%d",ilayer),128,-0.2,0.2);
	}
	TChain * chains[8];
	double slZ[8];
	double chi2[8];
	Long64_t N = 0;
	for (int ilayer = 1; ilayer<=7; ilayer++){
		if (ilayer == theLayer) continue;
		chains[ilayer] = new TChain("t");
		chains[ilayer]->Add(Form("../root/t_%d.layer%d",runNo,ilayer)+suffix+".root");
		chains[ilayer]->SetBranchAddress("slZ",&(slZ[ilayer]));
		chains[ilayer]->SetBranchAddress("chi2",&(chi2[ilayer]));
		N = chains[ilayer]->GetEntries();
	}

	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" events..."<<std::endl;
	for ( int i = 0 ; i<N; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/N*100);
		for (int ilayer = 1; ilayer<=7; ilayer++){
			if (ilayer == theLayer) continue;
			chains[ilayer]->GetEntry(i);
		}
		int nlayers = 0;
		double mean = 0;
		for (int ilayer = 1; ilayer<=7; ilayer++){
			if (ilayer == theLayer) continue;
			if (chi2[ilayer]<3){
				h_slz[ilayer]->Fill(slZ[ilayer]);
				mean+=slZ[ilayer];
				nlayers++;
			}
		}
		if (nlayers<=2) continue;
		mean/=nlayers;
		double rms = 0;
		for (int ilayer = 1; ilayer<=7; ilayer++){
			if (ilayer == theLayer) continue;
			if (chi2[ilayer]<3){
				rms+=pow(slZ[ilayer]-mean,2);
			}
		}
		rms/=nlayers;
		h_deltaslz->Fill(sqrt(rms));
	}

	TFile *ofile = new TFile(Form("../info/hists.%d",runNo)+suffix+".root","RECREATE");
	TTree * otree = new TTree("t","t");
	h_deltaslz->Write();
	TCanvas * canvas = new TCanvas("c","c",1024,768);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	bool notdrawn = true;
	double currentMaximum = 0;
	for (int ilayer = 1; ilayer<=7; ilayer++){
		if (ilayer == theLayer) continue;
		if (h_slz[ilayer]->GetMaximum()>currentMaximum) currentMaximum=h_slz[ilayer]->GetMaximum();
	}
	for (int ilayer = 1; ilayer<=7; ilayer++){
		if (ilayer == theLayer) continue;
		h_slz[ilayer]->Write();
		h_slz[ilayer]->SetLineColor(colors[ilayer]);
		if (notdrawn){
			h_slz[ilayer]->GetYaxis()->SetRangeUser(0,currentMaximum*1.2);
			h_slz[ilayer]->Draw();
			notdrawn = false;
		}
		else{
			h_slz[ilayer]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("hists.%d",runNo)+suffix+".png");
	otree->Write();
	ofile->Close();

	return 0;
}
