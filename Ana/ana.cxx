#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <sstream>
#include <iostream>

int main(int argc, char** argv){
	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.v2.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bd;
	int wp_ch;
	int wp_wid;
	int wp_lid;
	double wp_xc;
	double wp_yc;
	int map_lid[2][48];
	int map_wid[2][48];
	double map_xc[7][11];
	double map_yc[7][11];
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("b",&wp_bd);
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xc",&wp_xc);
	TTree_wirepos->SetBranchAddress("yc",&wp_yc);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<=7){
			map_xc[wp_lid-1][wp_wid] = wp_xc/10.;
			map_yc[wp_lid-1][wp_wid] = wp_yc/10.;
			map_wid[wp_bd][wp_ch] = wp_wid;
			map_lid[wp_bd][wp_ch] = wp_lid-1;
		}
	}

	//TFile * ifile = new TFile("run.002.try003.root");
	TFile * ifile = new TFile("../root/d_42.root");
	TTree * it = (TTree*) ifile->Get("t");

	int triggerNumber;
	int iHit[96];
	double driftTime[96];

	it->SetBranchAddress("triggerNumber",&triggerNumber);
	it->SetBranchAddress("i",&iHit);
	it->SetBranchAddress("dt",driftTime);

	TFile *ofile = new TFile("output.root","RECREATE");
	TTree * otree = new TTree("t","t");
	int Nhits;
	double meantime;
	otree->Branch("n",&Nhits);
	otree->Branch("mt",&meantime);
	TH1D * h_nhits = new TH1D("nhits","nhits",100,0,50);
	TH1D * h_meantime = new TH1D("meantime","meantime",256,0,1000);
	TH1D * h_ilayer = new TH1D("ilayer","ilayer",18,0,9);

	bool layerhit[9];
	int Nevents = it->GetEntries();
	for ( int i = 0 ; i<it->GetEntries(); i++){
		if (i%100==0) printf("%lf%...\n",(double)i/it->GetEntries()*100);
		it->GetEntry(i);
		Nhits = 0;
		meantime = 0;
		for (int ilayer = 0; ilayer<9; ilayer++){
			layerhit[ilayer] = false;
		}
		for (int ch = 0; ch<96; ch++){
			if (iHit[ch]>=0){
				meantime += driftTime[ch];
				Nhits++;
				int ilayer = map_lid[ch/48][ch%48]+1;
				layerhit[ilayer] = true;
			}
		}
		for (int ilayer = 0; ilayer<9; ilayer++){
			if (layerhit[ilayer]) h_ilayer->Fill(ilayer);
		}
		meantime/=Nhits;
		h_nhits->Fill(Nhits);
		h_meantime->Fill(meantime);
		otree->Fill();
	}
//	h_meantime->Scale(1./Nevents);
//	h_nhits->Scale(1./Nevents);
//	h_ilayer->Scale(1./Nevents);
	h_meantime->Write();
	h_nhits->Write();
	h_ilayer->Write();
	otree->Write();
	ofile->Close();

	return 0;
}
