#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF1.h"
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define NLAY 9

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

	TFile * ifile = new TFile(Form("../root/h_%d",runNo)+suffix+".root");
	TTree * it = (TTree*) ifile->Get("t");

	std::vector<int> * layerID = 0;
	std::vector<int> * wireID = 0;
	std::vector<double> * driftT = 0;

	it->SetBranchAddress("layerID",&layerID);
	it->SetBranchAddress("wireID",&wireID);
	it->SetBranchAddress("driftT",&driftT);

	TFile *ofile = new TFile("output.root","RECREATE");
	TTree * ot= new TTree("t","t");
	int nLayersWith1GoodHit;
	int nLayersWith1OKHit;
	int nLayersWith1Hit;
	int nHitsGood;
	int nHitsOK;
	int nHits;
	int nHitsAndNoise;
	ot->Branch("n1g",&nLayersWith1GoodHit);
	ot->Branch("n1o",&nLayersWith1OKHit);
	ot->Branch("n1",&nLayersWith1Hit);
	ot->Branch("N",&nHits);
	ot->Branch("No",&nHitsOK);
	ot->Branch("Ng",&nHitsGood);
	ot->Branch("Na",&nHitsAndNoise);

	int isgood[NLAY];
	int isok[NLAY];
	int isall[NLAY];

	Long64_t N = it->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" events..."<<std::endl;
	for ( int i = 0 ; i<N; i++){
		if (i%1000==0) printf("%lf%...\n",(double)i/N*100);
		it->GetEntry(i);
		nLayersWith1GoodHit=0;
		nLayersWith1OKHit=0;
		nLayersWith1Hit=0;
		nHitsGood=0;
		nHitsOK=0;
		nHits=0;
		nHitsAndNoise=0;
		for (int ilayer = 0; ilayer<NLAY; ilayer++){
			isgood[ilayer] = 0;
			isok[ilayer] = 0;
			isall[ilayer] = 0;
		}
		for (int ihit = 0; ihit<layerID->size(); ihit++){
			if ((*layerID)[ihit]==9||(*layerID)[ihit]==5) continue;
			if (fabs((*driftT)[ihit])>20&&fabs((*driftT)[ihit])<210){
				isgood[(*layerID)[ihit]]++;
				nHitsGood++;
			}
			if (fabs((*driftT)[ihit])>5&&fabs((*driftT)[ihit])<220){
				isok[(*layerID)[ihit]]++;
				nHitsOK++;
			}
			if (fabs((*driftT)[ihit])>-10&&fabs((*driftT)[ihit])<235.3){
				isall[(*layerID)[ihit]]++;
				nHits++;
			}
			nHitsAndNoise++;
		}
		for (int ilayer = 0; ilayer<NLAY; ilayer++){
			if (isgood[ilayer]==1) nLayersWith1GoodHit++;
			if (isall[ilayer]==1) nLayersWith1Hit++;
			if (isok[ilayer]==1) nLayersWith1OKHit++;
		}
		ot->Fill();
	}
	ot->Write();
	ofile->Close();
	return 0;
}
