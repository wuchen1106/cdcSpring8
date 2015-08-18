#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
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
		suffix=suffix+".";
	}

	TH1D * h_t[77][65];
	TF1 * f_t[77][65];
	for (int i = 0; i<65; i++){
		for (int j = 0; j<77; j++){
			h_t[j][i] = new TH1D(Form("h%d_%d_%d",j/11+1,j%11,i),Form("h%d",i),260,-10,250);
			f_t[j][i] = new TF1(Form("f%d_%d_%d",j/11+1,j%11,i),"gaus",-10,250);
		}
	}

	for (int theLayer = 1; theLayer<=7; theLayer++){
		std::cout<<"################################"<<std::endl;
		std::cout<<"       "<<theLayer<<std::endl;
		TFile * ifile = new TFile(Form("../root/t_%d.layer%d.",runNo,theLayer)+suffix+"root");
		if (!ifile) continue;
		TTree * it = (TTree*) ifile->Get("t");

		std::vector<int> * layerID = 0;
		std::vector<int> * wireID = 0;
		std::vector<double> * fitD = 0;
		std::vector<double> * driftT = 0;
		double chi2;

		it->SetBranchAddress("layerID",&layerID);
		it->SetBranchAddress("wireID",&wireID);
		it->SetBranchAddress("fitD",&fitD);
		it->SetBranchAddress("driftT",&driftT);
		it->SetBranchAddress("chi2",&chi2);

		Long64_t N = it->GetEntries();
		if (nEventMax&&nEventMax<N) N = nEventMax;
		std::cout<<"Processing "<<N<<" events..."<<std::endl;
		for ( int i = 0 ; i<N; i++){
			if (i%100==0) printf("%lf%...\n",(double)i/N*100);
			it->GetEntry(i);
			if (chi2>5) continue;
			for (int ihit = 0; ihit<layerID->size(); ihit++){
				if ((*layerID)[ihit]==theLayer){
					int idiv = (*fitD)[ihit]/0.025+32;
					if (idiv>=0&&idiv<65)
						h_t[(theLayer-1)*11+(*wireID)[ihit]][idiv]->Fill((*driftT)[ihit]);
				}
			}
		}
	}

	TFile *ofile = new TFile(Form("../info/xt.%d.",runNo)+suffix+".root","RECREATE");
	TTree * otree = new TTree("t","t");
	int lid;
	int wid;
	double x;
	double t;
	double sig;
	int nent;
	otree->Branch("x",&x);
	otree->Branch("t",&t);
	otree->Branch("lid",&lid);
	otree->Branch("wid",&wid);
	otree->Branch("sig",&sig);
	otree->Branch("n",&nent);

	TCanvas * c = new TCanvas();
	for (int il = 0; il<7; il++){
		for (int iw = 0; iw<11; iw++){
			std::cout<<"##"<<il<<" "<<iw<<std::endl;
			int index = il*11+iw;
			for (int i = 0; i<65; i++){
				nent = h_t[index][i]->Integral();
				if (nent<100) continue;
				int hmax = h_t[index][i]->GetMaximum();
				double tleft = 0;
				double tright = 0;
				for (int ibin = 1; ibin<=260; ibin++){
					if (!tleft&&h_t[index][i]->GetBinContent(ibin)>hmax/3.) tleft = h_t[index][i]->GetBinCenter(ibin);
					if (tleft&&!tright&&h_t[index][i]->GetBinContent(ibin)>hmax/3.){
						tright = h_t[index][i]->GetBinCenter(ibin);
						break;
					}
				}
				lid = il+1;
				wid = iw;
				h_t[index][i]->Fit(Form("f%d_%d_%d",lid,iw,i),"qN0","",tleft,tright);
				h_t[index][i]->Draw();
				f_t[index][i]->Draw("SAME");
				t = f_t[index][i]->GetParameter(1);
				sig = f_t[index][i]->GetParameter(2);
				x = (i-32)*0.025;
				otree->Fill();
				c->SaveAs(Form("dt.%d.%d.%d.%d.png",runNo,lid,iw,i));
			}
		}
	}

	for (int i = 0; i<65; i++){
		for (int j = 0; j<77; j++){
			h_t[j][i]->Write();
		}
	}
	std::cout<<"finish"<<std::endl;
	otree->Write();
	ofile->Close();

	return 0;
}
