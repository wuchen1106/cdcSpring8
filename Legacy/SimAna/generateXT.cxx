#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TRandom.h"
#include <iostream>

void smearx(double &x);
void smeart(double &t);
double x2t(double x);

TF1 * f_resx;
TF1 * f_rest;

int main(int argc, char** argv)
{
	TFile * ofile = new TFile("output.root","RECREATE");
	TTree * otree = new TTree("t","t");
	double x,t;
	otree->Branch("x",&x);
	otree->Branch("t",&t);

	f_resx = new TF1("f_resx","gaus",-2,2);
	f_resx->SetParameter(0,1);
	f_resx->SetParameter(1,0);
	f_resx->SetParameter(2,0.2);

	f_rest = new TF1("f_rest","gaus",-1,1);
	f_rest->SetParameter(0,1);
	f_rest->SetParameter(1,0);
	f_rest->SetParameter(2,10);

	TRandom random;

	std::cout<<"Start!"<<std::endl;
	for (int iev = 0; iev<1e6; iev++){
		if (iev%1000==0) std::cout<<iev/1.e6*100<<"%..."<<std::endl;
		x = random.Uniform(8);
		t = x2t(x);
		smearx(x);
//		smeart(t);
		t+=random.Uniform(2)-1;
		otree->Fill();
	}
	otree->Write();
	ofile->Close();

	return 0;
}

double x2t(double x){
	return x*250./8;
}

void smearx(double &x){
	x+=f_resx->GetRandom();
}

void smeart(double &t){
	t+=f_rest->GetRandom();
}
