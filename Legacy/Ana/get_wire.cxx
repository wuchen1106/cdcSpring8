#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraph.h"
#include "TChain.h"
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

	// Load pre XT
	double tres = 10;
	double driftTMax = 0; // ns
	std::vector<double> v_x_all_right;
	std::vector<double> v_t_all_right;
	std::vector<double> v_x_all_left;
	std::vector<double> v_t_all_left;
	TGraph * g_xt_all_right = 0;
	TGraph * g_xt_all_left = 0;
	TF1 * f_xt_all_right = 0;
	TF1 * f_xt_all_left = 0;
	TChain * c_xt = new TChain("t","t");
	c_xt->Add(Form("../info/xt.%d.root",runNo));
	double i_xt_x, i_xt_t;
	int i_xt_n;
	c_xt->SetBranchAddress("x",&i_xt_x);
	c_xt->SetBranchAddress("t",&i_xt_t);
	for ( int i = 0; i<c_xt->GetEntries(); i++ ){
		c_xt->GetEntry(i);
		if (i_xt_x>=-0.1){
			v_x_all_right.push_back(i_xt_x);
			v_t_all_right.push_back(i_xt_t);
		}
		if (i_xt_x<=0.1){
			v_x_all_left.push_back(i_xt_x);
			v_t_all_left.push_back(i_xt_t);
		}
		if (i_xt_t>driftTMax) driftTMax=i_xt_t;
	}
	double temp;
	for (int i = 0; i<v_x_all_left.size(); i++){
		for (int j = i+1; j<v_x_all_left.size(); j++){
			if (v_t_all_left[j]<v_t_all_left[i]){
				temp = v_t_all_left[i];
				v_t_all_left[i]=v_t_all_left[j];
				v_t_all_left[j]=temp;
				temp = v_x_all_left[i];
				v_x_all_left[i]=v_x_all_left[j];
				v_x_all_left[j]=temp;
			}
		}
	}
	for (int i = 0; i<v_x_all_right.size(); i++){
		for (int j = i+1; j<v_x_all_right.size(); j++){
			if (v_t_all_right[j]<v_t_all_right[i]){
				temp = v_t_all_right[i];
				v_t_all_right[i]=v_t_all_right[j];
				v_t_all_right[j]=temp;
				temp = v_x_all_right[i];
				v_x_all_right[i]=v_x_all_right[j];
				v_x_all_right[j]=temp;
			}
		}
	}
	f_xt_all_right = new TF1("f_xt_all_right","pol9",0-tres,driftTMax);
	f_xt_all_left = new TF1("f_xt_all_left","pol9",0-tres,driftTMax);
	g_xt_all_right = new TGraph(v_x_all_right.size(),&(v_t_all_right[0]),&(v_x_all_right[0]));
	g_xt_all_left = new TGraph(v_x_all_left.size(),&(v_t_all_left[0]),&(v_x_all_left[0]));
	g_xt_all_right->Fit("f_xt_all_right","qN0","");
	g_xt_all_left->Fit("f_xt_all_left","qN0","");
	TCanvas * c1 = new TCanvas("c1","c1");
	g_xt_all_left->Draw("ALP");
	f_xt_all_left->Draw("SAME");
	c1->SaveAs("left.pdf");
	g_xt_all_right->Draw("ALP");
	f_xt_all_right->Draw("SAME");
	c1->SaveAs("right.pdf");

	// Check post XT, find delta
	double offset[8] = {0};
	int npoints[8] = {0};
	TChain * chain_postXT = new TChain("t");
	chain_postXT->Add(Form("../info/xt.%d",runNo)+suffix+".root");
	double postXT_x;
	double postXT_t;
	int postXT_lid;
	int postXT_n;
	chain_postXT->SetBranchAddress("x",&postXT_x);
	chain_postXT->SetBranchAddress("t",&postXT_t);
	chain_postXT->SetBranchAddress("lid",&postXT_lid);
	chain_postXT->SetBranchAddress("n",&postXT_n);
	int N = chain_postXT->GetEntries();
	TF1 * f;
	for (int i = 0; i<=N; i++){
		chain_postXT->GetEntry(i);
		if (postXT_n<800) continue;
		if (postXT_x>=0)
			f=f_xt_all_right;
		else
			f=f_xt_all_left;
		double fitD = f->Eval(postXT_t);
		if (postXT_t>driftTMax+tres) continue;
		else if (postXT_t>driftTMax) fitD = postXT_x>=0?8:-8;
		offset[postXT_lid]+=postXT_x-fitD;
	//	std::cout<<postXT_lid<<": "<<postXT_x<<"-("<<postXT_t<<")"<<fitD<<"="<<postXT_x-fitD<<std::endl;
		npoints[postXT_lid]++;
	}
	for (int ilayer = 1; ilayer<=7; ilayer++){
		if (npoints[ilayer]) offset[ilayer]/=npoints[ilayer];
	}
	for (int ilayer = 1; ilayer<=7; ilayer++){
		if (ilayer==4) continue;
		offset[ilayer]-=offset[4];
//		std::cout<<ilayer<<": "<<offset[ilayer]<<std::endl;
		offset[ilayer]/=2;
	}
	offset[4] = 0;

	// Read pre wire position and save new one
	double xhv,yhv,xc,yc,xro,yro;
	int b,ch,l,h,w;
	TChain * chain_preWire = new TChain("t");
	chain_preWire->Add("../info/wire-position.root");
	chain_preWire->SetBranchAddress("xhv",&xhv);
	chain_preWire->SetBranchAddress("yhv",&yhv);
	chain_preWire->SetBranchAddress("xc",&xc);
	chain_preWire->SetBranchAddress("yc",&yc);
	chain_preWire->SetBranchAddress("xro",&xro);
	chain_preWire->SetBranchAddress("yro",&yro);
	chain_preWire->SetBranchAddress("b",&b);
	chain_preWire->SetBranchAddress("ch",&ch);
	chain_preWire->SetBranchAddress("l",&l);
	chain_preWire->SetBranchAddress("w",&w);
	chain_preWire->SetBranchAddress("h",&h);
	TFile * outFile =new TFile("../info/wire-position"+suffix+".root","RECREATE");
	TTree * outTree = new TTree("t","t");
	outTree->Branch("xhv",&xhv);
	outTree->Branch("yhv",&yhv);
	outTree->Branch("xc",&xc);
	outTree->Branch("yc",&yc);
	outTree->Branch("xro",&xro);
	outTree->Branch("yro",&yro);
	outTree->Branch("b",&b);
	outTree->Branch("ch",&ch);
	outTree->Branch("l",&l);
	outTree->Branch("w",&w);
	outTree->Branch("h",&h);
	N = chain_preWire->GetEntries();
	for (int i = 0; i<N; i++){
		chain_preWire->GetEntry(i);
		if (l<=3){
			if (w==3){
				xhv+=offset[l];
			}
		}
		else{
			if (w==4){
				xhv+=offset[l];
			}
		}
		std::cout<<l<<":"<<offset[l]<<std::endl;
		outTree->Fill();
	}
	outTree->Write();
	outFile->Close();

	return 0;
}
