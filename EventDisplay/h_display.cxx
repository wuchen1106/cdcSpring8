#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <sstream>
#include <iostream>
#include <stdlib.h>

int main(int argc, char** argv){
	if (argc<2) return -1;
	int runNo = (int)strtol(argv[1],NULL,10);

	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.v3.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_wid;
	int wp_lid;
	double wp_xro;
	double wp_yro;
	double wp_xhv;
	double wp_yhv;
	double map_xro[8][11];
	double map_yro[8][11];
	double map_xhv[8][11];
	double map_yhv[8][11];
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
	TTree_wirepos->SetBranchAddress("xro",&wp_xro);
	TTree_wirepos->SetBranchAddress("yro",&wp_yro);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<=7){
			map_xro[wp_lid][wp_wid] = wp_xro/10.;
			map_yro[wp_lid][wp_wid] = wp_yro/10.;
			map_xhv[wp_lid][wp_wid] = wp_xhv/10.;
			map_yhv[wp_lid][wp_wid] = wp_yhv/10.;
		}
	}
	double yup = 62.397007;
	double ydown = 52.760011;

	//===================Get ROOT File============================
	//TChain * c = new TChain("t","t");
	TChain * c = new TChain("t","t");
	std::stringstream buf;
	buf.str(""); buf.clear();
	buf<<"../root/i_"<<runNo<<".root";
	c->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	Long64_t triggerNumber;
	double xup, xdown;
	double zup, zdown;
	std::vector<double> * i_driftD = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_layerID = 0;
	c->SetBranchAddress("wireID",&i_wireID);
	c->SetBranchAddress("layerID",&i_layerID);
	c->SetBranchAddress("driftD",&i_driftD);
	c->SetBranchAddress("tx1",&xup);
	c->SetBranchAddress("tz1",&zup);
	c->SetBranchAddress("tx2",&xdown);
	c->SetBranchAddress("tz2",&zdown);
	c->SetBranchAddress("triggerNumber",&triggerNumber);

	TEllipse * ewiret[8][11];
	double wx,wy,wz,dd,ddt;
	double wxro,wyro,wzro;
	double wxhv,wyhv,wzhv;
	wzro = 59.917/2;
	wzhv = -59.917/2;
	int lid,wid;
	TCanvas * ca = new TCanvas("ca","ca",896,896);
	TH2D * h0 = new TH2D("h0","h0",128,-13,13,128,50,65);
	gStyle->SetOptStat(0);
	TLine * l = new TLine();
	l->SetLineColor(kRed);
	l->SetY1(yup);
	l->SetY2(ydown);
	TText * text[8][11];
	for (int lid = 0; lid<8; lid++){
		for (int wid = 0; wid<11; wid++){
			text[lid][wid] = 0;
			ewiret[lid][wid] = 0;
		}
	}

	h0->Draw();
//	TTree_wirepos->SetMarkerStyle(20);
//	TTree_wirepos->SetMarkerSize(0.5);
//	TTree_wirepos->Draw("yc/10.:xc/10.","","SAME");
	std::string suffix = ".pdf";
//	for ( int i = 0 ; i<c->GetEntries(); i++){
	for ( int i = 0 ; i<100; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/c->GetEntries()*100);
		c->GetEntry(i);
		buf.str("");
		buf.clear();
		//FIXME
		buf.str("");
		buf.clear();
		buf<<"Entry#"<<i<<", TriggerNumber#"<<triggerNumber;
		h0->SetTitle(buf.str().c_str());
		// FIXME
//		if (iHit[83]<0) continue;
		int nHits = 0;
		for (int lid = 0; lid<8; lid++){
			for (int wid = 0; wid<11; wid++){
				if (ewiret[lid][wid]){
					delete ewiret[lid][wid];
					ewiret[lid][wid] = 0;
				}
				if (text[lid][wid]){
					delete text[lid][wid];
					text[lid][wid] = 0;
				}
			}
		}
		for (int ihit = 0; ihit<i_driftD->size(); ihit++){
			lid = (*i_layerID)[ihit];
			wid = (*i_wireID)[ihit];
			if (lid<0) continue;
			wxro = map_xro[lid][wid];
			wyro = map_yro[lid][wid];
			wxhv = map_xhv[lid][wid];
			wyhv = map_yhv[lid][wid];
			wy = (wyro+wyhv)/2.;
			wz = ((yup-wy)*zdown+(wy-ydown)*zup)/(yup-ydown);
			wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
			wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
			dd = (*i_driftD)[ihit];
			ewiret[lid][wid] = new TEllipse(wx,wy,dd,dd);
			ewiret[lid][wid]->SetFillStyle(0);
			ewiret[lid][wid]->SetLineColor(kRed);
			ewiret[lid][wid]->Draw("SAME"); // Draw hits.
			text[lid][wid]= new TText(wx,wy,Form("%d,%d",lid,wid));
			text[lid][wid]->SetTextSize(0.02);
			text[lid][wid]->Draw("SAME");
		}
		l->SetX1(xup);
		l->SetX2(xdown);
		l->Draw("SAME");
		buf.str("");
		buf.clear();
		//FIXME
//		buf<<i<<"_before.pdf";
		buf<<i<<suffix;
//		buf<<i<<"_after.pdf";
		ca->SaveAs(buf.str().c_str());
//		ca->WaitPrimitive();
	//	ca->Update();
//		while(1){}
	}
	return 0;
}
