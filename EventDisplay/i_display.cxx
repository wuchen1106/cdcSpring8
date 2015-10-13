#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraph.h"
#include <sstream>
#include <iostream>
#include <stdlib.h>

int main(int argc, char** argv){
	if (argc<2) return -1;
	int runNo = (int)strtol(argv[1],NULL,10);
	std::string suffix = "";
	if (argc>=3){
		suffix  = argv[2];
		suffix="."+suffix;
	}
	int nEventMax = 0;
	if (argc>=4) nEventMax = (int)strtol(argv[3],NULL,10);

	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bid;
	int wp_ch;
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
	int map_ch[8][11];
	TTree_wirepos->SetBranchAddress("b",&wp_bid);
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
	TTree_wirepos->SetBranchAddress("xro",&wp_xro);
	TTree_wirepos->SetBranchAddress("yro",&wp_yro);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<=7){
			map_xro[wp_lid][wp_wid] = wp_xro;
			map_yro[wp_lid][wp_wid] = wp_yro;
			map_xhv[wp_lid][wp_wid] = wp_xhv;
			map_yhv[wp_lid][wp_wid] = wp_yhv;
			map_ch[wp_lid][wp_wid] = wp_ch+wp_bid*48;
		}
	}
	double yup = 623.97007;
	double ydown = 527.60011;

	//===================Get ROOT File============================
	//TChain * c = new TChain("t","t");
	TChain * c = new TChain("t","t");
	std::stringstream buf;
	buf.str(""); buf.clear();
	buf<<"../root/i_"<<runNo<<suffix<<".root";
	c->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	int triggerNumber;
	double xup;
	double xdown;
	std::vector<double> * i_driftD = 0;
	std::vector<int> * i_lr = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_layerID = 0;
	c->SetBranchAddress("lr",&i_lr);
	c->SetBranchAddress("wireID",&i_wireID);
	c->SetBranchAddress("layerID",&i_layerID);
	c->SetBranchAddress("driftD",&i_driftD);
	c->SetBranchAddress("tx1",&xup);
	c->SetBranchAddress("tx2",&xdown);
	c->SetBranchAddress("triggerNumber",&triggerNumber);

	//==================Get ADC==========================
	TChain * chain2 = new TChain("tree","tree");
	chain2->Add(Form("../root/run_%0.6d_built.root",runNo));
	int adc[99][32];
	int tdc[99][32];
	int clockNumberDriftTime[99][32];
	int tdcNhit[99];
	int triggerNumber2;
	chain2->SetBranchAddress("adc",adc);
	chain2->SetBranchAddress("driftTime",tdc);
	chain2->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	chain2->SetBranchAddress("tdcNhit",tdcNhit);
	chain2->SetBranchAddress("triggerNumber",&triggerNumber2);

	TGraph *gr_waveForm = 0;
	Int_t vSample[32];
	for (Int_t i=0; i<32; i++){
		vSample[i] = i;
	}
	TLatex *textTDC[32];
	TMarker *markerTDC[32];
	for (Int_t j=0; j<32; j++) {
		textTDC[j] = new TLatex(0,235,"");
		textTDC[j]->SetTextSize(0.04);
		textTDC[j]->SetTextColor(kRed);
		markerTDC[j] = new TMarker(0,0,20);
		markerTDC[j]->SetMarkerSize(0.55);
	}
	TEllipse * ewiret[8][11];
	double wx,wy,wz,dd,ddt;
	double wxro,wyro,wzro;
	double wxhv,wyhv,wzhv;
	wzro = 599.17/2;
	wzhv = -599.17/2;
	int lid,wid;
	TCanvas * ca = new TCanvas("ca","ca",900,900);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	TH2D * h0 = new TH2D("h0","h0",128,-130,130,128,500,650);
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

//	TTree_wirepos->SetMarkerStyle(20);
//	TTree_wirepos->SetMarkerSize(0.5);
//	TTree_wirepos->Draw("yc/10.:xc/10.","","SAME");
	std::string suf = ".png";
//	for ( int i = 0 ; i<c->GetEntries(); i++){
	for ( int i = 0; i<100; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/c->GetEntries()*100);
		c->GetEntry(i);

		buf.str("");
		buf.clear();
		//FIXME
		buf.str("");
		buf.clear();
		buf<<"Entry#"<<i<<", TriggerNumber#"<<triggerNumber;
		h0->SetTitle(buf.str().c_str());
		h0->Draw();
		// FIXME
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
			wxro = map_xro[lid][wid];
			wyro = map_yro[lid][wid];
			wxhv = map_xhv[lid][wid];
			wyhv = map_yhv[lid][wid];
			wy = (wyro+wyhv)/2.;
			wx = (wxro+wxhv)/2.;
			dd = (*i_driftD)[ihit];
			ewiret[lid][wid] = new TEllipse(wx,wy,dd,dd);
			ewiret[lid][wid]->SetFillStyle(0);
			if ((*i_lr)[ihit])
				ewiret[lid][wid]->SetLineColor(kRed);
			else
				ewiret[lid][wid]->SetLineColor(kBlue);
			ewiret[lid][wid]->Draw("SAME"); // Draw hits.
			text[lid][wid]= new TText(wx,wy,Form("%d,%d",lid,wid));
			text[lid][wid]->SetTextSize(0.02);
			text[lid][wid]->Draw("SAME");
		}
		l->SetX1(xup);
		l->SetX2(xdown);
		l->Draw("SAME");

//		for (int index2 = triggerNumber; ;index2--){
//			chain2->GetEntry(index2);
//			if (triggerNumber==triggerNumber2) break;
//		}

		buf.str("");
		buf.clear();
		//FIXME
//		buf<<i<<"_before.pdf";
		buf<<i<<suf;
//		buf<<i<<"_after.pdf";
		ca->SaveAs(buf.str().c_str());
//		ca->WaitPrimitive();
	//	ca->Update();
//		while(1){}
	}
	return 0;
}
