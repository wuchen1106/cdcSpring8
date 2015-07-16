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
#include <stdlib.h>

int main(int argc, char** argv){
	if (argc<2) return -1;
	int runNo = (int)strtol(argv[1],NULL,10);

	double dtmax = 550;
	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.v3.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bd;
	int wp_ch;
	int wp_wid;
	int wp_lid;
	double wp_xc;
	double wp_yc;
	int map_lid[2][48];
	int map_wid[2][48];
	double map_xc[8][11];
	double map_yc[8][11];
	for (int i = 0; i<2; i++){
		for ( int j = 0; j<48; j++){
			map_lid[i][j] = -1;
			map_wid[i][j] = -1;
		}
	}
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("b",&wp_bd);
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xc",&wp_xc);
	TTree_wirepos->SetBranchAddress("yc",&wp_yc);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<=7){
			map_xc[wp_lid][wp_wid] = wp_xc/10.;
			map_yc[wp_lid][wp_wid] = wp_yc/10.;
			map_wid[wp_bd][wp_ch] = wp_wid;
			map_lid[wp_bd][wp_ch] = wp_lid;
		}
	}

	TFile * ifile = new TFile(Form("../root/d_%d.root",runNo));
	TTree * it = (TTree*) ifile->Get("t");

	int triggerNumber;
	int iHit[96];
	double driftTime[96];

	it->SetBranchAddress("triggerNumber",&triggerNumber);
	it->SetBranchAddress("i",&iHit);
	it->SetBranchAddress("dt",driftTime);

	TEllipse * ewiret[96];
	double wx,wy,wz,dd,ddt;
	int chs,bd,lid,wid;
	std::stringstream buf;
	TCanvas * c = new TCanvas("c","c",896,896);
	c->SetGridx(1);
	c->SetGridy(1);
	TH2D * h0 = new TH2D("h0","h0",128,-13,13,128,50,65);
	gStyle->SetOptStat(0);
	TLine * l = new TLine();
	l->SetLineColor(kRed);
	TText * text[96];
	for (int ch = 0; ch<96; ch++){
		text[ch] = 0;
		ewiret[ch] = 0;
	}

	h0->Draw();
	TTree_wirepos->SetMarkerStyle(20);
	TTree_wirepos->SetMarkerSize(0.5);
	TTree_wirepos->Draw("yc/10.:xc/10.","","SAME");
	std::string suffix = ".pdf";
//	for ( int i = 0 ; i<it->GetEntries(); i++){
	for ( int i = 0 ; i<100; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/it->GetEntries()*100);
		it->GetEntry(i);
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
		for (int ch = 0; ch<96; ch++){
			bd = ch/48;
			chs = ch%48;
			lid = map_lid[bd][chs];
			if (lid<0) continue;
			if (iHit[ch]>=0&&driftTime[ch]<dtmax) nHits++;
		}
		if (nHits>9){
			suffix = ".l9.pdf";
		}
		else{
			suffix = ".s9.pdf";
		}
		for (int ch = 0; ch<96; ch++){
			if (ewiret[ch]){
				delete ewiret[ch];
				ewiret[ch] = 0;
			}
			if (text[ch]){
				delete text[ch];
				text[ch] = 0;
			}
			if (iHit[ch]<0||driftTime[ch]>dtmax) continue;
			chs = ch%48;
			bd = ch/48;
			lid = map_lid[bd][chs];
			wid = map_wid[bd][chs];
			if (lid<0) continue;
			wx = map_xc[lid][wid];
			wy = map_yc[lid][wid];
			dd = driftTime[ch]*0.8/dtmax;
			ewiret[ch] = new TEllipse(wx,wy,dd,dd);
			ewiret[ch]->SetFillStyle(0);
			ewiret[ch]->SetLineColor(kRed);
			ewiret[ch]->Draw("SAME"); // Draw hits.
			text[ch]= new TText(wx,wy,Form("%d,%d;%d",lid,wid,ch%48));
			text[ch]->SetTextSize(0.02);
			text[ch]->Draw("SAME");
		}
		buf.str("");
		buf.clear();
		//FIXME
//		buf<<i<<"_before.pdf";
		buf<<i<<suffix;
//		buf<<i<<"_after.pdf";
		c->SaveAs(buf.str().c_str());
//		c->WaitPrimitive();
	//	c->Update();
//		while(1){}
	}
	return 0;
}
