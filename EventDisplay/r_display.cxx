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
	TFile * TFile_wirepos = new TFile("../info/wire-position.root");
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
			map_lid[wp_bd][wp_ch] = wp_lid;
		}
	}

	//TFile * ifile = new TFile("run.002.try003.root");
	TFile * ifile = new TFile("../root/run_000117_built.root");
	TTree * it = (TTree*) ifile->Get("tree");

	int triggerNumber;
	int tdcNhit[96];
	int driftTime[96][32];

	it->SetBranchAddress("triggerNumber",&triggerNumber);
	it->SetBranchAddress("tdcNhit",tdcNhit);
	it->SetBranchAddress("driftTime",driftTime);

	TEllipse * ewiret[96];
	double wx,wy,wz,dd,ddt;
	int chs,bd,lid,wid;
	std::stringstream buf;
	TCanvas * c = new TCanvas("c","c",896,896);
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
//	for ( int i = 0 ; i<it->GetEntries(); i++){
	for ( int i = 0 ; i<100; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/it->GetEntries()*100);
		it->GetEntry(i);
//		if (i!=10) continue; // 0
		buf.str("");
		buf.clear();
		//FIXME
		buf.str("");
		buf.clear();
		buf<<"Entry#"<<i<<", TriggerNumber#"<<triggerNumber;
		h0->SetTitle(buf.str().c_str());
		for (int ch = 0; ch<96; ch++){
			if (ewiret[ch]){
				delete ewiret[ch];
				ewiret[ch] = 0;
			}
			if (text[ch]){
				delete text[ch];
				text[ch] = 0;
			}
			if (tdcNhit[ch]<=0) continue;
			chs = ch%48;
			bd = ch/48;
			lid = map_lid[bd][chs];
			wid = map_wid[bd][chs];
			if (lid>=7) continue;
			if (wid>=11) continue;
			wx = map_xc[lid][wid];
			wy = map_yc[lid][wid];
			dd = (driftTime[ch][0]+840)/0.96*0.8/200;
			ewiret[ch] = new TEllipse(wx,wy,dd,dd);
			ewiret[ch]->SetFillStyle(0);
			ewiret[ch]->SetLineColor(kRed);
			ewiret[ch]->Draw("SAME"); // Draw hits.
			text[ch]= new TText(wx,wy,Form("%d,%d;%d",lid,wid,ch));
			text[ch]->SetTextSize(0.02);
			text[ch]->Draw("SAME");
		}
		buf.str("");
		buf.clear();
		//FIXME
//		buf<<i<<"_before.pdf";
		buf<<i<<".pdf";
//		buf<<i<<"_after.pdf";
		c->SaveAs(buf.str().c_str());
//		c->WaitPrimitive();
	//	c->Update();
//		while(1){}
	}
	return 0;
}
