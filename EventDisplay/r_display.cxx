#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include <sstream>
#include <iostream>

int main(int argc, char** argv){
	int hmin = 13;
	int t0 = -849;
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
	double map_xc[9][11];
	double map_yc[9][11];
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("b",&wp_bd);
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xc",&wp_xc);
	TTree_wirepos->SetBranchAddress("yc",&wp_yc);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		map_wid[wp_bd][wp_ch] = wp_wid;
		map_lid[wp_bd][wp_ch] = wp_lid;
		map_xc[wp_lid][wp_wid] = wp_xc/10.;
		map_yc[wp_lid][wp_wid] = wp_yc/10.;
	}

	//TFile * ifile = new TFile("run.002.try003.root");
	gROOT->ProcessLine(".L loader.C+");
	TFile * ifile = new TFile("../root/p_117.root");
	TTree * it = (TTree*) ifile->Get("t");

	int triggerNumber;
	int i_nHits[50];
	double i_p[96];
	double i_aa[96];
	int i_n[96][50];
	std::vector<std::vector<int> > * i_tdc = 0;
	std::vector<std::vector<int> > * i_peak = 0;
	std::vector<std::vector<int> > * i_clk = 0;
	std::vector<std::vector<int> > * i_width = 0;
	std::vector<std::vector<double> > * i_sum = 0;

	it->SetBranchAddress("triggerNumber",&triggerNumber);
	it->SetBranchAddress("nh",i_nHits);
	it->SetBranchAddress("np",i_n);
	it->SetBranchAddress("ped",i_p);
	it->SetBranchAddress("aa",i_aa);
	it->SetBranchAddress("clk",&i_clk);
	it->SetBranchAddress("peak",&i_peak);
	it->SetBranchAddress("width",&i_width);
	it->SetBranchAddress("sum",&i_sum);
	it->SetBranchAddress("tdc",&i_tdc);

	TFile * iifile = new TFile("../root/run_000117_built.root");
	TTree * iit = (TTree*) iifile->Get("tree");
	iit->SetMarkerStyle(20);

	TString prefix = "";
	TEllipse * ewiret[96];
	double wx,wy,wz,dd,ddt;
	int chs,bd,lid,wid;
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
	TPad * pad[96];
	TCanvas * c2 = new TCanvas("c2","c2",896,896);
	c2->cd();
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}
	TCanvas * c3 = new TCanvas("c3","c3",896,896);
	c3->cd();
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i+48;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}

	TTree_wirepos->SetMarkerStyle(20);
	TTree_wirepos->SetMarkerSize(0.5);
	for ( int i = 0 ; i<it->GetEntries(); i++){
//	for ( int i = 0 ; i<100; i++){
//		if (i!=195
//		  &&i!=205
//		  )
//		continue;
		if (i%100==0) printf("%lf%...\n",(double)i/it->GetEntries()*100);
		it->GetEntry(i);
		if (i_nHits[hmin]!=10&&i_nHits[hmin]!=11) continue;
//		if (i!=10) continue; // 0
		h0->SetTitle(Form("Entry#%d, TriggerNumber#%d, %d hits",i,triggerNumber,i_nHits[hmin]));
		c->cd();
		h0->Draw();
		TTree_wirepos->Draw("yc/10.:xc/10.","","SAME");
		for (int ch = 0; ch<96; ch++){
			pad[ch]->cd();
			iit->Draw(Form("adc[%d]:Iteration$>>hist%d(32,0,32,400,180,580)",ch,ch),"","LP",1,i);
//			it->Draw(Form("adc[%d]:Iteration$",ch),"","LP",1,i);
			c->cd();
			if (ewiret[ch]){
				delete ewiret[ch];
				ewiret[ch] = 0;
			}
			if (text[ch]){
				delete text[ch];
				text[ch] = 0;
			}
			if (i_n[ch][hmin]<=0) continue;
			chs = ch%48;
			bd = ch/48;
			lid = map_lid[bd][chs];
			wid = map_wid[bd][chs];
			wx = map_xc[lid][wid];
			wy = map_yc[lid][wid];
			for (int ipeak = 0; ipeak<i_n[ch][0]; ipeak++){
				if ((*i_peak)[ch][ipeak]-i_p[ch]>=hmin){
					dd = ((*i_tdc)[ch][0]-t0)/0.96*0.8/200;
					break;
				}
			}
			ewiret[ch] = new TEllipse(wx,wy,dd,dd);
			ewiret[ch]->SetFillStyle(0);
			ewiret[ch]->SetLineColor(kRed);
			ewiret[ch]->Draw("SAME"); // Draw hits.
			text[ch]= new TText(wx,wy,Form("%d,%d;%d",lid,wid,ch));
			text[ch]->SetTextSize(0.02);
			text[ch]->Draw("SAME");
		}
		if (i_nHits[hmin]==10) prefix = "n10.";
		if (i_nHits[hmin]==11) prefix = "n11.";
//		else if (i_nHits[hmin]>9) prefix = "multi.";
//		else if (i_nHits[hmin]>7) prefix = "single.";
//		else prefix = "incom.";
		c->SaveAs(prefix+Form("ep.%d.pdf",i));
		c2->SaveAs(prefix+Form("wf.0.%d.pdf",i));
		c3->SaveAs(prefix+Form("wf.1.%d.pdf",i));
//		c->WaitPrimitive();
	//	c->Update();
//		while(1){}
	}
	return 0;
}
