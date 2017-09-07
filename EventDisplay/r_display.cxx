#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMarker.h"
#include "TLatex.h"
#include <stdlib.h>
#include <sstream>
#include <iostream>

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}

int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int hmin = 13;
	int smin = 10.7;
	int aamin = 40;
	int t0 = -849;
	int tmax = -635;
	int dtmax = tmax-t0;
	int tres = 2;
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
	TFile * ifile = new TFile(Form("../root/p_%d.root",runNo));
	TTree * it = (TTree*) ifile->Get("t");

	int triggerNumber;
	int i_nh[50];
	double i_p[96];
	double i_aa[96];
	int i_np[96][50];
	std::vector<std::vector<int> > * i_tdc = 0;
	std::vector<std::vector<int> > * i_peak = 0;
	std::vector<std::vector<int> > * i_clk = 0;
	std::vector<std::vector<int> > * i_width = 0;
	std::vector<std::vector<double> > * i_sum = 0;

	it->SetBranchAddress("triggerNumber",&triggerNumber);
	it->SetBranchAddress("nh",i_nh);
	it->SetBranchAddress("np",i_np);
	it->SetBranchAddress("ped",i_p);
	it->SetBranchAddress("aa",i_aa);
	it->SetBranchAddress("clk",&i_clk);
	it->SetBranchAddress("peak",&i_peak);
	it->SetBranchAddress("width",&i_width);
	it->SetBranchAddress("sum",&i_sum);
	it->SetBranchAddress("tdc",&i_tdc);

	TFile * iifile = new TFile(Form("../root/run_%0.6d_built.root",runNo));
	TTree * iit = (TTree*) iifile->Get("tree");
	iit->SetMarkerStyle(20);
	iit->SetMarkerSize(0.3);
	int i_adc[96][32];
	iit->SetBranchAddress("adc",i_adc);

	TMarker * marker[96][32];
	TLatex * text2[96][32];
	for ( int ch = 0; ch<96; ch++){
		for ( int i = 0; i<32; i++){
			marker[ch][i] = new TMarker(1,1,20);
			marker[ch][i]->SetMarkerColor(kRed);
			marker[ch][i]->SetMarkerSize(0.4);
			text2[ch][i] = new TLatex(1,1,"");
			text2[ch][i]->SetTextColor(kRed);
		}
	}

	TString prefix = "";
	TEllipse * ewiret[96];
	TEllipse * ewiret2[96];
	TEllipse * ewiret3[96];
	double wx,wy,wz,dd,ddt;
	int chs,bd,lid,wid;
	TCanvas * c = new TCanvas("c","c",896,896);
	TLatex * canvtitle = new TLatex(0.1,0.98,"");
	canvtitle->SetTextSize(0.02);
	TH2D * h0 = new TH2D("h0","h0",128,-13,13,128,50,65);
	gStyle->SetOptStat(0);
	TLine * l = new TLine();
	l->SetLineColor(kRed);
	TText * text[96];
	TLatex* text3[96];
	for (int ch = 0; ch<96; ch++){
		text[ch] = 0;
		text3[ch]= new TLatex(0.1,0.8,"");
		text3[ch]->SetTextColor(kRed);
		ewiret[ch] = 0;
		ewiret2[ch] = 0;
		ewiret3[ch] = 0;
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
//	for ( int i = 0 ; i<it->GetEntries(); i++){
	for ( int i = 0 ; i<100; i++){
//		if (i!=195
//		  &&i!=205
//		  )
//		continue;
		if (i%100==0) printf("%lf%...\n",(double)i/it->GetEntries()*100);
		it->GetEntry(i);
		prefix = "";

		// ************Cutting on channel*****************
//		if (i_nh[hmin]!=10&&i_nh[hmin]!=11) continue;
//		if (i!=10) continue; // 0
		if (i_np[2][0]!=0) prefix = "n2";
//		if (i_nh[hmin]<=7) prefix = "incom.";
//		else if (i_nh[hmin]==8) prefix = "single.";
//		else if (i_nh[hmin]==9) prefix = "single.";
//		else if (i_nh[hmin]==10) prefix = "n10.";
//		else if (i_nh[hmin]==11) prefix = "n11.";
//		else if (i_nh[hmin]==12) prefix = "n12.";
//		else if (i_nh[hmin]>=13) prefix = "multi.";

		// ************Cutting on peaks in ch15****************
//		int thepeak = 0;
//		for (; thepeak<(*i_peak)[15].size(); thepeak++){
//			if ((*i_sum)[15][thepeak]>=smin&&(*i_tdc)[15][thepeak]<=tmax+tres&&(*i_tdc)[15][thepeak]>=t0-tres) break;
//		}
//		if (thepeak==(*i_peak)[15].size()){
//			if (i_aa[15]>aamin) prefix = "missed.";
//		}
//		else{
//			if (thepeak<(*i_peak)[15].size()-1&&(*i_sum)[15][thepeak]<(*i_sum)[15][thepeak+1])
//				prefix = "shadowed.";
//			if (thepeak>0&&(*i_tdc)[15][thepeak-1]<=tmax+tres&&(*i_tdc)[15][thepeak-1]>=t0-tres)
//				prefix = "prepeak.";
//		}

//		if (prefix=="") continue;

		h0->SetTitle(Form("Entry#%d, TriggerNumber#%d, %d hits",i,triggerNumber,i_nh[hmin]));
		c->cd();
		h0->Draw();
		TTree_wirepos->Draw("yc/10.:xc/10.","","SAME");
		canvtitle->SetText(0.1,0.98,Form("Entry#%d, TriggerNumber#%d, %d hits",i,triggerNumber,i_nh[hmin]));
		for (int ch = 0; ch<96; ch++){
			if (ch==0) {
				c2->cd();
				canvtitle->Draw("SAME");
			}
			else if (ch==48) {
				c3->cd();
				canvtitle->Draw("SAME");
			}
			pad[ch]->cd();
			iit->Draw(Form("adc[%d]:Iteration$>>hist%d(32,0,32,700,50,750)",ch,ch),"","LP",1,i);
			for (int ipeak = 0; ipeak<(*i_peak)[ch].size(); ipeak++){
				marker[ch][ipeak]->SetX((*i_clk)[ch][ipeak]);
				marker[ch][ipeak]->SetY(i_adc[ch][(*i_clk)[ch][ipeak]]);
				text2[ch][ipeak]->SetText((*i_clk)[ch][ipeak],i_adc[ch][(*i_clk)[ch][ipeak]]+(0.5-ipeak%2)*100,Form("%d",(*i_tdc)[ch][ipeak]));
				marker[ch][ipeak]->Draw("SAME");
				text2[ch][ipeak]->Draw("SAME");
			}
			text3[ch]->SetText(2,700,Form("%.1lf",i_aa[ch]));
			text3[ch]->Draw("SAME");
//			it->Draw(Form("adc[%d]:Iteration$",ch),"","LP",1,i);
			c->cd();
			if (ewiret[ch]){
				delete ewiret[ch];
				ewiret[ch] = 0;
			}
			if (ewiret3[ch]){
				delete ewiret3[ch];
				ewiret3[ch] = 0;
			}
			if (ewiret2[ch]){
				delete ewiret2[ch];
				ewiret2[ch] = 0;
			}
			if (text[ch]){
				delete text[ch];
				text[ch] = 0;
			}
			chs = ch%48;
			bd = ch/48;
			lid = map_lid[bd][chs];
			wid = map_wid[bd][chs];
			wx = map_xc[lid][wid];
			wy = map_yc[lid][wid];
			int ipeak = 0;
			for (; ipeak<i_np[ch][0]; ipeak++){
				if ((*i_sum)[ch][ipeak]>=smin){
					break;
				}
			}
			if (ipeak<i_np[ch][0]){
				dd = ((*i_tdc)[ch][ipeak]-t0)/0.96*0.8/dtmax;
				ewiret[ch] = new TEllipse(wx,wy,dd,dd);
				ewiret[ch]->SetFillStyle(0);
				ewiret[ch]->SetLineColor(kRed);
				ewiret[ch]->Draw("SAME"); // Draw hits.
				text[ch]= new TText(wx,wy,Form("%d,%d;%d",lid,wid,ch));
				text[ch]->SetTextSize(0.02);
				text[ch]->Draw("SAME");
			}
			if (ipeak>0){
				dd = ((*i_tdc)[ch][0]-t0)/0.96*0.8/dtmax;
				ewiret2[ch] = new TEllipse(wx,wy,dd,dd);
				ewiret2[ch]->SetFillStyle(0);
				ewiret2[ch]->SetLineColor(kBlue);
				ewiret2[ch]->Draw("SAME"); // Draw hits.
			}
			int jpeak = ipeak+1;
			for (; jpeak<i_np[ch][0]; jpeak++){
				if ((*i_sum)[ch][jpeak]>(*i_sum)[ch][ipeak]){
					dd = ((*i_tdc)[ch][jpeak]-t0)/0.96*0.8/dtmax;
					break;
				}
			}
			if (jpeak<i_np[ch][0]){
				dd = ((*i_tdc)[ch][jpeak]-t0)/0.96*0.8/dtmax;
				ewiret3[ch] = new TEllipse(wx,wy,dd,dd);
				ewiret3[ch]->SetFillStyle(0);
				ewiret3[ch]->SetLineColor(kOrange);
				ewiret3[ch]->Draw("SAME"); // Draw hits.
			}
		}
		c->SaveAs(prefix+Form("ep.%d.pdf",i));
		c2->SaveAs(prefix+Form("wf.0.%d.pdf",i));
		c3->SaveAs(prefix+Form("wf.1.%d.pdf",i));
//		c->WaitPrimitive();
	//	c->Update();
//		while(1){}
	}
	return 0;
}
