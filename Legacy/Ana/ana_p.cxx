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
#include <sstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv){
	int tdcmin = -846;
	int ch = 15;
	int hmin = 13;

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

	TH1D * h_h_1p = new TH1D("h_1p","adc peak",700,50,750);
	TH1D * h_h_2p = new TH1D("h_2p","adc peak",700,50,750);

	TH1D * h_h_w1 = new TH1D("h_w1","adc peak",700,50,750);
	TH1D * h_h_w2 = new TH1D("h_w2","adc peak",700,50,750);
	TH1D * h_h_w3 = new TH1D("h_w3","adc peak",700,50,750);
	TH1D * h_h_w4 = new TH1D("h_w4","adc peak",700,50,750);

	TH2D * h_wf_1p = new TH2D("wf_1p","wf",42,-10,32,700,50,750);
	TH2D * h_wf_2p = new TH2D("wf_2p","wf",42,-10,32,700,50,750);

	TH2D * h_wf_w1 = new TH2D("wf_w1","wf",42,-10,32,700,50,750);
	TH2D * h_wf_w2 = new TH2D("wf_w2","wf",42,-10,32,700,50,750);
	TH2D * h_wf_w3 = new TH2D("wf_w3","wf",42,-10,32,700,50,750);
	TH2D * h_wf_w4 = new TH2D("wf_w4","wf",42,-10,32,700,50,750);

	TCanvas * c = new TCanvas();
	TMarker * marker[32];
	TLatex * text[32];
	for ( int i = 0; i<32; i++){
		marker[i] = new TMarker(1,1,20);
		marker[i]->SetMarkerColor(kRed);
		text[i] = new TLatex(1,1,"");
		text[i]->SetTextColor(kBlack);
	}

	TFile * iifile = new TFile("../root/run_000117_built.root");
	TTree * iit = (TTree*) iifile->Get("tree");
	int i_adc[96][32];
	iit->SetBranchAddress("adc",i_adc);

	TFile *ofile = new TFile("output.root","RECREATE");
	TTree * otree = new TTree("t","t");

int Nevents = it->GetEntries();
	for ( int i = 0 ; i<Nevents; i++){
//	for ( int i = 0 ; i<10000; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/Nevents*100);
		it->GetEntry(i);

		if (i_nHits[hmin]<8||i_nHits[hmin]>9) continue;

		double offset = 0;
		int npeaks = 0;
		int thepeak = -1;
		for ( int ipeak = 0; ipeak<i_n[ch][0]; ipeak++){
			if ((*i_tdc)[ch][ipeak]<tdcmin) continue;
			if ((*i_peak)[ch][ipeak]-i_p[ch]<hmin) continue;
			npeaks++;
			if (thepeak==-1) thepeak = ipeak;
		}

		if (npeaks<=0) continue;
		offset = ((*i_tdc)[ch][thepeak]+1000)/1000.*32;
		int peak = (*i_peak)[ch][thepeak];
		int width = (*i_width)[ch][thepeak];

		TString suffix = ".";

		if (npeaks==1){
			h_h_1p->Fill(peak);
			if (width==1){
				h_h_w1->Fill(peak);
				suffix = ".1p.w1.";
			}
			else if (width==2){
				h_h_w2->Fill(peak);
				suffix = ".1p.w2.";
			}
			else if (width==3){
				h_h_w3->Fill(peak);
				suffix = ".1p.w3.";
			}
			else if (width==4){
				h_h_w4->Fill(peak);
				suffix = ".1p.w4.";
			}
		}
		else if (npeaks==2){
			h_h_2p->Fill(peak);
			suffix = ".2p.";
		}

		iit->GetEntry(i);
		for ( int iadc = 0; iadc<32; iadc++ ){
			if (npeaks==1){
				h_wf_1p->Fill(iadc-offset,i_adc[ch][iadc]);
				if (width==1){
					h_wf_w1->Fill(iadc-offset,i_adc[ch][iadc]);
				}
				else if (width==2){
					h_wf_w2->Fill(iadc-offset,i_adc[ch][iadc]);
				}
				else if (width==3){
					h_wf_w3->Fill(iadc-offset,i_adc[ch][iadc]);
				}
				else if (width==4){
					h_wf_w4->Fill(iadc-offset,i_adc[ch][iadc]);
				}
			}
			else if (npeaks==2){
				h_wf_2p->Fill(iadc-offset,i_adc[ch][iadc]);
			}
		}


		if (i<500){
			iit->Draw(Form("adc[%d]:Iteration$",ch),"","LP",1,i);
			for (int ipeak = 0; ipeak<(*i_peak)[ch].size(); ipeak++){
				marker[ipeak]->SetX((*i_clk)[ch][ipeak]);
				marker[ipeak]->SetY(i_adc[ch][(*i_clk)[ch][ipeak]]);
				text[ipeak]->SetText((*i_clk)[ch][ipeak]+0.5,i_adc[ch][(*i_clk)[ch][ipeak]],Form("%d",(*i_tdc)[ch][ipeak]));
				if (ipeak == thepeak){
					marker[ipeak]->SetMarkerColor(kRed);
					text[ipeak]->SetTextColor(kRed);
				}
				else{
					marker[ipeak]->SetMarkerColor(kBlack);
					text[ipeak]->SetTextColor(kBlack);
				}
				marker[ipeak]->Draw("SAME");
				text[ipeak]->Draw("SAME");
			}
			c->SaveAs(Form("%d",i)+suffix+"pdf");
		}

		otree->Fill();
	}
	h_wf_1p->Write();
	h_wf_2p->Write();
	h_wf_w1->Write();
	h_wf_w2->Write();
	h_wf_w3->Write();
	h_wf_w4->Write();
	h_h_1p->Write();
	h_h_2p->Write();
	h_h_w1->Write();
	h_h_w2->Write();
	h_h_w3->Write();
	h_h_w4->Write();
	otree->Write();
	ofile->Close();

	return 0;
}
