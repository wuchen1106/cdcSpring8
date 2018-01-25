#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLatex.h"

#include "header.h"

double tscale = 0.96;

void print_usage(char* prog_name);
int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	if (argc>=3) nEventMax = (int)strtol(argv[2],NULL,10);
	TString suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix=suffix+".";
	}

    TString HOME=getenv("CDCS8WORKING_DIR");

	//===================Get channel map============================
	TFile * TFile_wirepos = new TFile(HOME+"/info/wire-position.v3.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bid;
	int wp_wid;
	int wp_lid;
	int wp_ch;
	int map_wid[NCHT];
	int map_lid[NCHT];
	int widmax[NLAY];
	for ( int il = 0; il<NLAY; il++ ){
		widmax[il] = 0;
	}
	for ( int ch = 0; ch<NCHT; ch++ ){
		map_lid[ch] = -1;
		map_wid[ch] = -1;
	}
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("b",&wp_bid);
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		map_wid[wp_ch+wp_bid*NCHS] = wp_wid;
		map_lid[wp_ch+wp_bid*NCHS] = wp_lid;
		if (widmax[wp_lid]<wp_wid) widmax[wp_lid] = wp_wid;
	}

	//===================Get run info============================
	TFile * if_run = new TFile(HOME+"/info/run-info.root");
	TTree * t_run = (TTree*) if_run->Get("t");
	int i_runNo, gasID, runGr, HV, THR;
	char runDu[128];
	double t00, t01, aacut, sumcut;
	t_run->SetBranchAddress("run_number",&i_runNo);
	t_run->SetBranchAddress("gas_mixture_id",&gasID);
	t_run->SetBranchAddress("hv_ch0",&HV);
	t_run->SetBranchAddress("recbe_th_input_bd0",&THR);
	t_run->SetBranchAddress("duration",&runDu);
	t_run->SetBranchAddress("run_grade",&runGr);
	t_run->SetBranchAddress("t00",&t00);
	t_run->SetBranchAddress("t01",&t01);
	t_run->SetBranchAddress("aa",&aacut);
	t_run->SetBranchAddress("sum",&sumcut);
	for(int i = 0; i<t_run->GetEntries(); i++){
		t_run->GetEntry(i);
		if (i_runNo == runNo) break;
	}
	double npair = 17.96;
	TString gastype = "He:C_{2}H_{4}(50:50)";
	if (gasID==1){
		gastype = "He:iC_{4}H_{10}(90:10)";
		npair = 27.96;
	}
	else if (gasID==2){
		gastype = "He:CH_{4}(80:20)";
		npair = 56.10;
	}
	TString duration = runDu;
	const char *sep = ":";
	char * durationSep = strtok(runDu,sep);
	double durationTime = 0;
	double timeunit = 3600;
	while(durationSep){
		durationTime += timeunit*strtol(durationSep,NULL,10);
		timeunit/=60;
		durationSep = strtok(NULL,sep);
	}
	double t0[NBRD];
	t0[0] = t00;
	t0[1] = t01;
	std::cout<<"runNo#"<<runNo<<": "<<gastype<<", "<<runGr<<", "<<duration<<", "<<HV<<" V, "<<THR<<" mV, "<<durationTime<<"sec"<<std::endl;

	//===================Get raw input ROOT file============================
	TChain * c_raw = new TChain("tree","tree");
	c_raw->Add(HOME+Form("/root/run_%0.6d_built.root",runNo));
	int i_adc[NCHT][NSAM];
	c_raw->SetBranchAddress("adc",i_adc);

	//===================Get peak input ROOT file============================
	TChain * c_peak = new TChain("t","t");
	c_peak->Add(HOME+Form("/root/p_%d.root",runNo));
	int triggerNumber;
	int i_nh;
	int i_np[NCHT];
	double i_ped[NCHT];
	double i_aa[NCHT];
	std::vector<std::vector<int> > * i_clk = 0;
	std::vector<std::vector<int> > * i_tdc = 0;
	std::vector<std::vector<int> > * i_peak = 0;
	std::vector<std::vector<double> > * i_sum = 0;
	std::vector<std::vector<int> > * i_width = 0;
	std::vector<std::vector<int> > * i_height = 0;
	std::vector<std::vector<int> > * i_mpn = 0;
	std::vector<std::vector<int> > * i_mpi = 0;
	c_peak->SetBranchAddress("triggerNumber",&triggerNumber);
	c_peak->SetBranchAddress("nh",&i_nh);
	c_peak->SetBranchAddress("np",i_np);
	c_peak->SetBranchAddress("ped",i_ped);
	c_peak->SetBranchAddress("aa",i_aa);
	c_peak->SetBranchAddress("clk",&i_clk);
	c_peak->SetBranchAddress("tdc",&i_tdc);
	c_peak->SetBranchAddress("peak",&i_peak);
	c_peak->SetBranchAddress("sum",&i_sum);
	c_peak->SetBranchAddress("width",&i_width);
	c_peak->SetBranchAddress("height",&i_height);
	c_peak->SetBranchAddress("mpn",&i_mpn);
	c_peak->SetBranchAddress("mpi",&i_mpi);

	//===================Prepare output ROOT file============================
	int o_nHits;
	int o_nLayers;
	std::vector<int> * o_layerID = 0;
	std::vector<int> * o_wireID = 0;
	std::vector<int> * o_type = 0;
	std::vector<int> * o_np = 0;
	std::vector<int> * o_ip = 0;
	std::vector<int> * o_clk = 0;
	std::vector<int> * o_width = 0;
	std::vector<int> * o_peak = 0;
	std::vector<int> * o_height = 0;
	std::vector<int> * o_mpn = 0;
	std::vector<int> * o_mpi = 0;
	std::vector<int> * o_rank = 0;
	std::vector<double> * o_ped = 0;
	std::vector<double> * o_sum = 0;
	std::vector<double> * o_aa = 0;
	std::vector<double> * o_driftT = 0;
	TFile * f = new TFile(HOME+Form("/root/h_%d.",runNo)+suffix+"root","RECREATE");
	TTree * t = new TTree("t","t");
	t->Branch("triggerNumber",&triggerNumber);
	t->Branch("nHits",&o_nHits);
	t->Branch("nLayers",&o_nLayers);
	t->Branch("layerID",&o_layerID);
	t->Branch("wireID",&o_wireID);
	t->Branch("type",&o_type);
	t->Branch("np",&o_np);
	t->Branch("ip",&o_ip);
	t->Branch("clk",&o_clk);
	t->Branch("width",&o_width);
	t->Branch("peak",&o_peak);
	t->Branch("height",&o_height);
	t->Branch("mpn",&o_mpn);
	t->Branch("mpi",&o_mpi);
	t->Branch("rank",&o_rank);
	t->Branch("ped",&o_ped);
	t->Branch("sum",&o_sum);
	t->Branch("aa",&o_aa);
	t->Branch("driftT",&o_driftT);
	o_layerID = new std::vector<int>;
	o_wireID = new std::vector<int>;
	o_type = new std::vector<int>;
	o_ip = new std::vector<int>;
	o_np = new std::vector<int>;
	o_clk = new std::vector<int>;
	o_width = new std::vector<int>;
	o_peak = new std::vector<int>;
	o_height = new std::vector<int>;
	o_mpn = new std::vector<int>;
	o_mpi = new std::vector<int>;
	o_rank = new std::vector<int>;
	o_ped = new std::vector<double>;
	o_sum = new std::vector<double>;
	o_aa = new std::vector<double>;
	o_driftT = new std::vector<double>;

	//===================Prepare Histograms============================
	TCanvas * canvas;
	TPad * pad[NCHS];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	TLine *line_ht[NCHS];
	TLatex *text_ht[NCHS];
	TLine *line_ht2[NCHS];
	TLatex *text_ht2[NCHS];
	for (Int_t i=0; i<NCHS; i++) {
		line_ht[i] = new TLine();
		line_ht[i]->SetLineColor(kRed);
		line_ht[i]->SetLineWidth(0.5);
		text_ht[i] = new TLatex(0,0,"");
		text_ht[i]->SetTextSize(0.04);
		text_ht[i]->SetTextColor(kRed);
		line_ht2[i] = new TLine();
		line_ht2[i]->SetLineColor(kBlack);
		line_ht2[i]->SetLineWidth(0.5);
		text_ht2[i] = new TLatex(0,0,"");
		text_ht2[i]->SetTextSize(0.04);
		text_ht2[i]->SetTextColor(kBlack);
	}
	TLatex * text = new TLatex();
	text->SetTextSize(0.02);
	TString name;
	TString title;

	TH2D* hwf[NCHT];
	TH2D* hst[NCHT];
	TH2D* hat[NCHT];
	TH1D* htdc[NCHT];
	TH1D* ha[NCHT];
	TH1D* hs[NCHT];

	for (int i = 0; i<NCHT; i++){
		if (map_lid[i]==-1){
			title = Form("ch#%d",i%NCHS);
		}
		else{
			title = Form("ch#%d, layer#%d, wire#%d",i%NCHS,map_lid[i],map_wid[i]);
		}
		//hst
		name = Form("hst_%d",i);
		hst[i] = new TH2D(name,title,NBINS,-1000,0,NBINS,MIN_ADCL,MAX_ADCL);
		hst[i]->GetXaxis()->SetTitle("TDC");
		hst[i]->GetYaxis()->SetTitle("ADC (-pedestal) Sum (of Single Peak) [ADC]");
		//hat
		name = Form("hat_%d",i);
		hat[i] = new TH2D(name,title,NBINS,-1000,0,NBINS,MIN_ADCL,MAX_ADCL);
		hat[i]->GetXaxis()->SetTitle("TDC");
		hat[i]->GetYaxis()->SetTitle("ADC (-pedestal) Sum (of All Peaks) [ADC]");
		//wf
		name = Form("hwf_%d",i);
		hwf[i] = new TH2D(name,title,NSAM+10,-10,NSAM,MAX_ADCL-MIN_ADCL,MIN_ADCL,MAX_ADCL);
		hwf[i]->GetXaxis()->SetTitle("SampleID - SampleID_{first hit}");
		hwf[i]->GetYaxis()->SetTitle("ADC");
		//tdc
		name = Form("htdc_%d",i);
		htdc[i] = new TH1D(name,title,100,-900,-800);
		htdc[i]->GetXaxis()->SetTitle("TDC");
		//adc sum
		name = Form("has_%d",i);
		ha[i] = new TH1D(name,title,1000,MIN_ADCL,MAX_ADCL);
		ha[i]->GetXaxis()->SetTitle("ADC (-pedestal) Sum (of All Peaks) [ADC]");
		//adc peak
		name = Form("hap_%d",i);
		hs[i] = new TH1D(name,title,(MAX_ADCL-MIN_ADCL)/2,MIN_ADCL,MAX_ADCL);
		hs[i]->GetXaxis()->SetTitle("ADC (-pedestal) Sum (of Single Peak) [ADC]");

		if (map_lid[i]==-1){
			hst[i]->GetXaxis()->SetAxisColor(kRed);
			hst[i]->GetYaxis()->SetAxisColor(kRed);
			hst[i]->GetXaxis()->SetLabelColor(kRed);
			hst[i]->GetYaxis()->SetLabelColor(kRed);

			hat[i]->GetXaxis()->SetAxisColor(kRed);
			hat[i]->GetYaxis()->SetAxisColor(kRed);
			hat[i]->GetXaxis()->SetLabelColor(kRed);
			hat[i]->GetYaxis()->SetLabelColor(kRed);

			hwf[i]->GetXaxis()->SetAxisColor(kRed);
			hwf[i]->GetYaxis()->SetAxisColor(kRed);
			hwf[i]->GetXaxis()->SetLabelColor(kRed);
			hwf[i]->GetYaxis()->SetLabelColor(kRed);

			htdc[i]->GetXaxis()->SetAxisColor(kRed);
			htdc[i]->GetYaxis()->SetAxisColor(kRed);
			htdc[i]->GetXaxis()->SetLabelColor(kRed);
			htdc[i]->GetYaxis()->SetLabelColor(kRed);

			ha[i]->GetXaxis()->SetAxisColor(kRed);
			ha[i]->GetYaxis()->SetAxisColor(kRed);
			ha[i]->GetXaxis()->SetLabelColor(kRed);
			ha[i]->GetYaxis()->SetLabelColor(kRed);

			hs[i]->GetXaxis()->SetAxisColor(kRed);
			hs[i]->GetYaxis()->SetAxisColor(kRed);
			hs[i]->GetXaxis()->SetLabelColor(kRed);
			hs[i]->GetYaxis()->SetLabelColor(kRed);
		}
		else if (map_lid[i]==0){
			hst[i]->GetXaxis()->SetAxisColor(kGreen);
			hst[i]->GetYaxis()->SetAxisColor(kGreen);
			hst[i]->GetXaxis()->SetLabelColor(kGreen);
			hst[i]->GetYaxis()->SetLabelColor(kGreen);

			hat[i]->GetXaxis()->SetAxisColor(kGreen);
			hat[i]->GetYaxis()->SetAxisColor(kGreen);
			hat[i]->GetXaxis()->SetLabelColor(kGreen);
			hat[i]->GetYaxis()->SetLabelColor(kGreen);

			hwf[i]->GetXaxis()->SetAxisColor(kGreen);
			hwf[i]->GetYaxis()->SetAxisColor(kGreen);
			hwf[i]->GetXaxis()->SetLabelColor(kGreen);
			hwf[i]->GetYaxis()->SetLabelColor(kGreen);

			htdc[i]->GetXaxis()->SetAxisColor(kGreen);
			htdc[i]->GetYaxis()->SetAxisColor(kGreen);
			htdc[i]->GetXaxis()->SetLabelColor(kGreen);
			htdc[i]->GetYaxis()->SetLabelColor(kGreen);

			ha[i]->GetXaxis()->SetAxisColor(kGreen);
			ha[i]->GetYaxis()->SetAxisColor(kGreen);
			ha[i]->GetXaxis()->SetLabelColor(kGreen);
			ha[i]->GetYaxis()->SetLabelColor(kGreen);

			hs[i]->GetXaxis()->SetAxisColor(kGreen);
			hs[i]->GetYaxis()->SetAxisColor(kGreen);
			hs[i]->GetXaxis()->SetLabelColor(kGreen);
			hs[i]->GetYaxis()->SetLabelColor(kGreen);
		}
	}

	// Loop in events
	Long64_t N = c_peak->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" events..."<<std::endl;
	int triggerNumberMax = 0;
	int nPeaks;
	bool layerhit[NLAY];
	int ip = -1;
	double avped[NCHT];
	for ( int ch = 0; ch<NCHT; ch++ ){
		avped[ch] = 0;
	}
	int rank[NSAM];
	//N=1;
	for (Long64_t i = 0;i<N; i++){
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c_peak->GetEntry(i);
		c_raw->GetEntry(i);
		if (triggerNumberMax<triggerNumber) triggerNumberMax = triggerNumber;
		o_nHits = 0;
		o_nLayers = 0;
        o_layerID->clear();
        o_wireID->clear();
        o_type->clear();
        o_ip->clear();
        o_np->clear();
        o_clk->clear();
        o_width->clear();
        o_peak->clear();
        o_height->clear();
        o_mpn->clear();
        o_mpi->clear();
        o_rank->clear();
        o_ped->clear();
        o_sum->clear();
        o_aa->clear();
        o_driftT->clear();
		for (int il = 0; il < NLAY; il++){
			layerhit[il] = false;
		}
		for(int ch = 0; ch<NCHT; ch++){
			nPeaks = i_np[ch];
			if (!isnan(i_ped[ch])&&!isinf(i_ped[ch]))
				avped[ch]+= i_ped[ch];

            // Fill histograms
            if (i_aa[ch]) ha[ch]->Fill(i_aa[ch]);
            int offset;
            if (nPeaks<=0) offset = 0;
            else offset = (*i_clk)[ch][0]; // FIXME: should consider about choosing a proper peak as the center of waveform
            for ( int clk = 0; clk<NSAM; clk++ ){
                hwf[ch]->Fill(clk-offset,i_adc[ch][clk]);
            }

            // sort peaks by sum
			for(int ip = 0;ip<nPeaks; ip++){
				rank[ip] = 0;
			}
			for(int ip = 0;ip<nPeaks; ip++){
				for(int jp = 0;jp<nPeaks; jp++){
					if (jp==ip) continue;
					if ((*i_sum)[ch][ip]<(*i_sum)[ch][jp]){
						rank[ip]++;
					}
				}
			}

			// get peaks;
			for(int ip = 0;ip<nPeaks; ip++){
                // Fill histograms
                if ((*i_mpi)[ch][ip]==0){ // do not consider the secondary peaks in the save packet, just take the first one
                    hs[ch]->Fill((*i_sum)[ch][ip]);
                    hst[ch]->Fill((*i_tdc)[ch][ip],(*i_sum)[ch][ip]);
                    hat[ch]->Fill((*i_tdc)[ch][ip],i_aa[ch]);
                    htdc[ch]->Fill((*i_tdc)[ch][ip]);
                }

                // save hits
                if (map_lid[ch]<0) continue; // only save the channel which is connected
                int bid = ch/NCHS;
                if (map_lid[ch]==0) o_type->push_back(4); // dummy
                else if (map_lid[ch]==NLAY-1) o_type->push_back(3); // guard
                else{
                    if (map_wid[ch]==0)  o_type->push_back(1); // left
                    else if (map_wid[ch]==widmax[ch])  o_type->push_back(2); // right
                    else{
                        o_type->push_back(0); // center
                        layerhit[map_lid[ch]] = true;
                    }
                }
                o_nHits++;
                o_ip->push_back(ip);
                o_driftT->push_back(((*i_tdc)[ch][ip]-t0[bid])/tscale);
                o_layerID->push_back(map_lid[ch]);
                o_wireID->push_back(map_wid[ch]);
                o_np->push_back(nPeaks);
                o_clk->push_back((*i_clk)[ch][ip]);
                o_width->push_back((*i_width)[ch][ip]);
                o_peak->push_back((*i_peak)[ch][ip]);
                o_height->push_back((*i_height)[ch][ip]);
                o_mpn->push_back((*i_mpn)[ch][ip]);
                o_mpi->push_back((*i_mpi)[ch][ip]);
                o_rank->push_back(rank[ip]);
                o_ped->push_back(i_ped[ch]);
                o_sum->push_back((*i_sum)[ch][ip]);
                o_aa->push_back(i_aa[ch]);
			}
		}
		for (int il = 0; il < NLAY; il++ ){
			if (layerhit[il]) o_nLayers++;
		}
		t->Fill();
	}
	t->Write();
	for ( int ch = 0; ch<NCHT; ch++ ){
		avped[ch] = avped[ch]/N;
	}

	// Prepare canvas
	canvas = new TCanvas("c","c",1024,768);
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index] = new TPad(Form("p%d_%d",i,j),Form("p%d_%d",i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
			pad[index]->Draw();
			pad[index]->SetGridx(1);
			pad[index]->SetGridy(1);
		}
	}

	// run summary
	TLatex * text2 = new TLatex();
	text2->SetTextSize(0.02);
	text2->SetText(0.1,0.96,Form("run#%d ",runNo)+gastype+Form(", %d V,%d mV, Grade#%d",HV,THR,runGr)+", "+duration+Form(", %d events, Eff_{daq} = %2.2lf%%, Rate_{tri} = %1.1lfkHz",N,((double)N)/(triggerNumberMax+1)*100,(triggerNumberMax+1)/durationTime/1000));

	// Draw sum VS dt for Board 0
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of Single Peak) VS TDC (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			hst[index]->Draw("COLZ");
			line_ht[index]->SetX1(-1000);
			line_ht[index]->SetY1(sumcut);
			line_ht[index]->SetX2(0);
			line_ht[index]->SetY2(sumcut);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(-100,sumcut+5,Form("%.2lf",sumcut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.st.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.st.b0.png",runNo));
	// Draw sum VS dt for Board 1
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of Single Peak) VS TDC (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			hst[index2]->Draw("COLZ");
			line_ht[index]->SetX1(-100);
			line_ht[index]->SetY1(sumcut);
			line_ht[index]->SetX2(900);
			line_ht[index]->SetY2(sumcut);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(-100,sumcut+5,Form("%.2lf",sumcut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.st.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.st.b1.png",runNo));

	// Draw aa VS dt for Board 0
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of All Peaks) VS TDC (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			hat[index]->Draw("COLZ");
			line_ht[index]->SetX1(-1000);
			line_ht[index]->SetY1(aacut);
			line_ht[index]->SetX2(0);
			line_ht[index]->SetY2(aacut);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(-100,aacut+5,Form("%.2lf",aacut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.at.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.at.b0.png",runNo));
	// Draw h VS dt for Board 1
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of All Peak1) VS TDC (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			hst[index2]->Draw("COLZ");
			line_ht[index]->SetX1(-1000);
			line_ht[index]->SetY1(aacut);
			line_ht[index]->SetX2(0);
			line_ht[index]->SetY2(aacut);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(-100,aacut+5,Form("%.2lf",aacut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.at.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.at.b1.png",runNo));

	// Draw waveform for board 0
	canvas->cd();
	text->SetText(0.1,0.98,"Waveform (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			pad[index]->SetLogz(1);
			hwf[index]->Draw("COLZ");
			double ped = avped[index];
			line_ht2[index]->SetX1(-10);
			line_ht2[index]->SetY1(ped);
			line_ht2[index]->SetX2(32);
			line_ht2[index]->SetY2(ped);
			line_ht2[index]->Draw("SAME");
			text_ht2[index]->SetText(25,ped+5,Form("%.2lf",ped));
			text_ht2[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.wf.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.wf.b0.png",runNo));
	// Draw waveform for board 1
	canvas->cd();
	text->SetText(0.1,0.98,"Waveform (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			hwf[index2]->Draw("COLZ");
			double ped = avped[index2];
			line_ht2[index]->SetX1(-10);
			line_ht2[index]->SetY1(ped);
			line_ht2[index]->SetX2(32);
			line_ht2[index]->SetY2(ped);
			line_ht2[index]->Draw("SAME");
			text_ht2[index]->SetText(25,ped+5,Form("%.2lf",ped));
			text_ht2[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.wf.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.wf.b1.png",runNo));

	// Draw TDC for board 0
	canvas->cd();
	text->SetText(0.1,0.98,"TDC (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			htdc[index]->Draw("");
			int bd = index/48;
			line_ht[index]->SetX1(t0[bd]);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(t0[bd]);
			int max = htdc[index]->GetBinContent(htdc[index]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(t0[bd],max*0.1,Form("%.2lf",t0[bd]));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.tdc.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.tdc.b0.png",runNo));
	// Draw TDC for board 1
	canvas->cd();
	text->SetText(0.1,0.98,"TDC (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			htdc[index2]->Draw("");
			int bd = index/48;
			line_ht[index]->SetX1(t0[bd]);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(t0[bd]);
			int max = htdc[index2]->GetBinContent(htdc[index2]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(t0[bd],max*0.1,Form("%.2lf",t0[bd]));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.tdc.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.tdc.b1.png",runNo));

	// Draw ADC (-pedestal) Sum (of Single Peak) for board 0
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of Single Peak) (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			hs[index]->Draw("");
			line_ht[index]->SetX1(sumcut);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(sumcut);
			int max = hs[index]->GetBinContent(hs[index]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(sumcut,max*0.2,Form("%.2lf",sumcut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.sum.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.sum.b0.png",runNo));
	// Draw ADC (-pedestal) Sum (of Single Peak) for board 1
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of Single Peak) (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			hs[index2]->Draw("");
			line_ht[index]->SetX1(sumcut);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(sumcut);
			int max = hs[index2]->GetBinContent(hs[index2]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(sumcut,max*0.2,Form("%.2lf",sumcut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.sum.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.sum.b1.png",runNo));

	// Draw ADC sum for board 0
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of All Peaks) (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			ha[index]->Draw("");
			line_ht[index]->SetX1(aacut);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(aacut);
			int max = hs[index]->GetBinContent(hs[index]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(aacut,max*0.2,Form("%.2lf",aacut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.aa.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.aa.b0.png",runNo));
	// Draw ADC sum for board 1
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) Sum (of All Peaks) (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			ha[index2]->Draw("");
			line_ht[index]->SetX1(aacut);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(aacut);
			int max = hs[index2]->GetBinContent(hs[index2]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(aacut,max*0.2,Form("%.2lf",aacut));
			text_ht[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.aa.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.aa.b1.png",runNo));

	for (int i = 0; i<NCHT; i++) hst[i]->Write();;
	for (int i = 0; i<NCHT; i++) hat[i]->Write();;
	for (int i = 0; i<NCHT; i++) hwf[i]->Write();;
	for (int i = 0; i<NCHT; i++) htdc[i]->Write();;
	for (int i = 0; i<NCHT; i++) ha[i]->Write();;
	for (int i = 0; i<NCHT; i++) hs[i]->Write();;

	return 0;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
