#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"

#define NSAM 32
#define NCHS 48
#define NCHT 96
#define NBD 2
#define MIN_ADC 50
#define MAX_ADC 750
#define NBINS  256

#define DEBUG

int power2_15 = pow(2,15);

bool connected(int ch);
void print_usage(char* prog_name);
int get_bid(int i);
int get_bid_core(int i);
int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	if (argc>=3) nEventMax = (int)strtol(argv[2],NULL,10);
	std::string suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix=suffix+".";
	}

	//===================Get wire position============================
	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.v3.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bid;
	int wp_wid;
	int wp_lid;
	int wp_ch;
	int map_wid[NCHT];
	int map_lid[NCHT];
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
	}

	//===================Get input ROOT file============================
	TChain * c = new TChain("tree","tree");
	char inputName[128];
	sprintf(inputName,"../root/run_%0.6d_built.root",runNo);
	c->Add(inputName);
	int tdcNhit[NCHT];
	int clockNumberDriftTime[NCHT][NSAM];
	int adc[NCHT][NSAM];
	int driftTime[NCHT][NSAM];
	int triggerNumber;
	int triggerNumberMax = 0;
	c->SetBranchAddress("triggerNumber",&triggerNumber);
	c->SetBranchAddress("tdcNhit",tdcNhit);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("driftTime",driftTime);
	c->SetBranchAddress("adc",adc);

	//===================Get run info============================
	TFile * if_run = new TFile("../info/run-info.root");
	TTree * t_run = (TTree*) if_run->Get("t");
	int i_runNo, gasID, runGr, HV, THR;
	char runDu[128];
	t_run->SetBranchAddress("run_number",&i_runNo);
	t_run->SetBranchAddress("gas_mixture_id",&gasID);
	t_run->SetBranchAddress("hv_ch0",&HV);
	t_run->SetBranchAddress("recbe_th_input_bd0",&THR);
	t_run->SetBranchAddress("duration",&runDu);
	t_run->SetBranchAddress("run_grade",&runGr);
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
	std::cout<<"runNo#"<<runNo<<": "<<gastype<<", "<<runGr<<", "<<duration<<", "<<HV<<" V, "<<THR<<" mV, "<<durationTime<<"sec"<<std::endl;

	//===================Prepare output ROOT file============================
	char outputName[128];
	int o_iHits[NCHT];
	int o_nHits[NCHT];
	int o_pedestalN[NCHT];
	double o_pedestal[NCHT];
	double o_pedestalChi2[NCHT];
	double o_areaall[NCHT];
	double o_area[NCHT];
	double o_time[NCHT];
	int o_height[NCHT];
	sprintf(outputName,("../root/d_%d."+suffix+"root").c_str(),runNo);
	TFile * f = new TFile(outputName,"RECREATE");
	TTree * t = new TTree("t","t");
	t->Branch("triggerNumber",&triggerNumber);
	t->Branch("i",o_iHits,Form("i[%d]/I",NCHT));
	t->Branch("n",o_nHits,Form("n[%d]/I",NCHT));
	t->Branch("p",o_pedestal,Form("p[%d]/D",NCHT));
	t->Branch("pn",o_pedestalN,Form("pn[%d]/I",NCHT));
	t->Branch("pchi2",o_pedestalChi2,Form("pchi2[%d]/D",NCHT));
	t->Branch("aa",o_areaall,Form("aa[%d]/D",NCHT));
	t->Branch("a",o_area,Form("a[%d]/D",NCHT));
	t->Branch("dt",o_time,Form("dt[%d]/D",NCHT));
	t->Branch("h",o_height,Form("h[%d]/I",NCHT));
	//FIXME better not to record
	//t->Branch("clockNumberDriftTime",clockNumberDriftTime,"clockNumberDriftTime[NCHT][NSAM]/I");
	//t->Branch("adc",adc,"adc[NCHT][NSAM]/I");

	//===================Prepare Histograms============================
	TCanvas * canvas;
	TPad * pad[NCHS];
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	TLine *line_ht[NCHS];
	TLine *line_ht2[NCHS];
	TLatex *text_ht[NCHS];
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
	TH2D* hht[NCHT];
	TH1D* htdc[NCHT];
	TH1D* has[NCHT];
	TH1D* hap[NCHT];

	for (int i = 0; i<NCHT; i++){
		if (map_lid[i]==-1){
			title = Form("ch#%d",i%NCHS);
		}
		else{
			title = Form("ch#%d, layer#%d, wire#%d",i%NCHS,map_lid[i],map_wid[i]);
		}
		//hht
		name = Form("hht_%d",i);
		hht[i] = new TH2D(name,title,NBINS,-100,900,NBINS,MIN_ADC,MAX_ADC);
		hht[i]->GetXaxis()->SetTitle("Time [ns]");
		hht[i]->GetYaxis()->SetTitle("Peak Height [ADC]");
		//wf
		name = Form("hwf_%d",i);
		hwf[i] = new TH2D(name,title,NSAM+10,-10,NSAM,MAX_ADC-MIN_ADC,MIN_ADC,MAX_ADC);
		hwf[i]->GetXaxis()->SetTitle("SampleID - SampleID_{first hit}");
		hwf[i]->GetYaxis()->SetTitle("ADC");
		//tdc
		name = Form("htdc_%d",i);
		htdc[i] = new TH1D(name,title,100,-900,-800);
		htdc[i]->GetXaxis()->SetTitle("TDC");
		//adc sum
		name = Form("has_%d",i);
		has[i] = new TH1D(name,title,1000,0,2000);
		has[i]->GetXaxis()->SetTitle("ADC (-pedestal) sum");
		//adc peak
		name = Form("hap_%d",i);
		hap[i] = new TH1D(name,title,(MAX_ADC-MIN_ADC)/2,MIN_ADC,MAX_ADC);
		hap[i]->GetXaxis()->SetTitle("ADC peak");

		if (map_lid[i]==-1){
			hht[i]->GetXaxis()->SetAxisColor(kRed);
			hht[i]->GetYaxis()->SetAxisColor(kRed);
			hht[i]->GetXaxis()->SetLabelColor(kRed);
			hht[i]->GetYaxis()->SetLabelColor(kRed);

			hwf[i]->GetXaxis()->SetAxisColor(kRed);
			hwf[i]->GetYaxis()->SetAxisColor(kRed);
			hwf[i]->GetXaxis()->SetLabelColor(kRed);
			hwf[i]->GetYaxis()->SetLabelColor(kRed);

			htdc[i]->GetXaxis()->SetAxisColor(kRed);
			htdc[i]->GetYaxis()->SetAxisColor(kRed);
			htdc[i]->GetXaxis()->SetLabelColor(kRed);
			htdc[i]->GetYaxis()->SetLabelColor(kRed);

			has[i]->GetXaxis()->SetAxisColor(kRed);
			has[i]->GetYaxis()->SetAxisColor(kRed);
			has[i]->GetXaxis()->SetLabelColor(kRed);
			has[i]->GetYaxis()->SetLabelColor(kRed);

			hap[i]->GetXaxis()->SetAxisColor(kRed);
			hap[i]->GetYaxis()->SetAxisColor(kRed);
			hap[i]->GetXaxis()->SetLabelColor(kRed);
			hap[i]->GetYaxis()->SetLabelColor(kRed);
		}
	}

	// other histograms
	TH1D * hh_noise[NCHT];
	//TH1D * hh_noise2[NCHT];
	TH1D * hh_sigall[NCHT];
	TH1D * hh_sighigh[NCHT];
	TH1D * hh_siglow[NCHT];
	double pedestal[NCHT]={0};
	for ( int i = 0; i<NCHT; i++){
		hh_noise[i] = new TH1D(Form("hh_noise_%d",i),"Peak Height of Noise Peaks",500,200,700);
		//hh_noise2[i] = new TH1D(Form("h_h2_%d",i),"hh_noise2",500,200,700);
		hh_sigall[i] = new TH1D(Form("hh_sigall_%d",i),"hh_sigall",300,200,500);
		hh_sighigh[i] = new TH1D(Form("hh_sighigh_%d",i),"hh_sighigh",50,250,700);
		hh_siglow[i] = new TH1D(Form("hh_siglow_%d",i),"hh_siglow",150,200,250);
	}

	// Loop in events
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" events..."<<std::endl;
	int tdcNhitwire = -1;
	//N=1;
	for (Long64_t i = 0;i<N; i++){
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c->GetEntry(i);
		if (triggerNumberMax<triggerNumber) triggerNumberMax = triggerNumber;
		for(int ch = 0; ch<NCHT; ch++){
			tdcNhitwire = tdcNhit[ch];
			for ( int ihit = 0; ihit<tdcNhitwire; ihit++){
				if (driftTime[ch][ihit]>0) driftTime[ch][ihit]-=power2_15;
				htdc[ch]->Fill(driftTime[ch][ihit]);
				int height=0;
				for(int clk = clockNumberDriftTime[ch][ihit]; clk<(ihit+1>=tdcNhitwire?NSAM:clockNumberDriftTime[ch][ihit+1]); clk++){
					if (adc[ch][clk]<height){
						break;
					}
					height=adc[ch][clk];
				}
				if (driftTime[ch][ihit]<-880)
					hh_noise[ch]->Fill(height);
				//else if (driftTime[ch][ihit]>-400)
				//	hh_noise2[ch]->Fill(height);
				if (driftTime[ch][ihit]<-725&&driftTime[ch][ihit]>-825){
					hh_sigall[ch]->Fill(height);
					hh_sighigh[ch]->Fill(height);
					hh_siglow[ch]->Fill(height);
				}
			}
			int clk;
			double ped = 0;
			for(clk = 0; clk<(tdcNhitwire<=0?NSAM:clockNumberDriftTime[ch][0]-1); clk++){
				ped+=adc[ch][clk];
			}
			if (clk>0) ped/=(double)clk;
			pedestal[ch] += ped;
		}
	}

	// Get t0 and Hmin
	double Hmin[NCHT];
	TF1 *f1 = 0;
	TF1 *f2 = 0;
	TF1 *f3 = 0;
	for ( int i = 0; i<NCHT; i++){
		Hmin[i] = pedestal[i]/N+15;
		/*
		if(hh_sighigh[i]->GetMaximumBin()>5){
			if (f2) delete f2;
			f2 = new TF1(Form("hmin2_%d",i),"gaus",200,300);
			double center = hh_siglow[i]->GetBinCenter(hh_siglow[i]->GetMaximumBin());
			hh_siglow[i]->Fit(Form("hmin2_%d",i),"qN0","",center-8,center+8);
			Hmin[i] = f2->GetParameter(1)+6*f2->GetParameter(2);
			printf("super large gas gain!\n");
			printf("Hmin[%d] = %.2lf; %d\n",i,Hmin[i],hh_sighigh[i]->GetMaximumBin());
		}
		else{
			double h_maxbin_300 = hh_siglow[i]->GetBinCenter(hh_siglow[i]->GetMaximumBin());
			double h_ratio_300 = ((double)hh_siglow[i]->GetBinContent(99))/hh_siglow[i]->GetBinContent(hh_siglow[i]->GetMaximumBin());
			printf("ratio_300 = %.2lf, %.2lf - %.2lf;\n",h_ratio_300,h_maxbin_300,pedestal[i]/N);
			if (h_ratio_300>0.01&&h_maxbin_300-pedestal[i]/N<8){
				if (f2) delete f2;
				f2 = new TF1(Form("hmin2_%d",i),"gaus",200,300);
				double center = hh_siglow[i]->GetBinCenter(hh_siglow[i]->GetMaximumBin());
				f2->SetParameters(hh_siglow[i]->Integral(),center,10);
				hh_siglow[i]->Fit(Form("hmin2_%d",i),"qN0","",center-8,center+8);
				if (f2->GetParameter(1)-center>2&&f2->GetParameter(2)>3)
					Hmin[i] = center+2;
				else
					Hmin[i] = f2->GetParameter(1)+5*f2->GetParameter(2);
				printf("Hmin[%d] = %.2lf; (%.2lf,%.2lf,%.2lf)\n",i,Hmin[i],f2->GetParameter(1),f2->GetParameter(2),center);
			}
			else{
				if (f1) delete f1;
				f1 = new TF1(Form("hmin_%d",i),"gaus",200,700);
				double center = hh_noise[i]->GetBinCenter(hh_noise[i]->GetMaximumBin());
				f1->SetParameters(hh_noise[i]->Integral(),center,10);
				hh_noise[i]->Fit(Form("hmin_%d",i),"qN0","",center-8,center+8);
				double hmin = f1->GetParameter(1)+3*f1->GetParameter(2);

				if (f3) delete f3;
				f3 = new TF1(Form("hmin3_%d",i),"landau",200,500);
				center = hh_sigall[i]->GetBinCenter(hh_sigall[i]->GetMaximumBin());
				f3->SetParameters(hh_sigall[i]->Integral(),center,10);
				hh_sigall[i]->Fit(Form("hmin3_%d",i),"qN0","",center-30,center+60);
				double hmin3 = f3->GetParameter(1)-2.5*f3->GetParameter(2);

				if ((double)hh_noise[i]->GetEntries()>N*8.e-3&&hmin>hmin3)
					Hmin[i] = hmin;
				else
					if (hmin3>600||f3->GetParameter(1)>f1->GetParameter(1)*1.5)
						Hmin[i] = f1->GetParameter(1)+f1->GetParameter(2);
					else
						Hmin[i] = hmin3;


				printf("Hmin[%d] = %.2lf; (%.2lf,%.2lf,%.2lf) (%.2lf,%.2lf,%.2lf)\n",i,Hmin[i],f3->GetParameter(1),f3->GetParameter(2),center,f1->GetParameter(1),f1->GetParameter(2),hh_noise[i]->GetEntries()/((double)N));
			}
		}
		*/
	}

	double t0[NBD];
	double t0_nmax[NBD];
	double t0_temp[NCHT];
	for(int i = 0; i <NBD; i++){
		t0[i] = 0;
		std::vector<int> indice;
		std::vector<int> indice2;
		for (int j = 0; j <NCHT; j++){
			if (connected(j)) indice.push_back(j);
		}
		for ( int j = 0; j<indice.size(); j++){
			int ch = indice[j];
			int maxbin = htdc[ch]->GetMaximumBin();
			double max = htdc[ch]->GetBinContent(maxbin);
			if (t0_nmax[i]<max) t0_nmax[i] = max;
			double min = 0;
			for (int k = 1; k<=5; k++){
				min += htdc[ch]->GetBinContent(k);
			}
			min /= 5.;
			double th = (max-min)/5+min;
			int ibin = 1;
			for (; ibin<=60; ibin++){
				int height = htdc[ch]->GetBinContent(ibin);
				if (height>th) break;
			}
			t0_temp[indice[j]] = htdc[ch]->GetBinCenter(ibin);
			if (max>30){
				t0[i] += t0_temp[indice[j]];
				indice2.push_back(ch);
			}
		}
		int N = indice2.size();
		double average = t0[i]/(N);
		for ( int j = 0; j<indice2.size(); j++){
			if (fabs(average-t0_temp[indice2[j]])>5){
				t0[i]-=t0_temp[indice2[j]];
				N--;
			}
		}
		t0[i]/=N;
	}

	// Loop in events
	for (Long64_t i = 0;i<N; i++){
		//FIXME
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c->GetEntry(i);
		for(int ch = 0; ch<NCHT; ch++){
			// reset
			tdcNhitwire = tdcNhit[ch];
			o_nHits[ch] = tdcNhitwire;
			o_iHits[ch]=-1;
			o_time[ch] = -1e9;
			o_height[ch] = -1e9;
			o_area[ch] = -1e9;

			// Get location
			int bid = get_bid(ch);

			// Get height and etc
			int o_height_all[NCHT][NSAM];
			double o_time_all[NCHT][NSAM];
			double o_area_all[NCHT][NSAM];
			o_pedestal[ch]=0;
			int clk;
			for(clk = 0; clk<(tdcNhitwire<=0?NSAM:clockNumberDriftTime[ch][0]-1); clk++){
				o_pedestal[ch]+=adc[ch][clk];
			}
			o_pedestal[ch]/=clk;
			o_pedestalN[ch] = clk;
			o_pedestalChi2[ch]=0;
			for(clk = 0; clk<(tdcNhitwire<=0?NSAM:clockNumberDriftTime[ch][0]); clk++){
				o_pedestalChi2[ch]+=pow(adc[ch][clk]-o_pedestal[ch],2);
			}
			o_areaall[ch] = 0;
			for ( int ihit = 0; ihit<tdcNhitwire; ihit++){
				if (driftTime[ch][ihit]>0) driftTime[ch][ihit]-=power2_15;
				o_time_all[ch][ihit]=(driftTime[ch][ihit]-t0[bid])/0.96;
				o_height_all[ch][ihit] = 0;
				for(int clk = clockNumberDriftTime[ch][ihit]; clk<(ihit+1>=tdcNhitwire?NSAM:clockNumberDriftTime[ch][ihit+1]); clk++){
					if (adc[ch][clk]<o_height_all[ch][ihit]){
						break;
					}
					o_height_all[ch][ihit]=adc[ch][clk];
				}

				o_area_all[ch][ihit] = 0;
				for(int clk = clockNumberDriftTime[ch][ihit]; clk<(ihit+1>=tdcNhitwire?NSAM:clockNumberDriftTime[ch][ihit+1]); clk++){
					if (clk!=clockNumberDriftTime[ch][ihit]&&adc[ch][clk]<o_pedestal[ch]) break;
					if (adc[ch][clk]>o_pedestal[ch]) o_area_all[ch][ihit] += adc[ch][clk]-o_pedestal[ch];
				}
				o_areaall[ch] += o_area_all[ch][ihit];
			}

			// find the peak
			for ( int ihit = 0; ihit<tdcNhit[ch]; ihit++){
				hap[ch]->Fill(o_height_all[ch][ihit]);
				if (o_height_all[ch][ihit]>Hmin[ch]){
					o_iHits[ch] = ihit;
					o_height[ch] = o_height_all[ch][ihit];
					o_area[ch] = o_area_all[ch][ihit];
					o_time[ch] = o_time_all[ch][ihit];
					break;
				}
			}

			// Fill histograms
			int offset;
			if (tdcNhitwire==0) offset = 0;
			else offset = clockNumberDriftTime[ch][0];
			for ( int clk = 0; clk<NSAM; clk++ ){
				hwf[ch]->Fill(clk-offset,adc[ch][clk]);
			}

			hht[ch]->Fill(o_time_all[ch][0],o_height_all[ch][0]);
			htdc[ch]->Fill(o_time_all[ch][0]);
			if (o_areaall[ch])
				has[ch]->Fill(o_areaall[ch]);
		}
		t->Fill();
	}
	t->Write();

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

	// Draw h VS dt for Board 0
	canvas->cd();
	text->SetText(0.1,0.98,"Peak Height VS Time (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			hht[index]->Draw("COLZ");
			double hmin = Hmin[index];
			line_ht[index]->SetX1(-100);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(900);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(800,hmin+5,Form("%.2lf",hmin));
			text_ht[index]->Draw("SAME");
			double ped = pedestal[index]/N;
			line_ht2[index]->SetX1(-100);
			line_ht2[index]->SetY1(ped);
			line_ht2[index]->SetX2(900);
			line_ht2[index]->SetY2(ped);
			line_ht2[index]->Draw("SAME");
			text_ht2[index]->SetText(600,ped+5,Form("%.2lf",ped));
			text_ht2[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.b0.png",runNo));
	// Draw h VS dt for Board 1
	canvas->cd();
	text->SetText(0.1,0.98,"Peak Height VS Time (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			hht[index2]->Draw("COLZ");
			double hmin = Hmin[index2];
			line_ht[index]->SetX1(-100);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(900);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(800,hmin+5,Form("%.2lf",hmin));
			text_ht[index]->Draw("SAME");
			double ped = pedestal[index2]/N;
			line_ht2[index]->SetX1(-100);
			line_ht2[index]->SetY1(ped);
			line_ht2[index]->SetX2(900);
			line_ht2[index]->SetY2(ped);
			line_ht2[index]->Draw("SAME");
			text_ht2[index]->SetText(600,ped+5,Form("%.2lf",ped));
			text_ht2[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ht.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.ht.b1.png",runNo));

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
			double hmin = Hmin[index];
			line_ht[index]->SetX1(-10);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(32);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(30,hmin+5,Form("%.2lf",hmin));
			text_ht[index]->Draw("SAME");
			double ped = pedestal[index]/N;
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
			double hmin = Hmin[index2];
			line_ht[index]->SetX1(-10);
			line_ht[index]->SetY1(hmin);
			line_ht[index]->SetX2(32);
			line_ht[index]->SetY2(hmin);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(30,hmin+5,Form("%.2lf",hmin));
			text_ht[index]->Draw("SAME");
			double ped = pedestal[index2]/N;
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

	// Draw ADC peak for board 0
	canvas->cd();
	text->SetText(0.1,0.98,"ADC peak (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			hap[index]->Draw("");
			line_ht[index]->SetX1(Hmin[index]);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(Hmin[index]);
			int max = hap[index]->GetBinContent(hap[index]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(Hmin[index],max*0.2,Form("%.2lf",Hmin[index]));
			text_ht[index]->Draw("SAME");
			double ped = pedestal[index]/N;
			line_ht2[index]->SetX1(ped);
			line_ht2[index]->SetY1(0);
			line_ht2[index]->SetX2(ped);
			line_ht2[index]->SetY2(max);
			line_ht2[index]->Draw("SAME");
			text_ht2[index]->SetText(ped,max*0.1,Form("%.2lf",ped));
			text_ht2[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ap.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.ap.b0.png",runNo));
	// Draw ADC peak for board 1
	canvas->cd();
	text->SetText(0.1,0.98,"ADC peak (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			hap[index2]->Draw("");
			line_ht[index]->SetX1(Hmin[index2]);
			line_ht[index]->SetY1(0);
			line_ht[index]->SetX2(Hmin[index2]);
			int max = hap[index2]->GetBinContent(hap[index2]->GetMaximumBin());
			line_ht[index]->SetY2(max);
			line_ht[index]->Draw("SAME");
			text_ht[index]->SetText(Hmin[index2],max*0.2,Form("%.2lf",Hmin[index2]));
			text_ht[index]->Draw("SAME");
			double ped = pedestal[index2]/N;
			line_ht2[index]->SetX1(ped);
			line_ht2[index]->SetY1(0);
			line_ht2[index]->SetX2(ped);
			line_ht2[index]->SetY2(max);
			line_ht2[index]->Draw("SAME");
			text_ht2[index]->SetText(ped,max*0.1,Form("%.2lf",ped));
			text_ht2[index]->Draw("SAME");
		}
	}
	canvas->SaveAs(Form("run%d.ap.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.ap.b1.png",runNo));

	// Draw ADC sum for board 0
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) sum (Board #0)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			pad[index]->cd();
			has[index]->Draw("");
		}
	}
	canvas->SaveAs(Form("run%d.as.b0.pdf",runNo));
	canvas->SaveAs(Form("run%d.as.b0.png",runNo));
	// Draw ADC sum for board 1
	canvas->cd();
	text->SetText(0.1,0.98,"ADC (-pedestal) sum (Board #1)");
	text->Draw("SAME");
	text2->Draw("SAME");
	for (int i = 0; i<8; i++){
		for (int j = 0; j<6; j++){
			int index = j*8+i;
			int index2 = index+NCHS;
			pad[index]->cd();
			has[index2]->Draw("");
		}
	}
	canvas->SaveAs(Form("run%d.as.b1.pdf",runNo));
	canvas->SaveAs(Form("run%d.as.b1.png",runNo));


	for (int i = 0; i<NCHT; i++) hht[i]->Write();;
	for (int i = 0; i<NCHT; i++) hwf[i]->Write();;
	for (int i = 0; i<NCHT; i++) htdc[i]->Write();;
	for (int i = 0; i<NCHT; i++) has[i]->Write();;
	for (int i = 0; i<NCHT; i++) hap[i]->Write();;
	for (int i = 0; i<NCHT; i++) hh_noise[i]->Write();;
	//for (int i = 0; i<NCHT; i++) hh_noise2[i]->Write();;
	for (int i = 0; i<NCHT; i++) hh_sigall[i]->Write();;
	for (int i = 0; i<NCHT; i++) hh_sighigh[i]->Write();;
	for (int i = 0; i<NCHT; i++) hh_siglow[i]->Write();;

	return 0;
}

int get_bid(int i){
	int bid = -1;
	if (i<NCHS) bid = 0;
	else bid = 1;
	return bid;
}

int get_bid_core(int i){
	int bid = -1;
	if (i>=0&&i<6) bid = 0;
	else if (i==15||i==19||i==24) bid = 1;
	else if (i==42||i==46||i==51) bid = 2;
	else if (i>=60&&i<96) bid = 3;
	return bid;
}

bool connected(int ch){
	int chs = ch%NCHS;
	if (chs==22||chs==23||chs==46) return false;
	else return true;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
