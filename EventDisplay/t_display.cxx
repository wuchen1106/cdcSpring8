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
	TFile * TFile_wirepos = new TFile("../info/wire-position.v3.root");
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
			map_xro[wp_lid][wp_wid] = wp_xro/10.;
			map_yro[wp_lid][wp_wid] = wp_yro/10.;
			map_xhv[wp_lid][wp_wid] = wp_xhv/10.;
			map_yhv[wp_lid][wp_wid] = wp_yhv/10.;
			map_ch[wp_lid][wp_wid] = wp_ch+wp_bid*48;
		}
	}
	double yup = 62.397007;
	double ydown = 52.760011;

	//===================Get ROOT File============================
	//TChain * c = new TChain("t","t");
	TChain * c = new TChain("t","t");
	std::stringstream buf;
	buf.str(""); buf.clear();
	buf<<"../root/t_"<<runNo<<suffix<<".root";
	c->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	int triggerNumber;
	double xup, zup, slx, slz;
	double xdown, zdown;
	double xdowni, zdowni;
	double xupi, zupi, slix, sliz;
	double chi2;
	std::vector<double> * i_driftD = 0;
	std::vector<double> * i_driftT = 0;
	std::vector<double> * i_fitD = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_layerID = 0;
	c->SetBranchAddress("wireID",&i_wireID);
	c->SetBranchAddress("layerID",&i_layerID);
	c->SetBranchAddress("driftD",&i_driftD);
	c->SetBranchAddress("driftT",&i_driftT);
	c->SetBranchAddress("fitD",&i_fitD);
	c->SetBranchAddress("inX",&xup);
	c->SetBranchAddress("inZ",&zup);
	c->SetBranchAddress("slX",&slx);
	c->SetBranchAddress("slZ",&slz);
	c->SetBranchAddress("iniX",&xupi);
	c->SetBranchAddress("iniZ",&zupi);
	c->SetBranchAddress("sliX",&slix);
	c->SetBranchAddress("sliZ",&sliz);
	c->SetBranchAddress("chi2",&chi2);
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
	TEllipse * ewiret2[8][11];
	double wx,wy,wz,dd,ddt;
	double wxro,wyro,wzro;
	double wxhv,wyhv,wzhv;
	wzro = 59.917/2;
	wzhv = -59.917/2;
	int lid,wid;
	TCanvas * ca = new TCanvas("ca","ca",896,1024);
	TPad * p1 = new TPad("p1","p1",0,0.2,1,1);
	TPad * p2 = new TPad("p2","p2",0,0,1,0.2);
	p1->Draw();
	p2->Draw();
	p1->SetGridx(1);
	p1->SetGridy(1);
	p2->SetGridx(1);
	p2->SetGridy(1);
	TH2D * h0 = new TH2D("h0","h0",128,-13,13,128,50,65);
	gStyle->SetOptStat(0);
	TLine * l = new TLine();
	l->SetLineColor(kRed);
	l->SetY1(yup);
	l->SetY2(ydown);
	TLine * l2 = new TLine();
	l2->SetLineColor(kBlue);
	l2->SetY1(yup);
	l2->SetY2(ydown);
	TText * text[8][11];
	for (int lid = 0; lid<8; lid++){
		for (int wid = 0; wid<11; wid++){
			text[lid][wid] = 0;
			ewiret[lid][wid] = 0;
			ewiret2[lid][wid] = 0;
		}
	}

//	TTree_wirepos->SetMarkerStyle(20);
//	TTree_wirepos->SetMarkerSize(0.5);
//	TTree_wirepos->Draw("yc/10.:xc/10.","","SAME");
	std::string suf = ".pdf";
//	for ( int i = 0 ; i<c->GetEntries(); i++){
	for ( int i = 0 ; i<100; i++){
		if (i%100==0) printf("%lf%...\n",(double)i/c->GetEntries()*100);
		c->GetEntry(i);
//		if (i!=352
//		  &&i!=367
//		  &&i!=758
//		  ) continue;
		int ihit = 0;
		int ch = -1;
		for (; ihit<i_driftD->size(); ihit++){
			lid = (*i_layerID)[ihit];
			wid = (*i_wireID)[ihit];
			if (lid==4){
				ch = map_ch[lid][wid];
				break;
			}
		}
		if (ch==-1) continue;
//		if ((*i_fitD)[ihit]>-0.52||(*i_fitD)[ihit]<-0.58||(*i_driftT)[ihit]<170||chi2>1) continue;

		buf.str("");
		buf.clear();
		//FIXME
		buf.str("");
		buf.clear();
		buf<<"Entry#"<<i<<", TriggerNumber#"<<triggerNumber<<", chi2 = "<<chi2;
		h0->SetTitle(buf.str().c_str());
		p1->cd();
		h0->Draw();
		// FIXME
		int nHits = 0;
		for (int lid = 0; lid<8; lid++){
			for (int wid = 0; wid<11; wid++){
				if (ewiret[lid][wid]){
					delete ewiret[lid][wid];
					ewiret[lid][wid] = 0;
				}
				if (ewiret2[lid][wid]){
					delete ewiret2[lid][wid];
					ewiret2[lid][wid] = 0;
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
			zdown = zup + (ydown-yup)*slz;
			xdown = xup + (ydown-yup)*slx;
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
			zdowni = zupi + (ydown-yup)*sliz;
			xdowni = xupi + (ydown-yup)*slix;
			wz = ((yup-wy)*zdowni+(wy-ydown)*zupi)/(yup-ydown);
			//wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
			//wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
			wx = (wxhv+wxro)/2.;
			wy = (wyhv+wyro)/2.;
			ewiret2[lid][wid] = new TEllipse(wx,wy,dd,dd);
			ewiret2[lid][wid]->SetFillStyle(0);
			ewiret2[lid][wid]->SetLineColor(kBlue);
			ewiret2[lid][wid]->Draw("SAME"); // Draw hits.
		}
		l->SetX1(xup);
		l->SetX2(xdown);
		l->Draw("SAME");
		l2->SetX1(xupi);
		l2->SetX2(xdowni);
		l2->Draw("SAME");

		for (int index2 = triggerNumber; ;index2--){
			chain2->GetEntry(index2);
			if (triggerNumber==triggerNumber2) break;
		}

		p2->cd();
		if (gr_waveForm) delete gr_waveForm; gr_waveForm = new TGraph(32,vSample,adc[ch]);
		gr_waveForm->Draw("APL");
		for (Int_t j=0; j<tdcNhit[ch]; j++) {
			textTDC[j]->SetText(clockNumberDriftTime[ch][j],adc[ch][clockNumberDriftTime[ch][j]],Form("%d",(int)(tdc[ch][j])));
			textTDC[j]->Draw();
			markerTDC[j]->SetX(clockNumberDriftTime[ch][j]);
			markerTDC[j]->SetY(adc[ch][clockNumberDriftTime[ch][j]]);
			//printf("markerX = clockNumberDriftTime[%d][%d] = %d\n",ch,j,clockNumberDriftTime[ch][j]);
			//printf("markerY = adc[%d][%d] = %d\n",ch,clockNumberDriftTime[ch][j],adc[ch][clockNumberDriftTime[ch][j]]);
			// FIXME
			markerTDC[j]->SetMarkerColor(kRed);
			//if (adc[ch][j]>Hmin[ch]) markerTDC[i][j]->SetMarkerColor(kRed);
			//else markerTDC[i][j]->SetMarkerColor(kBlue);
			markerTDC[j]->Draw();
		}

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
