#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TChain.h"
#include "TLine.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TF1.h"
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define TSTART -10
#define TMIN 5
#define TTURN 220
#define TMAX 455
#define NBIN1 250
#define STEP1 2
#define NBIN2 100
#define STEP2 4
#define NBIN 330
// SETP1*NBIN1+STEP2*NBIN2 = (TMAX-TMIN)*2

double n2t(int ibin){
	double t;
	if (ibin<(NBIN1+NBIN2)/2){// left
		if (ibin<NBIN2/2){// left end
			t = TMAX-(ibin+0.5)*STEP2;
		}
		else{// left normal
			t = TMAX-NBIN2/2*STEP2-(ibin-NBIN2/2+0.5)*STEP1;
		}
	}
	else{// right
		if (ibin<NBIN1+NBIN2/2){// right normal
			t = TMIN+(ibin-NBIN1/2-NBIN2/2+0.5)*STEP1;
		}
		else{// right end
			t = TMIN+NBIN1/2*STEP1+(ibin-NBIN1-NBIN2/2+0.5)*STEP2;
		}
	}
	return t;
}

int t2n(double t,bool isright){
	int ibin;
	if (t>TMAX||t<TMIN) return -1;
	if (t>NBIN1/2*STEP1+TMIN){ // end
		if (isright) ibin = NBIN1+NBIN2/2+(t-(NBIN1/2*STEP1+TMIN))/STEP2;
		else ibin = (TMAX-t)/STEP2;
	}
	else{ // normal
		if (isright) ibin = NBIN1/2+NBIN2/2+(t-TMIN)/STEP1;
		else ibin = NBIN2/2+((NBIN1/2*STEP1+TMIN)-t)/STEP1;
	}
	return ibin;
}

int main(int argc, char** argv){

	if (argc<2){
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	//if (argc>=4) nEventMax = (int)strtol(argv[3],NULL,10);
	TString suffix = "";
	if (argc>=3){
		suffix  = argv[2];
		suffix="."+suffix;
	}

	TH1D * h_x[88][NBIN];
	TF1 * f_t[88][NBIN];
	for (int i = 0; i<NBIN; i++){
		for (int j = 0; j<88; j++){
			h_x[j][i] = new TH1D(Form("h%d_%d_%d",j/11+1,j%11,i),Form("Layer%d Cell%d: Drift Distance (Drift Time = %lf)",j/11+1,j%11,n2t(i)),256,-10,10);
			f_t[j][i] = new TF1(Form("f%d_%d_%d",j/11+1,j%11,i),"gaus",TSTART,TMAX);
		}
	}
	TF1 * f_left[88];
	TF1 * f_right[88];
	TF1 * f_left_end[88];
	TF1 * f_right_end[88];
	for (int j = 0; j<88; j++){
		f_left[j] = new TF1(Form("f_left_%d_%d",j/11+1,j%11),"pol5",TSTART,TTURN+30);
		f_left_end[j] = new TF1(Form("f_left_end_%d_%d",j/11+1,j%11),"pol2",TTURN-30,TMAX);
		f_right[j] = new TF1(Form("f_right_%d_%d",j/11+1,j%11),"pol5",TSTART,TTURN+30);
		f_right_end[j] = new TF1(Form("f_right_end_%d_%d",j/11+1,j%11),"pol2",TTURN-30,TMAX);
		f_left[j]->SetLineWidth(0.3);
		f_left_end[j]->SetLineWidth(0.3);
		f_right[j]->SetLineWidth(0.3);
		f_right_end[j]->SetLineWidth(0.3);
	}

	TFile *ofile = new TFile(Form("../info/xt.%d",runNo)+suffix+".root","RECREATE");
	TTree * otree = new TTree("t","t");
	TTree * otree2 = new TTree("p","p");
	int tlel;
	int tler;
	int tlml;
	int tlmr;
	int trml;
	int trmr;
	int trel;
	int trer;
	int lid;
	int wid;
	double x;
	double t;
	double sig;
	int nent;
	otree2->Branch("lid",&lid);
	otree2->Branch("wid",&wid);
	otree2->Branch("tlel",&tlel);
	otree2->Branch("tler",&tler);
	otree2->Branch("tlml",&tlml);
	otree2->Branch("tlmr",&tlmr);
	otree2->Branch("trml",&trml);
	otree2->Branch("trmr",&trmr);
	otree2->Branch("trel",&trel);
	otree2->Branch("trer",&trer);
	otree->Branch("x",&x);
	otree->Branch("t",&t);
	otree->Branch("lid",&lid);
	otree->Branch("wid",&wid);
	otree->Branch("sig",&sig);
	otree->Branch("n",&nent);
	otree->SetMarkerStyle(7);
	otree->SetMarkerSize(0.8);
	otree->SetMarkerColor(kBlack);
	otree->SetLineColor(kBlack);

	TCanvas * canvas = new TCanvas("c","c",1440,600);
	TPad * pad1  =new TPad("p1","p1",0,0,0.4,1);
	TPad * pad2  =new TPad("p2","p2",0.4,0,0.7,1);
	TPad * pad3  =new TPad("p3","p3",0.7,0,1,1);
	pad1->Draw();
	pad2->Draw();
	pad3->Draw();
	pad1->SetGridx(1);
	pad1->SetGridy(1);
	pad2->SetGridx(1);
	pad2->SetGridy(1);
	pad3->SetGridx(1);
	pad3->SetGridy(1);
	gStyle->SetStatX(0.6);
	TLine * l0 = new TLine(50,-10,50,10);
	TLine * l1 = new TLine(0,-10,0,10);
	TString cutstr;

	TCanvas * c = new TCanvas("c2","c2",800,600);

	std::vector<double> vt_left;
	std::vector<double> vx_left;
	std::vector<double> vt_right;
	std::vector<double> vx_right;
	std::vector<double> vt_left_end;
	std::vector<double> vx_left_end;
	std::vector<double> vt_right_end;
	std::vector<double> vx_right_end;
	TGraph * g_left = 0;
	TGraph * g_right = 0;
	TGraph * g_left_end = 0;
	TGraph * g_right_end = 0;

	for (int il = 1; il<=8; il++){
		std::cout<<"################################"<<std::endl;
		std::cout<<"       "<<il<<std::endl;
		TChain * ichain = new TChain("t","t");
		ichain->SetMarkerStyle(1);
		ichain->SetMarkerColor(kMagenta);
		ichain->Add(Form("../root/t_%d.layer%d",runNo,il)+suffix+".root");

		std::vector<int> * layerID = 0;
		std::vector<int> * wireID = 0;
		std::vector<double> * fitD = 0;
		std::vector<double> * driftT = 0;
		double chi2;
		double slz;

		ichain->SetBranchAddress("layerID",&layerID);
		ichain->SetBranchAddress("wireID",&wireID);
		ichain->SetBranchAddress("fitD",&fitD);
		ichain->SetBranchAddress("driftT",&driftT);
		ichain->SetBranchAddress("chi2",&chi2);
		ichain->SetBranchAddress("slZ",&slz);

		Long64_t N = ichain->GetEntries();
		if (nEventMax&&nEventMax<N) N = nEventMax;
		std::cout<<"Processing "<<N<<" events..."<<std::endl;
		for ( int i = 0 ; i<N; i++){
			if (i%1000==0) printf("%lf%...\n",(double)i/N*100);
			ichain->GetEntry(i);
			// FIXME
			if (chi2>5) continue;
			if (fabs(slz-0.01)>0.1) continue;
			for (int ihit = 0; ihit<layerID->size(); ihit++){
				if ((*layerID)[ihit]==il){
					int idiv = t2n((*driftT)[ihit],(*fitD)[ihit]>0);
					if (idiv>=0){
						// FIXME
//						h_x[(il-1)*11+(*wireID)[ihit]][idiv]->Fill((*fitD)[ihit]);
						h_x[(il-1)*11][idiv]->Fill((*fitD)[ihit]);
					}
				}
			}
		}
		for (int iw = 0; iw<11; iw++){
			std::cout<<"##"<<il<<" "<<iw<<std::endl;
			lid=il;
			wid=iw;
			int index = (il-1)*11+iw;
			// FIXME
			if (iw>0){
				otree2->Fill();
				f_left[index] = f_left[index-iw];
				f_right[index] = f_right[index-iw];
				f_left_end[index] = f_left_end[index-iw];
				f_right_end[index] = f_right_end[index-iw];
				continue;
			}
			double avsig = 0;
			int nxtpoints = 0;
			double tturnl = 0;
			double tturnr = 0;
			for (int i = 0; i<NBIN; i++){
				nent = h_x[index][i]->Integral();
				int hmax = h_x[index][i]->GetMaximum();
				double tleft = 0;
				double tright = 0;
				for (int ibin = 1; ibin<=256; ibin++){
					if (!tleft&&h_x[index][i]->GetBinContent(ibin)>hmax/3.) tleft = h_x[index][i]->GetBinCenter(ibin);
					if (tleft&&!tright&&h_x[index][i]->GetBinContent(ibin)>hmax/3.){
						tright = h_x[index][i]->GetBinCenter(ibin);
						break;
					}
				}
				if (nent){
					h_x[index][i]->Fit(Form("f%d_%d_%d",il,iw,i),"qN0","",tleft,tright);
					c->cd();
					//h_x[index][i]->Draw();
					//f_t[index][i]->Draw("SAME");
					//c->SaveAs(Form("dt.%d.%d.%d.%d.png",runNo,il,iw,i));
				}
				x = f_t[index][i]->GetParameter(1);
				sig = f_t[index][i]->GetParameter(2);
				t = n2t(i);
				otree->Fill();
				if (nent>50&&fabs(x)<7&&fabs(x)>2){avsig+=sig;nxtpoints++;}
				if (!tturnl&&nent>50&&x>=-7.5) tturnl = t;
				if (!tturnr&&nent>50&&x>=7.5) tturnr = t;
			}
			if (!tturnl) tturnl = TTURN;
			if (!tturnr) tturnr = TTURN;
			avsig/=nxtpoints;
			std::cout<<"avsig="<<avsig<<", tturnl="<<tturnl<<", tturnr = "<<tturnr<<std::endl;
			int Nentries = otree->GetEntries();
			vx_left.clear();
			vt_left.clear();
			vx_right.clear();
			vt_right.clear();
			vx_left_end.clear();
			vt_left_end.clear();
			vx_right_end.clear();
			vt_right_end.clear();
			int nMax = 0;
			for (int i = 0; i<NBIN; i++){
				otree->GetEntry(Nentries-NBIN+i);
				if (nent>nMax) nMax = nent;
				if (nent<50) continue;
				if (fabs(x)<7&&fabs(x)>2&&sig>avsig*1.1) continue;
				if (x>0&&t<tturnr) {
					vx_right.push_back(x);
					vt_right.push_back(t);
				}
				else if (x<0&&t<tturnl){
					vx_left.push_back(x);
					vt_left.push_back(t);
				}
				else if (t>=tturnl&&x<0){
					vx_left_end.push_back(x);
					vt_left_end.push_back(t);
				}
				else if (t>=tturnr&&x>0){
					vx_right_end.push_back(x);
					vt_right_end.push_back(t);
				}
			}
			if (vt_left_end.size()){
				tlel=vt_left_end[0]+STEP1/2.;
				tler=vt_left_end[vt_left_end.size()-1]-STEP1/2.;
			}
			else{
				tlel=0;
				tler=0;
			}
			if (vt_left.size()){
				tlml=vt_left[0]+STEP1/2.;
				tlmr=vt_left[vt_left.size()-1]-STEP1/2.;
			}
			else{
				tlml=0;
				tlmr=0;
			}
			if (vt_right_end.size()){
				trel=vt_right_end[0]-STEP1/2.;
				trer=vt_right_end[vt_right_end.size()-1]+STEP1/2.;
			}
			else{
				trel=0;
				trer=0;
			}
			if (vt_right.size()){
				trml=vt_right[0]-STEP1/2.;
				trmr=vt_right[vt_right.size()-1]+STEP1/2.;
			}
			else{
				trml=0;
				trmr=0;
			}
			otree2->Fill();
			if (g_left) delete g_left;
			if (g_right) delete g_right;
			if (g_left_end) delete g_left_end;
			if (g_right_end) delete g_right_end;
			g_left = new TGraph(vx_left.size(),&(vt_left[0]),&(vx_left[0]));
			g_right = new TGraph(vx_right.size(),&(vt_right[0]),&(vx_right[0]));
			g_left_end = new TGraph(vx_left_end.size(),&(vt_left_end[0]),&(vx_left_end[0]));
			g_right_end = new TGraph(vx_right_end.size(),&(vt_right_end[0]),&(vx_right_end[0]));
			g_left->Fit(Form("f_left_%d_%d",il,iw),"qN0","");
			f_left[index]->SetRange(tlmr,tlml);
			g_right->Fit(Form("f_right_%d_%d",il,iw),"qN0","");
			f_right[index]->SetRange(trml,trmr);
			g_left_end->Fit(Form("f_left_end_%d_%d",il,iw),"qN0","");
			f_left_end[index]->SetRange(tler,tlel);
			g_right_end->Fit(Form("f_right_end_%d_%d",il,iw),"qN0","");
			f_right_end[index]->SetRange(trel,trer);
			g_left->SetMarkerColor(kRed);
			g_left->SetMarkerStyle(7);
			g_left->SetMarkerSize(0.5);
			g_right->SetMarkerColor(kRed);
			g_right->SetMarkerStyle(7);
			g_right->SetMarkerSize(0.5);
			g_left_end->SetMarkerColor(kRed);
			g_left_end->SetMarkerStyle(7);
			g_left_end->SetMarkerSize(0.5);
			g_right_end->SetMarkerColor(kRed);
			g_right_end->SetMarkerStyle(7);
			g_right_end->SetMarkerSize(0.5);
			//FIXME
			//cutstr = Form("chi2<10&&layerID==%d&&wireID==%d&&abs(slZ)<0.15&&abs(inZ)<19",il,iw);
			cutstr = Form("chi2<5&&layerID==%d&&abs(slZ-0.01)<0.1&&abs(inZ)<19",il);
			if (ichain->GetEntries(cutstr)>50){
				pad1->cd();
				ichain->Draw(Form("fitD:driftT>>h(%d,%d,%d,500,-10,10)",NBIN,TSTART,TMAX),cutstr,"COLZ");
				ichain->Draw("driftD:driftT",cutstr,"SAME");
				otree->Draw("x:t",Form("lid==%d&&wid==%d&&n>50",il,iw),"PSAME");
				f_left[(il-1)*11+iw]->Draw("SAME");
				f_right[(il-1)*11+iw]->Draw("SAME");
				f_left_end[(il-1)*11+iw]->Draw("SAME");
				f_right_end[(il-1)*11+iw]->Draw("SAME");
				g_left->Draw("PSAME");
				g_right->Draw("PSAME");
				g_left_end->Draw("PSAME");
				g_right_end->Draw("PSAME");
				pad2->cd();
				otree->Draw("x:sig>>htemp(1024,0,0.8,1024,-10,10)&&n>50",Form("lid==%d&&wid==%d",il,iw),"P");
				l1->SetX1(avsig*1.1),l1->SetX2(avsig*1.1);
				l1->Draw("SAME");
				pad3->cd();
				otree->Draw(Form("x:n>>htemp(1024,0,%d,1024,-10,10)",nMax),Form("lid==%d&&wid==%d",il,iw),"P");
				l0->Draw("SAME");
				canvas->SaveAs(Form("xt.%d.%d.%d.png",runNo,il,iw));
			}
		}
	}

//	for (int i = 0; i<NBIN; i++){
//		for (int j = 0; j<88; j++){
//			h_x[j][i]->Write();
//		}
//	}
	for (int j = 0; j<88; j++){
		f_left[j]->Write();
		f_right[j]->Write();
		f_left_end[j]->Write();
		f_right_end[j]->Write();
	}
	std::cout<<"finish"<<std::endl;
	otree->Write();
	otree2->Write();
	ofile->Close();

	return 0;
}
