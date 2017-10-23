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

#define NBIN1 125
#define NBIN2 50
#define NBIN 350 // NBIN = 2*NBIN1+2*NBIN2

#define NBINX 256

#define NLAY 9
#define NCEL 11
#define NCELA 99

int TSTART=-10;
int TMIN=0;
int TTURN=220;
int TMAX=450;
double STEP1=(TTURN-TMIN)/((double)NBIN1);
double STEP2=(TMAX-TTURN)/((double)NBIN2);

double n2t(int ibin){
	double t;
	if (ibin<NBIN1+NBIN2){// left
		if (ibin<NBIN2){// left end
			t = TMAX-(ibin+0.5)*STEP2;
		}
		else{// left normal
			t = TMAX-NBIN2*STEP2-(ibin-NBIN2+0.5)*STEP1;
		}
	}
	else{// right
		if (ibin<NBIN1*2+NBIN2){// right normal
			t = TMIN+(ibin-NBIN1-NBIN2+0.5)*STEP1;
		}
		else{// right end
			t = TMIN+NBIN1*STEP1+(ibin-NBIN1*2-NBIN2+0.5)*STEP2;
		}
	}
	return t;
}

int t2n(double t,bool isright){
	int ibin;
	if (t>TMAX||t<TMIN) return -1;
	if (t>NBIN1*STEP1+TMIN){ // end
		if (isright) ibin = NBIN1*2+NBIN2+(t-(NBIN1*STEP1+TMIN))/STEP2;
		else ibin = (TMAX-t)/STEP2;
	}
	else{ // normal
		if (isright) ibin = NBIN1+NBIN2+(t-TMIN)/STEP1;
		else ibin = NBIN2+((NBIN1*STEP1+TMIN)-t)/STEP1;
	}
	return ibin;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [suffix] [workMode: 1, i_XXX; 2, t_XXX] [xtType: 2, sym, thr 0; 1, sym; 0, no req] <[TMIN] [TTURN] [TMAX]>\n",name);
}

int main(int argc, char** argv){

	if (argc<5){
	    printUsage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
    TString suffix  = argv[2];
	int workMode = (int)strtol(argv[3],NULL,10);
	int xtType = (int)strtol(argv[4],NULL,10);
	if (argc>=8){
	    TMIN = (int)strtol(argv[5],NULL,10);
	    TTURN = (int)strtol(argv[6],NULL,10);
	    TMAX = (int)strtol(argv[7],NULL,10);
        STEP1=(TTURN-TMIN)/NBIN1;
        STEP2=(TMAX-TTURN)/NBIN2;
	}

	TH1D * h_x[NCELA][NBIN];
	TF1 * f_x = new TF1("f_x","gaus",-10,10);
	for (int i = 0; i<NBIN; i++){
		for (int j = 0; j<NCELA; j++){
			h_x[j][i] = new TH1D(Form("h%d_%d_%d",j/NCEL,j%NCEL,i),Form("Layer%d Cell%d: Drift Distance (Drift Time = %lf)",j/NCEL,j%NCEL,n2t(i)),NBINX,-10,10);
		}
	}
	TF1 * f_left[NCELA];
	TF1 * f_right[NCELA];
	TF1 * f_left_end[NCELA];
	TF1 * f_right_end[NCELA];
	for (int j = 0; j<NCELA; j++){
		f_left[j] = new TF1(Form("f_left_%d_%d",j/NCEL,j%NCEL),"pol5",TMIN,TMAX);
		f_left_end[j] = new TF1(Form("f_left_end_%d_%d",j/NCEL,j%NCEL),"pol2",TMIN,TMAX);
		f_right[j] = new TF1(Form("f_right_%d_%d",j/NCEL,j%NCEL),"pol5",TMIN,TMAX);
		f_right_end[j] = new TF1(Form("f_right_end_%d_%d",j/NCEL,j%NCEL),"pol2",TMIN,TMAX);
	}
	TF1 * fl[NCELA];
	TF1 * fr[NCELA];

	TFile * ofile = new TFile(Form("../info/xt.%d.",runNo)+suffix+".root","RECREATE");
	TTree * otree = new TTree("t","t");
	int lid;
	int wid;
	double x;
	double t;
	double sig;
	int nent;
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
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
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
	TLine * l0 = new TLine(0,-10,0,10);
	l0->SetLineStyle(1);
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

	TChain * ichain = 0;
	for (int il = 1; il<=8; il++){
		std::cout<<"################################"<<std::endl;
		std::cout<<"       "<<il<<std::endl;
		if (ichain) delete ichain;
		ichain = new TChain("t","t");
		ichain->SetMarkerStyle(1);
		ichain->SetMarkerColor(kMagenta);
		ichain->Add(Form("../root/t_%d.layer%d.",runNo,il)+suffix+".root");

		std::vector<int> * layerID = 0;
		std::vector<int> * wireID = 0;
		std::vector<double> * fitD = 0;
		std::vector<double> * driftT = 0;
		double chi2;
		double slz;

		ichain->SetBranchAddress("layerID",&layerID);
		ichain->SetBranchAddress("wireID",&wireID);
		ichain->SetBranchAddress("fitD0",&fitD);
		ichain->SetBranchAddress("driftT",&driftT);
		ichain->SetBranchAddress("chi20",&chi2);
		ichain->SetBranchAddress("slz0",&slz);

		Long64_t N = ichain->GetEntries();
		if (nEventMax&&nEventMax<N) N = nEventMax;
		std::cout<<"Processing "<<N<<" events..."<<std::endl;
		for ( int i = 0 ; i<N; i++){
			ichain->GetEntry(i);
			// FIXME cut for good event
			if (chi2>5) continue;
			if (fabs(slz-0.01)>0.1) continue;
			for (int ihit = 0; ihit<layerID->size(); ihit++){
				if ((*layerID)[ihit]==il){
                    if (xtType>0) (*fitD)[ihit] = fabs((*fitD)[ihit]);
					int idiv = t2n((*driftT)[ihit],(*fitD)[ihit]>0);
					if (idiv>=0){
						// FIXME merge the data from all wires in the same layer into wire 0
//						h_x[il*NCEL+(*wireID)[ihit]][idiv]->Fill((*fitD)[ihit]);
						h_x[il*NCEL][idiv]->Fill((*fitD)[ihit]);
					}
				}
			}
		}
		for (int iw = 0; iw<NCEL; iw++){
			std::cout<<"##"<<il<<" "<<iw<<std::endl;
			lid=il;
			wid=iw;
			int index = il*NCEL+iw;
            // FIXME only fit wire 0 and skip others
			if (iw>0){
				continue;
			}
			double avsig = 0;
			int nxtpoints = 0;
			for (int i = 0; i<NBIN; i++){
				nent = h_x[index][i]->Integral();
				int hmax = h_x[index][i]->GetMaximum();
				double xleft = 0;
				double xright = 0;
				for (int ibin = 1; ibin<=NBINX; ibin++){
					if (!xleft&&h_x[index][i]->GetBinContent(ibin)>hmax/3.) xleft = h_x[index][i]->GetBinCenter(ibin);
					if (xleft&&!xright&&h_x[index][i]->GetBinContent(ibin)>hmax/3.){
						xright = h_x[index][i]->GetBinCenter(ibin);
						break;
					}
				}
				if (nent){
					h_x[index][i]->Fit("f_x","qN0","",xleft,xright);
					//c->cd();
					//h_x[index][i]->Draw();
					//f_x->Draw("SAME");
					//c->SaveAs(Form("dt.%d.%d.%d.%d.png",runNo,il,iw,i));
				}
				x = f_x->GetParameter(1);
				sig = f_x->GetParameter(2);
				t = n2t(i);
				otree->Fill();
				// FIXME: cut for good xt point
				if (nent>50&&fabs(x)<7&&fabs(x)>2){avsig+=sig;nxtpoints++;}
			}
			avsig/=nxtpoints;
			std::cout<<"avsig="<<avsig<<std::endl;
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
				if (x>0&&t<TTURN) {
					vx_right.push_back(x);
					vt_right.push_back(t);
				}
				else if (x<0&&t<TTURN){
					vx_left.push_back(x);
					vt_left.push_back(t);
				}
				else if (t>=TTURN&&x<0){
					vx_left_end.push_back(x);
					vt_left_end.push_back(t);
				}
				else if (t>=TTURN&&x>0){
					vx_right_end.push_back(x);
					vt_right_end.push_back(t);
				}
			}
			if (g_left) delete g_left;
			if (g_right) delete g_right;
			if (g_left_end) delete g_left_end;
			if (g_right_end) delete g_right_end;
			g_left = new TGraph(vx_left.size(),&(vt_left[0]),&(vx_left[0]));
			g_right = new TGraph(vx_right.size(),&(vt_right[0]),&(vx_right[0]));
			g_left_end = new TGraph(vx_left_end.size(),&(vt_left_end[0]),&(vx_left_end[0]));
			g_right_end = new TGraph(vx_right_end.size(),&(vt_right_end[0]),&(vx_right_end[0]));
			g_left->Fit(Form("f_left_%d_%d",il,iw),"qN0","");
			g_right->Fit(Form("f_right_%d_%d",il,iw),"qN0","");
			g_left_end->Fit(Form("f_left_end_%d_%d",il,iw),"qN0","");
			g_right_end->Fit(Form("f_right_end_%d_%d",il,iw),"qN0","");
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
			g_left_end->SetName(Form("gle_%d_%d",il,iw));
			g_left_end->Write();
			g_left->SetName(Form("gl_%d_%d",il,iw));
			g_left->Write();
			g_right_end->SetName(Form("gre_%d_%d",il,iw));
			g_right_end->Write();
			g_right->SetName(Form("gr_%d_%d",il,iw));
			g_right->Write();
			double p0l_end = f_left_end[index]->GetParameter(0);
			double p1l_end = f_left_end[index]->GetParameter(1);
			double p2l_end = f_left_end[index]->GetParameter(2);
			double p0r_end = f_right_end[index]->GetParameter(0);
			double p1r_end = f_right_end[index]->GetParameter(1);
			double p2r_end = f_right_end[index]->GetParameter(2);
            TF1 * f_delta_l = new TF1(Form("dl_%d_%d",il,iw),Form("f_left_%d_%d-x*x*%.7e-x*%.7e-%.7e",il,iw,p2l_end,p1l_end,p0l_end),TMIN,TTURN);
            double tl = 0;
            for (int tmax = TTURN; tmax<TMAX; tmax+=5){
                f_delta_l->SetRange(TMIN,tmax);
                tl = f_delta_l->GetX(0);
                if (tl!=tmax) break;
            }
            TF1 * f_delta_r = new TF1(Form("dr_%d_%d",il,iw),Form("f_right_%d_%d-x*x*%.7e-x*%.7e-%.7e",il,iw,p2r_end,p1r_end,p0r_end),TMIN,TTURN);
            double tr = 0;
            for (int tmax = TTURN; tmax<TMAX; tmax+=5){
                f_delta_r->SetRange(TMIN,tmax);
                tr = f_delta_r->GetX(0);
                if (tr!=tmax) break;
            }
            printf("tl[%d][%d] = %.3e, tr[%d][%d] = %.3e\n",il,iw,tl,il,iw,tr);
            fr[index] = new TF1(Form("fr_%d_%d",il,iw),Form("(x<%.3e)*f_right_%d_%d+(x>=%.3e)*(x*x*%.7e+x*%.7e+%.7e)",tr,il,iw,tr,p2r_end,p1r_end,p0r_end),TMIN,TMAX);
            if (xtType>0)
                fl[index] = new TF1(Form("fl_%d_%d",il,iw),Form("-(x<%.3e)*f_right_%d_%d-(x>=%.3e)*(x*x*%.7e+x*%.7e+%.7e)",tr,il,iw,tr,p2r_end,p1r_end,p0r_end),TMIN,TMAX);
            else
                fl[index] = new TF1(Form("fl_%d_%d",il,iw),Form("(x<%.3e)*f_left_%d_%d+(x>=%.3e)*(x*x*%.7e+x*%.7e+%.7e)",tl,il,iw,tl,p2l_end,p1l_end,p0l_end),TMIN,TMAX);
            fr[index]->SetLineWidth(0.3);
            fl[index]->SetLineWidth(0.3);
            fr[index]->Write();
            fl[index]->Write();
            f_right_end[index]->Write();
            f_right[index]->Write();
            f_left_end[index]->Write();
            f_left[index]->Write();
			//FIXME get all wires in this layer
			//cutstr = Form("chi2<10&&layerID==%d&&wireID==%d&&abs(slZ)<0.15&&abs(inZ)<19",il,iw);
			cutstr = Form("chi20<5&&layerID==%d&&abs(slz0-0.01)<0.1&&abs(inz0)<19",il);
			if (ichain->GetEntries(cutstr)>50){
				pad1->cd();
				ichain->Draw(Form("fitD0:driftT>>h(%d,%d,%d,500,-10,10)",NBIN,TSTART,TMAX),cutstr,"COLZ");
				ichain->Draw("driftD0:driftT",cutstr,"SAME");
				otree->Draw("x:t",Form("lid==%d&&wid==%d&&n>50",il,iw),"PSAME");
				fl[index]->Draw("SAME");
				fr[index]->Draw("SAME");
				g_left->Draw("PSAME");
				g_right->Draw("PSAME");
				g_left_end->Draw("PSAME");
				g_right_end->Draw("PSAME");
				pad2->cd();
				otree->Draw("x:sig>>htemp(1024,0,0.8,1024,-10,10)&&n>50",Form("lid==%d&&wid==%d",il,iw),"P");
				l0->SetX1(avsig),l0->SetX2(avsig);
				l0->Draw("SAME");
				l1->SetX1(avsig*1.1),l1->SetX2(avsig*1.1);
				l1->Draw("SAME");
				pad3->cd();
				otree->Draw(Form("x:n>>htemp(1024,0,%d,1024,-10,10)",nMax),Form("lid==%d&&wid==%d",il,iw),"P");
				canvas->SaveAs(Form("xt.%d.%d.%d.png",runNo,il,iw));
			}
		}
	}

	std::cout<<"finish"<<std::endl;
	otree->Write();
	ofile->Close();

	return 0;
}
