#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "header.hxx"

#define NBINS    20
#define MAXTRANC 5

//===================Chamber Parameter============================
double  npair_per_cm = 0;
TF1 * f_left0 = 0;
TF1 * f_right0 = 0;
TF1 * f_left = 0;
TF1 * f_right = 0;
TF1 * f_left_mid = 0;
TF1 * f_right_mid = 0;
TF1 * f_left_cent = 0;
TF1 * f_right_cent = 0;
TF1 * f_left_end = 0;
TF1 * f_right_end = 0;
int xtType = 0;
double t7l = 0;
double t8l = 0;
double t7r = 0;
double t8r = 0;
double W = 0;

//==================About RECBE======================
TF1 * fADC2ChargeFunction = 0;

//==================About Scintillator======================
int geoSetup = 0; // 0: normal; 1: finger
// normal scintillator
double sciYup = 0;
double sciYdown = 0;
double sciHL = 0;
double sciHW = 0;

//===================for get dist============================
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

//===================for binning============================
double mTmin = -25-1/0.96/2; // t range for one x bin
double mTmax = 800+1/0.96/2;
int    mNbint = 792+1;
double mXmax = 10; // x range for one t bin
int    mNbinx = 256;
int    mNbinRes = 256;
double minchi2p = 1;
double closestchi2 = 1e9;
double maxRes = 2;

double ADC2Charge(double adc); // fC
double getGG(double aa, double slx, double slz);
int t2d(double t, double & d, bool isRight);
double findFirstX(TF1 * f, double val, double xmin, double xmax, double delta);
void doFit(TH1D * h,double leftRatio = 1/3., double rightRatio = 1/3., double leftEnd = -maxRes, double rightEnd = maxRes);
void printUsage(char * name);

int main(int argc, char** argv){

	//=================================================Get options========================================================
	int runNo = -1;
    int testLayer = 4;
	xtType = 0;
    if (argc>=5)
		xtType = (int)strtol(argv[4],NULL,10);
    double maxchi2 = 1;
    if (argc>=6)
        maxchi2 = (double)strtod(argv[5],NULL);
    int nHitsMax = 0;
    if (argc>=7)
        nHitsMax = (int)strtol(argv[6],NULL,10);
    int aaCut = 0;
    if (argc>=8)
        aaCut = (int)strtol(argv[7],NULL,10);
	int savehists = 0;
    if (argc>=9)
        savehists = (int)strtol(argv[8],NULL,10);
    int debugLevel = 0;
    if (argc>=10)
        debugLevel = (int)strtol(argv[9],NULL,10);
    int iEntryStart = 0;
    if (argc>=11)
        iEntryStart = (int)strtol(argv[10],NULL,10);
    int iEntryStop = 0;
    if (argc>=12)
        iEntryStop = (int)strtol(argv[11],NULL,10);
    int nHitsSmin = 7;
    printf("##############%s with %d Parameters##################\n",argv[0],argc);
    printf("xtType:       %d\n",xtType);
    printf("maxchi2     = %.3e\n",maxchi2);
    printf("test layer:   %d\n",testLayer);
    printf("nHits max   = %d\n",nHitsMax);
    printf("Q cut       = %d\n",aaCut);
    printf("savehists:    %s\n",savehists?"yes":"no");
    printf("debug       = %d\n",debugLevel);
    printf("Entries:     [%d~%d]\n",iEntryStart,iEntryStop);
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

	//=================================================Get related info========================================================
	// get run info
	int gasID = 1;
	int HV = 1800;
	int THR = 0;
	TString gasName;
	if (gasID==0){
		gasName = "C2H6";
	}
	else if (gasID==1){
		gasName = "C4H10";
	}
	else if (gasID==2){
		gasName = "CH4";
	}

    // get XT file
    TFile * XTFile = new TFile(Form("%s/info/xt.%s.%d.root",HOME.Data(),gasName.Data(),HV));
    if (!XTFile){
    	fprintf(stderr,"Cannot find XT file!\n");
    	return 0;
	}
	f_left0 = (TF1*) XTFile->Get("fl_0");
	f_right0 = (TF1*) XTFile->Get("fr_0");
	f_left = (TF1*) XTFile->Get(Form("flc_%d",testLayer));
	f_right = (TF1*) XTFile->Get(Form("frc_%d",testLayer));
	f_left_cent = (TF1*) XTFile->Get(Form("flce_%d",testLayer));
	f_right_cent = (TF1*) XTFile->Get(Form("frce_%d",testLayer));
	f_left_mid = (TF1*) XTFile->Get(Form("flm_%d",testLayer));
	f_right_mid = (TF1*) XTFile->Get(Form("frm_%d",testLayer));
	f_left_end = (TF1*) XTFile->Get(Form("fle_%d",testLayer));
	f_right_end = (TF1*) XTFile->Get(Form("fre_%d",testLayer));
    if (!f_left||!f_right|!f_left0||!f_right0||!f_left_mid||!f_right_mid||!f_left_cent||!f_right_cent||!f_left_end||!f_right_end){
    	fprintf(stderr,"Cannot find XT functions!\n");
    	return 0;
	}
	if (xtType%10==0||xtType%10==1){
		t7r = findFirstX(f_right,7,-10,800,10);
		t8r = findFirstX(f_right,8,-10,800,10);
		t7l = findFirstX(f_left,-7,-10,800,10);
		t8l = findFirstX(f_left,-8,-10,800,10);
	}
	else{
		t7r = findFirstX(f_right_mid,7,-10,800,10);
		t8r = findFirstX(f_right_mid,8,-10,800,10);
		t7l = findFirstX(f_left_mid,-7,-10,800,10);
		t8l = findFirstX(f_left_mid,-8,-10,800,10);
	}
    double minDT = 0;
    if (xtType%10==3)
		minDT = f_left_mid->GetXmin()>f_right_mid->GetXmin()?f_right_mid->GetXmin():f_left_mid->GetXmin();
    else
		minDT = f_left->GetXmin()>f_right->GetXmin()?f_right->GetXmin():f_left->GetXmin();
    double maxDT = 0;
	if (xtType/10==0){
		if (xtType%10==0)
			maxDT = f_left0->GetXmax()<f_right0->GetXmax()?f_right0->GetXmax():f_left0->GetXmax();
		else if (xtType%10==1)
			maxDT = f_left->GetXmax()<f_right->GetXmax()?f_right->GetXmax():f_left->GetXmax();
		else
			maxDT = f_left_mid->GetXmax()<f_right_mid->GetXmax()?f_right_mid->GetXmax():f_left_mid->GetXmax();
	}
	else
		maxDT = t8l<t8r?t8r:t8l;

	//==============================================Prepare input file & output variables=================================================
    // input file
	//----------------------------------Set input file--------------------------------------------
	double i_dca;
	double i_driftTime;
	TChain * ichain = new TChain("t","t");
	ichain->Add(Form("%s/Input/garfXT.%s.%d.root",HOME.Data(),gasName.Data(),HV));
	ichain->SetBranchAddress("dca",&i_dca);
	ichain->SetBranchAddress("driftTime",&i_driftTime);

	TFile * ofile = new TFile(Form("%s/root/res_%s.%d.root",HOME.Data(),gasName.Data(),HV),"RECREATE");
	TTree * otree = new TTree("t","t");
	int o_ibin;
	double o_xmin;
	double o_xmid;
	double o_xmax;
	double o_xres;
	double o_xreserr;
	double o_xrms;
	double o_xrmserr;
	double o_xeff;
	double o_xeff3sig;
	double o_xeff5sig;
	double o_xeff500um;
	double o_xeff1mm;
	double o_xoff;
	double o_dres;
	double o_dreserr;
	double o_drms;
	double o_drmserr;
	double o_deff3sig;
	double o_deff5sig;
	double o_deff500um;
	double o_deff1mm;
	double o_doff;
	double o_nx;
	double o_nxh;
	double o_nd;
	otree->Branch("ibin",&o_ibin);
	otree->Branch("xmin",&o_xmin);
	otree->Branch("xmid",&o_xmid);
	otree->Branch("xmax",&o_xmax);
	otree->Branch("xres",&o_xres);
	otree->Branch("xreserr",&o_xreserr);
	otree->Branch("xrms",&o_xrms);
	otree->Branch("xrmserr",&o_xrmserr);
	otree->Branch("xeff",&o_xeff);
	otree->Branch("xeff3sig",&o_xeff3sig);
	otree->Branch("xeff5sig",&o_xeff5sig);
	otree->Branch("xeff500um",&o_xeff500um);
	otree->Branch("xeff1mm",&o_xeff1mm);
	otree->Branch("xoff",&o_xoff);
	otree->Branch("dres",&o_dres);
	otree->Branch("dreserr",&o_dreserr);
	otree->Branch("drms",&o_drms);
	otree->Branch("drmserr",&o_drmserr);
	otree->Branch("deff3sig",&o_deff3sig);
	otree->Branch("deff5sig",&o_deff5sig);
	otree->Branch("deff500um",&o_deff500um);
	otree->Branch("deff1mm",&o_deff1mm);
	otree->Branch("doff",&o_doff);
	otree->Branch("nx",&o_nx);
	otree->Branch("nxh",&o_nxh);
	otree->Branch("nd",&o_nd);

	//=================================================Prepare histograms, function and counters================================================
	TF1 * f_res = new TF1("fres","gaus",-maxRes,maxRes);
	TH1I * h_nHits = new TH1I("hnHits","Number of TDC hits in each event",100,0,100);
	TH1I * h_DOF = new TH1I("hDOF","Number of DOF",5,0,5);
	TH1D * h_chi2 = new TH1D("hchi2","#chi^{2} of fitting",256,0,10);
	TH1D * h_chi2p = new TH1D("hchi2p","p value of fitting",256,0,1);
	TH1D * h_DOCA = new TH1D("hDOCA","DOCA with Left/Right",mNbinx,-mXmax,mXmax);
	TH1D * h_DOCAb = new TH1D("hDOCAb","DOCA without Left/Right",mNbinx/2,0,mXmax);
	TH1D * h_DriftD = new TH1D("hDriftD","Drift distance with Left/Right",mNbinx,-mXmax,mXmax);
	TH1D * h_DriftDb = new TH1D("hDriftDb","Drift distance without Left/Right",mNbinx/2,0,mXmax);
	TH2D * h_aaVST = new TH2D("haaVST","ADC sum VS driftT",mNbint,mTmin,mTmax,200,-50,550);
	TH2D * h_aaVSD = new TH2D("haaVSD","ADC sum VS driftD",mNbinx,0,mXmax,200,-50,550);
	TH2D * h_ggVSX = new TH2D("hggVSX","Gas gain VS DOCA",mNbinx,-mXmax,mXmax,256,0,2e5);
	TH2D * h_xt = new TH2D("hxt","Time space relation",mNbint,mTmin,mTmax,mNbinx,-mXmax,mXmax);
	TH2D * h_tx = new TH2D("htx","Space time relation",mNbinx,-mXmax,mXmax,mNbint,mTmin,mTmax);
	TH2D * h_resVSX = new TH2D("hresVSX","Residual VS DOCA",mNbinx,-mXmax,mXmax,mNbinRes,-maxRes,maxRes);
	TH2D * h_resVSD = new TH2D("hresVSD","Residual VS drift distance",mNbinx,-mXmax,mXmax,mNbinRes,-maxRes,maxRes);
	TH1D * h_resD[NBINS];
	TH1D * h_resX[NBINS];
	TH1D * h_dedx[MAXTRANC];
	for (int i = 0; i<NBINS; i++){
		double xmin = mXmax*(i)/NBINS;
		double xmax = mXmax*(i+1)/NBINS;
		h_resD[i] = new TH1D(Form("hresD%d",i),Form("Resolution with drift distance in [%.1f,%.1f] mm",xmin,xmax),mNbinRes,-maxRes,maxRes);
		h_resD[i]->GetXaxis()->SetTitle("Residual [mm]");
		h_resX[i] = new TH1D(Form("hresX%d",i),Form("Resolution with DOCA in [%.1f,%.1f] mm",xmin,xmax),mNbinRes,-maxRes,maxRes);
		h_resX[i]->GetXaxis()->SetTitle("Residual [mm]");
	}
	for (int i = 0; i<MAXTRANC; i++){
		h_dedx[i] = new TH1D(Form("hdedx%d",i),Form("dEdX with %d hits omitted",i),128,0,3);
		h_dedx[i]->GetXaxis()->SetTitle("dE/dX [keV/cm]");
	}
	h_nHits->GetXaxis()->SetTitle("Number of Hits");
	h_DOF->GetXaxis()->SetTitle("DOF");
	h_chi2->GetXaxis()->SetTitle("#chi^{2}");
	h_chi2p->GetXaxis()->SetTitle("p value");
	h_DOCA->GetXaxis()->SetTitle("DOCA [mm]");
	h_DOCAb->GetXaxis()->SetTitle("DOCA with left(-)/right(+) [mm]");
	h_DriftD->GetXaxis()->SetTitle("Drift distance [mm]");
	h_DriftDb->GetXaxis()->SetTitle("Drift distance with left(-)/right(+) [mm]");
	h_aaVST->GetXaxis()->SetTitle("Drift time [ns]");
	h_aaVST->GetYaxis()->SetTitle("ADC sum");
	h_aaVSD->GetXaxis()->SetTitle("Drift distance [mm]");
	h_aaVSD->GetYaxis()->SetTitle("ADC sum");
	h_ggVSX->GetXaxis()->SetTitle("DOCA [mm]");
	h_ggVSX->GetYaxis()->SetTitle("Gas gain");
	h_xt->GetXaxis()->SetTitle("Drift time [ns]");
	h_xt->GetYaxis()->SetTitle("DOCA [mm]");
	h_tx->GetXaxis()->SetTitle("DOCA [mm]");
	h_tx->GetYaxis()->SetTitle("Drift time [ns]");
	h_resVSX->GetXaxis()->SetTitle("DOCA [mm]");
	h_resVSX->GetYaxis()->SetTitle("Residual [mm]");
	h_resVSD->GetXaxis()->SetTitle("Drift distance [mm]");
	h_resVSD->GetYaxis()->SetTitle("Residual [mm]");
	std::vector<double> v_xx;
	std::vector<double> v_xxerr;
	std::vector<double> v_dx;
	std::vector<double> v_dxerr;
	std::vector<double> v_xeff;
	std::vector<double> v_xeff500um;
	std::vector<double> v_xres;
	std::vector<double> v_xrms;
	std::vector<double> v_xreserr;
	std::vector<double> v_xrmserr;
	std::vector<double> v_xoff;
	std::vector<double> v_deff;
	std::vector<double> v_deff500um;
	std::vector<double> v_dres;
	std::vector<double> v_drms;
	std::vector<double> v_dreserr;
	std::vector<double> v_drmserr;
	std::vector<double> v_doff;

	int N_ALL = 0;
	int N_CUT1 = 0;
	int N_CUT2 = 0;
	int N_CUT3 = 0;
	int N_CUT4 = 0;
	int N_BIN[NBINS] = {0};

	//=================================================Loop in events====================================================
	Long64_t N = ichain->GetEntries();
	if (!iEntryStart&&!iEntryStop){
		iEntryStart = 0;
		iEntryStop = N-1;
	}
	if (debugLevel>0) {printf("Processing %d events\n",N);fflush(stdout);}
	for ( int iEntry = iEntryStart ; iEntry<=iEntryStop; iEntry++){
		if (iEntry%10000==0) printf("%d\n",iEntry);
		if (debugLevel>=20) printf("Entry %d: \n",iEntry);
		ichain->GetEntry(iEntry);

		double fd = i_dca*10; // cm->mm
		double dt = i_driftTime;
		double dd = 0;
		int status = t2d(dt,dd,fd>0);
		if (!status) continue; // out of range
		h_DriftD->Fill(dd);
		h_DriftDb->Fill(fabs(dd));
		int ibx = fabs(fd)/mXmax*NBINS;
		int ib = fabs(dd)/mXmax*NBINS;
		//dd = (*i_driftD)[ihit];
		double res = fabs(dd) - fabs(fd);
		// xt
		h_tx->Fill(fd,dt);
		h_xt->Fill(dt,fd);
		// res VS x/d
		h_resVSX->Fill(fd,res);
		h_resVSD->Fill(dd,res);
		// res
		h_resD[ib]->Fill(res);
		h_resX[ibx]->Fill(res);
		N_BIN[ibx]++;
	}

	//=================================================Get bin by bin information====================================================
	int ibinl = 0;
	int ibinr = 0;
	TCanvas * canv_bin = 0;
	if (savehists){
		canv_bin = new TCanvas("canv_bin","canv_bin",1024,768);
		gStyle->SetPalette(1);
		gStyle->SetOptStat(0);
		gStyle->SetPadTickX(1);
		gStyle->SetPadTickY(1);
		gStyle->SetOptFit(1);
		gPad->SetGridx(1);
		gPad->SetGridy(1);
	}
	int NusedD = 0;
	int NusedX = 0;
	double averageEffD = 0;
	double averageResD = 0;
	double averageRMSD = 0;
	double averageEffX = 0;
	double averageResX = 0;
	double averageRMSX = 0;
	double bestEffD = 0;
	double bestResD = 1e6;
	double bestRMSD = 1e6;
	double bestEffX = 0;
	double bestResX = 1e6;
	double bestRMSX = 1e6;
	for (int ibin = 0; ibin<NBINS; ibin++){
		o_ibin = ibin;
		o_xmin = mXmax*(ibin)/NBINS;
		o_xmid = mXmax*(ibin+0.5)/NBINS;
		o_xmax = mXmax*(ibin+1)/NBINS;
		o_nx = N_BIN[ibin];
		o_nxh = h_resX[ibin]->GetEntries();
		o_nd = h_resD[ibin]->GetEntries();
		o_xrms = h_resX[ibin]->GetRMS();
		o_xrmserr = h_resX[ibin]->GetRMSError();
		o_drms = h_resD[ibin]->GetRMS();
		o_drmserr = h_resD[ibin]->GetRMSError();
		if (o_nx&&h_resX[ibin]->Integral()>0){
			o_xeff = (double)o_nxh/o_nx;
			if (o_xmid<0.5)
				doFit(h_resX[ibin],3/4.,1/4.,0);
			else if (o_xmid>7)
				doFit(h_resX[ibin],1/3.,2/3.);
			else
				doFit(h_resX[ibin],1/3.,1/3.);
			o_xres = f_res->GetParameter(2);
			o_xreserr = f_res->GetParError(2);
			o_xoff = f_res->GetParameter(1);
			ibinl = h_resX[ibin]->FindBin(o_xoff-o_xres*3); 
			ibinr = h_resX[ibin]->FindBin(o_xoff+o_xres*3); 
			o_xeff3sig = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
			ibinl = h_resX[ibin]->FindBin(o_xoff-o_xres*5); 
			ibinr = h_resX[ibin]->FindBin(o_xoff+o_xres*5); 
			o_xeff5sig = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
			ibinl = h_resX[ibin]->FindBin(o_xoff-0.5); 
			ibinr = h_resX[ibin]->FindBin(o_xoff+0.5); 
			o_xeff500um = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
			ibinl = h_resX[ibin]->FindBin(o_xoff-1); 
			ibinr = h_resX[ibin]->FindBin(o_xoff+1); 
			o_xeff1mm = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
			if (savehists){
				h_resX[ibin]->Draw();
				canv_bin->SaveAs(Form("resX%d_%s.%d.png",ibin,gasName.Data(),HV));
			}
		}
		else{
			o_xres = 0;
			o_xoff = 0;
			o_xeff = 0;
			o_xeff3sig  = 0;
			o_xeff5sig  = 0;
			o_xeff500um = 0;
			o_xeff1mm   = 0;
		}
		if (h_resD[ibin]->Integral()>0){
			if (o_xmid<0.5)
				doFit(h_resD[ibin],1/4.,3/4.);
			else if (o_xmid>7)
				doFit(h_resD[ibin],2/3.,1/3.);
			else
				doFit(h_resD[ibin],1/3.,1/3.);
			o_dres = f_res->GetParameter(2);
			o_dreserr = f_res->GetParError(2);
			o_doff = f_res->GetParameter(1);
			ibinl = h_resD[ibin]->FindBin(o_doff-o_dres*3); 
			ibinr = h_resD[ibin]->FindBin(o_doff+o_dres*3); 
			o_deff3sig = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
			ibinl = h_resD[ibin]->FindBin(o_doff-o_dres*5); 
			ibinr = h_resD[ibin]->FindBin(o_doff+o_dres*5); 
			o_deff5sig = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
			ibinl = h_resD[ibin]->FindBin(o_doff-0.5); 
			ibinr = h_resD[ibin]->FindBin(o_doff+0.5); 
			o_deff500um = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
			ibinl = h_resD[ibin]->FindBin(o_doff-1); 
			ibinr = h_resD[ibin]->FindBin(o_doff+1); 
			o_deff1mm = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
			if (savehists){
				h_resD[ibin]->Draw();
				canv_bin->SaveAs(Form("resD%d_%s.%d.png",ibin,gasName.Data(),HV));
			}
		}
		else{
			o_dres = 0;
			o_doff = 0;
			o_deff3sig  = 0;
			o_deff5sig  = 0;
			o_deff500um = 0;
			o_deff1mm   = 0;
		}
		otree->Fill();
		if (o_nd>100){
			if (o_xmid<8){
				NusedD++;
				averageResD+=o_dres;
				averageRMSD+=o_drms;
				averageEffD+=o_deff500um;
				if (bestResD>o_dres) bestResD = o_dres;
				if (bestRMSD>o_drms) bestRMSD = o_drms;
				if (bestEffD<o_deff500um) bestEffD = o_deff500um;
			}
			v_dx.push_back(o_xmid);
			v_dxerr.push_back((o_xmax-o_xmin)/2.);
			v_deff500um.push_back(o_deff500um);
			v_dres.push_back(o_dres);
			v_dreserr.push_back(o_dreserr);
			v_drms.push_back(o_drms);
			v_drmserr.push_back(o_drmserr);
			v_doff.push_back(o_doff);
		}
		if (o_nxh>100){
			if (o_xmid<8){
				NusedX++;
				averageResX+=o_xres;
				averageRMSX+=o_xrms;
				averageEffX+=o_xeff500um;
				if (bestResX>o_xres) bestResX = o_xres;
				if (bestRMSX>o_xrms) bestRMSX = o_xrms;
				if (bestEffX<o_xeff500um) bestEffX = o_xeff500um;
			}
			v_xx.push_back(o_xmid);
			v_xxerr.push_back((o_xmax-o_xmin)/2.);
			v_xeff.push_back(o_xeff);
			v_xeff500um.push_back(o_xeff500um);
			v_xres.push_back(o_xres);
			v_xreserr.push_back(o_xreserr);
			v_xrms.push_back(o_xrms);
			v_xrmserr.push_back(o_xrmserr);
			v_xoff.push_back(o_xoff);
		}
	}
	averageEffD/=NusedD;
	averageResD/=NusedD;
	averageRMSD/=NusedD;
	averageEffX/=NusedX;
	averageRMSX/=NusedX;
	printf("=>  %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n","runNo","testLayer","HV","THR","gasID","aaCut","averageEffD","averageResD","averageEffX","averageResX","averageRMSD","averageResX","bestEffD","bestResD","bestEffX","bestResX","bestRMSD","bestRMSX");
	printf("==> %d %d %d %d %d %d %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",runNo,testLayer,HV,THR,gasID,aaCut,averageEffD,averageResD,averageEffX,averageResX,averageRMSD,averageResX,bestEffD,bestResD,bestEffX,bestResX,bestRMSD,bestRMSX);

	//=================================================Draw====================================================
	TCanvas * canv_tracking = new TCanvas("canv_tracking","canv_tracking",1024,768);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	canv_tracking->Divide(2,2);
	canv_tracking->cd(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_nHits->Draw();
	TLine * line_nHits = new TLine(nHitsMax,0,nHitsMax,h_nHits->GetMaximum());
	line_nHits->SetLineColor(kRed);
	line_nHits->Draw("SAME");
	TLatex * text_nHits = new TLatex(nHitsMax,h_nHits->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT1,(double)N_CUT1/N_ALL*100));
	text_nHits->SetTextColor(kRed);
	text_nHits->SetTextSize(0.04);
	text_nHits->Draw("SAME");
	canv_tracking->cd(2);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_DOF->Draw();
	TLine * line_DOF = new TLine(nHitsSmin-4,0,nHitsSmin-4,h_DOF->GetMaximum());
	line_DOF->SetLineColor(kRed);
	line_DOF->Draw("SAME");
	TLatex * text_DOF = new TLatex(nHitsSmin-4,h_DOF->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT2,(double)N_CUT2/N_ALL*100));
	text_DOF->SetTextColor(kRed);
	text_DOF->SetTextSize(0.04);
	text_DOF->Draw("SAME");
	canv_tracking->cd(3);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_chi2->Draw();
	TLine * line_chi2 = new TLine(maxchi2,0,maxchi2,h_chi2->GetMaximum());
	line_chi2->SetLineColor(kRed);
	line_chi2->Draw("SAME");
	TLatex * text_chi2 = new TLatex(maxchi2,h_chi2->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT3,(double)N_CUT3/N_ALL*100));
	text_chi2->SetTextColor(kRed);
	text_chi2->SetTextSize(0.04);
	text_chi2->Draw("SAME");
	canv_tracking->cd(4);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_chi2p->Draw();
	TLine * line_chi2p = new TLine(minchi2p,0,minchi2p,h_chi2p->GetMaximum());
	line_chi2p->SetLineColor(kRed);
	line_chi2p->Draw("SAME");
	TLatex * text_chi2p = new TLatex(minchi2p,h_chi2p->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT3,(double)N_CUT3/N_ALL*100));
	text_chi2p->SetTextColor(kRed);
	text_chi2p->SetTextSize(0.04);
	text_chi2p->Draw("SAME");
	canv_tracking->SaveAs(Form("track_%s.%d.pdf",gasName.Data(),HV));
	canv_tracking->SaveAs(Form("track_%s.%d.png",gasName.Data(),HV));

	TCanvas * canv_DOCA = new TCanvas("canv_DOCA","canv_DOCA",600,800);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	canv_DOCA->Divide(1,2);
	canv_DOCA->cd(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_DOCA->GetYaxis()->SetRangeUser(0,h_DOCA->GetMaximum()*1.1);
	h_DOCA->Draw();
	canv_DOCA->cd(2);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_DOCAb->GetYaxis()->SetRangeUser(0,h_DOCAb->GetMaximum()*1.1);
	h_DOCAb->Draw();
	canv_DOCA->SaveAs(Form("DOCA_%s.%d.pdf",gasName.Data(),HV));
	canv_DOCA->SaveAs(Form("DOCA_%s.%d.png",gasName.Data(),HV));
	canv_DOCA->cd(1);
	h_DriftD->GetYaxis()->SetRangeUser(0,h_DriftD->GetMaximum()*1.1);
	h_DriftD->Draw();
	canv_DOCA->cd(2);
	h_DriftDb->GetYaxis()->SetRangeUser(0,h_DriftDb->GetMaximum()*1.1);
	h_DriftDb->Draw();
	canv_DOCA->SaveAs(Form("DriftD_%s.%d.pdf",gasName.Data(),HV));
	canv_DOCA->SaveAs(Form("DriftD_%s.%d.png",gasName.Data(),HV));

	TCanvas * canv_XT = new TCanvas("canv_XT","canv_XT",800,600);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_xt->Draw("COLZ");
	h_xt->GetXaxis()->SetRangeUser(minDT,maxDT);
	if (xtType%10==0){
		f_left0->Draw("SAME");
		f_right0->Draw("SAME");
	}
	else if (xtType%10==1){
		f_left->Draw("SAME");
		f_right->Draw("SAME");
	}
	else if (xtType%10==3){
		f_left_mid->Draw("SAME");
		f_right_mid->Draw("SAME");
	}
	else{ // fxc_4|fxm_4
		f_left_mid->SetRange(t7l,maxDT);
		f_left_mid->Draw("SAME");
		f_left->SetRange(minDT,t7l);
		f_left->Draw("SAME");
		f_right_mid->SetRange(t7r,maxDT);
		f_right_mid->Draw("SAME");
		f_right->SetRange(minDT,t7r);
		f_right->Draw("SAME");
	}
	canv_XT->SaveAs(Form("XT_%s.%d.pdf",gasName.Data(),HV));
	canv_XT->SaveAs(Form("XT_%s.%d.png",gasName.Data(),HV));

	TCanvas * canv_TX = new TCanvas("canv_TX","canv_TX",600,800);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_tx->Draw("COLZ");
	h_tx->Draw("COLZ");
	h_tx->GetYaxis()->SetRangeUser(minDT,maxDT);
	canv_TX->SaveAs(Form("TX_%s.%d.pdf",gasName.Data(),HV));
	canv_TX->SaveAs(Form("TX_%s.%d.png",gasName.Data(),HV));

	TCanvas * canv_general = new TCanvas("canv_general","canv_general",1024,768);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_resVSX->Draw("COLZ");
	canv_general->SaveAs(Form("resVSX_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("resVSX_%s.%d.png",gasName.Data(),HV));

	h_resVSD->Draw("COLZ");
	canv_general->SaveAs(Form("resVSD_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("resVSD_%s.%d.png",gasName.Data(),HV));

	h_aaVST->Draw("COLZ");
	TLine * line_aaVST = new TLine(mTmin,aaCut,mTmax,aaCut);
	line_aaVST->SetLineColor(kRed);
	line_aaVST->Draw("SAME");
	canv_general->SaveAs(Form("aaVST_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("aaVST_%s.%d.png",gasName.Data(),HV));

	h_aaVSD->Draw("COLZ");
	TLine * line_aaVSD = new TLine(0,aaCut,mXmax,aaCut);
	line_aaVSD->SetLineColor(kRed);
	line_aaVSD->Draw("SAME");
	canv_general->SaveAs(Form("aaVSD_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("aaVSD_%s.%d.png",gasName.Data(),HV));

	h_ggVSX->Draw("COLZ");
	canv_general->SaveAs(Form("ggVSX_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("ggVSX_%s.%d.png",gasName.Data(),HV));

	TLegend* leg_dedx = new TLegend(0.7,0.7,0.9,0.9);
	leg_dedx->AddEntry(h_dedx[0],Form("All hits used"));
	h_dedx[0]->SetLineColor(1);
	h_dedx[0]->GetYaxis()->SetRangeUser(0,h_dedx[MAXTRANC-1]->GetMaximum());
	h_dedx[0]->Draw();
	for (int i = 1; i<MAXTRANC; i++){
		h_dedx[i]->SetLineColor(i+1);
		leg_dedx->AddEntry(h_dedx[i],Form("Neglecting %d hits",i));
		h_dedx[i]->Draw("SAME");
	}
	leg_dedx->Draw();
	canv_general->SaveAs(Form("dedx_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("dedx_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_xeff = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xeff[0]));
	TGraphErrors * g_xeff500 = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xeff500um[0]));
	g_xeff->SetName("gxeff");
	g_xeff500->SetName("gxeff500um");
	g_xeff->SetTitle("Efficiency VS DOCA");
	g_xeff->GetXaxis()->SetTitle("DOCA [mm]");
	g_xeff->GetYaxis()->SetTitle("Efficiency");
	g_xeff->SetMarkerStyle(20);
	g_xeff->SetMarkerColor(kBlack);
	g_xeff->SetLineColor(kBlack);
	g_xeff->GetYaxis()->SetRangeUser(0,1.1);
	g_xeff->Draw("APL");
	g_xeff500->SetMarkerStyle(20);
	g_xeff500->SetMarkerColor(kRed);
	g_xeff500->SetLineColor(kRed);
	g_xeff500->Draw("PLSAME");
	TLegend * leg_xeff = new TLegend(0.1,0.1,0.4,0.4);
	leg_xeff->AddEntry(g_xeff,"Raw efficiency","PL");
	leg_xeff->AddEntry(g_xeff500,"Efficiency with 500 um cut","PL");
	leg_xeff->Draw("SAME");
	canv_general->SaveAs(Form("effx_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("effx_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_deff500 = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_deff500um[0]));
	g_deff500->SetName("gdeff500um");
	g_deff500->SetTitle("Efficiency VS drift distance");
	g_deff500->GetXaxis()->SetTitle("Drift distance [mm]");
	g_deff500->GetYaxis()->SetTitle("Efficiency");
	g_deff500->SetMarkerStyle(20);
	g_deff500->SetMarkerColor(kRed);
	g_deff500->SetLineColor(kRed);
	g_deff500->GetYaxis()->SetRangeUser(0,1.1);
	g_deff500->Draw("APL");
	canv_general->SaveAs(Form("effd_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("effd_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_xres = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xres[0]),&(v_xxerr[0]),&(v_xreserr[0]));
	g_xres->SetName("gxres");
	g_xres->SetTitle("Spatial resolution VS DOCA");
	g_xres->GetXaxis()->SetTitle("Drift distance [mm]");
	g_xres->GetYaxis()->SetTitle("Spatial resolution [mm]");
	g_xres->SetMarkerStyle(20);
	g_xres->SetMarkerColor(kBlack);
	g_xres->SetLineColor(kBlack);
	g_xres->Draw("APL");
	canv_general->SaveAs(Form("resx_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("resx_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_dres = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_dres[0]),&(v_dxerr[0]),&(v_dreserr[0]));
	g_dres->SetName("gdres");
	g_dres->SetTitle("Spatial resolution VS drift distance");
	g_dres->GetXaxis()->SetTitle("Drift distance [mm]");
	g_dres->GetYaxis()->SetTitle("Spatial resolution [mm]");
	g_dres->SetMarkerStyle(20);
	g_dres->SetMarkerColor(kBlack);
	g_dres->SetLineColor(kBlack);
	g_dres->Draw("APL");
	canv_general->SaveAs(Form("resd_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("resd_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_xrms = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xrms[0]),&(v_xxerr[0]),&(v_xrmserr[0]));
	g_xrms->SetName("gxrms");
	g_xrms->SetTitle("Spatial rmsolution VS DOCA");
	g_xrms->GetXaxis()->SetTitle("Drift distance [mm]");
	g_xrms->GetYaxis()->SetTitle("Spatial rmsolution [mm]");
	g_xrms->SetMarkerStyle(20);
	g_xrms->SetMarkerColor(kBlack);
	g_xrms->SetLineColor(kBlack);
	g_xrms->Draw("APL");
	canv_general->SaveAs(Form("rmsx_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("rmsx_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_drms = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_drms[0]),&(v_dxerr[0]),&(v_drmserr[0]));
	g_drms->SetName("gdrms");
	g_drms->SetTitle("Spatial rmsolution VS drift distance");
	g_drms->GetXaxis()->SetTitle("Drift distance [mm]");
	g_drms->GetYaxis()->SetTitle("Spatial rmsolution [mm]");
	g_drms->SetMarkerStyle(20);
	g_drms->SetMarkerColor(kBlack);
	g_drms->SetLineColor(kBlack);
	g_drms->Draw("APL");
	canv_general->SaveAs(Form("rmsd_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("rmsd_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_xoff = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xoff[0]));
	g_xoff->SetName("gxoff");
	g_xoff->SetTitle("Offset VS DOCA");
	g_xoff->GetXaxis()->SetTitle("Drift distance [mm]");
	g_xoff->GetYaxis()->SetTitle("offset [mm]");
	g_xoff->SetMarkerStyle(20);
	g_xoff->SetMarkerColor(kBlack);
	g_xoff->SetLineColor(kBlack);
	g_xoff->Draw("APL");
	canv_general->SaveAs(Form("offx_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("offx_%s.%d.png",gasName.Data(),HV));

	TGraphErrors * g_doff = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_doff[0]));
	g_doff->SetName("gdoff");
	g_doff->SetTitle("Offset VS drift distance");
	g_doff->GetXaxis()->SetTitle("Drift distance [mm]");
	g_doff->GetYaxis()->SetTitle("offset [mm]");
	g_doff->SetMarkerStyle(20);
	g_doff->SetMarkerColor(kBlack);
	g_doff->SetLineColor(kBlack);
	g_doff->Draw("APL");
	canv_general->SaveAs(Form("offd_%s.%d.pdf",gasName.Data(),HV));
	canv_general->SaveAs(Form("offd_%s.%d.png",gasName.Data(),HV));

	//=================================================Save====================================================
	for (int i = 0; i<NBINS; i++){
		h_resD[i]->Write();
		h_resX[i]->Write();
	}
	for (int i = 0; i<MAXTRANC; i++){
		h_dedx[i]->Write();
	}
	h_resVSX->Write();
	h_resVSD->Write();
	h_xt->Write();
	h_tx->Write();
	h_aaVST->Write();
	h_aaVSD->Write();
	h_ggVSX->Write();
	h_nHits->Write();
	h_DOF->Write();
	h_chi2->Write();
	h_chi2p->Write();
	h_DOCA->Write();
	h_DOCAb->Write();
	h_DriftD->Write();
	h_DriftDb->Write();
	g_deff500->Write();
	g_xeff->Write();
	g_xeff500->Write();
	g_dres->Write();
	g_xres->Write();
	g_drms->Write();
	g_xrms->Write();
	g_doff->Write();
	g_xoff->Write();
	otree->Write();
	ofile->Close();

	printf("All events %d\n",N_ALL);
	printf("nHits<=%d: %d (%.1f%)\n",nHitsMax,N_CUT1,(double)N_CUT1/N_ALL*100);
	printf("nHitsS>=%d: %d (%.1f%)\n",nHitsSmin,N_CUT2,(double)N_CUT2/N_ALL*100);
	printf("chi2<%.1f (pvalue>%.1f): %d (%.1f%)\n",maxchi2,minchi2p,N_CUT3,(double)N_CUT3/N_ALL*100);
	printf("DOCA<=%.1f: %d (%.1f%)\n",mXmax,N_CUT4,(double)N_CUT4/N_ALL*100);
    return 0;
}

double ADC2Charge(double adc){ // fC
	double charge = 0;
	if (adc<735.346)
		charge = fADC2ChargeFunction->Eval(adc);
	else
		charge = fADC2ChargeFunction->Eval(735);
	return charge;
}

double getGG(double charge, double slx, double slz){
	double nPairs = CELLH*sqrt(1+slx*slx+slz*slz)*npair_per_cm/10;
	double gg = charge*1e-15/1.6e-19/nPairs; // FIXME: fake gas gain for now
	return gg;
}

int t2d(double t, double & d, bool isRight){
	int status = 0;
	TF1 * f = 0; // body

	if (xtType%10==0){ // fx_0
		if (isRight) f = f_right0;
		else f = f_left0;
	}
	else if (xtType%10==1){ // fxc_4
		if (isRight) f = f_right;
		else f = f_left;
	}
	else if (xtType%10==2){ // fxc_4|fxm_4
		if (isRight){
			if (t>t7r&&t<t8r)
				f = f_right_mid;
			else
				f = f_right;
		}
		else{
			if (t>t7l&&t<t8r)
				f = f_left_mid;
			else
				f = f_left;
		}
	}
	else{ // fxm_4
		if (isRight){
			if (t<t8r)
				f = f_right_mid;
			else
				f = f_right;
		}
		else{
			if (t<t8r)
				f = f_left_mid;
			else
				f = f_left;
		}
	}

	double tRight = f->GetXmax();
	double tmin = f->GetXmin();
	double tmax = 0;
	if (xtType/10==0){
		tmax = tRight;
	}
	else{
		if (isRight) tmax = t8r;
		else tmax = t8l;
	}

	if (t<tmin){
		status = 0;
		d = 0;
	}
	else if (t>tmax){
		status = 0;
		if (t>tRight) d = f->Eval(tRight);
		else d = f->Eval(t);
	}
	else{
		status = 1;
		d = f->Eval(t);
	}

	return status;
}

double findFirstX(TF1 * f, double val, double xmin, double xmax, double delta){
	double theX = 0;
	for (double x = xmin+delta; x<xmax; x+=delta){ // At least two solutions. Scan to find the smallest one
		theX = f->GetX(val,xmin,x);
		if (fabs(theX-x)>delta/10.&&fabs(theX-xmin)>delta/10.){
			break;
		}
	}
	return theX;
}

void doFit(TH1D * h,double leftRatio, double rightRatio, double leftEnd, double rightEnd){
	int bmax = h->GetMaximumBin();
	double max = h->GetBinContent(bmax)*leftRatio;
	int binl = bmax-1;
	for (;binl>=3; binl--){
		double height3bins = h->GetBinContent(binl);
		height3bins+=h->GetBinContent(binl-1);
		height3bins+=h->GetBinContent(binl-2);
		if (height3bins/3<max) break;
	}
	binl-=1;
	max = h->GetBinContent(bmax)*rightRatio;
	int binr = bmax+1;
	for (;binr<=mNbinRes-3; binr++){
		double height3bins = h->GetBinContent(binr);
		height3bins+=h->GetBinContent(binr+1);
		height3bins+=h->GetBinContent(binr+2);
		if (height3bins/3<max) break;
	}
	binr+=1;
	double left = h->GetBinCenter(binl);
	double right = h->GetBinCenter(binr);
	if (left<leftEnd) left = leftEnd;
	if (right>rightEnd) right = rightEnd;
	h->Fit("fres","QG","",left,right);
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [runname] <[testLayer (4)] [xtType: (0), fx_0; 1, fxc_4; 2, fxc_4|fxm_4; 3, fxm_4; 10, fx_0, 8mm; 11, fxc_4, 8mm; 12, fxc_4|fxm_4, 8mm; 13, fxm_4, 8mm] [maxchi2 (1)] [nHitsMax (0)] [aaCut (0)] [savehists: (0)] [debug: 0;...] [iEntryStart (0)] [iEntryStop (0)]>\n",name);
}
