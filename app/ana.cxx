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

#include "header.h"

#define NBINS 25

//===================Chamber Parameter============================
double U = 8; // mm
double chamberHL = 599.17/2; // mm
double chamberHH = 170.05/2; // mm
double chamberCY = 572; // mm
// map for wire position
double  map_x[NLAY][NCEL][2];
double  map_y[NLAY][NCEL][2];
double  map_off[NLAY][NCEL];
bool    map_has[NLAY][NCEL];

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

double get_dist(int lid, int wid, double slx, double inx, double slz, double inz);
double getGG(double aa, double slx, double slz);
void printUsage(char * name);

int main(int argc, char** argv){

	//=================================================Get options========================================================
	if (argc<3){
	    printUsage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
    TString runname = argv[2];
    int testLayer = 4;
    if (argc>=4)
        testLayer = (int)strtol(argv[3],NULL,10);
	int xtType = 6;
    if (argc>=5)
		xtType = (int)strtol(argv[4],NULL,10);
	geoSetup = 0; // 0: normal scintillator; 1: finger scintillator
    if (argc>=6)
        geoSetup = (int)strtol(argv[5],NULL,10);
    double maxchi2 = 1;
    if (argc>=7)
        maxchi2 = (double)strtod(argv[6],NULL);
    int nHitsMax = 0;
    if (argc>=8)
        nHitsMax = (int)strtol(argv[7],NULL,10);
    int debugLevel = 0;
    if (argc>=9)
        debugLevel = (int)strtol(argv[8],NULL,10);
    int iEntryStart = 0;
    if (argc>=10)
        iEntryStart = (int)strtol(argv[9],NULL,10);
    int iEntryStop = 0;
    if (argc>=11)
        iEntryStop = (int)strtol(argv[10],NULL,10);
    int nHitsSmin = 7;
    printf("##############Input %d Parameters##################\n",argc);
    printf("runNo       = %d\n",runNo);
    printf("runname     = \"%s\"\n",runname.Data());
    printf("geoSetup:     %s\n",geoSetup==0?"normal scintillator":"finger scintillator");
    printf("xtType:       %d: %s\n",xtType,xtType==0?"asymmetric":(xtType==1?"symmetric, with offset":(xtType==2?"symmetric, thru 0":(xtType==3?"symmetric with nLHits==0":(xtType==4?"symmetric with smallest chi2a":(xtType==5?"symmetric with smallest chi2":"others?"))))));
    printf("maxchi2     = %.3e\n",maxchi2);
    printf("test layer:   %d\n",testLayer);
    printf("maxNhits    = %d\n",nHitsMax);
    printf("debug       = %d\n",debugLevel);
    printf("Entries:     [%d~%d]\n",iEntryStart,iEntryStop);
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

	//=================================================Get related info========================================================
	// get run info
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
	printf("runNo#%d: %s, %d, %s, %d V, %d mV, %.0f sec\n",runNo,gastype.Data(),runGr,duration.Data(),HV,THR,durationTime);

    //Prepare Maps
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
			map_has[lid][wid] = false;
        }
    }

	// get offset
	if (xtType==1||xtType==6){
		TChain * iChain_off = new TChain("t","t");
		iChain_off->Add(Form("%s/info/offset.%d.%s.root",HOME.Data(),runNo,runname.Data()));
		double i_off_delta;
		int i_off_lid;
		int i_off_wid;
		iChain_off->SetBranchAddress("d",&i_off_delta);
		iChain_off->SetBranchAddress("wid",&i_off_wid);
		iChain_off->SetBranchAddress("lid",&i_off_lid);
		int N = iChain_off->GetEntries();
		for (int i = 0; i<N; i++){
			iChain_off->GetEntry(i);
			if (i_off_lid>=0&&i_off_lid<NLAY&&i_off_wid>=0&&i_off_wid<NCEL)
				map_off[i_off_lid][i_off_wid] = i_off_delta;
		}
	}

    // Get Wire Position
    TFile * TFile_wirepos = new TFile(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),runNo,runname.Data()));
    TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
    int     wp_bid;
    int     wp_ch;
    int     wp_wid;
    int     wp_lid;
    double  wp_xro;
    double  wp_yro;
    double  wp_xhv;
    double  wp_yhv;
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
        if (wp_lid>=0&&wp_lid<NLAY&&wp_wid>=0&&wp_wid<NCEL){
            map_x[wp_lid][wp_wid][0] = wp_xhv;
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_x[wp_lid][wp_wid][1] = wp_xro;
            map_y[wp_lid][wp_wid][1] = wp_yro;
            map_has[wp_lid][wp_wid] = true;
        }
    }
    TFile_wirepos->Close();


	//===================Set scintillator geometry============================
	if (geoSetup==0){
        // normal scintillator
        sciYup = chamberCY+chamberHH+180; // mm
        sciYdown = chamberCY-chamberHH-180; 
        sciHL = 300/2.;
        sciHW = 90/2.;
    }
    else{
        // finger scintillator
        sciYup = chamberCY+chamberHH+250; // mm
        sciYdown = chamberCY-chamberHH-195; 
        sciHL = 33/2.;
        sciHW = 33/2.;
    }
    printf("##############Geometry##################\n");
    printf("sciYup      = %.3e\n",sciYup);
    printf("sciYdown    = %.3e\n",sciYdown);
    printf("sciHL       = %.3e\n",sciHL);
    printf("sciHW       = %.3e\n",sciHW);


	//==============================================Prepare input file & output variables=================================================
    // input file
    int triggerNumber;
    int nHits;
    int nHitsG;
    int npairs;
    int isel;
    int icom;
    double iinx;
    double islx;
    double iinz;
    double islz;
    double chi2x;
    double chi2z;
    double chi2i;
    int nHitsS;
    double inx;
    double slx;
    double inz;
    double slz;
    double chi2;
    double chi2p;
    double chi2a;
	int nHitsSmallAll = 0;
	int nHitsSmallSASD = 0;
	int nSmallSumHits = 0;
	int nShadowedHits = 0;
	int nLateHits = 0;
	int nBoundaryHits = 0;
	int nSmallBoundaryHits = 0;
	// the closest peak to the track in the test layer
	double res = 1e9;
	double theDD = 1e9;
	double theDT = 1e9;
	int has = 0;
	int theWid = -1;
	double theSum = 0;
	double thePeak = 0;
	double theHeight = 0;
	double sum1st = 0;
	double dt1st = 0;
	int theIp = 0;
	int theMpi = 0;
	// the highest hit in this event
	int highBid = 0;
	int highCh = 0;
	int highLid = 0;
	int highWid = 0;
	int highIp = 0;
	double highSum = 0;
	double highAA = 0;
	double highDT = 0;
	// hit list
	std::vector<int> * i_layerID = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_ip = 0;
	std::vector<double> * i_aa = 0;
	std::vector<double> * i_driftT = 0;
	std::vector<double> * i_driftD = 0;
	std::vector<double> * i_fitD = 0;
	//----------------------------------Set input file--------------------------------------------
	TChain * ichain = new TChain("t","t");
	ichain->Add(Form("%s/root/ana_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),testLayer));
	ichain->SetBranchAddress("triggerNumber",&triggerNumber);
	ichain->SetBranchAddress("res",&res);
	ichain->SetBranchAddress("theDD",&theDD);
	ichain->SetBranchAddress("theDT",&theDT);
	ichain->SetBranchAddress("theWid",&theWid);
	ichain->SetBranchAddress("theSum",&theSum);
	ichain->SetBranchAddress("sum1st",&sum1st);
	ichain->SetBranchAddress("dt1st",&dt1st);
	ichain->SetBranchAddress("has",&has);
	ichain->SetBranchAddress("thePeak",&thePeak);
	ichain->SetBranchAddress("theHeight",&theHeight);
	ichain->SetBranchAddress("theIp",&theIp);
	ichain->SetBranchAddress("theMpi",&theMpi);
	ichain->SetBranchAddress("highBid",&highBid);
	ichain->SetBranchAddress("highCh",&highCh);
	ichain->SetBranchAddress("highLid",&highLid);
	ichain->SetBranchAddress("highWid",&highWid);
	ichain->SetBranchAddress("highIp",&highIp);
	ichain->SetBranchAddress("highSum",&highSum);
	ichain->SetBranchAddress("highAA",&highAA);
	ichain->SetBranchAddress("highDT",&highDT);
	ichain->SetBranchAddress("nHitsSmallSASD",&nHitsSmallSASD);
	ichain->SetBranchAddress("nHitsSmallAll",&nHitsSmallAll);
	ichain->SetBranchAddress("nSHits",&nShadowedHits);
	ichain->SetBranchAddress("nLHits",&nLateHits);
	ichain->SetBranchAddress("nSSHits",&nSmallSumHits);
	ichain->SetBranchAddress("nBHits",&nBoundaryHits);
	ichain->SetBranchAddress("nSBHits",&nSmallBoundaryHits);
	ichain->SetBranchAddress("nHits",&nHits);
	ichain->SetBranchAddress("nHitsG",&nHitsG);
	ichain->SetBranchAddress("npairs",&npairs);
	ichain->SetBranchAddress("isel",&isel);
	ichain->SetBranchAddress("icom",&icom);
	ichain->SetBranchAddress("islx",&islx);
	ichain->SetBranchAddress("islz",&islz);
	ichain->SetBranchAddress("iinx",&iinx);
	ichain->SetBranchAddress("iinz",&iinz);
	ichain->SetBranchAddress("chi2x",&chi2x);
	ichain->SetBranchAddress("chi2z",&chi2z);
	ichain->SetBranchAddress("chi2i",&chi2i);
	ichain->SetBranchAddress("nHitsS",&nHitsS);
	ichain->SetBranchAddress("slx",&slx);
	ichain->SetBranchAddress("slz",&slz);
	ichain->SetBranchAddress("inx",&inx);
	ichain->SetBranchAddress("inz",&inz);
	ichain->SetBranchAddress("chi2",&chi2);
	ichain->SetBranchAddress("chi2p",&chi2p);
	ichain->SetBranchAddress("chi2a",&chi2a);
	ichain->SetBranchAddress("layerID",&i_layerID);
	ichain->SetBranchAddress("wireID",&i_wireID);
	ichain->SetBranchAddress("ip",&i_ip);
	ichain->SetBranchAddress("aa",&i_aa);
	ichain->SetBranchAddress("driftT",&i_driftT);
	ichain->SetBranchAddress("driftD",&i_driftD);
	ichain->SetBranchAddress("fitD",&i_fitD);
	double theAA;
	double theFitD;
	double theGG;

	TFile * ofile = new TFile(Form("%s/root/res_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),testLayer),"RECREATE");
	TTree * otree = new TTree("t","t");
	int o_ibin;
	double o_xmin;
	double o_xmax;
	double o_res;
	double o_eff;
	double o_eff3sig;
	double o_eff5sig;
	double o_eff500um;
	double o_eff1mm;
	double o_off;
	otree->Branch("ibin",&o_ibin);
	otree->Branch("xmin",&o_xmin);
	otree->Branch("xmax",&o_xmax);
	otree->Branch("res",&o_res);
	otree->Branch("eff",&o_eff);
	otree->Branch("eff3sig",&o_eff3sig);
	otree->Branch("eff5sig",&o_eff5sig);
	otree->Branch("eff500um",&o_eff500um);
	otree->Branch("eff1mm",&o_eff1mm);
	otree->Branch("off",&o_off);

	//=================================================Prepare histograms, function and counters====================================================
	double mTmin = -25-1/0.96/2; // t range for one x bin
	double mTmax = 800+1/0.96/2;
	int    mNbint = 792+1;
	double mXmax = 10; // x range for one t bin
	int    mNbinx = 256;
	double minchi2p = 1;
	double closestchi2 = 1e9;
	double maxRes = 2;
	TF1 * f_res = new TF1("fres","gaus",-maxRes,maxRes);
	TH1I * h_nHits = new TH1I("hnHits","Number of TDC hits in each event",100,0,100);
	TH1I * h_DOF = new TH1I("hDOF","Number of DOF",5,0,5);
	TH1D * h_chi2 = new TH1D("hchi2","#chi^{2} of fitting",256,0,10);
	TH1D * h_chi2p = new TH1D("hchi2p","p value of fitting",256,0,1);
	TH1D * h_DOCA = new TH1D("hDOCA","DOCA with Left/Right",mNbinx,-mXmax,mXmax);
	TH1D * h_DOCAb = new TH1D("hDOCAb","DOCA without Left/Right",mNbinx/2,0,mXmax);
	TH2D * h_aaVST = new TH2D("haaVST","ADC sum VS driftT",mNbint,mTmin,mTmax,200,-50,550);
	TH2D * h_ggVSX = new TH2D("hggVSX","Gas gain VS DOCA",mNbinx,-mXmax,mXmax,256,0,2e5);
	TH2D * h_xt = new TH2D("hxt","X-T relation",mNbinx,-mXmax,mXmax,mNbint,mTmin,mTmax);
	TH2D * h_tx = new TH2D("htx","T-X relation",mNbint,mTmin,mTmax,mNbinx,-mXmax,mXmax);
	TH2D * h_resVSX = new TH2D("hresVSX","Residual VS DOCA",mNbinx,-mXmax,mXmax,256,-maxRes,maxRes);
	TH1D * h_res[NBINS];
	for (int i = 0; i<NBINS; i++){
		double xmin = mXmax*(i)/NBINS;
		double xmax = mXmax*(i+1)/NBINS;
		h_res[i] = new TH1D(Form("hresX%d",i),Form("Resolution with DOCA in [%.1f,%.1f] mm",xmin,xmax),256,-maxRes,maxRes);
	}

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

		// update minchi2p
		if (fabs(chi2-maxchi2)<fabs(closestchi2-maxchi2)){
			closestchi2 = chi2;
			minchi2p = chi2p;
		}

		// ignore events with bad fitting
		N_ALL++;
		h_nHits->Fill(nHits);
		if (nHitsMax&&nHits>nHitsMax) continue;
		N_CUT1++;
		h_DOF->Fill(nHitsS-4);
		if (nHitsS<nHitsSmin) continue;
		N_CUT2++;
		h_chi2->Fill(chi2);
		h_chi2p->Fill(chi2p);
		if (chi2>maxchi2) continue;
		N_CUT3++;

		// get closest wire
		theFitD = 1e3;
		theWid = 0;
		for (int wid = 0; wid<NCEL; wid++){
			double fitD = get_dist(testLayer,wid,slx,inx,slz,inz);
			if (fabs(theFitD)>fabs(fitD)){
				theWid = wid;
				theFitD = fitD;
			}
		}
		int ibinX = fabs(theFitD)/mXmax*NBINS;
		if (ibinX>=NBINS) continue; // too far away from any cell
		h_DOCA->Fill(theFitD);
		h_DOCAb->Fill(abs(theFitD));
		N_CUT4++;
		N_BIN[ibinX]++;

		// fill hit level histograms
		for (int ihit = 0; ihit<nHits; ihit++){
			if ((*i_layerID)[ihit]==testLayer&&(*i_ip)[ihit]==0){ // only consider first peaks in test layer hits
				theAA = (*i_aa)[ihit];
				theDT = (*i_driftT)[ihit];
				// aa
				h_aaVST->Fill(theDT,theAA); // first peak in all hits
				if ((*i_wireID)[ihit]==theWid){ // signal peak
					theDD = (*i_driftD)[ihit];
					res = fabs(theDD) - fabs(theFitD);
					theGG = getGG(theAA,slx,slz);
					// gas gain
					h_ggVSX->Fill(theFitD,theGG);
					// xt
					h_xt->Fill(theFitD,theDT);
					h_tx->Fill(theDT,theFitD);
					// res VS x
					h_resVSX->Fill(theFitD,res);
					// res
					h_res[ibinX]->Fill(res);
				}
			}
		}

		if (debugLevel>=20) printf("  Good Event! Looping in %d hits\n",nHits);
	}

	//=================================================Get bin by bin information====================================================
	int ibinl = 0;
	int ibinr = 0;
	for (int ibin = 0; ibin<NBINS; ibin++){
		o_ibin = ibin;
		o_xmin = mXmax*(ibin)/NBINS;
		o_xmax = mXmax*(ibin+1)/NBINS;
		if (h_res[ibin]->GetEntries()&&N_BIN[ibin]){
			h_res[ibin]->Fit("fres","QN0G","");
			o_res = f_res->GetParameter(2);
			o_off = f_res->GetParameter(1);
			o_eff = h_res[ibin]->GetEntries()/N_BIN[ibin];
			ibinl = h_res[ibin]->FindBin(o_off-o_res*3); 
			ibinr = h_res[ibin]->FindBin(o_off+o_res*3); 
			o_eff3sig = h_res[ibin]->Integral(ibinl,ibinr)/N_BIN[ibin];
			ibinl = h_res[ibin]->FindBin(o_off-o_res*5); 
			ibinr = h_res[ibin]->FindBin(o_off+o_res*5); 
			o_eff5sig = h_res[ibin]->Integral(ibinl,ibinr)/N_BIN[ibin];
			ibinl = h_res[ibin]->FindBin(o_off-0.5); 
			ibinr = h_res[ibin]->FindBin(o_off+0.5); 
			o_eff500um = h_res[ibin]->Integral(ibinl,ibinr)/N_BIN[ibin];
			ibinl = h_res[ibin]->FindBin(o_off-1); 
			ibinr = h_res[ibin]->FindBin(o_off+1); 
			o_eff1mm = h_res[ibin]->Integral(ibinl,ibinr)/N_BIN[ibin];
		}
		else{
			o_res = 0;
			o_off = 0;
			o_eff = 0;
			o_eff3sig  = 0;
			o_eff5sig  = 0;
			o_eff500um = 0;
			o_eff1mm   = 0;
		}
		otree->Fill();
	}

	//=================================================Draw====================================================
	TCanvas * canv_tracking = new TCanvas("canv_tracking","canv_tracking",1024,768);
	canv_tracking->Divide(2,2);
	canv_tracking->cd(1);
	h_nHits->Draw();
	canv_tracking->cd(2);
	h_DOF->Draw();
	canv_tracking->cd(3);
	h_chi2->Draw();
	canv_tracking->cd(4);
	h_chi2p->Draw();
	canv_tracking->SaveAs(Form("track_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_tracking->SaveAs(Form("track_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	TCanvas * canv_DOCA = new TCanvas("canv_DOCA","canv_DOCA",600,800);
	canv_DOCA->Divide(1,2);
	canv_DOCA->cd(1);
	h_DOCA->Draw();
	canv_DOCA->cd(2);
	h_DOCAb->Draw();
	canv_DOCA->SaveAs(Form("DOCA_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_DOCA->SaveAs(Form("DOCA_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	TCanvas * canv_XT = new TCanvas("canv_XT","canv_XT",600,800);
	h_xt->Draw("COLZ");
	canv_XT->SaveAs(Form("XT_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_XT->SaveAs(Form("XT_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	TCanvas * canv_TX = new TCanvas("canv_TX","canv_TX",800,600);
	h_tx->Draw("COLZ");
	canv_TX->SaveAs(Form("TX_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_TX->SaveAs(Form("TX_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	TCanvas * canv_res = new TCanvas("canv_res","canv_res",1024,768);
	h_resVSX->Draw("COLZ");
	canv_res->SaveAs(Form("res_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_res->SaveAs(Form("res_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	TCanvas * canv_aa = new TCanvas("canv_aa","canv_aa",1024,768);
	h_aaVST->Draw("COLZ");
	canv_aa->SaveAs(Form("aa_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_aa->SaveAs(Form("aa_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	TCanvas * canv_gg = new TCanvas("canv_gg","canv_gg",1024,768);
	h_ggVSX->Draw("COLZ");
	canv_gg->SaveAs(Form("gg_%d.%s.layer%d.pdf",runNo,runname.Data(),testLayer));
	canv_gg->SaveAs(Form("gg_%d.%s.layer%d.png",runNo,runname.Data(),testLayer));

	//=================================================Save====================================================
	for (int i = 0; i<NBINS; i++){
		h_res[i]->Write();
	}
	h_resVSX->Write();
	h_xt->Write();
	h_tx->Write();
	h_aaVST->Write();
	h_ggVSX->Write();
	h_nHits->Write();
	h_DOF->Write();
	h_chi2->Write();
	h_chi2p->Write();
	h_DOCA->Write();
	h_DOCAb->Write();
	otree->Write();
	ofile->Close();

	printf("All events %d\n",N_ALL);
	printf("nHits<=%d: %d (%.1f\%)\n",nHitsMax,N_CUT1,(double)N_CUT1/N_ALL*100);
	printf("nHitsS>=%d: %d (%.1f\%)\n",nHitsSmin,N_CUT2,(double)N_CUT2/N_ALL*100);
	printf("chi2<%.1f (pvalue>%.1f): %d (%.1f\%)\n",maxchi2,minchi2p,N_CUT3,(double)N_CUT3/N_ALL*100);
	printf("DOCA<=%.1f: %d (%.1f\%)\n",mXmax,N_CUT4,(double)N_CUT4/N_ALL*100);
    return 0;
}

double getGG(double aa, double slx, double slz){
	return aa/sqrt(1+slx*slx+slz*slz)*100; // FIXME: fake gas gain for now
}

double get_dist(int lid, int wid, double slx, double inx, double slz, double inz)
{
	if (!map_has[lid][wid]) return 1e3;
	double xdown = inx-slx*(sciYup-sciYdown);
	double zdown = inz-slz*(sciYup-sciYdown);
	vTrackU.SetXYZ(inx,sciYup,inz);
	vTrackD.SetXYZ(xdown,sciYdown,zdown);
	vWireHV.SetXYZ(map_x[lid][wid][0],map_y[lid][wid][0],-chamberHL);
	vWireRO.SetXYZ(map_x[lid][wid][1],map_y[lid][wid][1],chamberHL);
	vTrack = vTrackD-vTrackU;
	vWire = vWireRO-vWireHV;
	vDist = vWireHV-vTrackU;
	vAxis = vWire.Cross(vTrack);
	double value = -vDist*(vAxis.Unit());
	return value-map_off[lid][wid];
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [runname] <[testLayer (4)] [xtType: (6), sym+offset+1st peak; 3, sym with min nLHits, 2, sym, thr 0; 1, sym+offset; 0, no req] [geoSetup: (0), normal;1, finger] [maxchi2 (1)] [nHitsMax (0)] [debug: 0;...] [iEntryStart (0)] [iEntryStop (0)]>\n",name);
}
