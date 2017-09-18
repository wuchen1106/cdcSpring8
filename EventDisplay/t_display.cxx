#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TF1.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>

//#define PRINT_CROSSPOINTS

#define XMAX 130
#define ZMAX 350

#define NLAY 9
#define NZXP 8
#define NCEL 11

#define NBRD 2
#define NCHS 48
#define NCHT 96
#define NSMP 32
#define ADCMIN 125
#define ADCMAX 700

//===================About xt============================
double tres = 0;
TF1 * f_left_end[88];
TF1 * f_right_end[88];
TF1 * f_left[88];
TF1 * f_right[88];
int vtrel[88];
int vtrer[88];
int vtlel[88];
int vtler[88];
int vtrml[88];
int vtrmr[88];
int vtlml[88];
int vtlmr[88];

double t2x(double time, int lid, int wid, int lr, int & status);

int main(int argc, char** argv){

	//===================Get Arguments============================
	if (argc<3) return -1;
	int runNo = (int)strtol(argv[1],NULL,10);
	int workMode = (int)strtol(argv[2],NULL,10); // 0: h_XXX; 1: i_XXX; 2: t_XXX
	int thelayer = 5;
	if (argc>=4){
	    thelayer = atoi(argv[3]);
    }
	int thewire = -1; // negative value means just pick up the first hit, regardless of the target layer.
	if (argc>=5){
	    thewire = atoi(argv[4]);
    }
	TString suffix = "";
	if (argc>=6){
		suffix  = argv[5];
	}
    int iEntryStart = 0;
    int iEntryStop = 9;
	if (argc>=8){
	    iEntryStart = (int)strtol(argv[6],NULL,10);
	    iEntryStop = (int)strtol(argv[7],NULL,10);
    }
    else if (argc>=7){
	    iEntryStart=0;
	    iEntryStop=(int)strtol(argv[6],NULL,10)-1;
	}
	printf("runNo:              %d\n",runNo);
	printf("workMode:           %d\n",workMode);
	printf("test wire:          [%d,%d]\n",thelayer,thewire);
	printf("suffix:             %s\n",suffix.Data());
	printf("Entries:            %d~%d\n",iEntryStart,iEntryStop);

	//===================Chamber Parameter============================
	double ytop = 642.;         // top of the chamber
	double yup = 623.97007;     // up plane used as the incident plane
	double ydown = 527.60011;   // bottom of the chamber
	double chamberHL = 599.17/2.;   // half of the chamber length;

	//===================Prepare Maps============================
	// map for wire position
    double  map_x[NLAY][NCEL][2];
    double  map_y[NLAY][NCEL][2];
    double  map_z[NLAY][NCEL][2];
	int     map_ch[NLAY][NCEL];
	int     map_bid[NLAY][NCEL];
	int     map_lid[NBRD][NCHS];
	int     map_wid[NBRD][NCHS];
    int     map_check[NLAY][NCEL];
	// mcp for cross points
    double mcp_xc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
    double mcp_zc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
	for(int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
		    map_x[lid][wid][0] = 0;
		    map_y[lid][wid][0] = 0;
		    map_z[lid][wid][0] = 0;
		    map_x[lid][wid][1] = 0;
		    map_y[lid][wid][1] = 0;
		    map_z[lid][wid][1] = 0;
			map_check[lid][wid]=0;
			if (lid <NZXP){ // z-x planes corresponding to the layerID of the lower layer counting from 1 
                for (int wjd = 0; wjd<NCEL; wjd++){
                    mcp_xc[lid][wid][wjd] = 999;
                    mcp_zc[lid][wid][wjd] = 999;
                }
            }
		}
	}

	//===================Get Wire Position============================
	TFile * TFile_wirepos = new TFile("../info/wire-position.root");
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
    std::vector<double> v_wire_xc; // to get a graph of wire centers
    std::vector<double> v_wire_yc; // to get a graph of wire centers
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=0&&wp_lid<NLAY&&wp_wid>=0&&wp_wid<NCEL){
            map_x[wp_lid][wp_wid][0] = wp_xhv;
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_z[wp_lid][wp_wid][0] = -chamberHL;
            map_x[wp_lid][wp_wid][1] = wp_xro;
            map_y[wp_lid][wp_wid][1] = wp_yro;
            map_z[wp_lid][wp_wid][1] = chamberHL;
            map_ch[wp_lid][wp_wid] = wp_ch;
            map_bid[wp_lid][wp_wid] = wp_bid*48;
            if (wp_lid>=1){ // ignore the first layer (dummy)
                map_check[wp_lid][wp_wid] = 1;
                v_wire_xc.push_back((wp_xhv+wp_xro)/2.);
                v_wire_yc.push_back((wp_yhv+wp_yro)/2.);
            }
		}
        else{
            fprintf(stderr,"ERROR: Entry %d in wiremap file, lid = %d wid = %d out of range (%d,%d)!\n",i,wp_lid,wp_wid,NLAY,NCEL);
            return -1;
        }
		if (wp_bid>=0&&wp_bid<NBRD&&wp_ch>=0&&wp_ch<NCHS){
			map_lid[wp_bid][wp_ch] = wp_lid;
			map_wid[wp_bid][wp_ch] = wp_wid;
        }
        else{
            fprintf(stderr,"ERROR: Entry %d in wiremap file, bid = %d ch = %d out of range (%d,%d)!\n",i,wp_bid,wp_ch,NBRD,NCHS);
            return -1;
        }
	}
    TTree_wirepos->SetMarkerStyle(20);
    TTree_wirepos->SetMarkerSize(0.5);
	TFile_wirepos->Close();

	//==================Get Crosspoints==========================
    TFile * TFile_crosspoint = new TFile("../info/crosspoint.root");
    TTree * TTree_crosspoint = (TTree*) TFile_crosspoint->Get("t");
    int     cp_l1;
	int     cp_l2;
	int     cp_w1;
	int     cp_w2;
    double  cp_zc;
	double  cp_xc;
    TTree_crosspoint->SetBranchAddress("l1",&cp_l1);
    TTree_crosspoint->SetBranchAddress("l2",&cp_l2);
    TTree_crosspoint->SetBranchAddress("w1",&cp_w1);
    TTree_crosspoint->SetBranchAddress("w2",&cp_w2);
    TTree_crosspoint->SetBranchAddress("z",&cp_zc);
    TTree_crosspoint->SetBranchAddress("x",&cp_xc);
    int nEntries_crosspoint = TTree_crosspoint->GetEntries();
    for (int iEntry = 0; iEntry<nEntries_crosspoint; iEntry++){
        TTree_crosspoint->GetEntry(iEntry);
        mcp_xc[cp_l1][cp_w1][cp_w2] = cp_xc;
        mcp_zc[cp_l1][cp_w1][cp_w2] = cp_zc;
    }
    TFile_crosspoint->Close();

    //===================Get XT============================
	TFile * i_xt = new TFile(Form("../info/xt.%d.root",runNo));
	for (int i = 0; i<88; i++){
		f_left_end[i] = (TF1*) i_xt->Get(Form("f_left_end_%d_%d",i/11+1,i%11));
		f_right_end[i] = (TF1*) i_xt->Get(Form("f_right_end_%d_%d",i/11+1,i%11));
		f_left[i] = (TF1*) i_xt->Get(Form("f_left_%d_%d",i/11+1,i%11));
		f_right[i] = (TF1*) i_xt->Get(Form("f_right_%d_%d",i/11+1,i%11));
	}
	TTree * itree_xt = (TTree*) i_xt->Get("p");
	int trel,trer,tlel,tler,tlml,tlmr,trml,trmr;
	int lid,wid;
	itree_xt->SetBranchAddress("lid",&lid);
	itree_xt->SetBranchAddress("wid",&wid);
	itree_xt->SetBranchAddress("tlel",&tlel);
	itree_xt->SetBranchAddress("tler",&tler);
	itree_xt->SetBranchAddress("tlml",&tlml);
	itree_xt->SetBranchAddress("tlmr",&tlmr);
	itree_xt->SetBranchAddress("trml",&trml);
	itree_xt->SetBranchAddress("trmr",&trmr);
	itree_xt->SetBranchAddress("trel",&trel);
	itree_xt->SetBranchAddress("trer",&trer);
	for (int i = 0; i<itree_xt->GetEntries(); i++){
		itree_xt->GetEntry(i);
		vtlel[(lid-1)*11+wid] = tlel;
		vtler[(lid-1)*11+wid] = tler;
		vtlml[(lid-1)*11+wid] = tlml;
		vtlmr[(lid-1)*11+wid] = tlmr;
		vtrel[(lid-1)*11+wid] = trel;
		vtrer[(lid-1)*11+wid] = trer;
		vtrml[(lid-1)*11+wid] = trml;
		vtrmr[(lid-1)*11+wid] = trmr;
	}

	//==================Get ADC==========================
	TChain * iChain_ADC = new TChain("tree","tree");
	iChain_ADC->Add(Form("../root/run_%0.6d_built.root",runNo));
	int adc[NCHT][NSMP];
	int tdc[NCHT][NSMP];
	int clockNumberDriftTime[NCHT][NSMP];
	int tdcNhit[NCHT];
	int triggerNumber_ADC;
	iChain_ADC->SetBranchAddress("adc",adc);
	iChain_ADC->SetBranchAddress("driftTime",tdc);
	iChain_ADC->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	iChain_ADC->SetBranchAddress("tdcNhit",tdcNhit);
	iChain_ADC->SetBranchAddress("triggerNumber",&triggerNumber_ADC);

	//===================Get ROOT File============================
	int triggerNumber;
	int nHits;
	double xup;
	double zup;
	double slx;
	double slz;
	double xupi;
	double zupi;
	double slix;
	double sliz;
	double chi2;
    std::vector<double> * i_fitD = 0;
	std::vector<double> * i_driftD = 0;
	std::vector<double> * i_driftT = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_layerID = 0;
	// file from hit picking
	TChain * iChain_h = new TChain("t","t");
	// file from track finding
	// FIXME: TBA
	// file from track fitting
	TChain * iChain_t = new TChain("t","t");
    // decide the main file to get hits
    TChain * iChain = 0;
    if (workMode==0){ // 0: h_XXX; 1: i_XXX; 2: t_XXX
        iChain_h->Add(Form("../root/h_%d.",runNo)+suffix+".root");
        iChain = iChain_h;
    }
    else if (workMode==1){
        fprintf(stderr,"workMode 1 (i_XXX) is not supported yet!\n");
        return -1;
    }
    else if (workMode==2){
        iChain_t->Add(Form("../root/t_%d.layer%d.",runNo,thelayer)+suffix+".root");
        iChain_t->SetBranchAddress("inX",&xup);
        iChain_t->SetBranchAddress("inZ",&zup);
        iChain_t->SetBranchAddress("slX",&slx);
        iChain_t->SetBranchAddress("slZ",&slz);
        iChain_t->SetBranchAddress("iniX",&xupi);
        iChain_t->SetBranchAddress("iniZ",&zupi);
        iChain_t->SetBranchAddress("sliX",&slix);
        iChain_t->SetBranchAddress("sliZ",&sliz);
        iChain_t->SetBranchAddress("chi2",&chi2);
        iChain_t->SetBranchAddress("fitD",&i_fitD);
        iChain_t->SetBranchAddress("driftD",&i_driftD);
        iChain = iChain_t;
    }
    else{
        fprintf(stderr,"workMode %d is not supported! please chose from 0,1,2\n",workMode);
        return -1;
    }
	iChain->SetBranchAddress("triggerNumber",&triggerNumber);
	iChain->SetBranchAddress("nHits",&nHits);
	iChain->SetBranchAddress("wireID",&i_wireID);
	iChain->SetBranchAddress("layerID",&i_layerID);
	iChain->SetBranchAddress("driftT",&i_driftT);

	//==================Prepare for drawing==========================
	//Prepare the Canvas for waveforms
	TLatex * canvtitle = new TLatex(0,0,"");
	canvtitle->SetTextSize(0.02);
	TCanvas * ca_WF[NBRD];
	TPad * pa_WF[NBRD][NCHS];
	for (int bid = 0; bid<NBRD; bid++){
        ca_WF[bid] = new TCanvas(Form("ca_WF_%d",bid),"ca_WF",896,896);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        for (int i = 0; i<8; i++){
            for (int j = 0; j<6; j++){
                int index = j*8+i;
                pa_WF[bid][index] = new TPad(Form("pa_%d_%d_%d",bid,i,j),Form("pa_%d_%d_%d",bid,i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
                pa_WF[bid][index]->Draw();
                pa_WF[bid][index]->SetGridx(1);
                pa_WF[bid][index]->SetGridy(1);
            }
        }
    }
	//Prepare the Canvas for x-y plane and target chanel ADC
	TCanvas * ca_xyADC = new TCanvas("ca_xyADC","ca_xyADC",896,1024);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
	TPad * p_xyADC[2];
	for (int ipad = 0; ipad<2; ipad++){
        p_xyADC[ipad] = new TPad(Form("p_xyADC_%i",ipad),"p_xyADC",0,ipad?0:0.2,1,ipad?0.2:1);
        p_xyADC[ipad]->Draw();
        p_xyADC[ipad]->SetGridx(1);
        p_xyADC[ipad]->SetGridy(1);
	}
	//Prepare the Canvas for z-x planes
	TCanvas* ca_zx[NZXP]; // z-x planes corresponding to the layerID of the lower layer counting from 1
	for (int lid = 1; lid<NZXP; lid++){
	    ca_zx[lid] = new TCanvas(Form("ca_zx_%d",lid),"ca_zx",1024,768);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
    }
	//Prepare for ADC
	TGraph * gr_waveForm[NBRD][NCHS];
	int vSample[NSMP];
	for (int i=0; i<NSMP; i++){
		vSample[i] = i;
	}
	TLatex *textTDC[NBRD][NCHS][NSMP];
	TMarker *markerTDC[NBRD][NCHS][NSMP];
    for (int bid = 0; bid<NBRD; bid++){
        for (int ch = 0; ch<NCHS; ch++){
            for (int sid=0; sid<NSMP; sid++) {
                textTDC[bid][ch][sid] = new TLatex(0,0,"");
                textTDC[bid][ch][sid]->SetTextSize(0.04);
                textTDC[bid][ch][sid]->SetTextColor(kRed);
                markerTDC[bid][ch][sid] = new TMarker(0,0,20);
                markerTDC[bid][ch][sid]->SetMarkerSize(0.3);
                markerTDC[bid][ch][sid]->SetMarkerColor(kRed);
            }
        } 
	}
	//Prepare for hit circles on x-y plane
	TText * text[NLAY][NCEL];
	TEllipse * ewiret[NLAY][NCEL];  // value from track finding/fitting;
	TEllipse * ewiret2[NLAY][NCEL]; // initial value for track fitting;
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			ewiret[lid][wid] = new TEllipse(0,0,1,1);
			ewiret2[lid][wid] = new TEllipse(0,0,1,1);
		}
	}
	//Prepare two tracks on x-y plane: l: track after fitting; l2: track before fitting, after track finding.
	TLine * l = new TLine();
	l->SetLineColor(kRed);
	l->SetY1(ytop);
	l->SetY2(ydown);
	TLine * l2 = new TLine();
	l2->SetLineColor(kBlue);
	l2->SetY1(ytop);
	l2->SetY2(ydown);
    // Background histogram for x-y plane
    TGraph * gr_wireCenter = new TGraph(v_wire_yc.size(),&(v_wire_xc[0]),&(v_wire_yc[0]));
    gr_wireCenter->SetMarkerStyle(20);
    gr_wireCenter->SetMarkerSize(0.55);
    gr_wireCenter->GetXaxis()->SetTitle("x [mm]");
    gr_wireCenter->GetYaxis()->SetTitle("y [mm]");
    // Prepare colors to be used for each wire
    int color[NCEL];
    int icolor = 0;
    color[icolor++] = kBlack;
    color[icolor++] = kMagenta+2;
    color[icolor++] = kMagenta;
    color[icolor++] = kBlue;
    color[icolor++] = kBlue-2;
    color[icolor++] = kCyan;
    color[icolor++] = kGreen;
    color[icolor++] = kYellow;
    color[icolor++] = kOrange;
    color[icolor++] = kRed;
    color[icolor++] = kRed+2;
    // Prepare driftT lines on z-x planes
    TLine * l_zx[NLAY][NCEL][4];
    for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
            for (int ilr = 0; ilr<2; ilr++){
                l_zx[lid][wid][ilr] = new TLine();
                l_zx[lid][wid][ilr]->SetLineColor(color[wid]);
                if (workMode>=1&&ilr)
                    l_zx[lid][wid][ilr]->SetLineStyle(2);
                else
                    l_zx[lid][wid][ilr]->SetLineStyle(3);
            }
        }
    }
    // Prepare cross points of driftT lines on z-x planes
    TMarker * point_cross_zx[NZXP][NCEL][NCEL][4];
    for (int lid = 1; lid<NZXP; lid++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
        for (int wid = 0; wid<NCEL; wid++){
            for (int wjd = 0; wjd<NCEL; wjd++){
                for (int icombi = 0; icombi<4; icombi++){
                    if (workMode>=1&&icombi==3) // reverse lid when icombi is 0 or 1, reverse lid+1 when icombi is 0 or 2. Keep them unchanged only when icombi is 3
                        point_cross_zx[lid][wid][wjd][icombi] = new TMarker(0,0,20);
                    else
                        point_cross_zx[lid][wid][wjd][icombi] = new TMarker(0,0,4);
                    point_cross_zx[lid][wid][wjd][icombi]->SetMarkerColor(kBlack);
                    point_cross_zx[lid][wid][wjd][icombi]->SetMarkerSize(0.55);
                }
            }
        }
    }
    // Prepare track points on z-x planes
    TMarker * point_track_zx[NZXP];
    if (workMode>=1){
        for (int lid = 1; lid<NZXP; lid++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
            point_track_zx[lid] = new TMarker(0,0,20);
            point_track_zx[lid]->SetMarkerColor(kRed);
            point_track_zx[lid]->SetMarkerSize(0.5);
        }
    }
    // Prepare wires in each layer
    //  and to record the position of hits in each layer
    int    nHits_zx[NLAY];
    bool   check_zx[NLAY][NCEL];
    double y_zx[NLAY][NCEL];
    double dd_zx[NLAY][NCEL];
    double fd_zx[NLAY][NCEL];
    double theta_zx[NLAY][NCEL];
    TGraph * gr_wire[NLAY][NCEL];
    for (int lid = 1; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            if (map_check[lid][wid]==1){
                gr_wire[lid][wid] = new TGraph(2,map_z[lid][wid],map_x[lid][wid]);
                gr_wire[lid][wid]->SetLineColor(color[wid]);
                gr_wire[lid][wid]->SetMarkerColor(color[wid]);
                y_zx[lid][wid] = (map_y[lid][wid][1]+map_y[lid][wid][0])/2.;
            }
        }
    }
    // Prepare texts and markers for wire on each z-x planes
	TMarker * point_cross_wire[NZXP][NCEL][NCEL];
	TLatex * text_cross_wire[NZXP][NCEL][NCEL];
	TLatex * text_cr_1[NZXP];
	TLatex * text_cr_1l[NZXP][NCEL];
	TLatex * text_cr_2[NZXP];
	TLatex * text_cr_2r[NZXP][NCEL];
    TGraph * gr_all[NZXP];
    double ar_cr_x[2] = {-ZMAX,ZMAX};
    double ar_cr_y[2] = {-XMAX,XMAX};
    for (int lid = 1; lid<NZXP; lid++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
        // Background graph for z-x planes
        gr_all[lid] = new TGraph(2,ar_cr_x,ar_cr_y);
        gr_all[lid]->SetTitle(Form("layer #%d and layer #%d",lid,lid+1));
        gr_all[lid]->SetMarkerColor(kWhite);
        gr_all[lid]->SetMarkerSize(0.1);
        gr_all[lid]->GetXaxis()->SetRangeUser(-ZMAX,ZMAX);
        gr_all[lid]->GetYaxis()->SetRangeUser(-XMAX,XMAX);
        gr_all[lid]->GetXaxis()->SetTitle("z [mm]");
        gr_all[lid]->GetYaxis()->SetTitle("x [mm]");
        text_cr_1[lid] = new TLatex(map_z[lid][0][0]-30,XMAX-20,Form("Layer %d",lid));
        text_cr_1[lid]->SetTextSize(0.025);
        text_cr_2[lid] = new TLatex(map_z[lid+1][0][1]-30,XMAX-20,Form("Layer %d",lid+1));
        text_cr_2[lid]->SetTextSize(0.025);
        for (int wid = 0; wid<NCEL; wid++){
            if (map_check[lid][wid]==1){
                text_cr_1l[lid][wid] = new TLatex(map_z[lid][wid][0]-10,map_x[lid][wid][0],Form("%d",wid));
                text_cr_1l[lid][wid]->SetTextColor(color[wid]);
                text_cr_1l[lid][wid]->SetTextSize(0.02);
            }
            if (map_check[lid+1][wid]==1){
                text_cr_2r[lid][wid] = new TLatex(map_z[lid+1][wid][1]+10,map_x[lid+1][wid][1],Form("%d",wid));
                text_cr_2r[lid][wid]->SetTextColor(color[wid]);
                text_cr_2r[lid][wid]->SetTextSize(0.02);
            }
#ifdef PRINT_CROSSPOINTS
            for (int wjd = 0; wjd < NCEL; wjd++){
                if (fabs(mcp_zc[lid][wid][wjd])>300) continue;
                point_cross_wire[lid][wid][wjd] = new TMarker(mcp_zc[lid][wid][wjd],mcp_xc[lid][wid][wjd],20);
                point_cross_wire[lid][wid][wjd]->SetMarkerColor(color[wid]);
                point_cross_wire[lid][wid][wjd]->SetMarkerSize(0.4);
                text_cross_wire[lid][wid][wjd] = new TLatex(mcp_zc[lid][wid][wjd],mcp_xc[lid][wid][wjd]+5,Form("%d,%d",wid,wjd));
                text_cross_wire[lid][wid][wjd]->SetTextColor(color[wid]);
                text_cross_wire[lid][wid][wjd]->SetTextSize(0.02);
            }
#endif
        }
    }

	//===================Loop in Events============================
	for ( int iEntry = iEntryStart; iEntry<=iEntryStop; iEntry++){
		if ((iEntry-iEntryStart)%1000==0) printf("%lf\n",iEntry);
		iChain->GetEntry(iEntry);
		if (nHits<=0) continue;

        // Find the target channel
		int the_ihit = 0;
		int the_bid = -1;
		int the_ch = -1;
        if (thelayer<0||thewire<0){ // negative value means just pick up the first hit, regardless of the target layer
            int lid = (*i_layerID)[0];
            int wid = (*i_wireID)[0];
            the_bid = map_bid[lid][wid];
            the_ch = map_ch[lid][wid];
        }
        else{
            for (; the_ihit<i_driftT->size(); the_ihit++){
                int lid = (*i_layerID)[the_ihit];
                int wid = (*i_wireID)[the_ihit];
                if (lid==thelayer&&wid==thewire){
                    the_ch = map_ch[lid][wid];
                    the_bid = map_bid[lid][wid];
                    break;
                }
            }
        }
		if (the_bid==-1||the_ch==-1) continue; // only show the events with hit in the target cell

        // Get ADC information
		for (int iEntry2 = triggerNumber; iEntry2>=0 ;iEntry2--){ // in case that the raw ROOT file missed some triggerNumbers
			iChain_ADC->GetEntry(iEntry2);
			if (triggerNumber==triggerNumber_ADC) break;
		}

//		if ((*i_fitD)[the_ihit]>-2.3||(*i_fitD)[the_ihit]<-3.2||fabs((*i_driftD)[the_ihit]-(*i_fitD)[the_ihit])<0.5||chi2>1||fabs(slz)>0.15) continue;
//		if ((*i_fitD)[the_ihit]>0.32||(*i_fitD)[the_ihit]<0.28||chi2>3) continue;
//		if ((*i_driftT)[the_ihit]>80||chi2>3) continue;

        // clear objects from the previous event
	    for (int bid = 0; bid<NBRD; bid++){
	        for (int ch = 0; ch<NCHS; ch++){
                if (gr_waveForm[bid][ch]) delete gr_waveForm[bid][ch]; gr_waveForm[bid][ch] = 0;
	        }
	    }

        // Reset counters
		for (int lid = 0; lid<NLAY; lid++){
		    nHits_zx[lid] = 0;
			for (int wid = 0; wid<NCEL; wid++){
                check_zx[lid][wid] = false;
			}
		}

        // Draw waveforms
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->cd();
            canvtitle->SetText(0.1,0.98,Form("Entry %d, Trigger Number %d, Board %d",iEntry,triggerNumber,bid));
            canvtitle->Draw();
            for (int ch = 0; ch<NCHS; ch++){
                int chg = bid*NCHS+ch;
                pa_WF[bid][ch]->cd();
                gr_waveForm[bid][ch] = new TGraph(NSMP,vSample,adc[chg]);
                gr_waveForm[bid][ch]->SetTitle(Form("Channel %d Layer %d Wire %d",ch,map_lid[bid][ch],map_wid[bid][ch]));
                gr_waveForm[bid][ch]->GetXaxis()->SetRangeUser(0,NSMP-1);
                gr_waveForm[bid][ch]->GetYaxis()->SetRangeUser(ADCMIN,ADCMAX);
                gr_waveForm[bid][ch]->GetXaxis()->SetTitle("Sample Index");
                gr_waveForm[bid][ch]->GetYaxis()->SetTitle("ADC");
                gr_waveForm[bid][ch]->SetMarkerStyle(20);
                gr_waveForm[bid][ch]->SetMarkerSize(0.3);
                gr_waveForm[bid][ch]->Draw("APL");
                for (int ihit=0; ihit<tdcNhit[chg]; ihit++) {
                    int clk = clockNumberDriftTime[chg][ihit];
                    markerTDC[bid][ch][ihit]->SetX(clk);
                    markerTDC[bid][ch][ihit]->SetY(adc[chg][clk]);
                    //if (adc[chg][clk]>Hmin[chg]) markerTDC[ihit]->SetMarkerColor(kRed);
                    //else markerTDC[ihit]->SetMarkerColor(kBlue);
                    markerTDC[bid][ch][ihit]->Draw();
                    textTDC[bid][ch][ihit]->SetText(clk,adc[chg][clk]+(0.5-ihit%2)*50,Form("%d",(int)(tdc[chg][ihit])));
                    textTDC[bid][ch][ihit]->Draw();
                }
            }
        }

        // draw the target channel ADC
		p_xyADC[1]->cd();
		gr_waveForm[the_bid][the_ch]->Draw("APL");
		int the_chg = the_bid*NCHS+the_ch;
		for (int ihit=0; ihit<tdcNhit[the_chg]; ihit++) {
			textTDC[the_bid][the_ch][ihit]->Draw();
			markerTDC[the_bid][the_ch][ihit]->Draw();
		}

        // get the hits and track information
        int nHitsAll = 0;
        double zdown  = 0; 
        double xdown  = 0; 
        double xtop   = 0; 
        double zdowni = 0; 
        double xdowni = 0; 
        double xtopi  = 0; 
		for (int ihit = 0; ihit<i_driftT->size(); ihit++){
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			if (lid<=0) continue; // ignore the first layer (dummy)
			nHitsAll++;

            // Get hit information
			double fitd = 0;
			double dd = 0;
			double dt = (*i_driftT)[ihit];
			double wxro = map_x[lid][wid][1];
			double wyro = map_y[lid][wid][1];
            double wzro = chamberHL;
			double wxhv = map_x[lid][wid][0];
			double wyhv = map_y[lid][wid][0];
            double wzhv = -chamberHL;
			double wy = (wyro+wyhv)/2.;
			double wz = 0;
			double wx = 0;
			if (workMode>=1){
                fitd = (*i_fitD)[ihit]; // !!! make sure it's mm
                dd = (*i_driftD)[ihit]; // !!! make sure it's mm
			    // correct wx wy wz according to the track position
                zdown = zup + (ydown-yup)*slz;
                xdown = xup + (ydown-yup)*slx;
                xtop = xup + (ytop-yup)*slx;
                if (workMode>=2){
                    zdowni = zupi + (ydown-yup)*sliz;
                    xdowni = xupi + (ydown-yup)*slix;
                    xtopi = xupi + (ytop-yup)*slix;
                }
                wz = ((yup-wy)*zdown+(wy-ydown)*zup)/(yup-ydown);
                wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
            }
            else{
                // choose the center point of the wire
                wz = 0;
                wx = (wxhv+wxro)/2.;
                int status = 0;
                dd = t2x(dt,lid,wid,0,status);// FIXME: should get x-t curve
            }

            // Get information for cross points
			double delta = dd*chamberHL*2/sqrt(chamberHL*chamberHL*4+(map_x[lid][wid][0]-map_x[lid][wid][1])*(map_x[lid][wid][0]-map_x[lid][wid][1]));
			for (int ilr = 0; ilr<2; ilr++){
                l_zx[lid][wid][ilr]->SetX1(-chamberHL);
                l_zx[lid][wid][ilr]->SetX2(chamberHL);
                l_zx[lid][wid][ilr]->SetY1(wxhv+(ilr?delta:-delta));
                l_zx[lid][wid][ilr]->SetY2(wxro+(ilr?delta:-delta));
            }
            nHits_zx[lid]++;
            check_zx[lid][wid] = true;
			y_zx[lid][wid] = wy;
            dd_zx[lid][wid] = dd;
            fd_zx[lid][wid] = fitd;
            theta_zx[lid][wid] = atan(-(map_x[lid][wid][0]-map_x[lid][wid][1])/chamberHL/2); // rotation angle w.r.t the dart plane: read out  plane, i.e. z>0, i.e. iplane = 1; positive rotation angle point to -x direction

            // set the hit circles
            ewiret[lid][wid]->SetX1(wx);
            ewiret[lid][wid]->SetY1(wy);
            ewiret[lid][wid]->SetR1(dd);
            ewiret[lid][wid]->SetR2(dd);
			ewiret[lid][wid]->SetFillStyle(0);
			ewiret[lid][wid]->SetLineColor(kRed);
            // set the text of each hit
			if (lid==thelayer&&wid==thewire) text[lid][wid]= new TText(wx,wy,Form("%d,%d,%.1lf,%.1lf",lid,wid,fitd,dt));
			else text[lid][wid]= new TText(wx,wy,Form("%d,%d",lid,wid));
			text[lid][wid]->SetTextSize(0.02);

            if (workMode>=2){
                // do the correction according to the input track parameters from track finding
                wz = ((yup-wy)*zdowni+(wy-ydown)*zupi)/(yup-ydown);
                wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
                ewiret2[lid][wid]->SetX1(wx);
                ewiret2[lid][wid]->SetY1(wy);
                ewiret2[lid][wid]->SetR1(dd);
                ewiret2[lid][wid]->SetR2(dd);
                ewiret2[lid][wid]->SetFillStyle(0);
                ewiret2[lid][wid]->SetLineColor(kBlue);
            }
		}

        // Draw the background histogram for x-y plane
		p_xyADC[0]->cd();
		gr_wireCenter->SetTitle(Form("Entry %d, Trigger Number %d, nHits = %d, chi2 = %.2f",iEntry,triggerNumber,nHitsAll,chi2));
		gr_wireCenter->Draw("AP");
        // Draw the hit circles
		for (int ihit = 0; ihit<i_driftT->size(); ihit++){
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			if (ewiret[lid][wid]&&text[lid][wid]){
                ewiret[lid][wid]->Draw(); // from track fitting/finding
                text[lid][wid]->Draw();
            }
            if (workMode>=2&&ewiret2[lid][wid]){
                ewiret2[lid][wid]->Draw(); // input values for track fitting.
            }
		}

        // draw the tracks on the x-y plane
        if (workMode>=1){ // from track finding or track fitting
            l->SetX1(xtop);
            l->SetX2(xdown);
            l->Draw();
            l2->SetX1(xtopi);
            l2->SetX2(xdowni);
            l2->Draw();
        }

        // draw the z-x planes
        for (int lid = 1; lid<NZXP; lid++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
            ca_zx[lid]->cd();
            // draw Background graph for z-x planes
            gr_all[lid]->Draw("AP");
            text_cr_1[lid]->Draw();
            text_cr_2[lid]->Draw();
            for (int wid = 0; wid<NCEL; wid++){
                if (map_check[lid][wid]==1){
                    gr_wire[lid][wid]->Draw("PLSAME");
                    text_cr_1l[lid][wid]->Draw();
                }
                if (map_check[lid+1][wid]==1){
                    gr_wire[lid+1][wid]->Draw("PLSAME");
                    text_cr_2r[lid][wid]->Draw();
                }
#ifdef PRINT_CROSSPOINTS
                for (int wjd = 0; wjd < NCEL; wjd++){
                    if (!point_cross_wire[lid][wid][wjd]
                       ||!text_cross_wire[lid][wid][wjd]) continue;
                    point_cross_wire[lid][wid][wjd]->Draw();
                    text_cross_wire[lid][wid][wjd]->Draw();
                }
#endif
            }
            // draw the track and driftT lines on the z-x planes
            for (int ilr = 0; ilr<2; ilr++){
                for (int wid = 0; wid<NCEL; wid++){
                    if (check_zx[lid][wid])
                        l_zx[lid][wid][ilr]->Draw();
                    if (check_zx[lid+1][wid])
                        l_zx[lid+1][wid][ilr]->Draw();
                }
            }
            // position of the track point
            double z_track = 0;
            double x_track = 0;
            // is there a cross point?
            if (nHits_zx[lid]>0&&nHits_zx[lid+1]>0){
                for (int wid = 0; wid<NCEL; wid++){
                    if (!check_zx[lid][wid]) continue;
                    for (int wjd = 0; wjd<NCEL; wjd++){
                        if (!check_zx[lid+1][wjd]) continue;
                        // position of the cross point
                        for (int icombi = 0; icombi<4; icombi++){
                            double dd1 = dd_zx[lid][wid];
                            double dd2 = dd_zx[lid+1][wjd];
                            if (icombi<2) dd1 = -dd1; // reverse lid when icombi is 0 or 1
                            if (icombi%2==0) dd2 = -dd2; // reverse lid+1 when icombi is 0 or 2
                            double theta1 = theta_zx[lid][wid];
                            double theta2 = theta_zx[lid+1][wjd];
                            double sintheta12 = sin(theta1-theta2);
                            double zc_fix_slx = 0;
                            if (workMode>=1){
                                double deltaY = y_zx[lid+1][wjd]-y_zx[lid][wid];
                                zc_fix_slx = deltaY*slx/(tan(theta2)-tan(theta1));
                            }
                            double xc = mcp_xc[lid][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
                            double zc = mcp_zc[lid][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
                            point_cross_zx[lid][wid][wjd][icombi]->SetX(zc);
                            point_cross_zx[lid][wid][wjd][icombi]->SetY(xc);
                            point_cross_zx[lid][wid][wjd][icombi]->Draw();
                            if (workMode>=1){
                                if (icombi==3){
                                    double y_track = (y_zx[lid][wid]+y_zx[lid+1][wjd])/2.;
                                    z_track = zup+(y_track-yup)*slz;
                                    x_track = xup+(y_track-yup)*slx;
                                    double fd1 = fd_zx[lid][wid];
                                    double fd2 = fd_zx[lid+1][wjd];
                                    double xcf = mcp_xc[lid][wid][wjd]+fd1*sin(theta2)/(-sintheta12)+fd2*sin(theta1)/sintheta12;
                                    double zcf = mcp_zc[lid][wid][wjd]+fd1*cos(theta2)/(-sintheta12)+fd2*cos(theta1)/sintheta12+zc_fix_slx;
                                    gr_all[lid]->SetTitle(Form("DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) #Delta_{x}: %.0lf #mum #Delta_{z}: %.0lf #mum",lid,wid,dd1,(fd1-dd1)*1000,lid+1,wjd,dd2,(fd2-dd2)*1000,(xc-x_track)*1000,(zc-z_track)*1000));
                                }
                            }
                            else{
                                gr_all[lid]->SetTitle(Form("Layer %d and Layer %d",lid,lid+1));
                            }
                        }
                    }
                }
            }
            else{
                int wid = 0;
                for (; wid<NCEL; wid++){
                    if (check_zx[lid][wid]) break;
                }
                if (wid==NCEL) wid = NCEL/2;
                int wjd = 0;
                for (; wjd<NCEL; wjd++){
                    if (check_zx[lid+1][wjd]) break;
                }
                if (wjd==NCEL) wjd = NCEL/2;
                if (workMode>=1){
                    double y_track = (y_zx[lid][wid]+y_zx[lid+1][wjd])/2.; // take the y value from a previous event
                    z_track = zup+(y_track-yup)*slz;
                    x_track = xup+(y_track-yup)*slx;
                }
                gr_all[lid]->SetTitle(Form("Layer %d and Layer %d",lid,lid+1));
            }
            // draw the track point
            if (workMode>=1){
                point_track_zx[lid]->SetX(z_track);
                point_track_zx[lid]->SetY(x_track);
                point_track_zx[lid]->Draw();
            }
        }

        // save the canvases
		ca_xyADC->SaveAs(Form("xyADC.%d.pdf",triggerNumber));
		ca_xyADC->SaveAs(Form("xyADC.%d.png",triggerNumber));
        for (int lid = 1; lid<NZXP; lid++){
            ca_zx[lid]->SaveAs(Form("zx.%d.l%d.pdf",triggerNumber,lid));
            ca_zx[lid]->SaveAs(Form("zx.%d.l%d.png",triggerNumber,lid));
        }
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->SaveAs(Form("wf.%d.b%d.pdf",triggerNumber,bid));
            ca_WF[bid]->SaveAs(Form("wf.%d.b%d.png",triggerNumber,bid));
        }

		// FIXME: in interactive mode
//		ca_xyADC->WaitPrimitive();
	//	ca_xyADC->Update();
//		while(1){}
	}
	return 0;
}

double t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
	TF1* fl=0;
	TF1* fr=0;
	// FIXME
	//int index = (lid-1)*11+wid;
	int index = (6-1)*11;
	status = 0;
	if (time<=vtlel[index]&&time>vtler[index]){
		fl = f_left_end[index];
		status = -2;
	}
	else if (time<=vtlml[index]&&time>=vtlmr[index]){
		fl = f_left[index];
		status = -1;
	}
	if (time>=vtrel[index]&&time<vtrer[index]){
		fr = f_right_end[index];
		status = 2;
	}
	else if (time>=vtrml[index]&&time<=vtrmr[index]){
		fr = f_right[index];
		status = 1;
	}
	double dd=0;
	if (lr>=0){
		if (fr) dd = fr->Eval(time);
	}
	else{
		if (fl) dd = fl->Eval(time);
	}
	//std::cout<<"t2x("<<time<<","<<lid<<","<<wid<<","<<lr<<") = "<<dd<<std::endl;
	return dd;
}

