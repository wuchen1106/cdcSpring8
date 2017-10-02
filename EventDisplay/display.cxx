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
#include <stdlib.h>

//#define PRINT_CROSSPOINTS

#define XMAX 130
#define ZMAX 350

#define NLAY 9
#define NZXP 8
#define NCEL 11
#define NCELA 88

#define NBRD 2
#define NCHS 48
#define NCHT 96
#define NSMP 32
#define ADCMIN 125
#define ADCMAX 700

#define NCUT 50

//===================About xt============================
double tres = 0;
TF1 * f_left_end[NCELA];
TF1 * f_right_end[NCELA];
TF1 * f_left[NCELA];
TF1 * f_right[NCELA];
int vtrel[NCELA];
int vtrer[NCELA];
int vtlel[NCELA];
int vtler[NCELA];
int vtrml[NCELA];
int vtrmr[NCELA];
int vtlml[NCELA];
int vtlmr[NCELA];

double tdc2t(int tdc);

int t0 = -849;
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
	int pk_sumcut = 30;
	int pk_aacut = 30;
	printf("runNo:              %d\n",runNo);
	printf("workMode:           %d\n",workMode);
	printf("test wire:          [%d,%d]\n",thelayer,thewire);
	printf("suffix:             %s\n",suffix.Data());
	printf("Entries:            %d~%d\n",iEntryStart,iEntryStop);
	printf("min ADC sum:        %d\n",pk_sumcut);
	printf("min ADC area:       %d\n",pk_aacut);
	printf("t0:                 %d\n",t0);

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
            map_bid[wp_lid][wp_wid] = wp_bid;
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
	for (int i = 0; i<NCELA; i++){
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
		vtlel[(lid)*11+wid] = tlel;
		vtler[(lid)*11+wid] = tler;
		vtlml[(lid)*11+wid] = tlml;
		vtlmr[(lid)*11+wid] = tlmr;
		vtrel[(lid)*11+wid] = trel;
		vtrer[(lid)*11+wid] = trer;
		vtrml[(lid)*11+wid] = trml;
		vtrmr[(lid)*11+wid] = trmr;
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

	//==================Get Peaks==========================
	TChain * iChain_p = new TChain("t","t");
	iChain_p->Add(Form("../root/p_%d.root",runNo));
	double pk_aa[NCHT];
	int pk_np[NCHT][NCUT];
	std::vector<std::vector<int> > * pk_peak = 0;
	std::vector<std::vector<int> > * pk_width = 0;
	std::vector<std::vector<int> > * pk_clk = 0;
	std::vector<std::vector<int> > * pk_tdc = 0;
	std::vector<std::vector<double> > * pk_sum = 0;
	iChain_p->SetBranchAddress("aa",pk_aa);
	iChain_p->SetBranchAddress("np",pk_np);
	iChain_p->SetBranchAddress("peak",&pk_peak);
	iChain_p->SetBranchAddress("width",&pk_width);
	iChain_p->SetBranchAddress("clk",&pk_clk);
	iChain_p->SetBranchAddress("tdc",&pk_tdc);
	iChain_p->SetBranchAddress("sum",&pk_sum);
	int pk_ip[NCHT]; // pick peak
	int pk_jp[NCHT]; // post peak

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
	TPad * pad_WF[NBRD][NCHS];
	for (int bid = 0; bid<NBRD; bid++){
        ca_WF[bid] = new TCanvas(Form("ca_WF_%d",bid),"ca_WF",896,896);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        for (int i = 0; i<8; i++){
            for (int j = 0; j<6; j++){
                int index = j*8+i;
                pad_WF[bid][index] = new TPad(Form("pad_%d_%d_%d",bid,i,j),Form("pad_%d_%d_%d",bid,i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
                pad_WF[bid][index]->Draw();
                pad_WF[bid][index]->SetGridx(1);
                pad_WF[bid][index]->SetGridy(1);
            }
        }
    }
	//Prepare the Canvas for x-y plane and target chanel ADC
	TCanvas * ca_xyADC = new TCanvas("ca_xyADC","ca_xyADC",896,1024);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
	TPad * pad_xyADC[2];
	for (int ipad = 0; ipad<2; ipad++){
        pad_xyADC[ipad] = new TPad(Form("pad_xyADC_%i",ipad),"pad_xyADC",0,ipad?0:0.2,1,ipad?0.2:1);
        pad_xyADC[ipad]->Draw();
        pad_xyADC[ipad]->SetGridx(1);
        pad_xyADC[ipad]->SetGridy(1);
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
	TLatex *textWF[NBRD][NCHS];
	TLatex *textTDC[NBRD][NCHS][NSMP];
	TMarker *markerTDC[NBRD][NCHS][NSMP];
    for (int bid = 0; bid<NBRD; bid++){
        for (int ch = 0; ch<NCHS; ch++){
            textWF[bid][ch] = new TLatex(0,0,"");
            textWF[bid][ch]->SetTextSize(0.04);
            for (int sid=0; sid<NSMP; sid++) {
                textTDC[bid][ch][sid] = new TLatex(0,0,"");
                textTDC[bid][ch][sid]->SetTextSize(0.04);
                markerTDC[bid][ch][sid] = new TMarker(0,0,20);
                markerTDC[bid][ch][sid]->SetMarkerSize(0.3);
            }
        } 
	}
	//Prepare for hit circles on x-y plane
	TText * text[NLAY][NCEL];
	TEllipse * ewiret_pick[NLAY][NCEL];  // value from track finding/fitting;
	TEllipse * ewiret_pick_ini[NLAY][NCEL]; // initial value for track fitting;
	TEllipse * ewiret_pre[NLAY][NCEL]; // peak before the picked one
	TEllipse * ewiret_post[NLAY][NCEL]; // peak after the picked one
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			ewiret_pick[lid][wid] = new TEllipse(0,0,1,1);
			ewiret_pick_ini[lid][wid] = new TEllipse(0,0,1,1);
			ewiret_pre[lid][wid] = new TEllipse(0,0,1,1);
			ewiret_post[lid][wid] = new TEllipse(0,0,1,1);
			ewiret_pick[lid][wid]->SetFillStyle(0);
			ewiret_pick[lid][wid]->SetLineColor(kRed);
            ewiret_pick_ini[lid][wid]->SetFillStyle(0);
            ewiret_pick_ini[lid][wid]->SetLineColor(kBlue);
            ewiret_pre[lid][wid]->SetFillStyle(0);
            ewiret_pre[lid][wid]->SetLineColor(kMagenta);
            ewiret_post[lid][wid]->SetFillStyle(0);
            ewiret_post[lid][wid]->SetLineColor(kOrange);
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
	TString prefix = "";
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

        // Get ADC and peak information
		for (int iEntry2 = triggerNumber; iEntry2>=0 ;iEntry2--){ // in case that the raw ROOT file missed some triggerNumbers
			iChain_ADC->GetEntry(iEntry2);
			if (triggerNumber==triggerNumber_ADC){
			    iChain_p->GetEntry(iEntry2);
			    break;
            }
		}
        int pk_nh = 0;
		for (int chg = 0; chg<NCHT; chg++){
            pk_ip[chg] = -1;
            int thepeak = 0;
			int tdcNhitwire = pk_np[chg][0];
			for(;thepeak<tdcNhitwire; thepeak++){ // find the first peak over sumcut
				if ((*pk_sum)[chg][thepeak]>=pk_sumcut) break;
			}
			if (thepeak==tdcNhitwire){// didn't find the peak, choose the first one;
				thepeak=0;
			}
			if (pk_aa[chg]>=pk_aacut){// didn't pass aa cut, consider no pick
			    pk_nh++;
            }
            else{
			    thepeak = -1;
			}
			pk_ip[chg] = thepeak;
            int postpeak = thepeak+1;
            for (; postpeak<tdcNhitwire; postpeak++){ // post peak;
                if ((*pk_sum)[chg][postpeak]>pk_sumcut) break;
            }
            if (postpeak==tdcNhitwire){
                // FIXME: should we consider small post peak?
//                postpeak = thepeak==tdcNhitwire-1?-1:thepeak+1;
                postpeak = -1;
            }
            pk_jp[chg] = postpeak;
            if (thelayer<0||thewire<0){ // in case we didn't specify a target cell, choose an interesting one
                if (postpeak>=0){ // post peak
                    the_bid = chg/NCHT;
                    the_ch = chg%NCHT;
                }
                else if (thepeak>0){ // pre peak
                    the_bid = chg/NCHT;
                    the_ch = chg%NCHT;
                }
            }
		}
		// ************Cutting on hits****************
		if (pk_nh<=7) prefix = "incom.";
		else if (pk_nh==8) prefix = "single.";
		else if (pk_nh==9) prefix = "single.";
		else if (pk_nh==10) prefix = "n10.";
		else if (pk_nh==11) prefix = "n11.";
		else if (pk_nh==12) prefix = "n12.";
		else if (pk_nh>=13) prefix = "multi.";
		// ************Cutting on peaks in ch15****************
//		int thepeak = 0;
//		for (; thepeak<(*pk_peak)[15].size(); thepeak++){
//			if ((*pk_sum)[15][thepeak]>=smin&&(*pk_tdc)[15][thepeak]<=tmax+tres&&(*pk_tdc)[15][thepeak]>=t0-tres) break;
//		}
//		if (thepeak==(*pk_peak)[15].size()){
//			if (pk_aa[15]>aamin) prefix = "missed.";
//		}
//		else{
//			if (thepeak<(*pk_peak)[15].size()-1&&(*pk_sum)[15][thepeak]<(*pk_sum)[15][thepeak+1])
//				prefix = "shadowed.";
//			if (thepeak>0&&(*pk_tdc)[15][thepeak-1]<=tmax+tres&&(*pk_tdc)[15][thepeak-1]>=t0-tres)
//				prefix = "prepeak.";
//		}


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
                pad_WF[bid][ch]->cd();
                gr_waveForm[bid][ch] = new TGraph(NSMP,vSample,adc[chg]);
                gr_waveForm[bid][ch]->SetTitle(Form("Channel %d Layer %d Wire %d",ch,map_lid[bid][ch],map_wid[bid][ch]));
                gr_waveForm[bid][ch]->GetXaxis()->SetRangeUser(0,NSMP-1);
                gr_waveForm[bid][ch]->GetYaxis()->SetRangeUser(ADCMIN,ADCMAX);
                gr_waveForm[bid][ch]->GetXaxis()->SetTitle("Sample Index");
                gr_waveForm[bid][ch]->GetYaxis()->SetTitle("ADC");
                gr_waveForm[bid][ch]->SetMarkerStyle(20);
                gr_waveForm[bid][ch]->SetMarkerSize(0.3);
                gr_waveForm[bid][ch]->Draw("APL");
                int ipk = pk_ip[chg];
                int jpk = pk_jp[chg];
                if (ipk>=0){
                    textWF[bid][ch]->SetText(1,ADCMAX-50,Form("%d peaks area %.0lf pick %d w %d s %.0lf p %d",pk_np[chg][0],pk_aa[chg],ipk,(*pk_width)[chg][ipk],(*pk_sum)[chg][ipk],(*pk_peak)[chg][ipk]));
                    if (jpk>=0) // post peak
                        textWF[bid][ch]->SetTextColor(kOrange);
                    else if (ipk>0) // pre peak
                        textWF[bid][ch]->SetTextColor(kMagenta);
                    else
                        textWF[bid][ch]->SetTextColor(kRed);
                }
                else{
                    textWF[bid][ch]->SetText(1,ADCMAX-50,Form("%d peaks area %.0lf no pick",pk_np[chg][0],pk_aa[chg]));
                    if (pk_np[chg][0]>0)
                        textWF[bid][ch]->SetTextColor(kGreen);
                    else
                        textWF[bid][ch]->SetTextColor(kBlue);
                }
                textWF[bid][ch]->Draw();
                for (int ihit=0; ihit<tdcNhit[chg]; ihit++) {
                    int clk = clockNumberDriftTime[chg][ihit];
                    markerTDC[bid][ch][ihit]->SetX(clk);
                    markerTDC[bid][ch][ihit]->SetY(adc[chg][clk]);
                    if (ihit==ipk){
                        markerTDC[bid][ch][ihit]->SetMarkerColor(kRed);
                        textTDC[bid][ch][ihit]->SetTextColor(kRed);
                    }
                    else if (ipk>0&&ihit==0){ // pre peak
                        markerTDC[bid][ch][ihit]->SetMarkerColor(kMagenta);
                        textTDC[bid][ch][ihit]->SetTextColor(kMagenta);
                    }
                    else if (ihit==jpk){ // post peak
                        markerTDC[bid][ch][ihit]->SetMarkerColor(kOrange);
                        textTDC[bid][ch][ihit]->SetTextColor(kOrange);
                    }
                    else{
                        markerTDC[bid][ch][ihit]->SetMarkerColor(kBlue);
                        textTDC[bid][ch][ihit]->SetTextColor(kBlue);
                    }
                    markerTDC[bid][ch][ihit]->Draw();
                    textTDC[bid][ch][ihit]->SetText(clk,adc[chg][clk]+(0.5-ihit%2)*50,Form("%d",(int)(tdc[chg][ihit])));
                    textTDC[bid][ch][ihit]->Draw();
                }
            }
        }

        // draw the target channel ADC
		pad_xyADC[1]->cd();
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
			int dd_status = 0;
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
                dd = t2x(dt,lid,wid,0,dd_status);
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
            ewiret_pick[lid][wid]->SetX1(wx);
            ewiret_pick[lid][wid]->SetY1(wy);
            ewiret_pick[lid][wid]->SetR1(dd);
            ewiret_pick[lid][wid]->SetR2(dd);
            if (abs(dd_status)>2) ewiret_pick[lid][wid]->SetLineStyle(2);
            else ewiret_pick[lid][wid]->SetLineStyle(1);
            ewiret_pre[lid][wid]->SetX1(wx);
            ewiret_pre[lid][wid]->SetY1(wy);
            ewiret_post[lid][wid]->SetX1(wx);
            ewiret_post[lid][wid]->SetY1(wy);
            // set the text of each hit
			if (lid==thelayer&&wid==thewire) text[lid][wid]= new TText(wx,wy,Form("%d,%d,%.1lf,%.1lf",lid,wid,fitd,dt));
			else text[lid][wid]= new TText(wx,wy,Form("%d,%d,%.2f",lid,wid,dd));
			text[lid][wid]->SetTextSize(0.02);

            if (workMode>=2){
                // do the correction according to the input track parameters from track finding
                wz = ((yup-wy)*zdowni+(wy-ydown)*zupi)/(yup-ydown);
                wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
                ewiret_pick_ini[lid][wid]->SetX1(wx);
                ewiret_pick_ini[lid][wid]->SetY1(wy);
                ewiret_pick_ini[lid][wid]->SetR1(dd);
                ewiret_pick_ini[lid][wid]->SetR2(dd);
                if (abs(dd_status)>2) ewiret_pick_ini[lid][wid]->SetLineStyle(2);
                else ewiret_pick_ini[lid][wid]->SetLineStyle(1);
            }
		}

        // Draw the background histogram for x-y plane
		pad_xyADC[0]->cd();
		gr_wireCenter->SetTitle(Form("Entry %d, Trigger Number %d, nHits = %d, chi2 = %.2f",iEntry,triggerNumber,nHitsAll,chi2));
		gr_wireCenter->Draw("AP");
        // Draw the hit circles
		for (int ihit = 0; ihit<i_driftT->size(); ihit++){
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			if (ewiret_pick[lid][wid]&&text[lid][wid]){
                ewiret_pick[lid][wid]->Draw(); // from track fitting/finding
                text[lid][wid]->Draw();
                int chg = map_bid[lid][wid]*NCHS+map_ch[lid][wid];
                int ipeak = pk_ip[chg];
                if (ipeak>0){
                    int status;
                    int ddp = t2x(tdc2t((*pk_tdc)[chg][0]),lid,wid,0,status);
                    ewiret_pre[lid][wid]->SetR1(ddp);
                    ewiret_pre[lid][wid]->SetR2(ddp);
                    ewiret_pre[lid][wid]->Draw();
                    if (abs(status)>2) ewiret_pre[lid][wid]->SetLineStyle(2);
                    else ewiret_pre[lid][wid]->SetLineStyle(1);
                }
                int jpk = pk_jp[chg];
                if (jpk<pk_np[chg][0]&&jpk>=0){
                    int status;
                    int ddp = t2x(tdc2t((*pk_tdc)[chg][jpk]),lid,wid,0,status);
                    ewiret_post[lid][wid]->SetR1(ddp);
                    ewiret_post[lid][wid]->SetR2(ddp);
                    ewiret_post[lid][wid]->Draw();
                    if (abs(status)>2) ewiret_post[lid][wid]->SetLineStyle(2);
                    else ewiret_post[lid][wid]->SetLineStyle(1);
                }
            }
            if (workMode>=2&&ewiret_pick_ini[lid][wid]){
                ewiret_pick_ini[lid][wid]->Draw(); // input values for track fitting.
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
                            if (zc>-chamberHL&&zc<chamberHL){
                                point_cross_zx[lid][wid][wjd][icombi]->SetX(zc);
                                point_cross_zx[lid][wid][wjd][icombi]->SetY(xc);
                                point_cross_zx[lid][wid][wjd][icombi]->Draw();
                            }
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
		ca_xyADC->SaveAs(prefix+Form("xyADC.%d.pdf",triggerNumber));
		ca_xyADC->SaveAs(prefix+Form("xyADC.%d.png",triggerNumber));
        for (int lid = 1; lid<NZXP; lid++){
            ca_zx[lid]->SaveAs(prefix+Form("zx.%d.l%d.pdf",triggerNumber,lid));
            ca_zx[lid]->SaveAs(prefix+Form("zx.%d.l%d.png",triggerNumber,lid));
        }
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->SaveAs(prefix+Form("wf.%d.b%d.pdf",triggerNumber,bid));
            ca_WF[bid]->SaveAs(prefix+Form("wf.%d.b%d.png",triggerNumber,bid));
        }

		// FIXME: in interactive mode
//		ca_xyADC->WaitPrimitive();
	//	ca_xyADC->Update();
//		while(1){}
	}
	return 0;
}

double t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
	TF1* f=0;
	// FIXME
	//int index = (lid)*NCEL+wid;
	int index = 4*NCEL;
	status = 0;
	if (lr>=0){
        if (time<=vtlel[index]&&time>vtler[index]){
            f = f_left_end[index];
            status = -2;
        }
        else if (time<=vtlml[index]&&time>=vtlmr[index]){
            f = f_left[index];
            status = -1;
        }
        else if (time>vtlel[index]){
            f = f_left_end[index];
            status = -3;
        }
    }
    else{
        if (time>=vtrel[index]&&time<vtrer[index]){
            f = f_right_end[index];
            status = 2;
        }
        else if (time>=vtrml[index]&&time<=vtrmr[index]){
            f = f_right[index];
            status = 1;
        }
        else if (time>vtrer[index]){
            f = f_right_end[index];
            status = 3;
        }
    }
	double dd=0;
    if (f) dd = f->Eval(time);
    else{
        if (!status){ // dt too small
            dd = 0;
        }
        else{ // dt too large
            dd = f->Eval(vtrer[index]);
        }
    }
	return dd;
}

double tdc2t(int tdc){
    return (tdc-t0)/0.96;
}
