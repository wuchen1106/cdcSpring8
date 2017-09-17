#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraph.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>

//#define PRINT_CROSSPOINTS

#define NLAY 9
#define NZXP  8
#define NCEL 11

int main(int argc, char** argv){

	//===================Get Arguments============================
	if (argc<2) return -1;
	int runNo = (int)strtol(argv[1],NULL,10);
	int thelayer = 5;
	if (argc>=3){
	    thelayer = atoi(argv[2]);
    }
	int thewire = 5;
	if (argc>=4){
	    thewire = atoi(argv[3]);
    }
	TString suffix = "";
	if (argc>=5){
		suffix  = argv[4];
	}
    int iEntryStart = 0;
    int iEntryStop = 9;
	if (argc>=7){
	    iEntryStart = (int)strtol(argv[5],NULL,10);
	    iEntryStop = (int)strtol(argv[6],NULL,10);
    }
    else if (argc>=6){
	    iEntryStart=0;
	    iEntryStop=(int)strtol(argv[5],NULL,10)-1;
	}

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
    int     map_check[NLAY][NCEL];
    double  map_y_hit[NLAY];
	// mcp for cross points
    double mcp_xc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
    double mcp_zc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
	for(int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
		    map_x[ilayer][iwire][0] = 0;
		    map_y[ilayer][iwire][0] = 0;
		    map_z[ilayer][iwire][0] = 0;
		    map_x[ilayer][iwire][1] = 0;
		    map_y[ilayer][iwire][1] = 0;
		    map_z[ilayer][iwire][1] = 0;
			map_check[ilayer][iwire]=0;
			if (ilayer <NZXP){ // z-x planes corresponding to the layerID of the lower layer counting from 1 
                for (int jwire = 0; jwire<NCEL; jwire++){
                    mcp_xc[ilayer][iwire][jwire] = 999;
                    mcp_zc[ilayer][iwire][jwire] = 999;
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
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<NLAY){
			map_x[wp_lid][wp_wid][0] = wp_xhv;
			map_y[wp_lid][wp_wid][0] = wp_yhv;
			map_z[wp_lid][wp_wid][0] = -chamberHL;
			map_x[wp_lid][wp_wid][1] = wp_xro;
			map_y[wp_lid][wp_wid][1] = wp_yro;
			map_z[wp_lid][wp_wid][1] = chamberHL;
            map_y_hit[wp_lid] = (wp_yro+wp_yhv)/2.;
			map_ch[wp_lid][wp_wid] = wp_ch+wp_bid*48;
			map_check[wp_lid][wp_wid] = 1;
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

	//==================Get ADC==========================
	TChain * iChain_ADC = new TChain("tree","tree");
	iChain_ADC->Add(Form("../root/run_%0.6d_built.root",runNo));
	int adc[99][32];
	int tdc[99][32];
	int clockNumberDriftTime[99][32];
	int tdcNhit[99];
	int triggerNumber_ADC;
	iChain_ADC->SetBranchAddress("adc",adc);
	iChain_ADC->SetBranchAddress("driftTime",tdc);
	iChain_ADC->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	iChain_ADC->SetBranchAddress("tdcNhit",tdcNhit);
	iChain_ADC->SetBranchAddress("triggerNumber",&triggerNumber_ADC);

	//===================Get ROOT File============================
	TChain * iChain = new TChain("t","t");
	iChain->Add(Form("../root/t_%d.layer%d.",runNo,thelayer)+suffix+".root");
	int triggerNumber;
	double xup;
	double zup;
	double slx;
	double slz;
	double xupi;
	double zupi;
	double slix;
	double sliz;
	double chi2;
	std::vector<double> * i_driftD = 0;
	std::vector<double> * i_driftT = 0;
	std::vector<double> * i_fitD = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_layerID = 0;
	iChain->SetBranchAddress("wireID",&i_wireID);
	iChain->SetBranchAddress("layerID",&i_layerID);
	iChain->SetBranchAddress("driftD",&i_driftD);
	iChain->SetBranchAddress("driftT",&i_driftT);
	iChain->SetBranchAddress("fitD",&i_fitD);
	iChain->SetBranchAddress("inX",&xup);
	iChain->SetBranchAddress("inZ",&zup);
	iChain->SetBranchAddress("slX",&slx);
	iChain->SetBranchAddress("slZ",&slz);
	iChain->SetBranchAddress("iniX",&xupi);
	iChain->SetBranchAddress("iniZ",&zupi);
	iChain->SetBranchAddress("sliX",&slix);
	iChain->SetBranchAddress("sliZ",&sliz);
	iChain->SetBranchAddress("chi2",&chi2);
	iChain->SetBranchAddress("triggerNumber",&triggerNumber);

	//==================Prepare for drawing==========================
	//Prepare the Canvas for x-y plane and ADC
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
    TCanvas * ca_zx = new TCanvas("ca_zx","ca_zx",768,1024);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
	TPad * p_zx[NZXP]; // z-x planes corresponding to the layerID of the lower layer counting from 1
	for (int ilayer = 0; ilayer<NZXP; ilayer++){
	    int ipad = ilayer-1;
	    int ir = ipad/2;
	    int il = ipad%2;
	    double x1 = il/2.;
	    double x2 = (il+1)/2.;
	    double y1 = (4-ir)/4.;
	    double y2 = (3-ir)/4.;
	    p_zx[ilayer] = new TPad(Form("p_zx_%d",ilayer),"p_zx",x1,y1,x2,y2);
        p_zx[ilayer]->Draw();
        p_zx[ilayer]->SetGridx(1);
        p_zx[ilayer]->SetGridy(1);
	}
	//Prepare for ADC
	TGraph *gr_waveForm = 0;
	int vSample[32];
	for (int i=0; i<32; i++){
		vSample[i] = i;
	}
	TLatex *textTDC[32];
	TMarker *markerTDC[32];
	for (int j=0; j<32; j++) {
		textTDC[j] = new TLatex(0,235,"");
		textTDC[j]->SetTextSize(0.04);
		textTDC[j]->SetTextColor(kRed);
		markerTDC[j] = new TMarker(0,0,20);
		markerTDC[j]->SetMarkerSize(0.55);
	}
	//Prepare for hit circles on x-y plane
	TEllipse * ewiret[NLAY][NCEL];
	TEllipse * ewiret2[NLAY][NCEL];
	TText * text[NLAY][NCEL];
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			text[lid][wid] = 0;
			ewiret[lid][wid] = 0;
			ewiret2[lid][wid] = 0;
		}
	}
	//Prepare two tracks: l: track after fitting; l2: track before fitting, after track finding.
	TLine * l = new TLine();
	l->SetLineColor(kRed);
	l->SetY1(ytop);
	l->SetY2(ydown);
	TLine * l2 = new TLine();
	l2->SetLineColor(kBlue);
	l2->SetY1(ytop);
	l2->SetY2(ydown);
    // Background histogram for x-y plane
	TH2D * h0 = new TH2D("h0","h0",512,-130,130,512,500,650);
    // Prepare wires on z-x planes
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
    TGraph * gr_wire[NLAY][NCEL];
    for (int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
		    if (map_check[ilayer][iwire]==0) continue;
			gr_wire[ilayer][iwire] = new TGraph(2,map_z[ilayer][iwire],map_x[ilayer][iwire]);
			gr_wire[ilayer][iwire]->SetLineColor(color[iwire]);
			gr_wire[ilayer][iwire]->SetMarkerColor(color[iwire]);
        }
    }
    // Prepare driftT lines on z-x planes
    TLine * l_zx[NLAY][NCEL][4];
    for (int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
            for (int ilr = 0; ilr<2; ilr++){
                l_zx[ilayer][iwire][ilr] = 0;
            }
        }
    }
    // Prepare cross points of driftT lines on z-x planes
    TMarker * point_cross_zx[NZXP][4];
    for (int ilayer = 1; ilayer<NZXP; ilayer++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
        for (int icombi = 0; icombi<4; icombi++){
            if (icombi==3) // reverse lid when icombi is 0 or 1, reverse lid+1 when icombi is 0 or 2. Keep them unchanged only when icombi is 3
                point_cross_zx[ilayer][icombi] = new TMarker(0,0,20);
            else
                point_cross_zx[ilayer][icombi] = new TMarker(0,0,4);
            point_cross_zx[ilayer][icombi]->SetMarkerColor(kBlack);
            point_cross_zx[ilayer][icombi]->SetMarkerSize(0.3);
        }
    }
    // to record the position of hits in each layer
    // FIXME: should consider about how to support more or less than one hit in one layer.
    //        now we only record one hit per layer: the last one in that layer.
    int    wid_zx[NLAY];
    double fd_zx[NLAY];
    double dd_zx[NLAY];
    double theta_zx[NLAY];
    // Prepare track points on z-x planes
    TMarker * point_track_zx[NZXP];
    for (int ilayer = 1; ilayer<NZXP; ilayer++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
        point_track_zx[ilayer] = new TMarker(0,0,20);
        point_track_zx[ilayer]->SetMarkerColor(kRed);
        point_track_zx[ilayer]->SetMarkerSize(0.25);
    }
    // Draw wires on the z-x planes
    TGraph * gr_all[NZXP];
    for (int ilayer = 1; ilayer<NZXP; ilayer++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
        p_zx[ilayer]->cd();
        // Background graph for z-x planes
        gr_all[ilayer] = new TGraph(NLAY*NCEL*2,&(map_z[0][0][0]),&(map_x[0][0][0]));
        gr_all[ilayer]->SetTitle(Form("layer #%d and layer #%d",ilayer,ilayer+1));
        gr_all[ilayer]->SetMarkerColor(kWhite);
        gr_all[ilayer]->SetMarkerSize(0.1);
        gr_all[ilayer]->GetXaxis()->SetTitle("z [mm]");
        gr_all[ilayer]->GetYaxis()->SetTitle("x [mm]");
        gr_all[ilayer]->Draw("AP");
        for (int iwire = 0; iwire<NCEL; iwire++){
            TLatex * text1 = new TLatex(map_z[ilayer][0][0]-20,map_x[ilayer][0][0]+10,Form("# %d",ilayer));
            text1->SetTextSize(0.02);
            text1->Draw("SAME");
            TLatex * text2 = new TLatex(map_z[ilayer+1][0][1]+10,map_x[ilayer+1][0][1]+10,Form("# %d",ilayer+1));
            text2->SetTextSize(0.02);
            text2->Draw("SAME");
            if (map_check[ilayer][iwire]==1){
                gr_wire[ilayer][iwire]->Draw("PLSAME");
                TLatex * textl = new TLatex(map_z[ilayer][iwire][0]-10,map_x[ilayer][iwire][0],Form("%d",iwire));
                TLatex * textr = new TLatex(map_z[ilayer][iwire][1]-10,map_x[ilayer][iwire][1],Form("%d",iwire));
                textl->SetTextColor(color[iwire]);
                textr->SetTextColor(color[iwire]);
                textl->SetTextSize(0.02);
                textr->SetTextSize(0.02);
                textl->Draw("SAME");
                textr->Draw("SAME");
            }
            if (map_check[ilayer+1][iwire]==1){
                gr_wire[ilayer+1][iwire]->Draw("PLSAME");
                TLatex * textl = new TLatex(map_z[ilayer+1][iwire][0]+10,map_x[ilayer+1][iwire][0],Form("%d",iwire));
                TLatex * textr = new TLatex(map_z[ilayer+1][iwire][1]+10,map_x[ilayer+1][iwire][1],Form("%d",iwire));
                textl->SetTextColor(color[iwire]);
                textr->SetTextColor(color[iwire]);
                textl->SetTextSize(0.02);
                textr->SetTextSize(0.02);
                textl->Draw("SAME");
                textr->Draw("SAME");
            }
#ifdef PRINT_CROSSPOINTS
            for (int jwire = 0; jwire < NCEL; jwire++){
                if (fabs(mcp_zc[ilayer][iwire][jwire])>300) continue;
                TMarker * p = new TMarker(mcp_zc[ilayer][iwire][jwire],mcp_xc[ilayer][iwire][jwire],20);
                p->SetMarkerColor(color[iwire]);
                p->SetMarkerSize(0.4);
                p->Draw("SAME");
                TLatex * text = new TLatex(mcp_zc[ilayer][iwire][jwire],mcp_xc[ilayer][iwire][jwire]+5,Form("%d,%d",iwire,jwire));
                text->SetTextColor(color[iwire]);
                text->SetTextSize(0.02);
                text->Draw("SAME");
            }
#endif
        }
    }

	//===================Loop in Events============================
	for ( int iEntry = iEntryStart; iEntry<iEntryStop; iEntry++){
		if (iEntry%1000==0) printf("%lf\n",iEntry);
		iChain->GetEntry(iEntry);
		int the_ihit = 0;
		int the_ch = -1;
		for (; the_ihit<i_driftD->size(); the_ihit++){
			int lid = (*i_layerID)[the_ihit];
			int wid = (*i_wireID)[the_ihit];
			if (lid==thelayer&&wid==thewire){
				the_ch = map_ch[lid][wid];
				break;
			}
		}
		if (the_ch==-1) continue; // only show the events with hit in the target cell

//		if ((*i_fitD)[the_ihit]>-2.3||(*i_fitD)[the_ihit]<-3.2||fabs((*i_driftD)[the_ihit]-(*i_fitD)[the_ihit])<0.5||chi2>1||fabs(slz)>0.15) continue;
//		if ((*i_fitD)[the_ihit]>0.32||(*i_fitD)[the_ihit]<0.28||chi2>3) continue;
//		if ((*i_driftT)[the_ihit]>80||chi2>3) continue;

        // clear objects from the previous event
		for (int lid = 0; lid<NLAY; lid++){
            wid_zx[lid] = -1;
            dd_zx[lid] = 0;
            fd_zx[lid] = 0;
			for (int wid = 0; wid<NCEL; wid++){
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
				for (int ilr = 0; ilr<2; ilr++){
                    if (l_zx[lid][wid][ilr]){
                        delete l_zx[lid][wid][ilr];
                        l_zx[lid][wid][ilr] = 0;
                    }
                }
			}
		}
        if (gr_waveForm) delete gr_waveForm; gr_waveForm = 0;

        // Draw the background histogram for x-y plane
		p_xyADC[0]->cd();
		h0->SetTitle(Form("Entry#%d, TriggerNumber#%d, chi2 = %.2f",iEntry,triggerNumber,chi2));
		h0->Draw();

        // get the tracks and draw them on the x-y plane
        double zdown = zup + (ydown-yup)*slz;
        double xdown = xup + (ydown-yup)*slx;
        double xtop = xup + (ytop-yup)*slx;
        double zdowni = zupi + (ydown-yup)*sliz;
        double xdowni = xupi + (ydown-yup)*slix;
        double xtopi = xupi + (ytop-yup)*slix;
		l->SetX1(xtop);
		l->SetX2(xdown);
		l->Draw("SAME");
		l2->SetX1(xtopi);
		l2->SetX2(xdowni);
		l2->Draw("SAME");

        // get the hits and draw them on the x-y plane
		for (int ihit = 0; ihit<i_driftD->size(); ihit++){
			int lid = (*i_layerID)[ihit];
			if (lid<=0) continue;
			int wid = (*i_wireID)[ihit];
			double fitd = (*i_fitD)[ihit]; // !!! make sure it's mm
			double dd = (*i_driftD)[ihit]; // !!! make sure it's mm
			double dt = (*i_driftT)[ihit];
			double wxro = map_x[lid][wid][1];
			double wyro = map_y[lid][wid][1];
            double wzro = chamberHL;
			double wxhv = map_x[lid][wid][0];
			double wyhv = map_y[lid][wid][0];
            double wzhv = -chamberHL;
			double wy = (wyro+wyhv)/2.;
			double wz = ((yup-wy)*zdown+(wy-ydown)*zup)/(yup-ydown);
			double wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
			wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
			map_y_hit[lid] = wy;
			double delta = dd*chamberHL*2/sqrt(chamberHL*chamberHL*4+(map_x[lid][wid][0]-map_x[lid][wid][1])*(map_x[lid][wid][0]-map_x[lid][wid][1]));
			for (int ilr = 0; ilr<2; ilr++){
                l_zx[lid][wid][ilr] = new TLine();
                l_zx[lid][wid][ilr]->SetLineColor(color[wid]);
                if (ilr)
                    l_zx[lid][wid][ilr]->SetLineStyle(2);
                else
                    l_zx[lid][wid][ilr]->SetLineStyle(3);
                l_zx[lid][wid][ilr]->SetX1(-chamberHL);
                l_zx[lid][wid][ilr]->SetX2(chamberHL);
                l_zx[lid][wid][ilr]->SetY1(wxhv+(ilr?delta:-delta));
                l_zx[lid][wid][ilr]->SetY2(wxro+(ilr?delta:-delta));
            }
            wid_zx[lid] = wid;
            dd_zx[lid] = dd;
            fd_zx[lid] = fitd;
            theta_zx[lid] = atan(-(map_x[lid][wid][0]-map_x[lid][wid][1])/chamberHL/2); // rotation angle w.r.t the dart plane: read out  plane, i.e. z>0, i.e. iplane = 1; positive rotation angle point to -x direction

			ewiret[lid][wid] = new TEllipse(wx,wy,dd,dd);
			ewiret[lid][wid]->SetFillStyle(0);
			ewiret[lid][wid]->SetLineColor(kRed);
			ewiret[lid][wid]->Draw("SAME"); // Draw hits.

			if (lid==thelayer&&wid==thewire) text[lid][wid]= new TText(wx,wy,Form("%d,%d,%.1lf,%.1lf",lid,wid,fitd,dt));
			else text[lid][wid]= new TText(wx,wy,Form("%d,%d",lid,wid));
			text[lid][wid]->SetTextSize(0.02);
			text[lid][wid]->Draw("SAME");

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

        // get ADC information and draw the waveform
		p_xyADC[1]->cd();
		for (int iEntry2 = triggerNumber; iEntry2>=0 ;iEntry2--){
			iChain_ADC->GetEntry(iEntry2);
			if (triggerNumber==triggerNumber_ADC) break;
		}
		gr_waveForm = new TGraph(32,vSample,adc[the_ch]);
		gr_waveForm->Draw("APL");
		for (int ihit=0; ihit<tdcNhit[the_ch]; ihit++) {
		    int clk = clockNumberDriftTime[the_ch][ihit];
			textTDC[ihit]->SetText(clk,adc[the_ch][clk],Form("%d",(int)(tdc[the_ch][ihit])));
			textTDC[ihit]->Draw();
			markerTDC[ihit]->SetX(clk);
			markerTDC[ihit]->SetY(adc[the_ch][clk]);
			markerTDC[ihit]->SetMarkerColor(kRed);
			//if (adc[the_ch][clk]>Hmin[the_ch]) markerTDC[ihit]->SetMarkerColor(kRed);
			//else markerTDC[ihit]->SetMarkerColor(kBlue);
			markerTDC[ihit]->Draw();
		}

        // draw the track and driftT lines on the z-x planes
        for (int lid = 1; lid<NZXP; lid++){
            p_zx[lid]->cd();
            for (int ilr = 0; ilr<2; ilr++){
                for (int wid = 0; wid<NCEL; wid++){
                    if (l_zx[lid][wid][ilr])
                        l_zx[lid][wid][ilr]->Draw("SAME");
                    if (l_zx[lid+1][wid][ilr])
                        l_zx[lid+1][wid][ilr]->Draw("SAME");
                }
            }
            // position of the track point
            double y = (map_y_hit[lid]+map_y_hit[lid+1])/2.;
            double z = zup+(y-yup)*slz;
            double x = xup+(y-yup)*slx;
            // is there a cross point?
            if (wid_zx[lid]>=0&&wid_zx[lid+1]>=0){
                // position of the cross point
                for (int icombi = 0; icombi<4; icombi++){
                    double fd1 = fd_zx[lid];
                    double fd2 = fd_zx[lid+1];
                    double dd1 = dd_zx[lid];
                    double dd2 = dd_zx[lid+1];
                    if (icombi<2) dd1 = -dd1; // reverse lid when icombi is 0 or 1
                    if (icombi%2==0) dd2 = -dd2; // reverse lid+1 when icombi is 0 or 2
                    double theta1 = theta_zx[lid];
                    double theta2 = theta_zx[lid+1];
                    double sintheta12 = sin(theta1-theta2);
                    double deltaY = map_y_hit[lid+1]-map_y_hit[lid];
                    double zc_fix_slx = deltaY*slx/(tan(theta2)-tan(theta1));
                    double xc = mcp_xc[lid][wid_zx[lid]][wid_zx[lid+1]]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
                    double zc = mcp_zc[lid][wid_zx[lid]][wid_zx[lid+1]]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
//                  point_cross_zx[lid][icombi]->SetMarkerColor(color[wid_zx[lid]]);
                    point_cross_zx[lid][icombi]->SetX(zc);
                    point_cross_zx[lid][icombi]->SetY(xc);
                    point_cross_zx[lid][icombi]->Draw("SAME");
                    if (icombi==3){
                        double xcf = mcp_xc[lid][wid_zx[lid]][wid_zx[lid+1]]+fd1*sin(theta2)/(-sintheta12)+fd2*sin(theta1)/sintheta12;
                        double zcf = mcp_zc[lid][wid_zx[lid]][wid_zx[lid+1]]+fd1*cos(theta2)/(-sintheta12)+fd2*cos(theta1)/sintheta12+zc_fix_slx;
                        gr_all[lid]->SetTitle(Form("DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) #Delta_{x}: %.0lf #mum #Delta_{z}: %.0lf #mum",lid,wid_zx[lid],dd1,(fd1-dd1)*1000,lid+1,wid_zx[lid+1],dd2,(fd2-dd2)*1000,(xc-x)*1000,(zc-z)*1000));
//                        gr_all[lid]->SetTitle(Form("#Wire: %d,%d D_{meas}: %.2lf,%.2lf mm #Delta_{D}: %.0lf,%.0lf #mum #Delta_{x,xf}: %.0lf,%.0lf #mum #Delta_{z,zf}: %.0lf,%.0lf #mum",wid_zx[lid],wid_zx[lid+1],dd1,dd2,(fd1-dd1)*1000,(fd2-dd2)*1000,(xc-x)*1000,(xcf-x)*1000,(zc-z)*1000,(zcf-z)*1000));
//                        printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",iEntry,lid,slx,slz,fd1-dd1,fd2-dd2,x,xc,xcf,z,zc,zcf);
                    }
                }
            }
            else{
                for (int icombi = 0; icombi<4; icombi++){
                    point_cross_zx[lid][icombi]->SetX(chamberHL*2); // move the point out of the canvas
                    point_cross_zx[lid][icombi]->Draw("SAME");
                }
                gr_all[lid]->SetTitle("No crosspoint");
            }
            // draw the track point
            point_track_zx[lid]->SetX(z);
            point_track_zx[lid]->SetY(x);
            point_track_zx[lid]->Draw("SAME");
        }

        // save the canvases
		ca_xyADC->SaveAs(Form("%d.xyADC.pdf",iEntry));
		ca_xyADC->SaveAs(Form("%d.xyADC.png",iEntry));
        ca_zx->SaveAs(Form("%d.zx.pdf",iEntry));
//        ca_zx->SaveAs(Form("%d.zx.png",iEntry)); // too vague to draw...

		// FIXME: in interactive mode
//		ca_xyADC->WaitPrimitive();
	//	ca_xyADC->Update();
//		while(1){}
	}
	return 0;
}
