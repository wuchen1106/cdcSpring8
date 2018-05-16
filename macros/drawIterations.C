#include "header.h"
int drawIterations(){
//	int runNo = 1012;
//	TString runname="1012.0122";
//	TString original="1012.original";
//	bool isMC = false;
//	int Niters = 9;
//	int runNo = 117;
//	TString runname="117.0124";
//	TString original="117.original";
//	bool isMC = false;
//	int Niters = 17;
	int runNo = 117;
	TString runname="117.0126";
	TString original="117.original";
	bool isMC = false;
	int Niters = 50;
//	int runNo = 100007;
//	TString runname="100007.0125.old";
//	TString original="100007.Garfield";
//	bool isMC = true;
//	int Niters = 40;
//	int runNo = 100007;
//	TString runname="100007.0125";
//	TString original="100007.Garfield";
//	bool isMC = true;
//	int Niters = 70;
//	int runNo = 100008;
//	TString runname="100008.0125";
//	TString original="100008.Garfield";
//	bool isMC = true;
//	int Niters = 74;
//	int runNo = 100009;
//	TString runname="100009.0125";
//	TString original="100009.Garfield";
//	bool isMC = true;
//	int Niters = 70;
//	int runNo = 100010;
//	TString runname="100010.0125";
//	TString original="100010.Garfield";
//	bool isMC = true;
//	int Niters = 25;
//	int runNo = 100011;
//	TString runname="100011.0125";
//	TString original="100011.Garfield";
//	bool isMC = true;
//	int Niters = 25;
//	int runNo = 100012;
//	TString runname="100012.0125";
//	TString original="100012.Garfield";
//	bool isMC = true;
//	int Niters = 25;

	bool drawXTs = false;
	int theLayer = 4;

	int color[NITERSMAX];
	for (int icolor = 0; icolor<NITERSMAX; icolor++){
		color[icolor] = kGray;
	}
	color[0] = kBlack;
	color[1] = kBlue-2;
	color[2] = kBlue+2;
	color[3] = kGreen-2;
	color[4] = kGreen+2;
	color[Niters-2] = kMagenta-2;
	color[Niters-1] = kMagenta+2;
	color[Niters] = kRed;

	int colorLayers[NLAY];
	icolor = 1;
	if ( icolor<NLAY) colorLayers[icolor] = kBlack; icolor++;
	if ( icolor<NLAY) colorLayers[icolor] = kGreen; icolor++;
	if ( icolor<NLAY) colorLayers[icolor] = kBlue; icolor++;
	if ( icolor<NLAY) colorLayers[icolor] = kRed; icolor++;
	if ( icolor<NLAY) colorLayers[icolor] = kMagenta; icolor++;
	if ( icolor<NLAY) colorLayers[icolor] = kMagenta+3; icolor++;
	if ( icolor<NLAY) colorLayers[icolor] = kCyan; icolor++;
	for (;icolor<NLAY; icolor++){
		colorLayers[icolor] = icolor;
	}
	int style[NCEL];
	int istyle = 0;
	if (istyle<NCEL) style[istyle] = 2; istyle++;
	if (istyle<NCEL) style[istyle] = 3; istyle++;
	if (istyle<NCEL) style[istyle] = 5; istyle++;
	if (istyle<NCEL) style[istyle] = 21; istyle++;
	if (istyle<NCEL) style[istyle] = 20; istyle++;
	if (istyle<NCEL) style[istyle] = 22; istyle++;
	if (istyle<NCEL) style[istyle] = 23; istyle++;
	if (istyle<NCEL) style[istyle] = 24; istyle++;
	if (istyle<NCEL) style[istyle] = 25; istyle++;
	if (istyle<NCEL) style[istyle] = 26; istyle++;
	if (istyle<NCEL) style[istyle] = 27; istyle++;
	if (istyle<NCEL) style[istyle] = 28; istyle++;

	if (drawXTs){
		TChain * ichain = new TChain("t","t");
		ichain->Add(Form("root/ana_%s.i%d.layer%d.root",runname.Data(),Niters,theLayer));
		//ichain->Draw("driftDmc:driftT>>h(600,-25,600,512,-10,10)",Form("chi2<2&&nHitsS>=7&&layerID==%d",theLayer),"COLZ");
		//ichain->Draw("res+theDD:theDT>>h(600,-25,600,512,-10,10)","chi2<0.5&&nHitsS>=6","COLZ");
		ichain->Draw("abs(driftDmc):driftT>>h(600,-25,600,256,0,10)",Form("chi2<2&&nHitsS>=7&&layerID==%d",theLayer),"COLZ");
		//ichain->Draw("abs(res+theDD):theDT>>h(600,-25,600,256,0,10)","chi2<0.5&&nHitsS>=6","COLZ");

		TFile * files[NITERSMAX];
		TF1 * fl[NITERSMAX];
		TF1 * fr[NITERSMAX];
		files[0] = new TFile(Form("info/xt.%s.root",original.Data()));
		fl[0] = (TF1*) files[0]->Get("fl_0");
		fr[0] = (TF1*) files[0]->Get("fr_0");
		fl[0]->SetLineColor(color[0]);
		fr[0]->SetLineColor(color[0]);
		if (!fl[0]||!fr[0]) {printf("Wrong with i%d\n",0); return 0;}
		for (int i = 1; i<=Niters; i++){
			files[i] = new TFile(Form("info/xt.%s.i%d.root",runname.Data(),i));
			fl[i] = (TF1*) files[i]->Get("fl_0");
			fr[i] = (TF1*) files[i]->Get("fr_0");
			if (!fl[i]||!fr[i]) {printf("Wrong with i%d\n",i); return 0;}
			fl[i]->SetLineColor(color[i]);
			fr[i]->SetLineColor(color[i]);
		}
		TGraph * gr_l_ce = (TGraph*) files[Niters]->Get(Form("gr_xt_lce_%d",theLayer));
		TGraph * gr_l_m = (TGraph*) files[Niters]->Get(Form("gr_xt_lm_%d",theLayer));
		TGraph * gr_l_e = (TGraph*) files[Niters]->Get(Form("gr_xt_le_%d",theLayer));
		TGraph * gr_r_ce = (TGraph*) files[Niters]->Get(Form("gr_xt_rce_%d",theLayer));
		TGraph * gr_r_m = (TGraph*) files[Niters]->Get(Form("gr_xt_rm_%d",theLayer));
		TGraph * gr_r_e = (TGraph*) files[Niters]->Get(Form("gr_xt_re_%d",theLayer));
		TGraph * gr_b_ce = (TGraph*) files[Niters]->Get(Form("gr_xt_bce_%d",theLayer));
		TGraph * gr_b_m = (TGraph*) files[Niters]->Get(Form("gr_xt_bm_%d",theLayer));
		TGraph * gr_b_e = (TGraph*) files[Niters]->Get(Form("gr_xt_be_%d",theLayer));
		gr_l_ce->SetMarkerStyle(20);
		gr_l_m->SetMarkerStyle(20);
		gr_l_e->SetMarkerStyle(20);
		gr_r_ce->SetMarkerStyle(20);
		gr_r_m->SetMarkerStyle(20);
		gr_r_e->SetMarkerStyle(20);
		gr_b_ce->SetMarkerStyle(20);
		gr_b_m->SetMarkerStyle(20);
		gr_b_e->SetMarkerStyle(20);
		//gr_l_ce->Draw("PSAME");
		//gr_l_m->Draw("PSAME");
		//gr_l_e->Draw("PSAME");
		//gr_r_ce->Draw("PSAME");
		//gr_r_m->Draw("PSAME");
		//gr_r_e->Draw("PSAME");
		gr_b_ce->Draw("PSAME");
		gr_b_m->Draw("PSAME");
		gr_b_e->Draw("PSAME");
		for (int i = 0; i<=Niters; i++){
			fl[i]->Draw("same");
			fr[i]->Draw("same");
		}
	}

	//====================================Wire calibration============================
	// prepare lists
	double deltaX[NCEL][NLAY][NITERSMAX]; // position difference with the true (designed) value
	double deltaX_max[NCEL];
	double deltaX_min[NCEL];
	double x[NCEL][NLAY]; // x positions of wires by design
	double off[NCEL][NLAY]; // offset set in MC
	TGraph * gr_wp[NCEL][NITERSMAX];

	double Iterations[NITERSMAX]; // iteration IDs
	for (int iter = 0; iter<=Niters; iter++){
		Iterations[iter] = iter;
	}
	bool   changed[NCEL][NLAY]; // is this wire calibrated?
	bool   changedW[NCEL]; // is this wire calibrated?
	for (int wid = 0; wid<NCEL; wid++){
		changedW[wid] = false;
		deltaX_max[wid] = -1e9;
		deltaX_min[wid] = 1e9;
		for (int lid = 0; lid<NLAY; lid++){
			changed[wid][lid] = false;
		}
	}

	// get offsets from MC
	if (isMC){
		TChain * ichain_off = new TChain("t","t");
		ichain_off->Add(Form("info/wire-offset.%d.root",runNo));
		double off_delta;
		int off_lid;
		int off_wid;
		ichain_off->SetBranchAddress("l",&off_lid);
		ichain_off->SetBranchAddress("w",&off_wid);
		ichain_off->SetBranchAddress("d",&off_delta);
		for (int i = 0; i<ichain_off->GetEntries(); i++){
			ichain_off->GetEntry(i);
			if (off_wid<0||off_wid>=NCEL||off_lid<0||off_lid>=NLAY) continue;
			off[off_wid][off_lid] = off_delta;
		}
	}

	// get wire positions
	TChain * ichain_wp = new TChain("t","t");
	ichain_wp->Add(Form("info/wire-position.%s.root",original.Data()));
	double wp_x;
	int wp_lid;
	int wp_wid;
	ichain_wp->SetBranchAddress("xro",&wp_x);
	ichain_wp->SetBranchAddress("l",&wp_lid);
	ichain_wp->SetBranchAddress("w",&wp_wid);
	for (int i = 0; i<ichain_wp->GetEntries(); i++){
		ichain_wp->GetEntry(i);
		if (wp_wid>=NCEL||wp_lid>=NLAY) continue;
		x[wp_wid][wp_lid] = wp_x;
		deltaX[wp_wid][wp_lid][0] = -off[wp_wid][wp_lid];
	}

	// check which wire is changed
	for (int iter = 1; iter<=Niters; iter++){
		ichain_wp = new TChain("t","t");
		ichain_wp->Add(Form("info/wire-position.%s.i%d.root",runname.Data(),iter));
		ichain_wp->SetBranchAddress("xro",&wp_x);
		ichain_wp->SetBranchAddress("l",&wp_lid);
		ichain_wp->SetBranchAddress("w",&wp_wid);
		for (int i = 0; i<ichain_wp->GetEntries(); i++){
			ichain_wp->GetEntry(i);
			if (wp_wid>=NCEL||wp_lid>=NLAY) continue;
			if (wp_x-x[wp_wid][wp_lid]){
				if (!changed[wp_wid][wp_lid]) printf("%d %d %d changed!\n",iter,wp_lid,wp_wid);
				changedW[wp_wid] = true;
				changed[wp_wid][wp_lid] = true;
			}
		}
	}
	// get delta X
	for (int iter = 1; iter<=Niters; iter++){
		ichain_wp = new TChain("t","t");
		ichain_wp->Add(Form("info/wire-position.%s.i%d.root",runname.Data(),iter));
		if (!ichain_wp->GetEntries()){
			printf("Cannot find info/wire-position.%s.i%d.root\n",runname.Data(),iter);
			ichain_wp->Add(Form("info/wire-position.%s.root",original.Data()));
		}
		ichain_wp->SetBranchAddress("xro",&wp_x);
		ichain_wp->SetBranchAddress("l",&wp_lid);
		ichain_wp->SetBranchAddress("w",&wp_wid);
		for (int i = 0; i<ichain_wp->GetEntries(); i++){
			ichain_wp->GetEntry(i);
			if (wp_wid>=NCEL||wp_lid>=NLAY) continue;
			if (changed[wp_wid][wp_lid]){
				deltaX[wp_wid][wp_lid][iter] = wp_x-x[wp_wid][wp_lid]-off[wp_wid][wp_lid];
				printf("deltaX[%d][%d][%d] = %.7e-%.7e-%3e = %.3e\n",wp_wid,wp_lid,iter,wp_x,x[wp_wid][wp_lid],off[wp_wid][wp_lid],deltaX[wp_wid][wp_lid][iter]);
				if (deltaX_max[wp_wid]<deltaX[wp_wid][wp_lid][iter]) deltaX_max[wp_wid] = deltaX[wp_wid][wp_lid][iter];
				if (deltaX_min[wp_wid]>deltaX[wp_wid][wp_lid][iter]) deltaX_min[wp_wid] = deltaX[wp_wid][wp_lid][iter];
			}
		}
	}

	TCanvas * canv_wp = new TCanvas();
	gPad->SetGridx(1);
	gPad->SetGridy(1);
//	TLegend * legend_lid = new TLegend(0.7,0.15,0.9,0.35);
//	TLegend * legend_wid = new TLegend(0.7,0.65,0.9,0.85);
	TLegend * legend_lid = new TLegend(0.1,0.15,0.3,0.35);
	TLegend * legend_wid = new TLegend(0.1,0.65,0.3,0.85);
	bool drawn = false;
	for (int wid = 0; wid<NCEL; wid++){
		if (!changedW[wid]) continue;
		printf("wire %d: %.3e~%.3e\n",wid,deltaX_min[wid],deltaX_max[wid]);
		for (int lid = 1; lid<NLAY; lid++){
			if (!changed[wid][lid]) continue;
			printf("  layer %d changed\n",lid);
			gr_wp[wid][lid] = new TGraph(Niters+1,Iterations,deltaX[wid][lid]);
			gr_wp[wid][lid]->SetTitle(Form("Offset of wire %d in each layer",wid));
			gr_wp[wid][lid]->SetMarkerStyle(style[wid]);gr_wp[wid][lid]->SetMarkerColor(colorLayers[lid]);gr_wp[wid][lid]->SetLineColor(colorLayers[lid]);
			if (wid==4) legend_lid->AddEntry(gr_wp[wid][lid],Form("layer %d",lid),"LP");
			if (!drawn){
				gr_wp[wid][lid]->Draw("APL");
//				gr_wp[wid][lid]->GetYaxis()->SetRangeUser(deltaX_min[wid]-0.05*(deltaX_max[wid]-deltaX_min[wid]),deltaX_max[wid]+0.05*(deltaX_max[wid]-deltaX_min[wid]));
				gr_wp[wid][lid]->GetYaxis()->SetRangeUser(-0.3,0.3);
				//gr_wp[wid][lid]->SetTitle("Offset of each wire w.r.t simulation position");
				gr_wp[wid][lid]->SetTitle("Offset of each wire w.r.t designed position");
				gr_wp[wid][lid]->GetXaxis()->SetTitle("# Iteration");
				gr_wp[wid][lid]->GetYaxis()->SetTitle("Offset [mm]");
				drawn = true;
			}
			else{
				gr_wp[wid][lid]->Draw("PLSAME");
			}
		}
		if (gr_wp[wid][4]) legend_wid->AddEntry(gr_wp[wid][4],Form("wire %d",wid),"LP");
		else if (gr_wp[wid][1]) legend_wid->AddEntry(gr_wp[wid][1],Form("wire %d",wid),"LP");
		else if (gr_wp[wid][8]) legend_wid->AddEntry(gr_wp[wid][8],Form("wire %d",wid),"LP");
	}
	legend_wid->Draw("SAME");
	legend_lid->Draw("SAME");
	canv_wp->SaveAs(Form("wireoff.%s.png",runname.Data()));

//	for (int lid = 1; lid<NLAY; lid++){
//		TLine * l = new TLine(0,off2[lid],Niters,off2[lid]);
//		l->SetLineColor(colorLayers[lid]);
//		canv_wp2->cd();
//		l->Draw("SAME");
//		l = new TLine(0,off3[lid],Niters,off3[lid]);
//		l->SetLineColor(colorLayers[lid]);
//		canv_wp3->cd();
//		l->Draw("SAME");
//		l = new TLine(0,off4[lid],Niters,off4[lid]);
//		l->SetLineColor(colorLayers[lid]);
//		canv_wp4->cd();
//		l->Draw("SAME");
//		l = new TLine(0,off5[lid],Niters,off5[lid]);
//		l->SetLineColor(colorLayers[lid]);
//		canv_wp5->cd();
//		l->Draw("SAME");
//		l = new TLine(0,off6[lid],Niters,off6[lid]);
//		l->SetLineColor(colorLayers[lid]);
//		canv_wp6->cd();
//		l->Draw("SAME");
//	}
}
