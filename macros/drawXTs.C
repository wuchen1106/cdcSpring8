#include "src/header.h"
int drawXTs(){
	std::vector<TString> runNames;
	runNames.push_back("1010.0226.sm10a30n35.i15");
	runNames.push_back("1012.0226.sm10a30n35.i15");
	runNames.push_back("1013.0226.sm10a30n35.i15");
	runNames.push_back("1014.0226.sm10a30n35.i15");
	runNames.push_back("1035.0226.sm10a30n35.i15");

	int nRuns = runNames.size();
	int theLayer = 4;

	int color[NITERSMAX];
	for (int icolor = 0; icolor<NITERSMAX; icolor++){
		color[icolor] = icolor+1;
	}

	TFile * files[NITERSMAX];
	TF1 * fl[NITERSMAX];
	TF1 * fr[NITERSMAX];
	TGraph * gr_l_ce[NITERSMAX];
	TGraph * gr_l_m[NITERSMAX];
	TGraph * gr_l_e[NITERSMAX];
	TGraph * gr_r_ce[NITERSMAX];
	TGraph * gr_r_m[NITERSMAX];
	TGraph * gr_r_e[NITERSMAX];
	TGraph * gr_b_ce[NITERSMAX];
	TGraph * gr_b_m[NITERSMAX];
	TGraph * gr_b_e[NITERSMAX];
	for (int i = 0; i<nRuns; i++){
		files[i] = new TFile(Form("info/xt.%s.root",runNames[i].Data()));
		if (!files[i]) {printf("Cannot find info/xt.%s.root",runNames[i].Data()); return 0;}
		fl[i] = (TF1*) files[i]->Get("fl_0");
		fr[i] = (TF1*) files[i]->Get("fr_0");
		if (!fl[i]||!fr[i]) {printf("Wrong with i%d\n",i); return 0;}
		fl[i]->SetLineColor(color[i]);
		fr[i]->SetLineColor(color[i]);
		gr_l_ce[i] = (TGraph*) files[i]->Get(Form("gr_xt_lce_%d",theLayer));
		gr_l_m[i] = (TGraph*) files[i]->Get(Form("gr_xt_lm_%d",theLayer));
		gr_l_e[i] = (TGraph*) files[i]->Get(Form("gr_xt_le_%d",theLayer));
		gr_r_ce[i] = (TGraph*) files[i]->Get(Form("gr_xt_rce_%d",theLayer));
		gr_r_m[i] = (TGraph*) files[i]->Get(Form("gr_xt_rm_%d",theLayer));
		gr_r_e[i] = (TGraph*) files[i]->Get(Form("gr_xt_re_%d",theLayer));
		gr_b_ce[i] = (TGraph*) files[i]->Get(Form("gr_xt_bce_%d",theLayer));
		gr_b_m[i] = (TGraph*) files[i]->Get(Form("gr_xt_bm_%d",theLayer));
		gr_b_e[i] = (TGraph*) files[i]->Get(Form("gr_xt_be_%d",theLayer));
		if (!gr_l_ce[i]||!gr_l_m[i]||!gr_l_e[i]||!gr_r_ce[i]||!gr_r_m[i]||!gr_r_e[i]||!gr_b_ce[i]||!gr_b_m[i]||!gr_b_e[i]) {printf("Wrong with i%d\n",i); return 0;}
		gr_l_ce[i]->SetMarkerStyle(20);
		gr_l_m[i]->SetMarkerStyle(20);
		gr_l_e[i]->SetMarkerStyle(20);
		gr_r_ce[i]->SetMarkerStyle(20);
		gr_r_m[i]->SetMarkerStyle(20);
		gr_r_e[i]->SetMarkerStyle(20);
		gr_b_ce[i]->SetMarkerStyle(20);
		gr_b_m[i]->SetMarkerStyle(20);
		gr_b_e[i]->SetMarkerStyle(20);
		gr_l_ce[i]->SetMarkerColor(color[i]);
		gr_l_m[i]->SetMarkerColor(color[i]);
		gr_l_e[i]->SetMarkerColor(color[i]);
		gr_r_ce[i]->SetMarkerColor(color[i]);
		gr_r_m[i]->SetMarkerColor(color[i]);
		gr_r_e[i]->SetMarkerColor(color[i]);
		gr_b_ce[i]->SetMarkerColor(color[i]);
		gr_b_m[i]->SetMarkerColor(color[i]);
		gr_b_e[i]->SetMarkerColor(color[i]);
	}
	TH2D * h2_xt = (TH2D*) files[0]->Get(Form("h2_xt_%d",theLayer));
	TH2D * h2_xtn = (TH2D*) files[0]->Get(Form("h2_xtn_%d",theLayer));

	h2_xtn->Draw("COLZ");
	for (int i = 0; i<nRuns; i++){
//		fl[i]->Draw("same");
		fr[i]->Draw("same");
		gr_b_ce[i]->Draw("PSAME");
		gr_b_m[i]->Draw("PSAME");
		gr_b_e[i]->Draw("PSAME");
	}
	TLegend * legend = new TLegend(0.7,0.7,0.9,0.9);
	for ( int i = 0; i<nRuns; i++){
		legend->AddEntry(fr[i],runNames[i]);
	}
	legend->Draw("SAME");
}
