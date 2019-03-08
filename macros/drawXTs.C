#include "src/header.hxx"
int drawXTs(){
	std::vector<TString> runNames; std::vector<int> theLayers;
	runNames.push_back("105.0110xt053t400l4.i3"); theLayers.push_back(4);
//	runNames.push_back("105.0110xt053t400l4.i25"); theLayers.push_back(4);
	runNames.push_back("105.0110xt053t400l5.i3"); theLayers.push_back(5);
//	runNames.push_back("105.0110xt053t400l5.i23"); theLayers.push_back(5);
	runNames.push_back("105.0116xt053t400l4.i3"); theLayers.push_back(4);
//	runNames.push_back("105.0116xt053t400l4.i25"); theLayers.push_back(4);
	runNames.push_back("105.1224xt053t400.i3"); theLayers.push_back(4);
//	runNames.push_back("105.1224xt053t400.i30"); theLayers.push_back(4);
	runNames.push_back("105.1224xt053t400l5.i3"); theLayers.push_back(5);
//	runNames.push_back("105.1224xt053t400l5.i26"); theLayers.push_back(5);

	int nRuns = runNames.size();

	int color[NITERSMAX];
	for (int icolor = 0; icolor<NITERSMAX; icolor++){
		color[icolor] = icolor+1;
	}

	TFile * files[NITERSMAX];
	TF1 * fl[NITERSMAX];
	TF1 * fr[NITERSMAX];
	TF1 * fbl[NITERSMAX];
	TF1 * fbr[NITERSMAX];
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
		fl[i] = (TF1*) files[i]->Get(Form("flc_%d",theLayers[i]));
		fr[i] = (TF1*) files[i]->Get(Form("frc_%d",theLayers[i]));
		fbl[i] = (TF1*) files[i]->Get(Form("fl_%d",theLayers[i]));
		fbr[i] = (TF1*) files[i]->Get(Form("fr_%d",theLayers[i]));
		if (!fl[i]||!fr[i]) {printf("Wrong with %s\n",runNames[i].Data()); return 0;}
		fl[i]->SetLineColor(color[i]);
		fr[i]->SetLineColor(color[i]);
		fbl[i]->SetLineColor(color[i]);
		fbr[i]->SetLineColor(color[i]);
		fbl[i]->SetLineStyle(2);
		fbr[i]->SetLineStyle(2);
		fl[i]->SetName(Form("flc_%d_%d",theLayers[i],i));
		fr[i]->SetName(Form("frc_%d_%d",theLayers[i],i));
		fbl[i]->SetName(Form("fl_%d_%d",theLayers[i],i));
		fbr[i]->SetName(Form("fr_%d_%d",theLayers[i],i));
		gr_l_ce[i] = (TGraph*) files[i]->Get(Form("gr_xt_lce_%d",theLayers[i]));
		gr_l_m[i] = (TGraph*) files[i]->Get(Form("gr_xt_lm_%d",theLayers[i]));
		gr_l_e[i] = (TGraph*) files[i]->Get(Form("gr_xt_le_%d",theLayers[i]));
		gr_r_ce[i] = (TGraph*) files[i]->Get(Form("gr_xt_rce_%d",theLayers[i]));
		gr_r_m[i] = (TGraph*) files[i]->Get(Form("gr_xt_rm_%d",theLayers[i]));
		gr_r_e[i] = (TGraph*) files[i]->Get(Form("gr_xt_re_%d",theLayers[i]));
		gr_b_ce[i] = (TGraph*) files[i]->Get(Form("gr_xt_bce_%d",theLayers[i]));
		gr_b_m[i] = (TGraph*) files[i]->Get(Form("gr_xt_bm_%d",theLayers[i]));
		gr_b_e[i] = (TGraph*) files[i]->Get(Form("gr_xt_be_%d",theLayers[i]));
		if (gr_l_ce[i]) gr_l_ce[i]->SetMarkerStyle(7);
		if (gr_l_m[i]) gr_l_m[i]->SetMarkerStyle(7);
		if (gr_l_e[i]) gr_l_e[i]->SetMarkerStyle(7);
		if (gr_r_ce[i]) gr_r_ce[i]->SetMarkerStyle(7);
		if (gr_r_m[i]) gr_r_m[i]->SetMarkerStyle(7);
		if (gr_r_e[i]) gr_r_e[i]->SetMarkerStyle(7);
		if (gr_b_ce[i]) gr_b_ce[i]->SetMarkerStyle(7);
		if (gr_b_m[i]) gr_b_m[i]->SetMarkerStyle(7);
		if (gr_b_e[i]) gr_b_e[i]->SetMarkerStyle(7);
		if (gr_l_ce[i]) gr_l_ce[i]->SetMarkerColor(color[i]);
		if (gr_l_m[i]) gr_l_m[i]->SetMarkerColor(color[i]);
		if (gr_l_e[i]) gr_l_e[i]->SetMarkerColor(color[i]);
		if (gr_r_ce[i]) gr_r_ce[i]->SetMarkerColor(color[i]);
		if (gr_r_m[i]) gr_r_m[i]->SetMarkerColor(color[i]);
		if (gr_r_e[i]) gr_r_e[i]->SetMarkerColor(color[i]);
		if (gr_b_ce[i]) gr_b_ce[i]->SetMarkerColor(color[i]);
		if (gr_b_m[i]) gr_b_m[i]->SetMarkerColor(color[i]);
		if (gr_b_e[i]) gr_b_e[i]->SetMarkerColor(color[i]);
	}
	TH2D * h2_xt = (TH2D*) files[0]->Get(Form("h2_xt_%d",theLayers[0]));
	TH2D * h2_xtn = (TH2D*) files[0]->Get(Form("h2_xtn_%d",theLayers[0]));
	int nBins = h2_xt->GetXaxis()->GetNbins();
	double xMin = h2_xt->GetXaxis()->GetXmin();
	double xMax = h2_xt->GetXaxis()->GetXmax();

    // Compare different runs, for l/r combined xts
    TCanvas * canv = new TCanvas("canv1","canv1",768,768);
    canv->Divide(1,2);
    canv->cd(1);
	h2_xtn->Draw("COLZ");
	TLegend * legend = new TLegend(0.7,0.7,1,1);
	for (int i = 0; i<nRuns; i++){
		fbr[i]->Draw("same");
//		if (gr_b_ce[i]) gr_b_ce[i]->Draw("PSAME");
//		if (gr_b_m[i]) gr_b_m[i]->Draw("PSAME");
//		if (gr_b_e[i]) gr_b_e[i]->Draw("PSAME");
		legend->AddEntry(fr[i],runNames[i]);
	}
	legend->Draw("SAME");

    canv->cd(2);
	TH2D * h2_diff = new TH2D("diffb","Difference among XTs from both sides",nBins,xMin,xMax,256,-0.5,0.5);
    h2_diff->Draw();
	legend = new TLegend(0.7,0.7,1,1);
	for (int i = 1; i<nRuns; i++){
	    TF1 * ftemp = new TF1(Form("fbdiff_%d_0",i),Form("%s-(%s)",fbr[i]->GetName(),fbr[0]->GetName()),-20,600);
	    ftemp->SetLineColor(color[i]);
	    ftemp->Draw("SAME");
        legend->AddEntry(ftemp,runNames[i]);
    }
	legend->Draw("SAME");

    // Compare different runs, for l/r separated xts
    canv = new TCanvas("canv2","canv2",768,1024);
    canv->Divide(1,3);
    canv->cd(1);
	h2_xt->Draw("COLZ");
	for (int i = 0; i<nRuns; i++){
		fl[i]->Draw("same");
		fr[i]->Draw("same");
//		if (gr_r_ce[i]) gr_r_ce[i]->Draw("PSAME");
//		if (gr_r_m[i]) gr_r_m[i]->Draw("PSAME");
//		if (gr_r_e[i]) gr_r_e[i]->Draw("PSAME");
//		if (gr_l_ce[i]) gr_l_ce[i]->Draw("PSAME");
//		if (gr_l_m[i]) gr_l_m[i]->Draw("PSAME");
//		if (gr_l_e[i]) gr_l_e[i]->Draw("PSAME");
	}
	legend = new TLegend(0.7,0.7,1,1);
	for ( int i = 0; i<nRuns; i++){
		legend->AddEntry(fr[i],runNames[i]);
	}
	legend->Draw("SAME");

    canv->cd(2);
	h2_diff = new TH2D("diffl","Difference among XTs from left sides",nBins,xMin,xMax,256,-0.5,0.5);
    h2_diff->Draw();
	legend = new TLegend(0.7,0.7,1,1);
	for (int i = 1; i<nRuns; i++){
	    TF1 * fltemp = new TF1(Form("fldiff_%d_0",i),Form("%s-(%s)",fl[i]->GetName(),fl[0]->GetName()),-20,600);
	    fltemp->SetLineColor(color[i]);
	    fltemp->Draw("SAME");
        legend->AddEntry(fltemp,runNames[i]);
    }
	legend->Draw("SAME");

    canv->cd(3);
	h2_diff = new TH2D("diffr","Difference among XTs from right sides",nBins,xMin,xMax,256,-0.5,0.5);
    h2_diff->Draw();
	legend = new TLegend(0.7,0.7,1,1);
	for (int i = 1; i<nRuns; i++){
	    TF1 * frtemp = new TF1(Form("frdiff_%d_0",i),Form("%s-(%s)",fr[i]->GetName(),fr[0]->GetName()),-20,600);
	    frtemp->SetLineColor(color[i]);
	    frtemp->Draw("SAME");
        legend->AddEntry(frtemp,runNames[i]);
    }
	legend->Draw("SAME");

    // Compare in the first run, for l/r separated and combined xts
    for (int i = 0; i<nRuns; i++){
        canv = new TCanvas(Form("canv_lrb_%d",i),"canv",768,768);
        h2_diff = new TH2D(Form("diff_lrb_%d",i),Form("Difference among XTs from separated sides for %s",runNames[i].Data()),nBins,xMin,xMax,256,-0.5,0.5);
        canv->Divide(1,2);
        canv->cd(1);
        h2_xt->Draw("COLZ");
        fl[i]->Draw("same");
        fr[i]->Draw("same");
        fbl[i]->Draw("same");
        fbr[i]->Draw("same");
        if (gr_r_ce[i]) gr_r_ce[i]->Draw("PSAME");
        if (gr_r_m[i]) gr_r_m[i]->Draw("PSAME");
        if (gr_r_e[i]) gr_r_e[i]->Draw("PSAME");
        if (gr_l_ce[i]) gr_l_ce[i]->Draw("PSAME");
        if (gr_l_m[i]) gr_l_m[i]->Draw("PSAME");
        if (gr_l_e[i]) gr_l_e[i]->Draw("PSAME");
        //if (gr_b_ce[i]) gr_b_ce[i]->Draw("PSAME");
        //if (gr_b_m[i]) gr_b_m[i]->Draw("PSAME");
        //if (gr_b_e[i]) gr_b_e[i]->Draw("PSAME");
        legend = new TLegend(0.7,0.7,1,1);
        legend->AddEntry(fr[i],runNames[i]);
        legend->Draw("SAME");
        canv->cd(2);
        h2_diff->Draw();
        TF1 * fltemp = new TF1(Form("flbdiff_%d",i),Form("%s-(%s)",fl[i]->GetName(),fbl[i]->GetName()),-20,600);
        TF1 * frtemp = new TF1(Form("frbdiff_%d",i),Form("%s-(%s)",fr[i]->GetName(),fbr[i]->GetName()),-20,600);
        fltemp->SetLineColor(kBlue);
        frtemp->SetLineColor(kRed);
        fltemp->Draw("SAME");
        frtemp->Draw("SAME");
    }
}
