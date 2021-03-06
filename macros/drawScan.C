{
    double m_maxEff[3];
    double m_minEff[3];
    double m_maxRes[3];
    double m_minRes[3];
    m_maxEff[0] = 1.05;
    m_minEff[0] = 0.75;
    m_maxRes[0] = 0.6;
    m_minRes[0] = 0.;
    m_maxEff[1] = 1.02;
    m_minEff[1] = 0.86;
    m_maxRes[1] = 0.35;
    m_minRes[1] = 0.1;
    m_maxEff[2] = 1.0;
    m_minEff[2] = 0.8;
    m_maxRes[2] = 0.6;
    m_minRes[2] = 0.2;

	TChain * ichain = new TChain("t","t");
	ichain->Add("result/ResEff/anabest.0614.resxAV.root");
	double gg;
	double eff;
	double rms;
	double res;
	double res0;
	double etrack;
	double etrack0;
	double ggerr = 0;
	double efferr = 0;
	double rmserr = 0;
	int hv;
	int gas;
	int N = 3e5;
	ichain->SetBranchAddress("averageGG",&gg);
	ichain->SetBranchAddress("averageGGErr",&ggerr);
//	ichain->SetBranchAddress("averageEffX",&eff);
	ichain->SetBranchAddress("averageEff3sX",&eff);
	ichain->SetBranchAddress("averageEff3sXErr",&efferr);
	ichain->SetBranchAddress("averageRMSX",&rms);
	ichain->SetBranchAddress("averageRMSXErr",&rmserr);
	ichain->SetBranchAddress("avresInix",&res);
	ichain->SetBranchAddress("avresInix0",&res0);
	ichain->SetBranchAddress("avEtrackx",&etrack);
	ichain->SetBranchAddress("avEtrackx0",&etrack0);
	ichain->SetBranchAddress("HV",&hv);
	ichain->SetBranchAddress("gasID",&gas);
	int nEntries = ichain->GetEntries();

	int maxNp = 10;
	int nPoints[3] = {0};
	double maxRMS = 1;
	TString gasName[3];
	gasName[0]="He-C_{2}H_{6} (50/50)";
	gasName[1]="He-iC_{4}H_{10} (90/10)";
	gasName[2]="He-CH_{4} (80/20)";
	TGraphErrors * gr_eff[3];
	TGraphErrors * gr_rms[3];
	TGraphErrors * gr_res[3];
	TGraphErrors * gr_etrack[3];
	TGraphErrors * gr_gg[3];
	for (int i = 0; i<3; i++){
		gr_eff[i] = new TGraphErrors();
		gr_rms[i] = new TGraphErrors();
		gr_res[i] = new TGraphErrors();
		gr_etrack[i] = new TGraphErrors();
		gr_gg[i] = new TGraphErrors();
		gr_eff[i]->Set(maxNp);
		gr_rms[i]->Set(maxNp);
		gr_res[i]->Set(maxNp);
		gr_etrack[i]->Set(maxNp);
		gr_gg[i]->Set(maxNp);
	}
	for (int iEntry = 0; iEntry<nEntries; iEntry++){
		ichain->GetEntry(iEntry);
		if (!eff) continue;
		if (gas>=3||gas<0) continue;
		if (gas==2&&hv<1900) continue;
		double hverr = 25;
		if (gas==0 && hv>=2450) hverr = 50;
		else if (gas==1 && hv>=1800) hverr = 15;
		else if (gas==2 && hv>=2000) hverr = 5;
		hverr = 0; // don't draw HV error bar
		double efferr_t = efferr;
		if (!efferr_t) efferr_t = sqrt(eff*(1-eff)/N);
		gr_eff[gas]->SetPoint(nPoints[gas],hv,eff);
		gr_eff[gas]->SetPointError(nPoints[gas],hverr,efferr_t);
		gr_rms[gas]->SetPoint(nPoints[gas],hv,(rms-m_minRes[gas])*(m_maxEff[gas]-m_minEff[gas])/(m_maxRes[gas]-m_minRes[gas])+m_minEff[gas]);
		gr_rms[gas]->SetPointError(nPoints[gas],hverr,rmserr*(m_maxEff[gas]-m_minEff[gas])/(m_maxRes[gas]-m_minRes[gas]));
		gr_res[gas]->SetPoint(nPoints[gas],hv,(res-m_minRes[gas])*(m_maxEff[gas]-m_minEff[gas])/(m_maxRes[gas]-m_minRes[gas])+m_minEff[gas]);
		gr_res[gas]->SetPointError(nPoints[gas],hverr,fabs(res-res0)*(m_maxEff[gas]-m_minEff[gas])/(m_maxRes[gas]-m_minRes[gas]));
		gr_etrack[gas]->SetPoint(nPoints[gas],hv,(etrack-m_minRes[gas])*(m_maxEff[gas]-m_minEff[gas])/(m_maxRes[gas]-m_minRes[gas])+m_minEff[gas]);
		gr_etrack[gas]->SetPointError(nPoints[gas],hverr,fabs(etrack-etrack0)*(m_maxEff[gas]-m_minEff[gas])/(m_maxRes[gas]-m_minRes[gas]));
		gr_gg[gas]->SetPoint(nPoints[gas],hv,gg);
		gr_gg[gas]->SetPointError(nPoints[gas],hverr,ggerr/sqrt(nEntries));
		nPoints[gas]++;
	}

	TCanvas * canvas = 0;
	for (int i = 0; i<3; i++){
		gr_eff[i]->Set(nPoints[i]);
		gr_rms[i]->Set(nPoints[i]);
		gr_res[i]->Set(nPoints[i]);
		gr_etrack[i]->Set(nPoints[i]);
		canvas = new TCanvas();
		gr_eff[i]->SetMarkerStyle(4);
		gr_rms[i]->SetMarkerStyle(26);
		gr_res[i]->SetMarkerStyle(26);
		gr_etrack[i]->SetMarkerStyle(26);

		gr_eff[i]->GetYaxis()->SetRangeUser(m_minEff[i],m_maxEff[i]);
		//gr_eff[i]->GetYaxis()->SetTitle("Hit efficiency (500 um cut)");
		gr_eff[i]->GetYaxis()->SetTitle("Hit efficiency");
		gr_eff[i]->GetXaxis()->SetTitle("HV [V]");
		gr_eff[i]->SetTitle(gasName[i]);
		gr_eff[i]->Draw("AP");

		double xmax = gr_eff[i]->GetXaxis()->GetXmax();
		double ymax = m_maxEff[i];
		double ymin = m_minEff[i];
		TGaxis * ax = new TGaxis(xmax,ymin,xmax,ymax,m_minRes[i],m_maxRes[i],510,"+L");
		ax->SetTitle("#sigma [mm]");
		int opt;
		TAxis * yax = gr_eff[i]->GetYaxis();
		opt = yax->GetTitleFont(); ax->SetTitleFont(opt);
		opt = yax->GetTickLength(); ax->SetTickLength(opt);
		opt = yax->GetLabelFont(); ax->SetLabelFont(opt);
		double optd;
		optd = yax->GetTitleSize(); ax->SetTitleSize(optd);
		optd = yax->GetLabelSize(); ax->SetLabelSize(optd);
		ax->Draw();

		//TLegend * leg = new TLegend(0.65,0.4,0.85,0.6);
		//leg->AddEntry(gr_eff[i],"Hit Efficiency","P");
		//leg->AddEntry(gr_rms[i],"Spatial Resolution","P");
		//leg->Draw();

		gr_res[i]->Draw("PSAME");
		canvas->SaveAs(Form("scan_res_eff_%d.eps",i));

		gr_res[i]->SetMarkerColor(kRed);
		gr_etrack[i]->SetMarkerColor(kBlue);
		gr_res[i]->SetLineColor(kRed);
		gr_etrack[i]->SetLineColor(kBlue);
		gr_rms[i]->Draw("PSAME");
		gr_etrack[i]->Draw("PSAME");
		canvas->SaveAs(Form("scan_resall_eff_%d.eps",i));
	}
	canvas = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetLogy(1);
	int style[3] = {20,21,22};
	TH2D * hbk = new TH2D("hbk","",1024,1550,2650,1024,1e3,2e5);
	hbk->GetYaxis()->SetTitle("Gas gain");
	hbk->GetXaxis()->SetTitle("HV [V]");
	hbk->Draw();
	//TLegend * legend = new TLegend(0.55,0.15,0.875,0.425);
	for (int i = 0; i<3; i++){
		gr_gg[i]->Set(nPoints[i]);
		gr_gg[i]->SetMarkerStyle(style[i]);
		gr_gg[i]->GetYaxis()->SetTitle("Gas gain");
		gr_gg[i]->GetXaxis()->SetTitle("HV [V]");
		gr_gg[i]->SetTitle(gasName[i]);
		gr_gg[i]->Draw("PSAME");
		//legend->AddEntry(gr_gg[i],gasName[i],"P");
	}
	//legend->Draw("SAME");
	canvas->SaveAs("scan_gg.eps");
}
