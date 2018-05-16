{
	int color[6];
	int icolor = 0;
	color[icolor] = kBlack;icolor++; 
	color[icolor] = kGreen;icolor++;
	color[icolor] = kCyan;icolor++;
	color[icolor] = kBlue;icolor++;
	color[icolor] = kMagenta;icolor++;
	color[icolor] = kRed;icolor++;

	TFile * ifile[6];
	TGraph * gr_resIni[6];
	TGraph * gr_resIniOld[6];
	for (int i = 0; i<6; i++){
		ifile[i] = new TFile(Form("info/res.2030.0426resxAV.i%d.root",i));
		gr_resIni[i] = (TGraph*) ifile[i]->Get("gr_resIni");
		gr_resIniOld[i] = (TGraph*) ifile[i]->Get("gr_resIniOldx");
		gr_resIni[i]->SetMarkerStyle(20);
		gr_resIni[i]->SetMarkerColor(color[i]);
		gr_resIni[i]->SetLineColor(color[i]);
		gr_resIniOld[i]->SetMarkerStyle(4);
		gr_resIniOld[i]->SetLineStyle(2);
		gr_resIniOld[i]->SetMarkerColor(color[i]);
		gr_resIniOld[i]->SetLineColor(color[i]);
		if (i==0) gr_resIni[i]->Draw("APL");
		else  gr_resIni[i]->Draw("PLSAME");
		gr_resIniOld[i]->Draw("PLSAME");
	}
}
