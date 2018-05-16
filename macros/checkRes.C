{
    TString runname = "0424resx";
//    TString runname = "0423resxAV";
    TString histname = "gr_Etrackx";
//    TString histname = "gr_resTotnew";
//    TString histname = "gr_resIni";
    TFile * file0 = new TFile("info/res.100014."+runname+".i0.root");
    TGraph * gr0 = (TGraph*) file0->Get(histname);
    TFile * file1 = new TFile("info/res.100014."+runname+".i1.root");
    TGraph * gr1 = (TGraph*) file1->Get(histname);
    TFile * file2 = new TFile("info/res.100014."+runname+".i2.root");
    TGraph * gr2 = (TGraph*) file2->Get(histname);
    TFile * file3 = new TFile("info/res.100014."+runname+".i3.root");
    TGraph * gr3 = (TGraph*) file3->Get(histname);
    TFile * file4 = new TFile("info/res.100014."+runname+".i4.root");
    TGraph * gr4 = (TGraph*) file4->Get(histname);
    TFile * file5 = new TFile("info/res.100014."+runname+".i5.root");
    TGraph * gr5 = (TGraph*) file5->Get(histname);

    if (gr0) gr0->SetLineColor(kBlack);
    if (gr0) gr0->SetMarkerColor(kBlack);
    if (gr0) gr0->SetMarkerStyle(20);
    if (gr1) gr1->SetLineColor(kBlue);
    if (gr1) gr1->SetMarkerColor(kBlue);
    if (gr1) gr1->SetMarkerStyle(20);
    if (gr2) gr2->SetLineColor(kGreen);
    if (gr2) gr2->SetMarkerColor(kGreen);
    if (gr2) gr2->SetMarkerStyle(20);
    if (gr3) gr3->SetLineColor(kYellow);
    if (gr3) gr3->SetMarkerColor(kYellow);
    if (gr3) gr3->SetMarkerStyle(20);
    if (gr4) gr4->SetLineColor(kMagenta);
    if (gr4) gr4->SetMarkerColor(kMagenta);
    if (gr4) gr4->SetMarkerStyle(20);
    if (gr5) gr5->SetLineColor(kRed);
    if (gr5) gr5->SetMarkerColor(kRed);
    if (gr5) gr5->SetMarkerStyle(20);

    if (gr0) gr0->GetYaxis()->SetRangeUser(0,0.4);
    if (gr0) gr0->GetYaxis()->SetTitle("Tracking Error [mm]");
    if (gr0) gr0->GetXaxis()->SetTitle("DOCA [mm]");
    if (gr0) gr0->Draw("APL");
    if (gr1) gr1->Draw("PLSAME");
    if (gr2) gr2->Draw("PLSAME");
    if (gr3) gr3->Draw("PLSAME");
    if (gr4) gr4->Draw("PLSAME");
    if (gr5) gr5->Draw("PLSAME");

    TLegend * l = new TLegend(0.2,0.2,0.5,0.5);
    if (gr0) l->AddEntry(gr0,"Start","LP");
    if (gr1) l->AddEntry(gr1,"Iteration 1","LP");
    if (gr2) l->AddEntry(gr2,"Iteration 2","LP");
    if (gr3) l->AddEntry(gr3,"Iteration 3","LP");
    if (gr4) l->AddEntry(gr4,"Iteration 4","LP");
    if (gr5) l->AddEntry(gr5,"Iteration 5","LP");
    l->Draw("SAME");
}
