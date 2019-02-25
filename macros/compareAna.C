{
    const int NF = 10;
    TChain * ichain[NF] = {0};
    int isGood[NF] = {0};
    int nHitsS[NF] = {0};
    int theWid[NF] = {0};
    double theFD[NF] = {0};
    double theDD[NF] = {0};
    double chi2p[NF] = {0};
    double slz[NF] = {0};
    double inz[NF] = {0};
    double slx[NF] = {0};
    double inx[NF] = {0};

    double zl7[NF] = {0};
    double xl7[NF] = {0};

    TFile * ofile = new TFile("compareAna.root","RECREATE");
    TTree * otree = new TTree("t","t");

    std::vector<TString> filenames;
    std::vector<TString> tagnames;
    filenames.push_back("105.1224xt053t400.i4.layer4"); tagnames.push_back("l4");
    filenames.push_back("105.1224xt053t400.i4.layer5"); tagnames.push_back("l5");
    filenames.push_back("105.1224xt053t400.i4.layer7"); tagnames.push_back("l7");

    for (int i = 0; i<filenames.size(); i++){
        ichain[i] = new TChain("t","t");
        ichain[i]->Add(Form("root/ana/ana_%s.root",filenames[i].Data()));
        ichain[i]->SetBranchAddress("isGood",&(isGood[i]));
        ichain[i]->SetBranchAddress("nHitsS",&(nHitsS[i]));
        ichain[i]->SetBranchAddress("theWid0",&(theWid[i]));
        ichain[i]->SetBranchAddress("theFD0",&(theFD[i]));
        ichain[i]->SetBranchAddress("theDD0",&(theDD[i]));
        ichain[i]->SetBranchAddress("chi2p",&(chi2p[i]));
        ichain[i]->SetBranchAddress("slz",&(slz[i]));
        ichain[i]->SetBranchAddress("slx",&(slx[i]));
        ichain[i]->SetBranchAddress("inz",&(inz[i]));
        ichain[i]->SetBranchAddress("inx",&(inx[i]));

        otree->Branch(Form("isGood%s",tagnames[i].Data()),&(isGood[i]));
        otree->Branch(Form("nHitsS%s",tagnames[i].Data()),&(nHitsS[i]));
        otree->Branch(Form("theWid%s",tagnames[i].Data()),&(theWid[i]));
        otree->Branch(Form("theFD%s",tagnames[i].Data()),&(theFD[i]));
        otree->Branch(Form("theDD%s",tagnames[i].Data()),&(theDD[i]));
        otree->Branch(Form("chi2p%s",tagnames[i].Data()),&(chi2p[i]));
        otree->Branch(Form("slz%s",tagnames[i].Data()),&(slz[i]));
        otree->Branch(Form("slx%s",tagnames[i].Data()),&(slx[i]));
        otree->Branch(Form("inz%s",tagnames[i].Data()),&(inz[i]));
        otree->Branch(Form("inx%s",tagnames[i].Data()),&(inx[i]));
        otree->Branch(Form("xl7%s",tagnames[i].Data()),&(xl7[i]));
        otree->Branch(Form("zl7%s",tagnames[i].Data()),&(zl7[i]));
    }

    TH2D * hchi2p = new TH2D("hchi2p","hchi2p",128,0,1,128,0,1);
    TH2D * hslz = new TH2D("hslz","hslz",128,-0.15,0.15,128,-0.05,0.05);
    TH2D * hslx = new TH2D("hslx","hslx",128,-0.1,0.1,128,-0.005,0.005);
    TH2D * hinz = new TH2D("hinz","hinz",128,-200,200,128,-25,25);
    TH2D * hzl7 = new TH2D("hzl7","hzl7",128,-200,200,128,-5,5);
    TH2D * hinx = new TH2D("hinx","hinx",128,-50,50,128,-1,1);
    TH2D * hxl7 = new TH2D("hxl7","hxl7",128,-50,50,128,-1,1);

    double chamberHH = 170.05/2; // mm
    double chamberCY = 572; // mm
    double sciYup = chamberCY+chamberHH+180; // mm
    double yl7 = 624;

    for (int iEntry = 0; iEntry<ichain[0]->GetEntries(); iEntry++){
        for (int i = 0; i<filenames.size(); i++){
            ichain[i]->GetEntry(iEntry);
            double zl7[i] = inz[i]+slz[i]*(yl7-sciYup);
            double xl7[i] = inx[i]+slx[i]*(yl7-sciYup);
        }
        if (nHitsS[2]>=7&&theFD[2]<-3&&theFD[2]>-5&&isGood[2]){ // file 2 is layer 7
            // compare with file 0 which is layer 4
            hchi2p->Fill(chi2p[0],chi2p[2]);
            hslz->Fill(slz[0],slz[2]-slz[0]);
            hslx->Fill(slx[0],slx[2]-slx[0]);
            hinz->Fill(inz[0],inz[2]-inz[0]);
            hzl7->Fill(zl7[0],zl7[2]-zl7[0]);
            hxl7->Fill(xl7[0],xl7[2]-xl7[0]);
            hinx->Fill(inx[0],inx[2]-inx[0]);
        }
        otree->Fill();
    }

    new TCanvas();
    hchi2p->Draw("COLZ");
    new TCanvas();
    hslz->Draw("COLZ");
    new TCanvas();
    hslx->Draw("COLZ");
    new TCanvas();
    hinz->Draw("COLZ");
    new TCanvas();
    hzl7->Draw("COLZ");
    new TCanvas();
    hinx->Draw("COLZ");
    new TCanvas();
    hxl7->Draw("COLZ");

    hchi2p->Write();
    hslz->Write();
    hslx->Write();
    hinz->Write();
    hzl7->Write();
    hinx->Write();
    hxl7->Write();
    otree->Write();
}
