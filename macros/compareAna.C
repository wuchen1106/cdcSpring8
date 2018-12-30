{
    TChain * ichain0 = new TChain("t","t");
    TChain * ichain1 = new TChain("t","t");
    ichain0->Add("root/ana/ana_105.1224xt053t400.i4.layer4.root");
    ichain1->Add("root/ana/ana_105.1224xt053t400l5.i4.layer4.root");
    int isGood0;
    int nHitsS0;
    int theWid0;
    double theFD0;
    double chi2p0;
    double slz0;
    double inz0;
    double slx0;
    double inx0;
    int isGood1;
    int nHitsS1;
    int theWid1;
    double theFD1;
    double chi2p1;
    double slz1;
    double inz1;
    double slx1;
    double inx1;
    ichain0->SetBranchAddress("isGood",&isGood0);
    ichain0->SetBranchAddress("nHitsS",&nHitsS0);
    ichain0->SetBranchAddress("theWid0",&theWid0);
    ichain0->SetBranchAddress("theFD0",&theFD0);
    ichain0->SetBranchAddress("chi2p",&chi2p0);
    ichain0->SetBranchAddress("slz",&slz0);
    ichain0->SetBranchAddress("slx",&slx0);
    ichain0->SetBranchAddress("inz",&inz0);
    ichain0->SetBranchAddress("inx",&inx0);
    ichain1->SetBranchAddress("isGood",&isGood1);
    ichain1->SetBranchAddress("nHitsS",&nHitsS1);
    ichain1->SetBranchAddress("theWid0",&theWid1);
    ichain1->SetBranchAddress("theFD0",&theFD1);
    ichain1->SetBranchAddress("chi2p",&chi2p1);
    ichain1->SetBranchAddress("slz",&slz1);
    ichain1->SetBranchAddress("slx",&slx1);
    ichain1->SetBranchAddress("inz",&inz1);
    ichain1->SetBranchAddress("inx",&inx1);

    TFile * ofile = new TFile("compareAna.root","RECREATE");
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

    double zl70 = 0;
    double zl71 = 0;
    double xl70 = 0;
    double xl71 = 0;

    TTree * otree = new TTree("t","t");
    otree->Branch("isGood0",&isGood0);
    otree->Branch("nHitsS0",&nHitsS0);
    otree->Branch("theWid0",&theWid0);
    otree->Branch("theFD0",&theFD0);
    otree->Branch("chi2p0",&chi2p0);
    otree->Branch("slz0",&slz0);
    otree->Branch("slx0",&slx0);
    otree->Branch("inz0",&inz0);
    otree->Branch("inx0",&inx0);
    otree->Branch("xl70",&xl70);
    otree->Branch("zl70",&zl70);
    otree->Branch("isGood1",&isGood1);
    otree->Branch("theFD1",&theFD1);
    otree->Branch("nHitsS1",&nHitsS1);
    otree->Branch("theWid1",&theWid1);
    otree->Branch("chi2p1",&chi2p1);
    otree->Branch("slz1",&slz1);
    otree->Branch("slx1",&slx1);
    otree->Branch("inz1",&inz1);
    otree->Branch("inx1",&inx1);
    otree->Branch("xl71",&xl71);
    otree->Branch("zl71",&zl71);

    for (int iEntry = 0; iEntry<ichain0->GetEntries(); iEntry++){
        ichain0->GetEntry(iEntry);
        ichain1->GetEntry(iEntry);
        double zl70 = inz0+slz0*(yl7-sciYup);
        double zl71 = inz1+slz1*(yl7-sciYup);
        double xl70 = inx0+slx0*(yl7-sciYup);
        double xl71 = inx1+slx1*(yl7-sciYup);
        if (nHitsS0>=7&&nHitsS1>=7&&abs(theFD1)<6&&abs(theFD1)>2&&isGood1){
            hchi2p->Fill(chi2p0,chi2p1);
            hslz->Fill(slz0,slz1-slz0);
            hslx->Fill(slx0,slx1-slx0);
            hinz->Fill(inz0,inz1-inz0);
            hzl7->Fill(zl70,zl71-zl70);
            hxl7->Fill(xl70,xl71-xl70);
            hinx->Fill(inx0,inx1-inx0);
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
