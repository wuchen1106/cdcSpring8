{
    int NLAY = 9;

    double Xmin = -45;
    double Xmax = 45;
    double Zmin = -150;
    double Zmax = 150;
    double Ymin = 306.975;
    double Ymax = 837.025;
    double Yref = Ymax;

    double FDmin[9] = {0,0,-3,-3,0,0,-3,-3,0};
    double FDmax[9] = {3,3, 0, 0,3,3, 0, 0,3};
    int    WireID[9] = {0,4,4,4,4,5,5,6,5};
    int    color[9] = {1,2,3,4,5,6,7,8,9};

    int NBinY = 256;
    int NBinXZ = 256;
    TH2D * hYZ[NLAY];
    TH2D * hYX[NLAY];
    TH3D * hYXZ = new TH3D("hYXZ","",NBinXZ,Xmin,Xmax,NBinXZ,Zmin,Zmax,NBinY,Ymin,Ymax);
    hYXZ->GetXaxis()->SetTitle("X");
    hYXZ->GetYaxis()->SetTitle("Y");
    hYXZ->GetZaxis()->SetTitle("Z");
    for (int i = 0; i<NLAY; i++){
        hYZ[i] = new TH2D(Form("hYZ%d",i),Form("Beam in Layer %d with DOCA %.0f~%.0f mm",i,FDmin[i],FDmax[i]),NBinXZ,Zmin,Zmax,NBinY,Ymin,Ymax);
        hYX[i] = new TH2D(Form("hYX%d",i),Form("Beam in Layer %d with DOCA %.0f~%.0f mm",i,FDmin[i],FDmax[i]),NBinXZ,Xmin,Xmax,NBinY,Ymin,Ymax);
    }

    int Iteration = 5;
    TCanvas * canv3D = new TCanvas("canv3D","",512,512);
    canv3D->cd();
    hYXZ->Draw();
    for (int iLayer = 4; iLayer<=5; iLayer++){
        TChain * ichain = new TChain("t","t");
        ichain->Add(Form("root/ana/ana_58.0902xt053l4.i%d.layer%d.root",Iteration,iLayer));
        int    nHitsS;
        double chi2;
        int    theWid0;
        double theFD0;
        double slx;
        double inx;
        double slz;
        double inz;
        ichain->SetBranchAddress("chi2",&chi2);
        ichain->SetBranchAddress("nHitsS",&nHitsS);
        ichain->SetBranchAddress("theWid0",&theWid0);
        ichain->SetBranchAddress("theFD0",&theFD0);
        ichain->SetBranchAddress("slx",&slx);
        ichain->SetBranchAddress("inx",&inx);
        ichain->SetBranchAddress("slz",&slz);
        ichain->SetBranchAddress("inz",&inz);
        canv3D->cd();
        Long64_t nEntries = ichain->GetEntries();
        for (Long64_t iEntry = 0; iEntry<10000; iEntry++){
            ichain->GetEntry(iEntry);
            if (chi2>2||nHitsS<7||theFD0<FDmin[iLayer]||theFD0>FDmax[iLayer]||theWid0!=WireID[iLayer]) continue;
            for (int iBin = 1; iBin<NBinY; iBin++){
                double y = Ymin+(Ymax-Ymin)/NBinY*(iBin-0.5);
                double x = inx+(y-Yref)*slx;
                double z = inz+(y-Yref)*slz;
                hYZ[iLayer]->Fill(z,y);
                hYX[iLayer]->Fill(x,y);
            }
            TPolyLine3D * line = new TPolyLine3D(2);
            line->SetPoint(0,inx,inz,Yref);
            line->SetPoint(1,inx+(Ymin-Yref)*slx,inz+(Ymin-Yref)*slz,Ymin);
            line->SetLineColor(color[iLayer]);
            line->Draw();
        }
        TCanvas * canv = new TCanvas(Form("canv%d",iLayer),"",1024,512);
        canv->Divide(2,1);
        canv->cd(1); hYX[iLayer]->Draw("COLZ");
        canv->cd(2); hYZ[iLayer]->Draw("COLZ");
    }
}
