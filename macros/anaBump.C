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
    TH2D * hXZ[NLAY];
    TH2D * hLD[NLAY];
    for (int i = 0; i<NLAY; i++){
        hYZ[i] = new TH2D(Form("hYZ%d",i),Form("Beam in Layer %d with DOCA %.0f~%.0f mm",i,FDmin[i],FDmax[i]),NBinXZ,Zmin,Zmax,NBinY,Ymin,Ymax);
        hYX[i] = new TH2D(Form("hYX%d",i),Form("Beam in Layer %d with DOCA %.0f~%.0f mm",i,FDmin[i],FDmax[i]),NBinXZ,Xmin,Xmax,NBinY,Ymin,Ymax);
        hXZ[i] = new TH2D(Form("hXZ%d",i),Form("Beam in Layer %d with DOCA %.0f~%.0f mm",i,FDmin[i],FDmax[i]),NBinXZ,Xmin,Xmax,NBinXZ,Zmin,Zmax);
        hLD[i] = new TH2D(Form("hLD%d",i),Form("Beam in Layer %d with DOCA %.0f~%.0f mm",i,FDmin[i],FDmax[i]),NBinXZ,-10,10,NLAY,0,NLAY);
    }

    int Iteration = 5;
    TChain * ichain = new TChain("t","t");
    ichain->Add("root/ana/ana_58.0902xt053l4.i5.layer4.root");
    int    nHitsS;
    double chi2;
    double slx;
    double inx;
    double slz;
    double inz;
    std::vector<double> * fitD = 0;
    std::vector<int> * layerID = 0;
    std::vector<int> * wireID = 0;
    ichain->SetBranchAddress("chi2",&chi2);
    ichain->SetBranchAddress("nHitsS",&nHitsS);
    ichain->SetBranchAddress("fitD",&fitD);
    ichain->SetBranchAddress("layerID",&layerID);
    ichain->SetBranchAddress("wireID",&wireID);
    ichain->SetBranchAddress("slx",&slx);
    ichain->SetBranchAddress("inx",&inx);
    ichain->SetBranchAddress("slz",&slz);
    ichain->SetBranchAddress("inz",&inz);
    double minFD[NLAY];
    Long64_t nEntries = ichain->GetEntries();
    for (Long64_t iEntry = 0; iEntry<nEntries; iEntry++){
        ichain->GetEntry(iEntry);
        for (int iLayer = 1; iLayer<NLAY; iLayer++){
            for (int jLayer = 1; jLayer<NLAY; jLayer++){
                minFD[jLayer] = 10;
            }
            bool found = false;
            for (int iHit = 0; iHit<fitD->size(); iHit++){
                //if (chi2>2||nHitsS<7||fitD->at(iHit)<FDmin[iLayer]||fitD->at(iHit)>FDmax[iLayer]||layerID->at(iHit)!=iLayer||wireID->at(iHit)!=WireID[iLayer]) continue;
                if (chi2>2||nHitsS<7) continue;
                if (fabs(minFD[layerID->at(iHit)])>fabs(fitD->at(iHit))) minFD[layerID->at(iHit)] = fitD->at(iHit);
                if(fitD->at(iHit)<FDmin[iLayer]||fitD->at(iHit)>FDmax[iLayer]||layerID->at(iHit)!=iLayer) continue;
                int preBinXZ = -1;
                for (int iBin = 1; iBin<NBinY; iBin++){
                    double y = Ymin+(Ymax-Ymin)/NBinY*(iBin-0.5);
                    double x = inx+(y-Yref)*slx;
                    double z = inz+(y-Yref)*slz;
                    int binXZ = hXZ[iLayer]->FindBin(x,z);
                    hYZ[iLayer]->Fill(z,y);
                    hYX[iLayer]->Fill(x,y);
                    if (binXZ!=preBinXZ) hXZ[iLayer]->Fill(x,z);
                    preBinXZ = binXZ;
                    found = true;
                }
            }
            if (found){
                for (int jLayer = 1; jLayer<NLAY; jLayer++){
                    hLD[iLayer]->Fill(minFD[jLayer],jLayer);
                }
            }
        }
    }
    TCanvas * canvYX = new TCanvas("canvYX","",1024,512);
    TCanvas * canvYZ = new TCanvas("canvYZ","",1024,512);
    TCanvas * canvXZ = new TCanvas("canvXZ","",1024,512);
    TCanvas * canvLD = new TCanvas("canvLD","",1024,512);
    canvYX->Divide(4,2);
    canvYZ->Divide(4,2);
    canvXZ->Divide(4,2);
    canvLD->Divide(4,2);
    for (int iLayer = 1; iLayer<NLAY; iLayer++){
        int i = ((iLayer+1)%2)*4+(iLayer+1)/2;
        canvYX->cd(i); gPad->SetGridx(1); gPad->SetGridy(1); hYX[iLayer]->Draw("COLZ");
        canvYZ->cd(i); gPad->SetGridx(1); gPad->SetGridy(1); hYZ[iLayer]->Draw("COLZ");
        canvXZ->cd(i); gPad->SetGridx(1); gPad->SetGridy(1); hXZ[iLayer]->Draw("COLZ");
        canvLD->cd(i); gPad->SetGridx(1); gPad->SetGridy(1); hLD[iLayer]->Draw("COLZ");
    }
}
