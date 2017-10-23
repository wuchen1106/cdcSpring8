{
    TFile * ifile = new TFile("../root/t_117.layer4.t0fit.1xt.18.root");
    TTree * itree = (TTree*) ifile->Get("t");
    TString cut1 = "layerID==4&&chi2<5&&abs(inZ)<19&&abs(slZ-0.01)<0.1&&abs(slX)<0.09";
    TString cut2 = "abs(fitD-driftD)<0.5";
    int N[18];
    int n[18];
    int x0 = -9;
    int dx = 1;
    for ( int i = 0; i<18; i++){
        int x = x0+i*dx;
        TString cut3 = Form("fitD>=%d&&fitD<%d",x,x+dx);
        N[i] = itree->GetEntries(cut1+"&&"+cut3);
        n[i] = itree->GetEntries(cut1+"&&"+cut2+"&&"+cut3);
        std::cout<<i<<" "<<x+dx/2.<<" "<<N[i]<<" "<<n[i]<<" "<<((double)n[i])/N[i]<<std::endl;
    }
}
