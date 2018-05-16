{
	std::vector<TString> runNames;
	runNames.push_back("sm20a20first");
	runNames.push_back("s10a20first");
	runNames.push_back("s20a20first");
	runNames.push_back("s30a46first");
	runNames.push_back("s40a46first");
	runNames.push_back("s10a20all");
	runNames.push_back("s20a20all");
	runNames.push_back("s30a46all");
	runNames.push_back("s40a46all");
	TChain * iChains[10];
	TH1F * hists4[10];
	TH1F * hists5[10];
	TCanvas * c4 = new TCanvas();
	TCanvas * c5 = new TCanvas();
	for ( int i = 0; i<runNames.size(); i++){
		iChains[i] = new TChain("t","t");
		iChains[i]->Add(Form("root/ana_1012.0131.%s.i10.layer4.root",runNames[i].Data()));
		iChains[i]->SetLineColor(i+1);
	}

	TFile * ofile = new TFile("output.root","RECREATE");
	c4->cd();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TLegend * l4 = new TLegend(0.7,0.7,0.9,0.9);
	for ( int i = 0; i<runNames.size(); i++){
		if (i==0) iChains[i]->Draw(Form("theDD+res>>h%s4(512,-10,10)",runNames[i].Data()),"chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&theWid==4");
		else      iChains[i]->Draw(Form("theDD+res>>h%s4(512,-10,10)",runNames[i].Data()),"chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&theWid==4","SAME");
		if (ofile->Get(Form("h%s4",runNames[i].Data()))){
			hists4[i] = (TH1F*) ofile->Get(Form("h%s4",runNames[i].Data()));
			l4->AddEntry(hists4[i],runNames[i]);
			printf("Tracks in cell[4,4] run %s: %d\n",runNames[i].Data(),hists4[i]->Integral());
		}
		else printf("Cannot find h%s4\n",runNames[i].Data());
	}
	l4->Draw();
	c5->cd();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	TLegend * l5 = new TLegend(0.7,0.7,0.9,0.9);
	for ( int i = 0; i<runNames.size(); i++){
		if (i==0) iChains[i]->Draw(Form("theDD+res>>h%s5(512,-10,10)",runNames[i].Data()),"chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&theWid==5");
		else      iChains[i]->Draw(Form("theDD+res>>h%s5(512,-10,10)",runNames[i].Data()),"chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&theWid==5","SAME");
		if (ofile->Get(Form("h%s5",runNames[i].Data()))){
			hists5[i] = (TH1F*) ofile->Get(Form("h%s5",runNames[i].Data()));
			l5->AddEntry(hists5[i],runNames[i]);
			printf("Tracks in cell[4,5] run %s: %d\n",runNames[i].Data(),hists5[i]->Integral());
		}
		else printf("Cannot find h%s5\n",runNames[i].Data());
	}
	l5->Draw();
}
