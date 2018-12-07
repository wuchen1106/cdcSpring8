{
//	int runNo = 105;
//	TString runName = "katayamapol1.0226.sm10a25n35.i15";
//	//TString runName = "katayama.0226.sm10a25n35.i15";
	int runNo = 1009;
	TString runName = "katayamapol1.0312.sm10a30n35.i3";
	//TString runName = "katayama.0312.sm10a30n35.i3";
//	int runNo = 2047;
//	TString runName = "katayamapol1.0312.sm10a10n0.i10";
//	//TString runName = "katayama.0312.sm10a10n0.i10";

	int testLayer = 4;
	const int nLayers4GG = 9; // number of layers
	const int MAXNTRUNC = 100;
	const double maxTruncRatio = 1/2.; // 50 % at most
	const int nLayersTruncMax = (nLayers4GG-1)*maxTruncRatio; // ignore the first layer
	std::vector<int> nLayersNeededLT;
	nLayersNeededLT.push_back(10);
	nLayersNeededLT.push_back(20);
	nLayersNeededLT.push_back(30);
	nLayersNeededLT.push_back(40);
	nLayersNeededLT.push_back(50);
	nLayersNeededLT.push_back(60);
	nLayersNeededLT.push_back(70);
	nLayersNeededLT.push_back(80);
	nLayersNeededLT.push_back(90);
	nLayersNeededLT.push_back(100);
	nLayersNeededLT.push_back(110);
	nLayersNeededLT.push_back(120);
	nLayersNeededLT.push_back(130);
	nLayersNeededLT.push_back(140);
	nLayersNeededLT.push_back(150);
//	nLayersNeededLT.push_back(1000);
//	nLayersNeededLT.push_back(5000);
//	nLayersNeededLT.push_back(10000);
	std::vector<int> nLayersTruncMaxLT;
	std::vector<int> colors;
	for (int i = 0; i<nLayersNeededLT.size(); i++){
		nLayersTruncMaxLT.push_back(nLayersNeededLT[i]*maxTruncRatio);
		colors.push_back(i+1);
	}

	double cellH = 1.6; // cm
	double W = 39; // eV C4H10
	//double W = 39; // eV CH4
	//double W = 32; // eV C2H6

	TFile * ifile = new TFile(Form("root/res_%d.%s.layer%d.root",runNo,runName.Data(),testLayer));
	TH2D * hggVSX = (TH2D*) ifile->Get("hggVSX");
	hggVSX->SetTitle("");
	hggVSX->SetLineColor(kBlack);
	hggVSX->RebinY(2);

	TCanvas * canvas = new TCanvas();
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	TH1D * hgg = (TH1D*) ifile->Get("hggall");
	//hgg->GetXaxis()->SetRangeUser(0,1.2e5);
	hgg->Draw();
	double averageGG = hgg->GetMean();
	canvas->SaveAs("gasgain.eps");

	hggVSX->RebinX(4);
	hggVSX->RebinY(4);
	hggVSX->GetXaxis()->SetRangeUser(-8.5,8.5);
	hggVSX->GetYaxis()->SetRangeUser(0,1.2e5);
//	hggVSX->Scale(0.1);
	hggVSX->Draw();
	canvas->SaveAs("gasgainVSX.eps");

	TGraphErrors * gr_dedxRes = new TGraphErrors();
	gr_dedxRes->Set(nLayersTruncMax);
	TF1 * fgaus = new TF1("fgaus","gaus",0,3);
	TH1D * h_dedx[nLayersTruncMax];
	printf("dE/dX resolution for prototype-4 (%d layers used, maximum truncated layers %d)\n",nLayers4GG-1,nLayersTruncMax);
	printf("itrunc mean sigma sigma/mean chi2\n");
	double maxReso = 0;
	for (int itrunc = 0; itrunc<nLayersTruncMax; itrunc++){
		h_dedx[itrunc] = (TH1D*) ifile->Get(Form("hdedx%d",itrunc));
		if (!h_dedx[itrunc]) continue;
//		if (itrunc==0) h_dedx[itrunc]->Draw();
//		else h_dedx[itrunc]->Draw("SAME");
//		h_dedx[itrunc]->Fit("fgaus","q","");
		h_dedx[itrunc]->Print();
		h_dedx[itrunc]->Fit("fgaus","qN0","");
		double mean = fgaus->GetParameter(1);
		double sigma = fgaus->GetParameter(2);
		double sigmaErr = fgaus->GetParError(2);
		double chi2 = fgaus->GetChisquare();
		// FIMXE: now use histogram's mean and rms
        sigma = h_dedx[itrunc]->GetRMS();
        sigmaErr = h_dedx[itrunc]->GetRMSError();
//        mean = h_dedx[itrunc]->GetMean();
		gr_dedxRes->SetPoint(itrunc,100.*(nLayers4GG-1-itrunc)/(nLayers4GG-1),sigma/mean*100);
		gr_dedxRes->SetPointError(itrunc,0,sigmaErr/mean*100);
		if (sigma/mean*100>maxReso) maxReso = sigma/mean*100;
	}
	gr_dedxRes->SetName("grdedxLT8");
	gr_dedxRes->GetXaxis()->SetTitle("Fraction of Used Layers [%]");
	gr_dedxRes->GetYaxis()->SetTitle("dE/dX Resolution [%]");
	gr_dedxRes->GetYaxis()->SetRangeUser(0,maxReso*1.5);
	gr_dedxRes->SetMarkerStyle(20);
	gr_dedxRes->Draw("AP");
	gr_dedxRes->Print();
	canvas->SaveAs("dedxreso.eps");
	ifile->Close();

	// Get simulated dE/dX resolution
	TChain * ichain  = new TChain("tree","tree");
	ichain->Add(Form("root/eres_%d.%s.layer%d.root",runNo,runName.Data(),testLayer));
	int    triggerNumber = 0;
	bool   isGood;
	double slx = 0;
	double slz = 0;
	double inx = 0;
	double inz = 0;
	double chargeOnTrack[nLayers4GG];
	double adcsumOnTrack[nLayers4GG];
	int    chargeOnTrackIndex[nLayers4GG];
	ichain->SetBranchAddress("triggerNumber",&triggerNumber);
	ichain->SetBranchAddress("isGood",&isGood);
	for (int lid = 0; lid<nLayers4GG; lid++){ // assuming we kept the same number of truncation options as the number of layers
		ichain->SetBranchAddress(Form("chargeOnTrack%d",lid),&(chargeOnTrack[lid]));
		ichain->SetBranchAddress(Form("adcsumOnTrack%d",lid),&(adcsumOnTrack[lid]));
		ichain->SetBranchAddress(Form("chargeOnTrackIndex%d",lid),&(chargeOnTrackIndex[lid]));
	}
	ichain->SetBranchAddress("slx",&slx);
	ichain->SetBranchAddress("slz",&slz);
	ichain->SetBranchAddress("inx",&inx);
	ichain->SetBranchAddress("inz",&inz);

	// prepare the output ROOT file
	TFile * ofile = new TFile(Form("root/dedx.%d.%s.layer%i.root",runNo,runName.Data(),testLayer),"RECREATE");
	double chargeCount[MAXNTRUNC] = {0};
	double adcsumCount[MAXNTRUNC] = {0};
	double lengthCount[MAXNTRUNC] = {0};
	TH1D * h_dedxLT[MAXNTRUNC] = {0};
	int triggerStart = 0;
	int triggerStop = 0;
	std::vector<double> * chargeList = 0;
	std::vector<double> * adcsumList = 0;
	std::vector<double> * lengthList = 0;
	std::vector<int>    * triggerList = 0;
	std::vector<int>    * indexList = 0;
	TGraphErrors * gr_dedxResVSN = new TGraphErrors();
	gr_dedxResVSN->Set(nLayersNeededLT.size());
	for (int iOpt = 0; iOpt<nLayersNeededLT.size(); iOpt++){
		int nLayersNeeded = nLayersNeededLT[iOpt];
		int nLayersTruncMax = nLayersTruncMaxLT[iOpt];
		if (chargeList ) delete chargeList ; chargeList  = 0;
		if (lengthList ) delete lengthList ; lengthList  = 0;
		if (triggerList) delete triggerList; triggerList = 0;
		if (indexList) delete indexList; indexList = 0;
		TString treeName = Form("t%d",nLayersNeeded);
		TTree * otree = new TTree(treeName,treeName);
		for (int itrunc = 0; itrunc<nLayersTruncMax; itrunc++){
			otree->Branch(Form("charge%d",itrunc),&chargeCount[itrunc]);
			otree->Branch(Form("adcsum%d",itrunc),&adcsumCount[itrunc]);
			otree->Branch(Form("length%d",itrunc),&lengthCount[itrunc]);
		}
		otree->Branch("chargeList",&chargeList);
		otree->Branch("adcsumList",&adcsumList);
		otree->Branch("lengthList",&lengthList);
		otree->Branch("triggerList",&triggerList);
		otree->Branch("indexList",&indexList);
		chargeList = new std::vector<double>(nLayersNeeded);
		adcsumList = new std::vector<double>(nLayersNeeded);
		lengthList = new std::vector<double>(nLayersNeeded);
		triggerList = new std::vector<int>(nLayersNeeded);
		indexList = new std::vector<int>(nLayersNeeded);
		for (int itrunc = 0; itrunc<nLayersTruncMax; itrunc++){
			// don't delete the histogram made in the previous optionLT: still need to keep all of them in the output ROOT file
			// if (h_dedxLT[itrunc]) delete h_dedxLT[itrunc];
			h_dedxLT[itrunc] = new TH1D(Form("hdedxLT%d_%d",nLayersNeeded,itrunc),Form("dE/dX of simulated track with %d hits (%.1f%% truncation)",nLayersNeeded,100.*(nLayersNeeded-itrunc)/nLayersNeeded),512,0,3);
		}

		int nLayersGot = 0;
		int nEntries = ichain->GetEntries();
		for (int iEntry = 0; iEntry<nEntries; iEntry++){
			ichain->GetEntry(iEntry);
			if (!isGood) continue;
			for (int lid = 1; lid<nLayers4GG; lid++){ // ignore layer 0
                (*chargeList)[nLayersGot] = chargeOnTrack[lid];
                (*adcsumList)[nLayersGot] = adcsumOnTrack[lid];
				(*lengthList)[nLayersGot] = cellH*sqrt(1+slx*slx+slz*slz);
				(*triggerList)[nLayersGot] = triggerNumber;
				(*indexList)[nLayersGot] = chargeOnTrackIndex[lid];
				nLayersGot++;
				if (nLayersGot==nLayersNeeded){ // finished one simulated track
					// sort from small to large
					for (int i=0; i<nLayersNeeded; i++){
						for (int j=i+1; j<nLayersNeeded; j++){
							if ((*chargeList)[i]>(*chargeList)[j]){
								double temp = (*chargeList)[i];
								(*chargeList)[i] = (*chargeList)[j];
								(*chargeList)[j] = temp;
								temp = (*adcsumList)[i];
								(*adcsumList)[i] = (*adcsumList)[j];
								(*adcsumList)[j] = temp;
								temp = (*lengthList)[i];
								(*lengthList)[i] = (*lengthList)[j];
								(*lengthList)[j] = temp;
								int itemp = (*triggerList)[i];
								(*triggerList)[i] = (*triggerList)[j];
								(*triggerList)[j] = itemp;
								itemp = (*indexList)[i];
								(*indexList)[i] = (*indexList)[j];
								(*indexList)[j] = itemp;
							}
						}
					}
					// Add them together to form a list of truncated tracks
					double chargeSum = 0;
					double adcsumSum = 0;
					double lengthSum = 0;
					for (int i = 0; i<nLayersNeeded; i++){
						chargeSum+=(*chargeList)[i];
						adcsumSum+=(*adcsumList)[i];
						lengthSum+=(*lengthList)[i];
						if (nLayersNeeded-(i+1)<nLayersTruncMax){
							int itrunc = nLayersNeeded-(i+1); // number of hits left
							chargeCount[itrunc] = chargeSum;
							adcsumCount[itrunc] = adcsumSum;
							lengthCount[itrunc] = lengthSum;
						}
					}
					// fill the histograms
					for (int itrunc = 0; itrunc<nLayersTruncMax; itrunc++){
						double theDE = chargeCount[itrunc]*1e-15/averageGG/1.6e-19/1000*W; // keV
						double theDX = lengthCount[itrunc]; // cm
						double dedx = theDE/theDX;
						h_dedxLT[itrunc]->Fill(dedx);
					}
					nLayersGot = 0;
					otree->Fill();
				}
			}
		}
		TGraphErrors * gr_dedxResLT = new TGraphErrors();
		gr_dedxResLT->Set(nLayersTruncMax);
		printf("dE/dX resolution for simulated track (%d layers used)\n",nLayersNeeded);
		printf("itrunc mean sigma sigma/mean chi2\n");
		double minRelSigma = 0;
		double minRelSigmaErr = 0;
		for (int itrunc = 0; itrunc<nLayersTruncMax; itrunc++){
			h_dedxLT[itrunc]->Write();
			h_dedxLT[itrunc]->Fit("fgaus","qN0","");
			double mean = fgaus->GetParameter(1);
			double sigma = fgaus->GetParameter(2);
			double chi2 = fgaus->GetChisquare();
			double sigmaErr = fgaus->GetParError(2);
            // FIMXE: now use histogram's mean and rms
			sigma = h_dedxLT[itrunc]->GetRMS();
			sigmaErr = h_dedxLT[itrunc]->GetRMSError();
			mean = h_dedxLT[itrunc]->GetMean();
			printf("%d %.3e %.3e %.3e %.3e\n",itrunc,mean,sigma,sigma/mean,chi2);
			gr_dedxResLT->SetPoint(itrunc,100.*(nLayersNeeded-itrunc)/nLayersNeeded,sigma/mean*100);
			gr_dedxResLT->SetPointError(itrunc,0,sigmaErr/mean*100);
			//if (sigma/mean*100<minRelSigma) minRelSigma = sigma/mean*100;
			// FIXME: now we fix it at 90%
			if (!minRelSigma&&1.*(nLayersNeeded-itrunc)/nLayersNeeded<=0.9){
				minRelSigma = sigma/mean*100;
				minRelSigmaErr = sigmaErr/mean*100;
			}
		}
		gr_dedxResLT->SetName(Form("grdedxLT%d",nLayersNeeded));
		gr_dedxResLT->GetXaxis()->SetTitle("Fraction of Used Layers [%]");
		gr_dedxResLT->GetYaxis()->SetTitle("dE/dX Resolution [%]");
		gr_dedxResLT->SetMarkerStyle(20);
		gr_dedxResLT->SetMarkerColor(colors[iOpt]);
		gr_dedxResLT->Write();
		if (iOpt==0)
			gr_dedxResLT->Draw("AP");
		else
			gr_dedxResLT->Draw("PSAME");
		otree->Write();
		gr_dedxResVSN->SetPoint(iOpt,nLayersNeeded,minRelSigma);
		gr_dedxResVSN->SetPointError(iOpt,0,minRelSigmaErr);
	}
	gr_dedxResVSN->SetName("dedxResVSN");
	gr_dedxResVSN->GetXaxis()->SetTitle("Number of used hits (90% accepted)");
	gr_dedxResVSN->GetYaxis()->SetTitle("dE/dX Resolution [%]");
	gr_dedxResVSN->Write();
	gr_dedxRes->Write();
	ofile->Close();
}
