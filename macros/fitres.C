{
	TChain * ichain = new TChain("t","t");
	ichain->Add("root/ana_1009.0312.sm10a30n35.i3.layer4.root");
	ichain->Draw("driftD-fitD>>h(256,-2,2)","layerID==4&&wireID==theWid&&ip==theIp&&has==1&&chi2<2&&nHitsS>=7&&driftT<15");
	TF1 * f = new TF1("f1","gaus(0)+gaus(3)",-2,2);
	double heightMax = h->GetMaximum();
	f->SetParameter(0,heightMax);
	f->SetParameter(1,0);
	f->SetParameter(2,0.15);
	f->SetParameter(3,heightMax);
	f->SetParameter(4,0);
	f->SetParameter(5,0.5);
	f->SetParLimits(0,0,heightMax*1.5);
	f->SetParLimits(1,-1,1);
	f->SetParLimits(2,0,0.5);
	f->SetParLimits(3,0,heightMax*1.5);
	f->SetParLimits(4,-1,1);
	f->SetParLimits(5,0,2);
	h->Fit("f1");
	double sigma1 = f->GetParameter(2);
	double sigma2 = f->GetParameter(5);
	double N1 = f->GetParameter(0);
	double N2 = f->GetParameter(3);
	double sigma = sqrt(N1/(N1+N2)*sigma1*sigma1+N2/(N1+N2)*sigma2*sigma2);
	printf("sigma = %.3e\n",sigma);
}
