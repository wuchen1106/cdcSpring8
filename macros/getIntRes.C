{
	TFile * ifile = new TFile("root/res_1012.0226.sm10a30n35.i15.layer4.root");
	float source[256];
	float response[256];
	TF1 * fres = new TF1("fres","gaus",-2,2);
	fres->SetParameter(0,1);
	fres->SetParameter(1,0);
	fres->SetParameter(2,0.05);
	for (int ibin = 0; ibin<256; ibin++){
		double x = 4/256.*(ibin+0.5)-2;
		response[ibin] = fres->Eval(x);
	}
//	for (int islice = 0; islice <20; islice++){
	for (int islice = 0; islice <1; islice++){
		TH1D * hres = (TH1D*) ifile->Get(Form("hresX%d",islice));
		TH1D * d = new TH1D("d","",256,-2,2);
		hres->Draw();
		for (int ibin = 0; ibin<256; ibin++){
			source[ibin] = hres->GetBinContent(ibin+1);
		}
		TSpectrum * s = new TSpectrum();
		s->Deconvolution(source,response,256,100,1,1);
		for (int ibin = 0; ibin<256; ibin++){
			d->SetBinContent(ibin+1,source[ibin]);
		}
		d->SetLineColor(kRed);
		d->Draw("SAME");
	}
}
