{
	TChain * iChain = new TChain("t","t");
	iChain->Add("xt.a0.m0.z0.iso.1800.root");
	double i_x;
	double i_t;
	iChain->SetBranchAddress("t",&i_t);
	iChain->SetBranchAddress("x",&i_x);

	std::vector<double> v_t_ce;
	std::vector<double> v_x_ce;
	std::vector<double> v_t_m;
	std::vector<double> v_x_m;

	int N = iChain->GetEntries();
	for (int i = 0; i<N; i++){
		iChain->GetEntry(i);
		i_x*=10;
		i_t*=1000;
		if (i_x<0) continue;
		if (i_t>360) continue;
		if (i_t<30){
			v_t_ce.push_back(i_t);
			v_x_ce.push_back(i_x);
		}
		if (i_t>25){
			v_t_m.push_back(i_t);
			v_x_m.push_back(i_x);
		}
	}

	TGraph * gr_ce = new TGraph(v_t_ce.size(),&(v_t_ce[0]),&(v_x_ce[0]));
	TGraph * gr_m = new TGraph(v_t_m.size(),&(v_t_m[0]),&(v_x_m[0]));
	gr_ce->SetMarkerStyle(20);
	gr_m->SetMarkerStyle(20);

	TF1 * fb_ce = new TF1("fb_ce","pol5",0,30);
	TF1 * fb_m = new TF1("fb_m","pol5",25,360);

	gr_ce->Fit("fb_ce","","",0,30);
	gr_m->Fit("fb_m","","",25,360);

	TF1 * fb_delta = new TF1("fb_delta","pol5",0,60);
	for (int i = 0; i<=5; i++){
		double p1 = fb_ce->GetParameter(i);
		double p2 = fb_m->GetParameter(i);
		fb_delta->SetParameter(i,p1-p2);
	}
	TCanvas * c0 = new TCanvas();
	fb_delta->Draw();
	double t_c2m = fb_delta->GetX(0,25,60);
	printf("t_c2m = %.7e\n",t_c2m)

	TCanvas * c1 = new TCanvas();
	gr_ce->Draw("APL");
	fb_m->Draw("SAME");

	TCanvas * c2 = new TCanvas();
	gr_m->Draw("APL");

	TF1 * fl_0 = new TF1("fl_0","(x<=[12]&&x>0)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>[12])*([6]+[7]*x+[8]*x*x+[9]*x*x*x+[10]*x*x*x*x+[11]*x*x*x*x*x)",-25,360);
	TF1 * fr_0 = new TF1("fr_0","(x<=[12]&&x>0)*([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x)+(x>[12])*([6]+[7]*x+[8]*x*x+[9]*x*x*x+[10]*x*x*x*x+[11]*x*x*x*x*x)",-25,360);
	fl_0->SetParameter(12,t_c2m);
	fr_0->SetParameter(12,t_c2m);
	for (int i = 0; i<=5; i++){
		double p1 = fb_ce->GetParameter(i);
		double p2 = fb_m->GetParameter(i);
		fr_0->SetParameter(i,p1);
		fr_0->SetParameter(i+6,p2);
		fl_0->SetParameter(i,-p1);
		fl_0->SetParameter(i+6,-p2);
	}

	TFile * f = new TFile("xt_test.root","RECREATE");
	fl_0->Write();
	fr_0->Write();
	f->Close();
}
