{
	TFile * f1 = new TFile("../info/xt.117.t0fit.1xt.7.root");
	TFile * f2 = new TFile("../info/xt.117.t0fit.1xt.8.root");
	TF1* f1_left[9];
	TF1* f1_left_end[9];
	TF1* f1_right[9];
	TF1* f1_right_end[9];
	TF1* f2_left[9];
	TF1* f2_left_end[9];
	TF1* f2_right[9];
	TF1* f2_right_end[9];
	TF1* fd_left[9];
	TF1* fd_left_end[9];
	TF1* fd_right[9];
	TF1* fd_right_end[9];
	double xtemp[2]; xtemp[0] = 0; xtemp[1] = 460;
	double ytemp[2]; ytemp[0] = -0.5; ytemp[1] = 0.5;
	TGraph * gtemp = new TGraph(2,xtemp,ytemp);
	gtemp->Draw("AP");
	std::vector<int >color;
	color.push_back(kBlack);
	color.push_back(kRed);
	color.push_back(kOrange);
	color.push_back(kYellow);
	color.push_back(kGreen);
	color.push_back(kCyan);
	color.push_back(kBlue);
	color.push_back(kMagenta);
	color.push_back(kBlack);

	std::vector<double> * vdx_left[9];
	std::vector<double> * vdx_right[9];
	std::vector<double> * vt_left[9];
	std::vector<double> * vt_right[9];
	TGraph * gd_left[9];
	TGraph * gd_right[9];

	TLegend * legend = new TLegend(0.7,0.5,0.9,0.9);
	for (int il = 1; il<=8; il++){
		f2_left[il] = (TF1*) f2->Get(Form("f_left_%d_0",il));
		f2_left_end[il] = (TF1*) f2->Get(Form("f_left_end_%d_0",il));
		f2_right[il] = (TF1*) f2->Get(Form("f_right_%d_0",il));
		f2_right_end[il] = (TF1*) f2->Get(Form("f_right_end_%d_0",il));
		f1_left[il] = (TF1*) f1->Get(Form("f_left_%d_0",il));
		f1_left_end[il] = (TF1*) f1->Get(Form("f_left_end_%d_0",il));
		f1_right[il] = (TF1*) f1->Get(Form("f_right_%d_0",il));
		f1_right_end[il] = (TF1*) f1->Get(Form("f_right_end_%d_0",il));

		vdx_right[il] = new std::vector<double>;
		vdx_left[il] = new std::vector<double>;
		vt_right[il] = new std::vector<double>;
		vt_left[il] = new std::vector<double>;

		double av = 0;
		for (double it = f2_left[il]->GetXmin(); it<=f2_left[il]->GetXmax(); it+=2){
			vt_left[il]->push_back(it);
			vdx_left[il]->push_back(f2_left[il]->Eval(it)-f1_left[il]->Eval(it));
			av+=fabs(f2_left[il]->Eval(it)-f1_left[il]->Eval(it));
		}
		//for (double it = f2_left_end[il]->GetXmin(); it<=f2_left_end[il]->GetXmax(); it+=2){
		//	vt_left[il]->push_back(it);
		//	vdx_left[il]->push_back(f2_left_end[il]->Eval(it)-f1_left_end[il]->Eval(it));
		//	av+=fabs(f2_left_end[il]->Eval(it)-f1_left_end[il]->Eval(it));
		//}
		for (double it = f2_right[il]->GetXmin(); it<=f2_right[il]->GetXmax(); it+=2){
			vt_right[il]->push_back(it);
			vdx_right[il]->push_back(f2_right[il]->Eval(it)-f1_right[il]->Eval(it));
			av+=fabs(f2_right[il]->Eval(it)-f1_right[il]->Eval(it));
		}
		//for (double it = f2_right_end[il]->GetXmin(); it<=f2_right_end[il]->GetXmax(); it+=2){
		//	vt_right[il]->push_back(it);
		//	vdx_right[il]->push_back(f2_right_end[il]->Eval(it)-f1_right_end[il]->Eval(it));
		//	av+=fabs(f2_right_end[il]->Eval(it)-f1_right_end[il]->Eval(it));
		//}
		av/=(vt_right[il]->size()+vt_left[il]->size());
		gd_left[il]=new TGraph(vt_left[il]->size(),&((*vt_left[il])[0]),&((*vdx_left[il])[0]));
		gd_right[il]=new TGraph(vt_right[il]->size(),&((*vt_right[il])[0]),&((*vdx_right[il])[0]));
		gd_left[il]->SetLineColor(color[il]);
		gd_right[il]->SetLineColor(color[il]);
		gd_left[il]->Draw("LSAME");
		gd_right[il]->Draw("LSAME");
		legend->AddEntry(gd_left[il],Form("Layer %d: %.2lf mm",il,av),"lp");
	}
	legend->Draw("SAME");
	new TCanvas();
	ytemp[0] = -10; ytemp[1] = 10;
	gtemp = new TGraph(2,xtemp,ytemp);
	gtemp->Draw("AP");
	for (int il = 1; il<=8; il++){
		f2_left[il]->SetLineColor(color[il]);
		f2_left_end[il]->SetLineColor(color[il]);
		f2_right[il]->SetLineColor(color[il]);
		f2_right_end[il]->SetLineColor(color[il]);
		f2_left[il]->Draw("LSAME");
		f2_left_end[il]->Draw("LSAME");
		f2_right[il]->Draw("LSAME");
		f2_right_end[il]->Draw("LSAME");
	}
}
