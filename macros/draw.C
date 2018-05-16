{
//	TFile * ifile = new TFile("root/ana_117.0126.i8.layer4.root");
	TFile * ifile = new TFile("root/ana_1012.0131.s30a46allt0_0_m1.i10.layer4.root");
	int lid = 4;
	int sumCut = -30;
	TTree * t = (TTree*) ifile->Get("t");
	t->Draw("fitD-driftD>>hsig(512,-10,10)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip==theIp",lid));
	t->Draw("fitD-driftD>>hearly(512,-10,10)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip<theIp",lid),"SAME");
	t->Draw("fitD-driftD>>hlate(512,-10,10)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip>theIp",lid),"SAME");
	new TCanvas();
	t->Draw("sum:driftT-theDT>>hnoiseT(1200,-625,625,700,-100,600)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip!=theIp",lid),"COLZ");
	new TCanvas();
	t->Draw("sum:driftT>>hearlyDT(600,-25,600,700,-100,600)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip<theIp",lid),"COLZ");
	new TCanvas();
	t->Draw("sum:driftT>>hlateDT(600,-25,600,700,-100,600)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip>theIp",lid),"COLZ");
	new TCanvas();
	t->Draw("theSum:theDT>>hsigDT(600,-25,600,700,-100,600)","chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1","COLZ");
	new TCanvas();
	t->Draw("driftT>>hdTearlyL(600,-25,600)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip<theIp&&sum>%d",lid,sumCut));
	t->Draw("driftT>>hdTearly(600,-25,600)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip<theIp",lid),"SAME");
	new TCanvas();
	t->Draw("fitD:driftT>>hXTearlyL(600,-25,600,512,-10,10)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip<theIp&&sum>%d",lid,sumCut),"COLZ");
	new TCanvas();
	t->Draw("fitD:driftT>>hXTearly(600,-25,600,512,-10,10)",Form("chi2<1&&nHitsS>=6&&abs(slz)<0.15&&has==1&&layerID==%d&&wireID==theWid&&ip<theIp",lid),"COLZ");
}
