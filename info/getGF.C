{
	TFile * f = new TFile("../../Simulate/proto/wirepos.cdc.pt.root");
	TTree * it = (TTree*)f->Get("t");
	double xu, yu, xd, yd;
	double l,w,t;
	it->SetBranchAddress("t",&t);
	it->SetBranchAddress("w",&w);
	it->SetBranchAddress("l",&l);
	it->SetBranchAddress("xd",&xd);
	it->SetBranchAddress("yd",&yd);
	it->SetBranchAddress("xu",&xu);
	it->SetBranchAddress("yu",&yu);
	for (int i = 0; i<it->GetEntries(); i++){
		it->GetEntry(i);
		if (t==0){
			printf("F 1 {FD} (%.4lf*(1-k)+%.4lf*(1+k))/2*COS+(%.4lf*(1-k)+%.4lf*(1+k))/2*SIN (%.4lf*(1-k)+%.4lf*(1+k))/2*COS-(%.4lf*(1-k)+%.4lf*(1+k))/2*SIN {FV} pi*({FD}/2)*({FD}/2)*l l FDE\n",xd/10.,xu/10.,yd/10.,yu/10.,yd/10.,yu/10.,xd/10.,xu/10.);
		}
		else{
			printf("S 1 {SD} (%.4lf*(1-k)+%.4lf*(1+k))/2*COS+(%.4lf*(1-k)+%.4lf*(1+k))/2*SIN (%.4lf*(1-k)+%.4lf*(1+k))/2*COS-(%.4lf*(1-k)+%.4lf*(1+k))/2*SIN {SV} pi*({SD}/2)*({SD}/2)*l l SDE\n",xd/10.,xu/10.,yd/10.,yu/10.,yd/10.,yu/10.,xd/10.,xu/10.);
		}
	}
}
