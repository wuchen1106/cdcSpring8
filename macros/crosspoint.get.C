{
    double chamberHL = 599.17/2; // mm
    int NLAY = 9;
    int NCEL = 11;
    double U = 8;
    double angle = 18.4*3.14/180;
    double offset = tan(angle)*U*2;

    double map_xc[NLAY][NCEL];
    double map_yc[NLAY][NCEL];
    double map_xhv[NLAY][NCEL];
    double map_yhv[NLAY][NCEL];
    double map_xro[NLAY][NCEL];
    double map_yro[NLAY][NCEL];
    double map_k[NLAY][NCEL][NCEL];
    int map_check[NLAY][NCEL];
	TFile * TFile_wirepos = new TFile("Input/wire-position.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bid;
	int wp_wid;
	int wp_lid;
	int wp_ch;
	double wp_xhv;
	double wp_yhv;
	double wp_xro;
	double wp_yro;
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("b",&wp_bid);
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xro",&wp_xro);
	TTree_wirepos->SetBranchAddress("yro",&wp_yro);
	TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
	for(int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
			map_check[ilayer][iwire]=0;
		}
	}
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<NLAY){
			map_xhv[wp_lid][wp_wid] = wp_xhv;
			map_yhv[wp_lid][wp_wid] = wp_yhv;
			map_xro[wp_lid][wp_wid] = wp_xro;
			map_yro[wp_lid][wp_wid] = wp_yro;
			map_xc[wp_lid][wp_wid] = (wp_xhv+wp_xro)/2.;
			map_yc[wp_lid][wp_wid] = (wp_yhv+wp_yro)/2.;
			map_check[wp_lid][wp_wid] = 1;
		}
	}
	for(int ilayer = 0; ilayer<NLAY-1; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
			if (!map_check[ilayer][iwire]) continue;
			for (int jwire = 0; jwire<NCEL; jwire++){
				if (!map_check[ilayer+1][jwire]) continue;
				double xro_u = map_xro[ilayer+1][jwire];
				double xro_d = map_xro[ilayer][iwire];
				double xhv_u = map_xhv[ilayer+1][jwire];
				double xhv_d = map_xhv[ilayer][iwire];
				if ((-(xro_u-xhv_u)+(xro_d-xhv_d))){
					map_k[ilayer][iwire][jwire] = ((xro_u+xhv_u)-(xro_d+xhv_d)+offset)/((xro_d-xhv_d)-(xro_u-xhv_u));
				}
				else{
					map_k[ilayer][iwire][jwire] = 10;
				}
				printf("%d %d %d %d %lf %lf\n",ilayer,ilayer+1,iwire,jwire,map_k[ilayer][iwire][jwire]*chamberHL,(xro_d*(1+map_k[ilayer][iwire][jwire])+xhv_d*(1-map_k[ilayer][iwire][jwire]))/2.);
			}
		}
	}
}
