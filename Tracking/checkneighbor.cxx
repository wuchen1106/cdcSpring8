	//===================Get wire position============================
	// For wire position
	TFile * TFile_wirepos = new TFile("wirepos.151225.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_wid;
	int wp_lid;
	int isSenseWire;
	double wp_xhv;
	double wp_yhv;
	double wp_xro;
	double wp_yro;
	TTree_wirepos->SetBranchAddress("isSenseWire",&isSenseWire);
	TTree_wirepos->SetBranchAddress("LayerID",&wp_lid);
	TTree_wirepos->SetBranchAddress("CellID",&wp_wid);
	TTree_wirepos->SetBranchAddress("xd",&wp_xro);
	TTree_wirepos->SetBranchAddress("yd",&wp_yro);
	TTree_wirepos->SetBranchAddress("xu",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yu",&wp_yhv);
	for (int i = 0; i<NLAY; i++){
		wmax[i] = -1;
	}

	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (!isSenseWire) continue;
		if (wp_lid>=1&&wp_lid<=8){
			map_xhv[wp_lid][wp_wid] = wp_xhv;
			map_yhv[wp_lid][wp_wid] = wp_yhv;
			map_xro[wp_lid][wp_wid] = wp_xro;
			map_yro[wp_lid][wp_wid] = wp_yro;
			errord[wp_lid][wp_wid] = 0.2;
			if (wmax[wp_lid]<wp_wid) wmax[wp_lid] = wp_wid;
		}
	}
	for (int k = 1; k<NLAY-1; k++){
		printf("layer %d & %d\n",k,k+1);
		for (int i = 0; i<=wmax[k]; i++){
			for (int j = 0; j<=wmax[k+1]; j++){
				if ((-(map_xro[k+1][j]-map_xhv[k+1][j])+(map_xro[k][i]-map_xhv[k][i]))){
					map_k[k][i][j] = ((map_xro[k+1][j]+map_xhv[k+1][j])-(map_xro[k][i]+map_xhv[k][i]))/(-(map_xro[k+1][j]-map_xhv[k+1][j])+(map_xro[k][i]-map_xhv[k][i]));
				}
				else{
					map_k[k][i][j] = 10;
				}
			}
		}
	}
