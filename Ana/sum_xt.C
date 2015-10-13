{
	int Off[8];
	Off[1]=0.0;
	Off[2]=0.0;
	Off[3]=0.0;
	Off[4]=0.0;
	Off[5]=0.0;
	Off[6]=0.0;
	Off[7]=0.0;
	TChain *chain0 = new TChain("t","t");
	chain0->Add("../info/xt.117.root");
	TChain *chain1 = new TChain("t","t");
	chain1->Add("../info/xt.117.8.new7.root");
	chain0->SetLineColor(kRed);
	chain0->SetMarkerColor(kRed);
	chain0->Draw("t:x","","LP");
	chain1->SetMarkerStyle(7);
	for (int ilayer = 1; ilayer<=7; ilayer++){
		chain1->Draw(Form("t:x+%lf",Off[ilayer]),Form("lid==%d&&wid==%d",ilayer,ilayer>=4?5:4),"LPSAME");
	}
}
