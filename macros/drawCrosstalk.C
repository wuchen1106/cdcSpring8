{
	TChain * ichain = new TChain("t","t");
	ichain->Add("root/ana_1012.0219.sp10a30n35.i10.layer4.root");
	TCanvas * canv_All = new TCanvas();
	ichain->Draw("(driftT-highDT)*0.96>>hdifDTAll(2000,-1000,1000)","(wireID!=highWid||layerID!=highLid||ip!=highIp)&&aa<10");
	TCanvas * canv_ASD = new TCanvas();
	ichain->Draw("(driftT-highDT)*0.96>>hdifDTASD(2000,-1000,1000)","(wireID!=highWid||layerID!=highLid||ip!=highIp)&&boardID==highBid&&fmod(channelID/8,1)==fmod(highCh/8,1)&&aa<10");
	TCanvas * canv_NeighbourSL = new TCanvas();
	ichain->Draw("(driftT-highDT)*0.96>>hdifDTNeighbourSL(2000,-1000,1000)","(wireID!=highWid||layerID!=highLid||ip!=highIp)&&layerID==highLid&&abs(wireID-highWid)==1&&aa<10");
	TCanvas * canv_NeighbourNL = new TCanvas();
	ichain->Draw("(driftT-highDT)*0.96>>hdifDTNeighbourNL(2000,-1000,1000)","(wireID!=highWid||layerID!=highLid||ip!=highIp)&&abs(layerID-highLid)==1&&abs(wireID-highWid)<=1&&aa<10");
	TCanvas * canv_NotNeighbour = new TCanvas();
	ichain->Draw("(driftT-highDT)*0.96>>hdifDTNotNeighbour(2000,-1000,1000)","(wireID!=highWid||layerID!=highLid||ip!=highIp)&&((layerID==highLid&&abs(wireID-highWid)>1)||(abs(layerID-highLid)==1&&abs(wireID-highWid)>1)||abs(layerID-highLid)>1)&&aa<10");
	TCanvas * canv_nHitsSmall = new TCanvas();
	ichain->Draw("nHitsSmallSASD:nHitsSmallAll","","COLZ");
}
