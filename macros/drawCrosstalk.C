{
	TChain * ichain = new TChain("t","t");
	ichain->Add("root/ana_1009.0312.sm10a30n35.i3.layer4.root");
	TString cutG = "chi2<2&&nHitsS>=7&&abs(slz)<0.1";
	gStyle->SetOptStat(0);

	TCanvas * canv_peakVSdt = new TCanvas();
	gPad->SetLogz(1);
	ichain->Draw("peak-ped:driftT>>hpeakVSdt(1200,-250,1000,512,-100,500)",cutG,"COLZ");
	hpeakVSdt->SetTitle("ADC peak VS drift time (all peaks)");
	hpeakVSdt->GetXaxis()->SetTitle("drift time [ns]");
	hpeakVSdt->GetYaxis()->SetTitle("ADC sum");
	canv_peakVSdt->SaveAs("peakVSdt.png");
	canv_peakVSdt->SaveAs("peakVSdt.pdf");

	TCanvas * canv_peakVSdtN = new TCanvas();
	gPad->SetLogz(1);
	ichain->Draw("peak-ped:driftT>>hpeakVSdtN(1200,-250,1000,800,-100,500)",cutG+"&&abs(fitD)>10","COLZ");
	hpeakVSdtN->SetTitle("ADC peak VS drift time (noise hit)");
	hpeakVSdtN->GetXaxis()->SetTitle("drift time [ns]");
	hpeakVSdtN->GetYaxis()->SetTitle("ADC sum");
	canv_peakVSdtN->SaveAs("peakVSdtN.png");
	canv_peakVSdtN->SaveAs("peakVSdtN.pdf");

	TCanvas * canv_geo = new TCanvas();
	ichain->Draw("layerID-highLid:wireID-highWid>>hgeo(21,-10.5,10.5,17,-8.5,8.5)",cutG+"&&(wireID!=highWid||layerID!=highLid)&&abs(fitD)>10","COLZ");
	hgeo->SetTitle("#Delta_{layerID} VS #Delta_{cellID}");
	hgeo->GetXaxis()->SetTitle("#Delta_{cellID}");
	hgeo->GetYaxis()->SetTitle("#Delta_{layerID}");
	canv_geo->SaveAs("geo.png");
	canv_geo->SaveAs("geo.pdf");

	TCanvas * canv_ch = new TCanvas();
	ichain->Draw("boardID*6+(channelID-(channelID%8))/8:highBid*6+(highCh-(highCh%8))/8>>hch(12,-0.5,11.5,12,-0.5,11.5)",cutG+"&&(wireID!=highWid||layerID!=highLid)&&abs(fitD)>10","COLZ");
	hch->SetTitle("ASD chip ID of noise hit and the source hit");
	hch->GetXaxis()->SetTitle("ASD_{ID} of source hit");
	hch->GetYaxis()->SetTitle("ASD_{ID} of noise hit");
	canv_ch->SaveAs("ch.png");
	canv_ch->SaveAs("ch.pdf");

	TCanvas * canv_chr = new TCanvas();
	ichain->Draw("boardID*6+(channelID-(channelID%8))/8:highBid*6+(highCh-(highCh%8))/8>>hchr(12,-0.5,11.5,12,-0.5,11.5)",cutG+"&&(wireID!=highWid||layerID!=highLid)&&abs(fitD)>10&&(abs(layerID-highLid)>3||abs(wireID-highWid)>5)","COLZ");
	hch->SetTitle("ASD chip ID of noise hit and the source hit (remote in chamber)");
	hch->GetXaxis()->SetTitle("ASD_{ID} of source hit");
	hch->GetYaxis()->SetTitle("ASD_{ID} of noise hit");
	canv_chr->SaveAs("chr.png");
	canv_chr->SaveAs("chr.pdf");

	TCanvas * canv_prVSddtNA = new TCanvas();
	gPad->SetLogz(1);
	ichain->Draw("(peak-ped)/highSum*100:driftT-highDT>>hprVSddtNA(1200,-625,625,256,-5,5)",cutG+"&&(wireID!=highWid||layerID!=highLid)&&abs(fitD)>10&&boardID*6+(channelID-(channelID%8))/8!=highBid*6+(highCh-(highCh%8))/8","COLZ");
	hprVSddtNA->SetTitle("ADC ratio VS drift time difference (different ASD chips)");
	hprVSddtNA->GetXaxis()->SetTitle("#Delta_{driftT} [ns]");
	hprVSddtNA->GetYaxis()->SetTitle("ADC ratio [%]");
	canv_prVSddtNA->SaveAs("prVSddtNA.png");
	canv_prVSddtNA->SaveAs("prVSddtNA.pdf");

	TCanvas * canv_prVSddtNH = new TCanvas();
	gPad->SetLogz(1);
	ichain->Draw("(peak-ped)/highSum*100:driftT-highDT>>hprVSddtNH(1200,-625,625,256,-5,5)",cutG+"&&(wireID!=highWid||layerID!=highLid)&&abs(fitD)>10&&boardID*6+(channelID-(channelID%8))/8==highBid*6+(highCh-(highCh%8))/8&&(abs(layerID-highLid)>1||abs(wireID-highWid)>3)","COLZ");
	hprVSddtNH->SetTitle("ADC ratio VS drift time difference (same ASD chip, remote in chamber)");
	hprVSddtNH->GetXaxis()->SetTitle("#Delta_{driftT} [ns]");
	hprVSddtNH->GetYaxis()->SetTitle("ADC ratio [%]");
	canv_prVSddtNH->SaveAs("prVSddtNH.png");
	canv_prVSddtNH->SaveAs("prVSddtNH.pdf");

	TCanvas * canv_prVSddtNAW1L0 = new TCanvas();
	gPad->SetLogz(1);
	ichain->Draw("(peak-ped)/highSum*100:driftT-highDT>>hprVSddtNAW1L0(1200,-625,625,256,-5,5)",cutG+"&&abs(wireID-highWid)==1&&abs(layerID-highLid)==0&&abs(fitD)>10&&boardID*6+(channelID-(channelID%8))/8!=highBid*6+(highCh-(highCh%8))/8","COLZ");
	hprVSddtNAW1L0->SetTitle("ADC ratio VS drift time difference (different ASD chip, same layer, neighbour cell)");
	hprVSddtNAW1L0->GetXaxis()->SetTitle("#Delta_{driftT} [ns]");
	hprVSddtNAW1L0->GetYaxis()->SetTitle("ADC ratio [%]");
	canv_prVSddtNAW1L0->SaveAs("prVSddtNAW1L0.png");
	canv_prVSddtNAW1L0->SaveAs("prVSddtNAW1L0.pdf");

	TCanvas * canv_prVSddtNAW2L0 = new TCanvas();
	gPad->SetLogz(1);
	ichain->Draw("(peak-ped)/highSum*100:driftT-highDT>>hprVSddtNAW2L0(1200,-625,625,256,-5,5)",cutG+"&&abs(wireID-highWid)==2&&abs(layerID-highLid)==0&&abs(fitD)>10&&boardID*6+(channelID-(channelID%8))/8!=highBid*6+(highCh-(highCh%8))/8","COLZ");
	hprVSddtNAW2L0->SetTitle("ADC ratio VS drift time difference (different ASD chip, same layer, next neighbour cell)");
	hprVSddtNAW2L0->GetXaxis()->SetTitle("#Delta_{driftT} [ns]");
	hprVSddtNAW2L0->GetYaxis()->SetTitle("ADC ratio [%]");
	canv_prVSddtNAW2L0->SaveAs("prVSddtNAW2L0.png");
	canv_prVSddtNAW2L0->SaveAs("prVSddtNAW2L0.pdf");
}
