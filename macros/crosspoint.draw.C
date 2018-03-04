{
    TFile * ifile = new TFile("../../input/wirepos.proto4.sw.ver7.root");
    TTree * itree = (TTree*) ifile->Get("t");
    int l,w;
    double xro,yro,xhv,yhv;
    itree->SetBranchAddress("l",&l);
    itree->SetBranchAddress("w",&w);
    itree->SetBranchAddress("xro",&xro);
    itree->SetBranchAddress("yro",&yro);
    itree->SetBranchAddress("xhv",&xhv);
    itree->SetBranchAddress("yhv",&yhv);

    TFile * ifile2 = new TFile("crosspoint.root");
    TTree * itree2 = (TTree*) ifile2->Get("t");
    int l1,l2,w1,w2;
    double zcross,xcross;
    itree2->SetBranchAddress("l1",&l1);
    itree2->SetBranchAddress("l2",&l2);
    itree2->SetBranchAddress("w1",&w1);
    itree2->SetBranchAddress("w2",&w2);
    itree2->SetBranchAddress("z",&zcross);
    itree2->SetBranchAddress("x",&xcross);

    int check[9][11];
    double x[9][11][2];
    double z[9][11][2];
    double z0 = 599.17/2;
    double xc[8][11][11];
    double zc[8][11][11];
	for(int ilayer = 0; ilayer<9; ilayer++){
		for (int iwire = 0; iwire<11; iwire++){
		    x[ilayer][iwire][0] = 0;
		    z[ilayer][iwire][0] = 0;
		    x[ilayer][iwire][1] = 0;
		    z[ilayer][iwire][1] = 0;
			check[ilayer][iwire]=0;
			if (ilayer <8){
                for (int jwire = 0; jwire<11; jwire++){
                    xc[ilayer][iwire][jwire] = 999;
                    zc[ilayer][iwire][jwire] = 999;
                }
            }
		}
	}

    int nEntries = itree->GetEntries();
    for (int iEntry = 0; iEntry<nEntries; iEntry++){
        itree->GetEntry(iEntry);
        check[l][w] = 1;
        x[l][w][0] = xhv;
        z[l][w][0] = -z0;
        x[l][w][1] = xro;
        z[l][w][1] = z0;
    }

    nEntries = itree2->GetEntries();
    for (int iEntry = 0; iEntry<nEntries; iEntry++){
        itree2->GetEntry(iEntry);
        xc[l1][w1][w2] = xcross;
        zc[l1][w1][w2] = zcross;
    }

    int color[11];
    int icolor = 0;
    color[icolor++] = kBlack;
    color[icolor++] = kMagenta+2;
    color[icolor++] = kMagenta;
    color[icolor++] = kBlue;
    color[icolor++] = kBlue-2;
    color[icolor++] = kCyan;
    color[icolor++] = kGreen;
    color[icolor++] = kYellow;
    color[icolor++] = kOrange;
    color[icolor++] = kRed;
    color[icolor++] = kRed+2;
    TGraph * g_wire[9][11];
    for (int ilayer = 0; ilayer<9; ilayer++){
		for (int iwire = 0; iwire<11; iwire++){
		    if (check[ilayer][iwire]==0) continue;
			g_wire[ilayer][iwire] = new TGraph(2,z[ilayer][iwire],x[ilayer][iwire]);
			g_wire[ilayer][iwire]->SetLineColor(color[iwire]);
			g_wire[ilayer][iwire]->SetMarkerColor(color[iwire]);
        }
    }

    for (int ilayer = 0; ilayer<=7; ilayer++){
        TCanvas * c = new TCanvas(Form("c%d",ilayer),"c",1024,768);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
		gPad->SetGridx(1);
		gPad->SetGridy(1);
        TGraph * g_all = new TGraph(9*11*2,&(z[0][0][0]),&(x[0][0][0]));
        g_all->SetTitle(Form("layer #%d and layer #%d",ilayer,ilayer+1));
        g_all->SetMarkerColor(kWhite);
        g_all->GetXaxis()->SetTitle("z [mm]");
        g_all->GetYaxis()->SetTitle("x [mm]");
        g_all->Draw("AP");
		for (int iwire = 0; iwire<11; iwire++){
            TLatex * text1 = new TLatex(z[ilayer][0][0]-10,x[ilayer][0][0]+10,Form("# %d",ilayer));
            text1->SetTextSize(0.02);
            text1->Draw("SAME");
            TLatex * text2 = new TLatex(z[ilayer+1][0][1]+10,x[ilayer+1][0][1]+10,Form("# %d",ilayer+1));
            text2->SetTextSize(0.02);
            text2->Draw("SAME");
		    if (check[ilayer][iwire]==1){
		        g_wire[ilayer][iwire]->Draw("PLSAME");
		        TLatex * textl = new TLatex(z[ilayer][iwire][0]-10,x[ilayer][iwire][0],Form("%d",iwire));
		        TLatex * textr = new TLatex(z[ilayer][iwire][1]-10,x[ilayer][iwire][1],Form("%d",iwire));
		        textl->SetTextColor(color[iwire]);
		        textr->SetTextColor(color[iwire]);
		        textl->SetTextSize(0.02);
		        textr->SetTextSize(0.02);
		        textl->Draw("SAME");
		        textr->Draw("SAME");
            }
		    if (check[ilayer+1][iwire]==1){
		        g_wire[ilayer+1][iwire]->Draw("PLSAME");
		        TLatex * textl = new TLatex(z[ilayer+1][iwire][0]+10,x[ilayer+1][iwire][0],Form("%d",iwire));
		        TLatex * textr = new TLatex(z[ilayer+1][iwire][1]+10,x[ilayer+1][iwire][1],Form("%d",iwire));
		        textl->SetTextColor(color[iwire]);
		        textr->SetTextColor(color[iwire]);
		        textl->SetTextSize(0.02);
		        textr->SetTextSize(0.02);
		        textl->Draw("SAME");
		        textr->Draw("SAME");
            }
            for (int jwire = 0; jwire < 11; jwire++){
                if (fabs(zc[ilayer][iwire][jwire])>300) continue;
                TMarker * p = new TMarker(zc[ilayer][iwire][jwire],xc[ilayer][iwire][jwire],20);
                p->SetMarkerColor(color[iwire]);
                p->Draw("SAME");
		        TLatex * text = new TLatex(zc[ilayer][iwire][jwire],xc[ilayer][iwire][jwire]+5,Form("%d,%d",iwire,jwire));
		        text->SetTextColor(color[iwire]);
		        text->SetTextSize(0.02);
		        text->Draw("SAME");
            }
        }
        c->SaveAs(Form("layer%d_%d.pdf",ilayer,ilayer+1));
        c->SaveAs(Form("layer%d_%d.png",ilayer,ilayer+1));
    }
}
