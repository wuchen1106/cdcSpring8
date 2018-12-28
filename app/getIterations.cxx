#include "TString.h"
#include "TChain.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include <vector>
#include <math.h>

#include "Log.hxx"
#include "header.hxx"

void print_usage(char * prog_name);

int main(int argc, char** argv){
	int m_runNo = 0;
    TString m_runname = "currun";
    TString m_orirunname = "ori";
    int m_Niters = 0;
	int m_geoSetup = 0;
    bool m_isMC = false;

    int    opt_result;
	while((opt_result=getopt(argc,argv,"R:I:G:M"))!=-1){
		switch(opt_result){
			/* INPUTS */
			case 'R':
			    m_runNo = atoi(optarg);
                printf("Run number set to %d\n",m_runNo);
				break;
			case 'I':
			    m_Niters = atoi(optarg);
                printf("Number of iterations set to %d\n",m_Niters);
				break;
			case 'G':
			    m_geoSetup = atoi(optarg);
                printf("Geometry set to %d\n",m_geoSetup);
				break;
			case 'M':
			    m_isMC = true;
                printf("Turn to use MC mode\n");
				break;
			case '?':
				printf("Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
			case 'h':
			default:
				print_usage(argv[0]);
				return 1;
		}
	}

	if (argc-optind<2){
	    print_usage(argv[0]);
		return -1;
    }
    m_orirunname = argv[optind++];
    m_runname= argv[optind++];

	// what do we want
    double iter_slx[NITERSMAX][NLAY][NCEL]; // mean slx of chosen events of each iteration in each layer run of each wire 
    double iter_inx[NITERSMAX][NLAY][NCEL]; // mean inx of chosen events of each iteration in each layer run of each wire 
    double iter_slxoff[NITERSMAX][NLAY][NCEL]; // mean slx offset of chosen events of each iteration in each layer run of each wire 
    double iter_slxmc[NITERSMAX][NLAY][NCEL]; // mean slxmc of chosen events of each iteration in each layer run of each wire 
    double iter_inxmc[NITERSMAX][NLAY][NCEL]; // mean inxmc of chosen events of each iteration in each layer run of each wire 
    double iter_slz[NITERSMAX][NLAY][NCEL]; // mean slz of chosen events of each iteration in each layer run of each wire 
    double iter_inz[NITERSMAX][NLAY][NCEL]; // mean inz of chosen events of each iteration in each layer run of each wire 
    double iter_slzmc[NITERSMAX][NLAY][NCEL]; // mean slzmc of chosen events of each iteration in each layer run of each wire 
    double iter_inzmc[NITERSMAX][NLAY][NCEL]; // mean inzmc of chosen events of each iteration in each layer run of each wire 
    double iter_chi2[NITERSMAX][NLAY][NCEL]; // mean chi2 of chosen events of each iteration in each layer run of each wire 
    int iter_nGood[NITERSMAX][NLAY][NCEL]; // number of chosen events of each iteration in each layer run of each wire
    double iter_offset[NITERSMAX][NLAY][NCEL]; // offset according to residual distribution of each wire in each iteration
    double iter_deltaX[NITERSMAX][NLAY][NCEL]; // deltaX according to true value of each wire in each iteration
    double iter_deltaXothers[NITERSMAX][NLAY][NCEL]; // mean deltaX of other wires used of each iteration in each layer run of each wire 
    double iter_deltaXmc[NITERSMAX][NLAY][NCEL]; // mean deltaX of driftD-driftDmc of each iteration in each layer run of each wire 
    double iter_deltaXfit[NITERSMAX][NLAY][NCEL]; // mean deltaX of fitD-fitDmc of each iteration in each layer run of each wire 
    int iter_Nothers[NITERSMAX][NLAY][NCEL][NLAY][NCEL]; // number of entries from other wires used in this fitting in each iteration
	double off[NLAY][NCEL]; // real offset set in MC
	double x[NLAY][NCEL]; // designed wire position
	bool changed[NLAY][NCEL]; // is this wire changed during the iteration?

    // preset
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            changed[lid][wid] = false;
            off[lid][wid] = 0;
            x[lid][wid] = 0;
        }
    }
    for (int iter = 0; iter<=m_Niters; iter++){
        for (int lid = 0; lid<NLAY; lid++){
            for (int wid = 0; wid<NCEL; wid++){
                iter_slx[iter][lid][wid] = 0;
                iter_inx[iter][lid][wid] = 0;
                iter_slxoff[iter][lid][wid] = 0;
                iter_slxmc[iter][lid][wid] = 0;
                iter_inxmc[iter][lid][wid] = 0;
                iter_slz[iter][lid][wid] = 0;
                iter_inz[iter][lid][wid] = 0;
                iter_slzmc[iter][lid][wid] = 0;
                iter_inzmc[iter][lid][wid] = 0;
                iter_chi2[iter][lid][wid] = 0;
                iter_nGood[iter][lid][wid] = 0;
                iter_offset[iter][lid][wid] = 0;
                iter_deltaX[iter][lid][wid] = 0;
                iter_deltaXothers[iter][lid][wid] = 0;
                iter_deltaXmc[iter][lid][wid] = 0;
                iter_deltaXfit[iter][lid][wid] = 0;
                for (int ljd = 0; ljd<NLAY; ljd++){
                    for (int wjd = 0; wjd<NCEL; wjd++){
                        iter_Nothers[iter][lid][wid][ljd][wjd] = 0;
                    }
                }
            }
        }
    }

    TString HOME=getenv("CDCS8WORKING_DIR");
    TString filename = "";

	// get offsets from MC
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NLAY; wid++){
            off[lid][wid] = 0;
        }
	}
	if (m_isMC){
		TChain * ichain_off = new TChain("t","t");
		filename = HOME+Form("Input/wire-offset.%d.root",m_runNo);
		ichain_off->Add(filename);
		if (!ichain_off->GetEntries()) {
            MyError(Form("Cannot find \"%.s\"!\n",filename.Data()));
		    return -1;
        }
		double off_delta;
		int off_lid;
		int off_wid;
		ichain_off->SetBranchAddress("l",&off_lid);
		ichain_off->SetBranchAddress("w",&off_wid);
		ichain_off->SetBranchAddress("d",&off_delta);
		for (int i = 0; i<ichain_off->GetEntries(); i++){
			ichain_off->GetEntry(i);
            if (off_lid<0||off_lid>=NLAY||off_wid<0||off_wid>=NCEL) continue;
			off[off_lid][off_wid] = off_delta;
		}
	}

	// get wire positions
	TChain * ichain_wp = new TChain("t","t");
    filename = HOME+"Input/wire-position.root";
    ichain_wp->Add(filename);
    if (!ichain_wp->GetEntries()) {
        MyError(Form("Cannot find \"%.s\"!\n",filename.Data()));
        return -1;
    }
	double wp_x;
	int wp_lid;
	int wp_wid;
	ichain_wp->SetBranchAddress("xro",&wp_x);
	ichain_wp->SetBranchAddress("l",&wp_lid);
	ichain_wp->SetBranchAddress("w",&wp_wid);
	for (int i = 0; i<ichain_wp->GetEntries(); i++){
		ichain_wp->GetEntry(i);
        if (wp_lid<0||wp_lid>=NLAY||wp_wid<0||wp_wid>=NCEL) continue;
		x[wp_lid][wp_wid] = wp_x;
	}
	delete ichain_wp;

	// check which wire is changed
	for (int iter = 1; iter<=m_Niters; iter++){
		ichain_wp = new TChain("t","t");
        filename = HOME+Form("info/offset.%d.%s.i%d.root",m_runNo,m_runname.Data(),iter);
		ichain_wp->Add(filename);
		if (!ichain_wp->GetEntries()) {
            MyError(Form("Cannot find \"%.s\"!\n",filename.Data()));
            return -1;
        }
        double off_adjust;
		ichain_wp->SetBranchAddress("adjust",&off_adjust);
		ichain_wp->SetBranchAddress("l",&wp_lid);
		ichain_wp->SetBranchAddress("w",&wp_wid);
		for (int i = 0; i<ichain_wp->GetEntries(); i++){
			ichain_wp->GetEntry(i);
            if (wp_lid<0||wp_lid>=NLAY||wp_wid<0||wp_wid>=NCEL) continue;
			if (off_adjust){
				changed[wp_lid][wp_wid] = true;
			}
		}
        delete ichain_wp;
	}

	// get delta X
	for (int iter = 0; iter<=m_Niters; iter++){
		ichain_wp = new TChain("t","t");
        filename = HOME+Form("/info/offset.%d.%s.i%d.root",m_runNo,m_runname.Data(),iter);
		ichain_wp->Add(filename);
        if (!ichain_wp->GetEntries()) {
            continue;
        }
        double off_adjust;
        double off_d;
		ichain_wp->SetBranchAddress("d",&off_d);
		ichain_wp->SetBranchAddress("adjust",&off_adjust);
		ichain_wp->SetBranchAddress("l",&wp_lid);
		ichain_wp->SetBranchAddress("w",&wp_wid);
		for (int i = 0; i<ichain_wp->GetEntries(); i++){
			ichain_wp->GetEntry(i);
            if (wp_lid<0||wp_lid>=NLAY||wp_wid<0||wp_wid>=NCEL) continue;
			if (changed[wp_lid][wp_wid]){
				iter_offset[iter][wp_lid][wp_wid] = off_d;
				iter_deltaX[iter][wp_lid][wp_wid] = off_adjust;
			}
		}
        delete ichain_wp;
	}

    // prepare graph and function to fit deltaX:layerID
    TGraph * gr_dxLid = new TGraph();
    TF1    * f_dxLid = new TF1("f_dxLid","pol1",0,NLAY);
    // get chi2, efficiency and Nothers
    if (m_isMC){
        for (int iter = 1; iter<=m_Niters; iter++){
            for (int lid = 1; lid<NLAY; lid++){
                TChain * ichain_ana = new TChain("t","t");
                ichain_ana->Add(Form("root/tracks/t_%d.%s.i%d.layer%d.root",m_runNo,m_runname.Data(),iter,lid));
                if (!ichain_ana->GetEntries()){
                    printf("Cannot find root/tracks/t_%d.%s.i%d.layer%d.root\n",m_runNo,m_runname.Data(),iter,lid);
                    continue;
                }
                double chi2;
                int nHitsS;
                double slx;
                double inx;
                double slxmc = 0;
                double inxmc = 0;
                double slz;
                double inz;
                double slzmc = 0;
                double inzmc = 0;
                std::vector<int> *    layerID = 0;
                std::vector<int> *    wireID = 0;
                std::vector<int> *    sel = 0;
                std::vector<double> * driftDmc = 0;
                std::vector<double> * driftD = 0;
                std::vector<double> * fitD = 0;
                ichain_ana->SetBranchAddress("chi20",&chi2);
                ichain_ana->SetBranchAddress("nHitsS0",&nHitsS);
                ichain_ana->SetBranchAddress("slz0",&slz);
                ichain_ana->SetBranchAddress("inz0",&inz);
                if (m_isMC) ichain_ana->SetBranchAddress("slzmc",&slzmc);
                if (m_isMC) ichain_ana->SetBranchAddress("inzmc",&inzmc);
                ichain_ana->SetBranchAddress("slx0",&slx);
                ichain_ana->SetBranchAddress("inx0",&inx);
                if (m_isMC) ichain_ana->SetBranchAddress("slxmc",&slxmc);
                if (m_isMC) ichain_ana->SetBranchAddress("inxmc",&inxmc);
                ichain_ana->SetBranchAddress("layerID",&layerID);
                ichain_ana->SetBranchAddress("wireID",&wireID);
                ichain_ana->SetBranchAddress("sel0",&sel);
                if (m_isMC) ichain_ana->SetBranchAddress("driftDmc",&driftDmc);
                ichain_ana->SetBranchAddress("driftD0",&driftD);
                ichain_ana->SetBranchAddress("fitD0",&fitD);
                int nEntries = ichain_ana->GetEntries();
                printf("Loading file root/tracks/t_%d.%s.i%d.layer%d.root, %d Entries\n",m_runNo,m_runname.Data(),iter,lid,nEntries);
                for (int iEntry = 0; iEntry<nEntries; iEntry++){
                    ichain_ana->GetEntry(iEntry);
                    int nHits = layerID->size();
                    if (iEntry%10000==0){
                        printf("%d, %d hits, chi2 = %.3e, nHitsS = %d, slz = %.3e\n",iEntry,nHitsS,chi2,nHitsS,slz);
                        fflush(stdout);
                    }
                    if (chi2>0.5||nHitsS<6) continue;
                    if (m_geoSetup==0&&fabs(slz)>0.15) continue;
                    if (m_geoSetup==1&&fabs(inz)>24) continue;
                    double minres = 1e9;
                    int has = 0;
                    int theWid = 0;
                    double theDD = 0;
                    double theDDmc = 0;
                    double theFitD = 0;
                    for (int ihit = 0; ihit<nHits; ihit++){
                        int tlayerID = (*layerID)[ihit];
                        if (tlayerID!=lid) continue;
                        int twireID = (*wireID)[ihit];
                        double tfitD = (*fitD)[ihit];
                        double tdriftD = (*driftD)[ihit];
                        double tdriftDmc = 0;
                        if (m_isMC) tdriftDmc = (*driftDmc)[ihit];
                        if (fabs(tfitD-tdriftD)<fabs(minres)){ // no cut for test layer!
                            minres = tfitD-tdriftD;
                            theWid = twireID;
                            theDD = tdriftD;
                            theFitD = tfitD;
                            theDDmc = tdriftDmc;
                            has = 1;
                        }
                    }
                    if (has==0||theWid<0||theWid>=NCEL) continue;
                    if (fabs(theDD)<2||fabs(theDD)>6) continue;
                    if (fabs(theFitD-theDD)>1) continue;
                    double tdeltaXothers = 0;
                    int nHitsUsed = 0;
                    for (int ihit = 0; ihit<nHits; ihit++){
                        if (!(*sel)[ihit]) continue;
                        int tlid = (*layerID)[ihit];
                        int twid = (*wireID)[ihit];
                        if (tlid<0||tlid>=NLAY||twid<0||twid>=NCEL) continue;
                        iter_Nothers[iter][lid][theWid][tlid][twid]++;
                        gr_dxLid->SetPoint(nHitsUsed,tlid,iter_deltaX[iter-1][tlid][twid]);
                        nHitsUsed++;
                    }
                    gr_dxLid->Set(nHitsUsed);
                    gr_dxLid->Fit("f_dxLid","QN0G","");
                    tdeltaXothers+=f_dxLid->Eval(lid);
                    iter_deltaXmc[iter][lid][theWid]+=theDD-theDDmc;
                    iter_deltaXfit[iter][lid][theWid]+=theFitD-(theDDmc-iter_deltaX[iter-1][lid][theWid])-tdeltaXothers;
                    iter_deltaXothers[iter-1][lid][theWid]+=tdeltaXothers;
                    iter_chi2[iter][lid][theWid]+=chi2;
                    iter_slx[iter][lid][theWid]+=slx;
                    iter_inx[iter][lid][theWid]+=inx;
                    iter_slxoff[iter][lid][theWid]+=f_dxLid->GetParameter(1);
                    iter_slxmc[iter][lid][theWid]+=slxmc;
                    iter_inxmc[iter][lid][theWid]+=inxmc;
                    iter_slz[iter][lid][theWid]+=slz;
                    iter_inz[iter][lid][theWid]+=inz;
                    iter_slzmc[iter][lid][theWid]+=slzmc;
                    iter_inzmc[iter][lid][theWid]+=inzmc;
                    iter_nGood[iter][lid][theWid]++;
                }
                for (int wid = 0; wid<NCEL; wid++){
                    iter_slx[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_inx[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_slxoff[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_slxmc[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_inxmc[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_slz[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_inz[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_slzmc[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_inzmc[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_chi2[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_deltaXothers[iter-1][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_deltaXmc[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                    iter_deltaXfit[iter][lid][wid]/=iter_nGood[iter][lid][wid];
                }
                delete ichain_ana;
            }
        }
    }

    // Output
    TFile * ofile = new TFile(HOME+Form("/info/iter.%d.%s.root",m_runNo,m_runname.Data()));

    // get histograms
    TH2D * hOffsets[NITERSMAX] = {0};
    TH2D * hAdjusts[NITERSMAX] = {0};
    for (int iter = 1; iter<=m_Niters; iter++){
        TCanvas * canv = new TCanvas(Form("Iter%d",iter),"canv",1024,800);
        canv->Divide(2,1);
        hOffsets[iter] = new TH2D(Form("hOffsets%d",iter),"Adjustments",11,-0.5,10.5,8,0.5,8.5);
        hAdjusts[iter] = new TH2D(Form("hAdjusts%d",iter),"Adjustments",11,-0.5,10.5,8,0.5,8.5);
        for (int lid = 0; lid<NLAY; lid++){
            for (int wid = 0; wid<NCEL; wid++){
                hOffsets[iter]->Fill(wid,lid,iter_offset[iter][lid][wid]);
                hAdjusts[iter]->Fill(wid,lid,iter_deltaX[iter][lid][wid]);
            }
        }
        canv->cd(1);
        hOffsets[iter]->Draw("COLZTEXT");
        canv->cd(2);
        hAdjusts[iter]->Draw("COLZTEXT");
        canv->SaveAs(Form("Wire_%d.%s.iter%d.png",m_runNo,m_runname.Data(),iter));
        hOffsets[iter]->Write();
        hAdjusts[iter]->Write();
    }

    int    o_lid;
    int    o_wid;
    double o_deltaX;
    double o_deltaXmc;
    double o_deltaXfit;
    double o_deltaXothers;
    double o_offset;
    double o_chi2;
    double o_slx;
    double o_inx;
    double o_slxoff;
    double o_slxmc;
    double o_slz;
    double o_inz;
    double o_slzmc;
    double o_inzmc;
    int    o_nGood;
    TTree * otree = new TTree("t","t");
    otree->Branch("lid",&o_lid);
    otree->Branch("wid",&o_wid);
    otree->Branch("deltaX",&o_deltaX);
    otree->Branch("deltaXmc",&o_deltaXmc);
    otree->Branch("deltaXfit",&o_deltaXfit);
    otree->Branch("deltaXothers",&o_deltaXothers);
    otree->Branch("offset",&o_offset);
    otree->Branch("chi2",&o_chi2);
    otree->Branch("slx",&o_slx);
    otree->Branch("inx",&o_inx);
    otree->Branch("slxoff",&o_slxoff);
    otree->Branch("slxmc",&o_slxmc);
    otree->Branch("slz",&o_slz);
    otree->Branch("inz",&o_inz);
    otree->Branch("slzmc",&o_slzmc);
    otree->Branch("inzmc",&o_inzmc);
    otree->Branch("nGood",&o_nGood);
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            if (!changed[lid][wid]) continue;
            for (int iter = 1; iter<=m_Niters; iter++){
                 o_lid = lid;
                 o_wid = wid;
                 o_deltaX = iter_deltaX[iter][lid][wid];
                 o_deltaXmc = iter_deltaXmc[iter][lid][wid];
                 o_deltaXfit = iter_deltaXfit[iter][lid][wid];
                 o_deltaXfit = iter_deltaXfit[iter][lid][wid];
                 o_offset = iter_offset[iter][lid][wid];
                 o_chi2 = iter_chi2[iter][lid][wid];
                 o_slx = iter_slx[iter][lid][wid];
                 o_inx = iter_inx[iter][lid][wid];
                 o_slxoff = iter_slxoff[iter][lid][wid];
                 o_slxmc = iter_slxmc[iter][lid][wid];
                 o_slz = iter_slz[iter][lid][wid];
                 o_inz = iter_inz[iter][lid][wid];
                 o_slzmc = iter_slzmc[iter][lid][wid];
                 o_inzmc = iter_inzmc[iter][lid][wid];
                 o_nGood = iter_nGood[iter][lid][wid];
                 otree->Fill();
                /*
                double deltaXothers = 0;
                for (int ljd = 0; ljd<NLAY; ljd++){
                    for (int wjd = 0; wjd<NCEL; wjd++){
                        int nEntries = iter_Nothers[iter][lid][wid][ljd][wjd];
                        deltaXothers += iter_deltaX[iter-1][ljd][wjd]/nEntries;
                        if (!nEntries) continue;
                    }
                }
                */
            }
        }
    }
    otree->Write();
    ofile->Close();

    return 0;
}

void print_usage(char * prog_name){
	fprintf(stderr,"Usage %s [options] orirunname runname\n",prog_name);
	fprintf(stderr,"[options]\n");
	fprintf(stderr,"\t -R <run>\n");
	fprintf(stderr,"\t\t Run number set to run\n");
	fprintf(stderr,"\t -I <Niter>\n");
	fprintf(stderr,"\t\t Number of iterations set to Niter\n");
	fprintf(stderr,"\t -G <geo>\n");
	fprintf(stderr,"\t\t Geometry set to geo\n");
	fprintf(stderr,"\t -M\n");
	fprintf(stderr,"\t\t Turn to MC mode (by default it's data mode)\n");
}
