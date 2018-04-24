#include "TString.h"
#include "TChain.h"
#include "TGraph.h"
#include "TF1.h"
#include <vector>

#include "header.h"

void printUsage(char * name);

int main(int argc, char** argv){
	if (argc<7){
	    printUsage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
    TString runname = argv[2];
    TString original = argv[3];
    int Niters = (int)strtol(argv[4],NULL,10);
	int geoSetup = (int)strtol(argv[5],NULL,10);
    int isMC = (int)strtol(argv[6],NULL,10);

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
    for (int iter = 0; iter<=Niters; iter++){
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

	// get offsets from MC
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NLAY; wid++){
            off[lid][wid] = 0;
        }
	}
	if (isMC){
		TChain * ichain_off = new TChain("t","t");
		ichain_off->Add(Form("Input/wire-offset.%d.root",runNo));
		if (!ichain_off->GetEntries()) {
		    printf("Cannot find Input/wire-offset.%d.root\n",runNo);
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
	ichain_wp->Add(Form("info/wire-position.%d.%s.root",runNo,original.Data()));
    if (!ichain_wp->GetEntries()) {
        printf("Cannot find info/wire-position.%d.%s.root\n",runNo,original.Data());
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
	for (int iter = 1; iter<=Niters; iter++){
		ichain_wp = new TChain("t","t");
		ichain_wp->Add(Form("info/wire-position.%d.%s.i%d.root",runNo,runname.Data(),iter));
		if (!ichain_wp->GetEntries()) {
            printf("Cannot find info/wire-position.%d.%s.i%d.root, will use original one\n",runNo,runname.Data(),iter);
		    ichain_wp->Add(Form("info/wire-position.%d.%s.root",runNo,original.Data()));
            if (!ichain_wp->GetEntries()) {
                printf("Cannot find info/wire-position.%d.%s.root\n",runNo,original.Data());
                return -1;
            }
        }
		ichain_wp->SetBranchAddress("xro",&wp_x);
		ichain_wp->SetBranchAddress("l",&wp_lid);
		ichain_wp->SetBranchAddress("w",&wp_wid);
		for (int i = 0; i<ichain_wp->GetEntries(); i++){
			ichain_wp->GetEntry(i);
            if (wp_lid<0||wp_lid>=NLAY||wp_wid<0||wp_wid>=NCEL) continue;
			printf("iter %d, [%d,%d], %.6e-%.6e=%.6e==0?%d\n",iter,wp_lid,wp_wid,wp_x,x[wp_lid][wp_wid],wp_x-x[wp_lid][wp_wid],wp_x-x[wp_lid][wp_wid]==0);
			if (wp_x-x[wp_lid][wp_wid]&&!changed[wp_lid][wp_wid]){
				changed[wp_lid][wp_wid] = true;
				printf(" check wire [%d,%d]\n",wp_lid,wp_wid);
			}
		}
        delete ichain_wp;
	}

	// get delta X
	for (int iter = 0; iter<=Niters; iter++){
		ichain_wp = new TChain("t","t");
		if (iter==0)
            ichain_wp->Add(Form("info/wire-position.%d.%s.root",runNo,original.Data()));
        else
            ichain_wp->Add(Form("info/wire-position.%d.%s.i%d.root",runNo,runname.Data(),iter));
		if (!ichain_wp->GetEntries()) {
            printf("Cannot find info/wire-position.%d.%s.i%d.root, will use original one\n",runNo,runname.Data(),iter);
		    ichain_wp->Add(Form("info/wire-position.%d.%s.root",runNo,original.Data()));
            if (!ichain_wp->GetEntries()) {
                printf("Cannot find info/wire-position.%d.%s.root\n",runNo,original.Data());
                return -1;
            }
        }
		ichain_wp->SetBranchAddress("xro",&wp_x);
		ichain_wp->SetBranchAddress("l",&wp_lid);
		ichain_wp->SetBranchAddress("w",&wp_wid);
		for (int i = 0; i<ichain_wp->GetEntries(); i++){
			ichain_wp->GetEntry(i);
            if (wp_lid<0||wp_lid>=NLAY||wp_wid<0||wp_wid>=NCEL) continue;
			if (changed[wp_lid][wp_wid]){
				iter_deltaX[iter][wp_lid][wp_wid] = wp_x-x[wp_lid][wp_wid]-off[wp_lid][wp_wid];
				printf("iter_deltaX[%d][%d][%d] = %.7e-%.7e-%.3e = %.3e\n",iter,wp_lid,wp_wid,wp_x,x[wp_lid][wp_wid],off[wp_lid][wp_wid],iter_deltaX[iter][wp_lid][wp_wid]);
			}
		}
        delete ichain_wp;
	}

	// get offset
	for (int iter = 1; iter<=Niters; iter++){
		TChain * ichain_off = new TChain("t","t");
		ichain_off->Add(Form("info/offset.%d.%s.i%d.root",runNo,runname.Data(),iter));
		if (!ichain_off->GetEntries()){
		    printf("Cannot find info/offset.%d.%s.i%d.root, will assume 0\n",runNo,runname.Data(),iter);
		    continue;
		}
		double off_d;
		int off_lid;
		int off_wid;
		ichain_off->SetBranchAddress("d",&off_d);
		ichain_off->SetBranchAddress("lid",&off_lid);
		ichain_off->SetBranchAddress("wid",&off_wid);
		for (int i = 0; i<ichain_off->GetEntries(); i++){
			ichain_off->GetEntry(i);
            if (off_lid<0||off_lid>=NLAY||off_wid<0||off_wid>=NCEL) continue;
			if (changed[off_lid][off_wid]){
				iter_offset[iter][off_lid][off_wid]= off_d;
			}
		}
        delete ichain_off;
	}

    // prepare graph and function to fit deltaX:layerID
    TGraph * gr_dxLid = new TGraph();
    TF1    * f_dxLid = new TF1("f_dxLid","pol1",0,NLAY);
    // get chi2, efficiency and Nothers
	for (int iter = 1; iter<=Niters; iter++){
	    for (int lid = 1; lid<NLAY; lid++){
            TChain * ichain_ana = new TChain("t","t");
            ichain_ana->Add(Form("root/t_%d.%s.i%d.layer%d.root",runNo,runname.Data(),iter,lid));
            if (!ichain_ana->GetEntries()){
                printf("Cannot find root/t_%d.%s.i%d.layer%d.root\n",runNo,runname.Data(),iter,lid);
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
            if (isMC) ichain_ana->SetBranchAddress("slzmc",&slzmc);
            if (isMC) ichain_ana->SetBranchAddress("inzmc",&inzmc);
            ichain_ana->SetBranchAddress("slx0",&slx);
            ichain_ana->SetBranchAddress("inx0",&inx);
            if (isMC) ichain_ana->SetBranchAddress("slxmc",&slxmc);
            if (isMC) ichain_ana->SetBranchAddress("inxmc",&inxmc);
            ichain_ana->SetBranchAddress("layerID",&layerID);
            ichain_ana->SetBranchAddress("wireID",&wireID);
            ichain_ana->SetBranchAddress("sel0",&sel);
            if (isMC) ichain_ana->SetBranchAddress("driftDmc",&driftDmc);
            ichain_ana->SetBranchAddress("driftD0",&driftD);
            ichain_ana->SetBranchAddress("fitD0",&fitD);
            int nEntries = ichain_ana->GetEntries();
            printf("Loading file root/t_%d.%s.i%d.layer%d.root, %d Entries\n",runNo,runname.Data(),iter,lid,nEntries);
            for (int iEntry = 0; iEntry<nEntries; iEntry++){
                ichain_ana->GetEntry(iEntry);
                int nHits = layerID->size();
                if (iEntry%10000==0){
                    printf("%d, %d hits, chi2 = %.3e, nHitsS = %d, slz = %.3e\n",iEntry,nHitsS,chi2,nHitsS,slz);
                    fflush(stdout);
                }
                if (chi2>0.5||nHitsS<6) continue;
                if (geoSetup==0&&fabs(slz)>0.15) continue;
                if (geoSetup==1&&fabs(inz)>24) continue;
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
                    if (isMC) tdriftDmc = (*driftDmc)[ihit];
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

    // print out the result
    printf("    iter/I lid/I wid/I dX dXo dXmc dXfit offset chi2 slx inx slxoff slxmc inxmc slz inz slzmc inzmc n\n");
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            if (!changed[lid][wid]) continue;
            printf("[%d,%d]:\n",lid,wid);
            for (int iter = 1; iter<=Niters; iter++){
                printf("  =>%d %d %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %d\n",
                        iter,
                        lid,wid,
                        iter_deltaX[iter-1][lid][wid],
                        iter_deltaXothers[iter-1][lid][wid],
                        iter_deltaXmc[iter][lid][wid],
                        iter_deltaXfit[iter][lid][wid],
                        iter_offset[iter][lid][wid],
                        iter_chi2[iter][lid][wid],
                        iter_slx[iter][lid][wid],
                        iter_inx[iter][lid][wid],
                        iter_slxoff[iter][lid][wid],
                        iter_slxmc[iter][lid][wid],
                        iter_inxmc[iter][lid][wid],
                        iter_slz[iter][lid][wid],
                        iter_inz[iter][lid][wid],
                        iter_slzmc[iter][lid][wid],
                        iter_inzmc[iter][lid][wid],
                        iter_nGood[iter][lid][wid]
                        );
                double deltaXothers = 0;
                for (int ljd = 0; ljd<NLAY; ljd++){
                    for (int wjd = 0; wjd<NCEL; wjd++){
                        int nEntries = iter_Nothers[iter][lid][wid][ljd][wjd];
                        deltaXothers += iter_deltaX[iter-1][ljd][wjd]/nEntries;
                        if (!nEntries) continue;
                        printf("            %d,%d: %d entries, deltaX = %.3e\n",ljd,wjd,nEntries,iter_deltaX[iter-1][ljd][wjd]);
                    }
                }
            }
        }
    }

    return 0;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [runname] [originalName] [Niter] [geoSetup: 0, normal;1, finger] [isMC]\n",name);
}
