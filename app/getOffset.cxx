#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

#include "XTAnalyzer.h"

#define NLAY  9
#define NCEL  11
#define NCAND 4

void printUsage(char * name);
int getHitType(int type,bool isRight);

int main(int argc, char** argv){

	if (argc<4){
	    printUsage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
    TString prerunname  = argv[2];
    TString runname = argv[3];
	int xtType = 2;
    if (argc>=5)
		xtType = (int)strtol(argv[4],NULL,10);
	int geoSetup = 0; // 0: normal scintillator; 1: finger scintillator
    if (argc>=6)
        geoSetup = (int)strtol(argv[5],NULL,10);
    int saveHists = 0;
    if (argc>=7)
        saveHists = (int)strtol(argv[6],NULL,10);
    int inputType = 0; // by defualt it's data
    if (argc>=8)
        inputType = (int)strtol(argv[7],NULL,10);
    TString offsetFile = "";
    if (argc>=9)
        offsetFile = argv[8];
    int debugLevel = 0;
    if (argc>=10)
        debugLevel = (int)strtol(argv[9],NULL,10);
    printf("##############Input %d Parameters##################\n",argc);
    printf("runNo       = %d\n",runNo);
    printf("prerunname  = \"%s\"\n",prerunname.Data());
    printf("runname     = \"%s\"\n",runname.Data());
    printf("geoSetup:     %s\n",geoSetup==0?"normal scintillator":"finger scintillator");
    printf("xtType:       %s\n",xtType==0?"asymmetric":(xtType==1?"symmetric, with offset":(xtType==2?"symmetric, thru 0":(xtType==3?"symmetric with nLHits==0":(xtType==4?"symmetric with smallest chi2a":(xtType==5?"symmetric with smallest chi2":"others?"))))));
    printf("save slice fittings? \"%s\"\n",saveHists?"yes":"no");
    printf("inputType   = %d, %s\n",inputType,inputType==0?"Real Data":"MC");
    printf("debug       = %d\n",debugLevel);
    printf("Offset file : \"%s\"\n",offsetFile.Data());
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

	// prepare offset
	TH1D * h_off[NLAY][NCAND];
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			h_off[lid][wid] = new TH1D(Form("h_off_%d_%d",lid,wid),Form("Offset of wire [%d,%d]",lid,wid),128,-1,1);
		}
	}

    // input file
    int triggerNumber;
    int nHits;
    int nHitsG;
    std::vector<int> *    i_layerID = 0;
    std::vector<int> *    i_wireID = 0;
    std::vector<double> * i_driftT = 0;
    std::vector<double> * i_driftDmc = 0;
    std::vector<int> *    i_type = 0;
    std::vector<int> *    i_np = 0;
    std::vector<int> *    i_ip = 0;
    std::vector<int> *    i_clk = 0;
    std::vector<int> *    i_width = 0;
    std::vector<int> *    i_peak = 0;
    std::vector<int> *    i_height = 0;
    std::vector<int> *    i_mpn = 0;
    std::vector<int> *    i_mpi = 0;
    std::vector<int> *    i_rank = 0;
    bool has_rank = false;
    std::vector<double> * i_aa = 0;
    std::vector<double> * i_ped = 0;
    bool has_ped = false;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_driftD[NCAND] = {0};
    std::vector<double> * i_calD[NCAND] = {0};
    int npairs[NCAND];
    int isel[NCAND];
    int icom[NCAND];
    double iinx[NCAND];
    double islx[NCAND];
    double iinz[NCAND];
    double islz[NCAND];
    double chi2x[NCAND];
    double chi2z[NCAND];
    double chi2i[NCAND];
    int nHitsS[NCAND];
    double inx[NCAND];
    double slx[NCAND];
    double inz[NCAND];
    double slz[NCAND];
    double chi2[NCAND];
    double chi2p[NCAND];
    double chi2a[NCAND];
    double chi2mc[NCAND];
    double chi2pmc[NCAND];
    double chi2amc[NCAND];
    double inxmc;
    double inzmc;
    double slxmc;
    double slzmc;
    std::vector<double> * i_fitD[NCAND] = {0};
    std::vector<int> * i_sel[NCAND] = {0};

    // Loop in layers
	for (int lid = 1; lid<=8; lid++){
        if (debugLevel>0) {printf("In Layer %d: preparing input TChain\n",lid);fflush(stdout);}
		TChain * ichain = new TChain("t","t");
		ichain->Add(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),lid));
        ichain->SetBranchAddress("triggerNumber",&triggerNumber);
        ichain->SetBranchAddress("nHits",&nHits);
        ichain->SetBranchAddress("nHitsG",&nHitsG);
        ichain->SetBranchAddress("layerID",&i_layerID);
        ichain->SetBranchAddress("wireID",&i_wireID);
        ichain->SetBranchAddress("driftT",&i_driftT);
        if (inputType) ichain->SetBranchAddress("driftDmc",&i_driftDmc);
        ichain->SetBranchAddress("type",&i_type);
        ichain->SetBranchAddress("np",&i_np);
        ichain->SetBranchAddress("ip",&i_ip);
        ichain->SetBranchAddress("clk",&i_clk);
        ichain->SetBranchAddress("width",&i_width);
        ichain->SetBranchAddress("peak",&i_peak);
        ichain->SetBranchAddress("height",&i_height);
        ichain->SetBranchAddress("mpn",&i_mpn);
        ichain->SetBranchAddress("mpi",&i_mpi);
        has_rank = (ichain->SetBranchAddress("rank",&i_rank)==0);
        ichain->SetBranchAddress("aa",&i_aa);
        has_ped = (ichain->SetBranchAddress("ped",&i_ped)==0);
		ichain->SetBranchAddress("sum",&i_sum);
		for (int iCand = 0; iCand<NCAND; iCand++){
			ichain->SetBranchAddress(Form("driftD%d",iCand),&(i_driftD[iCand]));
			ichain->SetBranchAddress(Form("npairs%d",iCand),&(npairs[iCand]));
			ichain->SetBranchAddress(Form("isel%d",iCand),&(isel[iCand]));
			ichain->SetBranchAddress(Form("icom%d",iCand),&(icom[iCand]));
			ichain->SetBranchAddress(Form("islx%d",iCand),&(islx[iCand]));
			ichain->SetBranchAddress(Form("islz%d",iCand),&(islz[iCand]));
			ichain->SetBranchAddress(Form("iinx%d",iCand),&(iinx[iCand]));
			ichain->SetBranchAddress(Form("iinz%d",iCand),&(iinz[iCand]));
			ichain->SetBranchAddress(Form("chi2x%d",iCand),&(chi2x[iCand]));
			ichain->SetBranchAddress(Form("chi2z%d",iCand),&(chi2z[iCand]));
			ichain->SetBranchAddress(Form("chi2i%d",iCand),&(chi2i[iCand]));
			ichain->SetBranchAddress(Form("calD%d",iCand),&(i_calD[iCand]));
			ichain->SetBranchAddress(Form("nHitsS%d",iCand),&(nHitsS[iCand]));
			ichain->SetBranchAddress(Form("slx%d",iCand),&(slx[iCand]));
			ichain->SetBranchAddress(Form("slz%d",iCand),&(slz[iCand]));
			ichain->SetBranchAddress(Form("inx%d",iCand),&(inx[iCand]));
			ichain->SetBranchAddress(Form("inz%d",iCand),&(inz[iCand]));
			ichain->SetBranchAddress(Form("chi2%d",iCand),&(chi2[iCand]));
			ichain->SetBranchAddress(Form("chi2p%d",iCand),&(chi2p[iCand]));
			ichain->SetBranchAddress(Form("chi2a%d",iCand),&(chi2a[iCand]));
			if (inputType){
                ichain->SetBranchAddress(Form("chi2mc%d",iCand),&(chi2mc[iCand]));
                ichain->SetBranchAddress(Form("chi2pmc%d",iCand),&(chi2pmc[iCand]));
                ichain->SetBranchAddress(Form("chi2amc%d",iCand),&(chi2amc[iCand]));
            }
			ichain->SetBranchAddress(Form("fitD%d",iCand),&(i_fitD[iCand]));
			ichain->SetBranchAddress(Form("sel%d",iCand),&(i_sel[iCand]));
		}
		if (inputType){
			ichain->SetBranchAddress("slxmc",&slxmc);
			ichain->SetBranchAddress("slzmc",&slzmc);
			ichain->SetBranchAddress("inxmc",&inxmc);
			ichain->SetBranchAddress("inzmc",&inzmc);
		}

        // Loop in events
        Long64_t N = ichain->GetEntries();
        if (N==0){
            fprintf(stderr,"WARNING: \"%s/root/t_%d.%s.layer%d.root\" is empty! Will ignore this layer.\n",HOME.Data(),runNo,runname.Data(),lid);
            continue;
        }
		int nSmallSumHits = 0;
		int nShadowedHits = 0;
		int nLateHits = 0;
		int nBoundaryHits = 0;
		int nSmallBoundaryHits = 0;
        if (debugLevel>0) {printf("Processing %d events\n",N);fflush(stdout);}
        for ( int iEntry = 0 ; iEntry<N; iEntry++){
            if (N%1000==0) printf("%d\n",N);
            if (debugLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

			// decide which candidate to use
			int theCand = 0;
			if (xtType==3){
				int nLateHitsMin = 1e9;
				for (int iCand = 0; iCand<NCAND; iCand++){
					nLateHits = 0;
					for (int ihit = 0; ihit<nHits; ihit++){
						int ip = 0;
						for (int jhit = ihit-1; jhit>0; jhit--){
							if ((*i_layerID)[jhit]!=(*i_layerID)[ihit]) break;
							int type = getHitType((*i_type)[jhit],(*i_fitD[iCand])[jhit]>=0);
							if (type<100) ip++;
						}
						if ((*i_sel[iCand])[ihit]==1){
							if(ip!=0)
								nLateHits++;
						}
					}
					if (nLateHits<nLateHitsMin){
						nLateHitsMin = nLateHits;
						theCand = iCand;
					}
				}
			}
            else if (xtType==4||xtType==5){
                double minchi2 = 1e9;
                int minNhitsS = 0;
				for (int iCand = 0; iCand<NCAND; iCand++){
                    if (xtType==4){
                        if ((minchi2>chi2a[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                            theCand = iCand;
                            minchi2 = chi2a[iCand];
                            minNhitsS = nHitsS[iCand];
                        }
                    }
                    else if (xtType==5){
                        if ((minchi2>chi2[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                            theCand = iCand;
                            minchi2 = chi2[iCand];
                            minNhitsS = nHitsS[iCand];
                        }
                    }
                }
            }

            // ignore events with bad fitting
            if (nHitsS[theCand]<7) continue;
			if (chi2[theCand]>2) continue;
            //if (nHitsG>nHitsS[theCand]) continue;
            if (geoSetup==1){
                if (fabs(inz[theCand])>24) continue;
            }
            else{
                if (fabs(slz[theCand])>0.15) continue;
            }

            if (debugLevel>=20) printf("  Good Event! Looping in %d hits\n",nHits);
            // find the closest hit in the test layer
            double minres = 1e9;
            bool has = false;
            int wireID;
            double driftD, driftT, fitD;
            // FIXME: test more cut
            bool hasBadHit = false;
            for (int ihit = 0; ihit<nHits; ihit++){
                int tlayerID = (*i_layerID)[ihit];
                int twireID = (*i_wireID)[ihit];
                double tfitD = (*i_fitD[theCand])[ihit];
                double tdriftD = (*i_driftD[theCand])[ihit];
            	if ((*i_sel[theCand])[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                if (tlayerID!=lid) continue;
				int ttype = getHitType((*i_type)[ihit],tfitD>=0);
                if (fabs(tfitD-tdriftD)<fabs(minres)){ // no cut for test layer!
                    minres = tfitD-tdriftD;
                    wireID = (*i_wireID)[ihit];
                    fitD = tfitD;
                    driftT = (*i_driftT)[ihit];
                    has = true;
                }
            }
            if (!has) continue; // no hits found in test layer
            //if (hasBadHit) continue;

            if (debugLevel>=20) printf("  Found hit! pushing to XTAnalyzer\n");
			// tell analyzer a new data point
            h_off[lid][wireID]->Fill(fitD-driftD);
        }
	}

	// output
	TFile * ofile = new TFile(offsetFile,"RECREATE");
	TTree * otree = new TTree("t","t");
	double o_off_delta;
	int o_off_lid;
	int o_off_wid;
	otree->Branch("d",&o_off_delta);
	otree->Branch("wid",&o_off_wid);
	otree->Branch("lid",&o_off_lid);
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			h_off[lid][wid]->Write();
			if (h_off[lid][wid]->GetEntries()<1000) continue;
			o_off_wid = wid;
			o_off_lid = lid;
			o_off_delta = h_off[lid][wid]->GetMean();
			otree->Fill();
		}
	}
	otree->Write();
	ofile->Close();

    return 0;
}

int getHitType(int type,bool isRight){
	int ttype = (type/10)%10;
	if (isRight){
		if (ttype==1||ttype==4) type-=ttype*10; // l- or l+
	}
	else{
		if (ttype==2||ttype==5) type-=ttype*10; // r- or r+
	}
	return type;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [runname] [geoSetup: 0, normal;1, finger] [offsetfile] [debug: 0;...]>\n",name);
}
