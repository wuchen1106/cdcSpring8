#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"

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
    int debugLevel = 0;
    if (argc>=9)
        debugLevel = (int)strtol(argv[8],NULL,10);
    printf("##############Input %d Parameters##################\n",argc);
    printf("runNo       = %d\n",runNo);
    printf("prerunname  = \"%s\"\n",prerunname.Data());
    printf("runname     = \"%s\"\n",runname.Data());
    printf("geoSetup:     %s\n",geoSetup==0?"normal scintillator":"finger scintillator");
    printf("xtType:       %s\n",xtType==0?"asymmetric":(xtType==1?"symmetric, with offset":(xtType==2?"symmetric, thru 0":(xtType==3?"symmetric with nLHits==0":(xtType==4?"symmetric with smallest chi2a":(xtType==5?"symmetric with smallest chi2":"others?"))))));
    printf("save slice fittings? \"%s\"\n",saveHists?"yes":"no");
    printf("inputType   = %d, %s\n",inputType,inputType==0?"Real Data":"MC");
    printf("debug       = %d\n",debugLevel);
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

	//===================Get run info============================
	TFile * if_run = new TFile(HOME+"/info/run-info.root");
	TTree * t_run = (TTree*) if_run->Get("t");
	int i_runNo, gasID, runGr, HV, THR;
	char runDu[128];
	double t00, t01, aacut, sumcut;
	t_run->SetBranchAddress("run_number",&i_runNo);
	t_run->SetBranchAddress("gas_mixture_id",&gasID);
	t_run->SetBranchAddress("hv_ch0",&HV);
	t_run->SetBranchAddress("recbe_th_input_bd0",&THR);
	t_run->SetBranchAddress("duration",&runDu);
	t_run->SetBranchAddress("run_grade",&runGr);
	t_run->SetBranchAddress("t00",&t00);
	t_run->SetBranchAddress("t01",&t01);
	t_run->SetBranchAddress("aa",&aacut);
	t_run->SetBranchAddress("sum",&sumcut);
	for(int i = 0; i<t_run->GetEntries(); i++){
		t_run->GetEntry(i);
		if (i_runNo == runNo) break;
	}
	double npair = 17.96;
	TString gastype = "He:C_{2}H_{4}(50:50)";
	if (gasID==1){
		gastype = "He:iC_{4}H_{10}(90:10)";
		npair = 27.96;
	}
	else if (gasID==2){
		gastype = "He:CH_{4}(80:20)";
		npair = 56.10;
	}
	TString duration = runDu;
	const char *sep = ":";
	char * durationSep = strtok(runDu,sep);
	double durationTime = 0;
	double timeunit = 3600;
	while(durationSep){
		durationTime += timeunit*strtol(durationSep,NULL,10);
		timeunit/=60;
		durationSep = strtok(NULL,sep);
	}
	printf("runNo#%d: %s, %d, %s, %d V, %d mV, %.0f sec\n",runNo,gastype.Data(),runGr,duration.Data(),HV,THR,durationTime);

	// get offset
	double off[NLAY][NCEL];
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			off[lid][wid] = 0;
		}
	}
	if (xtType==1){
		TChain * iChain_off = new TChain("t","t");
		iChain_off->Add(Form("%s/info/offset.%d.%s.root",HOME.Data(),runNo,runname.Data()));
		double i_off_delta;
		int i_off_lid;
		int i_off_wid;
		iChain_off->SetBranchAddress("d",&i_off_delta);
		iChain_off->SetBranchAddress("wid",&i_off_wid);
		iChain_off->SetBranchAddress("lid",&i_off_lid);
		int N = iChain_off->GetEntries();
		for (int i = 0; i<N; i++){
			iChain_off->GetEntry(i);
			if (i_off_lid>=0&&i_off_lid<NLAY&&i_off_wid>=0&&i_off_wid<NCEL)
				off[i_off_lid][i_off_wid] = i_off_delta;
		}
	}

    // get XT file of the previous run
    TFile * preXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),runNo,prerunname.Data()));

    // prepare new XT file for this run
    TFile * newXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),runNo,runname.Data()),"RECREATE");
    TTree * newXTTree = new TTree("t","t");
	double mX;
	double mT;
	int mLayerID;
	double mSig;
	double mChi2;
	int mEntries;
	int mType;
	newXTTree->Branch("x",&mX);
	newXTTree->Branch("t",&mT);
	newXTTree->Branch("lid",&mLayerID);
	newXTTree->Branch("sig",&mSig);
	newXTTree->Branch("chi2",&mChi2);
	newXTTree->Branch("n",&mEntries);
	newXTTree->Branch("type",&mType);

	// Prepare XTAnalyzer
	XTAnalyzer * fXTAnalyzer = new XTAnalyzer(gasID,debugLevel);

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

	// output file
    std::vector<double> * o_driftD = 0;
    std::vector<int>    * o_driftDs = 0;

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

		int saveEvenOdd = 0; if (lid==4) saveEvenOdd = 1; else if (lid==5) saveEvenOdd = -1;
		int statusInitialize = fXTAnalyzer->Initialize(Form("%d.%s.layer%d",runNo,runname.Data(),lid),lid,preXTFile,newXTFile,newXTTree,xtType,saveHists, lid==4, saveEvenOdd); // take the XT from layer 4 as default output XT (0). 4 as default even and 5 as default odd.
		if (statusInitialize){
			fprintf(stderr,"WARNING: something wrong with initializing XTAnalyzer for layer[%d], will ignore this layer!\n",lid);
			continue;
		}

        // Loop in events
        ichain->GetEntries();
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
            if (iEntry%10000==0) printf("%d\n",iEntry);
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
            if (nHitsS[theCand]<6) continue;
			if (chi2[theCand]>0.8) continue;
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
                double tfitD = (*i_fitD[theCand])[ihit]-off[tlayerID][twireID];
                double tdriftD = (*i_driftD[theCand])[ihit];
            	if ((*i_sel[theCand])[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                if (tlayerID!=lid) continue;
				int ttype = getHitType((*i_type)[ihit],tfitD>=0);
                //if (ttype<100&&fabs(tfitD-tdriftD)<fabs(minres)){
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
            fXTAnalyzer->Push(driftT,fitD);
        }
        if (debugLevel>0) printf("Starting XT analysis\n");
        // fit histograms/graphs, make plots, and save new xt file
        fXTAnalyzer->Process();

        // prepare for output ROOT file
        TFile * ofile = new TFile(Form("%s/root/ana_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),lid),"RECREATE");
        TTree * otree = new TTree("t","t");
		double minres = 1e9;
		double theDD = 1e9;
		double theDT = 1e9;
		int has = 0;
        int theWid = -1;
		double theSum = 0;
		double thePeak = 0;
		double theHeight = 0;
		double sum1st = 0;
		double dt1st = 0;
		int theIp = 0;
		int theMpi = 0;
		int theCand = 0;
        otree->Branch("triggerNumber",&triggerNumber);
        otree->Branch("res",&minres);
        otree->Branch("theDD",&theDD);
        otree->Branch("theDT",&theDT);
        otree->Branch("theWid",&theWid);
        otree->Branch("theSum",&theSum);
        otree->Branch("sum1st",&sum1st);
        otree->Branch("dt1st",&dt1st);
        otree->Branch("has",&has);
        otree->Branch("thePeak",&thePeak);
        otree->Branch("theHeight",&theHeight);
        otree->Branch("theIp",&theIp);
        otree->Branch("theMpi",&theMpi);
        otree->Branch("theCand",&theCand);
        otree->Branch("nSHits",&nShadowedHits);
        otree->Branch("nLHits",&nLateHits);
        otree->Branch("nSSHits",&nSmallSumHits);
        otree->Branch("nBHits",&nBoundaryHits);
        otree->Branch("nSBHits",&nSmallBoundaryHits);
        otree->Branch("nHits",&nHits);
        otree->Branch("nHitsG",&nHitsG);
        otree->Branch("layerID",&i_layerID);
        otree->Branch("wireID",&i_wireID);
        otree->Branch("driftT",&i_driftT);
        if (inputType) otree->Branch("driftDmc",&i_driftDmc);
        otree->Branch("type",&i_type);
        otree->Branch("np",&i_np);
        otree->Branch("ip",&i_ip);
        otree->Branch("clk",&i_clk);
        otree->Branch("width",&i_width);
        otree->Branch("peak",&i_peak);
        otree->Branch("height",&i_height);
        otree->Branch("mpn",&i_mpn);
        otree->Branch("mpi",&i_mpi);
        if (has_rank) otree->Branch("rank",&i_rank);
        otree->Branch("aa",&i_aa);
        if (has_ped) otree->Branch("ped",&i_ped);
        otree->Branch("sum",&i_sum);
        otree->Branch("driftD",&o_driftD);
        otree->Branch("driftDs",&o_driftDs);
        otree->Branch("driftD0",&(i_driftD[0]));
        otree->Branch("npairs",&(npairs[0]));
        otree->Branch("isel",&(isel[0]));
        otree->Branch("icom",&(icom[0]));
        otree->Branch("islx",&(islx[0]));
        otree->Branch("islz",&(islz[0]));
        otree->Branch("iinx",&(iinx[0]));
        otree->Branch("iinz",&(iinz[0]));
        otree->Branch("chi2x",&(chi2x[0]));
        otree->Branch("chi2z",&(chi2z[0]));
        otree->Branch("chi2i",&(chi2i[0]));
        otree->Branch("calD",&(i_calD[0]));
        otree->Branch("nHitsS",&(nHitsS[0]));
        otree->Branch("slx",&(slx[0]));
        otree->Branch("slz",&(slz[0]));
        otree->Branch("inx",&(inx[0]));
        otree->Branch("inz",&(inz[0]));
        otree->Branch("chi2",&(chi2[0]));
        otree->Branch("chi2p",&(chi2p[0]));
        otree->Branch("chi2a",&(chi2a[0]));
        if (inputType){
            otree->Branch("chi2mc",&(chi2mc[0]));
            otree->Branch("chi2pmc",&(chi2pmc[0]));
            otree->Branch("chi2amc",&(chi2amc[0]));
			otree->Branch("slxmc",&slxmc);
			otree->Branch("slzmc",&slzmc);
			otree->Branch("inxmc",&inxmc);
			otree->Branch("inzmc",&inzmc);
        }
        otree->Branch("fitD",&(i_fitD[0]));
        otree->Branch("sel",&(i_sel[0]));
        o_driftD = new std::vector<double>;
        o_driftDs = new std::vector<int>;

		// Get new driftD
        for ( int iEntry = 0 ; iEntry<N; iEntry++){
            if (iEntry%10000==0) printf("%d\n",iEntry);
            if (debugLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

			// decide which candidate to use
			theCand = 0;
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
            if (xtType>=3){
                otree->SetBranchAddress("driftD",&(i_driftD[theCand]));
                otree->SetBranchAddress("npairs",&(npairs[theCand]));
                otree->SetBranchAddress("isel",&(isel[theCand]));
                otree->SetBranchAddress("icom",&(icom[theCand]));
                otree->SetBranchAddress("islx",&(islx[theCand]));
                otree->SetBranchAddress("islz",&(islz[theCand]));
                otree->SetBranchAddress("iinx",&(iinx[theCand]));
                otree->SetBranchAddress("iinz",&(iinz[theCand]));
                otree->SetBranchAddress("chi2x",&(chi2x[theCand]));
                otree->SetBranchAddress("chi2z",&(chi2z[theCand]));
                otree->SetBranchAddress("chi2i",&(chi2i[theCand]));
                otree->SetBranchAddress("calD",&(i_calD[theCand]));
                otree->SetBranchAddress("nHitsS",&(nHitsS[theCand]));
                otree->SetBranchAddress("slx",&(slx[theCand]));
                otree->SetBranchAddress("slz",&(slz[theCand]));
                otree->SetBranchAddress("inx",&(inx[theCand]));
                otree->SetBranchAddress("inz",&(inz[theCand]));
                otree->SetBranchAddress("chi2",&(chi2[theCand]));
                otree->SetBranchAddress("chi2p",&(chi2p[theCand]));
                otree->SetBranchAddress("chi2a",&(chi2a[theCand]));
                if (inputType){
                    otree->SetBranchAddress("chi2mc",&(chi2mc[theCand]));
                    otree->SetBranchAddress("chi2pmc",&(chi2pmc[theCand]));
                    otree->SetBranchAddress("chi2amc",&(chi2amc[theCand]));
                    otree->SetBranchAddress("slxmc",&slxmc);
                    otree->SetBranchAddress("slzmc",&slzmc);
                    otree->SetBranchAddress("inxmc",&inxmc);
                    otree->SetBranchAddress("inzmc",&inzmc);
                }
                otree->SetBranchAddress("fitD",&(i_fitD[theCand]));
                otree->SetBranchAddress("sel",&(i_sel[theCand]));
            }

            // set driftD
            o_driftD->clear();
            o_driftDs->clear();
            minres = 1e9;
            theDD = 1e9;
            theDT = 1e9;
            has = 0;
            theWid = -1;
			theSum = 0;
			sum1st = 0;
			dt1st = 1e9;
			thePeak = 0;
			theHeight = 0;
			theIp = 0;
			theMpi = 0;
			nSmallSumHits = 0;
			nShadowedHits = 0;
			nLateHits = 0;
            nSmallBoundaryHits = 0;
            nBoundaryHits = 0;
            for (int ihit = 0; ihit<nHits; ihit++){
            	double dt = (*i_driftT)[ihit];
            	double dd0 = (*i_driftD[theCand])[ihit];
				int ip = 0;
				for (int jhit = ihit-1; jhit>0; jhit--){
					if ((*i_layerID)[jhit]!=(*i_layerID)[ihit]||(*i_wireID)[jhit]!=(*i_wireID)[ihit]) break;
					//int type = getHitType((*i_type)[jhit],(*i_fitD[theCand])[jhit]>=0);
					//if (type<100) ip++;
					ip++; // FIXME: ignore type for this moment
				}
            	if ((*i_sel[theCand])[ihit]==1){
            		if((fabs(dd0)<0.5||fabs(dd0)>7.5))
						nBoundaryHits++;
            		if((fabs(dd0)<0.25||fabs(dd0)>7.75))
						nSmallBoundaryHits++;
					if(ip!=0){
						nLateHits++;
						for (int jhit = ihit-1; jhit>0; jhit--){
							if ((*i_ip)[jhit]==0){
								if (dt1st>(*i_driftT)[jhit]){
									sum1st = (*i_sum)[jhit];
									dt1st = (*i_driftT)[jhit];
								}
								break;
							}
						}
					}
					if((*i_mpi)[ihit]!=0)
						nShadowedHits++;
					if(has_rank&&(*i_rank)[ihit]!=0)
						nSmallSumHits++;
				}
            	double dd;
            	int status = fXTAnalyzer->t2d(dt,dd,dd0>0);
            	o_driftD->push_back(dd);
            	o_driftDs->push_back(status);
                int tlid = (*i_layerID)[ihit];
                int twid = (*i_wireID)[ihit];
				(*i_fitD[theCand])[ihit]-=off[tlid][twid];
                double tfitD = (*i_fitD[theCand])[ihit];
				int ttype = getHitType((*i_type)[ihit],tfitD>=0);
                //if (tlid==lid&&ttype<100&&status==0&&fabs(tfitD-dd)<fabs(minres)){
                if (tlid==lid&&status==0&&fabs(tfitD-dd)<fabs(minres)){ // no cut for test layer
                    minres = tfitD-dd;
                    theDD = dd;
                    theDT = dt;
                    theWid = (*i_wireID)[ihit];
                    theSum = (*i_sum)[ihit];
                    double ped = 220;
                    if (has_ped) ped = (*i_ped)[ihit];
                    thePeak = (*i_peak)[ihit]-ped;
                    theHeight = (*i_height)[ihit]-ped;
                    theIp = ip;
                    theMpi = (*i_mpi)[ihit];
                    has = 1;
                }
			}
			otree->Fill();
		}
		otree->Write();
        ofile->Close();

        if (debugLevel>=20) printf("Finished!\n");
	}

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
    fprintf(stderr,"%s [runNo] [prerunname] [runname] <[xtType: 3, sym with min nLHits, 2, sym, thr 0; 1, sym+offset; 0, no req] [geoSetup: 0, normal;1, finger] [saveHists: 0;1] [inputType: 0, Real data; 1, MC] [debug: 0;...]>\n",name);
}
