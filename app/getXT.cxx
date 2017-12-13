#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"

#include "XTAnalyzer.h"

void printUsage(char * name);

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
    int saveHists = 0;
    if (argc>=6)
        saveHists = (int)strtol(argv[5],NULL,10);
    int debugLevel = 0;
    if (argc>=7)
        debugLevel = (int)strtol(argv[6],NULL,10);

    TString HOME=getenv("CDCS8WORKING_DIR");

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
	XTAnalyzer * fXTAnalyzer = new XTAnalyzer(debugLevel);

    // input file
    int nHits;
    int nHitsG;
    std::vector<int> *    i_layerID = 0;
    std::vector<int> *    i_wireID = 0;
    std::vector<double> * i_driftT = 0;
    std::vector<double> * i_driftD = 0;
    std::vector<int> *    i_type = 0;
    std::vector<int> *    i_np = 0;
    std::vector<int> *    i_ip = 0;
    std::vector<int> *    i_width = 0;
    std::vector<int> *    i_peak = 0;
    std::vector<int> *    i_height = 0;
    std::vector<int> *    i_mpn = 0;
    std::vector<int> *    i_mpi = 0;
    std::vector<double> * i_aa = 0;
    std::vector<double> * i_sum = 0;
    int npairs;
    int isel;
    int icom;
    double iinx;
    double islx;
    double iinz;
    double islz;
    double chi2x;
    double chi2z;
    double chi2i;
    int nHitsS;
    double inx;
    double slx;
    double inz;
    double slz;
    double chi2;
    std::vector<double> * i_fitD = 0;
    std::vector<int> * i_sel = 0;

	// output file
    std::vector<double> * o_driftD = 0;
    std::vector<int>    * o_driftDs = 0;

    // Loop in layers
	for (int lid = 1; lid<=8; lid++){
        if (debugLevel>0) printf("In Layer %d: preparing input TChain\n",lid);
		TChain * ichain = new TChain("t","t");
		ichain->Add(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),lid));
        ichain->SetBranchAddress("nHits",&nHits);
        ichain->SetBranchAddress("nHitsG",&nHitsG);
        ichain->SetBranchAddress("layerID",&i_layerID);
        ichain->SetBranchAddress("wireID",&i_wireID);
        ichain->SetBranchAddress("driftT",&i_driftT);
        ichain->SetBranchAddress("type",&i_type);
        ichain->SetBranchAddress("np",&i_np);
        ichain->SetBranchAddress("ip",&i_ip);
        ichain->SetBranchAddress("width",&i_width);
        ichain->SetBranchAddress("peak",&i_peak);
        ichain->SetBranchAddress("height",&i_height);
        ichain->SetBranchAddress("mpn",&i_mpn);
        ichain->SetBranchAddress("mpi",&i_mpi);
        ichain->SetBranchAddress("aa",&i_aa);
        ichain->SetBranchAddress("sum",&i_sum);
        ichain->SetBranchAddress("driftD0",&i_driftD);
        ichain->SetBranchAddress("npairs0",&npairs);
        ichain->SetBranchAddress("isel0",&isel);
        ichain->SetBranchAddress("icom0",&icom);
        ichain->SetBranchAddress("islx0",&islx);
        ichain->SetBranchAddress("islz0",&islz);
        ichain->SetBranchAddress("iinx0",&iinx);
        ichain->SetBranchAddress("iinz0",&iinz);
        ichain->SetBranchAddress("nHitsS0",&nHitsS);
        ichain->SetBranchAddress("slx0",&slx);
        ichain->SetBranchAddress("slz0",&slz);
        ichain->SetBranchAddress("inx0",&inx);
        ichain->SetBranchAddress("inz0",&inz);
        ichain->SetBranchAddress("chi20",&chi2);
        ichain->SetBranchAddress("fitD0",&i_fitD);
        ichain->SetBranchAddress("sel0",&i_sel);

		int statusInitialize = fXTAnalyzer->Initialize(Form("%d.%s.layer%d",runNo,runname.Data(),lid),lid,preXTFile,newXTFile,newXTTree,xtType,saveHists, lid==4); // take the XT from layer 4 as default output XT (0)
		if (statusInitialize){
			fprintf(stderr,"WARNING: something wrong with initializing XTAnalyzer for layer[%d], will ignore this layer!\n",lid);
			continue;
		}

        // Loop in events
        Long64_t N = ichain->GetEntries();
        if (N==0){
            fprintf(stderr,"WARNING: \"%s/root/t_%d.%s.layer%d.root\" is empty! Will ignore this layer.\n",HOME.Data(),runNo,runname.Data(),lid);
            continue;
        }
        if (debugLevel>0) printf("Processing %d events\n",N);
        for ( int iEntry = 0 ; iEntry<N; iEntry++){
            if (N%1000==0) printf("%d\n",N);
            if (debugLevel>=10) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

            // ignore events with bad fitting
            if (nHitsS<7) continue;
            //if (nHitsG>nHitsS) continue;
            if (chi2>1) continue;
            if (abs(inz)>24) continue;

            if (debugLevel>=11) printf("  Good Event! Looping in %d hits\n",nHits);
            // find the closest hit in the test layer
            double minres = 1e9;
            bool has = false;
            int wireID;
            double driftD, driftT, fitD;
            // FIXME: test more cut
            bool hasBadHit = false;
            for (int ihit = 0; ihit<nHits; ihit++){
                double tfitD = (*i_fitD)[ihit];
                double tdriftD = (*i_driftD)[ihit];
            	if ((*i_sel)[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                int tlayerID = (*i_layerID)[ihit];
                if (tlayerID!=lid) continue;
                int ttype = (*i_type)[ihit];
                if (ttype<=3&&fabs(tfitD-tdriftD)<fabs(minres)){
                    minres = tfitD-tdriftD;
                    wireID = (*i_wireID)[ihit];
                    fitD = (*i_fitD)[ihit];
                    driftT = (*i_driftT)[ihit];
                    has = true;
                }
            }
            if (!has) continue; // no hits found in test layer
            //if (hasBadHit) continue;

            if (debugLevel>=11) printf("  Found hit! pushing to XTAnalyzer\n");
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
        int theWid = -1;
		int has = 0;
		int nBad = 0;
		int nBadS = 0;
        otree->Branch("res",&minres);
        otree->Branch("theDD",&theDD);
        otree->Branch("theDT",&theDT);
        otree->Branch("theWid",&theWid);
        otree->Branch("has",&has);
        otree->Branch("nBad",&nBad);
        otree->Branch("nBadS",&nBadS);
        otree->Branch("nHits",&nHits);
        otree->Branch("nHitsG",&nHitsG);
        otree->Branch("layerID",&i_layerID);
        otree->Branch("wireID",&i_wireID);
        otree->Branch("driftT",&i_driftT);
        otree->Branch("type",&i_type);
        otree->Branch("np",&i_np);
        otree->Branch("ip",&i_ip);
        otree->Branch("width",&i_width);
        otree->Branch("peak",&i_peak);
        otree->Branch("height",&i_height);
        otree->Branch("mpn",&i_mpn);
        otree->Branch("mpi",&i_mpi);
        otree->Branch("aa",&i_aa);
        otree->Branch("sum",&i_sum);
        otree->Branch("driftD",&o_driftD);
        otree->Branch("driftDs",&o_driftDs);
        otree->Branch("driftD0",&i_driftD);
        otree->Branch("npairs0",&npairs);
        otree->Branch("isel0",&isel);
        otree->Branch("icom0",&icom);
        otree->Branch("islx0",&islx);
        otree->Branch("islz0",&islz);
        otree->Branch("iinx0",&iinx);
        otree->Branch("iinz0",&iinz);
        otree->Branch("nHitsS0",&nHitsS);
        otree->Branch("slx0",&slx);
        otree->Branch("slz0",&slz);
        otree->Branch("inx0",&inx);
        otree->Branch("inz0",&inz);
        otree->Branch("chi20",&chi2);
        otree->Branch("fitD0",&i_fitD);
        o_driftD = new std::vector<double>;
        o_driftDs = new std::vector<int>;

		// Get new driftD
        for ( int iEntry = 0 ; iEntry<N; iEntry++){
            if (N%1000==0) printf("%d\n",N);
            if (debugLevel>=10) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

            // set driftD
            o_driftD->clear();
            o_driftDs->clear();
            minres = 1e9;
            theDD = 1e9;
            theDT = 1e9;
            theWid = -1;
            has = 0;
            nBadS = 0;
            nBad = 0;
            for (int ihit = 0; ihit<nHits; ihit++){
            	double dt = (*i_driftT)[ihit];
            	double dd0 = (*i_driftD)[ihit];
            	if ((*i_sel)[ihit]==1){
            		if((fabs(dd0)<0.5||fabs(dd0)>7.5))
						nBad++;
            		if((fabs(dd0)<0.25||fabs(dd0)>7.75))
						nBadS++;
				}
            	double dd;
            	int status = fXTAnalyzer->t2d(dt,dd,dd0>0);
            	o_driftD->push_back(dd);
            	o_driftDs->push_back(status);
                double tfitD = (*i_fitD)[ihit];
                int ttype = (*i_type)[ihit];
                int tlid = (*i_layerID)[ihit];
                if (tlid==lid&&ttype<100&&status==0&&fabs(tfitD-dd)<fabs(minres)){
                    minres = tfitD-dd;
                    theDD = dd;
                    theDT = dt;
                    theWid = (*i_wireID)[ihit];
                    has = 1;
                }
			}
			otree->Fill();
		}
		otree->Write();
        ofile->Close();

        if (debugLevel>0) printf("Finished!\n");
	}

    return 0;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [prerunname] [runname] <[xtType: 2, sym, thr 0; 1, sym; 0, no req] [saveHists: 0;1] [debug: 0;...]>\n",name);
}
