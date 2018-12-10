#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"

#include "header.hxx"

#define NBIND 17
#define STEPD 0.5

double MAXFD = 6; 

void printUsage(char * name);
int doca2i(double d);
double i2doca(int i);
int getHitType(int type,bool isRight);

int main(int argc, char** argv){
    if (argc<4){
        printUsage(argv[0]);
        return -1;
    }
    int iArg = 1;
	int runNo = (int)strtol(argv[iArg],NULL,10); iArg++;
    TString originalname = argv[iArg]; iArg++;
    TString runname = argv[iArg]; iArg++;
    int isInitial=0;
    if (argc>=iArg+1)
		{isInitial = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int averageEtrack = 0;
    if (argc>=iArg+1)
		{averageEtrack = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int testLayer = 4;
    if (argc>=iArg+1)
        {testLayer = (int)strtol(argv[iArg],NULL,10);iArg++;}
	int geoSetup = 0; // 0: normal scintillator; 1: finger scintillator
    if (argc>=iArg+1)
        {geoSetup = (int)strtol(argv[iArg],NULL,10);iArg++;}
	int xtType = 2;
    if (argc>=iArg+1)
		{xtType = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int saveHists = 0;
    if (argc>=iArg+1)
        {saveHists = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int inputType = 2; // by defualt it's MC using X
    if (argc>=iArg+1)
        {inputType = (int)strtol(argv[iArg],NULL,10);iArg++;}
    double maxchi2 = 2;
    if (argc>=iArg+1)
        {maxchi2 = (double)strtod(argv[iArg],NULL);iArg++;}
    int nHitsMax = 0;
    if (argc>=iArg+1)
        {nHitsMax = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int debugLevel = 0;
    if (argc>=iArg+1)
        {debugLevel = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int iEntryStart = 0;
    if (argc>=iArg+1)
        {iEntryStart = (int)strtol(argv[iArg],NULL,10);iArg++;}
    int iEntryStop = 0;
    if (argc>=iArg+1)
        {iEntryStop = (int)strtol(argv[iArg],NULL,10);iArg++;}
    printf("##############%s with %d Parameters##################\n",argv[0],argc);
    printf("runNo        = %d\n",runNo);
    printf("originalname = \"%s\"\n",originalname.Data());
    printf("runname      = \"%s\"\n",runname.Data());
    printf("Is initial?     %s\n",isInitial?"yes":"no");
    printf("Average Etrack? %s\n",averageEtrack?"yes":"no");
    printf("geoSetup:      %s\n",geoSetup==0?"normal scintillator":"finger scintillator");
    printf("xtType:        %d\n",xtType);
    printf("save slice fittings? \"%s\"\n",saveHists?"yes":"no");
    printf("inputType    = %d, %s\n",inputType,inputType==0?"Real Data":"MC");
    printf("maxchi2      = %.3e\n",maxchi2);
    printf("maxNhits     = %d\n",nHitsMax);
    printf("debug        = %d\n",debugLevel);
    printf("Entries:      [%d~%d]\n",iEntryStart,iEntryStop);
    printf("Test layer   = %d\n",testLayer);
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

    // get XT file
    TFile * XTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),runNo,originalname.Data()));
    if (!XTFile){
        fprintf(stderr,"ERROR: Cannot find the original XT file!\n");
        return 1;
    }
    TF1 * f_l = (TF1*) XTFile->Get(Form("fl_%d",testLayer));
    TF1 * f_r = (TF1*) XTFile->Get(Form("fr_%d",testLayer));
    if (!f_l||!f_r){
        fprintf(stderr,"ERROR: Cannot find flc/frc in the original XT file!\n");
        return 1;
    }
    f_l->SetName("fl_0");
    f_r->SetName("fr_0");

    // get res file of the original run
    TFile * oriResoFile = new TFile(Form("%s/root/res_%d.%s.layer%d.root",HOME.Data(),runNo,originalname.Data(),testLayer));
    if (!oriResoFile){
        fprintf(stderr,"ERROR: Cannot the original resolution file!\n");
        return 1;
    }
    TGraph * gr_resTotx = (TGraph*) oriResoFile->Get("gxrms");
    TGraph * gr_resTotd = (TGraph*) oriResoFile->Get("gdrms");
    if (!gr_resTotd){
        fprintf(stderr,"ERROR: Cannot find error graph in the original resolution file!\n");
        return 1;
    }
    gr_resTotx->SetName("gr_resTotx");
    gr_resTotd->SetName("gr_resTotd");

    // prepare output graph
    TGraph * gr_Etrackx = new TGraph();
    gr_Etrackx->Set(NBIND);
    gr_Etrackx->SetName("gr_Etrackx");
    TGraph * gr_resTotnewx = new TGraph();
    gr_resTotnewx->Set(NBIND);
    gr_resTotnewx->SetName("gr_resTotnewx");
    TGraph * gr_resInix = new TGraph();
    gr_resInix->Set(NBIND);
    gr_resInix->SetName("gr_resInix");
    TGraph * gr_resIniOldx = new TGraph();
    gr_resIniOldx->Set(NBIND);
    gr_resIniOldx->SetName("gr_resIniOldx");
    TGraph * gr_Etrackd = new TGraph();
    gr_Etrackd->Set(NBIND);
    gr_Etrackd->SetName("gr_Etrackd");
    TGraph * gr_resTotnewd = new TGraph();
    gr_resTotnewd->Set(NBIND);
    gr_resTotnewd->SetName("gr_resTotnewd");
    TGraph * gr_resInid = new TGraph();
    gr_resInid->Set(NBIND);
    gr_resInid->SetName("gr_resInid");
    TGraph * gr_resIniOldd = new TGraph();
    gr_resIniOldd->Set(NBIND);
    gr_resIniOldd->SetName("gr_resIniOldd");

    if (!isInitial){ // not the initial step. Do the data analyzing
        // get the input file from tracking
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
        TChain * ichain = new TChain("t","t");
        ichain->Add(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),testLayer));
        ichain->GetEntries();
        Long64_t nEntries = ichain->GetEntries();
        if (nEntries==0){
            //        fprintf(stderr,"WARNING: \"%s/root/t_%d.%s.layer%d.root\" is empty! Will ignore this layer.\n",HOME.Data(),runNo,runname.Data(),testLayer);
            fprintf(stderr,"ERROR: \"%s/root/t_%d.%s.layer%d.root\" is empty!\n",HOME.Data(),runNo,runname.Data(),testLayer);
            return 0;
        }
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

        //prepare for output ROOT file
        std::vector<double> * o_driftD = 0;
        std::vector<int>    * o_driftDs = 0;
        std::vector<int> * o_channelID = 0;
        std::vector<int> * o_boardID = 0;
        int nHitsSmallAll = 0;
        int nHitsSmallSASD = 0;
        int nSmallSumHits = 0;
        int nShadowedHits = 0;
        int nLateHits = 0;
        int nBoundaryHits = 0;
        int nSmallBoundaryHits = 0;
        // the closest peak to the track in the test layer
        double minres = 1e9;
        double theFD = 1e9;
        double theDD = 1e9;
        double theDDmc = 1e9;
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
        // the highest hit in this event
        int highBid = 0;
        int highCh = 0;
        int highLid = 0;
        int highWid = 0;
        int highIp = 0;
        double highSum = 0;
        double highAA = 0;
        double highDT = 0;
        TFile * ofile = new TFile(Form("%s/root/anamc_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),testLayer),"RECREATE");
        TTree * otree = new TTree("t","t");
        otree->Branch("triggerNumber",&triggerNumber);
        otree->Branch("res",&minres);
        otree->Branch("theFD",&theFD);
        otree->Branch("theDD",&theDD);
        otree->Branch("theDDmc",&theDDmc);
        otree->Branch("theDT",&theDT);
        otree->Branch("theWid",&theWid);
        //    otree->Branch("theSum",&theSum);
        //    otree->Branch("sum1st",&sum1st);
        //    otree->Branch("dt1st",&dt1st);
        otree->Branch("has",&has);
        //    otree->Branch("thePeak",&thePeak);
        //    otree->Branch("theHeight",&theHeight);
        //    otree->Branch("theIp",&theIp);
        //    otree->Branch("theMpi",&theMpi);
        //    otree->Branch("theCand",&theCand);
        //    otree->Branch("highBid",&highBid);
        //    otree->Branch("highCh",&highCh);
        //    otree->Branch("highLid",&highLid);
        //    otree->Branch("highWid",&highWid);
        //    otree->Branch("highIp",&highIp);
        //    otree->Branch("highSum",&highSum);
        //    otree->Branch("highAA",&highAA);
        //    otree->Branch("highDT",&highDT);
        //    otree->Branch("nHitsSmallSASD",&nHitsSmallSASD);
        //    otree->Branch("nHitsSmallAll",&nHitsSmallAll);
        //    otree->Branch("nSHits",&nShadowedHits);
        //    otree->Branch("nLHits",&nLateHits);
        //    otree->Branch("nSSHits",&nSmallSumHits);
        //    otree->Branch("nBHits",&nBoundaryHits);
        //    otree->Branch("nSBHits",&nSmallBoundaryHits);
        otree->Branch("nHits",&nHits);
        otree->Branch("nHitsG",&nHitsG);
        otree->Branch("layerID",&i_layerID);
        otree->Branch("wireID",&i_wireID);
        //    otree->Branch("channelID",&o_channelID);
        //    otree->Branch("boardID",&o_boardID);
        otree->Branch("driftT",&i_driftT);
        if (inputType) otree->Branch("driftDmc",&i_driftDmc);
        if (inputType==2) otree->Branch("driftD0",&(i_driftD[0]));
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
        o_channelID = new std::vector<int>;
        o_boardID = new std::vector<int>;

        // prepare histograms and function
        TF1 * fres = new TF1("fres","gaus",-1,1);
        TH1D * h_Etrackx[NBIND];
        TH1D * h_resTotx[NBIND];
        TH1D * h_resIniOldx[NBIND];
        TH1D * h_Etrackd[NBIND];
        TH1D * h_resTotd[NBIND];
        TH1D * h_resIniOldd[NBIND];
        for (int i = 0; i<NBIND; i++){
            h_Etrackx[i] = new TH1D(Form("hex%d",i),"",256,-1,1);
            h_resTotx[i] = new TH1D(Form("hrx%d",i),"",256,-1,1);
            h_resIniOldx[i] = new TH1D(Form("hrix%d",i),"",256,-1,1);
            h_Etrackd[i] = new TH1D(Form("hed%d",i),"",256,-1,1);
            h_resTotd[i] = new TH1D(Form("hrd%d",i),"",256,-1,1);
            h_resIniOldd[i] = new TH1D(Form("hrid%d",i),"",256,-1,1);
        }

        // analyze the input data
        if (!iEntryStart&&!iEntryStop){
            iEntryStart = 0;
            iEntryStop = nEntries-1;
        }
        if (debugLevel>0) {printf("Processing %d events\n",nEntries);fflush(stdout);}
        for ( int iEntry = iEntryStart ; iEntry<=iEntryStop; iEntry++){
            if (iEntry%10000==0) printf("%d\n",iEntry);
            if (debugLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

            // decide which candidate to use
            int theCand = 0;
            // initialize
            minres = 1e9;
            theWid = -1;
            theFD = 0;
            theDD = 0;
            theDDmc = 0;
            theDT = 0;
            has = 0;

            // check quality
            bool isGoodEvent = true;
            // ignore events with bad fitting
            if (nHitsS[theCand]<7) isGoodEvent = false;
            if (chi2[theCand]>maxchi2) isGoodEvent = false;
            //if (nHitsG>nHitsS[theCand]) isGoodEvent = false;
            if (geoSetup==1){
                if (fabs(inz[theCand])>24) isGoodEvent = false;
            }
            else{
                if (fabs(slz[theCand])>0.11) isGoodEvent = false;
            }
            if (nHitsMax&&nHits>nHitsMax) isGoodEvent = false;

            // find the closest hit in the test layer
            if (isGoodEvent){
                for (int ihit = 0; ihit<nHits; ihit++){
                    int tlayerID = (*i_layerID)[ihit];
                    if (tlayerID!=testLayer) continue;
                    int twireID = (*i_wireID)[ihit];
                    double tdriftT = (*i_driftT)[ihit];
                    //                double tfitD = (*i_fitD[theCand])[ihit]-off[tlayerID][twireID];
                    double tfitD = (*i_fitD[theCand])[ihit]; // FIXME: at this moment don't consider about wire offset
                    double tdriftD = (*i_driftD[theCand])[ihit];
                    double tdriftDmc = 0; if (inputType) tdriftDmc = (*i_driftDmc)[ihit];
                    int ttype = getHitType((*i_type)[ihit],tfitD>=0);
                    if ((ttype<100&&(xtType==7||xtType==6)||(xtType!=6&&xtType!=7))&&fabs(tfitD-tdriftD)<fabs(minres)){ // Should have cut for test layer! otherwise XT will not be well tuned
                        //if (fabs(tfitD-tdriftD)<fabs(minres)) // no cut for test layer!
                        minres = tfitD-tdriftD;
                        theWid = twireID;
                        theFD = tfitD;
                        theDT = tdriftT;
                        theDD = tdriftD;
                        theDDmc = tdriftDmc;
                        has++;
                    }
                }
                if (!has) isGoodEvent = false; // didn't find any good hit
            }

            // fill the tree
            otree->Fill();

            // fill this event to the corresponding histogram
            if (isGoodEvent){
                int ibin = doca2i(theFD); // pairing with xrms
                h_Etrackx[ibin]->Fill(theFD-theDDmc);
                h_resTotx[ibin]->Fill(theDD-theFD);
                h_resIniOldx[ibin]->Fill(theDD-theDDmc);
                ibin = doca2i(theDD); // kind of driftD but more true-ish, so pairing with drms
                h_Etrackd[ibin]->Fill(theFD-theDDmc);
                h_resTotd[ibin]->Fill(theDD-theFD);
                h_resIniOldd[ibin]->Fill(theDD-theDDmc);
            }
        }

        // analyze and get tracking error for each DOCA slice
        for (int i = 0; i<NBIND; i++){
            double doca = i2doca(i);

            h_Etrackx[i]->Fit("fres","QG","");
            double sigma = fres->GetParameter(2);
            gr_Etrackx->SetPoint(i,doca,sigma);

            h_resTotx[i]->Fit("fres","QG","");
            sigma = fres->GetParameter(2);
            gr_resTotnewx->SetPoint(i,doca,sigma);

            h_resIniOldx[i]->Fit("fres","QG","");
            sigma = fres->GetParameter(2);
            gr_resIniOldx->SetPoint(i,doca,sigma);

            h_Etrackd[i]->Fit("fres","QG","");
            sigma = fres->GetParameter(2);
            gr_Etrackd->SetPoint(i,doca,sigma);

            h_resTotd[i]->Fit("fres","QG","");
            sigma = fres->GetParameter(2);
            gr_resTotnewd->SetPoint(i,doca,sigma);

            h_resIniOldd[i]->Fit("fres","QG","");
            sigma = fres->GetParameter(2);
            gr_resIniOldd->SetPoint(i,doca,sigma);

            h_Etrackx[i]->Write();
            h_resTotx[i]->Write();
            h_resIniOldx[i]->Write();
            h_Etrackd[i]->Write();
            h_resTotd[i]->Write();
            h_resIniOldd[i]->Write();
        }

        // save the tree
        otree->Write();
        ofile->Close();
    }
    else{ // the initial step. Make up a tracking resolution.
        for (int i = 0; i<NBIND; i++){
            double doca = i2doca(i);
            gr_Etrackx->SetPoint(i,doca,0);
            gr_Etrackd->SetPoint(i,doca,0);
        }
    }
    double avEtrackx = 0;
    double avEtrackd = 0;
    int nCount = 0;
    for (int i = 0; i<NBIND; i++){
        double doca,Etrack;
        gr_Etrackx->GetPoint(i,doca,Etrack);
        if (doca>MAXFD) continue;
        nCount++;
        avEtrackx+=Etrack;
        gr_Etrackd->GetPoint(i,doca,Etrack);
        avEtrackd+=Etrack;
    }
    avEtrackd/=nCount;
    avEtrackx/=nCount;
    printf("avEtrackx: %.3e, avEtrackd: %.3e\n",avEtrackx,avEtrackd);

    for (int i = 0; i<NBIND; i++){
        double resTot,resTrack,doca;
        gr_resTotx->GetPoint(i,doca,resTot);
        gr_Etrackx->GetPoint(i,doca,resTrack);
        if (averageEtrack) resTrack = avEtrackx;
        double reso = resTot>resTrack?sqrt(resTot*resTot-resTrack*resTrack):0;
        gr_resInix->SetPoint(i,doca,reso);

        gr_resTotd->GetPoint(i,doca,resTot);
        gr_Etrackd->GetPoint(i,doca,resTrack);
        if (averageEtrack) resTrack = avEtrackd;
        reso = resTot>resTrack?sqrt(resTot*resTot-resTrack*resTrack):0;
        gr_resInid->SetPoint(i,doca,reso);
    }
    double avresInix = 0;
    double avresInid = 0;
    nCount = 0;
    for (int i = 0; i<NBIND; i++){
        double doca,resIni;
        gr_resInix->GetPoint(i,doca,resIni);
        if (doca>MAXFD) continue;
        nCount++;
        avresInix+=resIni;
        gr_resInid->GetPoint(i,doca,resIni);
        avresInid+=resIni;
    }
    avresInid/=nCount;
    avresInix/=nCount;
    printf("avresInix: %.3e, avresInid: %.3e\n",avresInix,avresInid);

    // save the new tracking resolution
    TFile * ofileRes = new TFile(Form("%s/info/res.%d.layer%d.%s.root",HOME.Data(),runNo,testLayer,runname.Data()),"RECREATE");
    f_r->Write();
    f_l->Write();
    // FIXME: choose ini VS d or ini VS x for smearing?
//    gr_resInid->SetName("gr_resIni");
    gr_resInix->SetName("gr_resIni");
    gr_resTotx->Write();
    gr_resTotnewx->Write();
    gr_Etrackx->Write();
    gr_resInix->Write();
    gr_resIniOldx->Write();
    gr_resTotd->Write();
    gr_resTotnewd->Write();
    gr_Etrackd->Write();
    gr_resInid->Write();
    gr_resIniOldd->Write();
    ofileRes->Close();

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

int doca2i(double d){
    int n = d/STEPD;
    if (n<0) n = fabs(n);
    if (n>=NBIND) n = NBIND-1;
    return n;
}

double i2doca(int i){
    return (i+0.5)*STEPD;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [originalname] [runname] <[isInitial (0)] [averageEtrack (0)] [testLayer (4)] [geoSetup: (0), normal;1, finger] [xtType: (2) sym; 1 sym+offset; 0 no; 6 sym+offset+first OT peak; 7 sym+offset+first OT peak+2segments] [saveHists: (0);1] [inputType: 0, Real data; 1, MC T; (2), MC X] [maxchi2 (2)] [nHitsMax (0)] [debug: 0;...] [iEntryStart (0)] [iEntryStop (0)]>\n",name);
}
