#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"

#include "XTAnalyzer.h"

#include "header.h"

void printUsage(char * name);

int main(int argc, char** argv){

	//=================================================Get options========================================================
	if (argc<3){
	    printUsage(argv[0]);
		return 1;
	}
    TString runname  = argv[1];
    int gasID = (int)strtol(argv[2],NULL,10);
    int saveHists = 0;
    if (argc>=4)
        saveHists = (int)strtol(argv[3],NULL,10);
    int debugLevel = 0;
    if (argc>=5)
        debugLevel = (int)strtol(argv[4],NULL,10);
    printf("##############Input %d Parameters##################\n",argc);
    printf("runname     = \"%s\"\n",runname.Data());
    printf("gasID       = %d\n",gasID);
    printf("save slice fittings? \"%s\"\n",saveHists?"yes":"no");
    printf("debugLevel  = %d\n",debugLevel);
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

	//=================================================Get related info========================================================
    // prepare new XT file for this run
    TFile * newXTFile = new TFile(Form("%s/info/xt.%s.root",HOME.Data(),runname.Data()),"RECREATE");
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

	//==============================================Prepare input file & output variables=================================================
    // input file
    double i_driftTime = 0;
    double i_doca = 0;

    //=================================================Start to get XT====================================================
    // Prepare XTAnalyzer
    XTAnalyzer * fXTAnalyzer = new XTAnalyzer(gasID,debugLevel);
    //----------------------------------Set input file--------------------------------------------
    TChain * ichain = new TChain("t","t");
    ichain->Add(Form("%s/Input/garfXT.%s.root",HOME.Data(),runname.Data()));
    ichain->GetEntries();
    Long64_t N = ichain->GetEntries();
    if (N==0){
        fprintf(stderr,"ERROR: \"%s/Input/xt.%s.root\" is empty!\n",HOME.Data(),runname.Data());
        return -1;
    }
    ichain->SetBranchAddress("driftTime",&i_driftTime);
    ichain->SetBranchAddress("dca",&i_doca);

    //----------------------------------Initialize the analyzer--------------------------------------------
    int saveEvenOdd = 0;
    TFile * preXTFile = 0;
    int xtType = -1; // -1 stands for garfield simulation
    int testLayer=4;
    bool updateXT = true;
    bool isTheLayer = true;
    int statusInitialize = fXTAnalyzer->Initialize(runname.Data(),testLayer,preXTFile,newXTFile,newXTTree,xtType,saveHists, isTheLayer, saveEvenOdd, updateXT);
    if (statusInitialize){
        fprintf(stderr,"ERROR: something wrong with initializing XTAnalyzer!\n");
        return -1;
    }

    //----------------------------------Loop in events--------------------------------------------
    for ( int iEntry = 0; iEntry<N; iEntry++){
        ichain->GetEntry(iEntry);
        fXTAnalyzer->Push(i_driftTime,i_doca*10); // cm -> mm
    }
    // fit histograms/graphs, make plots, and save new xt file
    fXTAnalyzer->Process();

    return 0;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runname] [gasID: 0, C2H6; 1, C4H10;] <[saveHists: (0);1] [debugLevel: (0)]>\n",name);
}
