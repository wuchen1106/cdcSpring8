#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"

#include "XTAnalyzer.hxx"
#include "MyProcessManager.hxx"
#include "MyRuntimeParameters.hxx"
#include "Log.hxx"
#include "header.hxx"

void print_usage(char * prog_name);

MyProcessManager * pMyProcessManager;

int main(int argc, char** argv){

    TString m_runname = "currun";
    int m_iEntryStart = 0;
    int m_iEntryStop = 0;
    int m_gasID = 0; // 0, C2H6; 1, C4H10; 2, CH4
    int m_saveHists = 0;
    int m_verboseLevel = 0;
    int m_modulo = 10000;
    bool m_memdebug = false;
    TString m_configureFile = "";

    std::map<std::string, Log::ErrorPriority> namedDebugLevel;
    std::map<std::string, Log::LogPriority> namedLogLevel;
    int    opt_result;
	while((opt_result=getopt(argc,argv,"M:B:E:H:G:C:D:V:"))!=-1){
		switch(opt_result){
			/* INPUTS */
			case 'M':
			    m_modulo = atoi(optarg);
                printf("Printing modulo set to %d\n",m_modulo);
			case 'B':
			    m_iEntryStart = atoi(optarg);
                printf("Starting entry index set to %d\n",m_iEntryStart);
			case 'E':
			    m_iEntryStop = atoi(optarg);
                printf("Stopping entry index set to %d\n",m_iEntryStop);
			case 'H':
			    m_saveHists = atoi(optarg);
                printf("Histogram saving level set to %d\n",m_saveHists);
			case 'G':
			    m_gasID = atoi(optarg);
                printf("Gas ID set to %d\n",m_gasID);
			case 'C':
				m_configureFile = optarg;
                printf("Using configure file \"%s\"\n",optarg);
            case 'D':
                {
                    // Set the debug level for a named trace.
                    std::string arg(optarg);
                    std::size_t sep = arg.find("=");
                    if (sep != std::string::npos) {
                        std::string name = arg.substr(0,sep);
                        if (name=="Memory"||name=="memory") m_memdebug = true;
                        std::string levelName = arg.substr(sep+1);
                        switch (levelName[0]) {
                            case 'e': case 'E':
                                namedDebugLevel[name.c_str()] = Log::ErrorLevel;
                                break;
                            case 's': case 'S':
                                namedDebugLevel[name.c_str()] = Log::SevereLevel;
                                break;
                            case 'w': case 'W':
                                namedDebugLevel[name.c_str()] = Log::WarnLevel;
                                break;
                            case 'd': case 'D':
                                namedDebugLevel[name.c_str()] = Log::DebugLevel;
                                break;
                            case 't': case 'T':
                                namedDebugLevel[name.c_str()] = Log::TraceLevel;
                                break;
                            default:
                                print_usage(argv[0]);
                        }
                    }
                    break;
                }
            case 'V':
                {
                    // Set the debug level for a named trace.
                    std::string arg(optarg);
                    std::size_t sep = arg.find("=");
                    if (sep != std::string::npos) {
                        std::string name = arg.substr(0,sep);
                        std::string levelName = arg.substr(sep+1);
                        switch (levelName[0]) {
                            case 'q': case 'Q':
                                namedLogLevel[name.c_str()] = Log::QuietLevel;
                                break;
                            case 'l': case 'L':
                                namedLogLevel[name.c_str()] = Log::LogLevel;
                                break;
                            case 'i': case 'I':
                                namedLogLevel[name.c_str()] = Log::InfoLevel;
                                break;
                            case 'v': case 'V':
                                namedLogLevel[name.c_str()] = Log::VerboseLevel;
                                break;
                            default:
                                print_usage(argv[0]);
                        }
                    }
                    else{
                        m_verboseLevel = atoi(optarg);
                    }
                    break;
                }
			case '?':
				printf("Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
				break;
			case 'h':
			default:
				print_usage(argv[0]);
				return 1;
		}
	}
    for (std::map<std::string,Log::ErrorPriority>::iterator i 
            = namedDebugLevel.begin();
            i != namedDebugLevel.end();
            ++i) {
        Log::SetDebugLevel(i->first.c_str(), i->second);
    }

    for (std::map<std::string,Log::LogPriority>::iterator i 
            = namedLogLevel.begin();
            i != namedLogLevel.end();
            ++i) {
        Log::SetLogLevel(i->first.c_str(), i->second);
    }

    if (m_configureFile!=""){
        MyRuntimeParameters::Get().ReadParamOverrideFile(m_configureFile);
    }

    printf("##############%s with %d Parameters##################\n",argv[0],argc);
    printf("runname     = \"%s\"\n",m_runname.Data());
    printf("gasID       = %d\n",m_gasID);
    printf("save slice fittings at level %d\n",m_saveHists);
    printf("verboseLevel  = %d\n",m_verboseLevel);
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

	//=================================================Get related info========================================================
    // prepare new XT file for this run
    TFile * newXTFile = new TFile(Form("%s/info/xt.%s.root",HOME.Data(),m_runname.Data()),"RECREATE");
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
    XTAnalyzer * fXTAnalyzer = new XTAnalyzer(m_gasID,m_verboseLevel);
    //----------------------------------Set input file--------------------------------------------
    TChain * ichain = new TChain("t","t");
    ichain->Add(Form("%s/Input/garfXT.%s.root",HOME.Data(),m_runname.Data()));
    ichain->GetEntries();
    Long64_t N = ichain->GetEntries();
    if (N==0){
        fprintf(stderr,"ERROR: \"%s/Input/xt.%s.root\" is empty!\n",HOME.Data(),m_runname.Data());
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
    int statusInitialize = fXTAnalyzer->Initialize(m_runname.Data(),testLayer,preXTFile,newXTFile,newXTTree,xtType,m_saveHists, isTheLayer, saveEvenOdd, updateXT);
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

void print_usage(char * prog_name){
	fprintf(stderr,"Usage %s [options] prerunname runname\n",prog_name);
	fprintf(stderr,"[options]\n");
	fprintf(stderr,"\t -D <name>=[error,severe,warn,debug,trace]\n");
	fprintf(stderr,"\t\t Change the named debug level\n");
	fprintf(stderr,"\t -V <name>=[quiet,log,info,verbose]\n");
	fprintf(stderr,"\t\t Change the named log level\n");
	fprintf(stderr,"\t\t If equal sign is not found, set verbose level to the given value\n");
	fprintf(stderr,"\t -C <file>\n");
	fprintf(stderr,"\t\t Set the configure file\n");
	fprintf(stderr,"\t -M <n>\n");
	fprintf(stderr,"\t\t Printing modulo set to n\n");
	fprintf(stderr,"\t -B <n>\n");
	fprintf(stderr,"\t\t Starting entry index set to n\n");
	fprintf(stderr,"\t -E <n>\n");
	fprintf(stderr,"\t\t Stopping entry index set to n\n");
	fprintf(stderr,"\t -H <h>\n");
	fprintf(stderr,"\t\t Histogram saving level set to h\n");
	fprintf(stderr,"\t -G <g>\n");
	fprintf(stderr,"\t\t Gas ID set to h\n");
	fprintf(stderr,"\t\t 0, C2H6; 1, C4H10; 2, CH4\n");
}
