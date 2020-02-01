#include <vector>
#include <map>
#include <math.h>
#include <unistd.h> /* getopt */
#include <stdlib.h> /* atoi, atof */

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include "MyProcessManager.hxx"
#include "Log.hxx"

#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"
#include "BeamManager.hxx"
#include "GeometryManager.hxx"
#include "XTManager.hxx"
#include "InputOutputManager.hxx"

#include "XTAnalyzer.hxx"

// for selecting good hits in tracking
double t_min = 0;
double t_max = 0;
double sumCut = 0;
double aaCut = 0;
// definition of golden hit
double gold_t_min = 0;
double gold_t_max = 0;

void print_usage(char * prog_name);
void getRunTimeParameters(TString configureFile);

MyProcessManager * pMyProcessManager;

enum SmearType{
    kSpace,
    kSpaceUniform,
    kTime,
    kTimeUniform
} m_smearType;

int main(int argc, char** argv){
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int m_runNo = 0;
    TString m_inputXTFile = "";
    TString m_runName = "";
    int m_testLayer = 4;
    int m_iEntryStart = -1;
    int m_iEntryStop = -1;
    int m_nEntries = 0;
    int m_modulo = 10000;
    bool m_memdebug = false;
    bool m_DrawDetails = false;
    bool m_SeparateWires = false;
    double m_spatialResolution = 0.2; // by default use 200 um
    double m_timeResolution = 1; // by default use 1 ns
    TString m_wireAdjustmentFile = "";

    // Load options
    int    opt_result;
    std::size_t opt_pos;
    std::string opt_name;
    std::string opt_arg;
    std::string opt_value = "";
    while((opt_result=getopt(argc,argv,"A:B:C:D:E:HL:MN:P:R:V:Wr:h"))!=-1){
        switch(opt_result){
            /* INPUTS */
            case 'M':
                m_memdebug = true;
                printf("Turning on memory debug\n");
                break;
            case 'P':
                m_modulo = atoi(optarg);
                printf("Printing modulo set to %d\n",m_modulo);
                break;
            case 'R':
                m_runNo = atoi(optarg);
                printf("Run number set to %d\n",m_runNo);
                break;
            case 'L':
                m_testLayer = atoi(optarg);
                printf("Test layer set to %d\n",m_testLayer);
                break;
            case 'B':
                m_iEntryStart = atoi(optarg);
                printf("Starting entry index set to %d\n",m_iEntryStart);
                break;
            case 'E':
                m_iEntryStop = atoi(optarg);
                printf("Stopping entry index set to %d\n",m_iEntryStop);
                break;
            case 'N':
                m_nEntries = atoi(optarg);
                printf("Number of entries set to %d\n",m_nEntries);
                break;
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'H':
                m_DrawDetails = true;
                printf("Draw bin-by-bin histogram \n");
                break;
            case 'A':
                m_wireAdjustmentFile = optarg;
                printf("Using wire adjustment file \"%s\"\n",optarg);
                break;
            case 'r':
                opt_arg = optarg;
                opt_name = opt_arg;
                opt_pos = opt_arg.find("=");
                if (opt_pos != std::string::npos){
                    opt_name = opt_arg.substr(0,opt_pos);
                    opt_value = opt_arg.substr(opt_pos+1);
                }
                if (opt_name=="space"){
                    m_smearType = kSpace;
                    if (opt_value!=""){
                        m_smearType = kSpaceUniform;
                        m_spatialResolution = atof(opt_value.c_str());
                    }
                }
                else{
                    m_smearType = kTime;
                    if (opt_value!=""){
                        m_smearType = kTimeUniform;
                        m_timeResolution = atof(opt_value.c_str());
                    }
                }
                printf("Smearing on %s with %s\n",
                        (m_smearType==kSpace||m_smearType==kSpaceUniform?"space":"time"),
                        (m_smearType==kSpace||m_smearType==kTime?"given input XT file":(m_smearType==kSpaceUniform?Form("%.2f mm Gaussian",m_spatialResolution):Form("%.1f ns Gaussian",m_timeResolution)))
                        );
                break;
            case 'D':
                    if (!Log::ConfigureD(optarg)) print_usage(argv[0]);
                    break;
            case 'V':
                    if (!Log::ConfigureV(optarg)) print_usage(argv[0]);
                    break;
            case 'W':
                    m_SeparateWires = true;
                    printf("Will separate wires\n");
                    break;
            case '?':
                printf("Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    if (m_nEntries>0){
        if (m_iEntryStart<0) m_iEntryStart = 0;
        m_iEntryStop = m_iEntryStart+m_nEntries-1;
    }

    if (argc-optind<2){
        print_usage(argv[0]);
        return -1;
    }
    m_inputXTFile = argv[optind++];
    m_runName= argv[optind++];

    if (ParameterManager::Get().inputHitType == InputOutputManager::kData){
        MyError("Should not use Data as hit type! will change to DriftT and continue.");
        ParameterManager::Get().inputHitType = InputOutputManager::kMCDriftT;
    }
    printf("##############%s##################\n",argv[0]);
    printf("runNo               = %d\n",m_runNo);
    printf("input XT File       = \"%s\"\n",m_inputXTFile.Data());
    printf("runName             = \"%s\"\n",m_runName.Data());
    printf("Test layer ID       = %d\n",m_testLayer);
    printf("Start Entry         = %d\n",m_iEntryStart);
    printf("Stop Entry          = %d\n",m_iEntryStop);
    ParameterManager::Get().Print();

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    // Prepare managers
    bool success = false;
    success = RunInfoManager::Get().Initialize(m_runNo,m_runName,m_testLayer);RunInfoManager::Get().Print();
    if (!success) {MyError("Cannot initialize RunInfoManager"); return 1;}
    success = BeamManager::Get().Initialize(ParameterManager::Get().beamType);BeamManager::Get().Print();
    if (!success) {MyError("Cannot initialize BeamManager"); return 1;}
    success = GeometryManager::Get().Initialize(ParameterManager::Get().geoSetup,ParameterManager::Get().connectionType,ParameterManager::Get().chamberType); GeometryManager::Get().Print();
    if (!success) {MyError("Cannot initialize GeometryManager"); return 1;}
    if (m_wireAdjustmentFile!=""){
        success = GeometryManager::Get().AdjustWirePosition(m_wireAdjustmentFile);
        if (!success) MyWarn("Cannot load offset file for wire adjustment. Will ignore this step.");
    }
    if (m_inputXTFile!=""){
        success = XTManager::Get().SetInputFileXT(m_inputXTFile);
        if (!success){MyError("Invalid input XT file"); return 1;}
    }
    success = XTManager::Get().Initialize();XTManager::Get().Print();
    if (!success) {MyError("Cannot initialize XTManager"); return 1;}
    XTManager::Get().Print();
    InputOutputManager::Get().readTrackFile = true;
    InputOutputManager::Get().writeHitFile = true;
    success = InputOutputManager::Get().Initialize();
    if (!success) {MyError("Cannot initialize InputOutputManager"); return 1;}

    // for event selection
    TString CandSelBy = ParameterManager::Get().XTAnalyzerParameters.CandSelBy;
    bool RequireInTriggerCounter = ParameterManager::Get().XTAnalyzerParameters.RequireInTriggerCounter;
    bool RequireAllGoldenHits = ParameterManager::Get().XTAnalyzerParameters.RequireAllGoldenHits;
    int nHitsGMax = ParameterManager::Get().TrackingParameters.nHitsGMax;
    int nHits_max = ParameterManager::Get().XTAnalyzerParameters.nHits_max;
    int nHitsS_min = ParameterManager::Get().XTAnalyzerParameters.nHitsS_min;
    double chi2_max = ParameterManager::Get().XTAnalyzerParameters.chi2_max;
    double pValue_min = ParameterManager::Get().XTAnalyzerParameters.pValue_min;
    double slz_min = ParameterManager::Get().XTAnalyzerParameters.slz_min;
    double slz_max = ParameterManager::Get().XTAnalyzerParameters.slz_max;

    //----------------------------------Do the analysis--------------------------------------------
    // Efficiency Counters
    int N_trigger = 0;
    int N_good = 0;
    MyNamedLog("GenerateMC","  Start analyzing");
    Long64_t N = InputOutputManager::Get().GetEntries();
    int nGoodPeaks[50] = {0};
    if (m_iEntryStop<0||m_iEntryStart<0){m_iEntryStart = 0; m_iEntryStop=N-1;}
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    for (Long64_t iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
        MyNamedInfo("GenerateMC","############ Entry "<<iEntry<<" #############");
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        if (iEntry%m_modulo == 0){
            MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
            std::cout<<iEntry<<std::endl;
        }
        InputOutputManager::Get().Reset();
        InputOutputManager::Get().GetEntry(iEntry);
        N_trigger++; // triggered event

        int    nHits = InputOutputManager::Get().nHits;
        int    nHitsG = InputOutputManager::Get().nHitsG;
        bool isGoodEvent = true;
        // should have result
        if (!nHitsGMax&&InputOutputManager::Get().nCandidatesFound) isGoodEvent = false;
        // nHits cut
        if (nHits_max&&nHits>nHits_max) isGoodEvent = false;
        // decide which candidate to use
        int theCand = 0;
        double slx = InputOutputManager::Get().slopeX[theCand];
        double inx = InputOutputManager::Get().interceptX[theCand];
        double slz = InputOutputManager::Get().slopeZ[theCand];
        double inz = InputOutputManager::Get().interceptZ[theCand];
        double chi2 = InputOutputManager::Get().chi2[theCand];
        double pValue = InputOutputManager::Get().pValue[theCand];
        int    nHitsS = InputOutputManager::Get().nHitsS[theCand];
        // ignore events with bad fitting
        if (nHitsS<nHitsS_min) isGoodEvent = false;
        if (chi2_max&&chi2>chi2_max) isGoodEvent = false;
        if (pValue<pValue_min) isGoodEvent = false;
        // slope distribution?
        if (slz<slz_min) isGoodEvent = false;
        if (slz>slz_max) isGoodEvent = false;
        // extra cuts
        if (RequireInTriggerCounter&&!GeometryManager::Get().IsInScinti(1.5,inx,slx,inz,slz)) isGoodEvent = false;

        MyNamedVerbose("GenerateMC","  Good Event! Looping in "<<nHits<<" hits");

        // set the output
        InputOutputManager::Get().interceptXmc = inx;
        InputOutputManager::Get().interceptZmc = inz;
        InputOutputManager::Get().slopeXmc = slx;
        InputOutputManager::Get().slopeZmc = slz;
        if (isGoodEvent){
            for (int lid = 1; lid<=8; lid++){ // fill histograms related with signal hits
                double minDOCA = 1e9;
                int theWid = -1;
                for (int wid = 0; wid<NCEL; wid++){
                    double DOCA = GeometryManager::Get().GetDOCA(lid,wid,slx,inx,slz,inz);
                    if (fabs(DOCA)<fabs(minDOCA)){
                        minDOCA = DOCA;
                        theWid = wid;
                    }
                }
                if (theWid>=0&&fabs(minDOCA)<GeometryManager::Get().GetChamber()->cellWidth/2){
                    double driftT, driftD;
                    int status;
                    if (m_smearType==kTime||m_smearType==kTimeUniform){
                        if (m_smearType==kTime){
                            driftT = XTManager::Get().RandomDriftT(minDOCA,lid,theWid);
                        }
                        else{
                            driftT = XTManager::Get().x2t(minDOCA,lid,theWid);
                            driftT += gRandom->Gaus(0,m_timeResolution);
                        }
                        driftD = XTManager::Get().t2x(driftT,lid,theWid,minDOCA,status);
                    }
                    else if (m_smearType==kSpace||m_smearType==kSpaceUniform){
                        double err = m_spatialResolution;
                        if (m_smearType==kSpace){
                            err = XTManager::Get().GetError(minDOCA);
                        }
                        driftD = minDOCA+gRandom->Gaus(0,err);
                        driftT = XTManager::Get().x2t(minDOCA,lid,theWid);
                    }
                    InputOutputManager::Get().PushHitMC(lid,theWid,driftT,driftD,minDOCA);
                }
            }
            N_good++;
        }

        InputOutputManager::Get().Fill();
        MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    }// end of event loop
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    InputOutputManager::Get().Write();
    InputOutputManager::Get().Close();
    printf("Triggered Events: %d\n",N_trigger);
    printf("Good Events: %d\n",N_good);
    return 0;
}

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] inputXTFile runName\n",prog_name);
    fprintf(stderr,"[options]\n");
    fprintf(stderr,"\t -D <name>=[error,severe,warn,debug,trace]\n");
    fprintf(stderr,"\t\t Change the named debug level\n");
    fprintf(stderr,"\t\t if the given name is \"general\" then set the default debug level\n");
    fprintf(stderr,"\t -V <name>=[quiet,log,info,verbose]\n");
    fprintf(stderr,"\t\t Change the named log level\n");
    fprintf(stderr,"\t\t if the given name is \"general\" then set the default log level\n");
    fprintf(stderr,"\t -M\n");
    fprintf(stderr,"\t\t Turning on memory debug mode\n");
    fprintf(stderr,"\t -C <file>\n");
    fprintf(stderr,"\t\t Set the configure file\n");
    fprintf(stderr,"\t -P <n>\n");
    fprintf(stderr,"\t\t Printing modulo set to n\n");
    fprintf(stderr,"\t -R <run>\n");
    fprintf(stderr,"\t\t Run number set to run\n");
    fprintf(stderr,"\t -L <lid>\n");
    fprintf(stderr,"\t\t Test layer set to lid\n");
    fprintf(stderr,"\t -B <n>\n");
    fprintf(stderr,"\t\t Starting entry index set to n\n");
    fprintf(stderr,"\t -E <n>\n");
    fprintf(stderr,"\t\t Stopping entry index set to n\n");
    fprintf(stderr,"\t -N <n>\n");
    fprintf(stderr,"\t\t Maximum number of entries set to n\n");
    fprintf(stderr,"\t -W \n");
    fprintf(stderr,"\t\t Seperate wires\n");
    fprintf(stderr,"\t -H \n");
    fprintf(stderr,"\t\t Draw bin-by-bin histograms\n");
}

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        ParameterManager::Get().ReadInputFile(configureFile,"",false,false);
        ParameterManager::Get().LoadParameters(ParameterManager::kTracking);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTManager);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTAnalyzer);
    }
}

