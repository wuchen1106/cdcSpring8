#include <vector>
#include <map>
#include <math.h>
#include <unistd.h> /* getopt */
#include <stdlib.h> /* atoi, atof */

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "MyProcessManager.hxx"
#include "Log.hxx"

#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"
#include "BeamManager.hxx"
#include "GeometryManager.hxx"
#include "XTManager.hxx"
#include "InputOutputManager.hxx"

#include "XTAnalyzer.hxx"

#define NBINS    20
#define MAXTRUNC 6

#define NWIRES4CALIB 4
#define FIRSTWIRE4CALIB 3

// for selecting good hits in tracking
double t_min = 0;
double t_max = 0;
double sumCut = 0;
double aaCut = 0;
// definition of golden hit
double gold_t_min = 0;
double gold_t_max = 0;

void print_usage(char * prog_name);
bool isGoodHit(int iHit);
bool isGoldenHit(int iHit);
int CountGoodHitBeforeIt(int iHit);
int CountNotGoldenHitSelected(int iCand);
int GetCandidate(TString & candSelBy);
void getRunTimeParameters(TString configureFile);

MyProcessManager * pMyProcessManager;

int main(int argc, char** argv){
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int m_runNo = 0;
    TString m_preRunName = "pre";
    TString m_runName = "cur";
    int m_iEntryStart = -1;
    int m_iEntryStop = -1;
    int m_nEntries = 0;
    int m_modulo = 100;
    bool m_memdebug = false;
    int m_defaultLayerID = 4;
    bool m_DrawDetails = false;
    bool m_SaveHists = false;

    // Load options
    int    opt_result;
    while((opt_result=getopt(argc,argv,"B:C:D:E:HL:MN:P:R:S:V:h"))!=-1){
        switch(opt_result){
            /* INPUTS */
            case 'M':
                m_memdebug = true;
                printf("Turning on memory debug\n");
            case 'P':
                m_modulo = atoi(optarg);
                printf("Printing modulo set to %d\n",m_modulo);
                break;
            case 'R':
                m_runNo = atoi(optarg);
                printf("Run number set to %d\n",m_runNo);
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
            case 'L':
                m_defaultLayerID = atoi(optarg);
                printf("Test layer set to %d\n",m_defaultLayerID);
                break;
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'S':
                m_SaveHists = atoi(optarg);
                printf("Histogram saving level set to %d\n",m_SaveHists);
                break;
            case 'H':
                m_DrawDetails = atoi(optarg);
                printf("Draw bin-by-bin histogram \n");
                break;
            case 'D':
                {
                    if (!Log::ConfigureD(optarg)) print_usage(argv[0]);
                    break;
                }
            case 'V':
                {
                    if (!Log::ConfigureV(optarg)) print_usage(argv[0]);
                    break;
                }
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
    m_preRunName = argv[optind++];
    m_runName= argv[optind++];

    printf("##############%s##################\n",argv[0]);
    printf("runNo               = %d\n",m_runNo);
    printf("preRunName          = \"%s\"\n",m_preRunName.Data());
    printf("runName             = \"%s\"\n",m_runName.Data());
    printf("default test layer  = %d\n",m_defaultLayerID);
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
    success = RunInfoManager::Get().Initialize(m_runNo,m_preRunName,m_runName,m_defaultLayerID);RunInfoManager::Get().Print(); // note that the default layerID here is not important. We will set test layerID in the loop below
    if (!success) {MyError("Cannot initialize RunInfoManager"); return 1;}
    success = BeamManager::Get().Initialize(ParameterManager::Get().beamType);BeamManager::Get().Print();
    if (!success) {MyError("Cannot initialize BeamManager"); return 1;}
    success = GeometryManager::Get().Initialize(ParameterManager::Get().geoSetup,ParameterManager::Get().connectionType,ParameterManager::Get().chamberType); GeometryManager::Get().Print();
    if (!success) {MyError("Cannot initialize GeometryManager"); return 1;}
    success = GeometryManager::Get().AdjustWirePosition(Form("%s/info/offset.%d.%s.root",HOME.Data(),m_runNo,m_preRunName.Data()));
    if (!success) MyWarn("Cannot load offset file for wire adjustment. Will ignore this step.");
    success = XTManager::Get().Initialize();
    if (!success) {MyError("Cannot initialize XTManager"); return 1;}
    InputOutputManager::Get().readHitFile = true;
    InputOutputManager::Get().readTrackFile = true;

    // for getting XT
    sumCut = ParameterManager::Get().TrackingParameters.sumCut;
    aaCut = ParameterManager::Get().TrackingParameters.aaCut;
    t_min = ParameterManager::Get().TrackingParameters.tmin;
    t_max = ParameterManager::Get().TrackingParameters.tmax;
    gold_t_min = ParameterManager::Get().XTAnalyzerParameters.gold_t_min;
    gold_t_max = ParameterManager::Get().XTAnalyzerParameters.gold_t_max;
    TString CandSelBy = ParameterManager::Get().XTAnalyzerParameters.CandSelBy;
    bool RequireInTriggerCounter = ParameterManager::Get().XTAnalyzerParameters.RequireInTriggerCounter;
    bool RequireAllGoldenHits = ParameterManager::Get().XTAnalyzerParameters.RequireAllGoldenHits;
    bool ClosestPeak = ParameterManager::Get().XTAnalyzerParameters.ClosestPeak;
    bool UseGoodHit = ParameterManager::Get().XTAnalyzerParameters.UseGoodHit;
    bool AllGoodHitsUsed = ParameterManager::Get().XTAnalyzerParameters.AllGoodHitsUsed;
    int nHits_max = ParameterManager::Get().XTAnalyzerParameters.nHits_max;
    int nHitsS_min = ParameterManager::Get().XTAnalyzerParameters.nHitsS_min;
    double chi2_max = ParameterManager::Get().XTAnalyzerParameters.chi2_max;
    double slz_min = ParameterManager::Get().XTAnalyzerParameters.slz_min;
    double slz_max = ParameterManager::Get().XTAnalyzerParameters.slz_max;

    // Get the previous XT file used to do this tracking
    // NOTE: this will be used by XTAnalyzer to compare with the new XT as an iteration feedback
    TFile * preXTFile = 0;
    preXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_preRunName.Data()));
    if (!preXTFile||preXTFile->IsZombie()){
        MyNamedWarn("GetXT","Cannot find xt file according to the given prerunname. Will use garfield xt instead.");
        preXTFile = new TFile(HOME+Form("/info/xt.%s.%d.root",RunInfoManager::Get().gasTypeShort.Data(),RunInfoManager::Get().HV));
        if (!preXTFile||preXTFile->IsZombie()){
            MyError("Cannot find the default garfield xt: "<<HOME+Form("/info/xt.%s.%d.root",RunInfoManager::Get().gasTypeShort.Data(),RunInfoManager::Get().HV));
            return -1;
        }
    }

    // prepare XT files
    TFile * newXTFile = 0;
    TTree * newXTTree = 0;
    newXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_runName.Data()),"RECREATE");
    newXTTree = new TTree("t","t");
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
    XTAnalyzer * fXTAnalyzer = new XTAnalyzer(Form("%d.%s",m_runNo,m_runName.Data()),RunInfoManager::Get().gasID,preXTFile,newXTFile,newXTTree,m_defaultLayerID,m_SaveHists,m_DrawDetails);

    //=================================================Loop in layers to get XTs====================================================
    // Loop in layers
    for (int testLayer = 0; testLayer<NLAY; testLayer++){
        //----------------------------------Set input file--------------------------------------------
        MyNamedLog("GetXT","In Layer "<<testLayer<<": preparing input TChain");
        RunInfoManager::Get().testLayer = testLayer;
        success = InputOutputManager::Get().Initialize();
        if (!success) {MyError("Cannot initialize InputOutputManager for layer "<<testLayer); continue;}
        Long64_t N = InputOutputManager::Get().GetEntries();
        if (N==0){
            MyError("Input file for layer "<<testLayer<<" is empty!");
            continue;
        }

        //----------------------------------Start to get XT--------------------------------------------
        //Initialize the analyzer
        fXTAnalyzer->SetTestLayer(testLayer);
        int statusInitialize = fXTAnalyzer->Initialize();
        if (statusInitialize){
            fprintf(stderr,"WARNING: something wrong with initializing XTAnalyzer for layer[%d], will ignore this layer!\n",testLayer);
            continue;
        }
        MyNamedInfo("GetXT",Form("##############Layer %d: loop to get new XTs#############",testLayer));
        if (m_iEntryStop<0||m_iEntryStart<0){m_iEntryStart = 0; m_iEntryStop=N-1;}
        MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        for (Long64_t iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
            MyNamedInfo("GetXT","############ Entry "<<iEntry<<" #############");
            MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
            if (iEntry%m_modulo == 0){
                MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
                std::cout<<iEntry<<std::endl;
            }
            InputOutputManager::Get().Reset();
            InputOutputManager::Get().GetEntry(iEntry);

            // decide which candidate to use
            int theCand = GetCandidate(CandSelBy);
            double slx = InputOutputManager::Get().slopeX[theCand];
            double inx = InputOutputManager::Get().interceptX[theCand];
            double slz = InputOutputManager::Get().slopeZ[theCand];
            double inz = InputOutputManager::Get().interceptZ[theCand];

            // ignore events with bad fitting
            if (InputOutputManager::Get().nHitsS[theCand]<nHitsS_min) continue;
            if (InputOutputManager::Get().chi2[theCand]>chi2_max) continue;
            if (InputOutputManager::Get().slopeZ[theCand]<slz_min) continue;
            if (InputOutputManager::Get().slopeZ[theCand]>slz_max) continue;
            if (AllGoodHitsUsed&&InputOutputManager::Get().nHitsG>InputOutputManager::Get().nHitsS[theCand]) continue;
            if (RequireInTriggerCounter&&!GeometryManager::Get().IsInScinti(1.5,inx,slx,inz,slz)) continue;
            if (nHits_max&&InputOutputManager::Get().nHits>nHits_max) continue;
            if (RequireAllGoldenHits){
                if (CountNotGoldenHitSelected(theCand)) continue;
            }

            MyNamedVerbose("GetXT","  Good Event! Looping in "<<InputOutputManager::Get().nHits<<" hits");
            // find the closest hit in the test layer
            double residual_min = 1e9;
            bool foundTheHitInTestLayer = false;
            double driftT, fitD;
            for (int iHit = 0; iHit<InputOutputManager::Get().nHits; iHit++){
                int layerID = InputOutputManager::Get().LayerID->at(iHit);
                if (layerID!=testLayer) continue;
                int cellID = InputOutputManager::Get().CellID->at(iHit);
                if (UseGoodHit&&!isGoodHit(iHit)) continue; // only use good hit in the test layer if required
                if (ClosestPeak&&CountGoodHitBeforeIt(iHit)) continue; // if there is a good hit before this hit in the same cell, then skip it
                double tfitD = GeometryManager::Get().GetDOCA(layerID,cellID,slx,inx,slz,inz);
                double tdriftT = InputOutputManager::Get().DriftT->at(iHit)-InputOutputManager::Get().t0Offset[theCand]; // consider the t0 offset suggested by this candidate
                int status;
                double tdriftD = XTManager::Get().t2x(tdriftT,layerID,cellID,tfitD>0?1:-1,status);
                if (fabs(tfitD-tdriftD)<fabs(residual_min)){ // Get the one with smallest residual
                    residual_min = tfitD-tdriftD;
                    fitD = tfitD;
                    driftT = tdriftT; // consider the t0 offset suggested by this candidate
                    foundTheHitInTestLayer = true;
                }
            }
            if (!foundTheHitInTestLayer) continue; // no good hits found in test layer

            MyNamedVerbose("GetXT","  Found hit! pushing to XTAnalyzer");
            // tell analyzer a new data point
            fXTAnalyzer->Push(driftT,fitD);
        }
        MyNamedVerbose("GetXT","Starting XT analysis");
        // fit histograms/graphs, make plots, and save new xt file
        fXTAnalyzer->Process();
    }

    return 0;
}

bool isGoodHit(int iHit){ // Here I neglected the cut on ipeak, layerID cause when I call this function I'm counting hits in a selected cell or a cell in test layer
    if (iHit<0||iHit>=InputOutputManager::Get().nHits) return false;
    if (InputOutputManager::Get().DriftT->at(iHit)<t_min||InputOutputManager::Get().DriftT->at(iHit)>t_max) return false;
    if (InputOutputManager::Get().ADCsumAll->at(iHit)<aaCut) return false;
    if (InputOutputManager::Get().ADCsumPacket->at(iHit)<sumCut) return false;
    return true;
}

bool isGoldenHit(int iHit){ // Here I neglected the cut on ipeak, layerID cause when I call this function I'm counting hits in a selected cell or a cell in test layer
    if (!isGoodHit(iHit)) return false;
    if (InputOutputManager::Get().DriftT->at(iHit)<gold_t_min||InputOutputManager::Get().DriftT->at(iHit)>gold_t_max) return false;
    return true;
}

int CountGoodHitBeforeIt(int iHit){
    int ip = 0;
    for (int jHit = iHit-1; jHit>0; jHit--){ // Here I'm assuming the hits in InputOutputManager is still kept in its original order in DAQ
        if (InputOutputManager::Get().LayerID->at(jHit)!=InputOutputManager::Get().LayerID->at(iHit)||InputOutputManager::Get().CellID->at(jHit)!=InputOutputManager::Get().CellID->at(iHit)) break; // stop checking when it leaves the current cell
        if (isGoodHit(jHit)) ip++; // only count peaks satisfying good hit requirement
    }
    return ip;
}

int CountLateHitSelected(int iCand){
    int nHits = InputOutputManager::Get().nHits;
    int nLateHits = 0;
    for (int lid = 0; lid<NLAY; lid++){
        int iHit = InputOutputManager::Get().hitIndexSelected[lid][iCand];
        if (iHit>=0&&iHit<nHits){
            nLateHits+=CountGoodHitBeforeIt(iHit);
        }
    }
    return nLateHits;
}

int CountNotGoldenHitSelected(int iCand){
    int nHits = InputOutputManager::Get().nHits;
    int nBadHits = 0;
    for (int lid = 0; lid<NLAY; lid++){
        int iHit = InputOutputManager::Get().hitIndexSelected[lid][iCand];
        if (iHit>=0&&iHit<nHits){
            if (!isGoldenHit(iHit)) nBadHits++;
        }
    }
    return nBadHits;
}

int GetCandidate(TString & candSelBy){
    int cand = 0;
    if (candSelBy=="Original"){
        cand = 0;
    }
    else if (candSelBy=="LeastLateHit"){
        int nLateHits_min = 1e9;
        for (int iCand = 0; iCand<InputOutputManager::Get().nCandidatesFound; iCand++){
            int nLateHits = CountLateHitSelected(iCand);
            if (nLateHits<nLateHits_min){
                nLateHits_min = nLateHits;
                cand = iCand;
            }
        }
    }
    else if (candSelBy=="FittingChi2"||candSelBy=="GlobalChi2"){
        double minchi2 = 1e9;
        int minNhitsS = 0;
        for (int iCand = 0; iCand<NCAND; iCand++){
            if (candSelBy=="GlobalChi2"){
                if ((minchi2>InputOutputManager::Get().chi2WithTestLayer[iCand]&&minNhitsS==InputOutputManager::Get().nHitsS[iCand])||minNhitsS<InputOutputManager::Get().nHitsS[iCand]){
                    cand = iCand;
                    minchi2 = InputOutputManager::Get().chi2WithTestLayer[iCand];
                    minNhitsS = InputOutputManager::Get().nHitsS[iCand];
                }
            }
            else if (candSelBy=="FittingChi2"){
                if ((minchi2>InputOutputManager::Get().chi2[iCand]&&minNhitsS==InputOutputManager::Get().nHitsS[iCand])||minNhitsS<InputOutputManager::Get().nHitsS[iCand]){
                    cand = iCand;
                    minchi2 = InputOutputManager::Get().chi2[iCand];
                    minNhitsS = InputOutputManager::Get().nHitsS[iCand];
                }
            }
        }
    }
    else{
        cand = 0;
    }
    return cand;
}

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] preRunName runName\n",prog_name);
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
    fprintf(stderr,"\t -B <n>\n");
    fprintf(stderr,"\t\t Starting entry index set to n\n");
    fprintf(stderr,"\t -E <n>\n");
    fprintf(stderr,"\t\t Stopping entry index set to n\n");
    fprintf(stderr,"\t -N <n>\n");
    fprintf(stderr,"\t\t Maximum number of entries set to n\n");
    fprintf(stderr,"\t -L <l>\n");
    fprintf(stderr,"\t\t Test layer set to l\n");
    fprintf(stderr,"\t -H \n");
    fprintf(stderr,"\t\t Draw bin-by-bin histograms\n");
    fprintf(stderr,"\t -S <h>\n");
    fprintf(stderr,"\t\t Histogram saving level set to h\n");
}

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        ParameterManager::Get().ReadInputFile(configureFile,"",false,false);
        ParameterManager::Get().LoadParameters(ParameterManager::kTracking);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTAnalyzer);
    }
}

