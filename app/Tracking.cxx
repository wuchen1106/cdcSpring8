#include <stdio.h>  /* printf, getenv */
#include <unistd.h> /* getopt */
#include <stdlib.h> /* atoi, atof */
#include <iostream> /* cout */
#include <vector>
#include <string>
#include <map>
#include <math.h>

#include <TF1.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include "Log.hxx"
#include "MyProcessManager.hxx"

#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"
#include "BeamManager.hxx"
#include "GeometryManager.hxx"
#include "XTManager.hxx"
#include "InputOutputManager.hxx"
#include "Tracker.hxx"

#include "Track.hxx"
#include "Hit.hxx"

/// This application takes charge of:
/// 1. Load user instruction from option and relative parameters from ParameterManager and
///    * only hit selection related ones and
///    * input output related parameters and 
///    * geometry
/// 2. Prepare geometry, and apply alignment if needed
/// 2. Loop in events - hits, select good hits for each layer and set to the tracker
/// 3. Ask the tracker to do the tracking given the selected good hits
///    * tracker doesn't write to InputOutputManager but tracker can read hit information from InputOutputManager
///    * all tracking-related parameters are loaded here through ParameterManager
/// 4. Set the track results to InputOutputManager and fill into the output file
///    * number of candidates is adjustable: 0 means to save all and without sorting. others mean to save up to that number with sorting
///    * sorting method is also adjustable

//============================================================
// Global controlers
MyProcessManager * pMyProcessManager;

//============================================================
// Functions
void getRunTimeParameters(TString configureFile);
void print_usage(char* prog_name);

//============================================================
// Main function
int main(int argc, char** argv){
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int m_runNo = 0;
    TString m_inputXTFile = "";
    TString m_runName = "";
    int m_iEntryStart = -1;
    int m_iEntryStop = -1;
    int m_nEntries = 0;
    int m_modulo = 100;
    bool m_memdebug = false;
    int m_testLayer = 4;
    int m_resultsToSave = 4;
    bool m_createTrivialBranches = false;
    TString m_suffixHitFile = "";
    TString m_wireAdjustmentFile = "";
    int m_t0OffsetRange = 0;
    Tracker::SortType m_sortType = Tracker::NDFchi2;
    double m_minChi2Input = 0;
    bool m_skipLayerAllowed = false;
    bool m_trackingWithDriftDmc = false;

    // Load options
    int    opt_result;
    while((opt_result=getopt(argc,argv,"A:B:C:D:E:H:L:MN:O:P:R:S:TV:ac:dhsw"))!=-1){
        switch(opt_result){
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
                break;
            case 'B':
                m_iEntryStart = atoi(optarg);
                break;
            case 'E':
                m_iEntryStop = atoi(optarg);
                break;
            case 'N':
                m_nEntries = atoi(optarg);
                break;
            case 'L':
                m_testLayer = atoi(optarg);
                break;
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'O':
                m_t0OffsetRange = atoi(optarg);
                break;
            case 'S':
                m_resultsToSave = atoi(optarg);
                printf("Save up to %d fitting results; 0 means to save all without sorting.\n",m_resultsToSave);
                break;
            case 'T':
                m_createTrivialBranches = true;
                printf("Create trivial branches\n");
                break;
            case 'H':
                m_suffixHitFile = optarg;
                printf("Added suffix \"%s\" to the output file\n",optarg);
                break;
            case 'A':
                m_wireAdjustmentFile = optarg;
                break;
            case 'c':
                m_minChi2Input = atof(optarg);
                break;
            case 'd':
                m_trackingWithDriftDmc = true;
                break;
            case 's':
                m_skipLayerAllowed = true;
                break;
            case 'a':
                m_sortType = Tracker::chi2a;
                break;
            case 'w':
                m_sortType = Tracker::chi2WithTestLayer;
                break;
            case 'D':
                if (!Log::ConfigureD(optarg)) print_usage(argv[0]);
                break;
            case 'V':
                if (!Log::ConfigureV(optarg)) print_usage(argv[0]);
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

    printf("##############%s##################\n",argv[0]);
    printf("runNo               = %d\n",m_runNo);
    printf("input XT File       = \"%s\"\n",m_inputXTFile.Data());
    printf("runName             = \"%s\"\n",m_runName.Data());
    printf("Test layer ID       = %d\n",m_testLayer);
    printf("Start Entry         = %d\n",m_iEntryStart);
    printf("Stop Entry          = %d\n",m_iEntryStop);
    printf("Using wire adjustment file \"%s\"\n",m_wireAdjustmentFile.Data());
    printf("Scan t0 offset within -+ %d ns range.\n",m_t0OffsetRange);
    printf("Set input chi2 cut at %.3e to save tracking time\n",m_minChi2Input);
    printf("Allow tracker to skip layer even if there is good hit in it? %s\n",m_skipLayerAllowed?"yes":"no");
    printf("Ask tracker to use driftDmc instead of driftT? %s\n",m_trackingWithDriftDmc?"yes":"no");
    printf("Ask tracker to sort result by %s\n",m_sortType==Tracker::NDFchi2?"NDF&chi2":(m_sortType==Tracker::chi2a?"chi2a":"chi2WithTestLayer"));
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
    success = XTManager::Get().SetInputFileXT(m_inputXTFile);
    if (!success){MyWarn("Invalid input XT file! will continue using the default one!");}
    success = XTManager::Get().Initialize();XTManager::Get().Print();
    if (!success) {MyError("Cannot initialize XTManager"); return 1;}
    InputOutputManager::Get().readHitFile = true;
    InputOutputManager::Get().writeTrackFile = true;
    if (m_trackingWithDriftDmc) InputOutputManager::Get().hitFileIsMC = true;
    InputOutputManager::Get().SetHitFileSuffix(m_suffixHitFile); // the output file name will be like h_100SUFFIX.root
    success = InputOutputManager::Get().Initialize(m_createTrivialBranches);
    if (!success) {MyError("Cannot initialize InputOutputManager"); return 1;}

    // for track finding
    int lidStart = ParameterManager::Get().TrackingParameters.lidStart;
    int lidStop = ParameterManager::Get().TrackingParameters.lidStop;
    int blindLayer = ParameterManager::Get().TrackingParameters.BlindLayer;
    double sumCut = ParameterManager::Get().TrackingParameters.sumCut;
    double aaCut = ParameterManager::Get().TrackingParameters.aaCut;
    double tmin = ParameterManager::Get().TrackingParameters.tmin;
    double tmax = ParameterManager::Get().TrackingParameters.tmax;
    TrackingPara::PeakType peakType = ParameterManager::Get().peakType;

    // Prepare Tracker
    Tracker * tracker = new Tracker();
    tracker->SetT0OffsetRange(m_t0OffsetRange);
    tracker->SetLayerSkipping(m_skipLayerAllowed);
    tracker->SetMinChi2Input(m_minChi2Input);
    tracker->SetSortType(m_sortType);
    tracker->SetUseDriftDmc(m_trackingWithDriftDmc);
    tracker->SetMaxResults(m_resultsToSave);

    //===================Tracking====================================
    // Efficiency Counters
    int N_trigger = 0;
    int N_found = 0;
    int N_good = 0;
    Long64_t N = InputOutputManager::Get().GetEntries();
    if (m_iEntryStop<0||m_iEntryStart<0){m_iEntryStart = 0; m_iEntryStop=N-1;}
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    for (Long64_t iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
        MyNamedInfo("Tracking","############ Entry "<<iEntry<<" #############");
        MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        if (iEntry%m_modulo == 0){
            std::cout<<iEntry<<std::endl;
        }
        InputOutputManager::Get().Reset();
        InputOutputManager::Get().GetEntry(iEntry);
        tracker->Reset();
        N_trigger++; // triggered event

        /// 1. Scan raw hits, select a list of them according to predefined cuts
        InputOutputManager::Get().nHitsG = 0; // number of good hits in this event
        int iPeakOverCut = 0;
        for (int iHit = 0; iHit<InputOutputManager::Get().nHits; iHit++){
            double driftT = InputOutputManager::Get().DriftT->at(iHit);
            int iPeak = InputOutputManager::Get().iPeakInChannel->at(iHit);
            if (iPeak == 0) iPeakOverCut = 0;
            int lid = InputOutputManager::Get().LayerID->at(iHit);
            if (lid<lidStart||lid>lidStop) continue; // skip the layers outside of the tracking range
            if (lid==blindLayer) continue; // skip the blinded layer
            if (peakType!= TrackingPara::kAllPeaks&&iPeakOverCut!=0) continue; // TODO: not considering about the highest peak option yet
            if (lid==m_testLayer){
                tracker->hitIndexInTestLayer->push_back(iHit);
                continue;
            }
            double aa = InputOutputManager::Get().ADCsumAll->at(iHit);
            double sum = InputOutputManager::Get().ADCsumPacket->at(iHit);
            if (aa!=0){ // TODO: this is a temporary way to check if the input hit file is MC data
                if (aa<aaCut) continue;
                if (sum<sumCut) continue;
            }
            iPeakOverCut++; // only count peaks over size cut
            if (driftT<tmin||driftT>tmax) continue;
            tracker->hitLayerIndexMap->at(lid)->push_back(iHit);
            InputOutputManager::Get().nHitsG++;
        }
        if (Log::GetLogLevel("Tracking")>=Log::VerboseLevel) InputOutputManager::Get().Print("h"); // print hit level information

        /// 2. Do tracking if it's good
        if (tracker->GoodForTracking()){
            if (Log::GetLogLevel("Tracking")>=Log::VerboseLevel) tracker->Print("h"); // print hit level information
            N_found++; // found a candidate event
            MyNamedVerbose("Tracking","Applicable event, start tracking");
            tracker->DoTracking();
            /// 3. Check and save tracking results.
            if (tracker->nGoodTracks>0){
                int nHitsS = tracker->track3Ds[0].hitIndexSelected.size();
                MyNamedVerbose("Tracking","Good event, after tracking, "<<nHitsS<<" hits selected in the first candidate");
                if (Log::GetLogLevel("Tracking")>=Log::VerboseLevel) tracker->Print("t"); // print track level information
                N_good++; // successfully reconstructed a track
            }
        }
        else{
            MyNamedVerbose("Tracking","Not applicable for tracking");
        }

        InputOutputManager::Get().nCandidatesFound = tracker->nGoodTracks;
        for (int iFound = 0; iFound<tracker->nGoodTracks&&iFound<NCAND; iFound++){ // NCAND is defined in InputOutputManager as the maximum capacity
            InputOutputManager::Get().SetTrack(iFound,&(tracker->track3Ds[iFound]),&(tracker->track2Ds[iFound]));
        }

        InputOutputManager::Get().Fill();
        MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    }// end of event loop
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    InputOutputManager::Get().Write();
    InputOutputManager::Get().Close();
    printf("Triggered Events: %d\n",N_trigger);
    printf("Found Events: %d\n",N_found);
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
    fprintf(stderr,"\t -A <file>\n");
    fprintf(stderr,"\t\t Wire adjustment file set to file\n");
    fprintf(stderr,"\t -O <n>\n");
    fprintf(stderr,"\t\t Set to scan t0 offset within +-n ns range\n");
    fprintf(stderr,"\t -S <nResults>\n");
    fprintf(stderr,"\t\t Save up to <nResults> fitting results; 0 means to save all without sorting.\n");
    fprintf(stderr,"\t -c <v>\n");
    fprintf(stderr,"\t\t Set input chi2 cut at <v> to save tracking time. By default it's 0 which means no cut\n");
    fprintf(stderr,"\t -s\n");
    fprintf(stderr,"\t\t Allow tracker to skip layer even if there is good hit in it\n");
    fprintf(stderr,"\t -a\n");
    fprintf(stderr,"\t\t Ask tracker to sort result by chi2a instead of NDF&chi2.\n");
    fprintf(stderr,"\t -w\n");
    fprintf(stderr,"\t\t Ask tracker to sort result by chi2WithTestLayer instead of NDF&chi2.\n");
    fprintf(stderr,"\t -d\n");
    fprintf(stderr,"\t\t Ask tracker to use driftDmc instead of driftT.\n");
    fprintf(stderr,"\t -T\n");
    fprintf(stderr,"\t\t Create trivial branches in the output file\n");
    fprintf(stderr,"\t -H <suf>\n");
    fprintf(stderr,"\t\t Add suffix to the hit file, like h_100SUFFIX.root\n");
}

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        ParameterManager::Get().ReadInputFile(configureFile,"",false,false);
        ParameterManager::Get().LoadParameters(ParameterManager::kTracking);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTManager);
    }
}
