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
    TString m_preRunName = "pre";
    TString m_runName = "cur";
    int m_iEntryStart = -1;
    int m_iEntryStop = -1;
    int m_nEntries = 0;
    int m_modulo = 100;
    bool m_memdebug = false;
    int m_testLayer = 4;
    int m_resultsToSave = 4;
    bool m_createTrivalBranches = false;
    TString m_wireAdjustmentFile = "";

    // Load options
    int    opt_result;
    while((opt_result=getopt(argc,argv,"A:B:C:D:E:L:MN:P:R:S:TV:h"))!=-1){
        switch(opt_result){
            case 'M':
                m_memdebug = true;
                Log::ConfigureD("Memory=Debug");
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
                m_testLayer = atoi(optarg);
                printf("Test layer set to %d\n",m_testLayer);
                break;
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'S':
                m_resultsToSave = atoi(optarg);
                printf("Save up to %d fitting results\n",m_resultsToSave);
                break;
            case 'T':
                m_createTrivalBranches = true;
                printf("Create trivial branches\n");
                break;
            case 'A':
                m_wireAdjustmentFile = optarg;
                printf("Using wire adjustment file \"%s\"\n",optarg);
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
    m_preRunName = argv[optind++];
    m_runName= argv[optind++];

    printf("##############%s##################\n",argv[0]);
    printf("runNo       = %d\n",m_runNo);
    printf("preRunName  = \"%s\"\n",m_preRunName.Data());
    printf("runName     = \"%s\"\n",m_runName.Data());
    printf("test layer  = %d\n",m_testLayer);
    printf("Start Entry = %d\n",m_iEntryStart);
    printf("Stop Entry  = %d\n",m_iEntryStop);
    ParameterManager::Get().Print();

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    // Prepare managers
    bool success = false;
    success = RunInfoManager::Get().Initialize(m_runNo,m_preRunName,m_runName,m_testLayer);RunInfoManager::Get().Print();
    if (!success) {MyError("Cannot initialize RunInfoManager"); return 1;}
    success = BeamManager::Get().Initialize(ParameterManager::Get().beamType);BeamManager::Get().Print();
    if (!success) {MyError("Cannot initialize BeamManager"); return 1;}
    success = GeometryManager::Get().Initialize(ParameterManager::Get().geoSetup,ParameterManager::Get().connectionType,ParameterManager::Get().chamberType); GeometryManager::Get().Print();
    if (!success) {MyError("Cannot initialize GeometryManager"); return 1;}
    if (m_wireAdjustmentFile=="") m_wireAdjustmentFile = Form("%s/info/offset.%d.%s.root",HOME.Data(),m_runNo,m_preRunName.Data());
    success = GeometryManager::Get().AdjustWirePosition(m_wireAdjustmentFile);
    if (!success) MyWarn("Cannot load offset file for wire adjustment. Will ignore this step.");
    success = XTManager::Get().Initialize();
    if (!success) {MyError("Cannot initialize XTManager"); return 1;}
    XTManager::Get().Print();
    InputOutputManager::Get().readHitFile = true;
    InputOutputManager::Get().writeTrackFile = true;
    success = InputOutputManager::Get().Initialize(m_createTrivalBranches);
    if (!success) {MyError("Cannot initialize InputOutputManager"); return 1;}

    // for track finding
    double sciYup = GeometryManager::Get().GetScintillator()->Yup;
    double sciYdown = GeometryManager::Get().GetScintillator()->Ydown;
    int lidStart = ParameterManager::Get().TrackingParameters.lidStart;
    int lidStop = ParameterManager::Get().TrackingParameters.lidStop;
    int blindLayer = ParameterManager::Get().TrackingParameters.BlindLayer;
    int nPairsMin = ParameterManager::Get().TrackingParameters.nPairsMin;
    int nHitsGMax = ParameterManager::Get().TrackingParameters.nHitsGMax;
    int nHitsSMin = ParameterManager::Get().TrackingParameters.nHitsSMin;
    double sumCut = ParameterManager::Get().TrackingParameters.sumCut;
    double aaCut = ParameterManager::Get().TrackingParameters.aaCut;
    double tmin = ParameterManager::Get().TrackingParameters.tmin;
    double tmax = ParameterManager::Get().TrackingParameters.tmax;
    TrackingPara::PeakType peakType = ParameterManager::Get().peakType;
    InputOutputManager::InputHitType inputHitType = ParameterManager::Get().inputHitType;

    // Prepare Tracker
    Tracker * tracker = new Tracker(inputHitType);
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
        int nHitsG = 0; // number of good hits in this event
        int iPeakOverCut = 0;
        for (int iHit = 0; iHit<InputOutputManager::Get().nHits; iHit++){
            double aa = InputOutputManager::Get().ADCsumAll->at(iHit);
            double sum = InputOutputManager::Get().ADCsumPacket->at(iHit);
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
            if (aa<aaCut) continue;
            if (sum<sumCut) continue;
            iPeakOverCut++; // only count peaks over size cut
            if (driftT<tmin||driftT>tmax) continue;
            tracker->hitLayerIndexMap->at(lid)->push_back(iHit);
            nHitsG++;
        }
        InputOutputManager::Get().nHitsG = nHitsG;
        /// 2. Loop in layers, see how many pairs we can get and create a list of layers to form pairs
        ///    In this way we are picking one hit per layer to form pairs and get initial track parameters.
        ///    This tracking scheme doesn't support the tracking of multiple co-existing tracks.
        int nPairs = 0;
        for (int lid = lidStart; lid <= lidStop; lid++){
            if (tracker->hitLayerIndexMap->at(lid)->size()==0) continue;
            if(lid+1<=lidStop && tracker->hitLayerIndexMap->at(lid+1)->size()>0){
                tracker->pairableLayers->push_back(lid);
                nPairs++;
            }
            else if (lid-1>=lidStart && tracker->hitLayerIndexMap->at(lid-1)->size()>0){
                tracker->pairableLayers->push_back(lid);
            }
        }
        if (Log::GetLogLevel("Tracking")>=Log::VerboseLevel) InputOutputManager::Get().Print("h"); // print hit level information
        int nCombinations = 1;
        int nPairableLayers = tracker->pairableLayers->size();
        for (int ipick = 0; ipick<nPairableLayers; ipick++){
            int lid = tracker->pairableLayers->at(ipick);
            int nhits = tracker->hitLayerIndexMap->at(lid)->size();
            MyNamedInfo("Tracking",Form("  pairable layer %d: %d good hits",lid,nhits));
            for (int i = 0; i<nhits; i++){
                int iHit = tracker->hitLayerIndexMap->at(lid)->at(i);
                int wid = InputOutputManager::Get().CellID->at(iHit);
                MyNamedInfo("Tracking","    "<<iHit<<" ["<<lid<<","<<wid<<"]");
            }
            nCombinations*=nhits;
        }
        MyNamedInfo("Tracking",Form("  => %d pairs from %d good hits in %d pairable layers with %d combinations X 2^%d L/R choices",nPairs,nHitsG,nPairableLayers,nCombinations,nPairableLayers));

        /// 3. Apply tracking after cuts
        if (nHitsG<=nHitsGMax&&nPairs>=nPairsMin){
            N_found++; // found a track
            tracker->DoTracking();

            /// 4. Check and save tracking results.
            int nHitsS = tracker->trackResults[0].hitIndexSelected.size();
            MyNamedVerbose("Tracking","Good event, after tracking, "<<nHitsS<<" hits selected in the first candidate");
            if (nHitsS>=nHitsSMin){
                N_good++; // successfully reconstructed a track
            }
        }

        InputOutputManager::Get().nCandidatesFound = tracker->nGoodTracks;
        for (int iFound = 0; iFound<tracker->nGoodTracks&&iFound<NCAND; iFound++){
            InputOutputManager::Get().SetTrack(iFound,&(tracker->trackResults[iFound]));
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

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        ParameterManager::Get().ReadInputFile(configureFile,"",false,false);
        ParameterManager::Get().LoadParameters(ParameterManager::kTracking);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTManager);
    }
}

void print_usage(char* prog_name)
{
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
    fprintf(stderr,"\t -A <file>\n");
    fprintf(stderr,"\t\t Wire adjustment file set to file\n");
    fprintf(stderr,"\t -S <nResults>\n");
    fprintf(stderr,"\t\t Save up to <nResults> fitting results\n");
    fprintf(stderr,"\t -T\n");
    fprintf(stderr,"\t\t Create trivial branches in the output file\n");
}
