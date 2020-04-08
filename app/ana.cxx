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
bool isGoodHit(int iHit);
bool isGoldenHit(int iHit);
int CountGoodHitBeforeIt(int iHit);
int CountNotGoldenHitSelected(int iCand);
int GetCandidate(TString candSelBy, int nHitsS = 0);
void getRunTimeParameters(TString configureFile);

MyProcessManager * pMyProcessManager;

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
    bool m_SeparateWires = false;
    TString m_suffixHitFile = "";
    TString m_wireAdjustmentFile = "";
    bool m_isMC = false;
    bool m_useSideXT = false;

    // Load options
    int    opt_result;
    std::string opt_name;
    std::size_t opt_pos;
    while((opt_result=getopt(argc,argv,"A:B:C:D:E:H:L:MN:P:R:SV:Wmh"))!=-1){
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
                break;
            case 'L':
                m_testLayer = atoi(optarg);
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
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'H':
                m_suffixHitFile = optarg;
                printf("Added suffix \"%s\" to the input hit file\n",optarg);
                break;
            case 'A':
                m_wireAdjustmentFile = optarg;
                break;
            case 'D':
                if (!Log::ConfigureD(optarg)) print_usage(argv[0]);
                break;
            case 'V':
                if (!Log::ConfigureV(optarg)) print_usage(argv[0]);
                break;
            case 'W':
                m_SeparateWires = true;
                break;
            case 'S':
                m_useSideXT = true;
                break;
            case 'm':
                m_isMC = true;
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
    printf("Will separate wires? %s\n",m_SeparateWires?"yes":"no");
    printf("Input is MC sample? %s\n",m_isMC?"yes":"no");
    printf("Use side-by-side XT? %s\n",m_useSideXT?"yes":"no");
    ParameterManager::Get().Print();

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    // Prepare managers
    bool success = false;
    success = RunInfoManager::Get().Initialize(m_runNo,m_runName,m_testLayer);RunInfoManager::Get().Print();
    if (!success) {MyNamedError("Analyze","Cannot initialize RunInfoManager"); return 1;}
    success = BeamManager::Get().Initialize(ParameterManager::Get().beamType);BeamManager::Get().Print();
    if (!success) {MyNamedError("Analyze","Cannot initialize BeamManager"); return 1;}
    success = GeometryManager::Get().Initialize(ParameterManager::Get().geoSetup,ParameterManager::Get().connectionType,ParameterManager::Get().chamberType); GeometryManager::Get().Print();
    if (!success) {MyNamedError("Analyze","Cannot initialize GeometryManager"); return 1;}
    if (m_wireAdjustmentFile!=""){
        success = GeometryManager::Get().AdjustWirePosition(m_wireAdjustmentFile);
        if (!success) MyWarn("Cannot load offset file for wire adjustment. Will ignore this step.");
    }
    success = XTManager::Get().SetInputFileXT(m_inputXTFile);
    if (!success){MyError("Invalid input XT file"); return 1;}
    if (m_useSideXT) XTManager::Get().UseSideXT();
    success = XTManager::Get().Initialize();XTManager::Get().Print();
    if (!success) {MyNamedError("Analyze","Cannot initialize XTManager"); return 1;}
    InputOutputManager::Get().readHitFile = true;
    InputOutputManager::Get().readTrackFile = true;
    InputOutputManager::Get().SetHitFileSuffix(m_suffixHitFile); // the output file name will be like h_100SUFFIX.root
    if (m_isMC) InputOutputManager::Get().hitFileIsMC = true;
    success = InputOutputManager::Get().Initialize();
    if (!success) {MyNamedError("Analyze","Cannot initialize InputOutputManager"); return 1;}

    // for hit classification
    sumCut = ParameterManager::Get().TrackingParameters.sumCut;
    aaCut = ParameterManager::Get().TrackingParameters.aaCut;
    t_min = ParameterManager::Get().TrackingParameters.tmin;
    t_max = ParameterManager::Get().TrackingParameters.tmax;
    gold_t_min = ParameterManager::Get().XTAnalyzerParameters.gold_t_min;
    gold_t_max = ParameterManager::Get().XTAnalyzerParameters.gold_t_max;
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
    // for plotting
    double draw_tmin = ParameterManager::Get().XTAnalyzerParameters.draw_tmin;
    double draw_tmax = ParameterManager::Get().XTAnalyzerParameters.draw_tmax;
    double draw_xmin = ParameterManager::Get().XTAnalyzerParameters.draw_xmin;
    double draw_xmax = ParameterManager::Get().XTAnalyzerParameters.draw_xmax;
    int    draw_tnum = (draw_tmax-draw_tmin)*0.96;

    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    //----------------------------------Prepare tree--------------------------------------------
    TFile * ofile = new TFile(Form("root/ana/ana.%d.%s.layer%d.root",m_runNo,m_runName.Data(),m_testLayer),"RECREATE");
    TTree * otree = new TTree("t","t");
    int    o_nHits;
    int    o_nHitsG;
    int    o_nHitsS;
    double o_inx;
    double o_inz;
    double o_slx;
    double o_slz;
    double o_t0Offset;
    double o_inxmc;
    double o_inzmc;
    double o_slxmc;
    double o_slzmc;
    double o_t0mc;
    double o_chi2;
    double o_chi2a;
    double o_chi2w;
    double o_pValue;
    bool   o_isGood;
    int    o_wid[NLAY];
    double o_driftD[NLAY];
    double o_driftDmc[NLAY];
    double o_driftT[NLAY];
    double o_DOCA[NLAY];
    double o_DOCAmc[NLAY];
    int    o_ipeak[NLAY];
    double o_ADCheight[NLAY];
    double o_ADCsumPkt[NLAY];
    double o_ADCsumAll[NLAY];
    bool   o_foundTestLayerHit;
    if (m_isMC){
        otree->Branch("inxmc",&o_inxmc);
        otree->Branch("inzmc",&o_inzmc);
        otree->Branch("slxmc",&o_slxmc);
        otree->Branch("slzmc",&o_slzmc);
        otree->Branch("t0mc",&o_t0mc);
    }
    for (int i = 1; i<NLAY; i++){
        otree->Branch(Form("w%d",i),&o_wid[i]);
        otree->Branch(Form("D%d",i),&o_DOCA[i]);
        otree->Branch(Form("d%d",i),&o_driftD[i]);
        otree->Branch(Form("t%d",i),&o_driftT[i]);
        if (m_isMC){
            otree->Branch(Form("dmc%d",i),&o_driftDmc[i]);
            otree->Branch(Form("Dmc%d",i),&o_DOCAmc[i]);
        }
        else{
            otree->Branch(Form("i%d",i),&o_ipeak[i]);
            otree->Branch(Form("height%d",i),&o_ADCheight[i]);
            otree->Branch(Form("sum%d",i),&o_ADCsumPkt[i]);
            otree->Branch(Form("all%d",i),&o_ADCsumAll[i]);
        }
    }
    otree->Branch(Form("nHits"),&o_nHits);
    otree->Branch(Form("nHitsG"),&o_nHitsG);
    otree->Branch(Form("nHitsS"),&o_nHitsS);
    otree->Branch(Form("inx"),&o_inx);
    otree->Branch(Form("inz"),&o_inz);
    otree->Branch(Form("slx"),&o_slx);
    otree->Branch(Form("slz"),&o_slz);
    otree->Branch(Form("t0Offset"),&o_t0Offset);
    otree->Branch(Form("chi2"),&o_chi2);
    otree->Branch(Form("chi2a"),&o_chi2a);
    otree->Branch(Form("chi2w"),&o_chi2w);
    otree->Branch(Form("pValue"),&o_pValue);
    otree->Branch(Form("isGood"),&o_isGood);
    otree->Branch(Form("found"),&o_foundTestLayerHit);
    
    //----------------------------------Prepare histograms--------------------------------------------
    // Prepare histograms for hit map
    TH2D * h_hitmap = new TH2D("h_hitmap","Hit Map",1024,-10,0,NLAY,-0.5,NLAY-0.5); h_hitmap->SetContour(100);
    h_hitmap->GetYaxis()->SetTitle("Layer ID");
    h_hitmap->GetXaxis()->SetLabelOffset(10);
    // Prepare histograms for ADC & drift time distribution
    std::vector<std::vector<int> > wiresToConsider;
    for (int i = 0; i<NLAY; i++){
        std::vector<int> wires;
        wiresToConsider.push_back(wires);
    }
    int lidTemp = 1;
    wiresToConsider[lidTemp].push_back(4);wiresToConsider[lidTemp].push_back(5);lidTemp++; // layer 1
    wiresToConsider[lidTemp].push_back(3);wiresToConsider[lidTemp].push_back(4);lidTemp++; // layer 2
    wiresToConsider[lidTemp].push_back(4);wiresToConsider[lidTemp].push_back(5);lidTemp++; // layer 3
    wiresToConsider[lidTemp].push_back(4);wiresToConsider[lidTemp].push_back(5);lidTemp++; // layer 4
    wiresToConsider[lidTemp].push_back(5);wiresToConsider[lidTemp].push_back(6);lidTemp++; // layer 5
    wiresToConsider[lidTemp].push_back(4);wiresToConsider[lidTemp].push_back(5);lidTemp++; // layer 6
    wiresToConsider[lidTemp].push_back(5);wiresToConsider[lidTemp].push_back(6);lidTemp++; // layer 7
    wiresToConsider[lidTemp].push_back(4);wiresToConsider[lidTemp].push_back(5);lidTemp++; // layer 8
    TH1D * h_DriftT_all[NLAY][NCEL];
    TH1D * h_DriftT_sig[NLAY][NCEL];
    TH1D * h_DriftT_bkg[NLAY][NCEL];
    TH1D * h_ADCsumPkt_all[NLAY][NCEL];
    TH1D * h_ADCsumPkt_sig[NLAY][NCEL];
    TH1D * h_ADCsumPkt_bkg[NLAY][NCEL];
    TH1D * h_ADCsumAll_all[NLAY][NCEL];
    TH1D * h_ADCsumAll_sig[NLAY][NCEL];
    TH1D * h_ADCsumAll_bkg[NLAY][NCEL];
    TH2D * h_AVT_sig[NLAY][NCEL];
    TH2D * h_AVT_siglate[NLAY][NCEL];
    TH2D * h_AVT_bkg[NLAY][NCEL];
    TH2D * h_SVT_sig[NLAY][NCEL];
    TH2D * h_SVT_siglate[NLAY][NCEL];
    TH2D * h_SVT_bkg[NLAY][NCEL];

    for (int lid = 0; lid<=8; lid++){
        for (int wid = 0; wid<=10; wid++){
            h_DriftT_all[lid][wid] = new TH1D(Form("h_DriftT_all_%d_%d",lid,wid),Form("DriftT of all hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax); h_DriftT_all[lid][wid]->SetLineColor(kBlack);
            h_DriftT_all[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_DriftT_all[lid][wid]->GetYaxis()->SetTitle("Count");
            h_DriftT_sig[lid][wid] = new TH1D(Form("h_DriftT_sig_%d_%d",lid,wid),Form("DriftT of signal hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax); h_DriftT_sig[lid][wid]->SetLineColor(kRed);
            h_DriftT_sig[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_DriftT_sig[lid][wid]->GetYaxis()->SetTitle("Count");
            h_DriftT_bkg[lid][wid] = new TH1D(Form("h_DriftT_bkg_%d_%d",lid,wid),Form("DriftT of noise hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax); h_DriftT_bkg[lid][wid]->SetLineColor(kBlue);
            h_DriftT_bkg[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_DriftT_bkg[lid][wid]->GetYaxis()->SetTitle("Count");
            h_ADCsumPkt_all[lid][wid] = new TH1D(Form("h_ADCsumPkt_all_%d_%d",lid,wid),Form("ADCsumPkt of all hits: layer %d wire %d",lid,wid),850,-50,800); h_ADCsumPkt_all[lid][wid]->SetLineColor(kBlack);
            h_ADCsumPkt_all[lid][wid]->GetXaxis()->SetTitle("ADC sum");
            h_ADCsumPkt_all[lid][wid]->GetYaxis()->SetTitle("Count");
            h_ADCsumPkt_sig[lid][wid] = new TH1D(Form("h_ADCsumPkt_sig_%d_%d",lid,wid),Form("ADCsumPkt of signal hits: layer %d wire %d",lid,wid),850,-50,800); h_ADCsumPkt_sig[lid][wid]->SetLineColor(kRed);
            h_ADCsumPkt_sig[lid][wid]->GetXaxis()->SetTitle("ADC sum");
            h_ADCsumPkt_sig[lid][wid]->GetYaxis()->SetTitle("Count");
            h_ADCsumPkt_bkg[lid][wid] = new TH1D(Form("h_ADCsumPkt_bkg_%d_%d",lid,wid),Form("ADCsumPkt of noise hits: layer %d wire %d",lid,wid),850,-50,800); h_ADCsumPkt_bkg[lid][wid]->SetLineColor(kBlue);
            h_ADCsumPkt_bkg[lid][wid]->GetXaxis()->SetTitle("ADC sum");
            h_ADCsumPkt_bkg[lid][wid]->GetYaxis()->SetTitle("Count");
            h_ADCsumAll_all[lid][wid] = new TH1D(Form("h_ADCsumAll_all_%d_%d",lid,wid),Form("ADCsumAll of all hits: layer %d wire %d",lid,wid),850,-50,800); h_ADCsumAll_all[lid][wid]->SetLineColor(kBlack);
            h_ADCsumAll_all[lid][wid]->GetXaxis()->SetTitle("ADC sum");
            h_ADCsumAll_all[lid][wid]->GetYaxis()->SetTitle("Count");
            h_ADCsumAll_sig[lid][wid] = new TH1D(Form("h_ADCsumAll_sig_%d_%d",lid,wid),Form("ADCsumAll of signal hits: layer %d wire %d",lid,wid),850,-50,800); h_ADCsumAll_sig[lid][wid]->SetLineColor(kRed);
            h_ADCsumAll_sig[lid][wid]->GetXaxis()->SetTitle("ADC sum");
            h_ADCsumAll_sig[lid][wid]->GetYaxis()->SetTitle("Count");
            h_ADCsumAll_bkg[lid][wid] = new TH1D(Form("h_ADCsumAll_bkg_%d_%d",lid,wid),Form("ADCsumAll of noise hits: layer %d wire %d",lid,wid),850,-50,800); h_ADCsumAll_bkg[lid][wid]->SetLineColor(kBlue);
            h_ADCsumAll_bkg[lid][wid]->GetXaxis()->SetTitle("ADC sum");
            h_ADCsumAll_bkg[lid][wid]->GetYaxis()->SetTitle("Count");
            h_AVT_sig[lid][wid] = new TH2D(Form("h_AVT_sig_%d_%d",lid,wid),Form("ADCsumAll VS DriftT of signal hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax,850,-50,800); h_AVT_sig[lid][wid]->SetContour(100);
            h_AVT_sig[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_AVT_sig[lid][wid]->GetYaxis()->SetTitle("ADC");
            h_AVT_siglate[lid][wid] = new TH2D(Form("h_AVT_siglate_%d_%d",lid,wid),Form("ADCsumAll VS DriftT of late signal hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax,850,-50,800); h_AVT_siglate[lid][wid]->SetContour(100);
            h_AVT_siglate[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_AVT_siglate[lid][wid]->GetYaxis()->SetTitle("ADC");
            h_AVT_bkg[lid][wid] = new TH2D(Form("h_AVT_bkg_%d_%d",lid,wid),Form("ADCsumAll VS DriftT of noise hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax,850,-50,800); h_AVT_sig[lid][wid]->SetContour(100);
            h_AVT_bkg[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_AVT_bkg[lid][wid]->GetYaxis()->SetTitle("ADC");
            h_SVT_sig[lid][wid] = new TH2D(Form("h_SVT_sig_%d_%d",lid,wid),Form("ADCsumPkt VS DriftT of signal hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax,850,-50,800); h_SVT_sig[lid][wid]->SetContour(100);
            h_SVT_sig[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_SVT_sig[lid][wid]->GetYaxis()->SetTitle("ADC");
            h_SVT_siglate[lid][wid] = new TH2D(Form("h_SVT_siglate_%d_%d",lid,wid),Form("ADCsumPkt VS DriftT of late signal hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax,850,-50,800); h_SVT_siglate[lid][wid]->SetContour(100);
            h_SVT_siglate[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_SVT_siglate[lid][wid]->GetYaxis()->SetTitle("ADC");
            h_SVT_bkg[lid][wid] = new TH2D(Form("h_SVT_bkg_%d_%d",lid,wid),Form("ADCsumPkt VS DriftT of noise hits: layer %d wire %d",lid,wid),draw_tnum,draw_tmin,draw_tmax,850,-50,800); h_SVT_sig[lid][wid]->SetContour(100);
            h_SVT_bkg[lid][wid]->GetXaxis()->SetTitle("Drift Time [ns]");
            h_SVT_bkg[lid][wid]->GetYaxis()->SetTitle("ADC");
        }
    }
    // Prepare histograms for nHitsG
    TH2I * h_nGoodHits = new TH2I("h_nGoodHits","Number of good hits VS ADC sum cut",50,0,50,100,0,100); h_nGoodHits->SetContour(100);
    h_nGoodHits->GetXaxis()->SetTitle("ADC sum cut");
    h_nGoodHits->GetYaxis()->SetTitle("Number of good hits");

    //----------------------------------Do the analysis--------------------------------------------
    MyNamedLog("Analyze","  Start analyzing");
    Long64_t N = InputOutputManager::Get().GetEntries();
    int nGoodPeaks[50] = {0};
    if (m_iEntryStop<0||m_iEntryStart<0){m_iEntryStart = 0; m_iEntryStop=N-1;}
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    for (Long64_t iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
        MyNamedInfo("Analyze","############ Entry "<<iEntry<<" #############");
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        if (iEntry%m_modulo == 0){
            MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
            std::cout<<iEntry<<std::endl;
        }
        InputOutputManager::Get().Reset();
        InputOutputManager::Get().GetEntry(iEntry);

        // reset the values
        int    nHits = InputOutputManager::Get().nHits;
        int    nHitsG = InputOutputManager::Get().nHitsG;
        double slx = 0;
        double inx = 0;
        double slz = 0;
        double inz = 0;
        double t0Offset = 0;
        double chi2 = 0;
        double chi2a = 0;
        double chi2w = 0;
        double pValue = 0;
        int    nHitsS = 0;
        bool isGoodEvent = true;
        bool foundTestLayerHit = false;
        for (int iLayer = 0; iLayer < NLAY; iLayer++){
            o_wid[iLayer] = -1;
            o_driftD[iLayer] = 0;
            o_driftDmc[iLayer] = 0;
            o_driftT[iLayer] = 0;
            o_DOCA[iLayer] = 0;
            o_DOCAmc[iLayer] = 0;
            o_ipeak[iLayer] = 0;
            o_ADCheight[iLayer] = 0;
            o_ADCsumPkt[iLayer] = 0;
            o_ADCsumAll[iLayer] = 0;
        }
        // should have result
        // decide which candidate to use
        int theCand = GetCandidate(ParameterManager::Get().XTAnalyzerParameters.CandSelBy);
        if (theCand<0){
            MyError("Cannot find a candidate for this event! nCandidatesFound = "<<InputOutputManager::Get().nCandidatesFound);
            MyError("This will cause a mis match of entries among different files, thus we will stop this run.");
            return -1;
        }
        MyNamedInfo("Analyze","    => theCand = "<<theCand);
        slx = InputOutputManager::Get().slopeX[theCand];
        inx = InputOutputManager::Get().interceptX[theCand];
        slz = InputOutputManager::Get().slopeZ[theCand];
        inz = InputOutputManager::Get().interceptZ[theCand];
        t0Offset = InputOutputManager::Get().t0Offset[theCand];
        chi2 = InputOutputManager::Get().chi2[theCand];
        chi2a = InputOutputManager::Get().chi2a[theCand];
        chi2w = InputOutputManager::Get().chi2WithTestLayer[theCand];
        pValue = InputOutputManager::Get().pValue[theCand];
        nHitsS = InputOutputManager::Get().nHitsS[theCand];

        // nHits cut
        if (nHits_max&&nHits>nHits_max) isGoodEvent = false;
        // ignore events with bad fitting
        if (nHitsS<nHitsS_min) isGoodEvent = false;
        if (chi2_max&&chi2>chi2_max) isGoodEvent = false;
        if (pValue<pValue_min) isGoodEvent = false;
        // slope distribution?
        if (slz<slz_min) isGoodEvent = false;
        if (slz>slz_max) isGoodEvent = false;
        // extra cuts
        if (RequireInTriggerCounter&&!GeometryManager::Get().IsInScinti(1.5,inx,slx,inz,slz)) isGoodEvent = false;
        if (RequireAllGoldenHits&&CountNotGoldenHitSelected(theCand)) isGoodEvent = false;

        // reset o_wid
        for (int lid = 0; lid<NLAY; lid++){
            o_wid[lid] = -1;
        }

        MyNamedVerbose("Analyze","  Good Event! Looping in "<<nHits<<" hits");
        for (int i = 0; i<50; i++){
            nGoodPeaks[i] = 0;
        }
        foundTestLayerHit = false;
        double residualMinimal = 1e9;
        for (int iHit = 0; iHit<nHits; iHit++){
            int layerID = InputOutputManager::Get().LayerID->at(iHit);
            if (layerID==0){continue;} // don't include the guard layer
            int wireID = InputOutputManager::Get().CellID->at(iHit);
            int iPeak = InputOutputManager::Get().iPeakInChannel->at(iHit);
            double DOCAmc = 0;
            double DriftTmc = 0;
            double DriftDmc = 0;
            if (m_isMC){
                DOCAmc = InputOutputManager::Get().DOCA->at(iHit);
                DriftTmc = InputOutputManager::Get().DriftT->at(iHit)-InputOutputManager::Get().t0mc; // consider the t0 offset suggested by this candidate
                int status;
                DriftDmc = XTManager::Get().t2x(DriftTmc,layerID,wireID,DOCAmc,status);
            }
            double DOCA = GeometryManager::Get().GetDOCA(layerID,wireID,slx,inx,slz,inz);
            double DriftT = InputOutputManager::Get().DriftT->at(iHit)-InputOutputManager::Get().t0Offset[theCand]; // consider the t0 offset suggested by this candidate
            int status;
            double DriftD = XTManager::Get().t2x(DriftT,layerID,wireID,DOCA,status);
            double ADCheight = InputOutputManager::Get().ADCheight->at(iHit);
            double ADCsumPkt = InputOutputManager::Get().ADCsumPacket->at(iHit);
            double ADCsumAll = InputOutputManager::Get().ADCsumAll->at(iHit);
            h_DriftT_all[layerID][wireID]->Fill(DriftT);
            h_ADCsumPkt_all[layerID][wireID]->Fill(ADCsumPkt);
            if (iPeak==0){
                h_ADCsumAll_all[layerID][wireID]->Fill(ADCsumAll);
            }
            for (int i = 0; i<50; i++){
                sumCut = i; // adjust the cut
                aaCut = -1e9;
                if (isGoldenHit(iHit)) nGoodPeaks[i]++;
            }
            sumCut = ParameterManager::Get().TrackingParameters.sumCut;
            aaCut = ParameterManager::Get().TrackingParameters.aaCut;

            //if (!isGoodEvent) continue; // set the values regardless to if the event is good or not
            if (layerID==m_testLayer){
                if (isGoodHit(iHit)&&fabs(DOCA)<20){ // in test layer, on track, check if it's closer
                    foundTestLayerHit = true;
                    double residual = DriftD-DOCA;
                    if (fabs(residual)<fabs(residualMinimal)){
                        residualMinimal = residual;
                        o_wid[layerID] = wireID;
                        o_DOCA[layerID] = DOCA;
                        if (m_isMC){
                            o_DOCAmc[layerID] = DOCAmc;
                            o_driftDmc[layerID] = DriftDmc;
                        }
                        o_driftD[layerID] = DriftD;
                        o_driftT[layerID] = DriftT;
                        o_ipeak[layerID] = iPeak;
                        o_ADCheight[layerID] = ADCheight;
                        o_ADCsumPkt[layerID] = ADCsumPkt;
                        o_ADCsumAll[layerID] = ADCsumAll;
                    }
                }
            }
            else{
                if (iHit == InputOutputManager::Get().hitIndexSelected[layerID][theCand]){ // this is the signal hit
                    o_wid[layerID] = wireID;
                    o_DOCA[layerID] = DOCA;
                    if (m_isMC){
                        o_DOCAmc[layerID] = DOCAmc;
                        o_driftDmc[layerID] = DriftDmc;
                    }
                    o_driftD[layerID] = DriftD;
                    o_driftT[layerID] = DriftT;
                    o_ipeak[layerID] = iPeak;
                    o_ADCheight[layerID] = ADCheight;
                    o_ADCsumPkt[layerID] = ADCsumPkt;
                    o_ADCsumAll[layerID] = ADCsumAll;
                }
                else{ // this is the noise hit
                    h_DriftT_bkg[layerID][wireID]->Fill(DriftT);
                    h_ADCsumPkt_bkg[layerID][wireID]->Fill(ADCsumPkt);
                    h_AVT_bkg[layerID][wireID]->Fill(DriftT,ADCsumAll);
                    h_SVT_bkg[layerID][wireID]->Fill(DriftT,ADCsumPkt);
                    if (iPeak==0){
                        h_ADCsumAll_bkg[layerID][wireID]->Fill(ADCsumAll);
                    }
                }
            }
        }
        for (int i = 0; i<50; i++){
            h_nGoodHits->Fill(i,nGoodPeaks[i]);
        }

        if (isGoodEvent){
            for (int lid = 1; lid<=8; lid++){ // fill histograms related with signal hits
                if (lid==m_testLayer&&(!foundTestLayerHit)){ continue; }
                int wid = o_wid[lid];
                if (wid==-1) continue;
                double cellCoordinate = -wid+0.5+o_DOCA[lid]/GeometryManager::Get().fChamber->cellWidth;
                h_hitmap->Fill(cellCoordinate,lid);
                h_DriftT_sig[lid][wid]->Fill(o_driftT[lid]);
                h_ADCsumPkt_sig[lid][wid]->Fill(o_ADCsumPkt[lid]);
                h_AVT_sig[lid][wid]->Fill(o_driftT[lid],o_ADCsumAll[lid]);
                h_SVT_sig[lid][wid]->Fill(o_driftT[lid],o_ADCsumPkt[lid]);
                if (o_ipeak[lid]>0) {
                    h_AVT_siglate[lid][wid]->Fill(o_driftT[lid],o_ADCsumAll[lid]);
                    h_SVT_siglate[lid][wid]->Fill(o_driftT[lid],o_ADCsumPkt[lid]);
                }
                h_ADCsumAll_sig[lid][wid]->Fill(o_ADCsumAll[lid]);
            }
            for (int iHit = 0; iHit<nHits; iHit++){ // add noise hits of the test layer
                int layerID = InputOutputManager::Get().LayerID->at(iHit);
                int wireID = InputOutputManager::Get().CellID->at(iHit);
                int iPeak = InputOutputManager::Get().iPeakInChannel->at(iHit);
                if (layerID!=m_testLayer||(wireID==o_wid[layerID]&&iPeak==o_ipeak[layerID])){continue;} // only test layer, only hits not picked as signal
                double DOCA = GeometryManager::Get().GetDOCA(layerID,wireID,slx,inx,slz,inz);
                double DriftT = InputOutputManager::Get().DriftT->at(iHit);
                int status;
                double DriftD = XTManager::Get().t2x(DriftT,layerID,wireID,DOCA,status);
                double ADCheight = InputOutputManager::Get().ADCheight->at(iHit);
                double ADCsumPkt = InputOutputManager::Get().ADCsumPacket->at(iHit);
                double ADCsumAll = InputOutputManager::Get().ADCsumAll->at(iHit);
                h_DriftT_bkg[layerID][wireID]->Fill(DriftT);
                h_ADCsumPkt_bkg[layerID][wireID]->Fill(ADCsumPkt);
                h_AVT_bkg[layerID][wireID]->Fill(DriftT,ADCsumAll);
                h_SVT_bkg[layerID][wireID]->Fill(DriftT,ADCsumPkt);
                if (iPeak==0){
                    h_ADCsumAll_bkg[layerID][wireID]->Fill(ADCsumAll);
                }
            }
        }

        o_isGood = isGoodEvent;
        o_nHits = nHits;
        o_nHitsG = nHitsG;
        o_nHitsS = nHitsS;
        o_inx = inx;
        o_inz = inz;
        o_slx = slx;
        o_slz = slz;
        o_t0Offset = t0Offset;
        o_chi2 = chi2;
        o_chi2a = chi2a;
        o_chi2w = chi2w;
        o_pValue = pValue;
        o_foundTestLayerHit = foundTestLayerHit;

        if (m_isMC){
            o_inxmc = InputOutputManager::Get().interceptXmc;
            o_inzmc = InputOutputManager::Get().interceptZmc;
            o_slxmc = InputOutputManager::Get().slopeXmc;
            o_slzmc = InputOutputManager::Get().slopeZmc;
            o_t0mc = InputOutputManager::Get().t0mc;
        }
        otree->Fill();
        MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    // Draw the plots of ADC & drift time distributions
    TCanvas * canv = new TCanvas("canv","",1024,800);
    canv->Divide(3,3);
    for (int lid = 1; lid<=8; lid++){
        for (int iwire = 0; iwire<wiresToConsider[lid].size(); iwire++){
            int wid = wiresToConsider[lid][iwire];
            canv->cd(1);gPad->SetGridx(1);gPad->SetGridy(1); gPad->SetLogy(1);
            h_ADCsumPkt_all[lid][wid]->Draw();
            h_ADCsumPkt_sig[lid][wid]->Draw("SAME");
            h_ADCsumPkt_bkg[lid][wid]->Draw("SAME");
            canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
            h_SVT_sig[lid][wid]->Draw("COLZ");
            canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
            h_SVT_bkg[lid][wid]->Draw("COLZ");
            canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1); gPad->SetLogy(1);
            h_ADCsumAll_all[lid][wid]->Draw();
            h_ADCsumAll_sig[lid][wid]->Draw("SAME");
            h_ADCsumAll_bkg[lid][wid]->Draw("SAME");
            canv->cd(5);gPad->SetGridx(1);gPad->SetGridy(1);
            h_AVT_sig[lid][wid]->Draw("COLZ");
            canv->cd(6);gPad->SetGridx(1);gPad->SetGridy(1);
            h_AVT_bkg[lid][wid]->Draw("COLZ");
            canv->cd(7);gPad->SetGridx(1);gPad->SetGridy(1);
            h_DriftT_all[lid][wid]->Draw();
            h_DriftT_sig[lid][wid]->Draw("SAME");
            h_DriftT_bkg[lid][wid]->Draw("SAME");
            canv->cd(8);gPad->SetGridx(1);gPad->SetGridy(1); gPad->SetLogz(1);
            h_SVT_siglate[lid][wid]->Draw("COLZ");
            canv->cd(9);gPad->SetGridx(1);gPad->SetGridy(1); gPad->SetLogz(1);
            h_AVT_siglate[lid][wid]->Draw("COLZ");
            canv->SaveAs(Form("hits.%d.%s.layer%d.l%dw%d.png",m_runNo,m_runName.Data(),m_testLayer,lid,wid));
        }
    }
    // Draw the nHits distribution
    canv = new TCanvas("canv","",1024,800);
    h_nGoodHits->Draw("COLZ");
    canv->SaveAs(Form("nHitsG.%d.%s.layer%d.png",m_runNo,m_runName.Data(),m_testLayer));
    // Draw the hit map
    canv = new TCanvas("canv","",1024,800);
    gPad->SetGridx(1);gPad->SetGridy(1);
    h_hitmap->Draw("COLZ");
    TGaxis * newXaxis = new TGaxis(h_hitmap->GetXaxis()->GetXmax(),h_hitmap->GetYaxis()->GetXmin(),h_hitmap->GetXaxis()->GetXmin(),h_hitmap->GetYaxis()->GetXmin()
                                 ,-h_hitmap->GetXaxis()->GetXmax(),-h_hitmap->GetXaxis()->GetXmin(),510,"-");
    newXaxis->SetLabelOffset(-0.03);
    newXaxis->SetLabelFont(42);
    newXaxis->SetTextFont(42);
    newXaxis->SetTitle("Cell ID + DOCA/CellSize");
    newXaxis->Draw();
    canv->SaveAs(Form("map.%d.%s.layer%d.png",m_runNo,m_runName.Data(),m_testLayer));

    // Save the histograms and tree
    for (int lid = 1; lid<=8; lid++){
        for (int wid = 0; wid<=10; wid++){
            h_DriftT_all[lid][wid]->Write();
            h_DriftT_bkg[lid][wid]->Write();
            h_DriftT_sig[lid][wid]->Write();
            h_ADCsumPkt_all[lid][wid]->Write();
            h_ADCsumPkt_bkg[lid][wid]->Write();
            h_ADCsumPkt_sig[lid][wid]->Write();
            h_ADCsumAll_all[lid][wid]->Write();
            h_ADCsumAll_bkg[lid][wid]->Write();
            h_ADCsumAll_sig[lid][wid]->Write();
            h_AVT_bkg[lid][wid]->Write();
            h_AVT_sig[lid][wid]->Write();
        }
    }
    h_hitmap->Write();
    h_nGoodHits->Write();
    otree->Write();;
    ofile->Close();

    return 0;
}

bool isGoodHit(int iHit){ // Here I neglected the cut on ipeak, layerID cause when I call this function I'm counting hits in a selected cell or a cell in test layer
    if (iHit<0||iHit>=InputOutputManager::Get().nHits) return false;
    if (InputOutputManager::Get().DriftT->at(iHit)<t_min||InputOutputManager::Get().DriftT->at(iHit)>t_max) return false;
    if (InputOutputManager::Get().ADCsumAll->at(iHit)!=0){ // TODO: this is a temporary way to check if the input hit file is MC data
        if (InputOutputManager::Get().ADCsumAll->at(iHit)<aaCut) return false;
        if (InputOutputManager::Get().ADCsumPacket->at(iHit)<sumCut) return false;
    }
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

int GetCandidate(TString candSelBy, int nHitsS){
    int cand = 0;
    if (candSelBy=="Original"){
        cand = 0;
        if (nHitsS>0){
            cand = -1;
            for (int iCand = 0; iCand<InputOutputManager::Get().nCandidatesFound; iCand++){
                if (InputOutputManager::Get().nHitsS[iCand] == nHitsS){
                    cand = iCand;
                    break;
                }
            }
        }
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
    else if (candSelBy=="FittingChi2"||candSelBy=="GlobalChi2"||candSelBy=="GlobalChi2WithTestLayer"){
        double minchi2 = 1e9;
        int minNhitsS = 0;
        for (int iCand = 0; iCand<InputOutputManager::Get().nCandidatesFound; iCand++){
            if (candSelBy=="GlobalChi2"){
                if ((minchi2>InputOutputManager::Get().chi2a[iCand]&&minNhitsS==InputOutputManager::Get().nHitsS[iCand])||minNhitsS<InputOutputManager::Get().nHitsS[iCand]){
                    cand = iCand;
                    minchi2 = InputOutputManager::Get().chi2a[iCand];
                    minNhitsS = InputOutputManager::Get().nHitsS[iCand];
                }
            }
            else if (candSelBy=="GlobalChi2WithTestLayer"){
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
    fprintf(stderr,"\t -W \n");
    fprintf(stderr,"\t\t Seperate wires\n");
    fprintf(stderr,"\t -S \n");
    fprintf(stderr,"\t\t Use side-by-side XT\n");
    fprintf(stderr,"\t -H <suf>\n");
    fprintf(stderr,"\t\t Add suffix to the hit file, like h_100SUFFIX.root\n");
    fprintf(stderr,"\t -m\n");
    fprintf(stderr,"\t\t Treat the input sample as MC sample (by default it's considered as real data)\n");
}

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        ParameterManager::Get().ReadInputFile(configureFile,"",false,false);
        ParameterManager::Get().LoadParameters(ParameterManager::kTracking);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTManager);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTAnalyzer);
    }
}
