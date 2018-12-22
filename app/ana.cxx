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
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "XTAnalyzer.hxx"
#include "MyProcessManager.hxx"
#include "MyRuntimeParameters.hxx"
#include "Log.hxx"
#include "header.hxx"

#define NBINS    20
#define MAXTRUNC 6

// run time parameters
int m_runNo = 0;
TString m_prerunname  = "prerun";
TString m_runname = "currun";
TString m_suffix = "";
int m_iEntryStart = 0;
int m_iEntryStop = 0;
int m_verboseLevel = 0;
int m_modulo = 10000;
bool m_memdebug = false;
int m_saveHists = 0;
bool m_outputEventTree = false;
int m_defaultLayerID = 4;
int m_geoSetup = 0; // 0: normal scintillator; 1: finger scintillator
int m_inputType = 0; // 1 for MC; 0 for data
int m_xtType = 55; // XYZ means polX for center, polY for middle and polZ for tail. If X is 0 then let middle function fit the center region.
bool m_AsymXT = false; // use asymmetric xt curve or not
TString m_CandSelBy = "Original"; // find the candidate with the smallest chi2 without the testlayer; otherwise use chi2 wit the testlayer (default order from tracking);
bool m_ClosestPeak = false; // To get XT: find the peak with the smallest residual. Otherwise choose the first one over threshold; NOTE: in analysis part, it's always the closest peak chosen to get efficiency and residual
bool m_UpdateWireMap = false; // by default don't update the wire map (i.e. no calibration of wire position)
TString m_ExternalXT = ""; // if this is not an empty string, then skip generating xt file but use this one as external XTs
// for updating wiremap
double m_scale = 1; // normalize the step size with this scale when updating wiremap
double m_stepSize = 0;
double m_minDeltaSlz = 0;
double m_maxDeltaSlz = 0;
double m_minDeltaInx = 0;
double m_maxDeltaInx = 0;
std::map<int,bool> m_wireIDWhiteList; // a list of wire IDs to be updated in new wire map
// for cutting
bool m_RequireInTriggerCounter = true;
bool m_RequireAllGoldenHits = false;
int m_nHitsMax = 30;
int m_nHitsSmin = 7;
int m_sumCut = -30;
int m_aaCut = 0;
double m_maxchi2 = 2;
double m_maxslz = 0.1;
double m_maxFD = 6; // only count hits with DOCA smaller than 6 mm
int    m_tmaxSet = 0;
//for binning
double m_tmin = -25-1/0.96/2; // t range for one x bin
double m_tmax = 800+1/0.96/2;
double m_xmax = 10; // x range for one t bin
int    m_NbinT = 792+1;
int    m_NbinX = 256;
int    m_NbinRes = 256;
double m_minchi2p = 1;
double m_maxRes = 2;
// for ADC2Charge
TString m_adc2charge = "2.0158464*x";

// map for wire position
double  map_x[NLAY][NCEL][2];
double  map_y[NLAY][NCEL][2];
double  map_off[NLAY][NCEL];
double  map_slz[NLAY][NCEL]; // for wire calibration cuts
double  map_inx[NLAY][NCEL]; // for wire calibration cuts
double  map_inxmc[NLAY][NCEL]; // for wire calibration cuts
double  map_slzmc[NLAY][NCEL]; // for wire calibration cuts
bool    map_has[NLAY][NCEL];
int     map_ch[NLAY][NCEL];
int     map_bid[NLAY][NCEL];
double  npair_per_cm = 0;
TF1 * f_left = 0;
TF1 * f_right = 0;

//==================About RECBE======================
TF1 * fADC2ChargeFunction = 0;

//==================About Scintillator======================
// normal scintillator
double sciYup = 0;
double sciYdown = 0;
double sciHL = 0;
double sciHW = 0;

//===================for get dist============================
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

double ADC2Charge(double adc); // fC
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz);
double getGG(double aa, double slx, double slz);
int t2d(double t, double & d, bool isRight, double tmaxSet = 0);
double findFirstX(TF1 * f, double val, double xmin, double xmax, double delta);
void doFit(TH1D * h,double leftRatio = 1/3., double rightRatio = 1/3., double leftEnd = -2, double rightEnd = 2);
void print_usage(char * prog_name);
int getHitType(int type,bool isRight);
int GetCandidate(TString & candSelBy, std::vector<int> * layerID, std::vector<int> * type, std::vector<double> * fitD[NCAND], std::vector<int> * sel[NCAND], int * nHitsS, double * chi2, double * chi2a);
bool isInTriggerCounter(int geoSetup, double tinz, double tslz);
void getRunTimeParameters(TString configureFile);

MyProcessManager * pMyProcessManager;

int main(int argc, char** argv){
    int temp_geoSetup = 0; bool set_geoSetup = false;
    int temp_inputType = 0; bool set_inputType = false;
    int temp_xtType = 0; bool set_xtType = false;
    //for cutting
    int temp_nHitsMax = 0; bool set_nHitsMax = false;
    int temp_nHitsSmin = 0; bool set_nHitsSmin = false;
    int temp_sumCut = 0; bool set_sumCut = false;
    int temp_aaCut = 0; bool set_aaCut = false;
    double temp_maxchi2 = 0; bool set_maxchi2 = false;
    double temp_minchi2p = 0; bool set_minchi2p = false;
    double temp_maxRes = 0; bool set_maxRes = false;
    double temp_maxslz = 0; bool set_maxslz = false;
    double temp_maxFD = 0; bool set_maxFD = false;
    int    temp_tmaxSet = 0; bool set_tmaxSet = false;
    //for binning
    double temp_tmin = 0; bool set_tmin = false;
    double temp_tmax = 0; bool set_tmax = false;
    double temp_xmax = 0; bool set_xmax = false;
    int    temp_NbinT = 0; bool set_NbinT = false;
    int    temp_NbinX = 0; bool set_NbinX = false;
    int    temp_NbinRes = 0; bool set_NbinRes = false;
    bool temp_AsymXT = false; bool set_AsymXT = false;
    TString temp_CandSelBy = "Original"; bool set_CandSelBy = false;
    bool temp_ClosestPeak = false; bool set_ClosestPeak = false;

    std::map<std::string, Log::ErrorPriority> namedDebugLevel;
    std::map<std::string, Log::LogPriority> namedLogLevel;
    int    opt_result;
    while((opt_result=getopt(argc,argv,"M:R:B:E:L:H:C:TGF:WX:n:f:c:v:r:z:d:o:s:a:l:u:t:m:y:g:i:x:AS:PD:V:"))!=-1){
        switch(opt_result){
            /* INPUTS */
            case 'M':
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
            case 'L':
                m_defaultLayerID = atoi(optarg);
                printf("Test layer set to %d\n",m_defaultLayerID);
                break;
            case 'H':
                m_saveHists = atoi(optarg);
                printf("Histogram saving level set to %d\n",m_saveHists);
                break;
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'T':
                m_RequireInTriggerCounter = false;
                printf("Disable in-trigger cut\n");
            case 'G':
                m_RequireAllGoldenHits = true;
                printf("In getting offset and new xt, require all hits used in tracking are golden hits\n");
            case 'F':
                m_suffix = optarg;
                printf("Adding suffix \"%s\"\n",optarg);
                break;
            case 'W':
                m_UpdateWireMap = true;
                printf("Will update wire map\n");
                break;
            case 'X':
                m_ExternalXT = optarg;
                printf("Will not generate new xt file but use \"%s\"\n",optarg);
                break;
            case 'n':
                temp_nHitsMax = atoi(optarg);set_nHitsMax = true;
                printf("Maximum number of hits cut set to %d\n",temp_nHitsMax);
                break;
            case 'f':
                temp_nHitsSmin = atoi(optarg);set_nHitsSmin = true;
                printf("Minimum number of selected hits cut set to %d\n",temp_nHitsSmin);
                break;
            case 'c':
                temp_maxchi2 = atof(optarg);set_maxchi2 = true;
                printf("Maximum chi2 cut set to %.3e\n",temp_maxchi2);
                break;
            case 'v':
                temp_minchi2p = atof(optarg);set_minchi2p = true;
                printf("Minimum p-value cut set to %.3e\n",temp_minchi2p);
                break;
            case 'r':
                temp_maxRes = atof(optarg);set_maxRes = true;
                printf("Maximum resolution cut set to %.3e\n",temp_maxRes);
                break;
            case 'z':
                temp_maxslz = atof(optarg);set_maxslz = true;
                printf("Maximum y-z slope cut set to %.3e\n",temp_maxslz);
                break;
            case 'd':
                temp_maxFD = atof(optarg);set_maxFD = true;
                printf("Maximum fitD cut set to %.3e\n",temp_maxFD);
                break;
            case 'o':
                temp_tmaxSet = atoi(optarg);set_tmaxSet = true;
                printf("Maximum time range set to %d\n",temp_tmaxSet);
                break;
            case 's':
                temp_sumCut = atoi(optarg);set_sumCut = true;
                printf("ADC sum over peak cut set to %d\n",temp_sumCut);
                break;
            case 'a':
                temp_aaCut = atoi(optarg);set_aaCut = true;
                printf("ADC sum over all cut set to %d\n",temp_aaCut);
                break;
            case 'l':
                temp_tmin = atof(optarg);set_tmin = true;
                printf("Minimum time on axis set to %.3e\n",temp_tmin);
                break;
            case 'u':
                temp_tmax = atof(optarg);set_tmax = true;
                printf("Maximum time on axis set to %.3e\n",temp_tmax);
                break;
            case 't':
                temp_NbinT = atoi(optarg);set_NbinT = true;
                printf("Number of bins on time axis set to %d\n",temp_NbinT);
                break;
            case 'm':
                temp_NbinX = atoi(optarg);set_NbinX = true;
                printf("Number of bins on space axis set to %d\n",temp_NbinX);
                break;
            case 'y':
                temp_NbinRes = atoi(optarg);set_NbinRes = true;
                printf("Number of bins on resolution axis set to %d\n",temp_NbinRes);
                break;
            case 'g':
                temp_geoSetup = atoi(optarg);set_geoSetup = true;
                printf("Geometry setup set to %d\n",temp_geoSetup);
                break;
            case 'i':
                temp_inputType = atoi(optarg);set_inputType = true;
                printf("Input type set to %d\n",temp_inputType);
                break;
            case 'x':
                temp_xtType = atoi(optarg);set_xtType = true;
                printf("XT type set to %d\n",temp_xtType);
                break;
            case 'A':
                temp_AsymXT = true; set_AsymXT = true;
                printf("Use asymmetric XT\n");
                break;
            case 'S':
                temp_CandSelBy = optarg; set_CandSelBy = true;
                if (temp_CandSelBy!="Original"&&temp_CandSelBy!="FittingChi2"&&temp_CandSelBy!="GlobalChi2"&&temp_CandSelBy!="LeastLatePeak"){
                    MyNamedWarn("XT","The option \""<<temp_CandSelBy<<" is not known. Will choose Original");
                    temp_CandSelBy = "Original";
                }
                printf("Choose the candidate by %s\n",temp_CandSelBy.Data());
                break;
            case 'P':
                temp_ClosestPeak = true; set_ClosestPeak = true;
                printf("Use the peak with smallest residual to get driftT. Otherwise just use the first one\n");
                break;
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

    if (set_geoSetup) m_geoSetup = temp_geoSetup;
    if (set_inputType) m_inputType = temp_inputType;
    if (set_xtType) m_xtType = temp_xtType;
    if (set_AsymXT) m_AsymXT = temp_AsymXT;
    if (set_CandSelBy) m_CandSelBy = temp_CandSelBy;
    if (set_ClosestPeak) m_ClosestPeak = temp_ClosestPeak;
    // for cutting
    if (set_nHitsMax) m_nHitsMax = temp_nHitsMax;
    if (set_nHitsSmin) m_nHitsSmin = temp_nHitsSmin;
    if (set_sumCut) m_sumCut = temp_sumCut;
    if (set_aaCut) m_aaCut = temp_aaCut;
    if (set_maxchi2) m_maxchi2 = temp_maxchi2;
    if (set_maxslz) m_maxslz = temp_maxslz;
    if (set_maxFD) m_maxFD = temp_maxFD;
    if (set_tmaxSet) m_tmaxSet = temp_tmaxSet;
    //for binning
    if (set_tmin) m_tmin = temp_tmin;
    if (set_tmax) m_tmax = temp_tmax;
    if (set_xmax) m_xmax = temp_xmax;
    if (set_NbinT) m_NbinT = temp_NbinT;
    if (set_NbinX) m_NbinX = temp_NbinX;
    if (set_NbinRes) m_NbinRes = temp_NbinRes;
    if (set_minchi2p) m_minchi2p = temp_minchi2p;
    if (set_maxRes) m_maxRes = temp_maxRes;

    if (argc-optind<2){
        print_usage(argv[0]);
        return -1;
    }
    m_prerunname = argv[optind++];
    m_runname= argv[optind++];
    TString m_runnameout = m_runname;
    if (m_suffix!=""){
        m_runnameout=m_runname+"."+m_suffix;
    }

    if (argc-optind>0){
        for (;optind<argc; optind++){
            m_wireIDWhiteList[atoi(argv[optind])] = true;
        }
    }
    else{
        for (int i = 0; i<NCEL; i++){
            m_wireIDWhiteList[i] = true;
        }
    }

    printf("##############%s##################\n",argv[0]);
    printf("runNo       = %d\n",m_runNo);
    printf("prerunname  = \"%s\"\n",m_prerunname.Data());
    printf("runname     = \"%s\"\n",m_runname.Data());
    printf("suffix      = \"%s\"\n",m_suffix.Data());
    printf("default layer: %d\n",m_defaultLayerID);
    printf("geoSetup:     %s\n",m_geoSetup==0?"normal scintillator":"finger scintillator");
    printf("Use external XT? %s\n",m_ExternalXT==""?"No, generate new one":("Yes"+m_ExternalXT).Data());
    printf("Update wire map? %s\n",m_UpdateWireMap?"Yes":"No");
    printf("For wire map update, stepSize    = %.3e\n",m_stepSize);
    printf("For wire map update, min delta slz = %.3e\n",m_minDeltaSlz);
    printf("For wire map update, max delta slz = %.3e\n",m_maxDeltaSlz);
    printf("For wire map update, min delta inx = %.3e\n",m_minDeltaInx);
    printf("For wire map update, max delta inx = %.3e\n",m_maxDeltaInx);
    printf("For wire map update, wires to be calibrated: ");
    for (std::map<int,bool>::iterator i = m_wireIDWhiteList.begin(); i != m_wireIDWhiteList.end(); ++i) {
        printf("%d ",i->first);
    }
    printf("Require in trigger? %s\n",m_RequireInTriggerCounter?"Yes":"No");
    printf("Require all golden hit? %s\n",m_RequireAllGoldenHits?"Yes":"No");
    printf("xtType:       %d\n",m_xtType);
    printf("Asymmetric XT? %s\n",m_AsymXT?"yes":"no");
    printf("Candidate selected by? %s\n",m_CandSelBy.Data());
    printf("Choose Closest peak? %s\n",m_ClosestPeak?"yes":"no");
    printf("inputType   = %d, %s\n",m_inputType,m_inputType==0?"Real Data":"MC");
    printf("maxchi2     = %.3e\n",m_maxchi2);
    printf("maxslz      = %.3e\n",m_maxslz);
    printf("nHits max   = %d\n",m_nHitsMax);
    printf("nHitsSmin   = %d\n",m_nHitsSmin);
    printf("Q cut       = %d\n",m_aaCut);
    printf("tmaxSet     = %d\n",m_tmaxSet);
    printf("debug       = %d\n",m_verboseLevel);
    printf("print modulo= %d\n",m_modulo);
    printf("save fitting histograms at level %d\n",m_saveHists);
    printf("output EventTree? %s\n",m_outputEventTree?"yes":"no");
    printf("Entries:     [%d~%d]\n",m_iEntryStart,m_iEntryStop);
    printf("ADC -> Charge: %s\n",m_adc2charge.Data());
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    //=================================================Get related info========================================================
    // get run info
    TFile * if_run = new TFile(HOME+"/Input/run-info.root");
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
        if (i_runNo == m_runNo) break;
    }
    npair_per_cm = 60;
    TString gastype = "He:C_{2}H_{6}(50:50)";
    double W = 32; // eV
    if (gasID==1){
        gastype = "He:iC_{4}H_{10}(90:10)";
        npair_per_cm = 29;
        W = 39; // eV
    }
    else if (gasID==2){
        gastype = "He:CH_{4}(80:20)";
        npair_per_cm = 17;
        W = 39; // eV
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
    printf("runNo#%d: %s, %d, %s, %d V, %d mV, %.0f sec\n",m_runNo,gastype.Data(),runGr,duration.Data(),HV,THR,durationTime);

    //Prepare Maps
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            map_has[lid][wid] = false;
            map_off[lid][wid] = 0;
            map_ch[lid][wid] = -1;
            map_bid[lid][wid] = -1;
            map_slz[lid][wid] = 0;
            map_inx[lid][wid] = 0;
            map_inxmc[lid][wid] = 0;
            map_slzmc[lid][wid] = 0;
        }
    }

    // prepare for offset
    TH1D * h_celloff[NLAY][NCEL];
    TH1D * h_cellslz[NLAY][NCEL];
    TH1D * h_cellinx[NLAY][NCEL];
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            h_celloff[lid][wid] = new TH1D(Form("h_celloff_%d_%d",lid,wid),Form("Offset of wire [%d,%d]",lid,wid),128,-1,1);
            h_cellslz[lid][wid] = new TH1D(Form("h_cellslz_%d_%d",lid,wid),Form("slope Z of hits on wire [%d,%d]",lid,wid),128,-0.16,0.16);
            h_cellinx[lid][wid] = new TH1D(Form("h_cellinx_%d_%d",lid,wid),Form("intercetp X of hits on wire [%d,%d]",lid,wid),512,-50,50);
        }
    }

    // get slzmc inxmc of each wire
    TChain * ichain_beam = new TChain("t","t");
    ichain_beam->Add(Form("%s/Input/beammap.%d.root",HOME.Data(),m_runNo));
    if (ichain_beam->GetEntries()){
        double beam_inxmc;
        double beam_slzmc;
        int    beam_lid;
        int    beam_wid;
        ichain_beam->SetBranchAddress("inxmc",&beam_inxmc);
        ichain_beam->SetBranchAddress("slzmc",&beam_slzmc);
        ichain_beam->SetBranchAddress("lid",&beam_lid);
        ichain_beam->SetBranchAddress("wid",&beam_wid);
        for (int i = 0; i<ichain_beam->GetEntries(); i++){
            ichain_beam->GetEntry(i);
            if (beam_lid<0||beam_lid>=NLAY||beam_wid<0||beam_wid>=NCEL) continue;
            map_inxmc[beam_lid][beam_wid] = beam_inxmc;
            map_slzmc[beam_lid][beam_wid] = beam_slzmc;
        }
    }
    else{
        printf("Cannot find%s/Input/beammap.%d.root, will assume slzmc/inxmc center at 0\n",HOME.Data(),m_runNo);
    }

    // get wire map
    TChain * TChain_wirepos = new TChain("t","t");
    TChain_wirepos->Add(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),m_runNo,m_runname.Data()));
    if (!TChain_wirepos->GetEntries()){
        fprintf(stderr,"Cannot find %s/info/wire-position.%d.%s.root, using default one\n",HOME.Data(),m_runNo,m_runname.Data());
        TChain_wirepos->Add(Form("%s/Input/wire-position.root",HOME.Data()));
    }
    int     wp_bid;
    int     wp_ch;
    int     wp_wid;
    int     wp_lid;
    double  wp_xro;
    double  wp_yro;
    double  wp_xhv;
    double  wp_yhv;
    TChain_wirepos->SetBranchAddress("b",&wp_bid);
    TChain_wirepos->SetBranchAddress("ch",&wp_ch);
    TChain_wirepos->SetBranchAddress("l",&wp_lid);
    TChain_wirepos->SetBranchAddress("w",&wp_wid);
    TChain_wirepos->SetBranchAddress("xhv",&wp_xhv);
    TChain_wirepos->SetBranchAddress("yhv",&wp_yhv);
    TChain_wirepos->SetBranchAddress("xro",&wp_xro);
    TChain_wirepos->SetBranchAddress("yro",&wp_yro);
    for (int i = 0; i<TChain_wirepos->GetEntries(); i++){
        TChain_wirepos->GetEntry(i);
        if (wp_lid>=0&&wp_lid<NLAY&&wp_wid>=0&&wp_wid<NCEL){
            map_x[wp_lid][wp_wid][0] = wp_xhv;
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_x[wp_lid][wp_wid][1] = wp_xro;
            map_y[wp_lid][wp_wid][1] = wp_yro;
            map_ch[wp_lid][wp_wid] = wp_ch;
            map_bid[wp_lid][wp_wid] = wp_bid;
            map_has[wp_lid][wp_wid] = true;
        }
    }

    // prepare XT files
    TFile * preXTFile = 0;
    TFile * newXTFile = 0;
    TTree * newXTTree = 0;
    if (m_ExternalXT==""){
        preXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_prerunname.Data()));
        newXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_runnameout.Data()),"RECREATE");
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
    }
    else{
        newXTFile = new TFile(m_ExternalXT);
        if (!newXTFile||newXTFile->IsZombie()){
            MyError("External XT file \""<<m_ExternalXT<<"\" doesn't exist!");
            return -1;
        }
    }

    //===================Set scintillator geometry============================
    if (m_geoSetup==0){
        // normal scintillator
        sciYup = chamberCY+chamberHH+180; // mm
        sciYdown = chamberCY-chamberHH-180; 
        sciHL = 300/2.;
        sciHW = 90/2.;
    }
    else{
        // finger scintillator
        sciYup = chamberCY+chamberHH+250; // mm
        sciYdown = chamberCY-chamberHH-195; 
        sciHL = 33/2.;
        sciHW = 33/2.;
    }
    printf("    # Geometry:\n");
    printf("    sciYup      = %.3e\n",sciYup);
    printf("    sciYdown    = %.3e\n",sciYdown);
    printf("    sciHL       = %.3e\n",sciHL);
    printf("    sciHW       = %.3e\n",sciHW);

    // set RECBE ADC function
    fADC2ChargeFunction = new TF1("a2c",m_adc2charge,-10,800);

    //==============================================Prepare input file & output variables=================================================
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
    std::vector<double> * i_ped = 0;
    bool has_ped = false;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_aa = 0;
    std::vector<double> * i_driftD[NCAND] = {0};
    std::vector<double> * i_calD[NCAND] = {0};
    std::vector<double> * i_fitD[NCAND] = {0};
    std::vector<int> * i_sel[NCAND] = {0};
    int npairs[NCAND];
    int isel[NCAND];
    int icom[NCAND];
    double iinx[NCAND];
    double iinz[NCAND];
    double islx[NCAND];
    double islz[NCAND];
    double chi2x[NCAND];
    double chi2z[NCAND];
    double chi2i[NCAND];
    int nHitsS[NCAND];
    double inx[NCAND];
    double inz[NCAND];
    double slx[NCAND];
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

    // output file
    // the closest peak to the track in the test layer
    bool isGood = false; // if this event meets all cuts
    double res = 1e9;
    int    nHitsT = 0; // number of hits in the test layer is by default 1
    double theDD[2] = {1e9};
    double theDT[2] = {1e9};
    double theFD[2] = {1e9};
    int    theST[2] = {0}; // status of t-x translation
    double theAA[2] = {0};
    double theCharge = 0; // charge in the test layer
    bool   hashit[2] = {false};
    int    theWid[2] = {-1};
    double theSum[2] = {0};
    double thePeak[2] = {0};
    double theHeight[2] = {0};
    double sum1st[2] = {0};
    double dt1st[2] = {0};
    int theIp[2] = {0};
    int theMpi[2] = {0};
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
    // dE/dX related
    int nLayersOnTrack = 0; // number of layers used for chage on track
    double chargeOnTrack[NLAY]; // charge along the track
    double adcsumOnTrack[NLAY]; // ADC sum along the track
    int    chargeOnTrackIndex[NLAY]; // index of the corresponding hit along the track
    double theGG = 0;
    double trackGG = 0;
    // event info
    int nHitsSmallAll = 0;
    int nHitsSmallSASD = 0;
    int nSmallSumHits = 0;
    int nShadowedHits = 0;
    int nLateHits = 0;
    int nBoundaryHits = 0;
    int nSmallBoundaryHits = 0;
    std::vector<double> * o_driftD = 0;
    std::vector<int>    * o_driftDs = 0;
    std::vector<int> * o_channelID = 0;
    std::vector<int> * o_boardID = 0;

    //=================================================Prepare histograms, function and counters================================================
    TF1 * f_res = new TF1("fres","gaus",-m_maxRes,m_maxRes);
    TH1I * h_nHits = new TH1I("hnHits","Number of TDC hits in each event",100,0,100);
    TH1I * h_DOF = new TH1I("hDOF","Number of DOF",5,0,5);
    TH1D * h_chi2 = new TH1D("hchi2","#chi^{2} of fitting",256,0,10);
    TH1D * h_slz = new TH1D("hslz","Slope on z direction",256,-0.3,0.3);
    TH1D * h_DOCA = new TH1D("hDOCA","DOCA with Left/Right",m_NbinX,-m_xmax,m_xmax);
    TH1D * h_DOCAb = new TH1D("hDOCAb","DOCA without Left/Right",m_NbinX/2,0,m_xmax);
    TH1D * h_DriftD = new TH1D("hDriftD","Drift distance with Left/Right",m_NbinX,-m_xmax,m_xmax);
    TH1D * h_DriftDb = new TH1D("hDriftDb","Drift distance without Left/Right",m_NbinX/2,0,m_xmax);
    TH2D * h_aaVST = new TH2D("haaVST","ADC sum VS driftT",m_NbinT,m_tmin,m_tmax,200,-50,550);
    TH2D * h_aaVSD = new TH2D("haaVSD","ADC sum VS driftD",m_NbinX,0,m_xmax,200,-50,550);
    TH2D * h_ggVSX = new TH2D("hggVSX","Gas gain VS DOCA",m_NbinX,-m_xmax,m_xmax,256,0,2e5);
    TH1D * h_ggall = new TH1D("hggall","Gas gain",256,0,2e5);
    TH2D * h_xt = new TH2D("hxt","Time space relation",m_NbinT,m_tmin,m_tmax,m_NbinX,-m_xmax,m_xmax);
    TH2D * h_tx = new TH2D("htx","Space time relation",m_NbinX,-m_xmax,m_xmax,m_NbinT,m_tmin,m_tmax);
    TH2D * h_resVSX = new TH2D("hresVSX","Residual VS DOCA",m_NbinX,-m_xmax,m_xmax,m_NbinRes,-m_maxRes,m_maxRes);
    TH2D * h_resVSD = new TH2D("hresVSD","Residual VS drift distance",m_NbinX,-m_xmax,m_xmax,m_NbinRes,-m_maxRes,m_maxRes);
    TH1D * h_resD[NBINS];
    TH1D * h_resX[NBINS];
    TH1D * h_dedx[MAXTRUNC];
    for (int i = 0; i<NBINS; i++){
        double xmin = m_xmax*(i)/NBINS;
        double xmax = m_xmax*(i+1)/NBINS;
        h_resD[i] = new TH1D(Form("hresD%d",i),Form("Residual with drift distance in [%.1f,%.1f] mm",xmin,xmax),m_NbinRes,-m_maxRes,m_maxRes);
        h_resD[i]->GetXaxis()->SetTitle("Residual [mm]");
        h_resX[i] = new TH1D(Form("hresX%d",i),Form("Residual with DOCA in [%.1f,%.1f] mm",xmin,xmax),m_NbinRes,-m_maxRes,m_maxRes);
        h_resX[i]->GetXaxis()->SetTitle("Residual [mm]");
    }
    for (int i = 0; i<MAXTRUNC; i++){
        h_dedx[i] = new TH1D(Form("hdedx%d",i),Form("dEdX with %d hits omitted",i),256,0,3);
        h_dedx[i]->GetXaxis()->SetTitle("dE/dX [keV/cm]");
    }
    h_nHits->GetXaxis()->SetTitle("Number of Hits");
    h_DOF->GetXaxis()->SetTitle("DOF");
    h_chi2->GetXaxis()->SetTitle("#chi^{2}");
    h_slz->GetXaxis()->SetTitle("tan#theta_{z}");
    h_DOCA->GetXaxis()->SetTitle("DOCA with left(-)/right(+) [mm]");
    h_DOCAb->GetXaxis()->SetTitle("DOCA [mm]");
    h_DriftD->GetXaxis()->SetTitle("Drift distance [mm]");
    h_DriftDb->GetXaxis()->SetTitle("Drift distance with left(-)/right(+) [mm]");
    h_aaVST->GetXaxis()->SetTitle("Drift time [ns]");
    h_aaVST->GetYaxis()->SetTitle("ADC sum");
    h_aaVSD->GetXaxis()->SetTitle("Drift distance [mm]");
    h_aaVSD->GetYaxis()->SetTitle("ADC sum");
    h_ggVSX->GetXaxis()->SetTitle("DOCA [mm]");
    h_ggVSX->GetYaxis()->SetTitle("Gas gain");
    h_ggall->GetXaxis()->SetTitle("Gas gain");
    h_xt->GetXaxis()->SetTitle("Drift time [ns]");
    h_xt->GetYaxis()->SetTitle("DOCA [mm]");
    h_tx->GetXaxis()->SetTitle("DOCA [mm]");
    h_tx->GetYaxis()->SetTitle("Drift time [ns]");
    h_resVSX->GetXaxis()->SetTitle("DOCA [mm]");
    h_resVSX->GetYaxis()->SetTitle("Residual [mm]");
    h_resVSD->GetXaxis()->SetTitle("Drift distance [mm]");
    h_resVSD->GetYaxis()->SetTitle("Residual [mm]");

    int N_ALL = 0;
    int N_CUT1 = 0;
    int N_CUT2 = 0;
    int N_CUT3 = 0;
    int N_CUT4 = 0;
    int N_CUT5 = 0;
    int N_BIN[NBINS] = {0};

    //=================================================Loop in layers====================================================
    // Prepare XTAnalyzer
    XTAnalyzer * fXTAnalyzer = new XTAnalyzer(gasID,m_verboseLevel);
    // Loop in layers
    for (int testLayer = 0; testLayer<NLAY; testLayer++){
        //----------------------------------Set input file--------------------------------------------
        if (m_verboseLevel>0) {printf("In Layer %d: preparing input TChain\n",testLayer);fflush(stdout);}
        TChain * ichain = new TChain("t","t");
        ichain->Add(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),testLayer));
        ichain->GetEntries();
        Long64_t N = ichain->GetEntries();
        if (N==0){
            fprintf(stderr,"WARNING: \"%s/root/t_%d.%s.layer%d.root\" is empty! Will ignore this layer.\n",HOME.Data(),m_runNo,m_runname.Data(),testLayer);
            continue;
        }
        ichain->SetBranchAddress("triggerNumber",&triggerNumber);
        ichain->SetBranchAddress("nHits",&nHits);
        ichain->SetBranchAddress("nHitsG",&nHitsG);
        ichain->SetBranchAddress("layerID",&i_layerID);
        ichain->SetBranchAddress("wireID",&i_wireID);
        ichain->SetBranchAddress("driftT",&i_driftT);
        if (m_inputType) ichain->SetBranchAddress("driftDmc",&i_driftDmc);
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
        has_ped = (ichain->SetBranchAddress("ped",&i_ped)==0);
        ichain->SetBranchAddress("sum",&i_sum);
        ichain->SetBranchAddress("aa",&i_aa);
        for (int iCand = 0; iCand<NCAND; iCand++){
            ichain->SetBranchAddress(Form("driftD%d",iCand),&(i_driftD[iCand]));
            ichain->SetBranchAddress(Form("calD%d",iCand),&(i_calD[iCand]));
            ichain->SetBranchAddress(Form("fitD%d",iCand),&(i_fitD[iCand]));
            ichain->SetBranchAddress(Form("sel%d",iCand),&(i_sel[iCand]));
            ichain->SetBranchAddress(Form("npairs%d",iCand),&(npairs[iCand]));
            ichain->SetBranchAddress(Form("isel%d",iCand),&(isel[iCand]));
            ichain->SetBranchAddress(Form("icom%d",iCand),&(icom[iCand]));
            ichain->SetBranchAddress(Form("iinx%d",iCand),&(iinx[iCand]));
            ichain->SetBranchAddress(Form("iinz%d",iCand),&(iinz[iCand]));
            ichain->SetBranchAddress(Form("islx%d",iCand),&(islx[iCand]));
            ichain->SetBranchAddress(Form("islz%d",iCand),&(islz[iCand]));
            ichain->SetBranchAddress(Form("chi2x%d",iCand),&(chi2x[iCand]));
            ichain->SetBranchAddress(Form("chi2z%d",iCand),&(chi2z[iCand]));
            ichain->SetBranchAddress(Form("chi2i%d",iCand),&(chi2i[iCand]));
            ichain->SetBranchAddress(Form("nHitsS%d",iCand),&(nHitsS[iCand]));
            ichain->SetBranchAddress(Form("inx%d",iCand),&(inx[iCand]));
            ichain->SetBranchAddress(Form("inz%d",iCand),&(inz[iCand]));
            ichain->SetBranchAddress(Form("slx%d",iCand),&(slx[iCand]));
            ichain->SetBranchAddress(Form("slz%d",iCand),&(slz[iCand]));
            ichain->SetBranchAddress(Form("chi2%d",iCand),&(chi2[iCand]));
            ichain->SetBranchAddress(Form("chi2p%d",iCand),&(chi2p[iCand]));
            ichain->SetBranchAddress(Form("chi2a%d",iCand),&(chi2a[iCand]));
            if (m_inputType){
                ichain->SetBranchAddress(Form("chi2mc%d",iCand),&(chi2mc[iCand]));
                ichain->SetBranchAddress(Form("chi2pmc%d",iCand),&(chi2pmc[iCand]));
                ichain->SetBranchAddress(Form("chi2amc%d",iCand),&(chi2amc[iCand]));
            }
        }
        if (m_inputType){
            ichain->SetBranchAddress("inxmc",&inxmc);
            ichain->SetBranchAddress("inzmc",&inzmc);
            ichain->SetBranchAddress("slxmc",&slxmc);
            ichain->SetBranchAddress("slzmc",&slzmc);
        }

        //----------------------------------Get the offset--------------------------------------------
        MyNamedVerbose("Ana","##############The initial loop starts#############");
        // to get the offsets
        for ( int iEntry = 0 ; iEntry<N; iEntry++){
            if (iEntry%10000==0) printf("%d\n",iEntry);
            if (m_verboseLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

            // decide which candidate to use
            theCand = GetCandidate(m_CandSelBy, i_layerID, i_type, i_fitD, i_sel, nHitsS, chi2, chi2a);

            // ignore events with bad fitting
            if (nHitsS[theCand]<m_nHitsSmin) continue;
            if (chi2[theCand]>m_maxchi2) continue;
            //if (nHitsG>nHitsS[theCand]) continue;
            if (m_RequireInTriggerCounter&&!isInTriggerCounter(m_geoSetup,inz[theCand],slz[theCand])) continue;
            if (m_nHitsMax&&nHits>m_nHitsMax) continue;

            if (m_verboseLevel>=20) printf("  Good Event! Looping in %d hits\n",nHits);
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
                double tfitD = (*i_fitD[theCand])[ihit];
                double tdriftD = (*i_driftD[theCand])[ihit];
                if ((*i_sel[theCand])[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                if (tlayerID!=testLayer) continue;
                if (fabs(tfitD-tdriftD)<fabs(minres)){ // no cut for test layer!
                    minres = tfitD-tdriftD;
                    wireID = (*i_wireID)[ihit];
                    fitD = tfitD;
                    driftD = tdriftD;
                    driftT = (*i_driftT)[ihit];
                    has = true;
                }
            }
            if (!has) continue; // no hits found in test layer
            if (hasBadHit&&m_RequireAllGoldenHits) continue;

            if (m_verboseLevel>=20) printf("  Found hit! pushing to XTAnalyzer\n");
            if (fabs(driftD)>2&&fabs(driftD)<6&&abs(fitD-driftD)<1){ // only trust the body part
                h_celloff[testLayer][wireID]->Fill(fitD-driftD);
                h_cellslz[testLayer][wireID]->Fill(slz[theCand]);
                h_cellinx[testLayer][wireID]->Fill(inx[theCand]);
            }
        }

        // get offset
        for (int wid = 0; wid<NCEL; wid++){
            if (h_celloff[testLayer][wid]->GetEntries()<100) continue;
            map_off[testLayer][wid] = h_celloff[testLayer][wid]->GetMean();
            map_slz[testLayer][wid] = h_cellslz[testLayer][wid]->GetMean();
            map_inx[testLayer][wid] = h_cellinx[testLayer][wid]->GetMean();
        }

        if (m_ExternalXT==""){
            //----------------------------------Start to get XT--------------------------------------------
            //Initialize the analyzer
            int saveEvenOdd = 0; if (testLayer==4) saveEvenOdd = 1; else if (testLayer==5) saveEvenOdd = -1;
            int statusInitialize = fXTAnalyzer->Initialize(Form("%d.%s.layer%d",m_runNo,m_runnameout.Data(),testLayer),testLayer,preXTFile,newXTFile,newXTTree,m_xtType,!m_AsymXT,m_saveHists, testLayer==m_defaultLayerID, saveEvenOdd, testLayer!=0);
            if (statusInitialize){
                fprintf(stderr,"WARNING: something wrong with initializing XTAnalyzer for layer[%d], will ignore this layer!\n",testLayer);
                continue;
            }
            MyNamedVerbose("Ana","##############The fisrt loop starts: "<<N<<" entries#############");
            //for ( int iEntry = (m_iEntryStart?m_iEntryStart:0); iEntry<(m_iEntryStop?m_iEntryStop-1:N); iEntry++){
            for ( int iEntry = 0; iEntry<N; iEntry++){ // better not to skip anything before we get the XT file
                if (iEntry%m_modulo==0) printf("%d\n",iEntry);
                if (m_verboseLevel>=20) printf("Entry%d: \n",iEntry);
                ichain->GetEntry(iEntry);

                // decide which candidate to use
                theCand = GetCandidate(m_CandSelBy, i_layerID, i_type, i_fitD, i_sel, nHitsS, chi2, chi2a);

                // ignore events with bad fitting
                if (nHitsS[theCand]<m_nHitsSmin) continue;
                if (chi2[theCand]>m_maxchi2) continue;
                //if (nHitsG>nHitsS[theCand]) continue;
                if (m_RequireInTriggerCounter&&!isInTriggerCounter(m_geoSetup,inz[theCand],slz[theCand])) continue;
                if (m_nHitsMax&&nHits>m_nHitsMax) continue;

                if (m_verboseLevel>=20) printf("  Good Event! Looping in %d hits\n",nHits);
                // find the closest hit in the test layer
                double res_temp = 1e9;
                bool foundhit = false;
                double driftT, fitD;
                bool hasBadHit = false;
                bool wireChecked[NCEL] = {false};
                for (int ihit = 0; ihit<nHits; ihit++){
                    int tlayerID = (*i_layerID)[ihit];
                    int twireID = (*i_wireID)[ihit];
                    double tfitD = (*i_fitD[theCand])[ihit]-map_off[tlayerID][twireID];
                    double tdriftD = (*i_driftD[theCand])[ihit];
                    double tsum = (*i_sum)[ihit];
                    double taa = (*i_aa)[ihit];
                    if ((*i_sel[theCand])[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                    if (tlayerID!=testLayer) continue;
                    if (tsum<m_sumCut||taa<m_aaCut) continue;
                    if (m_ClosestPeak&&twireID<NCEL&&wireChecked[twireID]) continue; // for each cell only get the first peak over the threshold
                    if (fabs(tfitD-tdriftD)<fabs(res_temp)){ // Get the one with smallest residual
                        res_temp = tfitD-tdriftD;
                        fitD = tfitD;
                        driftT = (*i_driftT)[ihit];
                        if (twireID<NCEL) wireChecked[twireID] = true;
                        foundhit = true;
                    }
                }
                if (!foundhit) continue; // no hits found in test layer
                if (hasBadHit&&m_RequireAllGoldenHits) continue;

                if (m_verboseLevel>=20) printf("  Found hit! pushing to XTAnalyzer\n");
                // tell analyzer a new data point
                fXTAnalyzer->Push(driftT,fitD);
            }
            if (m_verboseLevel>0) printf("Starting XT analysis\n");
            // fit histograms/graphs, make plots, and save new xt file
            fXTAnalyzer->Process();
        }

        // get XT file
        TFile * XTFile = newXTFile;
        f_left = (TF1*) XTFile->Get(Form("flc_%d",testLayer));
        f_right = (TF1*) XTFile->Get(Form("frc_%d",testLayer));
        if (!f_left||!f_right){
            MyWarn("Cannot find the xt curves of this layer. Will load the default ones.");
            f_left = (TF1*) XTFile->Get("fl_0");
            f_right = (TF1*) XTFile->Get("fr_0");
        }
        if (!f_left||!f_right){
            MyError("Cannot find the default xt curves!");
            return 1;
        }
        double minDT = f_left->GetXmin()>f_right->GetXmin()?f_right->GetXmin():f_left->GetXmin();
        double maxDT =  f_left->GetXmax()<f_right->GetXmax()?f_right->GetXmax():f_left->GetXmax();
        if (m_tmaxSet)
            maxDT = m_tmaxSet;

        //----------------------------------prepare for output ROOT file--------------------------------------------
        TFile * ofile = new TFile(Form("%s/root/ana_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runnameout.Data(),testLayer),"RECREATE");
        TTree * otree = new TTree("t","t");
        otree->Branch("triggerNumber",&triggerNumber);
        otree->Branch("isGood",&isGood);
        otree->Branch("res",&res);
        otree->Branch("nHitsT",&nHitsT);
        otree->Branch("hashit",hashit,"hashit[nHitsT]/O");
        otree->Branch("theFD",theFD,"theFD[nHitsT]/D");
        otree->Branch("theDD",theDD,"theDD[nHitsT]/D");
        otree->Branch("theDT",theDT,"theDT[nHitsT]/D");
        otree->Branch("theST",theST,"theST[nHitsT]/I");
        otree->Branch("theAA",theAA,"theAA[nHitsT]/D");
        otree->Branch("theWid",&theWid,"theWid[nHitsT]/I");
        otree->Branch("theSum",&theSum,"theSum[nHitsT]/D");
        otree->Branch("sum1st",&sum1st,"sum1st[nHitsT]/D");
        otree->Branch("dt1st",&dt1st,"dt1st[nHitsT]/D");
        otree->Branch("thePeak",&thePeak,"thePeak[nHitsT]/D");
        otree->Branch("theHeight",&theHeight,"theHeight[nHitsT]/D");
        otree->Branch("theIp",&theIp,"theIp[nHitsT]/I");
        otree->Branch("theMpi",&theMpi,"theMpi[nHitsT]/I");
        otree->Branch("theCand",&theCand);
        otree->Branch("highBid",&highBid);
        otree->Branch("highCh",&highCh);
        otree->Branch("highLid",&highLid);
        otree->Branch("highWid",&highWid);
        otree->Branch("highIp",&highIp);
        otree->Branch("highSum",&highSum);
        otree->Branch("highAA",&highAA);
        otree->Branch("highDT",&highDT);
        otree->Branch("theCharge",&theCharge);
        for (int i = 0; i<NLAY; i++){
            otree->Branch(Form("chargeOnTrack%d",i),&(chargeOnTrack[i]));
            otree->Branch(Form("adcsumOnTrack%d",i),&(adcsumOnTrack[i]));
            otree->Branch(Form("chargeOnTrackIndex%d",i),&(chargeOnTrackIndex[i]));
        }
        otree->Branch("theGG",&theGG);
        otree->Branch("trackGG",&trackGG);
        otree->Branch("nLayers",&nLayersOnTrack);
        otree->Branch("nHitsSmallAll",&nHitsSmallAll);
        otree->Branch("nHitsSmallSASD",&nHitsSmallSASD);
        otree->Branch("nSHits",&nShadowedHits);
        otree->Branch("nLHits",&nLateHits);
        otree->Branch("nSSHits",&nSmallSumHits);
        otree->Branch("nBHits",&nBoundaryHits);
        otree->Branch("nSBHits",&nSmallBoundaryHits);
        otree->Branch("nHits",&nHits);
        otree->Branch("nHitsG",&nHitsG);
        otree->Branch("layerID",&i_layerID);
        otree->Branch("wireID",&i_wireID);
        otree->Branch("channelID",&o_channelID);
        otree->Branch("boardID",&o_boardID);
        otree->Branch("driftT",&i_driftT);
        if (m_inputType) otree->Branch("driftDmc",&i_driftDmc);
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
        if (has_ped) otree->Branch("ped",&i_ped);
        otree->Branch("sum",&i_sum);
        otree->Branch("aa",&i_aa);
        otree->Branch("driftD",&o_driftD);
        otree->Branch("driftDs",&o_driftDs);
        otree->Branch("driftD0",&(i_driftD[0]));
        otree->Branch("calD",&(i_calD[0]));
        otree->Branch("fitD",&(i_fitD[0]));
        otree->Branch("sel",&(i_sel[0]));
        otree->Branch("icom",&(icom[0]));
        otree->Branch("isel",&(isel[0]));
        otree->Branch("npairs",&(npairs[0]));
        otree->Branch("iinx",&(iinx[0]));
        otree->Branch("iinz",&(iinz[0]));
        otree->Branch("islx",&(islx[0]));
        otree->Branch("islz",&(islz[0]));
        otree->Branch("chi2x",&(chi2x[0]));
        otree->Branch("chi2z",&(chi2z[0]));
        otree->Branch("chi2i",&(chi2i[0]));
        otree->Branch("nHitsS",&(nHitsS[0]));
        otree->Branch("inx",&(inx[0]));
        otree->Branch("inz",&(inz[0]));
        otree->Branch("slx",&(slx[0]));
        otree->Branch("slz",&(slz[0]));
        otree->Branch("chi2",&(chi2[0]));
        otree->Branch("chi2p",&(chi2p[0]));
        otree->Branch("chi2a",&(chi2a[0]));
        if (m_inputType){
            otree->Branch("chi2mc",&(chi2mc[0]));
            otree->Branch("chi2pmc",&(chi2pmc[0]));
            otree->Branch("chi2amc",&(chi2amc[0]));
            otree->Branch("inxmc",&inxmc);
            otree->Branch("inzmc",&inzmc);
            otree->Branch("slxmc",&slxmc);
            otree->Branch("slzmc",&slzmc);
        }
        o_driftD = new std::vector<double>;
        o_driftDs = new std::vector<int>;
        o_channelID = new std::vector<int>;
        o_boardID = new std::vector<int>;

        MyNamedVerbose("Ana","##############The Second loop starts#############");
        double closeFD;
        double closeFD2;
        int    closeWid;
        int    closeWid2;
        double averageGG = 0;
        double averageGGErr = 0;
        int prevTheCand = 0;
        //----------------------------------Start the analysis to get residual and etc--------------------------------------------
        double closestchi2 = 1e9;
        for ( int iEntry = (m_iEntryStart?m_iEntryStart:0); iEntry<(m_iEntryStop?m_iEntryStop-1:N); iEntry++){
            if (iEntry%10000==0) printf("%d\n",iEntry);
            if (m_verboseLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

            // decide which candidate to use
            theCand = GetCandidate(m_CandSelBy, i_layerID, i_type, i_fitD, i_sel, nHitsS, chi2, chi2a);
            if (theCand!=prevTheCand){
                prevTheCand = theCand;
                otree->SetBranchAddress("driftD",&(i_driftD[theCand]));
                otree->SetBranchAddress("npairs",&(npairs[theCand]));
                otree->SetBranchAddress("isel",&(isel[theCand]));
                otree->SetBranchAddress("icom",&(icom[theCand]));
                otree->SetBranchAddress("iinx",&(iinx[theCand]));
                otree->SetBranchAddress("iinz",&(iinz[theCand]));
                otree->SetBranchAddress("islx",&(islx[theCand]));
                otree->SetBranchAddress("islz",&(islz[theCand]));
                otree->SetBranchAddress("chi2x",&(chi2x[theCand]));
                otree->SetBranchAddress("chi2z",&(chi2z[theCand]));
                otree->SetBranchAddress("chi2i",&(chi2i[theCand]));
                otree->SetBranchAddress("calD",&(i_calD[theCand]));
                otree->SetBranchAddress("nHitsS",&(nHitsS[theCand]));
                otree->SetBranchAddress("inx",&(inx[theCand]));
                otree->SetBranchAddress("inz",&(inz[theCand]));
                otree->SetBranchAddress("slx",&(slx[theCand]));
                otree->SetBranchAddress("slz",&(slz[theCand]));
                otree->SetBranchAddress("chi2",&(chi2[theCand]));
                otree->SetBranchAddress("chi2p",&(chi2p[theCand]));
                otree->SetBranchAddress("chi2a",&(chi2a[theCand]));
                if (m_inputType){
                    otree->SetBranchAddress("chi2mc",&(chi2mc[theCand]));
                    otree->SetBranchAddress("chi2pmc",&(chi2pmc[theCand]));
                    otree->SetBranchAddress("chi2amc",&(chi2amc[theCand]));
                    otree->SetBranchAddress("inxmc",&inxmc);
                    otree->SetBranchAddress("inzmc",&inzmc);
                    otree->SetBranchAddress("slxmc",&slxmc);
                    otree->SetBranchAddress("slzmc",&slzmc);
                }
                otree->SetBranchAddress("fitD",&(i_fitD[theCand]));
                otree->SetBranchAddress("sel",&(i_sel[theCand]));
            }

            // update m_minchi2p
            if (fabs(chi2[theCand]-m_maxchi2)<fabs(closestchi2-m_maxchi2)){
                closestchi2 = chi2[theCand];
                m_minchi2p = chi2p[theCand];
            }

            // ignore events with bad fitting
            isGood = true;
            N_ALL++;
            h_nHits->Fill(nHits);
            if (m_nHitsMax&&nHits>m_nHitsMax) isGood = false;
            if (isGood) N_CUT1++;
            if (isGood) h_DOF->Fill(nHitsS[theCand]-4);
            if (nHitsS[theCand]<m_nHitsSmin) isGood = false;
            if (isGood) N_CUT2++;
            if (isGood) h_chi2->Fill(chi2[theCand]);
            if (chi2[theCand]>m_maxchi2) isGood = false;
            if (isGood) N_CUT3++;
            if (isGood) h_slz->Fill(slz[theCand]);
            if (fabs(slz[theCand])>m_maxslz) isGood = false;
            if (isGood) N_CUT4++;

            // get closest wires
            closeFD = 1e3;
            closeWid = 0;
            closeFD2 = 1e3;
            closeWid2 = 0;
            for (int wid = 0; wid<NCEL; wid++){
                double fitD = get_dist(testLayer,wid,slx[theCand],inx[theCand],slz[theCand],inz[theCand]);
                if (fabs(closeFD)>fabs(fitD)){
                    closeFD2 = closeFD;
                    closeWid2 = closeWid;
                    closeWid = wid;
                    closeFD = fitD;
                }
                else if (fabs(closeFD2)>fabs(fitD)){
                    closeWid2 = wid;
                    closeFD2 = fitD;
                }
            }
            MyNamedVerbose("Ana","  iEntry = "<<iEntry<<", theCand = "<<theCand<<", "<<(isGood?"Good":"Bad")<<" event, closeWid = "<<closeWid<<", closeWid2 = "<<closeWid2<<", nHits = "<<nHits);
            theFD[0] = closeFD;
            int ibinX = fabs(closeFD)/m_xmax*NBINS; // bin index for DOCA
            if (ibinX>=NBINS) isGood = false; // too far away from any cell
            if (isGood) N_CUT5++;
            nHitsT = 1; // number of hits in the test layer is by default 1
            int ibinX2 = fabs(closeFD2)/m_xmax*NBINS; // bin index for the next closest DOCA
            if (fabs(closeFD2+closeFD)<4&&ibinX2<NBINS){
                nHitsT = 2; // there is one more hit in the test layer
                theFD[1] = closeFD2;
            }
            if (isGood){
                h_DOCA->Fill(closeFD);
                h_DOCAb->Fill(fabs(closeFD));
                N_BIN[ibinX]++;
                if (nHitsT==2){
                    h_DOCA->Fill(closeFD2);
                    h_DOCAb->Fill(fabs(closeFD2));
                    N_BIN[ibinX2]++;
                }
            }

            // find the signal hits (possibly 2)
            // record drift and charge info
            double minRes[2] = {1e3}; // minimal residual in the test layer
            hashit[0] = false;
            hashit[1] = false;
            for (int lid = 0; lid<NLAY; lid++){
                chargeOnTrack[lid] = 0;
                adcsumOnTrack[lid] = 0;
                chargeOnTrackIndex[lid] = 0;
            }
            theCharge = 0; // charge in the test layer
            res = 1e9;
            highBid = -1;
            highCh = -1;
            highLid = -1;
            highWid = -1;
            highIp = -1;
            highAA=0;
            highSum=-1e9;
            highDT=0;
            nHitsSmallSASD = 0;
            nHitsSmallAll = 0;
            nSmallSumHits = 0;
            nShadowedHits = 0;
            nLateHits = 0;
            nSmallBoundaryHits = 0;
            nBoundaryHits = 0;
            for (int ihit = 0; ihit<nHits; ihit++){
                // get hit info updated
                int lid = (*i_layerID)[ihit];
                int wid = (*i_wireID)[ihit];
                int ch = map_ch[lid][wid];
                int bid = map_bid[lid][wid];
                double aa = (*i_aa)[ihit];
                double sum = (*i_sum)[ihit];
                double peak = (*i_peak)[ihit];
                double height = (*i_height)[ihit];
                double fitD = (*i_fitD[theCand])[ihit]-map_off[lid][wid];
                int type = getHitType((*i_type)[ihit],fitD>=0);
                double dt = (*i_driftT)[ihit];
                double dd;
                double dd0 = (*i_driftD[theCand])[ihit];
                int status = t2d(dt,dd,fitD>0);
                o_driftD->push_back(dd);
                o_driftDs->push_back(status);
                o_channelID->push_back(ch);
                o_boardID->push_back(bid);
                (*i_fitD[theCand])[ihit]=fitD;
                int ip = (*i_ip)[ihit];
                int mpi = (*i_mpi)[ihit];
                double charge = ADC2Charge(aa); // ADC -> charge = e*Nt*GG; while Nt = dE/W = dEdX*trackL/W.
                double fd = get_dist(lid,wid,slx[theCand],inx[theCand],slz[theCand],inz[theCand]);

                // accumulate charge info
                if ((*i_ip)[ihit]==0){ // to calculate the distance to track, only count first peaks
                    //if (fabs(fd)<CELLW/2+0.5){ // along the track: distance smaller than half cell size (plus safety margin 0.5 mm)
                    if (fabs(fd)<m_maxFD){ // along the track: distance smaller than half cell size (plus safety margin 0.5 mm)
                        // FIXME: should decide whether to include the boundary layers or not: slightly smaller ADC, why?
                        //if (lid>0&&lid<NLAY-1){ // don't count the last layer: guard layer
                        //if (lid>1&&lid<NLAY){ // don't count the first layer: guard layer
                        if (lid>1&&lid<NLAY-1){ // don't count the first layer and the last layer: guard layers
                            //if (lid>0&&lid<NLAY){ // count all layers
                            if (!chargeOnTrack[lid]||charge>chargeOnTrack[lid]) chargeOnTrackIndex[lid] = ihit;
                            chargeOnTrack[lid]+=charge;
                            adcsumOnTrack[lid]+=aa;
                        }
                        if (lid==testLayer){ // in test layer hits
                            theCharge+=charge;
                        }
                    }
                    if (isGood&&(*i_layerID)[ihit]==testLayer){
                        // aa
                        h_aaVST->Fill(dt,aa);
                        h_aaVSD->Fill(dd,aa);
                    }
                }

                // check the closest peak in the test layer
                if ((*i_layerID)[ihit]==testLayer){
                    int isig = -1;
                    if (wid==closeWid) isig = 0;
                    else if (wid==closeWid2&&nHitsT==2) isig = 1;
                    if (isig>=0){
                        double resi = fabs(dd)-fabs(fd);
                        if (status==0&&fabs(resi)<fabs(minRes[isig])){ // Should have cut for test layer! otherwise XT will not be well tuned
                            hashit[isig] = true;
                            theFD[isig] = fd;
                            theDD[isig] = dd;
                            theDT[isig] = dt;
                            theST[isig] = status;
                            theAA[isig] = aa;
                            minRes[isig] = resi;
                            theWid[isig] = wid;
                            theSum[isig] = sum;
                            thePeak[isig] = peak;
                            theHeight[isig] = height;
                            theIp[isig] = ip;
                            theMpi[isig] = mpi;
                            if (fabs(resi)<fabs(res)) res = resi;
                        }
                        if (ip==0){ // record the ADC sum and the drift time for the first peak in the on track cell in the test layer
                            sum1st[isig] = sum;
                            dt1st[isig] = dt;
                        }
                    }
                }

                // check the highest hit
                if (highSum<sum){
                    highBid = bid;
                    highCh = ch;
                    highLid = lid;
                    highWid = wid;
                    highIp = ip;
                    highAA=aa;
                    highSum=sum;
                    highDT=dt;
                }

                // get nXXXHits according to the original hit distance
                if ((*i_sel[theCand])[ihit]==1){
                    if((fabs(dd0)<0.5||fabs(dd0)>7.5))
                        nBoundaryHits++;
                    if((fabs(dd0)<0.25||fabs(dd0)>7.75))
                        nSmallBoundaryHits++;
                    if(ip!=0){
                        nLateHits++;
                    }
                    if((*i_mpi)[ihit]!=0)
                        nShadowedHits++;
                    if(has_rank&&(*i_rank)[ihit]!=0)
                        nSmallSumHits++;
                }

            }

            // set driftD and extra info
            o_driftD->clear();
            o_driftDs->clear();
            o_channelID->clear();
            o_boardID->clear();

            // get statistics relating to the ADC with the highest hit
            for (int ihit = 0; ihit<nHits; ihit++){
                int lid = (*i_layerID)[ihit];
                int wid = (*i_wireID)[ihit];
                int ch = map_ch[lid][wid];
                int bid = map_bid[lid][wid];
                double aa = (*i_aa)[ihit];
                double dt = (*i_driftT)[ihit];
                if (aa<35){
                    nHitsSmallAll++;
                    if (bid==highBid&&ch/8==highCh/8){
                        nHitsSmallSASD++;
                    }
                }
            }

            // Count number of layers used for charge on track
            nLayersOnTrack = 0;
            double totalCharge = 0;
            for (int lid = 1; lid<NLAY; lid++){
                if (chargeOnTrack[lid]) nLayersOnTrack++;
                totalCharge+=chargeOnTrack[lid];
            }

            // sort the layers by charge from small to large
            for (int lid = 1; lid<NLAY; lid++){ // ignore layer 0
                for (int ljd = lid+1; ljd<NLAY; ljd++){
                    if (chargeOnTrack[lid]>chargeOnTrack[ljd]){
                        double temp = chargeOnTrack[lid];
                        chargeOnTrack[lid] = chargeOnTrack[ljd];
                        chargeOnTrack[ljd] = temp;
                        temp = adcsumOnTrack[lid];
                        adcsumOnTrack[lid] = adcsumOnTrack[ljd];
                        adcsumOnTrack[ljd] = temp;
                        int tempi = chargeOnTrackIndex[lid];
                        chargeOnTrackIndex[lid] = chargeOnTrackIndex[ljd];
                        chargeOnTrackIndex[ljd] = tempi;
                    }
                }
            }

            // get gg
            theGG = 0;
            if (isGood&&theCharge){
                int isig = -1;
                if (hashit[0]) isig = 0;
                else if (hashit[1]) isig = 1;
                if (isig!=-1){
                    double gg = getGG(theCharge,slx[theCand],slz[theCand]);
                    h_ggVSX->Fill(theFD[isig],gg); // only record the gas gain VS x in the test layer
                    if (gg>theGG) theGG = gg;
                }
                double gg = getGG(totalCharge/nLayersOnTrack,slx[theCand],slz[theCand]);
                h_ggall->Fill(gg);
                trackGG = gg;
            }
            otree->Fill();
        }
        // get averageGG
        //averageGG = h_ggVSX->GetMean(2);
        //averageGGErr = h_ggVSX->GetRMS(2);
        averageGG = h_ggall->GetMean();
        averageGGErr = h_ggall->GetRMS();

        MyNamedVerbose("Ana","##############The Third loop starts#############");
        double trackCharge[MAXTRUNC] = {0};
        //----------------------------------loop again for filling histograms--------------------------------------------
        for ( int iEntry = (m_iEntryStart?m_iEntryStart:0); iEntry<(m_iEntryStop?m_iEntryStop-1:N); iEntry++){
            otree->GetEntry(iEntry);
            MyNamedVerbose("Ana","  iEntry = "<<iEntry<<", "<<(isGood?"Good":"Bad")<<" event, nHits = "<<nHits<<", "<<nHitsT<<" hits in the test layer on track path");
            if (!isGood) continue; // not successfully reconstructed
            // get the truncated charge
            double totalCharge = 0;
            for (int itrunc = 0; itrunc < MAXTRUNC; itrunc++){
                trackCharge[itrunc] = 0;
            }
            for (int lid = 1; lid<NLAY; lid++){
                totalCharge+=chargeOnTrack[lid];
                if (NLAY-1-lid>=0&&NLAY-1-lid<MAXTRUNC){
                    trackCharge[NLAY-1-lid] = totalCharge;
                }
            }
            // Fill histograms if needed;
            for (int itrunc = 0; itrunc < MAXTRUNC; itrunc++){
                if (!trackCharge[itrunc]) continue;
                double theDE = trackCharge[itrunc]*1e-15/averageGG/1.6e-19*W;
                double theDX = (nLayersOnTrack-itrunc)*CELLH*sqrt(1+slx[theCand]*slx[theCand]+slz[theCand]*slz[theCand]);
                h_dedx[itrunc]->Fill(theDE/1000/(theDX/10));
            }
            for (int i = 0; i<nHitsT; i++){
                if (!hashit[i]) continue;
                int status = theST[i];
                double fd = theFD[i];
                double dd = theDD[i];
                double dt = theDT[i];
                double aa = theAA[i];
                if (status) continue; // out of range
                // if (aa<m_aaCut) continue; // too small
                h_DriftD->Fill(dd);
                h_DriftDb->Fill(fabs(dd));
                //dd = (*i_driftD)[ihit];
                double resi = fabs(dd) - fabs(fd);
                // xt
                h_tx->Fill(fd,dt);
                h_xt->Fill(dt,fd);
                // resi VS x/d
                h_resVSX->Fill(fd,resi);
                h_resVSD->Fill(dd,resi);
                // resi
                int ibx = fabs(fd)/m_xmax*NBINS;
                int ib = fabs(dd)/m_xmax*NBINS;
                if (ib>=0&&ib<NBINS) h_resD[ib]->Fill(resi);
                if (ibx>=0&&ibx<NBINS) h_resX[ibx]->Fill(resi);
            }
        }

        otree->Write();
        ofile->Close();

        //=================================================Get bin by bin information====================================================
        ofile = new TFile(Form("%s/info/resi_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runnameout.Data(),testLayer),"RECREATE");
        otree = new TTree("t","t");
        int o_ibin;
        double o_xmin;
        double o_xmid;
        double o_xmax;
        double o_xres;
        double o_xreserr;
        double o_xrms;
        double o_xrmserr;
        double o_xeff;
        double o_xeff3sig;
        double o_xeff3sigerr;
        double o_xeff5sig;
        double o_xeff500um;
        double o_xeff500umerr;
        double o_xeff1mm;
        double o_xoff;
        double o_dres;
        double o_dreserr;
        double o_drms;
        double o_drmserr;
        double o_deff3sig;
        double o_deff5sig;
        double o_deff500um;
        double o_deff1mm;
        double o_doff;
        double o_nx;
        double o_nxh;
        double o_nd;
        otree->Branch("ibin",&o_ibin);
        otree->Branch("xmin",&o_xmin);
        otree->Branch("xmid",&o_xmid);
        otree->Branch("xmax",&o_xmax);
        otree->Branch("xres",&o_xres);
        otree->Branch("xreserr",&o_xreserr);
        otree->Branch("xrms",&o_xrms);
        otree->Branch("xrmserr",&o_xrmserr);
        otree->Branch("xeff",&o_xeff);
        otree->Branch("xeff3sig",&o_xeff3sig);
        otree->Branch("xeff3sigerr",&o_xeff3sigerr);
        otree->Branch("xeff5sig",&o_xeff5sig);
        otree->Branch("xeff500um",&o_xeff500um);
        otree->Branch("xeff500umerr",&o_xeff500umerr);
        otree->Branch("xeff1mm",&o_xeff1mm);
        otree->Branch("xoff",&o_xoff);
        otree->Branch("dres",&o_dres);
        otree->Branch("dreserr",&o_dreserr);
        otree->Branch("drms",&o_drms);
        otree->Branch("drmserr",&o_drmserr);
        otree->Branch("deff3sig",&o_deff3sig);
        otree->Branch("deff5sig",&o_deff5sig);
        otree->Branch("deff500um",&o_deff500um);
        otree->Branch("deff1mm",&o_deff1mm);
        otree->Branch("doff",&o_doff);
        otree->Branch("nx",&o_nx);
        otree->Branch("nxh",&o_nxh);
        otree->Branch("nd",&o_nd);

        std::vector<double> v_xx;
        std::vector<double> v_xxerr;
        std::vector<double> v_dx;
        std::vector<double> v_dxerr;
        std::vector<double> v_xeff;
        std::vector<double> v_xeff3s;
        std::vector<double> v_xeff500um;
        std::vector<double> v_xres;
        std::vector<double> v_xrms;
        std::vector<double> v_xreserr;
        std::vector<double> v_xrmserr;
        std::vector<double> v_xoff;
        std::vector<double> v_deff;
        std::vector<double> v_deff3s;
        std::vector<double> v_deff500um;
        std::vector<double> v_dres;
        std::vector<double> v_drms;
        std::vector<double> v_dreserr;
        std::vector<double> v_drmserr;
        std::vector<double> v_doff;

        int ibinl = 0;
        int ibinr = 0;
        TCanvas * canv_bin = 0;
        if (m_saveHists){
            canv_bin = new TCanvas("canv_bin","canv_bin",1024,768);
            gStyle->SetPalette(1);
            gStyle->SetOptStat(0);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
            gStyle->SetOptFit(1);
            gPad->SetGridx(1);
            gPad->SetGridy(1);
        }
        int NusedD = 0;
        int NusedX = 0;
        double averageEffD = 0;
        double averageRMSD = 0;
        double averageEff3sD = 0;
        double averageResD = 0;
        double averageEffX = 0;
        double averageEffXErr = 0;
        double averageRMSX = 0;
        double averageRMSXErr = 0;
        double averageEff3sX = 0;
        double averageEff3sXErr = 0;
        double averageResX = 0;
        double averageResXErr = 0;
        double bestEffD = 0;
        double bestResD = 1e6;
        double bestRMSD = 1e6;
        double bestEffX = 0;
        double bestResX = 1e6;
        double bestRMSX = 1e6;
        for (int ibin = 0; ibin<NBINS; ibin++){
            o_ibin = ibin;
            o_xmin = m_xmax*(ibin)/NBINS;
            o_xmid = m_xmax*(ibin+0.5)/NBINS;
            o_xmax = m_xmax*(ibin+1)/NBINS;
            o_nx = N_BIN[ibin];
            o_nxh = h_resX[ibin]->GetEntries();
            o_nd = h_resD[ibin]->GetEntries();
            o_xrms = h_resX[ibin]->GetRMS();
            o_xrmserr = h_resX[ibin]->GetRMSError();
            o_drms = h_resD[ibin]->GetRMS();
            o_drmserr = h_resD[ibin]->GetRMSError();
            if (o_nx&&h_resX[ibin]->Integral()>0){
                o_xeff = (double)o_nxh/o_nx;
                if (o_xmid<0.5)
                    doFit(h_resX[ibin],3/4.,1/4.,0,m_maxRes);
                else if (o_xmid>7)
                    doFit(h_resX[ibin],1/3.,2/3.,-m_maxRes,m_maxRes);
                else
                    doFit(h_resX[ibin],1/3.,1/3.,-m_maxRes,m_maxRes);
                o_xres = f_res->GetParameter(2);
                o_xreserr = f_res->GetParError(2);
                o_xoff = f_res->GetParameter(1);
                ibinl = h_resX[ibin]->FindBin(o_xoff-o_xrms*3); 
                ibinr = h_resX[ibin]->FindBin(o_xoff+o_xrms*3); 
                o_xeff3sig = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
                o_xeff3sigerr = sqrt(o_xeff3sig*(1-o_xeff3sig)/o_nx);
                ibinl = h_resX[ibin]->FindBin(o_xoff-o_xrms*5); 
                ibinr = h_resX[ibin]->FindBin(o_xoff+o_xrms*5); 
                o_xeff5sig = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
                ibinl = h_resX[ibin]->FindBin(o_xoff-0.5); 
                ibinr = h_resX[ibin]->FindBin(o_xoff+0.5); 
                o_xeff500um = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
                o_xeff500umerr = sqrt(o_xeff500um*(1-o_xeff500um)/o_nx);
                ibinl = h_resX[ibin]->FindBin(o_xoff-1); 
                ibinr = h_resX[ibin]->FindBin(o_xoff+1); 
                o_xeff1mm = h_resX[ibin]->Integral(ibinl,ibinr)/o_nx;
                if (m_saveHists){
                    h_resX[ibin]->Draw();
                    canv_bin->SaveAs(Form("resX%d_%d.%s.layer%d.png",ibin,m_runNo,m_runnameout.Data(),testLayer));
                }
            }
            else{
                o_xres = 0;
                o_xoff = 0;
                o_xeff = 0;
                o_xeff3sig  = 0;
                o_xeff5sig  = 0;
                o_xeff500um = 0;
                o_xeff1mm   = 0;
            }
            if (h_resD[ibin]->Integral()>0){
                if (o_xmid<0.5)
                    doFit(h_resD[ibin],1/4.,3/4.,-m_maxRes,m_maxRes);
                else if (o_xmid>7)
                    doFit(h_resD[ibin],2/3.,1/3.,-m_maxRes,m_maxRes);
                else
                    doFit(h_resD[ibin],1/3.,1/3.,-m_maxRes,m_maxRes);
                o_dres = f_res->GetParameter(2);
                o_dreserr = f_res->GetParError(2);
                o_doff = f_res->GetParameter(1);
                ibinl = h_resD[ibin]->FindBin(o_doff-o_drms*3); 
                ibinr = h_resD[ibin]->FindBin(o_doff+o_drms*3); 
                o_deff3sig = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
                ibinl = h_resD[ibin]->FindBin(o_doff-o_drms*5); 
                ibinr = h_resD[ibin]->FindBin(o_doff+o_drms*5); 
                o_deff5sig = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
                ibinl = h_resD[ibin]->FindBin(o_doff-0.5); 
                ibinr = h_resD[ibin]->FindBin(o_doff+0.5); 
                o_deff500um = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
                ibinl = h_resD[ibin]->FindBin(o_doff-1); 
                ibinr = h_resD[ibin]->FindBin(o_doff+1); 
                o_deff1mm = h_resD[ibin]->Integral(ibinl,ibinr)/o_nd;
                if (m_saveHists){
                    h_resD[ibin]->Draw();
                    canv_bin->SaveAs(Form("resD%d_%d.%s.layer%d.png",ibin,m_runNo,m_runnameout.Data(),testLayer));
                }
            }
            else{
                o_dres = 0;
                o_doff = 0;
                o_deff3sig  = 0;
                o_deff5sig  = 0;
                o_deff500um = 0;
                o_deff1mm   = 0;
            }
            otree->Fill();
            if (o_nd>100){
                if (o_xmid<m_maxFD){
                    NusedD++;
                    averageResD+=o_dres;
                    averageRMSD+=o_drms;
                    averageEffD+=o_deff500um;
                    averageEff3sD+=o_deff3sig;
                    if (bestResD>o_dres) bestResD = o_dres;
                    if (bestRMSD>o_drms) bestRMSD = o_drms;
                    if (bestEffD<o_deff500um) bestEffD = o_deff500um;
                }
                v_dx.push_back(o_xmid);
                v_dxerr.push_back((o_xmax-o_xmin)/2.);
                v_deff3s.push_back(o_deff3sig);
                v_deff500um.push_back(o_deff500um);
                v_dres.push_back(o_dres);
                v_dreserr.push_back(o_dreserr);
                v_drms.push_back(o_drms);
                v_drmserr.push_back(o_drmserr);
                v_doff.push_back(o_doff);
            }
            if (o_nxh>100){
                if (o_xmid<m_maxFD){
                    NusedX++;
                    averageResX+=o_xres;
                    averageResXErr+=o_xreserr;
                    averageRMSX+=o_xrms;
                    averageRMSXErr+=o_xrmserr;
                    averageEffX+=o_xeff500um;
                    averageEffXErr+=o_xeff500umerr;
                    averageEff3sX+=o_xeff3sig;
                    averageEff3sXErr+=o_xeff3sigerr;
                    if (bestResX>o_xres) bestResX = o_xres;
                    if (bestRMSX>o_xrms) bestRMSX = o_xrms;
                    if (bestEffX<o_xeff500um) bestEffX = o_xeff500um;
                }
                v_xx.push_back(o_xmid);
                v_xxerr.push_back((o_xmax-o_xmin)/2.);
                v_xeff.push_back(o_xeff);
                v_xeff3s.push_back(o_xeff3sig);
                v_xeff500um.push_back(o_xeff500um);
                v_xres.push_back(o_xres);
                v_xreserr.push_back(o_xreserr);
                v_xrms.push_back(o_xrms);
                v_xrmserr.push_back(o_xrmserr);
                v_xoff.push_back(o_xoff);
            }
        }
        NusedD?averageEffD/=NusedD:averageEffD=0;
        NusedD?averageResD/=NusedD:averageResD=0;
        NusedD?averageRMSD/=NusedD:averageRMSD=0;
        NusedX?averageEff3sD/=NusedD:averageEff3sD=0;
        NusedX?averageEffX/=NusedX:averageEffX=0;
        NusedX?averageEffXErr/=NusedX:averageEffXErr=0;
        NusedX?averageEff3sX/=NusedX:averageEff3sX=0;
        NusedX?averageEff3sXErr/=NusedX:averageEff3sXErr=0;
        NusedX?averageResX/=NusedX:averageResX=0;
        NusedX?averageRMSX/=NusedX:averageRMSX=0;
        NusedX?averageRMSXErr/=NusedX:averageRMSXErr=0;
        // print out the result
        printf("=>  %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s\n",
                    "runNo/I","testLayer/I","HV/I","THR/I","gasID/I","aaCut/I",
                    "averageGG","averageGGErr",
                    "averageEffX","averageEff3sX","averageRMSX","averageResX",
                    "averageEffXErr","averageEff3sXErr","averageRMSXErr","averageResXErr",
                    "averageEffD","averageEff3sD","averageRMSD","averageResD",
                    "bestEffD","bestRMSD","bestResD","bestEffX","bestRMSX","bestResX");
        printf("==> %4d  %d   %d   %2d  %d   %2d  %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
                    m_runNo,testLayer,HV,THR,gasID,m_aaCut,
                    averageGG,averageGGErr,
                    averageEffX,averageEff3sX,averageRMSX,averageResX,
                    averageEffXErr,averageEff3sXErr,averageRMSXErr,averageResXErr,
                    averageEffD,averageEff3sD,averageRMSD,averageResD,
                    bestEffD,bestRMSD,bestResD,bestEffX,bestRMSX,bestResX);

        //=================================================Draw====================================================
        TCanvas * canv_tracking = new TCanvas("canv_tracking","canv_tracking",1024,768);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        canv_tracking->Divide(2,2);
        canv_tracking->cd(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_nHits->Draw();
        TLine * line_nHits = new TLine(m_nHitsMax,0,m_nHitsMax,h_nHits->GetMaximum());
        line_nHits->SetLineColor(kRed);
        line_nHits->Draw("SAME");
        TLatex * text_nHits = new TLatex(m_nHitsMax,h_nHits->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT1,(double)N_CUT1/N_ALL*100));
        text_nHits->SetTextColor(kRed);
        text_nHits->SetTextSize(0.04);
        text_nHits->Draw("SAME");
        canv_tracking->cd(2);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_DOF->Draw();
        TLine * line_DOF = new TLine(m_nHitsSmin-4,0,m_nHitsSmin-4,h_DOF->GetMaximum());
        line_DOF->SetLineColor(kRed);
        line_DOF->Draw("SAME");
        TLatex * text_DOF = new TLatex(m_nHitsSmin-4,h_DOF->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT2,(double)N_CUT2/N_ALL*100));
        text_DOF->SetTextColor(kRed);
        text_DOF->SetTextSize(0.04);
        text_DOF->Draw("SAME");
        canv_tracking->cd(3);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_chi2->Draw();
        TLine * line_chi2 = new TLine(m_maxchi2,0,m_maxchi2,h_chi2->GetMaximum());
        line_chi2->SetLineColor(kRed);
        line_chi2->Draw("SAME");
        TLatex * text_chi2 = new TLatex(m_maxchi2,h_chi2->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT3,(double)N_CUT3/N_ALL*100));
        text_chi2->SetTextColor(kRed);
        text_chi2->SetTextSize(0.04);
        text_chi2->Draw("SAME");
        canv_tracking->cd(4);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_slz->Draw();
        TLine * line_slzl = new TLine(-m_maxslz,0,-m_maxslz,h_slz->GetMaximum());
        line_slzl->SetLineColor(kRed);
        line_slzl->Draw("SAME");
        TLine * line_slzr = new TLine(m_maxslz,0,m_maxslz,h_slz->GetMaximum());
        line_slzr->SetLineColor(kRed);
        line_slzr->Draw("SAME");
        TLatex * text_slz = new TLatex(m_maxslz,h_slz->GetMaximum()*0.7,Form("%d(%.1f%%)",N_CUT4,(double)N_CUT4/N_ALL*100));
        text_slz->SetTextColor(kRed);
        text_slz->SetTextSize(0.04);
        text_slz->Draw("SAME");
        canv_tracking->SaveAs(Form("track_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_tracking->SaveAs(Form("track_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TCanvas * canv_DOCA = new TCanvas("canv_DOCA","canv_DOCA",600,800);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        canv_DOCA->Divide(1,2);
        canv_DOCA->cd(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_DOCA->GetYaxis()->SetRangeUser(0,h_DOCA->GetMaximum()*1.1);
        h_DOCA->Draw();
        canv_DOCA->cd(2);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_DOCAb->GetYaxis()->SetRangeUser(0,h_DOCAb->GetMaximum()*1.1);
        h_DOCAb->Draw();
        canv_DOCA->SaveAs(Form("DOCA_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_DOCA->SaveAs(Form("DOCA_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));
        canv_DOCA->cd(1);
        h_DriftD->GetYaxis()->SetRangeUser(0,h_DriftD->GetMaximum()*1.1);
        h_DriftD->Draw();
        canv_DOCA->cd(2);
        h_DriftDb->GetYaxis()->SetRangeUser(0,h_DriftDb->GetMaximum()*1.1);
        h_DriftDb->Draw();
        canv_DOCA->SaveAs(Form("DriftD_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_DOCA->SaveAs(Form("DriftD_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TCanvas * canv_XT = new TCanvas("canv_XT","canv_XT",800,600);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_xt->Draw("COLZ");
        h_xt->GetXaxis()->SetRangeUser(minDT,maxDT);
        f_left->Draw("SAME");
        f_right->Draw("SAME");
        canv_XT->SaveAs(Form("XT_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_XT->SaveAs(Form("XT_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TCanvas * canv_TX = new TCanvas("canv_TX","canv_TX",600,800);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_tx->Draw("COLZ");
        h_tx->Draw("COLZ");
        h_tx->GetYaxis()->SetRangeUser(minDT,maxDT);
        canv_TX->SaveAs(Form("TX_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_TX->SaveAs(Form("TX_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TCanvas * canv_general = new TCanvas("canv_general","canv_general",1024,768);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
        h_resVSX->Draw("COLZ");
        canv_general->SaveAs(Form("resVSX_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("resVSX_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        h_resVSD->Draw("COLZ");
        canv_general->SaveAs(Form("resVSD_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("resVSD_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        h_aaVST->Draw("COLZ");
        TLine * line_aaVST = new TLine(m_tmin,m_aaCut,m_tmax,m_aaCut);
        line_aaVST->SetLineColor(kRed);
        line_aaVST->Draw("SAME");
        canv_general->SaveAs(Form("aaVST_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("aaVST_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        h_aaVSD->Draw("COLZ");
        TLine * line_aaVSD = new TLine(0,m_aaCut,m_xmax,m_aaCut);
        line_aaVSD->SetLineColor(kRed);
        line_aaVSD->Draw("SAME");
        canv_general->SaveAs(Form("aaVSD_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("aaVSD_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        h_ggVSX->Draw("COLZ");
        canv_general->SaveAs(Form("ggVSX_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("ggVSX_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TLegend* leg_dedx = new TLegend(0.7,0.7,0.9,0.9);
        leg_dedx->AddEntry(h_dedx[0],Form("All hits used"));
        int maxEntry = 0;
        for (int i = 0; i<MAXTRUNC; i++){
            int entries = h_dedx[i]->GetMaximum();
            if (maxEntry<entries) maxEntry = entries;
        }
        h_dedx[0]->SetLineColor(1);
        h_dedx[0]->GetYaxis()->SetRangeUser(0,maxEntry*1.1);
        h_dedx[0]->Draw();
        for (int i = 1; i<MAXTRUNC; i++){
            h_dedx[i]->SetLineColor(i+1);
            leg_dedx->AddEntry(h_dedx[i],Form("Neglecting %d hits",i));
            h_dedx[i]->Draw("SAME");
        }
        leg_dedx->Draw();
        canv_general->SaveAs(Form("dedx_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("dedx_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_xeff = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xeff[0]));
        TGraphErrors * g_xeff3s = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xeff3s[0]));
        TGraphErrors * g_xeff500 = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xeff500um[0]));
        g_xeff->SetName("gxeff");
        g_xeff3s->SetName("gxeff3s");
        g_xeff500->SetName("gxeff500um");
        g_xeff->SetTitle("Efficiency VS DOCA");
        g_xeff->GetXaxis()->SetTitle("DOCA [mm]");
        g_xeff->GetYaxis()->SetTitle("Efficiency");
        g_xeff->SetMarkerStyle(20);
        g_xeff->SetMarkerColor(kBlack);
        g_xeff->SetLineColor(kBlack);
        g_xeff->GetYaxis()->SetRangeUser(0,1.1);
        g_xeff->Draw("APL");
        g_xeff500->SetMarkerStyle(20);
        g_xeff500->SetMarkerColor(kRed);
        g_xeff500->SetLineColor(kRed);
        g_xeff500->Draw("PLSAME");
        g_xeff3s->SetMarkerStyle(20);
        g_xeff3s->SetMarkerColor(kBlue);
        g_xeff3s->SetLineColor(kBlue);
        g_xeff3s->Draw("PLSAME");
        TLegend * leg_xeff = new TLegend(0.1,0.1,0.5,0.5);
        leg_xeff->AddEntry(g_xeff,"Raw efficiency","PL");
        leg_xeff->AddEntry(g_xeff3s,Form("Efficiency with 3#sigma cut %.1f%%",averageEff3sX*100),"PL");
        leg_xeff->AddEntry(g_xeff500,Form("Efficiency with 500 um cut %.1f%%",averageEffX*100),"PL");
        leg_xeff->Draw("SAME");
        canv_general->SaveAs(Form("effx_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("effx_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_deff3s = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_deff3s[0]));
        TGraphErrors * g_deff500 = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_deff500um[0]));
        g_deff3s->SetName("gdeff3s");
        g_deff500->SetName("gdeff500um");
        g_deff500->SetTitle("Efficiency VS drift distance");
        g_deff500->GetXaxis()->SetTitle("Drift distance [mm]");
        g_deff500->GetYaxis()->SetTitle("Efficiency");
        g_deff500->SetMarkerStyle(20);
        g_deff500->SetMarkerColor(kRed);
        g_deff500->SetLineColor(kRed);
        g_deff500->GetYaxis()->SetRangeUser(0,1.1);
        g_deff500->Draw("APL");
        g_deff3s->SetMarkerStyle(20);
        g_deff3s->SetMarkerColor(kBlue);
        g_deff3s->SetLineColor(kBlue);
        g_deff3s->Draw("PLSAME");
        TLegend * leg_deff = new TLegend(0.1,0.1,0.5,0.5);
        leg_deff->AddEntry(g_deff3s,Form("Efficiency with 3#sigma cut %.1f%%",averageEff3sD*100),"PL");
        leg_deff->AddEntry(g_deff500,Form("Efficiency with 500 um cut %.1f%%",averageEffD*100),"PL");
        leg_deff->Draw("SAME");
        canv_general->SaveAs(Form("effd_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("effd_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_xres = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xres[0]),&(v_xxerr[0]),&(v_xreserr[0]));
        g_xres->SetName("gxres");
        g_xres->SetTitle("#sigma of residual VS DOCA");
        g_xres->GetXaxis()->SetTitle("Drift distance [mm]");
        g_xres->GetYaxis()->SetTitle("#sigma [mm]");
        g_xres->SetMarkerStyle(20);
        g_xres->SetMarkerColor(kBlack);
        g_xres->SetLineColor(kBlack);
        g_xres->Draw("APL");
        canv_general->SaveAs(Form("resx_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("resx_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_dres = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_dres[0]),&(v_dxerr[0]),&(v_dreserr[0]));
        g_dres->SetName("gdres");
        g_dres->SetTitle("#sigma of residual VS drift distance");
        g_dres->GetXaxis()->SetTitle("Drift distance [mm]");
        g_dres->GetYaxis()->SetTitle("#sigma [mm]");
        g_dres->SetMarkerStyle(20);
        g_dres->SetMarkerColor(kBlack);
        g_dres->SetLineColor(kBlack);
        g_dres->Draw("APL");
        canv_general->SaveAs(Form("resd_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("resd_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_xrms = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xrms[0]),&(v_xxerr[0]),&(v_xrmserr[0]));
        g_xrms->SetName("gxrms");
        g_xrms->SetTitle("RMS of residual VS DOCA");
        g_xrms->GetXaxis()->SetTitle("Drift distance [mm]");
        g_xrms->GetYaxis()->SetTitle("RMS [mm]");
        g_xrms->SetMarkerStyle(20);
        g_xrms->SetMarkerColor(kBlack);
        g_xrms->SetLineColor(kBlack);
        g_xrms->Draw("APL");
        canv_general->SaveAs(Form("rmsx_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("rmsx_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_drms = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_drms[0]),&(v_dxerr[0]),&(v_drmserr[0]));
        g_drms->SetName("gdrms");
        g_drms->SetTitle("RMS of residual VS drift distance");
        g_drms->GetXaxis()->SetTitle("Drift distance [mm]");
        g_drms->GetYaxis()->SetTitle("RMS [mm]");
        g_drms->SetMarkerStyle(20);
        g_drms->SetMarkerColor(kBlack);
        g_drms->SetLineColor(kBlack);
        g_drms->Draw("APL");
        canv_general->SaveAs(Form("rmsd_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("rmsd_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_xoff = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xoff[0]));
        g_xoff->SetName("gxoff");
        g_xoff->SetTitle("Offset VS DOCA");
        g_xoff->GetXaxis()->SetTitle("Drift distance [mm]");
        g_xoff->GetYaxis()->SetTitle("offset [mm]");
        g_xoff->SetMarkerStyle(20);
        g_xoff->SetMarkerColor(kBlack);
        g_xoff->SetLineColor(kBlack);
        g_xoff->Draw("APL");
        canv_general->SaveAs(Form("offx_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("offx_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        TGraphErrors * g_doff = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_doff[0]));
        g_doff->SetName("gdoff");
        g_doff->SetTitle("Offset VS drift distance");
        g_doff->GetXaxis()->SetTitle("Drift distance [mm]");
        g_doff->GetYaxis()->SetTitle("offset [mm]");
        g_doff->SetMarkerStyle(20);
        g_doff->SetMarkerColor(kBlack);
        g_doff->SetLineColor(kBlack);
        g_doff->Draw("APL");
        canv_general->SaveAs(Form("offd_%d.%s.layer%d.pdf",m_runNo,m_runnameout.Data(),testLayer));
        canv_general->SaveAs(Form("offd_%d.%s.layer%d.png",m_runNo,m_runnameout.Data(),testLayer));

        //=================================================Save====================================================
        for (int i = 0; i<NBINS; i++){
            h_resD[i]->Write();
            h_resX[i]->Write();
        }
        for (int i = 0; i<MAXTRUNC; i++){
            h_dedx[i]->Write();
        }
        h_resVSX->Write();
        h_resVSD->Write();
        h_xt->Write();
        h_tx->Write();
        h_aaVST->Write();
        h_aaVSD->Write();
        h_ggVSX->Write();
        h_ggall->Write();
        h_nHits->Write();
        h_DOF->Write();
        h_chi2->Write();
        h_slz->Write();
        h_DOCA->Write();
        h_DOCAb->Write();
        h_DriftD->Write();
        h_DriftDb->Write();
        g_deff500->Write();
        g_xeff->Write();
        g_xeff500->Write();
        g_dres->Write();
        g_xres->Write();
        g_drms->Write();
        g_xrms->Write();
        g_doff->Write();
        g_xoff->Write();
        otree->Write();
        ofile->Close();

        printf("All events %d\n",N_ALL);
        printf("nHits<=%d: %d (%.1f%%)\n",m_nHitsMax,N_CUT1,(double)N_CUT1/N_ALL*100);
        printf("nHitsS>=%d: %d (%.1f%%)\n",m_nHitsSmin,N_CUT2,(double)N_CUT2/N_ALL*100);
        printf("chi2<%.1f (pvalue>%.1f): %d (%.1f%%)\n",m_maxchi2,m_minchi2p,N_CUT3,(double)N_CUT3/N_ALL*100);
        printf("|slz|<=%.1f: %d (%.1f%%)\n",m_maxslz,N_CUT4,(double)N_CUT4/N_ALL*100);
        printf("DOCA<=%.1f: %d (%.1f%%)\n",m_xmax,N_CUT5,(double)N_CUT5/N_ALL*100);

        if (m_verboseLevel>=20) printf("Finished!\n");
    }

    //===================output offset file============================
    TFile * ofile = new TFile(Form("%s/info/offset.%d.%s.root",HOME.Data(),m_runNo,m_runnameout.Data()),"RECREATE");
    TTree * otree = new TTree("t","t");
    double o_off_delta;
    double o_off_slz;
    double o_off_inx;
    int o_off_lid;
    int o_off_wid;
    otree->Branch("d",&o_off_delta);
    otree->Branch("slz",&o_off_slz);
    otree->Branch("inx",&o_off_inx);
    otree->Branch("wid",&o_off_wid);
    otree->Branch("lid",&o_off_lid);
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            h_celloff[lid][wid]->Write();
            h_cellslz[lid][wid]->Write();
            h_cellinx[lid][wid]->Write();
            if (!map_off[lid][wid]) continue;
            o_off_wid = wid;
            o_off_lid = lid;
            o_off_delta = map_off[lid][wid];
            o_off_slz = map_slz[lid][wid];
            o_off_inx = map_inx[lid][wid];
            otree->Fill();
        }
    }
    otree->Write();
    ofile->Close();

    //===================Set Wire Position============================
    if (m_UpdateWireMap){
        TFile * TFile_wirepos = new TFile(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),m_runNo,m_runnameout.Data()),"RECREATE");
        TTree * TTree_wirepos = new TTree("t","t");
        double wp_xc;
        double wp_yc;
        TTree_wirepos->Branch("b",&wp_bid);
        TTree_wirepos->Branch("ch",&wp_ch);
        TTree_wirepos->Branch("l",&wp_lid);
        TTree_wirepos->Branch("w",&wp_wid);
        TTree_wirepos->Branch("xhv",&wp_xhv);
        TTree_wirepos->Branch("yhv",&wp_yhv);
        TTree_wirepos->Branch("xc",&wp_xc);
        TTree_wirepos->Branch("yc",&wp_yc);
        TTree_wirepos->Branch("xro",&wp_xro);
        TTree_wirepos->Branch("yro",&wp_yro);
        for (wp_lid=0; wp_lid<NLAY; wp_lid++){
            for (wp_wid=0; wp_wid<NCEL; wp_wid++){
                if (!map_has[wp_lid][wp_wid]) continue;
                wp_xhv = map_x[wp_lid][wp_wid][0];
                wp_yhv = map_y[wp_lid][wp_wid][0];
                wp_xro = map_x[wp_lid][wp_wid][1];
                wp_yro = map_y[wp_lid][wp_wid][1];
                wp_xc = (wp_xro+wp_xhv)/2.;
                wp_yc = (wp_yro+wp_yhv)/2.;
                wp_ch = map_ch[wp_lid][wp_wid];
                wp_bid = map_bid[wp_lid][wp_wid];
                if (m_wireIDWhiteList[wp_wid]){
                    double theOff = map_off[wp_lid][wp_wid]*m_scale;
                    double deltaSlz = map_slz[wp_lid][wp_wid]-map_slzmc[wp_lid][wp_wid];
                    double deltaInx = map_inx[wp_lid][wp_wid]-map_inxmc[wp_lid][wp_wid];
                    if (m_stepSize){
                        if(fabs(theOff)>m_stepSize) theOff = theOff>0?m_stepSize:-m_stepSize;
                    }
                    if ((m_minDeltaSlz||m_maxDeltaSlz)&&(deltaSlz>m_maxDeltaSlz||deltaSlz<m_minDeltaSlz)){
                        theOff = 0;
                    }
                    if ((m_minDeltaInx||m_maxDeltaInx)&&(deltaInx>m_maxDeltaInx||deltaInx<m_minDeltaInx)){
                        theOff = 0;
                    }
                    wp_xro += theOff;
                    wp_xc += theOff;
                    wp_xhv += theOff;
                }
                TTree_wirepos->Fill();
            }
        }
        TTree_wirepos->Write();
        TFile_wirepos->Close();
    }

    return 0;
}

double ADC2Charge(double adc){ // fC
    double charge = 0;
    if (adc<735.346)
        charge = fADC2ChargeFunction->Eval(adc);
    else
        charge = fADC2ChargeFunction->Eval(735);
    return charge;
}

double getGG(double charge, double slx, double slz){
    double nPairs = CELLH*sqrt(1+slx*slx+slz*slz)*npair_per_cm/10;
    double gg = charge*1e-15/1.6e-19/nPairs;
    return gg;
}

double get_dist(int lid, int wid, double slx, double inx, double slz, double inz)
{
    if (!map_has[lid][wid]) return 1e3;
    double xdown = inx-slx*(sciYup-sciYdown);
    double zdown = inz-slz*(sciYup-sciYdown);
    vTrackU.SetXYZ(inx,sciYup,inz);
    vTrackD.SetXYZ(xdown,sciYdown,zdown);
    vWireHV.SetXYZ(map_x[lid][wid][0],map_y[lid][wid][0],-chamberHL);
    vWireRO.SetXYZ(map_x[lid][wid][1],map_y[lid][wid][1],chamberHL);
    vTrack = vTrackD-vTrackU;
    vWire = vWireRO-vWireHV;
    vDist = vWireHV-vTrackU;
    vAxis = vWire.Cross(vTrack);
    double value = -vDist*(vAxis.Unit());
    return value-map_off[lid][wid];
}

int t2d(double t, double & d, bool isRight, double tmaxSet){
    int status = 0;
    TF1 * f = 0; // body

    if (isRight) f = f_right;
    else f = f_left;

    double tRight = f->GetXmax();
    double tmax = tRight;
    if (tmaxSet&&tmax>tmaxSet) tmax = tmaxSet;
    double tmin = f->GetXmin();

    if (t<tmin){
        status = -1;
        d = 0;
    }
    else if (t>tmax){
        status = 1;
        if (t>tRight) d = f->Eval(tRight);
        else d = f->Eval(t);
    }
    else{
        status = 0;
        d = f->Eval(t);
    }

    return status;
}

double findFirstX(TF1 * f, double val, double xmin, double xmax, double delta){
    double theX = 0;
    for (double x = xmin+delta; x<xmax; x+=delta){ // At least two solutions. Scan to find the smallest one
        theX = f->GetX(val,xmin,x);
        if (fabs(theX-x)>delta/10.&&fabs(theX-xmin)>delta/10.){
            break;
        }
    }
    return theX;
}

void doFit(TH1D * h,double leftRatio, double rightRatio, double leftEnd, double rightEnd){
    int bmax = h->GetMaximumBin();
    double max = h->GetBinContent(bmax)*leftRatio;
    int binl = bmax-1;
    for (;binl>=3; binl--){
        double height3bins = h->GetBinContent(binl);
        height3bins+=h->GetBinContent(binl-1);
        height3bins+=h->GetBinContent(binl-2);
        if (height3bins/3<max) break;
    }
    binl-=1;
    max = h->GetBinContent(bmax)*rightRatio;
    int binr = bmax+1;
    for (;binr<=h->GetNbinsX()-3; binr++){
        double height3bins = h->GetBinContent(binr);
        height3bins+=h->GetBinContent(binr+1);
        height3bins+=h->GetBinContent(binr+2);
        if (height3bins/3<max) break;
    }
    binr+=1;
    double left = h->GetBinCenter(binl);
    double right = h->GetBinCenter(binr);
    if (left<leftEnd) left = leftEnd;
    if (right>rightEnd) right = rightEnd;
    h->Fit("fres","QG","",left,right);
}

int getHitType(int type,bool isRight){
    // Original type is in decimal: [IMASTR]. I: peak index (only counting peaks over m_sumCut); M: peak index in a packet; A: smaller than aa cut? S: smaller than sum cut? T: 0 good, 3/6/7 very bad, 1/2/4/5 partly bad; R: 0 center, 1 left, 2 right, 3 guard, 4 dummy
    int ttype = (type/10)%10;
    if (isRight){
        if (ttype==1||ttype==4) type-=ttype*10; // l- or l+
    }
    else{
        if (ttype==2||ttype==5) type-=ttype*10; // r- or r+
    }
    return type;
}

int GetCandidate(TString & candSelBy, std::vector<int> * layerID, std::vector<int> * type, std::vector<double> * fitD[NCAND], std::vector<int> * sel[NCAND], int * nHitsS, double * chi2, double * chi2a){
    int cand = 0;
    if (candSelBy=="Original"){
        cand = 0;
    }
    else if (candSelBy=="LeastLatePeak"){
        int nLateHitsMin = 1e9;
        for (int iCand = 0; iCand<NCAND; iCand++){
            int nLateHits = 0;
            for (int ihit = 0; ihit<layerID->size(); ihit++){
                int ip = 0;
                for (int jhit = ihit-1; jhit>0; jhit--){
                    if ((*layerID)[jhit]!=(*layerID)[ihit]) break;
                    int ttype = getHitType((*type)[jhit],(*fitD[iCand])[jhit]>=0);
                    if (ttype<100) ip++;
                }
                if ((*sel[iCand])[ihit]==1){
                    if(ip!=0)
                        nLateHits++;
                }
            }
            if (nLateHits<nLateHitsMin){
                nLateHitsMin = nLateHits;
                cand = iCand;
            }
        }
    }
    else if (candSelBy=="FittingChi2"||candSelBy=="GlobalChi2"){
        double minchi2 = 1e9;
        int minNhitsS = 0;
        for (int iCand = 0; iCand<NCAND; iCand++){
            if (candSelBy=="GlobalChi2"){
                if ((minchi2>chi2a[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                    cand = iCand;
                    minchi2 = chi2a[iCand];
                    minNhitsS = nHitsS[iCand];
                }
            }
            else if (candSelBy=="FittingChi2"){
                if ((minchi2>chi2[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                    cand = iCand;
                    minchi2 = chi2[iCand];
                    minNhitsS = nHitsS[iCand];
                }
            }
        }
    }
    else{
        cand = 0;
    }
    return cand;
}

bool isInTriggerCounter(int geoSetup, double tinz, double tslz){
    bool isIn = true;
    if (geoSetup==1){
        if (fabs(tinz)>24) isIn = false;
    }
    else{
        if (fabs(tslz)>0.11) isIn = false;
    }
    return isIn;
}

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        MyRuntimeParameters::Get().ReadParamOverrideFile(configureFile);
        if (MyRuntimeParameters::Get().HasParameter("geoSetup")) m_geoSetup = MyRuntimeParameters::Get().GetParameterI("geoSetup");
        if (MyRuntimeParameters::Get().HasParameter("inputType")) m_inputType = MyRuntimeParameters::Get().GetParameterI("inputType");
        if (MyRuntimeParameters::Get().HasParameter("xtType")) m_xtType = MyRuntimeParameters::Get().GetParameterI("xtType");
        if (MyRuntimeParameters::Get().HasParameter("ana.AsymXT")) m_AsymXT = MyRuntimeParameters::Get().GetParameterI("ana.AsymXT");
        if (MyRuntimeParameters::Get().HasParameter("ana.CandSelBy")) m_CandSelBy = MyRuntimeParameters::Get().GetParameterS("ana.CandSelBy");
        if (MyRuntimeParameters::Get().HasParameter("ana.ClosestPeak")) m_ClosestPeak = MyRuntimeParameters::Get().GetParameterI("ana.ClosestPeak");
        //for cutting
        if (MyRuntimeParameters::Get().HasParameter("ana.nHitsMax")) m_nHitsMax = MyRuntimeParameters::Get().GetParameterI("ana.nHitsMax");
        if (MyRuntimeParameters::Get().HasParameter("ana.nHitsSmin")) m_nHitsSmin = MyRuntimeParameters::Get().GetParameterI("ana.nHitsSmin");
        if (MyRuntimeParameters::Get().HasParameter("ana.sumCut")) m_sumCut = MyRuntimeParameters::Get().GetParameterI("ana.sumCut");
        if (MyRuntimeParameters::Get().HasParameter("ana.aaCut")) m_aaCut = MyRuntimeParameters::Get().GetParameterI("ana.aaCut");
        if (MyRuntimeParameters::Get().HasParameter("ana.maxchi2")) m_maxchi2 = MyRuntimeParameters::Get().GetParameterD("ana.maxchi2");
        if (MyRuntimeParameters::Get().HasParameter("ana.maxslz")) m_maxslz = MyRuntimeParameters::Get().GetParameterD("ana.maxslz");
        if (MyRuntimeParameters::Get().HasParameter("ana.maxFD")) m_maxFD = MyRuntimeParameters::Get().GetParameterD("ana.maxFD");
        if (MyRuntimeParameters::Get().HasParameter("ana.tmaxSet")) m_tmaxSet = MyRuntimeParameters::Get().GetParameterI("ana.tmaxSet");
        //for binning
        if (MyRuntimeParameters::Get().HasParameter("ana.tmin")) m_tmin = MyRuntimeParameters::Get().GetParameterD("ana.tmin");
        if (MyRuntimeParameters::Get().HasParameter("ana.tmax")) m_tmax = MyRuntimeParameters::Get().GetParameterD("ana.tmax");
        if (MyRuntimeParameters::Get().HasParameter("ana.xmax")) m_xmax = MyRuntimeParameters::Get().GetParameterD("ana.xmax");
        if (MyRuntimeParameters::Get().HasParameter("ana.NbinT")) m_NbinT = MyRuntimeParameters::Get().GetParameterI("ana.NbinT");
        if (MyRuntimeParameters::Get().HasParameter("ana.NbinX")) m_NbinX = MyRuntimeParameters::Get().GetParameterI("ana.NbinX");
        if (MyRuntimeParameters::Get().HasParameter("ana.NbinRes")) m_NbinRes = MyRuntimeParameters::Get().GetParameterI("ana.NbinRes");
        if (MyRuntimeParameters::Get().HasParameter("ana.minchi2p")) m_minchi2p = MyRuntimeParameters::Get().GetParameterD("ana.minchi2p");
        if (MyRuntimeParameters::Get().HasParameter("ana.maxRes")) m_maxRes = MyRuntimeParameters::Get().GetParameterD("ana.maxRes");
        //for wire position calibration
        if (MyRuntimeParameters::Get().HasParameter("ana.scale")) m_scale = MyRuntimeParameters::Get().GetParameterD("ana.scale");
        if (MyRuntimeParameters::Get().HasParameter("ana.stepSize")) m_stepSize = MyRuntimeParameters::Get().GetParameterD("ana.stepSize");
        if (MyRuntimeParameters::Get().HasParameter("ana.minDeltaSlz")) m_minDeltaSlz = MyRuntimeParameters::Get().GetParameterD("ana.minDeltaSlz");
        if (MyRuntimeParameters::Get().HasParameter("ana.maxDeltaSlz")) m_maxDeltaSlz = MyRuntimeParameters::Get().GetParameterD("ana.maxDeltaSlz");
        if (MyRuntimeParameters::Get().HasParameter("ana.minDeltaInx")) m_minDeltaInx = MyRuntimeParameters::Get().GetParameterD("ana.minDeltaInx");
        if (MyRuntimeParameters::Get().HasParameter("ana.maxDeltaInx")) m_maxDeltaInx = MyRuntimeParameters::Get().GetParameterD("ana.maxDeltaInx");
    }
}

void print_usage(char * prog_name){
    fprintf(stderr,"Usage %s [options] prerunname runname [wires to be calibrated (all)]\n",prog_name);
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
    fprintf(stderr,"\t -R <run>\n");
    fprintf(stderr,"\t\t Run number set to run\n");
    fprintf(stderr,"\t -B <n>\n");
    fprintf(stderr,"\t\t Starting entry index set to n\n");
    fprintf(stderr,"\t -E <n>\n");
    fprintf(stderr,"\t\t Stopping entry index set to n\n");
    fprintf(stderr,"\t -L <l>\n");
    fprintf(stderr,"\t\t Default layer set to l\n");
    fprintf(stderr,"\t -H <h>\n");
    fprintf(stderr,"\t\t Histogram saving level set to h\n");
    fprintf(stderr,"\t -T\n");
    fprintf(stderr,"\t\t Disable in-trigger cut\n");
    fprintf(stderr,"\t -G\n");
    fprintf(stderr,"\t\t In getting offset and new XT, require all hits used in tracking are golden hits\n");
    fprintf(stderr,"\t -F <suffix>\n");
    fprintf(stderr,"\t\t Adding suffix to output files\n");
    fprintf(stderr,"\t -W\n");
    fprintf(stderr,"\t\t (false) Update wire map\n");
    fprintf(stderr,"\t -X <xt file>\n");
    fprintf(stderr,"\t\t Use given xt file instead of generating a new one\n");
    fprintf(stderr,"\t -n <n>\n");
    fprintf(stderr,"\t\t Maximum number of hits cut set to n\n");
    fprintf(stderr,"\t -f <f>\n");
    fprintf(stderr,"\t\t Minimum number of selected hits cut set to f\n");
    fprintf(stderr,"\t -c <c>\n");
    fprintf(stderr,"\t\t Maximum chi2 cut set to c\n");
    fprintf(stderr,"\t -v <v>\n");
    fprintf(stderr,"\t\t Minimum p-value cut set to v\n");
    fprintf(stderr,"\t -r <r>\n");
    fprintf(stderr,"\t\t Maximum resolution cut set to r\n");
    fprintf(stderr,"\t -z <z>\n");
    fprintf(stderr,"\t\t Maximum y-z slope cut set to z\n");
    fprintf(stderr,"\t -d <d>\n");
    fprintf(stderr,"\t\t Maximum fitD cut set to d\n");
    fprintf(stderr,"\t -o <o>\n");
    fprintf(stderr,"\t\t Maximum time range set to o\n");
    fprintf(stderr,"\t -s <s>\n");
    fprintf(stderr,"\t\t ADC sum over peak cut set to s\n");
    fprintf(stderr,"\t -a <a>\n");
    fprintf(stderr,"\t\t ADC sum over all cut set to a\n");
    fprintf(stderr,"\t -l <l>\n");
    fprintf(stderr,"\t\t Minimum time on axis set to l\n");
    fprintf(stderr,"\t -u <u>\n");
    fprintf(stderr,"\t\t Maximum time on axis set to u\n");
    fprintf(stderr,"\t -t <t>\n");
    fprintf(stderr,"\t\t Number of bins on time axis set to t\n");
    fprintf(stderr,"\t -m <m>\n");
    fprintf(stderr,"\t\t Number of bins on space axis set to m\n");
    fprintf(stderr,"\t -y <y>\n");
    fprintf(stderr,"\t\t Number of bins on resolution axis set to y\n");
    fprintf(stderr,"\t -g <g>\n");
    fprintf(stderr,"\t\t Geometry setup set to g\n");
    fprintf(stderr,"\t\t (0): normal; 1: finger\n");
    fprintf(stderr,"\t -i <i>\n");
    fprintf(stderr,"\t\t Input type set to i\n");
    fprintf(stderr,"\t\t (0) for data; 1 for MC\n");
    fprintf(stderr,"\t -x <x>\n");
    fprintf(stderr,"\t\t xt type set to x\n");
    fprintf(stderr,"\t\t XYZ (055) means polX for center, polY for middle and polZ for tail. If X is 0 then let middle function fit the center region.\n");
    fprintf(stderr,"\t -A\n");
    fprintf(stderr,"\t\t Use asymmetric XT\n");
    fprintf(stderr,"\t -S\n");
    fprintf(stderr,"\t\t Select candidate by:\n");
    fprintf(stderr,"\t\t ((O)riginal): the first one given by tracking (global chi2)\n");
    fprintf(stderr,"\t\t (G)lobalChi2: find the candidate with the smallest global chi2 including the test layer\n");
    fprintf(stderr,"\t\t (F)ittingChi2: find the candidate with the smallest fitting chi2\n");
    fprintf(stderr,"\t\t (L)eastLatePeak: find the candidate with the least late arrival peaks\n");
    fprintf(stderr,"\t -P\n");
    fprintf(stderr,"\t\t Use the peak with smallest residual to get driftT. By default just use the first one\n");
}
