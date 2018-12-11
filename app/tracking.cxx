#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <stdio.h>  /* printf, fgets */
#include <unistd.h> /* getopt */
#include <stdlib.h> /* atoi, atof */
#include <iostream> /* cout */

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TVector3.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TMath.h"

#include "MyProcessManager.hxx"
#include "MyRuntimeParameters.hxx"
#include "Log.hxx"
#include "header.hxx"

#define MAXPICK 8

//===================Configure Parameter============================
int m_runNo = 0;
TString m_prerunname = "prerun";
TString m_runname = "currun";
int m_iEntryStart = 0;
int m_iEntryStop = 0;
int m_modulo = 100;
bool m_memdebug = false;
int m_testlayer = 4;
TString m_configureFile = "";
int m_nHitsMax = 30;
int m_t0shift0 = 0;
int m_t0shift1 = 0;
int m_tmin = -20;
int m_tmax = 800;
double m_sumCut = -10;
double m_aaCut = 0;
int m_geoSetup = 0; // 0: normal; 1: finger
int m_inputType = 0; // 1 for MC; 0 for data
int m_peakType = 0; // 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
int m_workType = 0; // fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers

//===================Chamber Parameter============================
// map for wire position
double  map_x[NLAY][NCEL][2];
double  map_y[NLAY][NCEL][2];
double  map_z[NLAY][NCEL][2];
int     map_ch[NLAY][NCEL];
int     map_bid[NLAY][NCEL];
double  map_theta[NLAY][NCEL];
int     map_lid[NBRD][NCHS];
int     map_wid[NBRD][NCHS];
// mcp for cross points
double mcp_xc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
double mcp_zc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
// for resolution
double errord[NLAY][NCEL];

//==================About Scintillator======================
// normal scintillator
double sciYup = 0;
double sciYdown = 0;
double sciHL = 0;
double sciHW = 0;

//==================About Beam======================
double beamSlzMax = 0;
double beamSlxMax = 0;
double beamInzMax = 0;
double beamInxMax = 0;

//===================Input Output============================
int i_nHits;
std::vector<int> * i_layerID = 0;
std::vector<int> * i_wireID = 0;
std::vector<int> * i_type = 0;
std::vector<double> * i_driftT = 0;
std::vector<double> * i_driftDmc = 0;
std::vector<double> * i_driftD = 0;
std::vector<double> * o_dxl = 0;
std::vector<double> * o_dxr = 0;
int o_nFind = 0;
int o_nFit = 0;
// for finding
std::vector<double> * o_driftD[NCAND] = {0};
int o_icombi[NCAND];
int o_iselec[NCAND];
int o_npairs[NCAND];
double o_islx[NCAND];
double o_iinx[NCAND];
double o_islz[NCAND];
double o_iinz[NCAND];
double o_chi2x[NCAND];
double o_chi2z[NCAND];
double o_chi2i[NCAND];
double o_chi2pi[NCAND];
double o_chi2ai[NCAND];
std::vector<double> * o_calD[NCAND] = {0};
// for fitting
std::vector<int> * o_sel[NCAND] = {0};
std::vector<double> * o_fitD[NCAND] = {0};
int o_nHitsS[NCAND];
double o_chi2[NCAND];
double o_chi2p[NCAND];
double o_chi2a[NCAND];
double o_slx[NCAND];
double o_inx[NCAND];
double o_slz[NCAND];
double o_inz[NCAND];
double o_chi2mc[NCAND];
double o_chi2pmc[NCAND];
double o_chi2amc[NCAND];
double i_inxmc;
double i_inzmc;
double i_slxmc;
double i_slzmc;

//===================About tracking============================
// for track finding
TF1 * f_x = new TF1("f_x","pol1",sciYdown,sciYup); // x VS y
TGraphErrors * g_x = 0; // x VS y
TF1 * f_z = new TF1("f_z","pol1",sciYdown,sciYup); // z VS y
TGraphErrors * g_z = 0; // z VS y
double chi2z = 0;
double iinz = 0;
double islz = 0;
double chi2x = 0;
double iinx = 0;
double islx = 0;
double chi2i = 0;
double chi2pi = 0;
double chi2ai = 0;
double chi2mc = 0;
double chi2pmc = 0;
double chi2amc = 0;
std::vector<double> * t_calD = new std::vector<double>;
std::vector<double> * t_fitD = new std::vector<double>;
std::vector<double> * t_driftD = new std::vector<double>;
std::vector<int> * t_sel = new std::vector<int>;
std::vector<int> * t_lr = new std::vector<int>;
// for track fitting
TMinuit *gMinuit = 0;
double arglist[10];
int ierflg = 0;
double amin,edm,errdef;
int nvpar,nparx,icstat;
double chi2 = 0;
double chi2p = 0;
double chi2a = 0;
double inz = 0;
double slz = 0;
double inx = 0;
double slx = 0;
// for get dist
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

//===================Hit Handlers============================
std::vector<std::vector<int> > v_layer_ihit; // vectors of indices of hits in each layer
std::vector<int> v_pick_lid; // layer index of each picked layer (for pair check)
std::vector<int> pick_ihit; // hit index of each picked hit (for pair check)
std::vector<double> pick_wy; // y position of each picked hit (for pair check)
std::vector<double> pair_wx; // x position of each picked hit pair (center position)
std::vector<double> pair_wy; // y position of each picked hit pair (center position)
std::vector<double> pair_wz; // z position of each picked hit pair (center position)

//===================Functions============================
void print_usage(char* prog_name);
double t2x(double time, int lid, int wid, int lr, int & status);
int Tracking(int ipick,int & iselection,int iEntry=0); //pick up one hit from each layer, and iterate in all combinations including left/right ambiguity
int doFitting(int nPicks,int iEntry=0,int iselection = 0);
void setLRdriftD(int nPicks,int icombi);
int updatePairPositions(int nPicks,int & nPairs);
int updateHitPositions(int nPicks);
int setErrors(int nPairs, bool noError = false);
int getChi2XZ(int nPairs, double & chi2x, double & chi2z);
int fityx(int nPairs);
int fityz(int nPairs);
bool checkScintillator(double saftyFactor,double inx, double slx, double inz, double slz);
bool checkChi2(int nHitsSel,int nPairs,int icombi, int iselection);
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz);
void getchi2(double &f, double & cp, double & ca, double slx, double inx, double slz, double inz,bool all = false);
void fcn(int &npar, double *gin, double &f, double *par, int iflag);
void do_fit(double slix, double inix,double sliz, double iniz);
int getHitIndex(int lid, int nHits);
int getHitType(int type,bool isRight);
bool isSame(int iCand);
double getError(int lid,double dt, bool isR);

MyProcessManager * pMyProcessManager;

//===================About xt============================
TF1 * f_left[NLAY+2];
TF1 * f_right[NLAY+2];
// for error function
TF1 * funcErr;
TGraph * gr_error_left[NLAY+2];
TGraph * gr_error_right[NLAY+2];

int main(int argc, char** argv){
    int temp_nHitsMax = 0; bool set_nHitsMax = false;
    int temp_t0shift0 = 0; bool set_t0shift0 = false;
    int temp_t0shift1 = 0; bool set_t0shift1 = false;
    int temp_tmin = 0; bool set_tmin = false;
    int temp_tmax = 0; bool set_tmax = false;
    double temp_sumCut = 0; bool set_sumCut = false;
    double temp_aaCut = 0; bool set_aaCut = false;
    int temp_geoSetup = 0; bool set_geoSetup = false;
    int temp_inputType = 0; bool set_inputType = false;
    int temp_peakType = 0; bool set_peakType = false;
    int temp_workType = 0; bool set_workType = false;

    std::map<std::string, Log::ErrorPriority> namedDebugLevel;
    std::map<std::string, Log::LogPriority> namedLogLevel;
    int    opt_result;
	while((opt_result=getopt(argc,argv,"M:R:B:E:L:C:n:x:y:l:u:s:a:g:i:p:w:D:V:"))!=-1){
		switch(opt_result){
			/* INPUTS */
			case 'M':
			    m_modulo = atoi(optarg);
                printf("Printing modulo set to %d\n",m_modulo);
			case 'R':
			    m_runNo = atoi(optarg);
                printf("Run number set to %d\n",m_runNo);
			case 'B':
			    m_iEntryStart = atoi(optarg);
                printf("Starting entry index set to %d\n",m_iEntryStart);
			case 'E':
			    m_iEntryStop = atoi(optarg);
                printf("Stopping entry index set to %d\n",m_iEntryStop);
			case 'L':
			    m_testlayer = atoi(optarg);
                printf("Test layer set to %d\n",m_testlayer);
			case 'C':
				m_configureFile = optarg;
                printf("Using configure file \"%s\"\n",optarg);
                break;
			case 'n':
			    temp_nHitsMax = atoi(optarg);set_nHitsMax = true;
                printf("Maximum number of hits cut set to %d\n",temp_nHitsMax);
			case 'x':
			    temp_t0shift0 = atoi(optarg);set_t0shift0 = true;
                printf("T0 shift on board 0 set to %d\n",temp_t0shift0);
			case 'y':
			    temp_t0shift1 = atoi(optarg);set_t0shift1 = true;
                printf("T0 shift on board 1 set to %d\n",temp_t0shift1);
			case 'l':
			    temp_tmin = atoi(optarg);set_tmin = true;
                printf("Minimum time cut set to %d\n",temp_tmin);
			case 'u':
			    temp_tmax = atoi(optarg);set_tmax = true;
                printf("Maximum time cut set to %d\n",temp_tmax);
			case 's':
			    temp_sumCut = atoi(optarg);set_sumCut = true;
                printf("ADC sum over peak cut set to %d\n",temp_sumCut);
			case 'a':
			    temp_aaCut = atoi(optarg);set_aaCut = true;
                printf("ADC sum over all cut set to %d\n",temp_aaCut);
			case 'g':
			    temp_geoSetup = atoi(optarg);set_geoSetup = true;
                printf("Geometry setup set to %d\n",temp_geoSetup);
			case 'i':
			    temp_inputType = atoi(optarg);set_inputType = true;
                printf("Input type set to %d\n",temp_inputType);
			case 'p':
			    temp_peakType = atoi(optarg);set_peakType = true;
                printf("Peak type set to %d\n",temp_peakType);
			case 'w':
			    temp_workType = atoi(optarg);set_workType = true;
                printf("Work type set to %d\n",temp_workType);
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
        m_inputType = MyRuntimeParameters::Get().GetParameterI("inputType");
        m_geoSetup = MyRuntimeParameters::Get().GetParameterI("geoSetup");
        m_peakType = MyRuntimeParameters::Get().GetParameterI("peakType");
        m_nHitsMax = MyRuntimeParameters::Get().GetParameterI("tracking.nHitsMax");
        m_t0shift0 = MyRuntimeParameters::Get().GetParameterI("tracking.t0shift0");
        m_t0shift1 = MyRuntimeParameters::Get().GetParameterI("tracking.t0shift1");
        m_tmin = MyRuntimeParameters::Get().GetParameterI("tracking.tmin");
        m_tmax = MyRuntimeParameters::Get().GetParameterI("tracking.tmax");
        m_sumCut = MyRuntimeParameters::Get().GetParameterD("tracking.sumCut");
        m_aaCut = MyRuntimeParameters::Get().GetParameterD("tracking.aaCut");
        m_workType = MyRuntimeParameters::Get().GetParameterI("tracking.workType");;
    }
    if (set_geoSetup) m_geoSetup = temp_geoSetup;
    if (set_nHitsMax) m_nHitsMax = temp_nHitsMax;
    if (set_t0shift0) m_t0shift0 = temp_t0shift0;
    if (set_t0shift1) m_t0shift1 = temp_t0shift1;
    if (set_tmin) m_tmin = temp_tmin;
    if (set_tmax) m_tmax = temp_tmax;
    if (set_sumCut) m_sumCut = temp_sumCut;
    if (set_aaCut) m_aaCut = temp_aaCut;
    if (set_workType) m_workType = temp_workType;
    if (set_inputType) m_inputType = temp_inputType;
    if (set_peakType) m_peakType = temp_peakType;

	if (argc-optind<2){
	    print_usage(argv[0]);
		return -1;
    }
    m_prerunname = argv[optind++];
    m_runname= argv[optind++];

    printf("##############%s##################\n",argv[0]);
    printf("runNo       = %d\n",m_runNo);
    printf("prerunname  = \"%s\"\n",m_prerunname.Data());
    printf("runname     = \"%s\"\n",m_runname.Data());
    printf("test layer  = %d\n",m_testlayer);
    printf("nHitsMax    = %d\n",m_nHitsMax);
    printf("t0shift b0  = %d\n",m_t0shift0);
    printf("t0shift b1  = %d\n",m_t0shift1);
    printf("tmin        = %d\n",m_tmin);
    printf("tmax        = %d\n",m_tmax);
    printf("sumCut      = %f\n",m_sumCut);
    printf("aaCut       = %f\n",m_aaCut);
    printf("geoSetup:     %s\n",m_geoSetup==0?"normal scintillator":"finger scintillator");
    printf("workType    = %d, %s\n",m_workType,m_workType==0?"all as 0":(m_workType==1?"even/odd":(m_workType==-1?"even/odd reversed":"all layers")));
    printf("inputType   = %d, %s\n",m_inputType,m_inputType==0?"Real Data":(m_inputType==2?"MC X":"MC T"));
    printf("peakType    = %d, %s\n",m_peakType,m_peakType==0?"First peak over threshold":"All peaks over threshold");
    printf("Start Entry = %d\n",m_iEntryStart);
    printf("Stop Entry  = %d\n",m_iEntryStop);

    TString HOME=getenv("CDCS8WORKING_DIR");

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

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
    printf("##############Geometry##################\n");
    printf("sciYup      = %.3e\n",sciYup);
    printf("sciYdown    = %.3e\n",sciYdown);
    printf("sciHL       = %.3e\n",sciHL);
    printf("sciHW       = %.3e\n",sciHW);

	//===================Set beam property============================
    // FIXME: currently set a broader range. Need further investigation
    beamSlzMax = 0.3;
    beamSlxMax = 0.1;
    beamInzMax = sciHL*1.5; // mm
    beamInxMax = sciHW*1.5;
    printf("##############Beam##################\n");
    printf("beamSlzMax  = %.3e\n",beamSlzMax);
    printf("beamSlxMax  = %.3e\n",beamSlxMax);
    printf("beamInzMax  = %.3e\n",beamInzMax);
    printf("beamInxMax  = %.3e\n",beamInxMax);

	//===================Get error============================
	funcErr = new TF1("funcErr","0.346904-0.221775*x+0.080226*x*x-0.0128037*x*x*x+0.000755738*x*x*x*x",0,10);

    //===================Prepare Maps============================
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            map_x[lid][wid][0] = 0;
            map_y[lid][wid][0] = 0;
            map_z[lid][wid][0] = 0;
            map_x[lid][wid][1] = 0;
            map_y[lid][wid][1] = 0;
            map_z[lid][wid][1] = 0;
        	map_ch[lid][wid] = -1;
        	map_bid[lid][wid] = -1;
            if (lid <NZXP){ // z-x planes corresponding to the layerID of the lower layer counting from 1 
                for (int wjd = 0; wjd<NCEL; wjd++){
                    mcp_xc[lid][wid][wjd] = 999;
                    mcp_zc[lid][wid][wjd] = 999;
                }
            }
        }
    }

    //===================Get Wire Position============================
    TFile * TFile_wirepos = new TFile(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),m_runNo,m_prerunname.Data()));
    TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
    int     wp_bid;
    int     wp_ch;
    int     wp_wid;
    int     wp_lid;
    double  wp_xro;
    double  wp_yro;
    double  wp_xhv;
    double  wp_yhv;
    TTree_wirepos->SetBranchAddress("b",&wp_bid);
    TTree_wirepos->SetBranchAddress("ch",&wp_ch);
    TTree_wirepos->SetBranchAddress("l",&wp_lid);
    TTree_wirepos->SetBranchAddress("w",&wp_wid);
    TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
    TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
    TTree_wirepos->SetBranchAddress("xro",&wp_xro);
    TTree_wirepos->SetBranchAddress("yro",&wp_yro);
    for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
        TTree_wirepos->GetEntry(i);
        if (wp_lid>=0&&wp_lid<NLAY&&wp_wid>=0&&wp_wid<NCEL){
            map_x[wp_lid][wp_wid][0] = wp_xhv;
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_z[wp_lid][wp_wid][0] = -chamberHL;
            map_x[wp_lid][wp_wid][1] = wp_xro;
            map_y[wp_lid][wp_wid][1] = wp_yro;
            map_z[wp_lid][wp_wid][1] = chamberHL;
            map_ch[wp_lid][wp_wid] = wp_ch;
            map_bid[wp_lid][wp_wid] = wp_bid;
            map_theta[wp_lid][wp_wid] = atan(-(wp_xhv-wp_xro)/chamberHL/2); // rotation angle w.r.t the dart plane: read out  plane; positive rotation angle point to -x direction
            MyNamedInfo("WireMap",Form("map_theta[%d][%d] = atan(-(%.3e-%.3e)/%.3e/2) = %.3e\n",wp_lid,wp_wid,wp_xhv,wp_xro,chamberHL,map_theta[wp_lid][wp_wid]));
        }
        else{
            fprintf(stderr,"WARNING: Entry %d in wiremap file, lid = %d wid = %d out of range (%d,%d)!\n",i,wp_lid,wp_wid,NLAY,NCEL);
        }
        if (wp_bid>=0&&wp_bid<NBRD&&wp_ch>=0&&wp_ch<NCHS){
            map_lid[wp_bid][wp_ch] = wp_lid;
            map_wid[wp_bid][wp_ch] = wp_wid;
        }
        else{
            fprintf(stderr,"WARNING: Entry %d in wiremap file, bid = %d ch = %d out of range (%d,%d)!\n",i,wp_bid,wp_ch,NBRD,NCHS);
        }
    }
    TFile_wirepos->Close();

    //==================Get Crosspoints==========================
    TFile * TFile_crosspoint = new TFile(HOME+"/Input/crosspoint.root");
    TTree * TTree_crosspoint = (TTree*) TFile_crosspoint->Get("t");
    int     cp_l1;
    int     cp_l2;
    int     cp_w1;
    int     cp_w2;
    double  cp_zc;
    double  cp_xc;
    TTree_crosspoint->SetBranchAddress("l1",&cp_l1);
    TTree_crosspoint->SetBranchAddress("l2",&cp_l2);
    TTree_crosspoint->SetBranchAddress("w1",&cp_w1);
    TTree_crosspoint->SetBranchAddress("w2",&cp_w2);
    TTree_crosspoint->SetBranchAddress("z",&cp_zc);
    TTree_crosspoint->SetBranchAddress("x",&cp_xc);
    int nEntries_crosspoint = TTree_crosspoint->GetEntries();
    for (int iEntry = 0; iEntry<nEntries_crosspoint; iEntry++){
        TTree_crosspoint->GetEntry(iEntry);
        if (cp_l1>=0&&cp_l1<NLAY&&cp_w1>=0&&cp_w1<NCEL&&cp_l2>=0&&cp_l2<NLAY&&cp_w2>=0&&cp_w2<NCEL){
            mcp_xc[cp_l1][cp_w1][cp_w2] = cp_xc;
            mcp_zc[cp_l1][cp_w1][cp_w2] = cp_zc;
        }
    }
    TFile_crosspoint->Close();

    //===================Prepare XT curves==============================
    printf("##############XT##################\n");
    printf("Reading from %s/info/xt.%d.%s.root\n",HOME.Data(),m_runNo,m_prerunname.Data());
    TFile * i_xt = new TFile(HOME+Form("/info/xt.%d.",m_runNo)+m_prerunname+".root");
    for (int i = 0; i<NLAY; i++){
        f_left[i] = (TF1*) i_xt->Get(Form("fl_%d",i));
        f_right[i] = (TF1*) i_xt->Get(Form("fr_%d",i));
		gr_error_left[i] = (TGraph*)i_xt->Get(Form("gr_sigts_slicetl_%d",i));
		gr_error_right[i] = (TGraph*)i_xt->Get(Form("gr_sigts_slicetr_%d",i));
        double tmaxl = 0;
        double tmaxr = 0;
        double tminl = 0;
        double tminr = 0;
        double xmaxl = 0;
        double xmaxr = 0;
        double xminl = 0;
        double xminr = 0;
        if (f_left[i]){
            tmaxl = f_left[i]->GetXmax();
            tminl = f_left[i]->GetXmin();
            xmaxl = f_left[i]->Eval(tmaxl);
            xminl = f_left[i]->Eval(tminl);
        }
        if (f_right[i]){
            tmaxr = f_right[i]->GetXmax();
            tminr = f_right[i]->GetXmin();
            xmaxr = f_right[i]->Eval(tmaxr);
            xminr = f_right[i]->Eval(tminr);
        }
        if (f_left[i]||f_right[i]){
        	printf("  XT in layer[%d]: (%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)-(%.3e,%.3e)\n",i,tmaxl,xmaxl,tminl,xminl,tminr,xminr,tmaxr,xmaxr);
			if (!gr_error_left[i]||!gr_error_right[i]) fprintf(stderr,"Cannot find gr_error_l/r[%d]! Would assume default error 0.2 mm\n",i);
		}
    }
	f_left[NLAY] = (TF1*) i_xt->Get("fl_even");
	f_right[NLAY] = (TF1*) i_xt->Get("fr_even");
	gr_error_left[NLAY] = (TGraph*)i_xt->Get("gr_sigts_slicetl_even");
	gr_error_right[NLAY] = (TGraph*)i_xt->Get("gr_sigts_slicetr_even");
	f_left[NLAY+1] = (TF1*) i_xt->Get("fl_odd");
	f_right[NLAY+1] = (TF1*) i_xt->Get("fr_odd");
	gr_error_left[NLAY+1] = (TGraph*)i_xt->Get("gr_sigts_slicetl_odd");
	gr_error_right[NLAY+1] = (TGraph*)i_xt->Get("gr_sigts_slicetr_odd");

    //===================Get input ROOT file============================
    TChain * c = new TChain("t","t");
    if (m_inputType==3)
        c->Add(HOME+Form("/root/h_%d.layer%d.MC.root",m_runNo,m_testlayer));
    else if (m_inputType==2)
        c->Add(HOME+Form("/root/h_%d.MC.root",m_runNo));
    else
        c->Add(HOME+Form("/root/h_%d.root",m_runNo));
    int triggerNumber;
    std::vector<int> * i_np = 0;
    std::vector<int> * i_ip = 0;
    std::vector<int> * i_clk = 0;
    std::vector<int> * i_width = 0;
    std::vector<int> * i_peak = 0;
    std::vector<int> * i_height = 0;
    std::vector<int> * i_mpn = 0;
    std::vector<int> * i_mpi = 0;
    std::vector<int> * i_rank = 0;
    std::vector<double> * i_ped = 0;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_aa = 0;
    c->SetBranchAddress("triggerNumber",&triggerNumber);
    c->SetBranchAddress("nHits",&i_nHits);
    c->SetBranchAddress("driftT",&i_driftT);
    if (m_inputType) c->SetBranchAddress("driftDmc",&i_driftDmc);
    if (m_inputType==2||m_inputType==3) c->SetBranchAddress("driftD",&i_driftD);
    c->SetBranchAddress("layerID",&i_layerID);
    c->SetBranchAddress("wireID",&i_wireID);
    c->SetBranchAddress("type",&i_type); // 0 center, 1 left, 2 right, 3 guard, 4 dummy
    c->SetBranchAddress("np",&i_np);
    c->SetBranchAddress("ip",&i_ip);
    c->SetBranchAddress("clk",&i_clk);
    c->SetBranchAddress("width",&i_width);
    c->SetBranchAddress("peak",&i_peak);
    c->SetBranchAddress("height",&i_height);
    c->SetBranchAddress("mpn",&i_mpn);
    c->SetBranchAddress("mpi",&i_mpi);
    c->SetBranchAddress("rank",&i_rank);
    c->SetBranchAddress("ped",&i_ped);
    c->SetBranchAddress("sum",&i_sum);
    c->SetBranchAddress("aa",&i_aa);
    if (m_inputType){
		c->SetBranchAddress("inxmc",&i_inxmc);
		c->SetBranchAddress("inzmc",&i_inzmc);
		c->SetBranchAddress("slxmc",&i_slxmc);
		c->SetBranchAddress("slzmc",&i_slzmc);
    }

    //===================Prepare output ROOT file============================
    printf("Output file: %s/root/t_%d.%s.layer%d.root\n",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer);
    TFile * of = new TFile(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer),"RECREATE"); 
    TTree * ot = new TTree("t","t");
    // from h_XXX
    ot->Branch("triggerNumber",&triggerNumber);
    ot->Branch("driftT",&i_driftT);
    if (m_inputType) ot->Branch("driftDmc",&i_driftDmc);
    ot->Branch("nHits",&i_nHits);
    ot->Branch("layerID",&i_layerID);
    ot->Branch("wireID",&i_wireID);
    ot->Branch("type",&i_type); // in dec, [IMASTR]. I: peak index (only counting peaks over m_sumCut); M: peak index in a packet; A: smaller than aa cut? S: smaller than sum cut? T: -1 <m_tmin, 0 good, 1 >m_tmax; R: 0 center, 1 left, 2 right, 3 guard, 4 dummy
    ot->Branch("np",&i_np);
    ot->Branch("ip",&i_ip);
    ot->Branch("clk",&i_clk);
    ot->Branch("width",&i_width);
    ot->Branch("peak",&i_peak);
    ot->Branch("height",&i_height);
    ot->Branch("mpn",&i_mpn);
    ot->Branch("mpi",&i_mpi);
    ot->Branch("rank",&i_rank);
    ot->Branch("ped",&i_ped);
    ot->Branch("sum",&i_sum);
    ot->Branch("aa",&i_aa);
    // basic
    int nHitsG;
    ot->Branch("nHitsG",&nHitsG); // number of good hits in layers other than the test one: in t region and with good peak quality
    ot->Branch("dxl",&o_dxl);
    ot->Branch("dxr",&o_dxr);
    // track finding/fitting with different candidates;
    ot->Branch("nFit",&o_nFit);
    ot->Branch("nFind",&o_nFind);
    for (int iCand = 0; iCand<NCAND; iCand++){
        ot->Branch(Form("driftD%d",iCand),&(o_driftD[iCand]));
        ot->Branch(Form("icom%d",iCand),&(o_icombi[iCand]));
        ot->Branch(Form("isel%d",iCand),&(o_iselec[iCand]));
        ot->Branch(Form("npairs%d",iCand),&(o_npairs[iCand]));
        ot->Branch(Form("islx%d",iCand),&(o_islx[iCand]));
        ot->Branch(Form("iinx%d",iCand),&(o_iinx[iCand]));
        ot->Branch(Form("islz%d",iCand),&(o_islz[iCand]));
        ot->Branch(Form("iinz%d",iCand),&(o_iinz[iCand]));
        ot->Branch(Form("chi2x%d",iCand),&(o_chi2x[iCand]));
        ot->Branch(Form("chi2z%d",iCand),&(o_chi2z[iCand]));
        ot->Branch(Form("chi2i%d",iCand),&(o_chi2i[iCand]));
        ot->Branch(Form("chi2pi%d",iCand),&(o_chi2pi[iCand]));
        ot->Branch(Form("chi2ai%d",iCand),&(o_chi2ai[iCand]));
        ot->Branch(Form("sel%d",iCand),&(o_sel[iCand]));
        ot->Branch(Form("calD%d",iCand),&(o_calD[iCand]));
        ot->Branch(Form("slx%d",iCand),&(o_slx[iCand]));
        ot->Branch(Form("inx%d",iCand),&(o_inx[iCand]));
        ot->Branch(Form("slz%d",iCand),&(o_slz[iCand]));
        ot->Branch(Form("inz%d",iCand),&(o_inz[iCand]));
        ot->Branch(Form("chi2%d",iCand),&(o_chi2[iCand]));
        ot->Branch(Form("chi2p%d",iCand),&(o_chi2p[iCand]));
        ot->Branch(Form("chi2a%d",iCand),&(o_chi2a[iCand]));
        ot->Branch(Form("fitD%d",iCand),&(o_fitD[iCand]));
        ot->Branch(Form("nHitsS%d",iCand),&(o_nHitsS[iCand])); // number of hits selected from finding and fed to fitting
		if (m_inputType){
			ot->Branch(Form("chi2mc%d",iCand), &(o_chi2mc[iCand]));
			ot->Branch(Form("chi2pmc%d",iCand),&(o_chi2pmc[iCand]));
			ot->Branch(Form("chi2amc%d",iCand),&(o_chi2amc[iCand]));
		}
    }
	if (m_inputType){
		ot->Branch("inxmc",&i_inxmc);
		ot->Branch("inzmc",&i_inzmc);
		ot->Branch("slxmc",&i_slxmc);
		ot->Branch("slzmc",&i_slzmc);
	}
    // initialize vectors
    o_dxl = new std::vector<double>;
    o_dxr = new std::vector<double>;
    for(int iCand = 0; iCand<NCAND; iCand++){
        o_driftD[iCand] = new std::vector<double>;
        o_calD[iCand] = new std::vector<double>;
        o_fitD[iCand] = new std::vector<double>;
        o_sel[iCand] = new std::vector<int>;
    }

    //===================Efficiency Counter============================
    int N_trigger = 0;
    int N_found = 0;
    int N_good = 0;

    //===================Set Hit handlers============================
    for (int i = 0; i<NLAY; i++){
        std::vector<int> hits;
        v_layer_ihit.push_back(hits);
    }
    pick_ihit.resize(NLAY);
    pick_wy.resize(NLAY);
    pair_wx.resize(NLAY);
    pair_wy.resize(NLAY);
    pair_wz.resize(NLAY);
    g_x = new TGraphErrors(NLAY,&(pair_wy[0]),&(pair_wx[0]),0,0);
    g_z = new TGraphErrors(NLAY,&(pair_wy[0]),&(pair_wz[0]),0,0);

    //===================Tracking====================================
    Long64_t N = c->GetEntries();
    if (!m_iEntryStop&&!m_iEntryStart){m_iEntryStart = 0; m_iEntryStop=N-1;}
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    for (Long64_t iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
        MyNamedInfo("Tracking","############ Entry "<<iEntry<<" #############");
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        if (iEntry%m_modulo == 0){
            MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
            std::cout<<iEntry<<std::endl;
        }
        c->GetEntry(iEntry);
        N_trigger++; // triggered event

        // prepare
        o_nFind = 0;
        o_nFit = 0;
        o_dxl->resize(i_nHits);
        o_dxr->resize(i_nHits);
        for (int iCand = 0; iCand<NCAND; iCand++){
            o_islx[iCand] = 0;
            o_iinx[iCand] = 0;
            o_islz[iCand] = 0;
            o_iinz[iCand] = 0;
            o_chi2x[iCand] = 1e9;
            o_chi2z[iCand] = 1e9;
            o_chi2i[iCand] = 1e9;
            o_chi2pi[iCand] = 1e9;
            o_chi2ai[iCand] = 1e9;
            o_npairs[iCand] = 0;
            o_icombi[iCand] = 0;
            o_iselec[iCand] = 0;
            o_calD[iCand]->resize(i_nHits);
            o_fitD[iCand]->resize(i_nHits);
            o_sel[iCand]->resize(i_nHits);
            o_driftD[iCand]->resize(i_nHits);
            o_nHitsS[iCand] = 0;
            o_slx[iCand] = 0;
            o_inx[iCand] = 0;
            o_slz[iCand] = 0;
            o_inz[iCand] = 0;
            o_chi2[iCand] = 1e9;
            o_chi2p[iCand] = 1e9;
            o_chi2a[iCand] = 1e9;
        }
        t_calD->resize(i_nHits);
        t_fitD->resize(i_nHits);
        t_driftD->resize(i_nHits);
        t_sel->resize(i_nHits);
        t_lr->resize(i_nHits);
        nHitsG = 0;
        v_pick_lid.clear();
        for (int i = 0; i<NLAY; i++){
            v_layer_ihit[i].clear();
        }

        // get basical cdc hit information
        int prevch = -1; // the mych ID of the previous hit
        int npoc = 0; // number of peaks over sum cut in this chanel
        for (int ihit = 0; ihit<i_nHits; ihit++){
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int bid = map_bid[lid][wid];
            (*i_driftT)[ihit]+=(bid==0?m_t0shift0:m_t0shift1); // fix driftT according to t0shift
            double dt = (*i_driftT)[ihit];
            int mych = lid*1000+wid;
            if (mych!=prevch) // new channel
                npoc = 0; // reset npoc
            prevch = mych; // record current channel ID
            int statusl,statusr; // 1:  large t; -1: small t; 0: good t
            if (m_inputType!=2&&m_inputType!=3){
                (*o_dxl)[ihit] = t2x(dt,lid,wid,-1,statusl);
                (*o_dxr)[ihit] = t2x(dt,lid,wid,1,statusr);
            }
            else{
                (*o_dxl)[ihit] = -fabs((*i_driftD)[ihit]);
                (*o_dxr)[ihit] = fabs((*i_driftD)[ihit]);
            }
            int type = (*i_type)[ihit]; // IMASTR
            // R: region
            // keep the original defination
            // T: time
            if (m_inputType!=2&&m_inputType!=3){
                if (dt<m_tmin) type+=3*10;
                else if (dt>m_tmax) type+=6*10;
                else{
                    if (statusl==-1&&statusr==0) type+=1*10;
                    else if (statusl==0&&statusr==-1) type+=2*10;
                    else if (statusl==-1&&statusr==-1) type+=3*10;
                    else if (statusl==1&&statusr==0) type+=4*10;
                    else if (statusl==0&&statusr==1) type+=5*10;
                    else if (statusl==1&&statusr==1) type+=6*10;
                    else if (statusl==0&&statusr==0) type+=0;
                    else type+=7*10;
                }
            }
            type+=npoc*100000; // number of peaks above threshold before this peak in this channel
            // S: sum of wave packet
            if ((*i_sum)[ihit]<m_sumCut) type+=1*100;
            else npoc++; // over sum cut, then increment npoc
            // A: sum of full waveform
            if ((*i_aa)[ihit]<m_aaCut) type+=1*1000;
            // M: index of peak in a multi-peak wave packet
            type+=(*i_mpi)[ihit]*10000;
            // I: index of peak over sum cut in the whole waveform
            // added above // type+=npoc*100000; 
            (*i_type)[ihit] = type;
            int ttype = getHitType(type,true);
            if (lid != m_testlayer&&ttype<=3){ // good hit
                MyNamedVerbose("Tracking",Form("  Entry %d: dxl[%d][%d] = dxl[%d] = t2x(%.3e) = %.3e\n",iEntry,lid,wid,ihit,dt,(*o_dxl)[ihit]));
                MyNamedVerbose("Tracking",Form("  Entry %d: dxr[%d][%d] = dxr[%d] = t2x(%.3e) = %.3e\n",iEntry,lid,wid,ihit,dt,(*o_dxr)[ihit]));
                v_layer_ihit[lid].push_back(ihit);
                nHitsG++;
            }
        }

        // To find pair candidates
        int npairs = 0;
        int prelid = -1;
        for(int lid = 1; lid<NLAY-1; lid++){
            if (lid==m_testlayer||lid+1==m_testlayer) continue;
            if (v_layer_ihit[lid].size()>0&&v_layer_ihit[lid+1].size()>0){// both layers have hits
                if (prelid+1 != lid)v_pick_lid.push_back(lid);
                v_pick_lid.push_back(lid+1);
                prelid = lid;
                npairs++;
                if (v_pick_lid.size()>=MAXPICK) break;
            }
        }
        for (int ipick = 0; ipick<v_pick_lid.size(); ipick++){
            MyNamedVerbose("Tracking",Form(" pick layer %d\n",v_pick_lid[ipick]));
        }

        // Do tracking
        MyNamedVerbose("Tracking",Form("nHitsG = %d, npairs = %d\n",nHitsG,npairs));
        if (nHitsG<=m_nHitsMax&&npairs>=3){ // need at least 3 pairs to do the fitting; number of hits should be small to control the time cost
            N_found++;

            int nSelections = 0;
            Tracking(0,nSelections,iEntry); // 0 means starting from the 1st pick; nSelections is the number of possible choices by selecting one hit per layer;
            if (o_nHitsS[0]>=5){ // at least 5 hits to fit: NDF of 3-D track without field is 4
                N_good++;
            }
        }
        ot->Fill();
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    }// end of event loop
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    ot->Write();
    of->Close();
    printf("Triggered Events: %d\n",N_trigger);
    printf("Found Events: %d\n",N_found);
    printf("Good Events: %d\n",N_good);
    return 0;
}

int Tracking(int ipick,int & iselection,int iEntry){
    if (ipick == v_pick_lid.size()){ // finished picking hits
        MyNamedVerbose("Tracking",Form(" Finished picking selection %d:\n",iselection));
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        doFitting(v_pick_lid.size(),iEntry,iselection);
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        iselection++;
    }
    else{
        int lid = v_pick_lid[ipick];
        for (int i = 0; i<v_layer_ihit[lid].size(); i++){
            int ihit = v_layer_ihit[lid][i];
            int wid = (*i_wireID)[ihit];
            MyNamedVerbose("Tracking",Form(" => pick # %d, layer %d, wire %d, hit[%d], ihit = %d\n",ipick,lid,wid,i,ihit));
            pick_ihit[ipick] = v_layer_ihit[lid][i];
            Tracking(ipick+1,iselection,iEntry);
        }
    }
    return 0;
}

int doFitting(int nPicks,int iEntry,int iselection){
    int ncombi = pow(2,nPicks);
    MyNamedVerbose("Tracking",Form("  %d picked layers -> %d combinations\n",nPicks,ncombi));
    for (int icombi = 0; icombi<ncombi; icombi++){ // each combination corresponds to a unique left/right selection set
    	o_nFind++;
        MyNamedVerbose("Tracking",Form("     combi %d\n",icombi));
        t_lr->clear();
        t_lr->resize(i_nHits,0); // 0 is used as a default value to indicate that this hit is not picked, thus left/right unfixed
        setLRdriftD(nPicks,icombi); // for picked hits
        f_x->SetParameters(0,0);
        f_z->SetParameters(0,0);
        int nPairs = 0;

        updateHitPositions(nPicks); // fix wy positions
        int result = updatePairPositions(nPicks,nPairs);
        if (result) continue;
        setErrors(nPairs,true);
        fityz(nPairs);
        fityx(nPairs);
        iinz = f_z->Eval(sciYup);
        islz = f_z->GetParameter(1);
        iinx = f_x->Eval(sciYup);
        islx = f_x->GetParameter(1);
        bool inScint = checkScintillator(2.5,iinx,islx,iinz,islz); // FIXME: need to tune
        bool fromSource = islx>-beamSlxMax&&islx<beamSlxMax&&islz>-beamSlzMax&&islz<beamSlzMax;
        int nGood = getChi2XZ(nPairs,chi2x,chi2z);
        MyNamedVerbose("Tracking",Form("       1st RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",islx,sciYup,iinx,chi2x,islz,sciYup,iinz,chi2z));
        if (!fromSource||!inScint) {f_x->SetParameters(0,0); f_z->SetParameters(0,0);}
        
        updateHitPositions(nPicks); // fix wy positions
        result = updatePairPositions(nPicks,nPairs);
        if (result) continue;
        setErrors(nPairs,false);
        fityz(nPairs);
        fityx(nPairs);
        iinz = f_z->Eval(sciYup);
        islz = f_z->GetParameter(1);
        iinx = f_x->Eval(sciYup);
        islx = f_x->GetParameter(1);
        inScint = checkScintillator(2.5,iinx,islx,iinz,islz); // FIXME: need to tune
        fromSource = islx>-beamSlxMax&&islx<beamSlxMax&&islz>-beamSlzMax&&islz<beamSlzMax;
        nGood = getChi2XZ(nPairs,chi2x,chi2z);
        MyNamedVerbose("Tracking",Form("       2nd RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",islx,sciYup,iinx,chi2x,islz,sciYup,iinz,chi2z));
        if (!fromSource||!inScint) {f_x->SetParameters(0,0); f_z->SetParameters(0,0);}
        
        updateHitPositions(nPicks); // fix wy positions
        result = updatePairPositions(nPicks,nPairs);
        if (result) continue;
        setErrors(nPairs,false);
        fityz(nPairs);
        fityx(nPairs);
        iinz = f_z->Eval(sciYup);
        islz = f_z->GetParameter(1);
        iinx = f_x->Eval(sciYup);
        islx = f_x->GetParameter(1);
        inScint = checkScintillator(2.5,iinx,islx,iinz,islz); // FIXME: need to tune
        fromSource = islx>-beamSlxMax&&islx<beamSlxMax&&islz>-beamSlzMax&&islz<beamSlzMax;
        nGood = getChi2XZ(nPairs,chi2x,chi2z);
        MyNamedVerbose("Tracking",Form("       3rd RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",islx,sciYup,iinx,chi2x,islz,sciYup,iinz,chi2z));

        if (inScint&&fromSource&&nGood>=3){ // good candidate
            // update calD for all hits and driftD for no-pick hits
            // get hit list
            int nHitsSel = 0;
            t_sel->clear();
            t_sel->resize(i_nHits,0);
            for (int ihit = 0; ihit<i_nHits; ihit++){
                int lid = (*i_layerID)[ihit];
                int wid = (*i_wireID)[ihit];
                double calD = get_dist(lid,wid,islx,iinx,islz,iinz);
                (*t_calD)[ihit] = calD;
                if (!(*t_lr)[ihit]){ // not picked
                    if (calD>0) (*t_driftD)[ihit] = (*o_dxr)[ihit];
                    else (*t_driftD)[ihit] = (*o_dxl)[ihit];
                }
                double dd = (*t_driftD)[ihit];
                int selected = 0;
                int type = getHitType((*i_type)[ihit],calD>=0);
                if (fabs(calD-dd)<2&&m_testlayer!=(*i_layerID)[ihit]&&type<=3){ // FIXME: should tune the error limit
                	int jhit = getHitIndex(lid,ihit);
                    if (jhit==-1){ // no previous chosen hit yet
                        selected = 1;
                        nHitsSel++;
                    }
                    else{
                    	if (fabs(calD-dd)<fabs((*t_calD)[jhit]-(*t_driftD)[jhit])){ // better than the previous hit, then replace it
                    		selected = 1;
                    		(*t_sel)[jhit] = 0;
                    	}
                    }
                }
                (*t_sel)[ihit] = selected;
            }
            MyNamedVerbose("Tracking",Form("       good! nHitsSel = %d\n",nHitsSel));
            if (nHitsSel>=5){ // at least 5 hits to fit: NDF of 3-D track without field is 4
            	o_nFit++;
                // fitting with TMinuit
                do_fit(islx,iinx,islz,iinz);
                double temp;
                gMinuit->GetParameter(0, slx, temp);
                gMinuit->GetParameter(1, inx, temp);
                gMinuit->GetParameter(2, slz, temp);
                gMinuit->GetParameter(3, inz, temp);
                // update fitD
                // reselect
                nHitsSel = 0;
                t_sel->clear();
                t_sel->resize(i_nHits,0);
                for (int ihit = 0; ihit<i_nHits; ihit++){
                    int lid = (*i_layerID)[ihit];
                    int wid = (*i_wireID)[ihit];
                    double fitD = get_dist(lid,wid,slx,inx,slz,inz);
                    (*t_fitD)[ihit] = fitD;
                    if (!(*t_lr)[ihit]){ // not picked
                        if (fitD>0) (*t_driftD)[ihit] = (*o_dxr)[ihit];
                        else (*t_driftD)[ihit] = (*o_dxl)[ihit];
                    }
                    double dd = (*t_driftD)[ihit];
                    int selected = 0;
					int type = getHitType((*i_type)[ihit],fitD>=0);
                    if (fabs(fitD-dd)<1&&m_testlayer!=(*i_layerID)[ihit]&&type<=3){ // FIXME: should tune the error limit
						int jhit = getHitIndex(lid,ihit);
						if (jhit==-1){ // no previous chosen hit yet
                            selected = 1;
                            nHitsSel++;
                        }
						else{
							if (fabs(fitD-dd)<fabs((*t_fitD)[jhit]-(*t_driftD)[jhit])){ // better than the previous hit, then replace it
								selected = 1;
								(*t_sel)[jhit] = 0;
							}
						}
                    }
                    (*t_sel)[ihit] = selected;
                }
                if (nHitsSel>=5){ // at least 5 hits to fit: NDF of 3-D track without field is 4
                    // fitting with TMinuit
                    do_fit(islx,iinx,islz,iinz);
                    double temp;
                    gMinuit->GetParameter(0, slx, temp);
                    gMinuit->GetParameter(1, inx, temp);
                    gMinuit->GetParameter(2, slz, temp);
                    gMinuit->GetParameter(3, inz, temp);
                    inScint = checkScintillator(1.5,inx,slx,inz,slz); // FIXME: error limit should be tuned
                    fromSource = slx>-beamSlxMax&&slx<beamSlxMax&&slz>-beamSlzMax&&slz<beamSlzMax;
                    if (inScint&&fromSource){
                        // update chi2
                        if (m_inputType)
							getchi2(chi2mc,chi2pmc,chi2amc,i_slxmc,i_inxmc,i_slzmc,i_inzmc,true);
                        getchi2(chi2i,chi2pi,chi2ai,islx,iinx,islz,iinz,true);
                        getchi2(chi2,chi2p,chi2a,slx,inx,slz,inz,true);
                        // check chi2 and see where the result fits
                        MyNamedVerbose("Tracking",Form("         final RESULT: x=%.3e*(y-%.3e)+%.3e, z=%.3e*(y-%.3e)+%.3e, chi2i = %.3e chi2 = %.3e\n",slx,sciYup,inx,slz,sciYup,inz,chi2i,chi2));
                        // at last, update driftD for unpick and unselected hits
                        for (int ihit = 0; ihit<i_nHits; ihit++){
                            int lid = (*i_layerID)[ihit];
                            int wid = (*i_wireID)[ihit];
                            double fitD = get_dist(lid,wid,slx,inx,slz,inz);
                            (*t_fitD)[ihit] = fitD;
                            if (!(*t_lr)[ihit]&&!(*t_sel)[ihit]){ // not picked
                                if (fitD>0) (*t_driftD)[ihit] = (*o_dxr)[ihit];
                                else (*t_driftD)[ihit] = (*o_dxl)[ihit];
                            }
                            int type = getHitType((*i_type)[ihit],fitD>=0);
                            if (type<=3)
                                MyNamedVerbose("Tracking",Form("        %d (%d,%d) dd %.3e fd %.3e res %.3e\n",ihit,lid,wid,(*t_driftD)[ihit],fitD,fitD-(*t_driftD)[ihit]));
                            else 
                                MyNamedVerbose("Tracking",Form("              # %d (%d,%d) dd %.3e fd %.3e res %.3e\n",ihit,lid,wid,(*t_driftD)[ihit],fitD,fitD-(*t_driftD)[ihit]));
                        }
                        checkChi2(nHitsSel,nGood,icombi,iselection);
                    }
                }
            }
        }
    }
}

void setLRdriftD(int nPicks,int icombi){
    for (int ipick = 0; ipick<nPicks; ipick++){
        int ilr = (icombi&(1<<ipick))>>ipick;
        if (ilr==0) ilr = -1;
        int ihit = pick_ihit[ipick];
        (*t_lr)[ihit] = ilr;
        if (ilr>0) (*t_driftD)[ihit] = (*o_dxr)[ihit];
        else       (*t_driftD)[ihit] = (*o_dxl)[ihit];
    }
}

int getHitIndex(int lid, int nHits){
    int theHit = -1;
    for (int ihit = 0; ihit<nHits; ihit++){
        int tlid = (*i_layerID)[ihit];
        int sel = (*t_sel)[ihit];
        if (tlid==lid&&sel==1){
            theHit = ihit;
            break;
        }
    }
    return theHit;
}

int getHitType(int type,bool isRight){ // see if the driftT is really out of range
	int ttype = (type/10)%10;
	if (isRight){
		if (ttype==1||ttype==4) type-=ttype*10; // l- or l+
	}
	else{
		if (ttype==2||ttype==5) type-=ttype*10; // r- or r+
	}
    if (m_peakType>1) type=type%10000; // ignoring npoc cut and mpi cut, leaving all the peaks over threshold competing
    else if (m_peakType) type=type%100000; // ignoring npoc cut, leaving all the peaks with mpi==0 over threshold competing
    if (m_inputType==2||m_inputType==3) type==0;
	return type;
}

bool isSame(int iCand){
	bool allTheSame = true;
	for (int ihit = 0; ihit<i_nHits; ihit++){
		if ((*o_sel[iCand])[ihit] != (*t_sel)[ihit] || (*o_fitD[iCand])[ihit]*(*t_fitD)[ihit] <0 ){
			allTheSame = false;
			break;
		}
	}
	return allTheSame;
}

int updateHitPositions(int nPicks){
    // calculate pick_wy
    for (int ipick = 0; ipick<nPicks; ipick++){
        // Get hit information
        int ihit = pick_ihit[ipick];
        int lid = (*i_layerID)[ihit];
        int wid = (*i_wireID)[ihit];
        double wyro = map_y[lid][wid][1];
        double wzro = chamberHL;
        double wyhv = map_y[lid][wid][0];
        double wzhv = -chamberHL;
        double wy = (wyro+wyhv)/2.;// assume wy
        double wz = f_z->Eval(wy);// get wz by extrapolating the track to wy
        // correct wy according to wz
        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
        wz = f_z->Eval(wy);
        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
        pick_wy[ipick] = wy;
    }
}

int updatePairPositions(int nPicks,int & nPairs){
    // calculate pair_wxyz
    nPairs = 0;
    int ipick = 0;
    for (; ipick<nPicks-1; ipick++){
        int ihit = pick_ihit[ipick];
        int jhit = pick_ihit[ipick+1];
        int lid = (*i_layerID)[ihit];
        int ljd = (*i_layerID)[jhit];
        if (lid+1!=ljd) continue; // not adjacent
        double deltaY = pick_wy[ipick+1]-pick_wy[ipick];
        double dd1 = (*t_driftD)[ihit];
        double dd2 = (*t_driftD)[jhit];
        int wid = (*i_wireID)[ihit];
        int wjd = (*i_wireID)[jhit];
        double theta1 = map_theta[lid][wid];
        double theta2 = map_theta[ljd][wjd];
        double sintheta12 = sin(theta1-theta2);
        double zc_fix_slx = deltaY*f_x->GetParameter(1)/(tan(theta2)-tan(theta1));
        double xc = mcp_xc[lid][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
        double zc = mcp_zc[lid][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
        pair_wx[nPairs] = xc;
        pair_wy[nPairs] = (pick_wy[ipick+1]+pick_wy[ipick])/2.;
        pair_wz[nPairs] = zc;
//      MyNamedVerbose("Tracking",Form("                  xc = %.3e+%.3e*sin(%.3e)/(-sin(%.3e-%.3e))+%.3e*sin(%.3e)/sin(%.3e-%.3e)\n",mcp_xc[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2));
 //     MyNamedVerbose("Tracking",Form("                  zc = %.3e+%.3e*cos(%.3e)/(-sin(%.3e-%.3e))+%.3e*cos(%.3e)/sin(%.3e-%.3e)+%.3e\n",mcp_zc[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2,zc_fix_slx,zc_fix_slx));
        MyNamedVerbose("Tracking",Form("       cp[%d,%d]: w(%d,%d) i(%d,%d) dd(%f,%f) xyz(%f,%f,%f)\n",lid,ljd,wid,wjd,ihit,jhit,dd1,dd2,xc,(pick_wy[ipick+1]+pick_wy[ipick])/2.,zc));
        if (zc<-chamberHL||zc>chamberHL){
            MyNamedVerbose("Tracking",Form("       bad combination!\n"));
            break;
        }
        nPairs++;
    }
    if (ipick==nPicks-1){
        MyNamedVerbose("Tracking",Form("       GOOD!\n"));
        return 0;
    }
    else{
        MyNamedVerbose("Tracking",Form("       BAD @ %d!\n",ipick));
        return 1;
    }
}

int setErrors(int nPairs, bool noError){
    // Should not delete the graph!
    if (g_z) delete g_z;
    if (g_x) delete g_x;
    g_z = new TGraphErrors(nPairs,&(pair_wy[0]),&(pair_wz[0]),0,0);
    g_x = new TGraphErrors(nPairs,&(pair_wy[0]),&(pair_wx[0]),0,0);
//    g_x->Set(nPairs);
//    g_z->Set(nPairs);
//    for (int i = 0; i<nPairs; i++){
//        g_x->SetPoint(i,pair_wy[i],pair_wx[i]);
//        g_z->SetPoint(i,pair_wy[i],pair_wz[i]);
//    }
    double errorzMax0 = 0;
    double errorzMax1 = 0;
    int errorzMax0_i = -1;
    int errorzMax1_i = -1;
    for (int ipair = 0; ipair<nPairs; ipair++){
        double errorz = 0;
        double errorx = 0;
        if (!noError){
            errorz = fabs(f_z->Eval(pair_wy[ipair])-pair_wz[ipair]);
            errorx = fabs(f_x->Eval(pair_wy[ipair])-pair_wx[ipair]);
        }
        if (errorzMax0<errorz){
            errorzMax0 = errorz;
            errorzMax0_i = ipair;
        }
        else if (errorzMax1<errorz){
            errorzMax1 = errorz;
            errorzMax1_i = ipair;
        }
	    g_z->SetPointError(ipair,0,0.1);
	    g_x->SetPointError(ipair,0,0.1);
    }
    if (errorzMax0>4) {
        g_z->SetPointError(errorzMax0_i,0,10);
        g_x->SetPointError(errorzMax0_i,0,10);
        MyNamedVerbose("Tracking",Form("           setErrors pair[%d]: error x = ->10, error z = %.3e->10\n",errorzMax0_i,errorzMax0));
    }
    if (errorzMax1>4) {
        g_z->SetPointError(errorzMax1_i,0,10);
        g_x->SetPointError(errorzMax1_i,0,10);
        MyNamedVerbose("Tracking",Form("           setErrors pair[%d]: error x = ->10, error z = %.3e->10\n",errorzMax1_i,errorzMax1));
    }
    return 0;
}

int getChi2XZ(int nPairs, double & chi2x, double & chi2z){
    // calculate pair_wxyz
    chi2x = 0;
    chi2z = 0;
    int nCount = 0;
    for (int ipair = 0; ipair<nPairs; ipair++){
        double tchi2z = pair_wz[ipair]-f_z->Eval(pair_wy[ipair]);
        double tchi2x = pair_wx[ipair]-f_x->Eval(pair_wy[ipair]);
        MyNamedVerbose("Tracking",Form("           getChi2XZ pair[%d]: error x = %.3e, error z = %.3e\n",ipair,tchi2x,tchi2z));
	    if (fabs(tchi2z)<4&&fabs(tchi2x)<1){ // FIXME: error limit should be tuned
            chi2z += pow(tchi2z,2);
            chi2x += pow(tchi2x,2);
            nCount++;
        }
    }
    if (nCount){
        chi2z/=nCount;
        chi2x/=nCount;
    }
    return nCount;
}

int fityz(int nPairs){
	g_z->Fit("f_z","QN0G","");
	return 0;
}

int fityx(int nPairs){
	g_x->Fit("f_x","QN0G","");
	return 0;
}

bool checkScintillator(double saftyFactor,double inx, double slx, double inz, double slz){
    double xtop = 1/saftyFactor*inx;
    double xbot = 1/saftyFactor*(inx+slx*(sciYdown-sciYup));
    double ztop = 1/saftyFactor*inz;
    double zbot = 1/saftyFactor*(inz+slz*(sciYdown-sciYup));
    MyNamedVerbose("Tracking",Form("              top xz(%.3e, %.3e) bottom xz(%.3e, %.3e)\n",xtop,ztop,xbot,zbot));
    if (xtop>sciHW||xtop<-sciHW||xbot>sciHW||xbot<-sciHW||ztop>sciHL||ztop<-sciHL||zbot>sciHL||zbot<-sciHL) return false;
    else return true;
}

bool checkChi2(int nHitsSel, int nPairs, int icombi, int iselection){
	bool issame = false;
	bool covered = false;
    for (int i = 0; i<NCAND; i++){
    	issame = isSame(i);
    	if (issame){ // yes, there is a candidate with the same hits
            MyNamedVerbose("Tracking",Form(" same with Cand#%d where chi2=%.3e\n",i,o_chi2[i]));
//    		if (chi2<o_chi2[i]){// better? then remove the old one
			// FIXME: WARNING, now we rely on total chi2 including test layer hit, a slight bias
    		if (chi2a<o_chi2a[i]){// better? then remove the old one 
                MyNamedVerbose("Tracking","   better than Cand#"<<i);
    		    covered = true;
				for (int j = i; j<NCAND-1; j++){
					o_iselec[j] = o_iselec[j+1];
					o_icombi[j] = o_icombi[j+1];
					o_npairs[j] = o_npairs[j+1];
					o_islx[j] = o_islx[j+1];
					o_iinx[j] = o_iinx[j+1];
					o_islz[j] = o_islz[j+1];
					o_iinz[j] = o_iinz[j+1];
					o_chi2x[j] = o_chi2x[j+1];
					o_chi2z[j] = o_chi2z[j+1];
					o_chi2i[j] = o_chi2i[j+1];
					o_chi2pi[j] = o_chi2pi[j+1];
					o_chi2ai[j] = o_chi2ai[j+1];
					o_nHitsS[j] = o_nHitsS[j+1];
					o_slx[j] = o_slx[j+1];
					o_inx[j] = o_inx[j+1];
					o_slz[j] = o_slz[j+1];
					o_inz[j] = o_inz[j+1];
					o_chi2[j] = o_chi2[j+1];
					o_chi2p[j] = o_chi2p[j+1];
					o_chi2a[j] = o_chi2a[j+1];
					o_chi2mc[j]  = o_chi2mc[j+1];
					o_chi2amc[j] = o_chi2amc[j+1];
					o_chi2pmc[j] = o_chi2pmc[j+1];
					for (int ihit = 0; ihit<i_nHits; ihit++){
						(*o_sel[j])[ihit] = (*o_sel[j+1])[ihit];
						(*o_calD[j])[ihit] = (*o_calD[j+1])[ihit];
						(*o_fitD[j])[ihit] = (*o_fitD[j+1])[ihit];
						(*o_driftD[j])[ihit] = (*o_driftD[j+1])[ihit];
					}
				}
				o_chi2[NCAND-1] = 1e9;
				o_nHitsS[NCAND-1] = 0;
			}
			break;
    	}
	}
	if (!issame||covered){ // didn't find a candidate with the same hits
		for (int i = 0; i<NCAND; i++){
//			if ((chi2<o_chi2[i]&&nHitsSel==o_nHitsS[i])||nHitsSel>o_nHitsS[i]){ // now we only pick up one hit per layer since the XT shape in the corener is very sensitive to position/angle thus less reliable
			// FIXME: WARNING, now we rely on total chi2 including test layer hit, a slight bias
			if ((chi2a<o_chi2a[i]&&nHitsSel==o_nHitsS[i])||nHitsSel>o_nHitsS[i]){ // now we only pick up one hit per layer since the XT shape in the corener is very sensitive to position/angle thus less reliable
                MyNamedVerbose("Tracking",Form("better than Cand#%d where chi2=%.3e, nHitsS=%d\n",i,o_chi2[i],o_nHitsS[i]));
				for (int j = NCAND-1; j>i; j--){
					o_iselec[j] = o_iselec[j-1];
					o_icombi[j] = o_icombi[j-1];
					o_npairs[j] = o_npairs[j-1];
					o_islx[j] = o_islx[j-1];
					o_iinx[j] = o_iinx[j-1];
					o_islz[j] = o_islz[j-1];
					o_iinz[j] = o_iinz[j-1];
					o_chi2x[j] = o_chi2x[j-1];
					o_chi2z[j] = o_chi2z[j-1];
					o_chi2i[j] = o_chi2i[j-1];
					o_chi2pi[j] = o_chi2pi[j-1];
					o_chi2ai[j] = o_chi2ai[j-1];
					o_nHitsS[j] = o_nHitsS[j-1];
					o_slx[j] = o_slx[j-1];
					o_inx[j] = o_inx[j-1];
					o_slz[j] = o_slz[j-1];
					o_inz[j] = o_inz[j-1];
					o_chi2[j] = o_chi2[j-1];
					o_chi2p[j] = o_chi2p[j-1];
					o_chi2a[j] = o_chi2a[j-1];
					for (int ihit = 0; ihit<i_nHits; ihit++){
						(*o_sel[j])[ihit] = (*o_sel[j-1])[ihit];
						(*o_calD[j])[ihit] = (*o_calD[j-1])[ihit];
						(*o_fitD[j])[ihit] = (*o_fitD[j-1])[ihit];
						(*o_driftD[j])[ihit] = (*o_driftD[j-1])[ihit];
					}
				}
				o_iselec[i] = iselection;
				o_icombi[i] = icombi;
				o_npairs[i] = nPairs;
				o_islx[i] = islx;
				o_iinx[i] = iinx;
				o_islz[i] = islz;
				o_iinz[i] = iinz;
				o_chi2z[i] = chi2z;
				o_chi2x[i] = chi2x;
				o_chi2i[i] = chi2i;
				o_chi2ai[i] = chi2ai;
				o_chi2pi[i] = chi2pi;
				o_nHitsS[i] = nHitsSel;
				o_slx[i] = slx;
				o_inx[i] = inx;
				o_slz[i] = slz;
				o_inz[i] = inz;
				o_chi2[i] = chi2;
				o_chi2p[i] = chi2p;
				o_chi2a[i] = chi2a;
				o_chi2mc[i] = chi2mc;
				o_chi2amc[i] = chi2amc;
				o_chi2pmc[i] = chi2pmc;
				for (int ihit = 0; ihit<i_nHits; ihit++){
					(*o_sel[i])[ihit] = (*t_sel)[ihit];
					(*o_calD[i])[ihit] = (*t_calD)[ihit];
					(*o_fitD[i])[ihit] = (*t_fitD)[ihit];
					(*o_driftD[i])[ihit] = (*t_driftD)[ihit];
				}
				break;
			}
		}
	}
}

double t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
    TF1* f=0;
    int theLayer = lid;
    if (m_workType==0) theLayer = 0;
	else if (m_workType==1) theLayer = (lid%2==0?NLAY:NLAY+1); // even/odd
	else if (m_workType==-1) theLayer = (lid%2==0?NLAY+1:NLAY); // even/odd reversed
    if (lr>=0){
        f = f_right[theLayer];
    }
    else {
        f = f_left[theLayer];
    }
    if (!f){
        fprintf(stderr,"Cannot get f[%d]!\n",theLayer);
        status = -2;
        return 0;
    }
    double m_tmax = f->GetXmax();
    double m_tmin = f->GetXmin();
    status = 0;
    double dd = 0;
    if (time>m_tmax){
        status = 1;
        dd = f->Eval(m_tmax);
    }
    else if (time<m_tmin){
        status = -1;
        dd = 0;
    }
    else {
        status = 0;
        dd = f->Eval(time);
    }
    return dd;
}

//______________________________________________________________________________
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz)
{
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
	return value;
}

//______________________________________________________________________________
void do_fit(double sliX, double iniX,double sliZ, double iniZ){
	if(gMinuit) delete gMinuit;
	gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);
	arglist[0] = 0;
	gMinuit->SetPrintLevel(-1); // no print
	gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	//if (m_inputType){
	//	iniZ = i_inzmc;
	//	sliZ = i_slzmc;
	//	iniX = i_inxmc;
	//	sliX = i_slxmc;
	//}

	// Set starting values and step sizes for parameters
	gMinuit->mnparm(0, "slopeX", sliX, beamSlxMax/1.e4, -beamSlxMax,beamSlxMax,ierflg);
	gMinuit->mnparm(1, "interceptX", iniX, beamInxMax/1.e4, -beamInxMax,beamInxMax,ierflg);
	gMinuit->mnparm(2, "slopeZ", sliZ, beamSlzMax/1.e4, -beamSlzMax,beamSlzMax,ierflg);
	gMinuit->mnparm(3, "interceptZ", iniZ, beamInzMax/1.e4, -beamInzMax,beamInzMax,ierflg);

	// Now ready for minimization step
	arglist[0] = 5000.0;
	arglist[1] = 1.0;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Print results
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//printf("====Rrestul====\n");
	//gMinuit->mnprin(3,amin);
}

//______________________________________________________________________________
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
	double cp,ca;
	getchi2(f,cp,ca,*par,*(par+1),*(par+2),*(par+3),false);
}

//______________________________________________________________________________
void getchi2(double &f, double & cp, double & ca, double slx, double inx, double slz, double inz,bool all)
{
	//calculate chisquare
	double chisq = 0;
	double delta;
	double dfit;

    int N = -4;
	for (int ihit=0;ihit<i_nHits; ihit++) {
		if ((*t_sel)[ihit]==0) continue;
		dfit = get_dist((*i_layerID)[ihit],(*i_wireID)[ihit],slx,inx,slz,inz);
		double dd = (*t_driftD)[ihit];
		double error = getError((*i_layerID)[ihit],(*i_driftT)[ihit],dd>0);
        delta  = (dfit-dd)/error;
		chisq += delta*delta;
		N++;
	}
	if (N>0) chisq/=N;
	f = chisq;
	if (all){ // should calculate chi2_pValue (cp) and chi2_all (ca)
		cp = TMath::Prob(f*N,N);
		// now find the closest hit in the test layer and add it into ca
		double minres = 1e9;
		bool found = false;
		for (int ihit=0;ihit<i_nHits; ihit++) {
			if ((*i_layerID)[ihit]!=m_testlayer) continue;
			dfit = get_dist((*i_layerID)[ihit],(*i_wireID)[ihit],slx,inx,slz,inz);
			int type = getHitType((*i_type)[ihit],dfit>=0);
			if (type>3) continue;
			double dd = dfit>0?(*o_dxr)[ihit]:(*o_dxl)[ihit];
			if (fabs(minres)>fabs(dfit-dd)){
				minres = dfit-dd;
				found = true;
			}
		}
		if (found){
			double error = 0.2;
			ca = f*N+minres*minres/error/error;
			ca/=(N+1);
		}
		else{
			ca = f;
		}
	}
}

double getError(int lid,double dt, bool isR){
	double error = 0.2; // default value 200 um
	TGraph * gr = 0;
    int theLayer = lid;
    if (m_workType==0) theLayer = 0;
	else if (m_workType==1) theLayer = (lid%2==0?NLAY:NLAY+1); // even/odd
	else if (m_workType==-1) theLayer = (lid%2==0?NLAY+1:NLAY); // even/odd reversed
	if (isR) // right side
		gr = gr_error_right[theLayer];
	else
		gr = gr_error_left[theLayer];
	if (!gr){
	    return error;
	}
	int N = gr->GetN();
	for (int i = 0; i<N-1; i++){
		double t1,sig1;
		double t2,sig2;
		gr->GetPoint(i,t1,sig1);
		gr->GetPoint(i+1,t2,sig2);
		if (t1<dt&&t2>=dt){
			error = (sig1+sig2)/2.;
		}
	}
	return error;
}

//______________________________________________________________________________
void print_usage(char* prog_name)
{
	fprintf(stderr,"Usage %s [options] prerunname runname\n",prog_name);
	fprintf(stderr,"[options]\n");
	fprintf(stderr,"\t -D <name>=[error,severe,warn,debug,trace]\n");
	fprintf(stderr,"\t\t Change the named debug level\n");
	fprintf(stderr,"\t -V <name>=[quiet,log,info,verbose]\n");
	fprintf(stderr,"\t\t Change the named log level\n");
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
    fprintf(stderr,"\t\t Test layer set to l\n");
    fprintf(stderr,"\t -n <n>\n");
    fprintf(stderr,"\t\t Maximum number of hits cut set to n\n");
    fprintf(stderr,"\t -x <x>\n");
    fprintf(stderr,"\t\t T0 shift on board 0 set to x\n");
    fprintf(stderr,"\t -y <y>\n");
    fprintf(stderr,"\t\t T0 shift on board 1 set to y\n");
    fprintf(stderr,"\t -l <l>\n");
    fprintf(stderr,"\t\t Minimum time cut set to l\n");
    fprintf(stderr,"\t -u <u>\n");
    fprintf(stderr,"\t\t Maximum time cut set to u\n");
    fprintf(stderr,"\t -s <s>\n");
    fprintf(stderr,"\t\t ADC sum over peak cut set to s\n");
    fprintf(stderr,"\t -a <a>\n");
    fprintf(stderr,"\t\t ADC sum over all cut set to a\n");
    fprintf(stderr,"\t -g <g>\n");
    fprintf(stderr,"\t\t Geometry setup set to g\n");
    fprintf(stderr,"\t\t (0): normal; 1: finger\n");
    fprintf(stderr,"\t -i <i>\n");
    fprintf(stderr,"\t\t Input type set to i\n");
    fprintf(stderr,"\t\t (0) for data; 1 for MC\n");
    fprintf(stderr,"\t -p <p>\n");
    fprintf(stderr,"\t\t Peak type set to p\n");
    fprintf(stderr,"\t\t (0) only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks\n");
    fprintf(stderr,"\t -w <w>\n");
    fprintf(stderr,"\t\t Work type set to w\n");
    fprintf(stderr,"\t\t 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers\n");
}
