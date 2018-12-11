#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>

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

#include "MyProcessManager.hxx"
#include "MyRuntimeParameters.hxx"
#include "Log.hxx"
#include "header.hxx"

#define NBINS    20
#define MAXTRUNC 6

//===================Configure Parameter============================
int m_runNo = 0;
TString m_runname = "currun";
int m_iEntryStart = 0;
int m_iEntryStop = 0;
int m_modulo = 10000;
bool m_memdebug = false;
bool m_savehists = false;
bool m_outputEventTree = false;
int m_testlayer = 4;
TString m_configureFile = "";
int m_geoSetup = 0; // 0: normal; 1: finger
// for cutting
int m_nHitsMax = 0;
int m_nHitsSmin = 7;
int m_sumCut = 0;
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

//===================Chamber Parameter============================
// map for wire position
double  map_x[NLAY][NCEL][2];
double  map_y[NLAY][NCEL][2];
double  map_off[NLAY][NCEL];
bool    map_has[NLAY][NCEL];
double  npair_per_cm = 0;
TF1 * f_left0 = 0;
TF1 * f_right0 = 0;
TF1 * f_left = 0;
TF1 * f_right = 0;
TF1 * f_left_mid = 0;
TF1 * f_right_mid = 0;
TF1 * f_left_cent = 0;
TF1 * f_right_cent = 0;
TF1 * f_left_end = 0;
TF1 * f_right_end = 0;
int m_xtType = 0;
double t7l = 0;
double t8l = 0;
double t7r = 0;
double t8r = 0;

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
int t2d(double t, double & d, bool isRight);
double findFirstX(TF1 * f, double val, double xmin, double xmax, double delta);
void doFit(TH1D * h,double leftRatio = 1/3., double rightRatio = 1/3., double leftEnd = -m_maxRes, double rightEnd = m_maxRes);
void print_usage(char * prog_name);

MyProcessManager * pMyProcessManager;

int main(int argc, char** argv){
    int temp_geoSetup = 0; bool set_geoSetup = false;
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

    std::map<std::string, Log::ErrorPriority> namedDebugLevel;
    std::map<std::string, Log::LogPriority> namedLogLevel;
    int    opt_result;
	while((opt_result=getopt(argc,argv,"M:R:B:E:L:C:n:f:c:p:r:z:d:o:s:a:l:u:t:x:y:g:D:V:"))!=-1){
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
			case 'f':
			    temp_nHitsSmin = atoi(optarg);set_nHitsSmin = true;
                printf("Minimum number of selected hits cut set to %d\n",temp_nHitsSmin);
			case 'c':
			    temp_maxchi2 = atof(optarg);set_maxchi2 = true;
                printf("Maximum chi2 cut set to %d\n",temp_maxchi2);
			case 'p':
			    temp_minchi2p = atof(optarg);set_minchi2p = true;
                printf("Minimum p-value cut set to %d\n",temp_minchi2p);
			case 'r':
			    temp_maxRes = atof(optarg);set_maxRes = true;
                printf("Maximum resolution cut set to %d\n",temp_maxRes);
			case 'z':
			    temp_maxslz = atof(optarg);set_maxslz = true;
                printf("Maximum y-z slope cut set to %d\n",temp_maxslz);
			case 'd':
			    temp_maxFD = atof(optarg);set_maxFD = true;
                printf("Maximum fitD cut set to %d\n",temp_maxFD);
			case 'o':
			    temp_tmaxSet = atoi(optarg);set_tmaxSet = true;
                printf("Maximum time range set to %d\n",temp_tmaxSet);
			case 's':
			    temp_sumCut = atoi(optarg);set_sumCut = true;
                printf("ADC sum over peak cut set to %d\n",temp_sumCut);
			case 'a':
			    temp_aaCut = atoi(optarg);set_aaCut = true;
                printf("ADC sum over all cut set to %d\n",temp_aaCut);
			case 'l':
			    temp_tmin = atof(optarg);set_tmin = true;
                printf("Minimum time on axis set to %d\n",temp_tmin);
			case 'u':
			    temp_tmax = atof(optarg);set_tmax = true;
                printf("Maximum time on axis set to %d\n",temp_tmax);
			case 't':
			    temp_NbinT = atoi(optarg);set_NbinT = true;
                printf("Number of bins on time axis set to %d\n",temp_NbinT);
			case 'x':
			    temp_NbinX = atoi(optarg);set_NbinX = true;
                printf("Number of bins on space axis set to %d\n",temp_NbinX);
			case 'y':
			    temp_NbinRes = atoi(optarg);set_NbinRes = true;
                printf("Number of bins on resolution axis set to %d\n",temp_NbinRes);
			case 'g':
			    temp_geoSetup = atoi(optarg);set_geoSetup = true;
                printf("Geometry setup set to %d\n",temp_geoSetup);
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
        m_geoSetup = MyRuntimeParameters::Get().GetParameterI("geoSetup");
        //for cutting
        m_nHitsMax = MyRuntimeParameters::Get().GetParameterI("ana.nHitsMax");
        m_nHitsSmin = MyRuntimeParameters::Get().GetParameterI("ana.nHitsSmin");
        m_sumCut = MyRuntimeParameters::Get().GetParameterI("ana.sumCut");
        m_aaCut = MyRuntimeParameters::Get().GetParameterI("ana.aaCut");
        m_maxchi2 = MyRuntimeParameters::Get().GetParameterD("ana.maxchi2");
        m_maxslz = MyRuntimeParameters::Get().GetParameterD("ana.maxslz");
        m_maxFD = MyRuntimeParameters::Get().GetParameterD("ana.maxFD");
        m_tmaxSet = MyRuntimeParameters::Get().GetParameterI("ana.tmaxSet");
        //for binning
        m_tmin = MyRuntimeParameters::Get().GetParameterD("ana.tmin");
        m_tmax = MyRuntimeParameters::Get().GetParameterD("ana.tmax");
        m_xmax = MyRuntimeParameters::Get().GetParameterD("ana.xmax");
        m_NbinT = MyRuntimeParameters::Get().GetParameterI("ana.NbinT");
        m_NbinX = MyRuntimeParameters::Get().GetParameterI("ana.NbinX");
        m_NbinRes = MyRuntimeParameters::Get().GetParameterI("ana.NbinRes");
        m_minchi2p = MyRuntimeParameters::Get().GetParameterD("ana.minchi2p");
        m_maxRes = MyRuntimeParameters::Get().GetParameterD("ana.maxRes");
    }
    if (set_geoSetup) m_geoSetup = temp_geoSetup;
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

	if (argc-optind<1){
	    print_usage(argv[0]);
		return -1;
    }
    m_runname= argv[optind++];

    printf("##############%s##################\n",argv[0]);
    printf("runNo       = %d\n",m_runNo);
    printf("runname     = \"%s\"\n",m_runname.Data());
    printf("test layer:   %d\n",m_testlayer);
    printf("xtType:       %d\n",m_xtType);
    printf("maxchi2     = %.3e\n",m_maxchi2);
    printf("maxslz      = %.3e\n",m_maxslz);
    printf("nHits max   = %d\n",m_nHitsMax);
    printf("nHitsSmin   = %d\n",m_nHitsSmin);
    printf("Q cut       = %d\n",m_aaCut);
    printf("geoSetup    = %d\n",m_geoSetup);
    printf("tmaxSet     = %d\n",m_tmaxSet);
    printf("savehists:    %s\n",m_savehists?"yes":"no");
    printf("output EventTree? %s\n",m_outputEventTree?"yes":"no");
    printf("Entries:     [%d~%d]\n",m_iEntryStart,m_iEntryStop);
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
	printf("m_runNo#%d: %s, %d, %s, %d V, %d mV, %.0f sec\n",m_runNo,gastype.Data(),runGr,duration.Data(),HV,THR,durationTime);

    //Prepare Maps
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
			map_has[lid][wid] = false;
			map_off[lid][wid] = 0;
        }
    }

	// get offset
	if (m_xtType!=-1){
		TChain * iChain_off = new TChain("t","t");
		iChain_off->Add(Form("%s/info/offset.%d.%s.root",HOME.Data(),m_runNo,m_runname.Data()));
		double i_off_delta;
		int i_off_lid;
		int i_off_wid;
		iChain_off->SetBranchAddress("d",&i_off_delta);
		iChain_off->SetBranchAddress("wid",&i_off_wid);
		iChain_off->SetBranchAddress("lid",&i_off_lid);
		int N = iChain_off->GetEntries();
		for (int i = 0; i<N; i++){
			iChain_off->GetEntry(i);
			if (i_off_lid>=0&&i_off_lid<NLAY&&i_off_wid>=0&&i_off_wid<NCEL)
				map_off[i_off_lid][i_off_wid] = i_off_delta;
		}
	}

    // Get Wire Position
    TChain * TChain_wirepos = new TChain("t","t");
    TChain_wirepos->Add(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),m_runNo,m_runname.Data()));
    if (!TChain_wirepos->GetEntries()){
    	fprintf(stderr,"Cannot find %s/info/wire-position.%d.%s.root, using default one\n",HOME.Data(),m_runNo,m_runname.Data());
		TChain_wirepos->Add(Form("%s/info/wire-position.root",HOME.Data()));
    }
    if (!TChain_wirepos->GetEntries()){
    	fprintf(stderr,"Cannot find default wire-position!\n");
    	return 1;
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
            map_has[wp_lid][wp_wid] = true;
        }
    }

    // get XT file
    TFile * XTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_runname.Data()));
	f_left0 = (TF1*) XTFile->Get("fl_0");
	f_right0 = (TF1*) XTFile->Get("fr_0");
	f_left = (TF1*) XTFile->Get(Form("flc_%d",m_testlayer));
	f_right = (TF1*) XTFile->Get(Form("frc_%d",m_testlayer));
	f_left_cent = (TF1*) XTFile->Get(Form("flce_%d",m_testlayer));
	f_right_cent = (TF1*) XTFile->Get(Form("frce_%d",m_testlayer));
	f_left_mid = (TF1*) XTFile->Get(Form("flm_%d",m_testlayer));
	f_right_mid = (TF1*) XTFile->Get(Form("frm_%d",m_testlayer));
	f_left_end = (TF1*) XTFile->Get(Form("fle_%d",m_testlayer));
	f_right_end = (TF1*) XTFile->Get(Form("fre_%d",m_testlayer));
    if (!f_left||!f_right|!f_left0||!f_right0||!f_left_mid||!f_right_mid||!f_left_cent||!f_right_cent||!f_left_end||!f_right_end){
    	fprintf(stderr,"Cannot find XT functions!\n");
    	return 0;
	}
	if (m_xtType%10==0||m_xtType%10==1){
		t7r = findFirstX(f_right,7,-10,800,10);
		t8r = findFirstX(f_right,8,-10,800,10);
		t7l = findFirstX(f_left,-7,-10,800,10);
		t8l = findFirstX(f_left,-8,-10,800,10);
	}
	else{
		t7r = findFirstX(f_right_mid,7,-10,800,10);
		t8r = findFirstX(f_right_mid,8,-10,800,10);
		t7l = findFirstX(f_left_mid,-7,-10,800,10);
		t8l = findFirstX(f_left_mid,-8,-10,800,10);
	}
    double minDT = 0;
    if (m_xtType%10==3)
		minDT = f_left_mid->GetXmin()>f_right_mid->GetXmin()?f_right_mid->GetXmin():f_left_mid->GetXmin();
    else
		minDT = f_left->GetXmin()>f_right->GetXmin()?f_right->GetXmin():f_left->GetXmin();
    double maxDT = 0;
	if (m_xtType/10==0){
		if (m_xtType%10==0)
			maxDT = f_left0->GetXmax()<f_right0->GetXmax()?f_right0->GetXmax():f_left0->GetXmax();
		else if (m_xtType%10==1)
			maxDT = f_left->GetXmax()<f_right->GetXmax()?f_right->GetXmax():f_left->GetXmax();
		else
			maxDT = f_left_mid->GetXmax()<f_right_mid->GetXmax()?f_right_mid->GetXmax():f_left_mid->GetXmax();
	}
	else{
	    if (m_tmaxSet)
	        maxDT = m_tmaxSet;
        else
            maxDT = t8l<t8r?t8r:t8l;
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
    printf("##############Geometry##################\n");
    printf("sciYup      = %.3e\n",sciYup);
    printf("sciYdown    = %.3e\n",sciYdown);
    printf("sciHL       = %.3e\n",sciHL);
    printf("sciHW       = %.3e\n",sciHW);

	// set RECBE ADC function
	//fADC2ChargeFunction = new TF1("a2c","5.98739+2.6652*x+0.000573394*x*x-5.21769e-05*x*x*x+3.05897e-07*x*x*x*x-7.54057e-10*x*x*x*x*x+8.60252e-13*    x*x*x*x*x*x-3.68603e-16*x*x*x*x*x*x*x",-10,800);
	//fADC2ChargeFunction = new TF1("a2c","5.98739+2.6652*x",-10,800);
	fADC2ChargeFunction = new TF1("a2c","2.0158464*x",-10,800);

	//==============================================Prepare input file & output variables=================================================
    // input file
    int triggerNumber;
    int nHits;
    int nHitsG;
    int npairs;
    int isel;
    int icom;
    double iinx;
    double islx;
    double iinz;
    double islz;
    double chi2x;
    double chi2z;
    double chi2i;
    int nHitsS;
    double inx;
    double slx;
    double inz;
    double slz;
    double chi2;
    double chi2p;
    double chi2a;
	int nHitsSmallAll = 0;
	int nHitsSmallSASD = 0;
	int nSmallSumHits = 0;
	int nShadowedHits = 0;
	int nLateHits = 0;
	int nBoundaryHits = 0;
	int nSmallBoundaryHits = 0;
	// the closest peak to the track in the test layer
	double res = 1e9;
	double theDD = 1e9;
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
	// the highest hit in this event
	int highBid = 0;
	int highCh = 0;
	int highLid = 0;
	int highWid = 0;
	int highIp = 0;
	double highSum = 0;
	double highAA = 0;
	double highDT = 0;
	// hit list
	std::vector<int> * i_layerID = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_ip = 0;
	std::vector<int> * i_sel = 0;
	std::vector<int> * i_peak = 0;
	std::vector<double> * i_ped = 0;
	std::vector<double> * i_aa = 0;
	std::vector<double> * i_driftT = 0;
	std::vector<double> * i_driftD = 0;
	std::vector<double> * i_fitD = 0;
	//----------------------------------Set input file--------------------------------------------
	TChain * ichain = new TChain("t","t");
	ichain->Add(Form("%s/root/ana_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer));
	ichain->SetBranchAddress("triggerNumber",&triggerNumber);
	ichain->SetBranchAddress("res",&res);
	ichain->SetBranchAddress("theDD",&theDD);
	ichain->SetBranchAddress("theDT",&theDT);
	ichain->SetBranchAddress("theWid",&theWid);
	ichain->SetBranchAddress("theSum",&theSum);
	ichain->SetBranchAddress("sum1st",&sum1st);
	ichain->SetBranchAddress("dt1st",&dt1st);
	ichain->SetBranchAddress("has",&has);
	ichain->SetBranchAddress("thePeak",&thePeak);
	ichain->SetBranchAddress("theHeight",&theHeight);
	ichain->SetBranchAddress("theIp",&theIp);
	ichain->SetBranchAddress("theMpi",&theMpi);
	ichain->SetBranchAddress("highBid",&highBid);
	ichain->SetBranchAddress("highCh",&highCh);
	ichain->SetBranchAddress("highLid",&highLid);
	ichain->SetBranchAddress("highWid",&highWid);
	ichain->SetBranchAddress("highIp",&highIp);
	ichain->SetBranchAddress("highSum",&highSum);
	ichain->SetBranchAddress("highAA",&highAA);
	ichain->SetBranchAddress("highDT",&highDT);
	ichain->SetBranchAddress("nHitsSmallSASD",&nHitsSmallSASD);
	ichain->SetBranchAddress("nHitsSmallAll",&nHitsSmallAll);
	ichain->SetBranchAddress("nSHits",&nShadowedHits);
	ichain->SetBranchAddress("nLHits",&nLateHits);
	ichain->SetBranchAddress("nSSHits",&nSmallSumHits);
	ichain->SetBranchAddress("nBHits",&nBoundaryHits);
	ichain->SetBranchAddress("nSBHits",&nSmallBoundaryHits);
	ichain->SetBranchAddress("nHits",&nHits);
	ichain->SetBranchAddress("nHitsG",&nHitsG);
	ichain->SetBranchAddress("npairs",&npairs);
	ichain->SetBranchAddress("isel",&isel);
	ichain->SetBranchAddress("icom",&icom);
	ichain->SetBranchAddress("islx",&islx);
	ichain->SetBranchAddress("islz",&islz);
	ichain->SetBranchAddress("iinx",&iinx);
	ichain->SetBranchAddress("iinz",&iinz);
	ichain->SetBranchAddress("chi2x",&chi2x);
	ichain->SetBranchAddress("chi2z",&chi2z);
	ichain->SetBranchAddress("chi2i",&chi2i);
	ichain->SetBranchAddress("nHitsS",&nHitsS);
	ichain->SetBranchAddress("slx",&slx);
	ichain->SetBranchAddress("slz",&slz);
	ichain->SetBranchAddress("inx",&inx);
	ichain->SetBranchAddress("inz",&inz);
	ichain->SetBranchAddress("chi2",&chi2);
	ichain->SetBranchAddress("chi2p",&chi2p);
	ichain->SetBranchAddress("chi2a",&chi2a);
	ichain->SetBranchAddress("layerID",&i_layerID);
	ichain->SetBranchAddress("wireID",&i_wireID);
	ichain->SetBranchAddress("ip",&i_ip);
	ichain->SetBranchAddress("aa",&i_aa);
	ichain->SetBranchAddress("peak",&i_peak);
	ichain->SetBranchAddress("ped",&i_ped);
	ichain->SetBranchAddress("sel",&i_sel);
	ichain->SetBranchAddress("driftT",&i_driftT);
	ichain->SetBranchAddress("driftD",&i_driftD);
	ichain->SetBranchAddress("fitD",&i_fitD);

    TFile * ofile_events = 0;
    TTree * otree_events = 0;
    bool isGood = false;
    double theCharge = 0;
	double chargeOnTrack[NLAY]; // charge along the track
	double adcsumOnTrack[NLAY]; // ADC sum along the track
	int    chargeOnTrackIndex[NLAY]; // index of the corresponding hit along the track
    int    nHitsT = 0;
    bool   o_hashit[2];
    double o_theFD[2];
    double o_theDD[2];
    double o_theDT[2];
    int    o_theST[2];
    double o_theAA[2];
    double theGG = 0;
    double trackGG = 0;
    int nLayersOnTrack = 0;
    ofile_events = new TFile(Form("%s/root/eres_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer),"RECREATE");
    otree_events = new TTree("tree","tree");
    otree_events->Branch("triggerNumber",&triggerNumber);
    otree_events->Branch("isGood",&isGood);
    otree_events->Branch("res",&res);
    otree_events->Branch("nHitsT",&nHitsT);
    otree_events->Branch("hashit",o_hashit,"hashit[nHitsT]/O");
    otree_events->Branch("theFD",o_theFD,"theFD[nHitsT]/D");
    otree_events->Branch("theDD",o_theDD,"theDD[nHitsT]/D");
    otree_events->Branch("theDT",o_theDT,"theDT[nHitsT]/D");
    otree_events->Branch("theST",o_theST,"theST[nHitsT]/I");
    otree_events->Branch("theAA",o_theAA,"theAA[nHitsT]/D");
    otree_events->Branch("theCharge",&theCharge);
    for (int i = 0; i<NLAY; i++){
        otree_events->Branch(Form("chargeOnTrack%d",i),&(chargeOnTrack[i]));
        otree_events->Branch(Form("adcsumOnTrack%d",i),&(adcsumOnTrack[i]));
        otree_events->Branch(Form("chargeOnTrackIndex%d",i),&(chargeOnTrackIndex[i]));
    }
    otree_events->Branch("theGG",&theGG);
    otree_events->Branch("trackGG",&trackGG);
    otree_events->Branch("nLayers",&nLayersOnTrack);
    otree_events->Branch("nBHits",&nBoundaryHits);
    otree_events->Branch("nSBHits",&nSmallBoundaryHits);
    otree_events->Branch("nHits",&nHits);
    otree_events->Branch("nHitsG",&nHitsG);
    otree_events->Branch("npairs",&npairs);
    otree_events->Branch("slx",&slx);
    otree_events->Branch("slz",&slz);
    otree_events->Branch("inx",&inx);
    otree_events->Branch("inz",&inz);
    otree_events->Branch("chi2",&chi2);
    otree_events->Branch("chi2p",&chi2p);
    otree_events->Branch("chi2a",&chi2a);

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

	//=================================================Loop in events====================================================
	double closeFD;
	double closeFD2;
	int    closeWid;
	int    closeWid2;
	double averageGG = 0;
	double averageGGErr = 0;
    double trackCharge[MAXTRUNC] = {0};
	Long64_t N = ichain->GetEntries();
	if (!m_iEntryStart&&!m_iEntryStop){
		m_iEntryStart = 0;
		m_iEntryStop = N-1;
	}
    double closestchi2 = 1e9;
    MyNamedLog("Ana","Processing "<<N<<" events");
	for ( int iEntry = m_iEntryStart ; iEntry<=m_iEntryStop; iEntry++){
		if (iEntry%10000==0) printf("%d\n",iEntry);
        MyNamedVerbose("Ana","Entry "<<iEntry);
		ichain->GetEntry(iEntry);

		// update m_minchi2p
		if (fabs(chi2-m_maxchi2)<fabs(closestchi2-m_maxchi2)){
			closestchi2 = chi2;
			m_minchi2p = chi2p;
		}

		// ignore events with bad fitting
		isGood = true;
		N_ALL++;
		h_nHits->Fill(nHits);
		if (m_nHitsMax&&nHits>m_nHitsMax) isGood = false;
		if (isGood) N_CUT1++;
		if (isGood) h_DOF->Fill(nHitsS-4);
		if (nHitsS<m_nHitsSmin) isGood = false;
		if (isGood) N_CUT2++;
		if (isGood) h_chi2->Fill(chi2);
		if (chi2>m_maxchi2) isGood = false;
		if (isGood) N_CUT3++;
		if (isGood) h_slz->Fill(slz);
		if (fabs(slz)>m_maxslz) isGood = false;
		if (isGood) N_CUT4++;
        MyNamedVerbose("Ana","  Good Event");

		// get closest wire
		closeFD = 1e3;
		closeWid = 0;
		closeFD2 = 1e3;
		closeWid2 = 0;
		for (int wid = 0; wid<NCEL; wid++){
			double fitD = get_dist(m_testlayer,wid,slx,inx,slz,inz);
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
		o_theFD[0] = closeFD;
		int ibinX = fabs(closeFD)/m_xmax*NBINS; // bin index for DOCA
		if (ibinX>=NBINS) isGood = false; // too far away from any cell
		if (isGood) N_CUT5++;
        nHitsT = 1; // number of hits in the test layer is by default 1
        int ibinX2 = fabs(closeFD2)/m_xmax*NBINS; // bin index for the next closest DOCA
        if (fabs(closeFD2+closeFD)<4&&ibinX2<NBINS){
            nHitsT = 2; // there is one more hit in the test layer
            o_theFD[1] = closeFD2;
        }
        // DOCA
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
		double minRes[2] = {1e3}; // minimal residual in the test layer
		o_hashit[0] = false;
		o_hashit[1] = false;
		for (int lid = 0; lid<NLAY; lid++){
            chargeOnTrack[lid] = 0;
            adcsumOnTrack[lid] = 0;
            chargeOnTrackIndex[lid] = 0;
        }
        theCharge = 0; // charge in the test layer
		for (int ihit = 0; ihit<nHits; ihit++){
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			double aa = (*i_aa)[ihit];
			double charge = ADC2Charge(aa); // ADC -> charge = e*Nt*GG; while Nt = dE/W = dEdX*trackL/W.
			double fd = get_dist(lid,wid,slx,inx,slz,inz);
			double dt = (*i_driftT)[ihit];
			double dd = 0;
			int status = t2d(dt,dd,fd>0);
            if ((*i_ip)[ihit]==0){ // to calculate the distance to track, only count first peaks
//                if (fabs(fd)<CELLW/2+0.5){ // along the track: distance smaller than half cell size (plus safety margin 0.5 mm)
                if (fabs(fd)<m_maxFD){ // along the track: distance smaller than half cell size (plus safety margin 0.5 mm)
                    // FIXME: should decide whether to include the boundary layers or not: slightly smaller ADC, why?
//                    if (lid>0&&lid<NLAY-1){ // don't count the last layer: guard layer
//                    if (lid>1&&lid<NLAY){ // don't count the first layer: guard layer
                    if (lid>1&&lid<NLAY-1){ // don't count the first layer and the last layer: guard layers
//                    if (lid>0&&lid<NLAY){ // count all layers
                        if (!chargeOnTrack[lid]||charge>chargeOnTrack[lid]) chargeOnTrackIndex[lid] = ihit;
                        chargeOnTrack[lid]+=charge;
                        adcsumOnTrack[lid]+=aa;
                    }
                    if (lid==m_testlayer){ // in test layer hits
                        theCharge+=charge;
                    }
				}
				if (isGood&&(*i_layerID)[ihit]==m_testlayer){
                    // aa
                    h_aaVST->Fill(dt,aa);
                    h_aaVSD->Fill(dd,aa);
                }
			}
			if ((*i_layerID)[ihit]==m_testlayer){
				int isig = -1;
				if (wid==closeWid) isig = 0;
				else if (wid==closeWid2&&nHitsT==2) isig = 1;
				if (isig>=0){
					double resi = fabs(dd)-fabs(fd);
					if (fabs(resi)<fabs(minRes[isig])){
						o_hashit[isig] = true;
						o_theFD[isig] = fd;
						o_theDD[isig] = dd;
						o_theDT[isig] = dt;
						o_theST[isig] = status;
						o_theAA[isig] = aa;
						minRes[isig] = resi;
					}
				}
			}
		}

        // Count number of layers used for charge on track
		nLayersOnTrack = 0;
		for (int lid = 1; lid<NLAY; lid++){
		    if (chargeOnTrack[lid]) nLayersOnTrack++;
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

        // get gg
        theGG = 0;
        if (isGood){
            int isig = -1;
            if (o_hashit[0]) isig = 0;
            else if (o_hashit[1]) isig = 1;
            if (isig!=-1){
                double gg = getGG(theCharge,slx,slz);
                h_ggVSX->Fill(o_theFD[isig],gg); // only record the gas gain VS x in the test layer
                if (gg>theGG) theGG = gg;
            }
            double gg = getGG(trackCharge[0]/nLayersOnTrack,slx,slz);
            h_ggall->Fill(gg);
            trackGG = gg;
        }

        otree_events->Fill();
	} // Finihsed event loop
	// get averageGG
//    averageGG = h_ggVSX->GetMean(2);
//    averageGGErr = h_ggVSX->GetRMS(2);
    averageGG = h_ggall->GetMean();
    averageGGErr = h_ggall->GetRMS();

    // loop again for filling histograms
	for ( int iEntry = m_iEntryStart ; iEntry<=m_iEntryStop; iEntry++){
	    otree_events->GetEntry(iEntry);
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
            double theDX = (nLayersOnTrack-itrunc)*CELLH*sqrt(1+slx*slx+slz*slz);
            h_dedx[itrunc]->Fill(theDE/1000/(theDX/10));
        }
        for (int i = 0; i<nHitsT; i++){
            if (!o_hashit[i]) continue;
            int status = o_theST[i];
            double fd = o_theFD[i];
            double dd = o_theDD[i];
            double dt = o_theDT[i];
            double aa = o_theAA[i];
            if (!status) continue; // out of range
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
    otree_events->Write();
//    ofile_events->Close();

	//=================================================Get bin by bin information====================================================
	TFile * ofile = new TFile(Form("%s/root/res_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer),"RECREATE");
	TTree * otree = new TTree("t","t");
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
	if (m_savehists){
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
				doFit(h_resX[ibin],3/4.,1/4.,0);
			else if (o_xmid>7)
				doFit(h_resX[ibin],1/3.,2/3.);
			else
				doFit(h_resX[ibin],1/3.,1/3.);
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
			if (m_savehists){
				h_resX[ibin]->Draw();
				canv_bin->SaveAs(Form("resX%d_%d.%s.layer%d.png",ibin,m_runNo,m_runname.Data(),m_testlayer));
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
				doFit(h_resD[ibin],1/4.,3/4.);
			else if (o_xmid>7)
				doFit(h_resD[ibin],2/3.,1/3.);
			else
				doFit(h_resD[ibin],1/3.,1/3.);
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
			if (m_savehists){
				h_resD[ibin]->Draw();
				canv_bin->SaveAs(Form("resD%d_%d.%s.layer%d.png",ibin,m_runNo,m_runname.Data(),m_testlayer));
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
	            "runNo/I","m_testlayer/I","HV/I","THR/I","gasID/I","aaCut/I",
	            "averageGG","averageGGErr",
	            "averageEffX","averageEff3sX","averageRMSX","averageResX",
	            "averageEffXErr","averageEff3sXErr","averageRMSXErr","averageResXErr",
	            "averageEffD","averageEff3sD","averageRMSD","averageResD",
	            "bestEffD","bestRMSD","bestResD","bestEffX","bestRMSX","bestResX");
	printf("==> %4d  %d   %d   %2d  %d   %2d  %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
	            m_runNo,m_testlayer,HV,THR,gasID,m_aaCut,
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
	canv_tracking->SaveAs(Form("track_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_tracking->SaveAs(Form("track_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

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
	canv_DOCA->SaveAs(Form("DOCA_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_DOCA->SaveAs(Form("DOCA_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));
	canv_DOCA->cd(1);
	h_DriftD->GetYaxis()->SetRangeUser(0,h_DriftD->GetMaximum()*1.1);
	h_DriftD->Draw();
	canv_DOCA->cd(2);
	h_DriftDb->GetYaxis()->SetRangeUser(0,h_DriftDb->GetMaximum()*1.1);
	h_DriftDb->Draw();
	canv_DOCA->SaveAs(Form("DriftD_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_DOCA->SaveAs(Form("DriftD_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TCanvas * canv_XT = new TCanvas("canv_XT","canv_XT",800,600);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_xt->Draw("COLZ");
	h_xt->GetXaxis()->SetRangeUser(minDT,maxDT);
	if (m_xtType%10==0){
		f_left0->Draw("SAME");
		f_right0->Draw("SAME");
	}
	else if (m_xtType%10==1){
		f_left->Draw("SAME");
		f_right->Draw("SAME");
	}
	else if (m_xtType%10==3){
		f_left_mid->Draw("SAME");
		f_right_mid->Draw("SAME");
	}
	else{ // fxc_4|fxm_4
		f_left_mid->SetRange(t7l,maxDT);
		f_left_mid->Draw("SAME");
		f_left->SetRange(minDT,t7l);
		f_left->Draw("SAME");
		f_right_mid->SetRange(t7r,maxDT);
		f_right_mid->Draw("SAME");
		f_right->SetRange(minDT,t7r);
		f_right->Draw("SAME");
	}
	canv_XT->SaveAs(Form("XT_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_XT->SaveAs(Form("XT_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

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
	canv_TX->SaveAs(Form("TX_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_TX->SaveAs(Form("TX_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TCanvas * canv_general = new TCanvas("canv_general","canv_general",1024,768);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h_resVSX->Draw("COLZ");
	canv_general->SaveAs(Form("resVSX_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("resVSX_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	h_resVSD->Draw("COLZ");
	canv_general->SaveAs(Form("resVSD_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("resVSD_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	h_aaVST->Draw("COLZ");
	TLine * line_aaVST = new TLine(m_tmin,m_aaCut,m_tmax,m_aaCut);
	line_aaVST->SetLineColor(kRed);
	line_aaVST->Draw("SAME");
	canv_general->SaveAs(Form("aaVST_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("aaVST_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	h_aaVSD->Draw("COLZ");
	TLine * line_aaVSD = new TLine(0,m_aaCut,m_xmax,m_aaCut);
	line_aaVSD->SetLineColor(kRed);
	line_aaVSD->Draw("SAME");
	canv_general->SaveAs(Form("aaVSD_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("aaVSD_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	h_ggVSX->Draw("COLZ");
	canv_general->SaveAs(Form("ggVSX_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("ggVSX_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

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
	canv_general->SaveAs(Form("dedx_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("dedx_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

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
	canv_general->SaveAs(Form("effx_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("effx_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

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
	canv_general->SaveAs(Form("effd_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("effd_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TGraphErrors * g_xres = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xres[0]),&(v_xxerr[0]),&(v_xreserr[0]));
	g_xres->SetName("gxres");
	g_xres->SetTitle("#sigma of residual VS DOCA");
	g_xres->GetXaxis()->SetTitle("Drift distance [mm]");
	g_xres->GetYaxis()->SetTitle("#sigma [mm]");
	g_xres->SetMarkerStyle(20);
	g_xres->SetMarkerColor(kBlack);
	g_xres->SetLineColor(kBlack);
	g_xres->Draw("APL");
	canv_general->SaveAs(Form("resx_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("resx_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TGraphErrors * g_dres = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_dres[0]),&(v_dxerr[0]),&(v_dreserr[0]));
	g_dres->SetName("gdres");
	g_dres->SetTitle("#sigma of residual VS drift distance");
	g_dres->GetXaxis()->SetTitle("Drift distance [mm]");
	g_dres->GetYaxis()->SetTitle("#sigma [mm]");
	g_dres->SetMarkerStyle(20);
	g_dres->SetMarkerColor(kBlack);
	g_dres->SetLineColor(kBlack);
	g_dres->Draw("APL");
	canv_general->SaveAs(Form("resd_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("resd_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TGraphErrors * g_xrms = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xrms[0]),&(v_xxerr[0]),&(v_xrmserr[0]));
	g_xrms->SetName("gxrms");
	g_xrms->SetTitle("RMS of residual VS DOCA");
	g_xrms->GetXaxis()->SetTitle("Drift distance [mm]");
	g_xrms->GetYaxis()->SetTitle("RMS [mm]");
	g_xrms->SetMarkerStyle(20);
	g_xrms->SetMarkerColor(kBlack);
	g_xrms->SetLineColor(kBlack);
	g_xrms->Draw("APL");
	canv_general->SaveAs(Form("rmsx_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("rmsx_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TGraphErrors * g_drms = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_drms[0]),&(v_dxerr[0]),&(v_drmserr[0]));
	g_drms->SetName("gdrms");
	g_drms->SetTitle("RMS of residual VS drift distance");
	g_drms->GetXaxis()->SetTitle("Drift distance [mm]");
	g_drms->GetYaxis()->SetTitle("RMS [mm]");
	g_drms->SetMarkerStyle(20);
	g_drms->SetMarkerColor(kBlack);
	g_drms->SetLineColor(kBlack);
	g_drms->Draw("APL");
	canv_general->SaveAs(Form("rmsd_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("rmsd_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TGraphErrors * g_xoff = new TGraphErrors(v_xx.size(),&(v_xx[0]),&(v_xoff[0]));
	g_xoff->SetName("gxoff");
	g_xoff->SetTitle("Offset VS DOCA");
	g_xoff->GetXaxis()->SetTitle("Drift distance [mm]");
	g_xoff->GetYaxis()->SetTitle("offset [mm]");
	g_xoff->SetMarkerStyle(20);
	g_xoff->SetMarkerColor(kBlack);
	g_xoff->SetLineColor(kBlack);
	g_xoff->Draw("APL");
	canv_general->SaveAs(Form("offx_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("offx_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

	TGraphErrors * g_doff = new TGraphErrors(v_dx.size(),&(v_dx[0]),&(v_doff[0]));
	g_doff->SetName("gdoff");
	g_doff->SetTitle("Offset VS drift distance");
	g_doff->GetXaxis()->SetTitle("Drift distance [mm]");
	g_doff->GetYaxis()->SetTitle("offset [mm]");
	g_doff->SetMarkerStyle(20);
	g_doff->SetMarkerColor(kBlack);
	g_doff->SetLineColor(kBlack);
	g_doff->Draw("APL");
	canv_general->SaveAs(Form("offd_%d.%s.layer%d.pdf",m_runNo,m_runname.Data(),m_testlayer));
	canv_general->SaveAs(Form("offd_%d.%s.layer%d.png",m_runNo,m_runname.Data(),m_testlayer));

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
	printf("nHits<=%d: %d (%.1f%)\n",m_nHitsMax,N_CUT1,(double)N_CUT1/N_ALL*100);
	printf("nHitsS>=%d: %d (%.1f%)\n",m_nHitsSmin,N_CUT2,(double)N_CUT2/N_ALL*100);
	printf("chi2<%.1f (pvalue>%.1f): %d (%.1f%)\n",m_maxchi2,m_minchi2p,N_CUT3,(double)N_CUT3/N_ALL*100);
	printf("|slz|<=%.1f: %d (%.1f%)\n",m_maxslz,N_CUT4,(double)N_CUT4/N_ALL*100);
	printf("DOCA<=%.1f: %d (%.1f%)\n",m_xmax,N_CUT5,(double)N_CUT5/N_ALL*100);
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

int t2d(double t, double & d, bool isRight){
	int status = 0;
	TF1 * f = 0; // body

	if (m_xtType%10==0){ // fx_0
		if (isRight) f = f_right0;
		else f = f_left0;
	}
	else if (m_xtType%10==1){ // fxc_4
		if (isRight) f = f_right;
		else f = f_left;
	}
	else if (m_xtType%10==2){ // fxc_4|fxm_4
		if (isRight){
			if (t>t7r&&t<t8r)
				f = f_right_mid;
			else
				f = f_right;
		}
		else{
			if (t>t7l&&t<t8r)
				f = f_left_mid;
			else
				f = f_left;
		}
	}
	else{ // fxm_4
		if (isRight){
			if (t<t8r)
				f = f_right_mid;
			else
				f = f_right;
		}
		else{
			if (t<t8r)
				f = f_left_mid;
			else
				f = f_left;
		}
	}

	double tRight = f->GetXmax();
	double tmin = f->GetXmin();
	double tmax = 0;
	if (m_xtType/10==0){
		tmax = tRight;
	}
	else{
		if (isRight) tmax = t8r;
		else tmax = t8l;
		if (m_tmaxSet) tmax = m_tmaxSet;
	}

	if (t<tmin){
		status = 0;
		d = 0;
	}
	else if (t>tmax){
		status = 0;
		if (t>tRight) d = f->Eval(tRight);
		else d = f->Eval(t);
	}
	else{
		status = 1;
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
	for (;binr<=m_NbinRes-3; binr++){
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

void print_usage(char * prog_name){
	fprintf(stderr,"Usage %s [options] runname\n",prog_name);
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
    fprintf(stderr,"\t -f <f>\n");
    fprintf(stderr,"\t\t Minimum number of selected hits cut set to f\n");
    fprintf(stderr,"\t -c <c>\n");
    fprintf(stderr,"\t\t Maximum chi2 cut set to c\n");
    fprintf(stderr,"\t -p <p>\n");
    fprintf(stderr,"\t\t Minimum p-value cut set to p\n");
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
    fprintf(stderr,"\t -x <x>\n");
    fprintf(stderr,"\t\t Number of bins on space axis set to x\n");
    fprintf(stderr,"\t -y <y>\n");
    fprintf(stderr,"\t\t Number of bins on resolution axis set to y\n");
    fprintf(stderr,"\t -g <g>\n");
    fprintf(stderr,"\t\t Geometry setup set to g\n");
    fprintf(stderr,"\t\t (0): normal; 1: finger\n");
}
