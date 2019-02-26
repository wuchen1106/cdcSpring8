#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TF1.h"

#include <math.h>
#include <stdlib.h>

#include "MyRuntimeParameters.hxx"
#include "Log.hxx"
#include "header.hxx"

//#define PRINT_CROSSPOINTS

#define MULTI 5 // maximum peak multiplicity of one cell to be considered in drawing

#define NLAY4DR 8
#define NWIRE4DR 8
#define NCELL4DR (NLAY4DR*NWIRE4DR)
#define IFWIRE4DR 1
#define IFLAY4DR 1

// for test 180404
//#define NLAY4DR 2
//#define NWIRE4DR 2
//#define NCELL4DR (NLAY4DR*NWIRE4DR)
//#define IFWIRE4DR 1
//#define IFLAY4DR 5

//===================About xt============================
TF1 * f_left[NCELA];
TF1 * f_right[NCELA];

double t0[NBRD];

double tdc2t(int tdc);
void cid4dr2lidwid(int cid, int & lid, int & wid);
bool lidwid2cid4dr(int lid, int wid, int & cid);

double t2x(double time, int lid, int wid, int lr, int & status);
void print_usage(char* prog_name);

int m_modulo = 0;
int m_runNo = 0;
int m_workMode = 0; // 0: h_XXX; 1: t_XXX; 2: ana_XXX; 10: h_XXX with zx; 11: t_XXX with zx
int m_testlayer = 4;
int m_thelayer = 4;
int m_thewire = -1; // negative value means just pick up the first hit, regardless of the target layer.
int m_iEntryStart = 0;
int m_iEntryStop = 9;
int m_geoSetup = 0; // 0: normal; 1: finger
TString m_xtFile = "";
TString m_runname = "";

int main(int argc, char** argv){
	//===================Get Arguments============================
    std::map<std::string, Log::ErrorPriority> namedDebugLevel;
    std::map<std::string, Log::LogPriority> namedLogLevel;
    int    opt_result;
	while((opt_result=getopt(argc,argv,"M:m:R:B:E:L:l:W:G:X:"))!=-1){
		switch(opt_result){
			/* INPUTS */
			case 'M':
			    m_modulo = atoi(optarg);
                printf("Printing modulo set to %d\n",m_modulo);
				break;
			case 'm':
			    m_workMode = atoi(optarg);
                printf("Work mode set to %d\n",m_workMode);
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
			    m_testlayer = atoi(optarg);
                printf("Test layer set to %d\n",m_testlayer);
				break;
			case 'l':
			    m_thelayer = atoi(optarg);
                printf("The layer set to %d\n",m_thelayer);
				break;
			case 'W':
			    m_thewire = atoi(optarg);
                printf("The wire set to %d\n",m_thelayer);
				break;
			case 'G':
			    m_geoSetup = atoi(optarg);
                printf("Geometry setup set to %d\n",m_geoSetup);
				break;
			case 'X':
			    m_xtFile = optarg;
                printf("Load xt curves from %s\n",m_xtFile.Data());
				break;
            case 'D':
                {
                    // Set the debug level for a named trace.
                    std::string arg(optarg);
                    std::size_t sep = arg.find("=");
                    if (sep != std::string::npos) {
                        std::string name = arg.substr(0,sep);
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

	if (argc-optind<1){
	    print_usage(argv[0]);
		return -1;
    }
    m_runname= argv[optind++];

    printf("##############%s##################\n",argv[0]);
	printf("runNo:              %d\n",m_runNo);
	printf("workMode:           %d\n",m_workMode);
	printf("test layer:         %d\n",m_testlayer);
	printf("check wire:         [%d,%d]\n",m_thelayer,m_thewire);
	printf("runname:             %s\n",m_runname.Data());
	printf("Entries:            %d~%d\n",m_iEntryStart,m_iEntryStop);
	printf("xt curves:          %s\n",m_xtFile==""?"self":m_xtFile.Data());
    printf("geoSetup:           %s\n",m_geoSetup==0?"normal scintillator":"finger scintillator");

    TString HOME=getenv("CDCS8WORKING_DIR");

    //===================Chamber Parameter============================
	double XMAX = 130; // range of x-z plane
	double ZMAX = 350; // range of x-z plane
	double yup = 645.;    // top of the chamber, for drawing x-y plane
	double ydown = 515;   // bottom of the chamber, for drawing x-y plane
    // normal scintillator
    double sciYup = 0;
    double sciYdown = 0;
    double sciHL = 0;
    double sciHW = 0;
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
    printf("sciYup      = %.3e\n",sciYup);
    printf("sciYdown    = %.3e\n",sciYdown);
    printf("sciHL       = %.3e\n",sciHL);
    printf("sciHW       = %.3e\n",sciHW);

	double iny = sciYup;     // incident y plane

	//===================Prepare Maps============================
	// map for wire position
    double  map_x[NLAY][NCEL][2];
    double  map_y[NLAY][NCEL][2];
    double  map_z[NLAY][NCEL][2];
	int     map_ch[NLAY][NCEL];
	int     map_bid[NLAY][NCEL];
	int     map_lid[NBRD][NCHS];
	int     map_wid[NBRD][NCHS];
    int     map_check[NLAY][NCEL];
    double  map_theta[NLAY][NCEL]; // rotation angle viewing from the dart plane (RO plane, z>0); positive rotation angle point to -x direction
	// mcp for cross points
    double mcp_xc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
    double mcp_zc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
    if (m_workMode/10==1){
		for(int lid = 0; lid<NLAY; lid++){
			for (int wid = 0; wid<NCEL; wid++){
				map_x[lid][wid][0] = 0;
				map_y[lid][wid][0] = 0;
				map_z[lid][wid][0] = 0;
				map_x[lid][wid][1] = 0;
				map_y[lid][wid][1] = 0;
				map_z[lid][wid][1] = 0;
				map_check[lid][wid]=0;
				if (lid <NZXP){ // z-x planes corresponding to the layerID of the lower layer counting from 1 
					for (int wjd = 0; wjd<NCEL; wjd++){
						mcp_xc[lid][wid][wjd] = 999;
						mcp_zc[lid][wid][wjd] = 999;
					}
				}
			}
		}
	}

	//===================Get Wire Position============================
	TFile * TFile_wirepos = new TFile(HOME+"/Input/wire-position.root");
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
    std::vector<double> v_wire_xc; // to get a graph of wire centers
    std::vector<double> v_wire_yc; // to get a graph of wire centers
    std::vector<double> v_wire_xhv; // to get a graph of wire centers
    std::vector<double> v_wire_yhv; // to get a graph of wire centers
    std::vector<double> v_wire_xro; // to get a graph of wire centers
    std::vector<double> v_wire_yro; // to get a graph of wire centers
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
            if (wp_lid>=1){ // ignore the first layer (dummy)
                map_check[wp_lid][wp_wid] = 1;
                v_wire_xc.push_back((wp_xhv+wp_xro)/2.);
                v_wire_yc.push_back((wp_yhv+wp_yro)/2.);
                v_wire_xhv.push_back(wp_xhv);
                v_wire_yhv.push_back(wp_yhv);
                v_wire_xro.push_back(wp_xro);
                v_wire_yro.push_back(wp_yro);
            }
            map_theta[wp_lid][wp_wid] = atan(-(map_x[wp_lid][wp_wid][0]-map_x[wp_lid][wp_wid][1])/chamberHL/2); // rotation angle viewing from the dart plane (RO plane, z>0); positive rotation angle point to -x direction
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
        mcp_xc[cp_l1][cp_w1][cp_w2] = cp_xc;
        mcp_zc[cp_l1][cp_w1][cp_w2] = cp_zc;
    }
    TFile_crosspoint->Close();

	//===================Get run info============================
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
	double npair = 17.96;
	TString gastype = "He:C_{2}H_{4}(50:50)";
	if (gasID==1){
		gastype = "He:iC_{4}H_{10}(90:10)";
		npair = 27.96;
	}
	else if (gasID==2){
		gastype = "He:CH_{4}(80:20)";
		npair = 56.10;
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
	t0[0] = t00;
	t0[1] = t01;
	std::cout<<"runNo#"<<m_runNo<<": "<<gastype<<", "<<runGr<<", "<<duration<<", "<<HV<<" V, "<<THR<<" mV, "<<durationTime<<"sec"<<std::endl;

    //===================Get XT============================
    TFile * i_xt = 0;
    if (m_xtFile==""){
        i_xt = new TFile(HOME+Form("/info/xt.%d.%s.root",m_runNo,m_runname.Data()));
    }
    else{
        i_xt = new TFile(m_xtFile);
    }
    if (!i_xt){
        fprintf(stderr,"Cannot load xt file!\n");
        return -1;
    }
    for (int i = 0; i<NCELA; i++){
        f_left[i] = (TF1*) i_xt->Get(Form("fl_%d",i/NCEL));
        f_right[i] = (TF1*) i_xt->Get(Form("fr_%d",i/NCEL));
    }

	//==================Get ADC==========================
	TChain * iChain_ADC = new TChain("tree","tree");
	iChain_ADC->Add(HOME+Form("/root/raw/run_%0.6d_built.root",m_runNo));
	int adc[NCHT][NSAM];
	int tdc[NCHT][NSAM];
	int clockNumberDriftTime[NCHT][NSAM];
	int tdcNhit[NCHT];
	iChain_ADC->SetBranchAddress("adc",adc);
	iChain_ADC->SetBranchAddress("driftTime",tdc);
	iChain_ADC->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	iChain_ADC->SetBranchAddress("tdcNhit",tdcNhit);
	if (!iChain_ADC->GetEntries()) {
	    fprintf(stderr,"Cannot find root/raw/run_%0.6d_built.root!\n",m_runNo);
	    return -1;
	}

	//==================Get Peaks==========================
	TChain * iChain_p = new TChain("t","t");
	iChain_p->Add(HOME+Form("/root/peaks/p_%d.root",m_runNo));
	double pk_aa[NCHT];
	iChain_p->SetBranchAddress("aa",pk_aa);

	//===================Get ROOT File============================
	// basic
	int triggerNumber;
	int nHits,nHitsG;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_layerID = 0;
	std::vector<double> * i_driftT = 0;
    std::vector<double> * i_dxl = 0;
    std::vector<double> * i_dxr = 0;
    std::vector<int> * i_np = 0;
    std::vector<int> * i_ip = 0;
    std::vector<int> * i_type = 0;
    std::vector<int> * i_clk = 0;
    std::vector<int> * i_width = 0;
    std::vector<int> * i_peak = 0;
    std::vector<int> * i_height = 0;
    std::vector<int> * i_mpn = 0;
    std::vector<int> * i_mpi = 0;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_aa = 0;
	// for candidates
    std::vector<double> * i_driftD[NCAND] = {0};
	int nHitsS[NCAND];
    int i_icombi[NCAND];
    int i_iselec[NCAND];
    int i_npairs[NCAND];
	double i_iinx[NCAND];
	double i_iinz[NCAND];
	double i_islx[NCAND];
	double i_islz[NCAND];
    double i_chi2x[NCAND];
    double i_chi2z[NCAND];
	double i_chi2i[NCAND];
    std::vector<double> * i_calD[NCAND] = {0};
    int i_nHitsS[NCAND];
	double i_inx[NCAND];
	double i_inz[NCAND];
	double i_slx[NCAND];
	double i_slz[NCAND];
	double i_chi2[NCAND];
    std::vector<double> * i_fitD[NCAND] = {0};
    std::vector<int> * i_sel[NCAND] = {0};
    double theRes;
    double theDD;
    int    theWid;
    int    has;
	TChain * iChain = new TChain("t","t");
    if (m_workMode%10==0){ // 0: h_XXX; 1: t_XXX
        iChain->Add(HOME+Form("/root/hits/h_%d.",m_runNo)+m_runname+".root");
    }
    else if (m_workMode%10==1){
        iChain->Add(Form("%s/root/tracks/t_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer));
    }
    else if (m_workMode%10==2){
        iChain->Add(Form("%s/root/ana/ana_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),m_testlayer));
    }
    int triggerNumberMax = 0;
    int nEntries = iChain->GetEntries();
	iChain->SetBranchAddress("triggerNumber",&triggerNumber);
	iChain->GetEntry(nEntries-1); triggerNumberMax = triggerNumber;
    if (m_workMode%10==0){
        i_driftD[0] = new std::vector<double>;
    }
    else if (m_workMode%10==1){
        iChain->SetBranchAddress("nHitsG",&nHitsG);
        iChain->SetBranchAddress("dxl",&i_dxl);
        iChain->SetBranchAddress("dxr",&i_dxr);
        for (int iCand = 0; iCand<NCAND; iCand++){
            iChain->SetBranchAddress(Form("driftD%d",iCand),&(i_driftD[iCand]));
            iChain->SetBranchAddress(Form("icom%d",iCand),&(i_icombi[iCand]));
            iChain->SetBranchAddress(Form("isel%d",iCand),&(i_iselec[iCand]));
            iChain->SetBranchAddress(Form("npairs%d",iCand),&(i_npairs[iCand]));
            iChain->SetBranchAddress(Form("islx%d",iCand),&(i_islx[iCand]));
            iChain->SetBranchAddress(Form("iinx%d",iCand),&(i_iinx[iCand]));
            iChain->SetBranchAddress(Form("islz%d",iCand),&(i_islz[iCand]));
            iChain->SetBranchAddress(Form("iinz%d",iCand),&(i_iinz[iCand]));
            iChain->SetBranchAddress(Form("chi2x%d",iCand),&(i_chi2x[iCand]));
            iChain->SetBranchAddress(Form("chi2z%d",iCand),&(i_chi2z[iCand]));
            iChain->SetBranchAddress(Form("chi2i%d",iCand),&(i_chi2i[iCand]));
            iChain->SetBranchAddress(Form("calD%d",iCand),&(i_calD[iCand]));
            iChain->SetBranchAddress(Form("nHitsS%d",iCand),&(nHitsS[iCand]));
            iChain->SetBranchAddress(Form("slx%d",iCand),&(i_slx[iCand]));
            iChain->SetBranchAddress(Form("inx%d",iCand),&(i_inx[iCand]));
            iChain->SetBranchAddress(Form("slz%d",iCand),&(i_slz[iCand]));
            iChain->SetBranchAddress(Form("inz%d",iCand),&(i_inz[iCand]));
            iChain->SetBranchAddress(Form("chi2%d",iCand),&(i_chi2[iCand]));
            iChain->SetBranchAddress(Form("fitD%d",iCand),&(i_fitD[iCand]));
            iChain->SetBranchAddress(Form("sel%d",iCand),&(i_sel[iCand]));
        }
    }
    else if (m_workMode%10==2){
        iChain->SetBranchAddress("res",&theRes);
        iChain->SetBranchAddress("theDD",&theDD);
        iChain->SetBranchAddress("theWid",&theWid);
        iChain->SetBranchAddress("has",&has);
        iChain->SetBranchAddress("nHitsG",&nHitsG);
		iChain->SetBranchAddress("driftD",&(i_driftD[0]));
		iChain->SetBranchAddress("icom",&(i_icombi[0]));
		iChain->SetBranchAddress("isel",&(i_iselec[0]));
		iChain->SetBranchAddress("npairs",&(i_npairs[0]));
		iChain->SetBranchAddress("islx",&(i_islx[0]));
		iChain->SetBranchAddress("iinx",&(i_iinx[0]));
		iChain->SetBranchAddress("islz",&(i_islz[0]));
		iChain->SetBranchAddress("iinz",&(i_iinz[0]));
		iChain->SetBranchAddress("chi2x",&(i_chi2x[0]));
		iChain->SetBranchAddress("chi2z",&(i_chi2z[0]));
		iChain->SetBranchAddress("chi2i",&(i_chi2i[0]));
		iChain->SetBranchAddress("calD",&(i_calD[0]));
		iChain->SetBranchAddress("nHitsS",&(nHitsS[0]));
		iChain->SetBranchAddress("slx",&(i_slx[0]));
		iChain->SetBranchAddress("inx",&(i_inx[0]));
		iChain->SetBranchAddress("slz",&(i_slz[0]));
		iChain->SetBranchAddress("inz",&(i_inz[0]));
		iChain->SetBranchAddress("chi2",&(i_chi2[0]));
		iChain->SetBranchAddress("fitD",&(i_fitD[0]));
		iChain->SetBranchAddress("sel",&(i_sel[0]));
    }
    else{
        fprintf(stderr,"workMode %d is not supported! please chose from 0,1\n",m_workMode);
        return -1;
    }
	iChain->SetBranchAddress("nHits",&nHits);
	iChain->SetBranchAddress("wireID",&i_wireID);
	iChain->SetBranchAddress("layerID",&i_layerID);
	iChain->SetBranchAddress("driftT",&i_driftT);
    iChain->SetBranchAddress("type",&i_type);
    iChain->SetBranchAddress("np",&i_np);
    iChain->SetBranchAddress("ip",&i_ip);
    iChain->SetBranchAddress("clk",&i_clk);
    iChain->SetBranchAddress("width",&i_width);
    iChain->SetBranchAddress("peak",&i_peak);
    iChain->SetBranchAddress("height",&i_height);
    iChain->SetBranchAddress("mpn",&i_mpn);
    iChain->SetBranchAddress("mpi",&i_mpi);
    iChain->SetBranchAddress("sum",&i_sum);
    iChain->SetBranchAddress("aa",&i_aa);

	//==================Prepare canvas for drawing==========================
	// run summary
	printf("Preparing canvas...\n");
	TLatex * text_runsum = new TLatex();
	text_runsum->SetTextSize(0.02);
	text_runsum->SetText(0.05,0.98,Form("run#%d ",m_runNo)+gastype+Form(", %d V,%d mV, Grade#%d",HV,THR,runGr)+", "+duration+Form(", %d events, Eff_{daq} = %2.2lf%%, Rate_{tri} = %1.1lfkHz",nEntries,((double)nEntries)/(triggerNumberMax+1)*100,(triggerNumberMax+1)/durationTime/1000));

	TLatex * text_title = new TLatex(0,0,"");
	text_title->SetTextSize(0.02);
	//Prepare the Canvas for waveforms (by board)
	TCanvas * ca_WF[NBRD];
	TPad * pad_WF[NBRD][NCHS];
	for (int bid = 0; bid<NBRD; bid++){
        ca_WF[bid] = new TCanvas(Form("ca_WF_%d",bid),"ca_WF",896,896);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        for (int i = 0; i<8; i++){
            for (int j = 0; j<6; j++){
                int index = j*8+i;
                pad_WF[bid][index] = new TPad(Form("pad_%d_%d_%d",bid,i,j),Form("pad_%d_%d_%d",bid,i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
                pad_WF[bid][index]->Draw();
                pad_WF[bid][index]->SetGridx(1);
                pad_WF[bid][index]->SetGridy(1);
            }
        }
    }
	//Prepare the Canvas for waveforms (by layer)
	TCanvas * ca_WFL;
	TPad * pad_WFL[NCELL4DR];
	ca_WFL = new TCanvas("ca_WFL","ca_WFL",1024,1024.*NLAY4DR/NWIRE4DR);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	for (int i = 0; i<NWIRE4DR; i++){
		for (int j = 0; j<NLAY4DR; j++){
			int index = j*NWIRE4DR+i;
			pad_WFL[index] = new TPad(Form("padL_%d_%d",i,j),Form("padL_%d_%d",i,j),1./NWIRE4DR*i,0.95/NLAY4DR*(NLAY4DR-1-j),1./NWIRE4DR*(i+1),0.95/NLAY4DR*(NLAY4DR-j));
			pad_WFL[index]->Draw();
//			pad_WFL[index]->SetGridx(1);
//			pad_WFL[index]->SetGridy(1);
		}
	}
	//Prepare the Canvas for x-y plane and target chanel ADC
	TCanvas * ca_xyADC = new TCanvas("ca_xyADC","ca_xyADC",896,1024);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
	TPad * pad_xyADC[2];
	for (int ipad = 0; ipad<2; ipad++){
        pad_xyADC[ipad] = new TPad(Form("pad_xyADC_%i",ipad),"pad_xyADC",0,ipad?0:0.2,1,ipad?0.2:0.96);
        pad_xyADC[ipad]->Draw();
        pad_xyADC[ipad]->SetGridx(1);
        pad_xyADC[ipad]->SetGridy(1);
	}
	//Prepare the Canvas for z-x planes
	TCanvas* ca_zx_all = new TCanvas("ca_zx_all","ca_zx_all",1024,768);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
	TCanvas* ca_zx[NZXP]; // z-x planes corresponding to the layerID of the lower layer counting from 1
	for (int izx = 1; izx<NZXP; izx++){
	    ca_zx[izx] = new TCanvas(Form("ca_zx_%d",izx),"ca_zx",1024,768);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gPad->SetGridx(1);
        gPad->SetGridy(1);
    }
    // Prepare colors to be used for each wire
    int color[NCEL];
    int icolor = 0;
    color[icolor++] = kBlack;
    color[icolor++] = kMagenta+2;
    color[icolor++] = kMagenta;
    color[icolor++] = kBlue;
    color[icolor++] = kBlue-2;
    color[icolor++] = kCyan;
    color[icolor++] = kGreen;
    color[icolor++] = kYellow;
    color[icolor++] = kOrange;
    color[icolor++] = kRed;
    color[icolor++] = kRed+2;

	//==================Prepare drawing objects==========================
	//Prepare for ADC
	TGraph * gr_waveForm[NBRD][NCHS] = {0};
	int vSample[NSAM];
	for (int i=0; i<NSAM; i++){
		vSample[i] = i;
	}
	TLatex *textWF[NBRD][NCHS];
	TLatex *textTDC[NBRD][NCHS][NSAM];
	TMarker *markerTDC[NBRD][NCHS][NSAM];
    for (int bid = 0; bid<NBRD; bid++){
        for (int ch = 0; ch<NCHS; ch++){
            textWF[bid][ch] = new TLatex(0,0,"");
            textWF[bid][ch]->SetTextSize(0.04);
            for (int sid=0; sid<NSAM; sid++) {
                textTDC[bid][ch][sid] = new TLatex(0,0,"");
                textTDC[bid][ch][sid]->SetTextSize(0.04);
                markerTDC[bid][ch][sid] = new TMarker(0,0,20);
                markerTDC[bid][ch][sid]->SetMarkerSize(0.3);
            }
        } 
	}
	//Prepare for hit circles on x-y plane
	TLatex * text_xyhit[NLAY][NCEL];
	TEllipse * circle_driftD[NLAY][NCEL][MULTI];  // value from track finding/fitting;
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
            text_xyhit[lid][wid] = new TLatex(0,0,"");
            text_xyhit[lid][wid]->SetTextSize(0.02);
		    for (int ip = 0; ip<MULTI; ip++){
                circle_driftD[lid][wid][ip] = new TEllipse(0,0,1,1);
                circle_driftD[lid][wid][ip]->SetFillStyle(0);
            }
		}
	}
	//Prepare two tracks on x-y plane: l_track: track after fitting; l_itrack: track before fitting, after track finding.
	TLine * l_track = new TLine();
	l_track->SetLineColor(kRed);
	l_track->SetY1(yup);
	l_track->SetY2(ydown);
	TLine * l_itrack = new TLine();
	l_itrack->SetLineColor(kBlue);
	l_itrack->SetY1(yup);
	l_itrack->SetY2(ydown);
    // Background graph for x-y plane
    TGraph * gr_wireCenter = new TGraph(v_wire_yc.size(),&(v_wire_xc[0]),&(v_wire_yc[0]));
    gr_wireCenter->SetMarkerStyle(20);
    gr_wireCenter->SetMarkerSize(0.55);
    gr_wireCenter->GetXaxis()->SetTitle("x [mm]");
    gr_wireCenter->GetYaxis()->SetTitle("y [mm]");
    // Prepare driftT lines on z-x planes
    TLine * l_zx[NLAY][NCEL][MULTI][2];
    if (m_workMode/10==1){
		for (int lid = 0; lid<NLAY; lid++){
			for (int wid = 0; wid<NCEL; wid++){
				for (int ip = 0; ip<MULTI; ip++){
					for (int ilr = 0; ilr<2; ilr++){
						l_zx[lid][wid][ip][ilr] = new TLine();
						l_zx[lid][wid][ip][ilr]->SetLineColor(color[wid]);
						l_zx[lid][wid][ip][ilr]->SetLineStyle(3);
						l_zx[lid][wid][ip][ilr]->SetX1(-chamberHL);
						l_zx[lid][wid][ip][ilr]->SetX2(chamberHL);
					}
				}
			}
		}
    }
    // Prepare cross points of driftT lines on z-x planes
    TMarker * point_cross_zx[NZXP][NCEL][NCEL][MULTI][MULTI][4];
    TText * text_cross_zx[NZXP][NCEL][NCEL][MULTI][MULTI][4];
    if (m_workMode/10==1){
		for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
			for (int wid = 0; wid<NCEL; wid++){
				for (int wjd = 0; wjd<NCEL; wjd++){
					for (int ip = 0; ip<MULTI; ip++){
						for (int jp = 0; jp<MULTI; jp++){
							for (int icombi = 0; icombi<4; icombi++){
								if (m_workMode%10>0&&icombi==3) // reverse izx when icombi is 0 or 1, reverse izx+1 when icombi is 0 or 2. Keep them unchanged only when icombi is 3
									point_cross_zx[izx][wid][wjd][ip][jp][icombi] = new TMarker(0,0,20);
								else
									point_cross_zx[izx][wid][wjd][ip][jp][icombi] = new TMarker(0,0,4);
								point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(kBlack);
								point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerSize(0.55);
								text_cross_zx[izx][wid][wjd][ip][jp][icombi] = new TText(0,0,Form("%d",izx));
								text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextColor(color[izx]);
								text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextSize(0.01);
							}
						}
					}
				}
			}
		}
    }
    // Prepare track points on z-x planes
    TMarker * point_track_zx[NZXP];
    TMarker * point_itrack_zx[NZXP];
    if (m_workMode%10>0&&m_workMode/10==1){
        for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
            point_track_zx[izx] = new TMarker(0,0,20);
            point_track_zx[izx]->SetMarkerColor(kRed);
            point_track_zx[izx]->SetMarkerSize(0.5);
            point_itrack_zx[izx] = new TMarker(0,0,4);
            point_itrack_zx[izx]->SetMarkerColor(kBlue);
            point_itrack_zx[izx]->SetMarkerSize(0.6);
        }
    }
    // Prepare wires in each layer on z-x planes
    TGraph * gr_wire[NLAY][NCEL];
    if (m_workMode/10==1){
		for (int lid = 1; lid<NLAY; lid++){
			for (int wid = 0; wid<NCEL; wid++){
				if (map_check[lid][wid]==1){
					gr_wire[lid][wid] = new TGraph(2,map_z[lid][wid],map_x[lid][wid]);
					gr_wire[lid][wid]->SetLineColor(color[wid]);
					gr_wire[lid][wid]->SetMarkerColor(color[wid]);
				}
			}
		}
	}
    // Prepare texts and markers for wire on each z-x planes
	TMarker * point_cross_wire[NZXP][NCEL][NCEL];
	TLatex * text_cross_wire[NZXP][NCEL][NCEL];
	TLatex * text_cr_1[NZXP];
	TLatex * text_cr_1l[NZXP][NCEL];
	TLatex * text_cr_2[NZXP];
	TLatex * text_cr_2r[NZXP][NCEL];
    TGraph * gr_all[NZXP];
    double ar_cr_x[2] = {-ZMAX,ZMAX};
    double ar_cr_y[2] = {-XMAX,XMAX};
    gr_all[0] = new TGraph(2,ar_cr_x,ar_cr_y);
    gr_all[0]->SetTitle(Form("Cross points of all layers"));
    gr_all[0]->SetMarkerColor(kWhite);
    gr_all[0]->SetMarkerSize(0.1);
    gr_all[0]->GetXaxis()->SetRangeUser(-ZMAX,ZMAX);
    gr_all[0]->GetYaxis()->SetRangeUser(-XMAX,XMAX);
    gr_all[0]->GetXaxis()->SetTitle("z [mm]");
    gr_all[0]->GetYaxis()->SetTitle("x [mm]");
    if (m_workMode/10==1){
		for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
			// Background graph for z-x planes
			gr_all[izx] = new TGraph(2,ar_cr_x,ar_cr_y);
			gr_all[izx]->SetTitle(Form("layer #%d and layer #%d",izx,izx+1));
			gr_all[izx]->SetMarkerColor(kWhite);
			gr_all[izx]->SetMarkerSize(0.1);
			gr_all[izx]->GetXaxis()->SetRangeUser(-ZMAX,ZMAX);
			gr_all[izx]->GetYaxis()->SetRangeUser(-XMAX,XMAX);
			gr_all[izx]->GetXaxis()->SetTitle("z [mm]");
			gr_all[izx]->GetYaxis()->SetTitle("x [mm]");
			text_cr_1[izx] = new TLatex(map_z[izx][0][0]-30,XMAX-20,Form("Layer %d",izx));
			text_cr_1[izx]->SetTextSize(0.025);
			text_cr_2[izx] = new TLatex(map_z[izx+1][0][1]-30,XMAX-20,Form("Layer %d",izx+1));
			text_cr_2[izx]->SetTextSize(0.025);
			for (int wid = 0; wid<NCEL; wid++){
				if (map_check[izx][wid]==1){
					text_cr_1l[izx][wid] = new TLatex(map_z[izx][wid][0]-10,map_x[izx][wid][0],Form("%d",wid));
					text_cr_1l[izx][wid]->SetTextColor(color[wid]);
					text_cr_1l[izx][wid]->SetTextSize(0.02);
				}
				if (map_check[izx+1][wid]==1){
					text_cr_2r[izx][wid] = new TLatex(map_z[izx+1][wid][1]+10,map_x[izx+1][wid][1],Form("%d",wid));
					text_cr_2r[izx][wid]->SetTextColor(color[wid]);
					text_cr_2r[izx][wid]->SetTextSize(0.02);
				}
#ifdef PRINT_CROSSPOINTS
				for (int wjd = 0; wjd < NCEL; wjd++){
					if (fabs(mcp_zc[izx][wid][wjd])>300) continue;
					point_cross_wire[izx][wid][wjd] = new TMarker(mcp_zc[izx][wid][wjd],mcp_xc[izx][wid][wjd],20);
					point_cross_wire[izx][wid][wjd]->SetMarkerColor(color[wid]);
					point_cross_wire[izx][wid][wjd]->SetMarkerSize(0.4);
					text_cross_wire[izx][wid][wjd] = new TLatex(mcp_zc[izx][wid][wjd],mcp_xc[izx][wid][wjd]+5,Form("%d,%d",wid,wjd));
					text_cross_wire[izx][wid][wjd]->SetTextColor(color[wid]);
					text_cross_wire[izx][wid][wjd]->SetTextSize(0.02);
				}
#endif
			}
		}
	}

	//===================Prepare counters============================
    //  and to record the position of hits in each layer
    int    nHits_layer[NLAY];
    int    nHits_cell[NLAY][NCEL];
    int    iHit_cell[NLAY][NCEL][MULTI];
    // to record each dd & fd
    double dd_cell[NLAY][NCEL][MULTI];
    double fd_cell[NLAY][NCEL][MULTI];
    // to record track position on each layer
    double y_cell[NLAY][NCEL];
    for (int lid = 1; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            if (map_check[lid][wid]==1){
                y_cell[lid][wid] = (map_y[lid][wid][1]+map_y[lid][wid][0])/2.;
            }
        }
    }

	//===================Loop in Events============================
	TString prefix = "";
	int the_bid = -1;
	int the_ch = -1;
	int the_ihit = -1;
	printf("the_ihit @ %p\n",(void*)(&the_ihit));
	printf("Looping in events %d~%d\n",m_iEntryStart,m_iEntryStop);
	for ( int iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
		iChain->GetEntry(iEntry);
		if (nHits<=0){
			printf("No hits in event %d, continue\n",iEntry);
			continue;
		}

        // Find the target channel
		the_bid = -1;
		the_ch = -1;
		the_ihit = -1;
        if (m_thelayer>=0&&m_thewire>=0){
            for (int ihit; ihit<i_driftT->size(); ihit++){
                int lid = (*i_layerID)[ihit];
                int wid = (*i_wireID)[ihit];
                if (lid==m_thelayer&&wid==m_thewire){
                    the_ch = map_ch[lid][wid];
                    the_bid = map_bid[lid][wid];
                    the_ihit = ihit;
                    break;
                }
            }
        }
		if (the_bid==-1||the_ch==-1){ // Didn't find the wire? just pick up the first hit, regardless of the target layer
			printf("Didn't find the wire to check! Will pick the first one\n");
            int lid = (*i_layerID)[0];
            int wid = (*i_wireID)[0];
            the_bid = map_bid[lid][wid];
            the_ch = map_ch[lid][wid];
            the_ihit = 0;
		}

        iChain_ADC->GetEntry(iEntry);
        iChain_p->GetEntry(iEntry);

		// set prefix
		if (m_workMode%10>0){
            if (nHitsG<=6) prefix = "incom.";
            else if (nHitsG==7 ) prefix = "single.";
            else if (nHitsG==8 ) prefix = "single.";
            else if (nHitsG==9 ) prefix = "single.";
            else if (nHitsG==10) prefix = "n10.";
            else if (nHitsG==11) prefix = "n11.";
            else if (nHitsG==12) prefix = "n12.";
            else if (nHitsG>=13) prefix = "multi.";
        }

        // Reset counters
		for (int lid = 0; lid<NLAY; lid++){
		    nHits_layer[lid] = 0;
			for (int wid = 0; wid<NCEL; wid++){
                nHits_cell[lid][wid] = 0;
			}
		}
		if (m_workMode%10==0) i_driftD[0]->clear();
        else if (m_workMode%10==1){
            for (int iCand = 0; iCand<(m_workMode%10==1?NCAND:1); iCand++){
                if (nHitsS[iCand]==0){// bad fitting, driftD nonsense.
                    i_driftD[iCand]->resize(nHits);
                }
            }
        }
		// count
        for (int ihit = 0; ihit<nHits; ihit++){
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int type = (*i_type)[ihit];
            if (m_workMode%10==0){// prepare dd and type
                double dt = (*i_driftT)[ihit];
                int dd_status = 0;
                double dd = t2x(dt,lid,wid,0,dd_status);
                if (dd_status==1) type+=20;
                else if (dd_status==-1) type+=10;
                (*i_type)[ihit] = type;
                i_driftD[0]->push_back(dd);
            }
			else if (m_workMode%10==1){
				for (int iCand = 0; iCand<(m_workMode%10==1?NCAND:1); iCand++){
					if (nHitsS[iCand]==0){// bad fitting, driftD nonsense.
						(*i_driftD[iCand])[ihit] = (*i_dxr)[ihit];
					}
				}
			}
            if (type<=3){ // good hit first
                if (nHits_cell[lid][wid]<MULTI){ // don't count the rest
                    iHit_cell[lid][wid][nHits_cell[lid][wid]] = ihit;
                    nHits_cell[lid][wid]++;
                }
                nHits_layer[lid]++;
            }
        }
        for (int ihit = 0; ihit<nHits; ihit++){// then bad hit for t_XXX
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int type = (*i_type)[ihit];
            if (type>3){
                if (nHits_cell[lid][wid]<MULTI){ // don't count the rest
                    iHit_cell[lid][wid][nHits_cell[lid][wid]] = ihit;
                    nHits_cell[lid][wid]++;
                }
                nHits_layer[lid]++;
            }
        }

        //===================Draw waveforms============================
        // clear objects from the previous event
	    for (int bid = 0; bid<NBRD; bid++){
	        for (int ch = 0; ch<NCHS; ch++){
                if (gr_waveForm[bid][ch]) delete gr_waveForm[bid][ch]; gr_waveForm[bid][ch] = 0;
	        }
	    }
        // Draw waveforms
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->cd();
			text_title->SetText(0.05,0.96,Form("Entry %d, Trigger Number %d, %d TDCs, %d chosen (first TDC in channel with ADC sum > %.0f)",iEntry,triggerNumber,nHits,nHitsG,aacut));
            text_title->Draw();
            text_runsum->Draw();
            for (int ch = 0; ch<NCHS; ch++){
                // get waveform and draw
                int chg = bid*NCHS+ch;
                int lid = map_lid[bid][ch];
                int wid = map_wid[bid][ch];
                pad_WF[bid][ch]->cd();
                gr_waveForm[bid][ch] = new TGraph(NSAM,vSample,adc[chg]);
				gr_waveForm[bid][ch]->SetTitle(Form("Layer %d Wire %d Board %d Channel %d ",lid,wid,bid,ch));
                gr_waveForm[bid][ch]->GetXaxis()->SetRangeUser(0,NSAM-1);
                gr_waveForm[bid][ch]->GetYaxis()->SetRangeUser(MIN_ADC,MAX_ADC);
                gr_waveForm[bid][ch]->GetXaxis()->SetTitle("Sample Index");
                gr_waveForm[bid][ch]->GetYaxis()->SetTitle("ADC");
                gr_waveForm[bid][ch]->SetMarkerStyle(20);
                gr_waveForm[bid][ch]->SetMarkerSize(0.3);
                gr_waveForm[bid][ch]->Draw("APL");
                // set title text
                textWF[bid][ch]->SetText(1,MAX_ADC-50,Form("%d peaks, ADC sum %.0lf",tdcNhit[chg],pk_aa[chg]));
                if (tdcNhit[chg]>0)
                    textWF[bid][ch]->SetTextColor(kBlue); // has hit
                else
                    textWF[bid][ch]->SetTextColor(kGreen); // no hit
                textWF[bid][ch]->Draw();
            }
        }
        for (int ihit = 0; ihit<nHits; ihit++){ // update with hit information
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int bid = map_bid[lid][wid];
            int ch = map_ch[lid][wid];
            int chg = bid*NCHS+ch;
            int type = (*i_type)[ihit];
            int ip = (*i_ip)[ihit];
            int clk = (*i_clk)[ihit];
            int height = (*i_height)[ihit];
            ca_WF[bid]->cd();
            pad_WF[bid][ch]->cd();
            markerTDC[bid][ch][ip]->SetX(clk);
            markerTDC[bid][ch][ip]->SetY(height);
            textTDC[bid][ch][ip]->SetText(clk,height,Form("%d",(int)(tdc[chg][ip])));
            if (type<=3){ // good hit
                markerTDC[bid][ch][ip]->SetMarkerColor(kRed);
                textTDC[bid][ch][ip]->SetTextColor(kRed);
                textWF[bid][ch]->SetTextColor(kRed);
                textWF[bid][ch]->Draw();
            }
            else{ // bad hit
                markerTDC[bid][ch][ip]->SetMarkerColor(14);
                textTDC[bid][ch][ip]->SetTextColor(14);
            }
            markerTDC[bid][ch][ip]->Draw();
            textTDC[bid][ch][ip]->Draw();
        }
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->SaveAs(prefix+Form("wf.%d.b%d.pdf",iEntry,bid));
            ca_WF[bid]->SaveAs(prefix+Form("wf.%d.b%d.png",iEntry,bid));
        }

        //===================Draw waveforms by layer============================
		ca_WFL->cd();
		text_title->Draw();
		text_runsum->Draw();
		for (int cell = 0; cell<NCELL4DR; cell++){
			// get waveform and draw
			int lid = 0;
			int wid = 0;
			cid4dr2lidwid(cell,lid,wid);
			int ch = map_ch[lid][wid];
			int bid = map_bid[lid][wid];
			pad_WFL[cell]->cd();
			gr_waveForm[bid][ch]->GetYaxis()->UnZoom();
			gr_waveForm[bid][ch]->Draw("APL");
			double maxadc  = gr_waveForm[bid][ch]->GetYaxis()->GetXmax();
			double minadc  = gr_waveForm[bid][ch]->GetYaxis()->GetXmin();
			double rangeadc = maxadc-minadc;
			textWF[bid][ch]->SetY(minadc+rangeadc/20);
			textWF[bid][ch]->Draw();
		}
        for (int ihit = 0; ihit<nHits; ihit++){ // update with hit information
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int cell = 0;
            bool notincluded = lidwid2cid4dr(lid,wid,cell);
            if (notincluded) continue;
            if (cell<0||cell>=NCELL4DR) continue;
            int bid = map_bid[lid][wid];
            int ch = map_ch[lid][wid];
            int chg = bid*NCHS+ch;
            int ip = (*i_ip)[ihit];
            pad_WFL[cell]->cd();
            int type = (*i_type)[ihit];
            if (type<=3){
				markerTDC[bid][ch][ip]->SetMarkerColor(kRed);
				textTDC[bid][ch][ip]->SetTextColor(kRed);
			}
			else{
				markerTDC[bid][ch][ip]->SetMarkerColor(kBlue);
				textTDC[bid][ch][ip]->SetTextColor(kBlue);
			}
			// for test 180404
			//markerTDC[bid][ch][ip]->SetMarkerSize(1);
            markerTDC[bid][ch][ip]->Draw();
            textTDC[bid][ch][ip]->Draw();
        }
		ca_WFL->SaveAs(prefix+Form("wfl.%d.pdf",iEntry));
		ca_WFL->SaveAs(prefix+Form("wfl.%d.png",iEntry));

        //===================Draw target waveform in xyADC plot============================
        // draw the target channel ADC
        ca_xyADC->cd();
        text_runsum->Draw();
		pad_xyADC[1]->cd();
		gr_waveForm[the_bid][the_ch]->Draw("APL");
		int the_chg = the_bid*NCHS+the_ch;
		for (int ip=0; ip<tdcNhit[the_chg]; ip++) {
            double sum = (*i_sum)[the_ihit+ip];
            double dd = (*i_driftD[0])[the_ihit+ip];
            double dt = (*i_driftT)[the_ihit+ip];
            int clk = (*i_clk)[the_ihit+ip];
            int height = (*i_height)[the_ihit+ip];
            textTDC[the_bid][the_ch][ip]->SetText(clk,height+(0.5-ip%2)*50,Form("s%.0f,%.0f,%.2f",sum,dt,dd));
			textTDC[the_bid][the_ch][ip]->Draw();
			markerTDC[the_bid][the_ch][ip]->Draw();
		}

        for (int iCand = 0; iCand<(m_workMode%10==1?NCAND:1); iCand++){
            //===================Draw XY in xyADC plot============================
            // update wire position
            for (int iwire = 0; iwire<v_wire_xc.size(); iwire++){
                double y = v_wire_yc[iwire];
                double x;
                if (m_workMode%10>0){ // according to z-y relation from tracking
                    double z = i_slz[iCand]*(y-iny)+i_inz[iCand];
                    x = ((chamberHL-z)*v_wire_xhv[iwire]+(chamberHL+z)*v_wire_xro[iwire])/chamberHL/2;
                    y = ((chamberHL-z)*v_wire_yhv[iwire]+(chamberHL+z)*v_wire_yro[iwire])/chamberHL/2;
                    z = i_slz[iCand]*(y-iny)+i_inz[iCand];
                    x = ((chamberHL-z)*v_wire_xhv[iwire]+(chamberHL+z)*v_wire_xro[iwire])/chamberHL/2;
                    y = ((chamberHL-z)*v_wire_yhv[iwire]+(chamberHL+z)*v_wire_yro[iwire])/chamberHL/2;
                }
                else{
                    x = v_wire_xc[iwire];
                }
                gr_wireCenter->SetPoint(iwire,x,y);
            }
            // Draw the background graph for x-y plane
            pad_xyADC[0]->cd();
            if (m_workMode%10>0)
                gr_wireCenter->SetTitle(Form("Ent %d, nHitsG (%d)%d(%d), icom %d, isel %d, sl_{z}: %.2e->%.2e, #chi^{2}: %.2e->%.2e",iEntry,nHits,nHitsG,nHitsS[iCand],i_icombi[iCand],i_iselec[iCand],i_islz[iCand],i_slz[iCand],i_chi2i[iCand],i_chi2[iCand]));
            else
                gr_wireCenter->SetTitle(Form("Entry %d nHits = %d",iEntry,nHits));
            gr_wireCenter->Draw("AP");
            // Draw the hit circles
            for (int lid = 0; lid<NLAY; lid++){
                for (int wid = 0; wid<NCEL; wid++){
                    // get wire position
                    double wxro = map_x[lid][wid][1];
                    double wyro = map_y[lid][wid][1];
                    double wzro = chamberHL;
                    double wxhv = map_x[lid][wid][0];
                    double wyhv = map_y[lid][wid][0];
                    double wzhv = -chamberHL;
                    double wy = (wyro+wyhv)/2.;
                    double wx = (wxro+wxhv)/2.;
                    if (m_workMode%10>0){
                        // correct wx wy wz according to the track position
                        double wz = i_slz[iCand]*(wy-iny)+i_inz[iCand];
                        wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
                        wz = i_slz[iCand]*(wy-iny)+i_inz[iCand];
                        wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
                    }
                    double resmin = 1e9;
                    double thefitd = 0;
                    for (int ip = 0; ip<nHits_cell[lid][wid]; ip++){
                        int ihit = iHit_cell[lid][wid][ip];
                        int type = (*i_type)[ihit];
                        // Get hit information
                        double fitd = 0;
                        double dd = (*i_driftD[iCand])[ihit];
                        if (m_workMode%10>0){
                            fitd = (*i_fitD[iCand])[ihit]; // !!! make sure it's mm
                        }
                        // set the hit circles
                        circle_driftD[lid][wid][ip]->SetX1(wx);
                        circle_driftD[lid][wid][ip]->SetY1(wy);
                        circle_driftD[lid][wid][ip]->SetR1(fabs(dd));
                        circle_driftD[lid][wid][ip]->SetR2(fabs(dd));
                        if (type>3){ // bad hit
                            circle_driftD[lid][wid][ip]->SetLineStyle(2);
                            circle_driftD[lid][wid][ip]->SetLineColor(14);
                        }
                        else{
                            circle_driftD[lid][wid][ip]->SetLineStyle(1);
                            if (m_workMode%10>0&&(*i_sel[iCand])[ihit]) // used for fitting
                                circle_driftD[lid][wid][ip]->SetLineColor(kRed);
                            else
                                circle_driftD[lid][wid][ip]->SetLineColor(kOrange);
                        }
                        // set the text of each hit
                        circle_driftD[lid][wid][ip]->Draw(); // from track fitting/finding

                        // get the min residual
                        if (fabs(fitd-dd)<fabs(resmin)&&(m_workMode%10>0&&type<=3)){ // only print when t_XXX and good hit
                            resmin = fitd-dd;
                            thefitd = fitd;
                        }

                        // Get information for cross points
                        dd_cell[lid][wid][ip] = dd;
                        fd_cell[lid][wid][ip] = fitd;
                        double delta = dd*chamberHL*2/sqrt(chamberHL*chamberHL*4+(map_x[lid][wid][0]-map_x[lid][wid][1])*(map_x[lid][wid][0]-map_x[lid][wid][1])); // correction for x position
                        if (m_workMode/10==1){
							for (int ilr = 0; ilr<2; ilr++){
								l_zx[lid][wid][ip][ilr]->SetY1(wxhv+(ilr?delta:-delta));
								l_zx[lid][wid][ip][ilr]->SetY2(wxro+(ilr?delta:-delta));
							}
						}
                    }
                    if (m_workMode%10==2&&wid==theWid&&has==1&&lid==m_testlayer){
                    	resmin = theRes;
                    	thefitd = theDD+theRes;
                    }
                    if (resmin<1e9){ // found a hit in this wire
                        text_xyhit[lid][wid]->SetText(wx,wy,Form("%d,%.3f,%.3f",wid,thefitd,resmin));// draw the one with the smallest res
                        text_xyhit[lid][wid]->Draw();
                    }
                    y_cell[lid][wid] = wy;
                }
            }
            // draw the tracks on the x-y plane
            double xdown  = i_inx[iCand] + (ydown-iny)*i_slx[iCand];
            double xup    = i_inx[iCand] + (yup-iny)*i_slx[iCand];
            double xdowni = i_iinx[iCand] + (ydown-iny)*i_islx[iCand];
            double xupi   = i_iinx[iCand] + (yup-iny)*i_islx[iCand];
            if (m_workMode%10>0){
                l_itrack->SetX1(xupi);
                l_itrack->SetX2(xdowni);
                l_itrack->Draw();
                l_track->SetX1(xup);
                l_track->SetX2(xdown);
                l_track->Draw();
            }
            if (m_workMode%10>0){
                ca_xyADC->SaveAs(prefix+Form("xyADC.%d.i%d.pdf",iEntry,iCand));
                ca_xyADC->SaveAs(prefix+Form("xyADC.%d.i%d.png",iEntry,iCand));
            }
            else{
                ca_xyADC->SaveAs(prefix+Form("xyADC.%d.pdf",iEntry));
                ca_xyADC->SaveAs(prefix+Form("xyADC.%d.png",iEntry));
            }

            //===================Draw ZX plots============================
            // draw the z-x planes
            if (m_workMode/10==1){
				for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
					ca_zx[izx]->cd();
					// draw Background graph for z-x planes
					gr_all[izx]->Draw("AP");
					text_cr_1[izx]->Draw();
					text_cr_2[izx]->Draw();
					for (int wid = 0; wid<NCEL; wid++){
						if (map_check[izx][wid]==1){
							gr_wire[izx][wid]->Draw("PLSAME");
							text_cr_1l[izx][wid]->Draw();
						}
						if (map_check[izx+1][wid]==1){
							gr_wire[izx+1][wid]->Draw("PLSAME");
							text_cr_2r[izx][wid]->Draw();
						}
#ifdef PRINT_CROSSPOINTS
						for (int wjd = 0; wjd < NCEL; wjd++){
							if (!point_cross_wire[izx][wid][wjd]
									||!text_cross_wire[izx][wid][wjd]) continue;
							point_cross_wire[izx][wid][wjd]->Draw();
							text_cross_wire[izx][wid][wjd]->Draw();
						}
#endif
					}
					// draw the track and driftT lines on the z-x planes
					for (int ilr = 0; ilr<2; ilr++){
						for (int wid = 0; wid<NCEL; wid++){
							for (int ip = 0; ip<nHits_cell[izx][wid]; ip++){
								int ihit = iHit_cell[izx][wid][ip];
								int itype = (*i_type)[ihit];
								if (itype<=3)
									l_zx[izx][wid][ip][ilr]->SetLineColor(color[wid]);
								else 
									l_zx[izx][wid][ip][ilr]->SetLineColor(14);
								l_zx[izx][wid][ip][ilr]->Draw();
							}
							for (int ip = 0; ip<nHits_cell[izx+1][wid]; ip++){
								int ihit = iHit_cell[izx+1][wid][ip];
								int itype = (*i_type)[ihit];
								if (itype<=3)
									l_zx[izx+1][wid][ip][ilr]->SetLineColor(color[wid]);
								else 
									l_zx[izx+1][wid][ip][ilr]->SetLineColor(14);
								l_zx[izx+1][wid][ip][ilr]->Draw();
							}
						}
					}
					// position of the track point
					double z_track = 0;
					double x_track = 0;
					double z_itrack = 0;
					double x_itrack = 0;
					// is there a cross point?
					if (nHits_layer[izx]>0&&nHits_layer[izx+1]>0){
						for (int wid = 0; wid<NCEL; wid++){
							for (int ip = 0; ip<nHits_cell[izx][wid]; ip++){
								int ihit = iHit_cell[izx][wid][ip];
								for (int wjd = 0; wjd<NCEL; wjd++){
									for (int jp = 0; jp<nHits_cell[izx+1][wjd]; jp++){
										int jhit = iHit_cell[izx+1][wjd][jp];
										int itype = (*i_type)[ihit];
										int jtype = (*i_type)[jhit];
										int isel = 0;
										int jsel = 0;
										if (m_workMode%10>0){
											isel = (*i_sel[iCand])[ihit];
											jsel = (*i_sel[iCand])[jhit];
										}
										// position of the cross point
										for (int icombi = 0; icombi<4; icombi++){
											if (itype>3||jtype>3){ // bad cross
												point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(14);
											}
											else if (isel&&jsel){ // selected cross
												point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(kBlack);
											}
											else{ // good cross
												point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(kGreen);
											}
											double dd1 = dd_cell[izx][wid][ip];
											double dd2 = dd_cell[izx+1][wjd][jp];
											if (icombi<2) dd1 = -dd1; // reverse izx when icombi is 0 or 1
											if (icombi%2==0) dd2 = -dd2; // reverse izx+1 when icombi is 0 or 2
											double theta1 = map_theta[izx][wid];
											double theta2 = map_theta[izx+1][wjd];
											double sintheta12 = sin(theta1-theta2);
											double zc_fix_slx = 0;
											if (m_workMode%10>0){
												double deltaY = y_cell[izx+1][wjd]-y_cell[izx][wid];
												zc_fix_slx = deltaY*i_slx[iCand]/(tan(theta2)-tan(theta1));
											}
											double xc = mcp_xc[izx][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
											double zc = mcp_zc[izx][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
											if (zc>-chamberHL&&zc<chamberHL){
												point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetX(zc);
												point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetY(xc);
												point_cross_zx[izx][wid][wjd][ip][jp][icombi]->Draw();
											}
											if (m_workMode%10>0){
												if (icombi==3){
													double y_track = (y_cell[izx][wid]+y_cell[izx+1][wjd])/2.;
													z_track = i_inz[iCand]+(y_track-iny)*i_slz[iCand];
													x_track = i_inx[iCand]+(y_track-iny)*i_slx[iCand];
													z_itrack = i_iinz[iCand]+(y_track-iny)*i_islz[iCand];
													x_itrack = i_iinx[iCand]+(y_track-iny)*i_islx[iCand];
													double fd1 = fd_cell[izx][wid][ip];
													double fd2 = fd_cell[izx+1][wjd][jp];
													double xcf = mcp_xc[izx][wid][wjd]+fd1*sin(theta2)/(-sintheta12)+fd2*sin(theta1)/sintheta12;
													double zcf = mcp_zc[izx][wid][wjd]+fd1*cos(theta2)/(-sintheta12)+fd2*cos(theta1)/sintheta12+zc_fix_slx;
													gr_all[izx]->SetTitle(Form("DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) #Delta_{x}: %.0lf #mum #Delta_{z}: %.0lf #mum",izx,wid,dd1,(fd1-dd1)*1000,izx+1,wjd,dd2,(fd2-dd2)*1000,(xc-x_track)*1000,(zc-z_track)*1000));
												}
											}
											else{
												gr_all[izx]->SetTitle(Form("Layer %d and Layer %d",izx,izx+1));
											}
										}
									}
								}
							}
						}
					}
					else{
						int wid = 0;
						for (; wid<NCEL; wid++){
							if (nHits_cell[izx][wid]) break;
						}
						if (wid==NCEL) wid = NCEL/2;
						int wjd = 0;
						for (; wjd<NCEL; wjd++){
							if (nHits_cell[izx+1][wjd]) break;
						}
						if (wjd==NCEL) wjd = NCEL/2;
						double y_track = (y_cell[izx][wid]+y_cell[izx+1][wjd])/2.; // take the y value from a previous event
						if (m_workMode%10>0){
							z_track = i_inz[iCand]+(y_track-iny)*i_slz[iCand];
							x_track = i_inx[iCand]+(y_track-iny)*i_slx[iCand];
							z_itrack = i_iinz[iCand]+(y_track-iny)*i_islz[iCand];
							x_itrack = i_iinx[iCand]+(y_track-iny)*i_islx[iCand];
						}
						gr_all[izx]->SetTitle(Form("Layer %d and Layer %d",izx,izx+1));
					}
					// draw the track point
					if (m_workMode%10>0){
						point_itrack_zx[izx]->SetX(z_itrack);
						point_itrack_zx[izx]->SetY(x_itrack);
						point_itrack_zx[izx]->Draw();
						point_track_zx[izx]->SetX(z_track);
						point_track_zx[izx]->SetY(x_track);
						point_track_zx[izx]->Draw();
					}
					if (m_workMode%10>0){
						ca_zx[izx]->SaveAs(prefix+Form("zx.%d.i%d.l%d.pdf",iEntry,iCand,izx));
						ca_zx[izx]->SaveAs(prefix+Form("zx.%d.i%d.l%d.png",iEntry,iCand,izx));
					}
					else{
						ca_zx[izx]->SaveAs(prefix+Form("zx.%d.l%d.pdf",iEntry,izx));
						ca_zx[izx]->SaveAs(prefix+Form("zx.%d.l%d.png",iEntry,izx));
					}
				}

				// create one more z-x plane with all points on it
				ca_zx_all->cd();
				gr_all[0]->Draw("AP");
				for (int izx = 1; izx<NZXP; izx++){
					if (m_workMode%10>0)
						point_track_zx[izx]->Draw();
					if (nHits_layer[izx]>0&&nHits_layer[izx+1]>0){
						for (int wid = 0; wid<NCEL; wid++){
							for (int ip = 0; ip<nHits_cell[izx][wid]; ip++){
								int ihit = iHit_cell[izx][wid][ip];
								for (int wjd = 0; wjd<NCEL; wjd++){
									for (int jp = 0; jp<nHits_cell[izx+1][wjd]; jp++){
										int jhit = iHit_cell[izx+1][wjd][jp];
										int itype = (*i_type)[ihit];
										int jtype = (*i_type)[jhit];
										// position of the cross point
										for (int icombi = 0; icombi<4; icombi++){
											if (itype>3||jtype>3){ // bad cross
												text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextColor(14);
											}
											double dd1 = dd_cell[izx][wid][ip];
											double dd2 = dd_cell[izx+1][wjd][jp];
											if (icombi<2) dd1 = -dd1; // reverse izx when icombi is 0 or 1
											if (icombi%2==0) dd2 = -dd2; // reverse izx+1 when icombi is 0 or 2
											double theta1 = map_theta[izx][wid];
											double theta2 = map_theta[izx+1][wjd];
											double sintheta12 = sin(theta1-theta2);
											double zc_fix_slx = 0;
											if (m_workMode%10>0){
												double deltaY = y_cell[izx+1][wjd]-y_cell[izx][wid];
												zc_fix_slx = deltaY*i_slx[iCand]/(tan(theta2)-tan(theta1));
											}
											double xc = mcp_xc[izx][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
											double zc = mcp_zc[izx][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
											if (zc>-chamberHL&&zc<chamberHL){
												text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetX(zc);
												text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetY(xc);
												text_cross_zx[izx][wid][wjd][ip][jp][icombi]->Draw();
											}
										}
									}
								}
							}
						}
					}
				}
				if (m_workMode%10>0){
					ca_zx_all->SaveAs(prefix+Form("zx.%d.i%d.all.png",iEntry,iCand));
					ca_zx_all->SaveAs(prefix+Form("zx.%d.i%d.all.pdf",iEntry,iCand));
				}
				else{
					ca_zx_all->SaveAs(prefix+Form("zx.%d.all.png",iEntry));
					ca_zx_all->SaveAs(prefix+Form("zx.%d.all.pdf",iEntry));
				}
			}
        }

		// FIXME: in interactive mode
//		ca_xyADC->WaitPrimitive();
	//	ca_xyADC->Update();
//		while(1){}
	}
	return 0;
}

double t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
    TF1* f=0;
    // FIXME: now we only take the default xt: layer 0 wire 0 (fake) 
    //int index = (lid-1)*NCEL+wid;
    int index = 0;
    if (lr>=0){
        f = f_right[index];
    }
    else {
        f = f_left[index];
    }
    if (!f){
        fprintf(stderr,"Cannot get f[%d]!\n",index);
        status = -2;
        return 0;
    }
    double tmax = f->GetXmax();
    double tmin = f->GetXmin();
    status = 0;
    double dd = 0;
    if (time>tmax){
        status = 1;
        dd = f->Eval(tmax);
    }
    else if (time<tmin){
        status = -1;
        dd = 0;
    }
    else {
        dd = f->Eval(time);
    }
    return dd;
}

double tdc2t(int deltaTDC){
    return (deltaTDC)/0.96;
}

void cid4dr2lidwid(int cid, int & lid, int & wid){
	lid = cid/NWIRE4DR+IFLAY4DR; // skip the first layer (dummy) by counting from 1.
	wid = cid%NWIRE4DR+IFWIRE4DR; // skip first wire (on boundary)  by counting from 1.

	// for test 180404
//	lid = 5;
//	if (cid==0) wid=6;
//	else if (cid==1) wid=7;
//	else if (cid==2) wid=1;
//	else if (cid==3) wid=2;
}

bool lidwid2cid4dr(int lid, int wid, int & cid){
    cid = (lid-IFLAY4DR)*NWIRE4DR+wid-IFWIRE4DR;
	if (lid<IFLAY4DR||lid-IFLAY4DR>=NLAY4DR||wid<IFWIRE4DR||wid-IFWIRE4DR>=NWIRE4DR) return true;
	else return false;

	// for test 180404
//	if (lid!=5) return true;
//	if (wid==6) cid=0;
//	else if (wid==7) cid=1;
//	else if (wid==1) cid=2;
//	else if (wid==2) cid=3;
//	else return true;
//	return false;
}

void print_usage(char* prog_name){
	fprintf(stderr,"Usage %s [options] runname\n",prog_name);
	fprintf(stderr,"[options]\n");
	fprintf(stderr,"\t -D <name>=[error,severe,warn,debug,trace]\n");
	fprintf(stderr,"\t\t Change the named debug level\n");
	fprintf(stderr,"\t -V <name>=[quiet,log,info,verbose]\n");
	fprintf(stderr,"\t\t Change the named log level\n");
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
    fprintf(stderr,"\t -l <l>\n");
    fprintf(stderr,"\t\t The layer set to l\n");
    fprintf(stderr,"\t -W <w>\n");
    fprintf(stderr,"\t\t The wire set to w\n");
    fprintf(stderr,"\t -G <g>\n");
    fprintf(stderr,"\t\t The geometry setup set to g\n");
    fprintf(stderr,"\t -X <xtfile>\n");
    fprintf(stderr,"\t\t Instead of searching for the xt file from this run, use the provided one\n");
    fprintf(stderr,"\t -m <mode>\n");
    fprintf(stderr,"\t\t Set the workmode (%d)\n",m_workMode);
    fprintf(stderr,"\t\t 0: h_XXX; 1: t_XXX; 2: ana_XXX; 10: h_XXX with zx; 11: t_XXX with zx\n");
	return;
}
