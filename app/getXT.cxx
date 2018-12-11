#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <stdio.h>  /* printf, fgets */
#include <unistd.h> /* getopt */
#include <stdlib.h> /* atoi, atof */
#include <iostream> /* cout */

#include "TString.h"
#include "TChain.h"
#include "TFile.h"

#include "XTAnalyzer.hxx"
#include "MyProcessManager.hxx"
#include "MyRuntimeParameters.hxx"
#include "Log.hxx"
#include "header.hxx"

void print_usage(char * prog_name);
int getHitType(int type,bool isRight);

MyProcessManager * pMyProcessManager;

int main(int argc, char** argv){
	int m_runNo = 0;
    TString m_prerunname  = "prerun";
    TString m_runname = "currun";
    int m_iEntryStart = 0;
    int m_iEntryStop = 0;
    int m_debugLevel = 0;
    int m_modulo = 10000;
    bool m_memdebug = false;
    int m_saveHists = 0;
    int m_defaultLayerID = 4;
    TString m_configureFile = "";
    int m_nHitsMax = 30;
    double m_maxchi2 = 2;
	int m_geoSetup = 0; // 0: normal scintillator; 1: finger scintillator
    int m_inputType = 0; // 1 for MC; 0 for data
	int m_peakType = 2; // 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
	int m_xtType = 2;

    double temp_maxchi2 = 0; bool set_maxchi2 = false;
    int temp_nHitsMax = 0; bool set_nHitsMax = false;
	int temp_geoSetup = 0; bool set_geoSetup = false;
    int temp_inputType = 0; bool set_inputType = false;
	int temp_peakType = 0; bool set_peakType = false;
	int temp_xtType = 0; bool set_xtType = false;

    std::map<std::string, Log::ErrorPriority> namedDebugLevel;
    std::map<std::string, Log::LogPriority> namedLogLevel;
    int    opt_result;
	while((opt_result=getopt(argc,argv,"M:R:B:E:L:C:n:c:g:i:p:x:D:V:"))!=-1){
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
			    m_defaultLayerID = atoi(optarg);
                printf("Test layer set to %d\n",m_defaultLayerID);
			case 'C':
				m_configureFile = optarg;
                printf("Using configure file \"%s\"\n",optarg);
                break;
			case 'n':
			    temp_nHitsMax = atoi(optarg);set_nHitsMax = true;
                printf("Maximum number of hits cut set to %d\n",temp_nHitsMax);
			case 'c':
			    temp_maxchi2 = atof(optarg);set_maxchi2 = true;
                printf("Maximum chi2 cut set to %d\n",temp_maxchi2);
			case 'g':
			    temp_geoSetup = atoi(optarg);set_geoSetup = true;
                printf("Geometry setup set to %d\n",temp_geoSetup);
			case 'i':
			    temp_inputType = atoi(optarg);set_inputType = true;
                printf("Input type set to %d\n",temp_inputType);
			case 'p':
			    temp_peakType = atoi(optarg);set_peakType = true;
                printf("Peak type set to %d\n",temp_peakType);
			case 'x':
			    temp_xtType = atoi(optarg);set_xtType = true;
                printf("XT type set to %d\n",temp_xtType);
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
        m_nHitsMax = MyRuntimeParameters::Get().GetParameterI("getXT.nHitsMax");
        m_maxchi2 = MyRuntimeParameters::Get().GetParameterD("getXT.maxchi2");
        m_geoSetup = MyRuntimeParameters::Get().GetParameterI("geoSetup");
        m_inputType = MyRuntimeParameters::Get().GetParameterI("inputType");
        m_peakType = MyRuntimeParameters::Get().GetParameterI("peakType");
        m_xtType = MyRuntimeParameters::Get().GetParameterI("getXT.xtType");;
    }
    if (set_nHitsMax) m_nHitsMax = temp_nHitsMax;
    if (set_maxchi2) m_maxchi2 = temp_maxchi2;
    if (set_geoSetup) m_geoSetup = temp_geoSetup;
    if (set_inputType) m_inputType = temp_inputType;
    if (set_peakType) m_peakType = temp_peakType;
    if (set_xtType) m_xtType = temp_xtType;

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
    printf("default layer: %d\n",m_defaultLayerID);
    printf("maxchi2     = %.3e\n",m_maxchi2);
    printf("maxNhits    = %d\n",m_nHitsMax);
    printf("geoSetup:     %s\n",m_geoSetup==0?"normal scintillator":"finger scintillator");
    printf("xtType:       %d\n",m_xtType);
    printf("inputType   = %d, %s\n",m_inputType,m_inputType==0?"Real Data":"MC");
    printf("debug       = %d\n",m_debugLevel);
    printf("print modulo= %d\n",m_modulo);
    printf("save slice fittings? \"%s\"\n",m_saveHists?"yes":"no");
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
	printf("runNo#%d: %s, %d, %s, %d V, %d mV, %.0f sec\n",m_runNo,gastype.Data(),runGr,duration.Data(),HV,THR,durationTime);

	// get offset
	double off[NLAY][NCEL];
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			off[lid][wid] = 0;
		}
	}
	if (m_xtType==1||m_xtType==6||m_xtType==7){
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
				off[i_off_lid][i_off_wid] = i_off_delta;
		}
	}

	// get wire map
	int map_ch[NLAY][NCEL];
	int map_bid[NLAY][NCEL];
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
        	map_ch[lid][wid] = -1;
        	map_bid[lid][wid] = -1;
		}
	}
    TFile * TFile_wirepos = new TFile(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),m_runNo,m_runname.Data()));
    TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
    int     wp_bid;
    int     wp_ch;
    int     wp_wid;
    int     wp_lid;
    TTree_wirepos->SetBranchAddress("b",&wp_bid);
    TTree_wirepos->SetBranchAddress("ch",&wp_ch);
    TTree_wirepos->SetBranchAddress("l",&wp_lid);
    TTree_wirepos->SetBranchAddress("w",&wp_wid);
    for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
        TTree_wirepos->GetEntry(i);
        if (wp_lid>=0&&wp_lid<NLAY&&wp_wid>=0&&wp_wid<NCEL){
            map_ch[wp_lid][wp_wid] = wp_ch;
            map_bid[wp_lid][wp_wid] = wp_bid;
        }
    }
    TFile_wirepos->Close();

    // get XT file of the previous run
    TFile * preXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_prerunname.Data()));

    // prepare new XT file for this run
    TFile * newXTFile = new TFile(Form("%s/info/xt.%d.%s.root",HOME.Data(),m_runNo,m_runname.Data()),"RECREATE");
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
    std::vector<double> * i_aa = 0;
    std::vector<double> * i_ped = 0;
    bool has_ped = false;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_driftD[NCAND] = {0};
    std::vector<double> * i_calD[NCAND] = {0};
    int npairs[NCAND];
    int isel[NCAND];
    int icom[NCAND];
    double iinx[NCAND];
    double islx[NCAND];
    double iinz[NCAND];
    double islz[NCAND];
    double chi2x[NCAND];
    double chi2z[NCAND];
    double chi2i[NCAND];
    int nHitsS[NCAND];
    double inx[NCAND];
    double slx[NCAND];
    double inz[NCAND];
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
    std::vector<double> * i_fitD[NCAND] = {0};
    std::vector<int> * i_sel[NCAND] = {0};

	// output file
    std::vector<double> * o_driftD = 0;
    std::vector<int>    * o_driftDs = 0;
    std::vector<int> * o_channelID = 0;
    std::vector<int> * o_boardID = 0;
	int nHitsSmallAll = 0;
	int nHitsSmallSASD = 0;
	int nSmallSumHits = 0;
	int nShadowedHits = 0;
	int nLateHits = 0;
	int nBoundaryHits = 0;
	int nSmallBoundaryHits = 0;
	// the closest peak to the track in the test layer
	double minres = 1e9;
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

	//=================================================Start to get XT====================================================
	// Prepare XTAnalyzer
	XTAnalyzer * fXTAnalyzer = new XTAnalyzer(gasID,m_debugLevel);
    // Loop in layers
	for (int testLayer = 0; testLayer<NLAY; testLayer++){
		//----------------------------------Set input file--------------------------------------------
        if (m_debugLevel>0) {printf("In Layer %d: preparing input TChain\n",testLayer);fflush(stdout);}
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
        ichain->SetBranchAddress("aa",&i_aa);
        has_ped = (ichain->SetBranchAddress("ped",&i_ped)==0);
		ichain->SetBranchAddress("sum",&i_sum);
		for (int iCand = 0; iCand<NCAND; iCand++){
			ichain->SetBranchAddress(Form("driftD%d",iCand),&(i_driftD[iCand]));
			ichain->SetBranchAddress(Form("npairs%d",iCand),&(npairs[iCand]));
			ichain->SetBranchAddress(Form("isel%d",iCand),&(isel[iCand]));
			ichain->SetBranchAddress(Form("icom%d",iCand),&(icom[iCand]));
			ichain->SetBranchAddress(Form("islx%d",iCand),&(islx[iCand]));
			ichain->SetBranchAddress(Form("islz%d",iCand),&(islz[iCand]));
			ichain->SetBranchAddress(Form("iinx%d",iCand),&(iinx[iCand]));
			ichain->SetBranchAddress(Form("iinz%d",iCand),&(iinz[iCand]));
			ichain->SetBranchAddress(Form("chi2x%d",iCand),&(chi2x[iCand]));
			ichain->SetBranchAddress(Form("chi2z%d",iCand),&(chi2z[iCand]));
			ichain->SetBranchAddress(Form("chi2i%d",iCand),&(chi2i[iCand]));
			ichain->SetBranchAddress(Form("calD%d",iCand),&(i_calD[iCand]));
			ichain->SetBranchAddress(Form("nHitsS%d",iCand),&(nHitsS[iCand]));
			ichain->SetBranchAddress(Form("slx%d",iCand),&(slx[iCand]));
			ichain->SetBranchAddress(Form("slz%d",iCand),&(slz[iCand]));
			ichain->SetBranchAddress(Form("inx%d",iCand),&(inx[iCand]));
			ichain->SetBranchAddress(Form("inz%d",iCand),&(inz[iCand]));
			ichain->SetBranchAddress(Form("chi2%d",iCand),&(chi2[iCand]));
			ichain->SetBranchAddress(Form("chi2p%d",iCand),&(chi2p[iCand]));
			ichain->SetBranchAddress(Form("chi2a%d",iCand),&(chi2a[iCand]));
			if (m_inputType){
                ichain->SetBranchAddress(Form("chi2mc%d",iCand),&(chi2mc[iCand]));
                ichain->SetBranchAddress(Form("chi2pmc%d",iCand),&(chi2pmc[iCand]));
                ichain->SetBranchAddress(Form("chi2amc%d",iCand),&(chi2amc[iCand]));
            }
			ichain->SetBranchAddress(Form("fitD%d",iCand),&(i_fitD[iCand]));
			ichain->SetBranchAddress(Form("sel%d",iCand),&(i_sel[iCand]));
		}
		if (m_inputType){
			ichain->SetBranchAddress("slxmc",&slxmc);
			ichain->SetBranchAddress("slzmc",&slzmc);
			ichain->SetBranchAddress("inxmc",&inxmc);
			ichain->SetBranchAddress("inzmc",&inzmc);
		}

		//----------------------------------Initialize the analyzer--------------------------------------------
		int saveEvenOdd = 0; if (testLayer==4) saveEvenOdd = 1; else if (testLayer==5) saveEvenOdd = -1;
		int statusInitialize = fXTAnalyzer->Initialize(Form("%d.%s.layer%d",m_runNo,m_runname.Data(),testLayer),testLayer,preXTFile,newXTFile,newXTTree,m_xtType,m_saveHists, testLayer==m_defaultLayerID, saveEvenOdd, testLayer!=0);
		if (statusInitialize){
			fprintf(stderr,"WARNING: something wrong with initializing XTAnalyzer for layer[%d], will ignore this layer!\n",testLayer);
			continue;
		}

		//----------------------------------Loop in events--------------------------------------------
        if (!m_iEntryStart&&!m_iEntryStop){
        	m_iEntryStart = 0;
        	m_iEntryStop = N-1;
        }
        if (m_debugLevel>0) {printf("Processing %d events\n",N);fflush(stdout);}
        for ( int iEntry = m_iEntryStart ; iEntry<=m_iEntryStop; iEntry++){
            if (iEntry%m_modulo==0) printf("%d\n",iEntry);
            if (m_debugLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

			// decide which candidate to use
			int theCand = 0;
			if (m_xtType==3){
				int nLateHitsMin = 1e9;
				for (int iCand = 0; iCand<NCAND; iCand++){
					nLateHits = 0;
					for (int ihit = 0; ihit<nHits; ihit++){
						int ip = 0;
						for (int jhit = ihit-1; jhit>0; jhit--){
							if ((*i_layerID)[jhit]!=(*i_layerID)[ihit]) break;
							int type = getHitType((*i_type)[jhit],(*i_fitD[iCand])[jhit]>=0);
							if (type<100) ip++;
						}
						if ((*i_sel[iCand])[ihit]==1){
							if(ip!=0)
								nLateHits++;
						}
					}
					if (nLateHits<nLateHitsMin){
						nLateHitsMin = nLateHits;
						theCand = iCand;
					}
				}
			}
            else if (m_xtType==4||m_xtType==5){
                double minchi2 = 1e9;
                int minNhitsS = 0;
				for (int iCand = 0; iCand<NCAND; iCand++){
                    if (m_xtType==4){
                        if ((minchi2>chi2a[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                            theCand = iCand;
                            minchi2 = chi2a[iCand];
                            minNhitsS = nHitsS[iCand];
                        }
                    }
                    else if (m_xtType==5){
                        if ((minchi2>chi2[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                            theCand = iCand;
                            minchi2 = chi2[iCand];
                            minNhitsS = nHitsS[iCand];
                        }
                    }
                }
            }

            // ignore events with bad fitting
            if (nHitsS[theCand]<7) continue;
			if (chi2[theCand]>m_maxchi2) continue;
            //if (nHitsG>nHitsS[theCand]) continue;
            if (m_geoSetup==1){
                if (fabs(inz[theCand])>24) continue;
            }
            else{
                if (fabs(slz[theCand])>0.11) continue;
            }
            if (m_nHitsMax&&nHits>m_nHitsMax) continue;

            if (m_debugLevel>=20) printf("  Good Event! Looping in %d hits\n",nHits);
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
                double tfitD = (*i_fitD[theCand])[ihit]-off[tlayerID][twireID];
                double tdriftD = (*i_driftD[theCand])[ihit];
            	if ((*i_sel[theCand])[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                if (tlayerID!=testLayer) continue;
				int ttype = getHitType((*i_type)[ihit],tfitD>=0);
                if ((ttype<100&&(m_xtType==7||m_xtType==6)||(m_xtType!=6&&m_xtType!=7))&&fabs(tfitD-tdriftD)<fabs(minres)){ // Should have cut for test layer! otherwise XT will not be well tuned
                //if (fabs(tfitD-tdriftD)<fabs(minres)){ // no cut for test layer!
                    minres = tfitD-tdriftD;
                    wireID = (*i_wireID)[ihit];
                    fitD = tfitD;
                    driftT = (*i_driftT)[ihit];
                    has = true;
                }
            }
            if (!has) continue; // no hits found in test layer
            //if (hasBadHit) continue;

            if (m_debugLevel>=20) printf("  Found hit! pushing to XTAnalyzer\n");
			// tell analyzer a new data point
            fXTAnalyzer->Push(driftT,fitD);
        }
        if (m_debugLevel>0) printf("Starting XT analysis\n");
        // fit histograms/graphs, make plots, and save new xt file
        fXTAnalyzer->Process();

		//----------------------------------prepare for output ROOT file--------------------------------------------
        TFile * ofile = new TFile(Form("%s/root/ana_%d.%s.layer%d.root",HOME.Data(),m_runNo,m_runname.Data(),testLayer),"RECREATE");
        TTree * otree = new TTree("t","t");
        otree->Branch("triggerNumber",&triggerNumber);
        otree->Branch("res",&minres);
        otree->Branch("theDD",&theDD);
        otree->Branch("theDT",&theDT);
        otree->Branch("theWid",&theWid);
        otree->Branch("theSum",&theSum);
        otree->Branch("sum1st",&sum1st);
        otree->Branch("dt1st",&dt1st);
        otree->Branch("has",&has);
        otree->Branch("thePeak",&thePeak);
        otree->Branch("theHeight",&theHeight);
        otree->Branch("theIp",&theIp);
        otree->Branch("theMpi",&theMpi);
        otree->Branch("theCand",&theCand);
        otree->Branch("highBid",&highBid);
        otree->Branch("highCh",&highCh);
        otree->Branch("highLid",&highLid);
        otree->Branch("highWid",&highWid);
        otree->Branch("highIp",&highIp);
        otree->Branch("highSum",&highSum);
        otree->Branch("highAA",&highAA);
        otree->Branch("highDT",&highDT);
        otree->Branch("nHitsSmallSASD",&nHitsSmallSASD);
        otree->Branch("nHitsSmallAll",&nHitsSmallAll);
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
        otree->Branch("aa",&i_aa);
        if (has_ped) otree->Branch("ped",&i_ped);
        otree->Branch("sum",&i_sum);
        otree->Branch("driftD",&o_driftD);
        otree->Branch("driftDs",&o_driftDs);
        otree->Branch("driftD0",&(i_driftD[0]));
        otree->Branch("npairs",&(npairs[0]));
        otree->Branch("isel",&(isel[0]));
        otree->Branch("icom",&(icom[0]));
        otree->Branch("islx",&(islx[0]));
        otree->Branch("islz",&(islz[0]));
        otree->Branch("iinx",&(iinx[0]));
        otree->Branch("iinz",&(iinz[0]));
        otree->Branch("chi2x",&(chi2x[0]));
        otree->Branch("chi2z",&(chi2z[0]));
        otree->Branch("chi2i",&(chi2i[0]));
        otree->Branch("calD",&(i_calD[0]));
        otree->Branch("nHitsS",&(nHitsS[0]));
        otree->Branch("slx",&(slx[0]));
        otree->Branch("slz",&(slz[0]));
        otree->Branch("inx",&(inx[0]));
        otree->Branch("inz",&(inz[0]));
        otree->Branch("chi2",&(chi2[0]));
        otree->Branch("chi2p",&(chi2p[0]));
        otree->Branch("chi2a",&(chi2a[0]));
        if (m_inputType){
            otree->Branch("chi2mc",&(chi2mc[0]));
            otree->Branch("chi2pmc",&(chi2pmc[0]));
            otree->Branch("chi2amc",&(chi2amc[0]));
			otree->Branch("slxmc",&slxmc);
			otree->Branch("slzmc",&slzmc);
			otree->Branch("inxmc",&inxmc);
			otree->Branch("inzmc",&inzmc);
        }
        otree->Branch("fitD",&(i_fitD[0]));
        otree->Branch("sel",&(i_sel[0]));
        o_driftD = new std::vector<double>;
        o_driftDs = new std::vector<int>;
        o_channelID = new std::vector<int>;
        o_boardID = new std::vector<int>;

		//----------------------------------Loop in events--------------------------------------------
        for ( int iEntry = m_iEntryStart ; iEntry<=m_iEntryStop; iEntry++){
            if (iEntry%10000==0) printf("%d\n",iEntry);
            if (m_debugLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

			// decide which candidate to use
			theCand = 0;
			if (m_xtType==3){
				int nLateHitsMin = 1e9;
				for (int iCand = 0; iCand<NCAND; iCand++){
					nLateHits = 0;
					for (int ihit = 0; ihit<nHits; ihit++){
						int ip = 0;
						for (int jhit = ihit-1; jhit>0; jhit--){
							if ((*i_layerID)[jhit]!=(*i_layerID)[ihit]) break;
							int type = getHitType((*i_type)[jhit],(*i_fitD[iCand])[jhit]>=0);
							if (type<100) ip++;
						}
						if ((*i_sel[iCand])[ihit]==1){
							if(ip!=0)
								nLateHits++;
						}
					}
					if (nLateHits<nLateHitsMin){
						nLateHitsMin = nLateHits;
						theCand = iCand;
					}
				}
                otree->SetBranchAddress("driftD",&(i_driftD[theCand]));
                otree->SetBranchAddress("npairs",&(npairs[theCand]));
                otree->SetBranchAddress("isel",&(isel[theCand]));
                otree->SetBranchAddress("icom",&(icom[theCand]));
                otree->SetBranchAddress("islx",&(islx[theCand]));
                otree->SetBranchAddress("islz",&(islz[theCand]));
                otree->SetBranchAddress("iinx",&(iinx[theCand]));
                otree->SetBranchAddress("iinz",&(iinz[theCand]));
                otree->SetBranchAddress("chi2x",&(chi2x[theCand]));
                otree->SetBranchAddress("chi2z",&(chi2z[theCand]));
                otree->SetBranchAddress("chi2i",&(chi2i[theCand]));
                otree->SetBranchAddress("calD",&(i_calD[theCand]));
                otree->SetBranchAddress("nHitsS",&(nHitsS[theCand]));
                otree->SetBranchAddress("slx",&(slx[theCand]));
                otree->SetBranchAddress("slz",&(slz[theCand]));
                otree->SetBranchAddress("inx",&(inx[theCand]));
                otree->SetBranchAddress("inz",&(inz[theCand]));
                otree->SetBranchAddress("chi2",&(chi2[theCand]));
                otree->SetBranchAddress("chi2p",&(chi2p[theCand]));
                otree->SetBranchAddress("chi2a",&(chi2a[theCand]));
                if (m_inputType){
                    otree->SetBranchAddress("chi2mc",&(chi2mc[theCand]));
                    otree->SetBranchAddress("chi2pmc",&(chi2pmc[theCand]));
                    otree->SetBranchAddress("chi2amc",&(chi2amc[theCand]));
                    otree->SetBranchAddress("slxmc",&slxmc);
                    otree->SetBranchAddress("slzmc",&slzmc);
                    otree->SetBranchAddress("inxmc",&inxmc);
                    otree->SetBranchAddress("inzmc",&inzmc);
                }
                otree->SetBranchAddress("fitD",&(i_fitD[theCand]));
                otree->SetBranchAddress("sel",&(i_sel[theCand]));
			}
            else if (m_xtType==4||m_xtType==5){
                double minchi2 = 1e9;
                int minNhitsS = 0;
				for (int iCand = 0; iCand<NCAND; iCand++){
                    if (m_xtType==4){
                        if ((minchi2>chi2a[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                            theCand = iCand;
                            minchi2 = chi2a[iCand];
                            minNhitsS = nHitsS[iCand];
                        }
                    }
                    else if (m_xtType==5){
                        if ((minchi2>chi2[iCand]&&minNhitsS==nHitsS[iCand])||minNhitsS<nHitsS[iCand]){
                            theCand = iCand;
                            minchi2 = chi2[iCand];
                            minNhitsS = nHitsS[iCand];
                        }
                    }
                }
            }

            // set driftD and extra info
            o_driftD->clear();
            o_driftDs->clear();
            o_channelID->clear();
            o_boardID->clear();
            minres = 1e9;
            theDD = 1e9;
            theDT = 1e9;
            has = 0;
            theWid = -1;
			theSum = 0;
			sum1st = 0;
			dt1st = 1e9;
			thePeak = 0;
			theHeight = 0;
			theIp = 0;
			theMpi = 0;
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
            	// get new ip // FIXME: now taking ip as its original meaning. Can redefine it by ignoring small peaks
				int ip = 0;
				for (int jhit = ihit-1; jhit>=0; jhit--){
					if ((*i_layerID)[jhit]!=(*i_layerID)[ihit]||(*i_wireID)[jhit]!=(*i_wireID)[ihit]) break;
					//int type = getHitType((*i_type)[jhit],(*i_fitD[theCand])[jhit]>=0);
					//if (type<100) ip++;
					ip++;
				}

				// get nXXXHits according to the original hit distance
            	double dd0 = (*i_driftD[theCand])[ihit];
            	if ((*i_sel[theCand])[ihit]==1){
            		if((fabs(dd0)<0.5||fabs(dd0)>7.5))
						nBoundaryHits++;
            		if((fabs(dd0)<0.25||fabs(dd0)>7.75))
						nSmallBoundaryHits++;
					if(ip!=0){
						nLateHits++;
						for (int jhit = ihit-1; jhit>0; jhit--){
							if ((*i_ip)[jhit]==0){
								if (dt1st>(*i_driftT)[jhit]){
									sum1st = (*i_sum)[jhit];
									dt1st = (*i_driftT)[jhit];
								}
								break;
							}
						}
					}
					if((*i_mpi)[ihit]!=0)
						nShadowedHits++;
					if(has_rank&&(*i_rank)[ihit]!=0)
						nSmallSumHits++;
				}

				// get other numbers updated
                int lid = (*i_layerID)[ihit];
                int wid = (*i_wireID)[ihit];
            	int ch = map_ch[lid][wid];
            	int bid = map_bid[lid][wid];
                double aa = (*i_aa)[ihit];
                double sum = (*i_sum)[ihit];
                double fitD = (*i_fitD[theCand])[ihit]-off[lid][wid];
				int type = getHitType((*i_type)[ihit],fitD>=0);
            	double dt = (*i_driftT)[ihit];
            	double dd;
            	int status = fXTAnalyzer->t2d(dt,dd,dd0>0);
            	o_driftD->push_back(dd);
            	o_driftDs->push_back(status);
            	o_channelID->push_back(ch);
            	o_boardID->push_back(bid);
				(*i_fitD[theCand])[ihit]=fitD;

				// check the closest peak
                if (lid==testLayer&&((type<100&&m_xtType==6)||m_xtType!=6)&&status==0&&fabs(fitD-dd)<fabs(minres)){ // Should have cut for test layer! otherwise XT will not be well tuned
                    minres = fitD-dd;
                    theDD = dd;
                    theDT = dt;
                    theWid = wid;
                    theSum = sum;
                    double ped = 220;
                    if (has_ped) ped = (*i_ped)[ihit];
                    thePeak = (*i_peak)[ihit]-ped;
                    theHeight = (*i_height)[ihit]-ped;
                    theIp = ip;
                    theMpi = (*i_mpi)[ihit];
                    has = 1;
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
			}

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

			otree->Fill();
		}
		otree->Write();
        ofile->Close();

        if (m_debugLevel>=20) printf("Finished!\n");
	}

    return 0;
}

int getHitType(int type,bool isRight){
	int ttype = (type/10)%10;
	if (isRight){
		if (ttype==1||ttype==4) type-=ttype*10; // l- or l+
	}
	else{
		if (ttype==2||ttype==5) type-=ttype*10; // r- or r+
	}
	return type;
}

void print_usage(char * prog_name){
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
    fprintf(stderr,"\t\t Default layer set to l\n");
    fprintf(stderr,"\t -n <n>\n");
    fprintf(stderr,"\t\t Maximum number of hits cut set to n\n");
    fprintf(stderr,"\t -c <c>\n");
    fprintf(stderr,"\t\t Maximum chi2 cut set to c\n");
    fprintf(stderr,"\t -g <g>\n");
    fprintf(stderr,"\t\t Geometry setup set to g\n");
    fprintf(stderr,"\t\t (0): normal; 1: finger\n");
    fprintf(stderr,"\t -i <i>\n");
    fprintf(stderr,"\t\t Input type set to i\n");
    fprintf(stderr,"\t\t (0) for data; 1 for MC\n");
    fprintf(stderr,"\t -p <p>\n");
    fprintf(stderr,"\t\t Peak type set to p\n");
    fprintf(stderr,"\t\t (0) only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks\n");
    fprintf(stderr,"\t -x <x>\n");
    fprintf(stderr,"\t\t xt type set to x\n");
    fprintf(stderr,"\t\t (2) sym; 1 sym+offset; 0 no; 6 sym+offset+first OT peak; 7 sym+offset+first OT peak+2segments\n");
}
