#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

#include "MyProcessManager.h"

#define NZXP 8
#define NLAY 9
#define NLAYG 8
#define NCEL 11

#define NBRD 2
#define NCHS 48

#define NCAND 4

int workType = 0;
int debug = 0;
int memdebug = 0;

//===================Chamber Parameter============================
double U = 8; // mm
double chamberHL = 599.17/2; // mm
double chamberHH = 170.05/2; // mm
double chamberCY = 572; // mm
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
int geoSetup = 0; // 0: normal; 1: finger
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

//===================About tracking============================
int testlayer = 0;
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

    if (argc<5){
        print_usage(argv[0]);
        return 1;
    }
    int runNo = (int)strtol(argv[1],NULL,10);
    testlayer = (int)strtol(argv[2],NULL,10);
    TString prerunname = argv[3];
    TString runname = argv[4];
    int nHitsMax = 10;
    if (argc>=6){
        nHitsMax = (int)strtol(argv[5],NULL,10);
    }
    int t0shift = 0;
    if (argc>=7) t0shift = (int)strtol(argv[6],NULL,10);
    int tmin = -10;
    if (argc>=8){
        tmin = (int)strtol(argv[7],NULL,10);
    }
    int tmax = 800;
    if (argc>=9){
        tmax = (int)strtol(argv[8],NULL,10);
    }
    if (argc>=10){
        geoSetup = (int)strtol(argv[9],NULL,10);
    }
    float sumCut = 0;
    if (argc>=11){
        sumCut = (float)atof(argv[10]);
    }
    float aaCut = 0;
    if (argc>=12){
        aaCut = (float)atof(argv[11]);
    }
    int iEntryStart = 0;
    int iEntryStop = 0;
    if (argc>=14){
        iEntryStart = (int)strtol(argv[12],NULL,10);
        iEntryStop = (int)strtol(argv[13],NULL,10);
    }
    if (argc>=15){
        workType = (int)strtol(argv[14],NULL,10);
    }
    if (argc>=16){
        debug = (int)strtol(argv[15],NULL,10);
    }
    if (argc>=17){
        memdebug = (int)strtol(argv[16],NULL,10);
    }
    printf("##############Input Parameters##################\n");
    printf("runNo       = %d\n",runNo);
    printf("test layer  = %d\n",testlayer);
    printf("prerunname  = \"%s\"\n",prerunname.Data());
    printf("runname     = \"%s\"\n",runname.Data());
    printf("nHitsMax    = %d\n",nHitsMax);
    printf("t0shift     = %d\n",t0shift);
    printf("tmin        = %d\n",tmin);
    printf("tmax        = %d\n",tmax);
    printf("geoSetup:     %s\n",geoSetup==0?"normal scintillator":"finger scintillator");
    printf("sumCut      = %f\n",sumCut);
    printf("aaCut       = %f\n",aaCut);
    printf("Start Entry = %d\n",iEntryStart);
    printf("Stop Entry  = %d\n",iEntryStop);
    printf("workType    = %s\n",workType==0?"all as 0":(workType==1?"even/odd":(workType==-1?"even/odd reversed":"all layers")));
    printf("debug       = %d\n",debug);
    printf("memdebug    = %d\n",memdebug);

    TString HOME=getenv("CDCS8WORKING_DIR");

    if (memdebug>0){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
        std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
    }

	//===================Set scintillator geometry============================
	if (geoSetup==0){
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
            if (lid <NZXP){ // z-x planes corresponding to the layerID of the lower layer counting from 1 
                for (int wjd = 0; wjd<NCEL; wjd++){
                    mcp_xc[lid][wid][wjd] = 999;
                    mcp_zc[lid][wid][wjd] = 999;
                }
            }
        }
    }

    //===================Get Wire Position============================
    TFile * TFile_wirepos = new TFile(HOME+"/info/wire-position.root");
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
            if(debug>0) printf("map_theta[%d][%d] = atan(-(%.3e-%.3e)/%.3e/2) = %.3e\n",wp_lid,wp_wid,wp_xhv,wp_xro,chamberHL,map_theta[wp_lid][wp_wid]);
        }
        else{
            fprintf(stderr,"ERROR: Entry %d in wiremap file, lid = %d wid = %d out of range (%d,%d)!\n",i,wp_lid,wp_wid,NLAY,NCEL);
            return -1;
        }
        if (wp_bid>=0&&wp_bid<NBRD&&wp_ch>=0&&wp_ch<NCHS){
            map_lid[wp_bid][wp_ch] = wp_lid;
            map_wid[wp_bid][wp_ch] = wp_wid;
        }
        else{
            fprintf(stderr,"ERROR: Entry %d in wiremap file, bid = %d ch = %d out of range (%d,%d)!\n",i,wp_bid,wp_ch,NBRD,NCHS);
            return -1;
        }
    }
    TFile_wirepos->Close();

    //==================Get Crosspoints==========================
    TFile * TFile_crosspoint = new TFile(HOME+"/info/crosspoint.root");
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

    //===================Prepare XT curves==============================
    printf("##############XT##################\n");
    printf("Reading from %s/info/xt.%d.%s.root\n",HOME.Data(),runNo,prerunname.Data());
    TFile * i_xt = new TFile(HOME+Form("/info/xt.%d.",runNo)+prerunname+".root");
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
    c->Add(HOME+Form("/root/h_%d.root",runNo));
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

    //===================Prepare output ROOT file============================
    printf("Output file: %s/root/t_%d.%s.layer%d.root\n",HOME.Data(),runNo,runname.Data(),testlayer);
    TFile * of = new TFile(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),testlayer),"RECREATE"); 
    TTree * ot = new TTree("t","t");
    // from h_XXX
    ot->Branch("triggerNumber",&triggerNumber);
    ot->Branch("driftT",&i_driftT);
    ot->Branch("nHits",&i_nHits);
    ot->Branch("layerID",&i_layerID);
    ot->Branch("wireID",&i_wireID);
    ot->Branch("type",&i_type); // in dec, [IMASTR]. I: peak index (FIXME: not good for noisy data); M: peak index in a packet; A: smaller than aa cut? S: smaller than sum cut? T: -1 <tmin, 0 good, 1 >tmax; R: 0 center, 1 left, 2 right, 3 guard, 4 dummy
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
    if (!iEntryStop&&iEntryStart){iEntryStart = 0; iEntryStop=N-1;}
    if (memdebug>0) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
    for (Long64_t iEntry = iEntryStart; iEntry<=iEntryStop; iEntry++){
        if (debug>10) printf("#####################################\n");
        if (debug>10) printf("Entry %d\n",iEntry);
        if (memdebug>10) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
        else if (iEntry%100==0&&memdebug>1) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
        if (iEntry%100==0) std::cout<<iEntry<<std::endl;
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
        for (int ihit = 0; ihit<i_nHits; ihit++){
            (*i_driftT)[ihit]+=t0shift; // fix driftT according to t0shift
            double dt = (*i_driftT)[ihit];
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int statusl,statusr; // 1:  large t; -1: small t; 0: good t
            (*o_dxl)[ihit] = t2x(dt,lid,wid,-1,statusl);
            (*o_dxr)[ihit] = t2x(dt,lid,wid,1,statusr);
            int type = (*i_type)[ihit]; // IMASTR
            // R: region
            // keep the original defination
            // T: time
            if (dt<tmin) type+=3*10;
            else if (dt>tmax) type+=6*10;
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
            // S: sum of wave packet
            if ((*i_sum)[ihit]<sumCut) type+=1*100;
            // A: sum of full waveform
            if ((*i_aa)[ihit]<aaCut) type+=1*1000;
            // M: index of peak in a multi-peak wave packet
            type+=(*i_mpi)[ihit]*10000;
            // I: index of peak in the hole waveform
            type+=(*i_ip)[ihit]*100000; // FIXME: this is shadowing all the late arriving peaks! Should modify when we have serious noise problem
            (*i_type)[ihit] = type;
            if (lid != testlayer&&type<=3){ // good hit
                if (debug>11) printf("  Entry %d: dxl[%d][%d] = dxl[%d] = t2x(%.3e) = %.3e\n",iEntry,lid,wid,ihit,dt,(*o_dxl)[ihit]);
                if (debug>11) printf("  Entry %d: dxr[%d][%d] = dxr[%d] = t2x(%.3e) = %.3e\n",iEntry,lid,wid,ihit,dt,(*o_dxr)[ihit]);
                v_layer_ihit[lid].push_back(ihit);
                nHitsG++;
            }
        }

        // To find pair candidates
        int npairs = 0;
        int prelid = -1;
        for(int lid = 1; lid<NLAY-1; lid++){
            if (lid==testlayer||lid+1==testlayer) continue;
            if (v_layer_ihit[lid].size()>0&&v_layer_ihit[lid+1].size()>0){// both layers have hits
                if (prelid+1 != lid)v_pick_lid.push_back(lid);
                v_pick_lid.push_back(lid+1);
                prelid = lid;
                npairs++;
            }
        }
        for (int ipick = 0; ipick<v_pick_lid.size(); ipick++){
            if (debug>11) printf(" pick layer %d\n",v_pick_lid[ipick]);
        }

        // Do tracking
        if (debug>11) printf("nHitsG = %d, npairs = %d\n",nHitsG,npairs);
        if (nHitsG<=nHitsMax&&npairs>=3){ // need at least 3 pairs to do the fitting; number of hits should be small to control the time cost
            N_found++;

            int nSelections = 0;
            Tracking(0,nSelections,iEntry); // 0 means starting from the 1st pick; nSelections is the number of possible choices by selecting one hit per layer;
            if (o_nHitsS[0]>=5){ // at least 5 hits to fit: NDF of 3-D track without field is 4
                N_good++;
            }
        }
        ot->Fill();
        if (memdebug>10) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
    }// end of event loop
    if (memdebug>0) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;

    ot->Write();
    of->Close();
    printf("Triggered Events: %d\n",N_trigger);
    printf("Found Events: %d\n",N_found);
    printf("Good Events: %d\n",N_good);
    return 0;
}

int Tracking(int ipick,int & iselection,int iEntry){
    if (ipick == v_pick_lid.size()){ // finished picking hits
        if (debug>11) printf(" Finished picking selection %d:\n",iselection);
        if (memdebug>11) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
        doFitting(v_pick_lid.size(),iEntry,iselection);
        if (memdebug>11) std::cout<<"Memory size @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize()<<std::endl;
        iselection++;
    }
    else{
        int lid = v_pick_lid[ipick];
        for (int i = 0; i<v_layer_ihit[lid].size(); i++){
            int ihit = v_layer_ihit[lid][i];
            int wid = (*i_wireID)[ihit];
            if (debug>11) printf(" => pick # %d, layer %d, wire %d, hit[%d], ihit = %d\n",ipick,lid,wid,i,ihit);
            pick_ihit[ipick] = v_layer_ihit[lid][i];
            Tracking(ipick+1,iselection,iEntry);
        }
    }
    return 0;
}

int doFitting(int nPicks,int iEntry,int iselection){
    int ncombi = pow(2,nPicks);
    if (debug>11) printf("  %d picked layers -> %d combinations\n",nPicks,ncombi);
    for (int icombi = 0; icombi<ncombi; icombi++){ // each combination corresponds to a unique left/right selection set
    	o_nFind++;
        if (debug>1) printf("     combi %d\n",icombi);
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
        if (debug>11) printf("       1st RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",islx,sciYup,iinx,chi2x,islz,sciYup,iinz,chi2z);
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
        if (debug>11) printf("       2nd RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",islx,sciYup,iinx,chi2x,islz,sciYup,iinz,chi2z);
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
        if (debug>11) printf("       3rd RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",islx,sciYup,iinx,chi2x,islz,sciYup,iinz,chi2z);

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
                if (fabs(calD-dd)<2&&testlayer!=(*i_layerID)[ihit]&&type<=3){ // FIXME: should tune the error limit
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
            if (debug>11) printf("       good! nHitsSel = %d\n",nHitsSel);
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
                    if (fabs(fitD-dd)<1&&testlayer!=(*i_layerID)[ihit]&&type<=3){ // FIXME: should tune the error limit
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
                    inScint = checkScintillator(1.5,inx,slx,inz,slz); // FIXME
                    fromSource = slx>-beamSlxMax&&slx<beamSlxMax&&slz>-beamSlzMax&&slz<beamSlzMax;
                    if (inScint&&fromSource){
                        // update chi2
                        getchi2(chi2i,chi2pi,chi2ai,islx,iinx,islz,iinz,true);
                        getchi2(chi2,chi2p,chi2a,slx,inx,slz,inz,true);
                        // check chi2 and see where the result fits
                        if (debug>11) printf("         final RESULT: x=%.3e*(y-%.3e)+%.3e, z=%.3e*(y-%.3e)+%.3e, chi2i = %.3e chi2 = %.3e\n",slx,sciYup,inx,slz,sciYup,inz,chi2i,chi2);
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
                            if (debug>12){
								int type = getHitType((*i_type)[ihit],fitD>=0);
                                if (type<=3)
                                    printf("        %d (%d,%d) dd %.3e fd %.3e res %.3e\n",ihit,lid,wid,(*t_driftD)[ihit],fitD,fitD-(*t_driftD)[ihit]);
                                else 
                                    printf("              # %d (%d,%d) dd %.3e fd %.3e res %.3e\n",ihit,lid,wid,(*t_driftD)[ihit],fitD,fitD-(*t_driftD)[ihit]);
							}
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
//        if (debug>11) printf("                  xc = %.3e+%.3e*sin(%.3e)/(-sin(%.3e-%.3e))+%.3e*sin(%.3e)/sin(%.3e-%.3e)\n",mcp_xc[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2);
 //       if (debug>11) printf("                  zc = %.3e+%.3e*cos(%.3e)/(-sin(%.3e-%.3e))+%.3e*cos(%.3e)/sin(%.3e-%.3e)+%.3e\n",mcp_zc[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2,zc_fix_slx,zc_fix_slx);
        if (debug>11) printf("       cp[%d,%d]: w(%d,%d) i(%d,%d) dd(%f,%f) xyz(%f,%f,%f)\n",lid,ljd,wid,wjd,ihit,jhit,dd1,dd2,xc,(pick_wy[ipick+1]+pick_wy[ipick])/2.,zc);
        if (zc<-chamberHL||zc>chamberHL){
            if (debug>11) printf("       bad combination!\n");
            break;
        }
        nPairs++;
    }
    if (ipick==nPicks-1){
        if (debug>11) printf("       GOOD!\n");
        return 0;
    }
    else{
        if (debug>11) printf("       BAD @ %d!\n",ipick);
        return 1;
    }
}

int setErrors(int nPairs, bool noError){
    if (g_z) delete g_z;
    if (g_x) delete g_x;
    g_z = new TGraphErrors(nPairs,&(pair_wy[0]),&(pair_wz[0]),0,0);
    g_x = new TGraphErrors(nPairs,&(pair_wy[0]),&(pair_wx[0]),0,0);
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
        if (debug>=12) printf("           setErrors pair[%d]: error x = ->10, error z = %.3e->10\n",errorzMax0_i,errorzMax0);
    }
    if (errorzMax1>4) {
        g_z->SetPointError(errorzMax1_i,0,10);
        g_x->SetPointError(errorzMax1_i,0,10);
        if (debug>=12) printf("           setErrors pair[%d]: error x = ->10, error z = %.3e->10\n",errorzMax1_i,errorzMax1);
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
        if (debug>=12) printf("           getChi2XZ pair[%d]: error x = %.3e, error z = %.3e\n",ipair,tchi2x,tchi2z);
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
    if (debug>12) printf("              top xz(%.3e, %.3e) bottom xz(%.3e, %.3e)\n",xtop,ztop,xbot,zbot);
    if (xtop>sciHW||xtop<-sciHW||xbot>sciHW||xbot<-sciHW||ztop>sciHL||ztop<-sciHL||zbot>sciHL||zbot<-sciHL) return false;
    else return true;
}

bool checkChi2(int nHitsSel, int nPairs, int icombi, int iselection){
	bool issame = false;
	bool covered = false;
    for (int i = 0; i<NCAND; i++){
    	issame = isSame(i);
    	if (issame){ // yes, there is a candidate with the same hits
            if (debug>11) printf(" same with Cand#%d where chi2=%.3e\n",i,o_chi2[i]);
//    		if (chi2<o_chi2[i]){// better? then remove the old one
			// FIXME: WARNING, now we rely on total chi2 including test layer hit, a slight bias
    		if (chi2a<o_chi2a[i]){// better? then remove the old one 
                if (debug>11) printf("   better than Cand#%d\n",i);
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
			    if (debug>11) printf("better than Cand#%d where chi2=%.3e, nHitsS=%d\n",i,o_chi2[i],o_nHitsS[i]);
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
    if (workType==0) theLayer = 0;
	else if (workType==1) theLayer = (lid%2==0?NLAY:NLAY+1); // even/odd
	else if (workType==-1) theLayer = (lid%2==0?NLAY+1:NLAY); // even/odd reversed
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

	// Set starting values and step sizes for parameters
	gMinuit->mnparm(0, "slopeX", sliX, beamSlxMax/1.e4, -beamSlxMax,beamSlxMax,ierflg);
	gMinuit->mnparm(1, "interceptX", iniX, beamInxMax/1.e4, -beamInxMax,beamInxMax,ierflg);
	gMinuit->mnparm(2, "slopeZ", sliZ, beamSlzMax/1.e4, -beamSlzMax,beamSlzMax,ierflg);
	gMinuit->mnparm(3, "interceptZ", iniZ, beamInzMax/1.e4, -beamInzMax,beamInzMax,ierflg);

	// Now ready for minimization step
	arglist[0] = 500.0;
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
		// FIXME: we should consider about the error
//		double error = errord[(*i_layerID)[ihit]][(*i_wireID)[ihit]]*(fabs(fabs(dfit)-4)+2/1.5)*1.5/4;
//		double error = errord[(*i_layerID)[ihit]][(*i_wireID)[ihit]];
//		double error = funcErr->Eval(fabs(dfit));
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
			if ((*i_layerID)[ihit]!=testlayer) continue;
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
			double error = 0.2; // FIXME: may consider non-constant error
			ca = f*N+minres*minres/error/error;
			ca/=(N+1);
		}
	}
}

double getError(int lid,double dt, bool isR){
	double error = 0.2; // default value 200 um
	TGraph * gr = 0;
    int theLayer = lid;
    if (workType==0) theLayer = 0;
	else if (workType==1) theLayer = (lid%2==0?NLAY:NLAY+1); // even/odd
	else if (workType==-1) theLayer = (lid%2==0?NLAY+1:NLAY); // even/odd reversed
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
    fprintf(stderr,"\t%s [runNo] [testlayer] [prerunname] [runname] <[nHitsMax] [t0shift] [tmin] [tmax] [geoSetup] [sumCut] [aaCut] [iEntryStart] [iEntryStop] [workType] [debug] [memdebug]>\n",prog_name);
}
