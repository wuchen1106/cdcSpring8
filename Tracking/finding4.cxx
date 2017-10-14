#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TString.h"
#include "TGraphErrors.h"

#define XMAX 130
#define ZMAX 350

#define NLAY 9
#define NZXP 8
#define NCEL 11
#define NCELA 88

#define NBRD 2
#define NCHS 48
#define NCHT 96
#define NSMP 32
#define ADCMIN 125
#define ADCMAX 700

#define NCUT 50

#define NCAND 4

int debug = 0;

std::vector<int> * i_layerID = 0;
std::vector<int> * i_wireID = 0;
std::vector<double> * o_dxl = 0;
std::vector<double> * o_dxr = 0;

//===================Chamber Parameter============================
double U = 8; // mm
double chamberHL = 599.17/2; // mm
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

//===================Hit Handlers============================
std::vector<std::vector<int> > v_layer_ihit; // vectors of indices of hits in each layer
std::vector<int> v_pick_lid; // layer index of each picked layer (for pair check)
std::vector<int> pick_ihit; // hit index of each picked hit (for pair check)
std::vector<double> pick_wy; // y position of each picked hit (for pair check)
std::vector<double> pair_wx; // x position of each picked hit pair (center position)
std::vector<double> pair_wy; // y position of each picked hit pair (center position)
std::vector<double> pair_wz; // z position of each picked hit pair (center position)

//===================About tracking============================
double yup = 623.97007;
double ydown = 527.60011;
TF1 * f_x = new TF1("f_x","pol1",500,640); // x VS y
TGraphErrors * g_x = 0; // x VS y
TF1 * f_z = new TF1("f_z","pol1",500,640); // z VS y
TGraphErrors * g_z = 0; // z VS y
double chi2_z = 0;
double inz = 0;
double slz = 0;
double chi2_x = 0;
double inx = 0;
double slx = 0;

//===================About beam============================
double slzMax = 0.2;
double slxMax = 0.05;

//===================About scintillator============================
double ytop = yup+50+8; // top scintillator
double ybot = ydown-50-8; // bottom scintillator
double hLtri = 150; // normal scintillator
double hHtri = 45; // normal scintillator

//===================Functions============================
void print_usage(char* prog_name);
double t2x(double time, int lid, int wid, int lr, int & status);
int ChooseHits(int ipick,int & iselection,int iEntry=0);
int checkCrossPoints(int nPicks,int iEntry=0,int iselection = 0);
int updatePairPositions(int icombi,int nPicks,int & nPairs);
int updateHitPositions(int nPicks,int icombi);
int setErrors(int nPairs, bool noError = false);
int getChi2XZ(int nPairs, double & chi2x, double & chi2z);
int fityx(int nPairs);
int fityz(int nPairs);
bool checkScintillator();

//===================About xt============================
double tres = 0;
TF1 * f_left_end[NCELA];
TF1 * f_right_end[NCELA];
TF1 * f_left[NCELA];
TF1 * f_right[NCELA];
int vtrel[NCELA];
int vtrer[NCELA];
int vtlel[NCELA];
int vtler[NCELA];
int vtrml[NCELA];
int vtrmr[NCELA];
int vtlml[NCELA];
int vtlmr[NCELA];

int main(int argc, char** argv){

    if (argc<3){
        print_usage(argv[0]);
        return 1;
    }
    int runNo = (int)strtol(argv[1],NULL,10);
    int testlayer = (int)strtol(argv[2],NULL,10);
    TString suffix = "";
    if (argc>=4){
        suffix  = argv[3];
        suffix="."+suffix;
    }
    int t0shift = 0;
    if (argc>=5) t0shift = (int)strtol(argv[4],NULL,10);
    int iEntryStart = 0;
    int iEntryStop = 0;
    if (argc>=6){
        iEntryStop = (int)strtol(argv[5],NULL,10)-1;
    }
    if (argc>=7){
        iEntryStart = (int)strtol(argv[5],NULL,10);
        iEntryStop = (int)strtol(argv[6],NULL,10);
    }
    int nHitsMax = 10;
    if (argc>=8){
        nHitsMax = (int)strtol(argv[7],NULL,10);
    }
    int tmin = -50;
    if (argc>=9){
        tmin = (int)strtol(argv[8],NULL,10);
    }
    int tmax = 1000; // FIXME: need to be properly set
    if (argc>=10){
        tmax = (int)strtol(argv[9],NULL,10);
    }
    if (argc>=11){
        debug = (int)strtol(argv[10],NULL,10);
    }
    printf("runNo       = %d\n",runNo);
    printf("test layer  = %d\n",testlayer);
    printf("suffix      = %s\n",suffix.Data());
    printf("t0shift     = %d\n",t0shift);
    printf("tmin        = %d\n",tmin);
    printf("tmax        = %d\n",tmax);
    printf("nHitsMax    = %d\n",nHitsMax);
    printf("Start Entry = %d\n",iEntryStart);
    printf("Stop Entry  = %d\n",iEntryStop);
    printf("debug       = %d\n",debug);

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
    TFile * TFile_wirepos = new TFile("../info/wire-position.root");
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
    TFile * TFile_crosspoint = new TFile("../info/crosspoint.root");
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
    TFile * i_xt = new TFile(Form("../info/xt.%d.root",runNo));
    for (int i = 0; i<NCELA; i++){
        f_left_end[i] = (TF1*) i_xt->Get(Form("f_left_end_%d_%d",i/NCEL+1,i%NCEL));
        f_right_end[i] = (TF1*) i_xt->Get(Form("f_right_end_%d_%d",i/NCEL+1,i%NCEL));
        f_left[i] = (TF1*) i_xt->Get(Form("f_left_%d_%d",i/NCEL+1,i%NCEL));
        f_right[i] = (TF1*) i_xt->Get(Form("f_right_%d_%d",i/NCEL+1,i%NCEL));
    }
    TTree * itree_xt = (TTree*) i_xt->Get("p");
    int trel,trer,tlel,tler,tlml,tlmr,trml,trmr;
    int lid,wid;
    itree_xt->SetBranchAddress("lid",&lid);
    itree_xt->SetBranchAddress("wid",&wid);
    itree_xt->SetBranchAddress("tlel",&tlel);
    itree_xt->SetBranchAddress("tler",&tler);
    itree_xt->SetBranchAddress("tlml",&tlml);
    itree_xt->SetBranchAddress("tlmr",&tlmr);
    itree_xt->SetBranchAddress("trml",&trml);
    itree_xt->SetBranchAddress("trmr",&trmr);
    itree_xt->SetBranchAddress("trel",&trel);
    itree_xt->SetBranchAddress("trer",&trer);
    if(debug>0) printf("lid,wid: tlel | tler | tlml | tlmr || trml | trmr | trel | trer\n",lid,wid,tlel,tler,tlml,tlmr,trml,trmr,trel,trer);
    for (int i = 0; i<itree_xt->GetEntries(); i++){
        itree_xt->GetEntry(i);
        vtlel[(lid)*NCEL+wid] = tlel;
        vtler[(lid)*NCEL+wid] = tler;
        vtlml[(lid)*NCEL+wid] = tlml;
        vtlmr[(lid)*NCEL+wid] = tlmr;
        vtrel[(lid)*NCEL+wid] = trel;
        vtrer[(lid)*NCEL+wid] = trer;
        vtrml[(lid)*NCEL+wid] = trml;
        vtrmr[(lid)*NCEL+wid] = trmr;
        if(debug>0) printf("%d,%d: %d | %d | %d | %d || %d | %d | %d | %d\n",lid,wid,tlel,tler,tlml,tlmr,trml,trmr,trel,trer);
    }

    //===================Get input ROOT file============================
    TChain * c = new TChain("t","t");
    c->Add(Form("../root/h_%d",runNo)+suffix+".root");
    int triggerNumber;
    int i_nHits;
    std::vector<double> * i_driftT = 0;
    std::vector<int> * i_type = 0;
    std::vector<int> * i_np = 0;
    std::vector<int> * i_ip = 0;
    std::vector<int> * i_clk = 0;
    std::vector<int> * i_width = 0;
    std::vector<int> * i_peak = 0;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_aa = 0;
    c->SetBranchAddress("triggerNumber",&triggerNumber);
    c->SetBranchAddress("nHits",&i_nHits);
    c->SetBranchAddress("driftT",&i_driftT);
    c->SetBranchAddress("layerID",&i_layerID);
    c->SetBranchAddress("wireID",&i_wireID);
    c->SetBranchAddress("type",&i_type); // -1: dummy layer; 1: guard layer; 2: left end; 3: right end; 0: center cell;
    c->SetBranchAddress("np",&i_np);
    c->SetBranchAddress("ip",&i_ip);
    c->SetBranchAddress("clk",&i_clk);
    c->SetBranchAddress("width",&i_width);
    c->SetBranchAddress("peak",&i_peak);
    c->SetBranchAddress("sum",&i_sum);
    c->SetBranchAddress("aa",&i_aa);

    //===================Prepare output ROOT file============================
    TFile * of = new TFile(Form("../root/i_%d.layer%d",runNo,testlayer)+suffix+".root","RECREATE"); 
    TTree * ot = new TTree("t","t");
    std::vector<double> * o_driftD = 0;
    std::vector<std::vector<double> > * o_fitD = 0;
    double o_tx1[NCAND];
    double o_tx2[NCAND];
    double o_tz1[NCAND];
    double o_tz2[NCAND];
    double o_chi2[NCAND];
    // from h_XXX
    ot->Branch("triggerNumber",&triggerNumber);
    ot->Branch("driftT",&i_driftT);
    ot->Branch("nHits",&i_nHits);
    ot->Branch("layerID",&i_layerID);
    ot->Branch("wireID",&i_wireID);
    ot->Branch("type",&i_type); // -1: dummy layer; 1(11): guard layer (out of t window); 2(12): left end (out of t window); 3(13): right end (out of t window); 0(10): center cell (out of t window);
    ot->Branch("np",&i_np);
    ot->Branch("ip",&i_ip);
    ot->Branch("clk",&i_clk);
    ot->Branch("width",&i_width);
    ot->Branch("peak",&i_peak);
    ot->Branch("sum",&i_sum);
    ot->Branch("aa",&i_aa);
    // 
    ot->Branch("driftD",&o_driftD);
    ot->Branch("dxl",&o_dxl);
    ot->Branch("dxr",&o_dxr);
    // after track finding, with different candidates;
    ot->Branch("fitD",&o_fitD); // [NCAND][nHits]
    ot->Branch("tx1",o_tx1);
    ot->Branch("tx2",o_tx2);
    ot->Branch("tz1",o_tz1);
    ot->Branch("tz2",o_tz2);
    ot->Branch("chi2",o_chi2);
    o_fitD = new std::vector<std::vector<double> >;
    for (int i = 0; i<NCAND; i++){
        std::vector<double> fitDs;
        o_fitD->push_back(fitDs);
    }
    o_driftD = new std::vector<double>;
    o_dxl = new std::vector<double>;
    o_dxr = new std::vector<double>;

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
    int nHitsgood; // number of good hits
    int nHitLayers; // number of layers with hit

    //===================Track Finding============================
    Long64_t N = c->GetEntries();
    for ( int iEntry = iEntryStart; iEntry<=iEntryStop; iEntry++){
        if (iEntry%1000==0) std::cout<<iEntry<<std::endl;
        c->GetEntry(iEntry);
        N_trigger++; // triggered event

        //========================================================================================================
        // prepare
        o_driftD->clear();
        o_dxl->clear();
        o_dxr->clear();
        for (int i = 0; i<NCAND; i++){
            (*o_fitD)[i].clear();
        }

        for (int i = 0; i<NLAY; i++){
            v_layer_ihit[i].clear();
        }
        nHitsgood = 0;
        nHitLayers = 0;
        v_pick_lid.clear();


        //========================================================================================================
        // get basical cdc hit information
        for (int ihit = 0; ihit<i_type->size(); ihit++){
            double dt = (*i_driftT)[ihit]+t0shift;
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            int status; // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
            o_dxl->push_back(t2x(dt,lid,wid,-1,status));
            o_dxr->push_back(t2x(dt,lid,wid,1,status));
            if (debug>0) printf("  Entry %d: dxl[%d][%d] = dxl[%d] = t2x(%.3e) = %.3e\n",iEntry,lid,wid,ihit,dt,(*o_dxl)[o_dxl->size()-1]);
            if (debug>0) printf("  Entry %d: dxr[%d][%d] = dxr[%d] = t2x(%.3e) = %.3e\n",iEntry,lid,wid,ihit,dt,(*o_dxr)[o_dxr->size()-1]);
            o_driftD->push_back((*o_dxr)[o_dxr->size()-1]);
            if (dt<tmin||dt>tmax) (*i_type)[ihit]+=10; // -1: dummy layer; 1(11): guard layer (out of t window); 2(12): left end (out of t window); 3(13): right end (out of t window); 0(10): center cell (out of t window);
            if (lid != testlayer&&(*i_type)[ihit]>=0&&(*i_type)[ihit]<10){ // not dummy layer or test layer; within time window
                v_layer_ihit[lid].push_back(ihit);
                nHitsgood++;
            }
        }

        // get number of layers with good hits;
        if (debug>0) printf("#####################################\n");
        if (debug>0) printf("Entry %d\n",iEntry);
        for(int lid = 1; lid<NLAY; lid++){
            if (debug>0)printf(" layer %d: %d hits\n",lid,v_layer_ihit[lid].size());
            if (v_layer_ihit[lid].size()>0)
                nHitLayers++;
        }

        if (nHitsgood>nHitsMax||nHitLayers<6) continue; // Need at least 6 layers with good hit and no more than nHitsMax good hits in total
        N_found++;

        //========================================================================================================
        // To find pair candidates
        int prelid = -1;
        for(int lid = 1; lid<NLAY-1; lid++){
            if (lid==testlayer||lid+1==testlayer) continue;
            if (v_layer_ihit[lid].size()>0&&v_layer_ihit[lid+1].size()>0){// both layers have hits
                if (prelid+1 != lid)v_pick_lid.push_back(lid);
                v_pick_lid.push_back(lid+1);
                prelid = lid;
            }
        }
        for (int ipick = 0; ipick<v_pick_lid.size(); ipick++){
            if (debug>0) printf(" pick layer %d\n",v_pick_lid[ipick]);
        }

        //========================================================================================================
        // pick up one hit from each layer, and iterate in all combinations
        int nSelections = 0;
        ChooseHits(0,nSelections,iEntry);

        ot->Fill();
    }// end of event loop

    ot->Write();
    of->Close();
    printf("Triggered Events: %d\n",N_trigger);
    printf("Found Events: %d\n",N_found);
    printf("Good Events: %d\n",N_good);
    return 0;
}

int ChooseHits(int ipick,int & iselection,int iEntry){
    if (ipick == v_pick_lid.size()){ // finished picking hits
        if (debug>0) printf(" Finished picking selection %d:\n",iselection);
        checkCrossPoints(v_pick_lid.size(),iEntry,iselection);
        iselection++;
    }
    else{
        int lid = v_pick_lid[ipick];
        for (int i = 0; i<v_layer_ihit[lid].size(); i++){
            int ihit = v_layer_ihit[lid][i];
            int wid = (*i_wireID)[ihit];
            if (debug>0) printf(" => pick # %d, layer %d, wire %d, hit[%d], ihit = %d\n",ipick,lid,wid,i,ihit);
            pick_ihit[ipick] = v_layer_ihit[lid][i];
            ChooseHits(ipick+1,iselection,iEntry);
        }
    }
    return 0;
}

int checkCrossPoints(int nPicks,int iEntry,int iselection){
    int ncombi = pow(2,nPicks);
    if (debug>0) printf("  %d picked layers -> %d combinations\n",nPicks,ncombi);
    for (int icombi = 0; icombi<ncombi; icombi++){ // each combination corresponds to a unique left/right selection set
        f_x->SetParameters(0,0);
        f_z->SetParameters(0,0);
        int nPairs = 0;

        updateHitPositions(nPicks,icombi); // fix wy positions
        int result = updatePairPositions(icombi,nPicks,nPairs);
        if (result) continue;
        setErrors(nPairs,true);
        fityz(nPairs);
        fityx(nPairs);
        bool inScint = checkScintillator();
        inz = f_z->Eval(yup);
        slz = f_z->GetParameter(1);
        inx = f_x->Eval(yup);
        slx = f_x->GetParameter(1);
        bool fromSource = slx>-slxMax&&slx<slxMax&&slz>-slzMax&&slz<slzMax;
        int nGood = getChi2XZ(nPairs,chi2_x,chi2_z);
        if (debug>0) printf("       1st RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",slx,yup,inx,chi2_x,slz,yup,inz,chi2_z);
        if (debug>=0&&inScint&&fromSource) printf("%d %d 0 %d %d %.4e %.4e\n",iEntry,icombi,iselection,nGood,chi2_x,chi2_z);
        if (!fromSource||!inScint) {f_x->SetParameters(0,0); f_z->SetParameters(0,0);}
        
        updateHitPositions(nPicks,icombi); // fix wy positions
        result = updatePairPositions(icombi,nPicks,nPairs);
        setErrors(nPairs,false);
        if (result) continue;
        fityz(nPairs);
        fityx(nPairs);
        inScint = checkScintillator();
        inz = f_z->Eval(yup);
        slz = f_z->GetParameter(1);
        inx = f_x->Eval(yup);
        slx = f_x->GetParameter(1);
        fromSource = slx>-slxMax&&slx<slxMax&&slz>-slzMax&&slz<slzMax;
        nGood = getChi2XZ(nPairs,chi2_x,chi2_z);
        if (debug>0) printf("       2nd RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",slx,yup,inx,chi2_x,slz,yup,inz,chi2_z);
        if (debug>=0&&inScint&&fromSource) printf("%d %d 1 %d %d %.4e %.4e\n",iEntry,icombi,iselection,nGood,chi2_x,chi2_z);
        if (!fromSource||!inScint) {f_x->SetParameters(0,0); f_z->SetParameters(0,0);}
        
        updateHitPositions(nPicks,icombi); // fix wy positions
        result = updatePairPositions(icombi,nPicks,nPairs);
        setErrors(nPairs,false);
        if (result) continue;
        fityz(nPairs);
        fityx(nPairs);
        inScint = checkScintillator();
        inz = f_z->Eval(yup);
        slz = f_z->GetParameter(1);
        inx = f_x->Eval(yup);
        slx = f_x->GetParameter(1);
        fromSource = slx>-slxMax&&slx<slxMax&&slz>-slzMax&&slz<slzMax;
        nGood = getChi2XZ(nPairs,chi2_x,chi2_z);
        if (debug>0) printf("       3rd RESULT: nGood = %d, inScint? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e\n",nGood,inScint?"yes":"no",slx,yup,inx,chi2_x,slz,yup,inz,chi2_z);
        if (debug>=0&&inScint&&fromSource) printf("%d %d 2 %d %d %.4e %.4e\n",iEntry,icombi,iselection,nGood,chi2_x,chi2_z);
    }
}

int updateHitPositions(int nPicks,int icombi){
    // calculate pick_wy
    for (int ipick = 0; ipick<nPicks; ipick++){
        // Get hit information
        int ilr = (icombi&(1<<ipick))>>ipick;
        int ihit = pick_ihit[ipick];
        int lid = (*i_layerID)[ihit];
        int wid = (*i_wireID)[ihit];
        double dd;
        if (ilr) dd = (*o_dxr)[ihit];
        else dd = (*o_dxl)[ihit];
        double wxro = map_x[lid][wid][1];
        double wyro = map_y[lid][wid][1];
        double wzro = chamberHL;
        double wxhv = map_x[lid][wid][0];
        double wyhv = map_y[lid][wid][0];
        double wzhv = -chamberHL;
        // assume wy
        double wy = (wyro+wyhv)/2.;
        // get wz by extrapolating the track to wy
        // FIXME: do we have to consider f_x here?
        double wz = f_z->Eval(wy);
        // correct wy according to wz
        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
        pick_wy[ipick] = wy;
    }
}

int updatePairPositions(int icombi,int nPicks,int & nPairs){
    // calculate pair_wxyz
    if (debug>1) printf("     combi %d\n",icombi);
    nPairs = 0;
    int ipick = 0;
    for (; ipick<nPicks-1; ipick++){
        int ihit = pick_ihit[ipick];
        int jhit = pick_ihit[ipick+1];
        int lid = (*i_layerID)[ihit];
        int ljd = (*i_layerID)[jhit];
        if (lid+1!=ljd) continue; // not adjacent
        int ilr = (icombi&(1<<ipick))>>ipick;
        int jlr = (icombi&(1<<(ipick+1)))>>(ipick+1);
        double deltaY = pick_wy[ipick+1]-pick_wy[ipick];
        double dd1, dd2;
        if (ilr) dd1 = (*o_dxr)[ihit]; // right
        else     dd1 = (*o_dxl)[ihit]; // left
        if (jlr) dd2 = (*o_dxr)[jhit]; // right
        else     dd2 = (*o_dxl)[jhit]; // left
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
//        if (debug>1) printf("                  xc = %.3e+%.3e*sin(%.3e)/(-sin(%.3e-%.3e))+%.3e*sin(%.3e)/sin(%.3e-%.3e)\n",mcp_xc[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2);
 //       if (debug>1) printf("                  zc = %.3e+%.3e*cos(%.3e)/(-sin(%.3e-%.3e))+%.3e*cos(%.3e)/sin(%.3e-%.3e)+%.3e\n",mcp_zc[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2,zc_fix_slx,zc_fix_slx);
        if (debug>1) printf("       cp[%d,%d]: lr(%d,%d) w(%d,%d) i(%d,%d) dd(%f,%f)] xyz(%f,%f,%f)\n",lid,ljd,ilr,jlr,wid,wjd,ihit,jhit,dd1,dd2,xc,(pick_wy[ipick+1]+pick_wy[ipick])/2.,zc);
        if (zc<-chamberHL||zc>chamberHL){
            if (debug>1) printf("       bad combination!\n");
            break;
        }
        nPairs++;
    }
    if (ipick==nPicks-1){
        if (debug>1) printf("       GOOD!\n");
        return 0;
    }
    else{
        if (debug>1) printf("       BAD @ %d!\n",ipick);
        return 1;
    }
}

int setErrors(int nPairs, bool noError){
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
        if (debug>=2) printf("           setErrors pair[%d]: error x = ->10, error z = %.3e->10\n",errorzMax0_i,errorzMax0);
    }
    if (errorzMax1>4) {
        g_z->SetPointError(errorzMax1_i,0,10);
        g_x->SetPointError(errorzMax1_i,0,10);
        if (debug>=2) printf("           setErrors pair[%d]: error x = ->10, error z = %.3e->10\n",errorzMax1_i,errorzMax1);
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
        if (debug>=2) printf("           getChi2XZ pair[%d]: error x = %.3e, error z = %.3e\n",ipair,tchi2x,tchi2z);
	    if (fabs(tchi2z)<4&&fabs(tchi2x)<1){
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
	for (int ipair = 0; ipair<nPairs; ipair++){
	    g_z->SetPoint(ipair,pair_wy[ipair],pair_wz[ipair]);
	}
	g_z->Set(nPairs);
	g_z->Fit("f_z","qN0F","");
	return 0;
}

int fityx(int nPairs){
	for (int ipair = 0; ipair<nPairs; ipair++){
	    g_x->SetPoint(ipair,pair_wy[ipair],pair_wx[ipair]);
	}
	g_x->Set(nPairs);
	g_x->Fit("f_x","qN0F","");
	return 0;
}

bool checkScintillator(){
    double xtop = f_x->Eval(ytop);
    double xbot = f_x->Eval(ybot);
    double ztop = f_z->Eval(ytop);
    double zbot = f_z->Eval(ybot);
    if (debug>=3) printf("              top xz(%.3e, %.3e) bottom xz(%.3e, %.3e)\n",xtop,ztop,xbot,zbot);
    if (xtop>hHtri||xtop<-hHtri||xbot>hHtri||xbot<-hHtri||ztop>hLtri||ztop<-hLtri||zbot>hLtri||zbot<-hLtri) return false;
    else return true;
}

double t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
    TF1* fl=0;
    TF1* fr=0;
    // FIXME: now we only take one xt: layer 5 cell 0 (fake)
    //int index = (lid-1)*NCEL+wid;
    int index = (6-1)*NCEL;
    status = 0;
    // FIXME: should we really set a boundary of f_left_end/f_right_end? driftT can be really long when it comes from corner...
//    if (time<=vtlel[index]&&time>vtler[index]){
    if (time>vtler[index]){
        fl = f_left_end[index];
        status = -2;
    }
    else if (time<=vtlml[index]&&time>=vtlmr[index]){
        fl = f_left[index];
        status = -1;
    }
//    if (time>=vtrel[index]&&time<vtrer[index]){
    if (time>vtrel[index]){
        fr = f_right_end[index];
        status = 2;
    }
    else if (time>=vtrml[index]&&time<=vtrmr[index]){
        fr = f_right[index];
        status = 1;
    }
    double dd=0;
    if (lr>=0){
        if (fr) dd = fr->Eval(time);
    }
    else{
        if (fl) dd = fl->Eval(time);
    }
    //std::cout<<"t2x("<<time<<","<<lid<<","<<wid<<","<<lr<<") = "<<dd<<std::endl;
    return dd;
}

//______________________________________________________________________________
void print_usage(char* prog_name)
{
    fprintf(stderr,"\t%s [runNo] [testlayer] <[suffix] [t0shift] [iEntryStart] [iEntryStop] [nHitsMax] [tmin] [tmax] [debug]>\n",prog_name);
}
