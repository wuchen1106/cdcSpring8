#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TString.h"
#include "TGraph.h"

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

#define NHITSMAX 10

std::vector<int> * i_layerID = 0;
std::vector<int> * i_wireID = 0;
std::vector<double> * o_driftD = 0;
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
int     map_lid[NBRD][NCHS];
int     map_wid[NBRD][NCHS];
int     map_theta[NBRD][NCHS];
// mcp for cross points
double mcp_xc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
double mcp_zc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 

//===================About tracking============================
double yup = 623.97007;
double ydown = 527.60011;
double zup;
double zdown;
double xup;
double xdown;
double slx;
double slz;
TF1 * f_x = new TF1("f_x","pol1",500,640); // x VS y
TGraph * g_x = 0; // x VS y
TF1 * f_z = new TF1("f_z","pol1",500,640); // z VS y
TGraph * g_z = 0; // z VS y

//===================Functions============================
void print_usage(char* prog_name);
double t2x(double time, int lid, int wid, int lr, int & status);
int ChooseHits(int ipick, std::vector<int> & v_pick_lid, std::vector<int> & pick_ihit, std::vector<double> & pick_wy, std::vector<std::vector<int> > & v_layer_ihit,int & iselection);
int checkCrossPoints(int nPicks,std::vector<int> & pick_ihit, std::vector<double> & pick_wy);
int getCrossPoints(int ihit,int jhit,int ilr,int jlr,double deltaY);
int updateHitPositions(int nPicks, std::vector<int> & pick_ihit, std::vector<double> & pick_wy);
int fityx();
int fityz();

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
    printf("runNo       = %d\n",runNo);
    printf("test layer  = %d\n",testlayer);
    printf("suffix      = %s\n",suffix.Data());
    printf("t0shift     = %d\n",t0shift);
    printf("Start Entry = %d\n",iEntryStart);
    printf("Stop Entry  = %d\n",iEntryStop);
    double tmin = -50;
    double tmax = 380; // FIXME: need to be properly set

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
    printf("lid,wid: tlel | tler | tlml | tlmr || trml | trmr | trel | trer\n",lid,wid,tlel,tler,tlml,tlmr,trml,trmr,trel,trer);
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
        printf("%d,%d: %d | %d | %d | %d || %d | %d | %d | %d\n",lid,wid,tlel,tler,tlml,tlmr,trml,trmr,trel,trer);
    }

    //===================Get input ROOT file============================
    TChain * c = new TChain("t","t");
    c->Add(Form("../root/h_%d.root",runNo));
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
    c->SetBranchAddress("type",&i_type);
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
    ot->Branch("type",&i_type);
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

    //===================Hit handlers============================
    std::vector<std::vector<int> > v_layer_ihit; // vectors of indices of hits in each layer
    for (int i = 0; i<NLAY; i++){
        std::vector<int> hits;
        v_layer_ihit.push_back(hits);
    }
    std::vector<int> v_pick_lid; // layer index of each picked layer (for pair check)
    std::vector<int> pick_ihit; // hit index of each picked hit (for pair check)
    std::vector<double> pick_wy; // y position of each picked hit pair (for pair check)
    pick_ihit.resize(NLAY);
    pick_wy.resize(NLAY);
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
            if ((*i_type)[ihit]!=0&&(*i_type)[ihit]!=1) continue; // including guard layer but without dummy layer
            double dt = (*i_driftT)[ihit]+t0shift;
            int lid = (*i_layerID)[ihit];
            int wid = (*i_wireID)[ihit];
            if (dt<tmin||dt>tmax) continue; // ignore the hits which is obviously out of the signal time window
            int status; // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
            o_dxl->push_back(t2x(dt,lid,wid,-1,status));
            o_dxr->push_back(t2x(dt,lid,wid,1,status));
            o_driftD->push_back((*o_dxr)[o_dxr->size()-1]);
            v_layer_ihit[lid].push_back(ihit);
            if (lid != testlayer) nHitsgood++;
        }

        // get number of layers with good hits;
        printf("#####################################\n");
        printf("Entry %d\n",iEntry);
        for(int lid = 1; lid<NLAY; lid++){
            printf(" layer %d: %d hits\n",lid,v_layer_ihit[lid].size());
            if (lid!=testlayer&&v_layer_ihit[lid].size()>0)
                nHitLayers++;
        }

        if (nHitsgood>NHITSMAX||nHitLayers<6) continue; // Need at least 6 layers with good hit and no more than NHITSMAX good hits in total
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
            printf(" pick layer %d\n",v_pick_lid[ipick]);
        }

        //========================================================================================================
        // pick up one hit from each layer, and iterate in all combinations
        int nSelections = 0;
        ChooseHits(0,v_pick_lid,pick_ihit,pick_wy,v_layer_ihit,nSelections);

        ot->Fill();
    }// end of event loop

    ot->Write();
    of->Close();
    printf("Triggered Events: %d\n",N_trigger);
    printf("Found Events: %d\n",N_found);
    printf("Good Events: %d\n",N_good);
    return 0;
}

int ChooseHits(int ipick, std::vector<int> & v_pick_lid, std::vector<int> & pick_ihit, std::vector<double> & pick_wy, std::vector<std::vector<int> > & v_layer_ihit,int & iselection){
    if (ipick == v_pick_lid.size()){ // finished picking hits
        printf(" Finished picking selection %d:\n",iselection);
        checkCrossPoints(v_pick_lid.size(),pick_ihit,pick_wy);
        iselection++;
    }
    else{
        int lid = v_pick_lid[ipick];
        for (int i = 0; i<v_layer_ihit[lid].size(); i++){
            printf(" => pick # %d, layer %d, hit[%d] %d\n",ipick,lid,i,v_layer_ihit[lid][i]);
            pick_ihit[ipick] = v_layer_ihit[lid][i];
            ChooseHits(ipick+1,v_pick_lid,pick_ihit,pick_wy,v_layer_ihit,iselection);
        }
    }
    return 0;
}

int checkCrossPoints(int nPicks,std::vector<int> & pick_ihit, std::vector<double> & pick_wy){
    int ncombi = pow(2,nPicks);
    printf("  %d picked layers -> %d combinations\n",nPicks,ncombi);
    updateHitPositions(nPicks, pick_ihit,pick_wy);
    for (int icombi = 0; icombi<ncombi; icombi++){ // each combination corresponds to a unique left/right selection set
        printf("     combi %d\n",icombi);
        for (int ipick = 0; ipick<nPicks-1; ipick++){
            int ihit = pick_ihit[ipick];
            int jhit = pick_ihit[ipick+1];
            int lid = (*i_layerID)[ihit];
            int ljd = (*i_layerID)[jhit];
            if (lid+1!=ljd) continue; // not adjacent
            int ilr = (icombi&(1<<ipick))>>ipick;
            int jlr = (icombi&(1<<(ipick+1)))>>(ipick+1);
            printf("      %d: %d,%d(%d) %d,%d(%d)\n",ipick,lid,ihit,ilr,ljd,jhit,jlr);
            getCrossPoints(ihit,jhit,ilr,jlr,pick_wy[ipick+1]-pick_wy[ipick]);
        }
    }
}

int updateHitPositions(int nPicks, std::vector<int> & pick_ihit, std::vector<double> & pick_wy){
    // calculate pick_wy
    for (int ipick = 0; ipick<nPicks; ipick++){
        // Get hit information
        int ihit = pick_ihit[ipick];
        int lid = (*i_layerID)[ihit];
        int wid = (*i_wireID)[ihit];
        double dd = (*o_driftD)[ihit]; // !!! make sure it's mm
        double wxro = map_x[lid][wid][1];
        double wyro = map_y[lid][wid][1];
        double wzro = chamberHL;
        double wxhv = map_x[lid][wid][0];
        double wyhv = map_y[lid][wid][0];
        double wzhv = -chamberHL;
        // assume wy
        double wy = (wyro+wyhv)/2.;
        // get wz by extrapolating the track to wy
        double zdown = zup + (ydown-yup)*slz;
        double xdown = xup + (ydown-yup)*slx;
        double wz = ((yup-wy)*zdown+(wy-ydown)*zup)/(yup-ydown);
        // get wx according to wz
        double wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
        // correct wy according to wz
        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
        pick_wy[ipick] = wy;
    }
}

int getCrossPoints(int ihit,int jhit,int ilr,int jlr,double deltaY){
    double dd1, dd2;
    if (ilr) dd1 = (*o_dxr)[ihit]; // right
    else     dd1 = (*o_dxl)[ihit]; // left
    if (jlr) dd2 = (*o_dxr)[jhit]; // right
    else     dd2 = (*o_dxl)[jhit]; // left
    int lid = (*i_layerID)[ihit];
    int ljd = (*i_layerID)[jhit];
    int wid = (*i_wireID)[ihit];
    int wjd = (*i_wireID)[jhit];
    double theta1 = map_theta[lid][wid];
    double theta2 = map_theta[ljd][wjd];
    double sintheta12 = sin(theta1-theta2);
    double zc_fix_slx = deltaY*slx/(tan(theta2)-tan(theta1));
    double xc = mcp_xc[lid][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
    double zc = mcp_zc[lid][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
    fityz();
}

int fityz(){
}

int fityx(){
}

double t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
    TF1* fl=0;
    TF1* fr=0;
    // FIXME
    //int index = (lid-1)*NCEL+wid;
    int index = (6-1)*NCEL;
    status = 0;
    if (time<=vtlel[index]&&time>vtler[index]){
        fl = f_left_end[index];
        status = -2;
    }
    else if (time<=vtlml[index]&&time>=vtlmr[index]){
        fl = f_left[index];
        status = -1;
    }
    if (time>=vtrel[index]&&time<vtrer[index]){
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
    fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
