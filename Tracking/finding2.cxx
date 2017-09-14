#include <vector>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "TMinuit.h"
#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TCanvas.h"

#define NLAY 9
#define NCEL 11
#define NCHT 96
#define NCHS 48

#define use_extra

int pow2i[NLAY];

//===================About tracking============================
std::vector<double> * i_driftT = 0;
std::vector<int> * i_layerID = 0;
std::vector<int> * i_wireID = 0;
std::vector<int> * i_type = 0;
std::vector<int> * O_wireID = 0;
std::vector<int> * O_layerID = 0;
std::vector<int> * O_lr = 0;
std::vector<double> * O_dxl= 0;
std::vector<double> * O_dxr= 0;
std::vector<double> * O_driftD = 0;
double ar_x[50]; // x VS y
double ar_y[50]; // x VS y
double yup = 623.97007;
double ydown = 527.60011;

TF1 * f_x = new TF1("f_x","pol1",500,640); // x VS y
TGraph * g_x = 0; // x VS y

int testlayer = 0;
int nHitsgood;
std::vector<int> v_hit_index;
std::vector<int> v_hit_flag;
std::vector<double> v_hit_dd;

//===================About chamber============================
double U = 8; // mm
double Lchamebr = 599.17; // mm
double zhv = -Lchamebr/2;
double zro = Lchamebr/2;

//===================About wireposition============================
double map_xc[NLAY][NCEL];
double map_yc[NLAY][NCEL];
double map_xhv[NLAY][NCEL];
double map_yhv[NLAY][NCEL];
double map_xro[NLAY][NCEL];
double map_yro[NLAY][NCEL];

//===================About xt============================
double tres = 0;
TF1 * f_left_end[88];
TF1 * f_right_end[88];
TF1 * f_left[88];
TF1 * f_right[88];
int vtrel[88];
int vtrer[88];
int vtlel[88];
int vtler[88];
int vtrml[88];
int vtrmr[88];
int vtlml[88];
int vtlmr[88];

//===================About hits============================
int a_nHitsLayer[NLAY];

//______________________________________________________________________________
void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}

int checkLR(){
	double	chi2min = 1e9;
	double	meanxmin = 1e9;
	int thecombi = -1;
	int ncombi = pow(2,nHitsgood);
	for (int icombi = 0; icombi<ncombi; icombi++){
		double chi2 = 0;
		double xc = 0;
		double avx = 0;
		int usedhits = 0;
//		std::cout<<"#################################################"<<std::endl;
		for (int index= 0; index<v_hit_index.size(); index++){
			if (v_hit_flag[index]) continue;
			int ihit = v_hit_index[index];
			int lr0 = (icombi/pow2i[usedhits])%2;
			usedhits++;

			double dd0 = v_hit_dd[index];
			xc = map_xc[(*i_layerID)[ihit]][(*i_wireID)[ihit]];
			if (lr0) xc+=dd0;
			else xc-=dd0;
			avx+=xc;
		}
		avx/=nHitsgood;
		usedhits=0;
		for (int index= 0; index<v_hit_index.size(); index++){
			if (v_hit_flag[index]) continue;
			int ihit = v_hit_index[index];
			int lr0 = (icombi/pow2i[usedhits])%2;
			usedhits++;

			double dd0 = v_hit_dd[index];
			xc = map_xc[(*i_layerID)[ihit]][(*i_wireID)[ihit]];
			if (lr0) xc+=dd0;
			else xc-=dd0;
			double dx = xc-avx;
			chi2+=dx*dx;
		}
//		std::cout<<icombi<<": chi2 = "<<chi2<<", avx = "<<avx<<std::endl;
		if (chi2<chi2min){
			chi2min=chi2;
			thecombi=icombi;
			meanxmin=avx;
		}
	}
	nHitsgood=0;
	for (int index = 0; index<v_hit_index.size(); index++){
		int ihit = v_hit_index[index];
		double xc = map_xc[(*i_layerID)[ihit]][(*i_wireID)[ihit]];
		double xr = xc+v_hit_dd[index];
		double xl = xc-v_hit_dd[index];
		if (fabs(xr-meanxmin)>=1&&fabs(xl-meanxmin)>=1) v_hit_flag[index] = -2;
		else {
			if ((*i_layerID)[ihit]!=testlayer) nHitsgood++;
			if (fabs(xr-meanxmin)<1&&fabs(xl-meanxmin)<1) v_hit_flag[index] = 0;
			else if (fabs(xr-meanxmin)<1) v_hit_flag[index] = 1;
			else v_hit_flag[index] = -1;
		}
	}
//	std::cout<<"meanxmin = "<<meanxmin<<", nHitsgood = "<<nHitsgood<<std::endl;
	if (nHitsgood<5) return 1;
	return 0;
}

int fitxy(){
	for (int ihit = 0; ihit<O_lr->size(); ihit++){
		int index = v_hit_index[ihit];
		ar_x[ihit] = map_xc[(*i_layerID)[index]][(*i_wireID)[index]]+v_hit_dd[ihit];
		ar_y[ihit] = map_yc[(*i_layerID)[index]][(*i_wireID)[index]];
	}
	if (g_x) delete g_x; g_x = new TGraph(O_lr->size(),&(ar_y[0]),&(ar_x[0]));
	g_x->Fit("f_x","qN0","");
	return 0;
}

double t2x(double time, int lid, int wid, int lr, int & status){
	TF1* fl=0;
	TF1* fr=0;
	// FIXME
	//int index = (lid-1)*11+wid;
	int index = (6-1)*11;
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

int main(int argc, char** argv){

	if (argc<3){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	testlayer = (int)strtol(argv[2],NULL,10);
	std::string suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix="."+suffix;
	}
	int t0shift = 0;
	if (argc>=5) t0shift = (int)strtol(argv[4],NULL,10);
	printf("t0shift = %d\n",t0shift);
	int nEventMax = 0;
	if (argc>=6) nEventMax = (int)strtol(argv[5],NULL,10);

	for (int i=0; i<NLAY; i++){
		pow2i[i] = pow(2,i);
	}

	//===================Get wire position============================
	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_bid;
	int wp_wid;
	int wp_lid;
	int wp_ch;
	double wp_xhv;
	double wp_yhv;
	double wp_xro;
	double wp_yro;
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("b",&wp_bid);
	TTree_wirepos->SetBranchAddress("ch",&wp_ch);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xro",&wp_xro);
	TTree_wirepos->SetBranchAddress("yro",&wp_yro);
	TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<NLAY){
			map_xhv[wp_lid][wp_wid] = wp_xhv;
			map_yhv[wp_lid][wp_wid] = wp_yhv;
			map_xro[wp_lid][wp_wid] = wp_xro;
			map_yro[wp_lid][wp_wid] = wp_yro;
		}
	}
	for(int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
			map_xc[ilayer][iwire] = (map_xhv[ilayer][iwire]+map_xro[ilayer][iwire])/2.;
			map_yc[ilayer][iwire] = (map_yhv[ilayer][iwire]+map_yro[ilayer][iwire])/2.;
		}
	}

	//===================Get run info============================
	TFile * if_run = new TFile("../info/run-info.root");
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
		if (i_runNo == runNo) break;
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
	std::cout<<"runNo#"<<runNo<<": "<<gastype<<", "<<runGr<<", "<<duration<<", "<<HV<<" V, "<<THR<<" mV, "<<durationTime<<"sec"<<std::endl;

	//===================Prepare XT curves==============================
	TFile * i_xt = new TFile(Form("../info/xt.%d.root",runNo));
	for (int i = 0; i<88; i++){
		f_left_end[i] = (TF1*) i_xt->Get(Form("f_left_end_%d_%d",i/11+1,i%11));
		f_right_end[i] = (TF1*) i_xt->Get(Form("f_right_end_%d_%d",i/11+1,i%11));
		f_left[i] = (TF1*) i_xt->Get(Form("f_left_%d_%d",i/11+1,i%11));
		f_right[i] = (TF1*) i_xt->Get(Form("f_right_%d_%d",i/11+1,i%11));
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
	for (int i = 0; i<itree_xt->GetEntries(); i++){
		itree_xt->GetEntry(i);
		vtlel[(lid-1)*11+wid] = tlel;
		vtler[(lid-1)*11+wid] = tler;
		vtlml[(lid-1)*11+wid] = tlml;
		vtlmr[(lid-1)*11+wid] = tlmr;
		vtrel[(lid-1)*11+wid] = trel;
		vtrer[(lid-1)*11+wid] = trer;
		vtrml[(lid-1)*11+wid] = trml;
		vtrmr[(lid-1)*11+wid] = trmr;
		printf("%d,%d: %d | %d | %d | %d || %d | %d | %d | %d\n",lid,wid,tlel,tler,tlml,tlmr,trml,trmr,trel,trer);
	}

	//===================Get input ROOT file============================
	TChain * c = new TChain("t","t");
	char fileName[128];
	sprintf(fileName,"../root/h_%d.root",runNo);
	c->Add(fileName);
	int triggerNumber;
	int i_nLayers, i_nHits;
	c->SetBranchAddress("triggerNumber",&triggerNumber);
	c->SetBranchAddress("nLayers",&i_nLayers);
	c->SetBranchAddress("nHits",&i_nHits);
	c->SetBranchAddress("driftT",&i_driftT);
	c->SetBranchAddress("layerID",&i_layerID);
	c->SetBranchAddress("wireID",&i_wireID);
	c->SetBranchAddress("type",&i_type);

	//===================Prepare output ROOT file============================
	std::vector<double> * O_driftT = 0;
	double O_tx1;
	double O_tx2;
	double O_tz1;
	double O_tz2;
	double O_chi2;

	std::stringstream buf;
	buf<<"../root/i_"<<runNo<<".layer"<<testlayer<<suffix<<".root";
	TFile * of = new TFile(buf.str().c_str(),"RECREATE"); 
	TTree * ot = new TTree("t","t");

	ot->Branch("nHits",&nHitsgood);
	ot->Branch("triggerNumber",&triggerNumber);
	ot->Branch("tx1",&O_tx1);
	ot->Branch("tx2",&O_tx2);
	ot->Branch("tz1",&O_tz1);
	ot->Branch("tz2",&O_tz2);
	ot->Branch("driftT",&O_driftT);
	ot->Branch("driftD",&O_driftD);
	ot->Branch("lr",&O_lr);
	ot->Branch("dxl",&O_dxl);
	ot->Branch("dxr",&O_dxr);
	ot->Branch("wireID",&O_wireID);
	ot->Branch("layerID",&O_layerID);
	ot->Branch("chi2",&O_chi2);

	//===================Efficiency Counter============================
	int N_trigger = 0;
	int N_found = 0;
	int N_good = 0;

	//===================Track Finding============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	for ( int i = 0; i<N; i++){
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		c->GetEntry(i);
//		printf("*****triggerNumber %d*******\n",triggerNumber);
		if(O_lr) delete O_lr; O_lr= new std::vector<int>;
		if(O_driftT) delete O_driftT; O_driftT = new std::vector<double>;
		if(O_driftD) delete O_driftD; O_driftD = new std::vector<double>;
		if(O_wireID) delete O_wireID; O_wireID = new std::vector<int>;
		if(O_layerID) delete O_layerID; O_layerID = new std::vector<int>;
		N_trigger++;

		// get basical cdc hit information
		nHitsgood = 0;
		v_hit_index.clear();
		v_hit_dd.clear();
		v_hit_flag.clear();
		int n1 = 0;
		for (int j = 0; j<NLAY; j++){
			a_nHitsLayer[j] = 0;
		}
		for (int ihit = 0; ihit<i_type->size(); ihit++){
			// FIXME: only use center cells
			//if ((*i_layerID)[ihit]<=3){
			//	if ((*i_wireID)[ihit]!=3){
			//		continue;
			//	}
			//}
			//else {
			//	if ((*i_wireID)[ihit]!=4){
			//		continue;
			//	}
			//}
			if ((*i_type)[ihit]!=0&&(*i_type)[ihit]!=1) continue; // including guard layer but without dummy layer
			double dt = (*i_driftT)[ihit]+t0shift;
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			int status;
			double dd = t2x(dt,lid,wid,0,status);
			// FIXME
			if (lid!=testlayer&&!status) continue; // not the test layer and beyond the known range.
			a_nHitsLayer[lid]++;
			v_hit_index.push_back(ihit);
			v_hit_dd.push_back(dd);
//			if (lid!=testlayer&&dd>1&&a_nHitsLayer[lid]==1){ // driftD larger than 1 mm and no good hit in this layer yet.
			if (lid!=testlayer&&dd>1){ // driftD larger than 1 mm
				v_hit_flag[ihit]=0; // golden hits
				nHitsgood++;
			}
            else{
                v_hit_flag.push_back(-2); // default: bad hit
            }
		}
//		std::cout<<"nHitsgood = "<<nHitsgood<<std::endl;
		if (nHitsgood<3) continue;
		N_found++;
//		std::cout<<"Found!"<<std::endl;

		if (checkLR()) continue;
		N_good++;
//		std::cout<<"Good!"<<std::endl;

		for(int index = 0; index<v_hit_index.size(); index++){
			int ihit = v_hit_index[index];
			if (v_hit_flag[index]>=-1&&v_hit_flag[index]<=1||(*i_layerID)[ihit]==testlayer){ // only keep test layer hits and good hits in other layers
				double dt = (*i_driftT)[ihit]+t0shift;
				int lid = (*i_layerID)[ihit];
				int wid = (*i_wireID)[ihit];
				O_wireID->push_back(wid);
				O_layerID->push_back(lid);
				O_driftT->push_back(dt);
				int status;
				O_driftD->push_back(t2x(dt,lid,wid,v_hit_flag[index],status));
				O_lr->push_back(v_hit_flag[index]);
			}
		}

		fitxy();
		O_tx1 = f_x->Eval(yup);
		O_tx2 = f_x->Eval(ydown);
		O_tz1 = 0;
		O_tz2 = 0;

//		for(int index = 0; index<v_hit_index.size(); index++){
//		}

		ot->Fill();
	}// end of event loop

	ot->Write();
	of->Close();
	printf("Triggered Events: %d\n",N_trigger);
	printf("Found Events: %d\n",N_found);
	printf("Good Events: %d\n",N_good);
	return 0;
}
