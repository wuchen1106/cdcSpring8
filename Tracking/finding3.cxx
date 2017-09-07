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

int pow2i[NLAY];

//===================Functions============================
double t2x(double time, int lid, int wid, int lr, int & status);
void print_usage(char* prog_name);
int checkHits();
int fityx();
int fityz();

//===================About tracking============================
std::vector<double> * i_driftT = 0;
std::vector<int> * i_layerID = 0;
std::vector<int> * i_wireID = 0;
std::vector<int> * i_type = 0;
double a_yx_x[50]; // x VS y
double a_yx_y[50]; // x VS y
int np_yx;
double a_yz_z[50]; // z VS y
double a_yz_y[50]; // z VS y
int np_yz;
double yup = 623.97007;
double ydown = 527.60011;
TF1 * f_x = new TF1("f_x","pol1",500,640); // x VS y
TGraph * g_x = 0; // x VS y
TF1 * f_z = new TF1("f_z","pol1",500,640); // z VS y
TGraph * g_z = 0; // x VS y

int testlayer = 0;
int nHitsgood;
std::vector<int> v_hit_index;
std::vector<int> v_hit_flag;
std::vector<double> v_hit_dd;
std::vector<int> v_pair_lid;
std::vector<int> v_pair_flag;
std::vector<int> v_pair_iup;
std::vector<int> v_pair_idown;
int nHits_ldd[NLAY];

//===================About chamber============================
double U = 8; // mm
double Lchamebr = 599.17; // mm
double zhv = -Lchamebr/2;
double zro = Lchamebr/2;

double map_xc[NLAY][NCEL];
double map_yc[NLAY][NCEL];
double map_xhv[NLAY][NCEL];
double map_yhv[NLAY][NCEL];
double map_xro[NLAY][NCEL];
double map_yro[NLAY][NCEL];
double map_k[NLAY][NCEL][NCEL];
int map_check[NLAY][NCEL];

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

int main(int argc, char** argv){

	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	std::string suffix = "";
	if (argc>=3){
		suffix  = argv[2];
		suffix="."+suffix;
	}
	if (argc>=4) testlayer = (int)strtol(argv[3],NULL,10);
	int t0shift = 0;
	if (argc>=5) t0shift = (int)strtol(argv[4],NULL,10);
	printf("t0shift = %d\n",t0shift);
	int nEventMax = 0;
	if (argc>=6) nEventMax = (int)strtol(argv[5],NULL,10);

	for (int i=0; i<NLAY; i++){
		pow2i[i] = pow(2,i);
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
	for(int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int iwire = 0; iwire<NCEL; iwire++){
			map_check[ilayer][iwire]=0;
		}
	}
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<NLAY){
			map_xhv[wp_lid][wp_wid] = wp_xhv;
			map_yhv[wp_lid][wp_wid] = wp_yhv;
			map_xro[wp_lid][wp_wid] = wp_xro;
			map_yro[wp_lid][wp_wid] = wp_yro;
			map_xc[wp_lid][wp_wid] = (wp_xhv+wp_xro)/2.;
			map_yc[wp_lid][wp_wid] = (wp_yhv+wp_yro)/2.;
			map_check[wp_lid][wp_wid] = 1;
		}
	}
	for(int ilayer = 0; ilayer<NLAY-1; ilayer++){
		printf("layer %d & %d\n",ilayer,ilayer+1);
		for (int iwire = 0; iwire<NCEL; iwire++){
			if (!map_check[ilayer][iwire]) continue;
			for (int jwire = 0; jwire<NCEL; jwire++){
				if (!map_check[ilayer+1][jwire]) continue;
				if ((-(map_xro[ilayer+1][jwire]-map_xhv[ilayer+1][jwire])+(map_xro[ilayer][iwire]-map_xhv[ilayer][iwire]))){
					map_k[ilayer][iwire][jwire] = ((map_xro[ilayer+1][jwire]+map_xhv[ilayer+1][jwire])-(map_xro[ilayer][iwire]+map_xhv[ilayer][iwire]))/(-(map_xro[ilayer+1][jwire]-map_xhv[ilayer+1][jwire])+(map_xro[ilayer][iwire]-map_xhv[ilayer][iwire]));
				}
				else{
					map_k[ilayer][iwire][jwire] = 10;
				}
				printf("  [%d,%d]: %lf\n",iwire,jwire,map_k[ilayer][iwire][jwire]*Lchamebr/2);
			}
		}
	}

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
	return 0;

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
	std::vector<int> * O_wireID = 0;
	std::vector<int> * O_layerID = 0;
	std::vector<int> * O_lr = 0;
	std::vector<double> * O_driftD = 0;

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

		//========================================================================================================
		// get basical cdc hit information
		nHitsgood = 0;
		v_hit_index.clear();
		v_hit_dd.clear();
		v_hit_flag.clear();
		int n1 = 0;
		for (int j = 0; j<NLAY; j++){
			nHits_ldd[j] = 0;
		}
		for (int ihit = 0; ihit<i_type->size(); ihit++){
			if ((*i_type)[ihit]!=0&&(*i_type)[ihit]!=1) continue; // including guard layer but without dummy layer
			double dt = (*i_driftT)[ihit]+t0shift;
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			int status;
			double dd = t2x(dt,lid,wid,0,status);
			// FIXME: maybe we want to put more strict limitations on this
			if (lid!=testlayer&&!status) continue; // beyond the known range.
			//if (lid!=testlayer&&abs(status)!=1) continue; // beyond the known range; without corner hits (tail)
			v_hit_index.push_back(ihit);
			v_hit_dd.push_back(dd);
			if (dd>1){
				nHits_ldd[lid]++;
				v_hit_flag.push_back(1); // default: right hit
			}
			else{
				v_hit_flag.push_back(0); // too small to tell left or right
			}
		}

		//========================================================================================================
		// To find pair candidates: two adjacent golden layers
		// golden layer: with only 1 hit and dd>1mm
		for(int ilayer = 0; ilayer<NLAY-1; ilayer++){
			if (ilayer==testlayer||ilayer+1==testlayer) continue;
			if (nHits_ldd[ilayer]==1&&nHits_ldd[ilayer+1]==1){// a new pair candidate
				v_pair_lid.push_back(ilayer);
				v_pair_flag.push_back(0); // default: no problem
				v_pair_idown.push_back(-1);
				v_pair_iup.push_back(-1);
			}
		}
//		std::cout<<"nPairs = "<<v_pair_lid.size()<<std::endl;
		if (v_pair_lid.size()<2) continue; // Need at least two points to determine the track on y-z projection
		N_found++;
		for(int ipair = 0; ipair<v_pair_lid.size(); ipair++){
			int thelayer=v_pair_lid[ipair];
			for(int ihit = 0; ihit<v_hit_index.size(); ihit++){
				if (!v_hit_flag[ihit]) continue; // golden hits should be with large drift distance (>1mm)
				if (lid==thelayer){
					v_pair_idown[ipair]=ihit;
				}
				else if (lid==thelayer+1){
					v_pair_iup[ipair]=ihit;
				}
			}
		}

		//========================================================================================================
		// To fit y-z projection
		if(fityz()) continue;
//		std::cout<<"y-z fitting succeeded!"<<std::endl;

		//========================================================================================================
		// To fit y-x projection
		if(fityx()) continue;
//		std::cout<<"y-x fitting succeeded!"<<std::endl;

		//========================================================================================================
		// To check hits
		if (checkHits()<5) continue;
		N_good++;
//		std::cout<<"Found enough hit candidates (>5)!"<<std::endl;

		O_tx1 = f_x->Eval(yup);
		O_tx2 = f_x->Eval(ydown);
		O_tz1 = f_z->Eval(yup);
		O_tz2 = f_z->Eval(ydown);

		for(int index = 0; index<v_hit_index.size(); index++){
			int ihit = v_hit_index[index];
			if (v_hit_flag[index]>=-1&&v_hit_flag[index]<=1){ // no bad hits
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

		ot->Fill();
	}// end of event loop

	ot->Write();
	of->Close();
	printf("Triggered Events: %d\n",N_trigger);
	printf("Found Events: %d\n",N_found);
	printf("Good Events: %d\n",N_good);
	return 0;
}

int checkHits(){
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

int fityz(){
	// loop in different left/right combinations
	for (int icombi = 0; icombi<ncombi; icombi++){
		// get y z arrays
		for (int ipair = 0; ipair<v_pair_iup.size(); ipair++){
			int lid = v_pair_lid[ipair];
			int ihit_up = v_pair_iup[ipair];
			int wid_up = (*i_wireID)[ihit_up];
			int ihit_down = v_pair_idown[ipair];
			int wid_down = (*i_wireID)[ihit_down];
			double k = map_k[lid][wid_down][wid_up];
		}
		// fit f_z
		// delete some points
		// fit f_z again
	}
}

int fityx(){
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

//______________________________________________________________________________
void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
