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
double driftTMax = 0; // ns
double driftTMin = 5; // ns
double driftgoodTMin = 20; // ns
std::vector<double> v_x_all_right;
std::vector<double> v_t_all_right;
std::vector<double> v_x_all_left;
std::vector<double> v_t_all_left;
TGraph * g_xt_all_right = 0;
TGraph * g_xt_all_left = 0;
TF1 * f_xt_all_right = 0;
TF1 * f_xt_all_left = 0;

//===================About hits============================
int c_layerhit[NLAY];

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
		int usedhits = 0;
		double chi2 = 0;
		double xc = 0;
		double xp = 0;
		double sumx = 0;
//		std::cout<<"#################################################"<<std::endl;
		for (int index= 0; index<v_hit_index.size(); index++){
			if (v_hit_flag[index]) continue;
			int ihit = v_hit_index[index];
			int lr0 = (icombi/pow2i[usedhits])%2;

			double dd0 = v_hit_dd[index];
			xc = map_xc[(*i_layerID)[ihit]][(*i_wireID)[ihit]];
			if (lr0) xc+=dd0;
			else xc-=dd0;

			if (usedhits!=0)
				chi2+= pow(xc-xp,2);

			xp = xc;
			sumx+=xc;

			usedhits++;
		}
//		std::cout<<"chi2 = "<<chi2<<std::endl;
		if (chi2<chi2min){
			chi2min=chi2;
			thecombi=icombi;
			meanxmin=sumx/nHitsgood;
		}
	}
	nHitsgood=0;
//	std::cout<<"meanxmin = "<<meanxmin<<std::endl;
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
//		if (v_hit_flag[index] == -2) std::cout<<"  ! x["<<(*i_layerID)[ihit]<<"]["<<(*i_wireID)[ihit]<<"] = "<<xc<<"+-"<<v_hit_dd[index]<<" = "<<xl<<" or "<<xr<<std::endl;
//		else if (v_hit_flag[index] == -1) std::cout<<"  x["<<(*i_layerID)[ihit]<<"]["<<(*i_wireID)[ihit]<<"] = "<<xc<<"-"<<v_hit_dd[index]<<" = "<<xl<<std::endl;
//		else if (v_hit_flag[index] == 1) std::cout<<"  x["<<(*i_layerID)[ihit]<<"]["<<(*i_wireID)[ihit]<<"] = "<<xc<<"+"<<v_hit_dd[index]<<" = "<<xr<<std::endl;
//		else if (v_hit_flag[index] == 0) std::cout<<"  x["<<(*i_layerID)[ihit]<<"]["<<(*i_wireID)[ihit]<<"] = "<<xc<<"+-"<<v_hit_dd[index]<<" = "<<xl<<" or "<<xr<<std::endl;
	}
	if (nHitsgood<5) return 1;
//	std::cout<<"chi2min = "<<chi2min<<std::endl;
	return 0;
}

int fitxy(){
	for (int ihit = 0; ihit<O_lr->size(); ihit++){
		ar_x[ihit] = map_xc[(*O_layerID)[ihit]][(*O_wireID)[ihit]]+(*O_driftD)[ihit];
		ar_y[ihit] = map_yc[(*O_layerID)[ihit]][(*O_wireID)[ihit]];
	}
	if (g_x) delete g_x; g_x = new TGraph(O_lr->size(),&(ar_y[0]),&(ar_x[0]));
	g_x->Fit("f_x","qN0","");
	return 0;
}

double t2x(double time, int lr){
	TF1 * f = 0;
	if (lr==-1) f = f_xt_all_left;
	else f = f_xt_all_right;
	double dd = f->Eval(time);
	if ((lr==1&&dd<0)||(lr==-1&&dd>0)) dd = 0;
	if (time>driftTMax) dd = lr==1?U:-U;
	return dd;
	// FIXME: may want to support multi-xt
//	double dd;
//	int index = -1;
//	std::vector<double> * vt;
//	std::vector<double> * vx;
//	if (sign>0){
//		vt = &v_t_all_right;
//		vx = &v_x_all_right;
//	}
//	else {
//		vt = &v_t_all_left;
//		vx = &v_x_all_left;
//	}
//	for (index = 0; index<vt->size(); index++){
//		if ((*vt)[index]>time) break;
//	}
//	if (index==vt->size()-1) dd = 0.8;
//	else if (index==0) dd = 0;
//	else dd = ((*vx)[index-1]*((*vt)[index]-time)+(*vx)[index]*(time-(*vt)[index-1]))/((*vt)[index]-(*vt)[index-1]);
//	return dd;
}

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
	TChain * c_xt = new TChain("t","t");
	c_xt->Add(Form("../info/xt.%d.root",runNo));
	double i_xt_x, i_xt_t;
	c_xt->SetBranchAddress("x",&i_xt_x);
	c_xt->SetBranchAddress("t",&i_xt_t);
	for ( int i = 0; i<c_xt->GetEntries(); i++ ){
		c_xt->GetEntry(i);
		if (i_xt_x>=-0.02){
			v_x_all_right.push_back(i_xt_x);
			v_t_all_right.push_back(i_xt_t);
		}
		if (i_xt_x<=0.02){
			v_x_all_left.push_back(i_xt_x);
			v_t_all_left.push_back(i_xt_t);
		}
		if (i_xt_t>driftTMax) driftTMax=i_xt_t;
	}
	double temp;
	for (int i = 0; i<v_x_all_left.size(); i++){
		for (int j = i+1; j<v_x_all_left.size(); j++){
			if (v_t_all_left[j]<v_t_all_left[i]){
				temp = v_t_all_left[i];
				v_t_all_left[i]=v_t_all_left[j];
				v_t_all_left[j]=temp;
				temp = v_x_all_left[i];
				v_x_all_left[i]=v_x_all_left[j];
				v_x_all_left[j]=temp;
			}
		}
	}
	for (int i = 0; i<v_x_all_right.size(); i++){
		for (int j = i+1; j<v_x_all_right.size(); j++){
			if (v_t_all_right[j]<v_t_all_right[i]){
				temp = v_t_all_right[i];
				v_t_all_right[i]=v_t_all_right[j];
				v_t_all_right[j]=temp;
				temp = v_x_all_right[i];
				v_x_all_right[i]=v_x_all_right[j];
				v_x_all_right[j]=temp;
			}
		}
	}
	f_xt_all_right = new TF1("f_xt_all_right","pol5",driftTMin-tres,driftTMax);
	f_xt_all_left = new TF1("f_xt_all_left","pol5",driftTMin-tres,driftTMax);
	g_xt_all_right = new TGraph(v_x_all_right.size(),&(v_t_all_right[0]),&(v_x_all_right[0]));
	g_xt_all_left = new TGraph(v_x_all_left.size(),&(v_t_all_left[0]),&(v_x_all_left[0]));
	g_xt_all_right->Fit("f_xt_all_right","qN0","");
	g_xt_all_left->Fit("f_xt_all_left","qN0","");
	TCanvas * c1 = new TCanvas("c1","c1");
	g_xt_all_left->Draw("ALP");
	f_xt_all_left->Draw("SAME");
	c1->SaveAs("left.pdf");
	g_xt_all_right->Draw("ALP");
	f_xt_all_right->Draw("SAME");
	c1->SaveAs("right.pdf");

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
	int nHits;

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
		nHits= 0;
		nHitsgood = 0;
		v_hit_index.clear();
		v_hit_dd.clear();
		v_hit_flag.clear();
		int n1 = 0;
		for (int j = 0; j<NLAY; j++){
			c_layerhit[j] = 0;
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
			if (lid!=testlayer&&(dt<driftTMin-tres||dt>driftTMax+tres)) continue;
			v_hit_index.push_back(ihit);
			v_hit_dd.push_back(t2x(dt,0));
			if (lid==testlayer){
				v_hit_flag.push_back(-10); // test hits
			}
			else{
				c_layerhit[lid]++;
				nHits++;
				if (dt>driftgoodTMin-tres){
					nHitsgood++;
					v_hit_flag.push_back(0); // good hits without left/right sign
				}
				else{
					v_hit_flag.push_back(-2); // default value: bad hits
				}
			}
		}
//		std::cout<<"nHits = "<<nHits<<", nHitsgood = "<<nHitsgood<<std::endl;
		if (nHits>7||nHitsgood<3) continue; // multiple tracks
		for(int ilayer = 1; ilayer<NLAY; ilayer++){
			if (ilayer==testlayer) continue;
			if (c_layerhit[ilayer]==1) n1++;
		}
//		std::cout<<"n1 = "<<n1<<std::endl;
		if (n1<5) continue;
		N_found++;

		if (checkLR()) continue;
		N_good++;
//		std::cout<<"Good!"<<std::endl;

		for(int index = 0; index<v_hit_index.size(); index++){
			int ihit = v_hit_index[index];
			if (v_hit_flag[index]>=-1&&v_hit_flag[index]<=1){
				double dt = (*i_driftT)[ihit]+t0shift;
				int lid = (*i_layerID)[ihit];
				int wid = (*i_wireID)[ihit];
				O_wireID->push_back(wid);
				O_layerID->push_back(lid);
				O_driftT->push_back(dt);
				O_driftD->push_back(t2x(dt,v_hit_flag[index]));
				O_lr->push_back(v_hit_flag[index]);
			}
		}

		fitxy();
		O_tx1 = f_x->Eval(yup);
		O_tx2 = f_x->Eval(ydown);
		O_tz1 = 0;
		O_tz2 = 0;

		ot->Fill();
	}// end of event loop

	ot->Write();
	of->Close();
	printf("Triggered Events: %d\n",N_trigger);
	printf("Found Events: %d\n",N_found);
	printf("Good Events: %d\n",N_good);
	return 0;
}
