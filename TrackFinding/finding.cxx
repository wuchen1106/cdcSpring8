#include <vector>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "TRandom.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TVector3.h"

#define NLAY 8
#define NCEL 11
#define NCHT 96
#define NCHS 48

#define use_extra

//===================About tracking============================
std::vector<double> vec_i; // k VS i: z/Lhalf VS ilayer
std::vector<double> vec_k; // k VS i: z/Lhalf VS ilayer
double ar_x[50]; // x VS y
double ar_y[50]; // x VS y
double aar_x[3][50]; // x VS y: all 3 possibilities
double aar_y[3][50]; // x VS y: all 3 possibilities
int n_xy = 0; // number of paris of adjacent hits; to make x-y point

TF1 * f_k = new TF1("f_k","pol1",-0.5,6.5); // k VS i: z/Lhalf VS ilayer
TGraph * g_k = 0; // k VS i: z/Lhalf VS ilayer
TF1 * f_x = new TF1("f_x","pol1",50,64); // x VS y
TGraph * g_x = 0; // x VS y

TF1 * f_xt_l = new TF1("f_xt_l","pol5",-10,250);
TF1 * f_xt_r = new TF1("f_xt_r","pol5",-10,250);

//===================About geometry============================
double map_wid[NCHT];
double map_lid[NCHT];
double map_xhv[NLAY][NCEL];
double map_yhv[NLAY][NCEL];
double map_xro[NLAY][NCEL];
double map_yro[NLAY][NCEL];
double map_k[NLAY][NCEL][NCEL];
int wmax[NLAY]; // maximum wireID
double U = 0.8; // cm

//===================About chamber============================
double Lchamebr = 59.917; // cm
double yup = 53; // cm; layer 1
double ydown = 62.6; // cm; layer 7
// FIXME: currently only use costantant drift velocity. Will implement an interface to load xt cureves later.
double driftTMax = 220; // ns

//===================About hits============================
int i_layerhit[NLAY][NCEL];
int c_layerhit[NLAY];

double t2x(double dt){
	double dd = 0;
	//dd = dt*U/driftTMax;
	return dd;
}

void print_usage(char* prog_name);
double getk(int lid, int wid, double dd){
	double xhv = map_xhv[lid][wid];
	double yhv = map_yhv[lid][wid];
	double xro = map_xro[lid][wid];
	double yro = map_yro[lid][wid];
	double y1 = ydown;
	double y2 = yup;
	double k1 = f_k->Eval(0);
	double k2 = f_k->Eval(6);
	double x1 = f_x->Eval(y1);
	double x2 = f_x->Eval(y2);
	double z1 = k1*Lchamebr/2;
	double z2 = k2*Lchamebr/2;
	TVector3 vt1t2(x2-x1,y2-y1,z2-z1);
	TVector3 vw1w2(xro-xhv,yro-yhv,Lchamebr);
	vt1t2 = vt1t2.Unit();
	vw1w2 = vw1w2.Unit();
	TVector3 vt1w2(xro-x1,yro-y1,Lchamebr/2-z1);
	TVector3 vdist = (vt1t2.Cross(vw1w2)).Unit();
	double d0 = fabs(vt1w2*vdist);
	double sina = vdist.y()/vdist.Perp();
	double cosa = vdist.x()/vdist.Perp();
	double xpro = xro*cosa+yro*sina;
	double xphv = xhv*cosa+yhv*sina;
	double xp1 = x1*cosa+y1*sina;
	double xp2 = x2*cosa+y2*sina;
	double ul = (xpro+xphv)/2.;
	double ur = (-k1*xp2+k2*xp1)/(k2-k1);
	double dl = (xpro-xphv)/2.;
	double dr = (xp2-xp1)/(k2-k1);
	double k = (-ul+ur)/(dl-dr);
	double z = k*Lchamebr/2;
	double deltaz = dd>d0?sqrt(dd*dd-d0*d0)/sina:0;
	double zini = f_k->Eval(lid)*Lchamebr/2;
//	printf("===>input:\n");
//	printf("    %lf,%lf -- %lf,%lf\n",xphv,-1.,xpro,1.);
//	printf("    %lf,%lf -- %lf,%lf\n",xp1,k1,xp2,k2);
//	printf("    theta = (%lf,%lf)\n",sina,cosa);
//	printf("    dd = %lf, d0 = %lf,vdist(%lf,%lf,%lf),vt1w2(%lf,%lf,%lf)\n",dd,d0,vdist.x(),vdist.y(),vdist.z(),vt1w2.x(),vt1w2.y(),vt1w2.z());
//	printf("  =>output:\n");
//	printf("    k = (-%lf+%lf)/(%lf-%lf) = %lf\n",ul,ur,dl,dr,k);
//	printf("    z = %lf+-%lf\n",z,deltaz);
	double zfinal = fabs(z+deltaz-zini)<fabs(z-deltaz-zini)?z+deltaz:z-deltaz;
	if (z2>z1){
		if (zfinal>z2) zfinal = z2;
		else if (zfinal<z1) zfinal = z1;
	}
	else{
		if (zfinal>z1) zfinal = z1;
		else if (zfinal<z2) zfinal = z2;
	}
	return (zfinal);
}

void getxy(double x1,double y1,double r1,double x2,double y2,double r2,double &ox1, double &ox2, double &ox3, double &oy1, double &oy2, double &oy3){
	double dist = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
	ox1 = (r1*x2+r2*x1)/(r1+r2);
	oy1 = (r1*y2+r2*y1)/(r1+r2);
	double ox0 = (x1+x2)/2.;
	double oy0 = (y1+y2)/2.;
	double sina = (y2-y1)/dist/dist;
	double cosa = (x2-x1)/dist/dist;
	double sinb = (r1-r2)/dist;
	double cosb = sqrt(1-sinb*sinb);
	double rbar = (r1+r2)/2.;
	ox2 = ox0+rbar*(sina*cosb-cosa*sinb);
	oy2 = oy0-rbar*(cosa*cosb+sina*sinb);
	ox3 = ox0-rbar*(sina*cosb+cosa*sinb);
	oy3 = oy0+rbar*(cosa*cosb-sina*sinb);
}

void fityz(std::vector<int> * wireID, std::vector<double> * driftD){
	vec_k.clear();
	vec_i.clear();
	for (int j = 0; j<7; j++){
		for (int hid = 0; hid<c_layerhit[j]; hid++){
			int index = i_layerhit[j][hid];
			vec_i.push_back(j);
			vec_k.push_back(getk(j,(*wireID)[index],(*driftD)[index]));
		}
	}
	if (g_k) delete g_k; g_k = new TGraph(vec_k.size(),&(vec_i[0]),&(vec_k[0]));
	g_k->Fit("f_k","qN0","");
}

int fityz(std::vector<int> * wireID){
	vec_k.clear();
	vec_i.clear();
	int wid1,wid2;
	double k;
	int nlayers = 0;
	for (int j = 1; j<7; j++){
		bool check = false;
		for (int hid1 = 0; hid1<c_layerhit[j]; hid1++){
			int index1 = i_layerhit[j][hid1];
			wid1 = (*wireID)[index1];
			for (int hid2 = 0; hid2<c_layerhit[j+1]; hid2++){
				int index2 = i_layerhit[j+1][hid2];
				wid2 = (*wireID)[index2];
				k = map_k[j][wid1][wid2];
				if (fabs(k)<1){
					vec_i.push_back(j+0.5);
					vec_k.push_back(k);
					if (!check) check = true;
				}
			}
		}
		if (check) nlayers++;
	}
	if (nlayers<3) return 1;
	if (g_k) delete g_k; g_k = new TGraph(vec_k.size(),&(vec_i[0]),&(vec_k[0]));
	g_k->Fit("f_k","qN0","");
	return 0;
	//g_k->Print();
	//f_k->Print();
}

int fitxy(std::vector<int> * wireID,std::vector<double> * driftD){
	int wid1,wid2;
	n_xy = 0;
	for (int j = 0; j<6; j++){
		for (int hid1 = 0; hid1<c_layerhit[j]; hid1++){
			int index1 = i_layerhit[j][hid1];
			wid1 = (*wireID)[index1];
			double k1 = f_k->Eval(j);
			double x1 = ((1+k1)*map_xro[j][wid1]+(1-k1)*map_xhv[j][wid1])/2.;
			double y1 = ((1+k1)*map_yro[j][wid1]+(1-k1)*map_yhv[j][wid1])/2.;
			double r1 = (*driftD)[index1];
			for (int hid2 = 0; hid2<c_layerhit[j+1]; hid2++){
				int index2 = i_layerhit[j+1][hid2];
				wid2 = (*wireID)[index2];
				double k2 = f_k->Eval(j+1);
				double x2 = ((1+k2)*map_xro[j+1][wid2]+(1-k2)*map_xhv[j+1][wid2])/2.;
				double y2 = ((1+k2)*map_yro[j+1][wid2]+(1-k2)*map_yhv[j+1][wid2])/2.;
				double r2 = (*driftD)[index2];
				getxy(x1,y1,r1,x2,y2,r2,aar_x[0][n_xy],aar_x[1][n_xy],aar_x[2][n_xy],aar_y[0][n_xy],aar_y[1][n_xy],aar_y[2][n_xy]);
				n_xy++;
			}
		}
	}
	if (n_xy==0) return 1;

	//double chi2min = 1e9;
	//int index = -1;
	//for (int i = 0; i<pow(3,n_xy); i++){
	//	int a = i;
	//	for (int j = 0; j<n_xy; j++){
	//		int pow3j = pow(3,n_xy-j-1);
	//		int b = a/pow3j;
	//		a-=pow3j*b;
	//		ar_x[j] = aar_x[b][j];
	//		ar_y[j] = aar_y[b][j];
	//	}
	//	if (g_x) delete g_x; g_x = new TGraph(n_xy,&(ar_x[0]),&(ar_y[0]));
	//	g_x->Fit("f_x","qN0","");
	//	double chi2 = f_x->GetChisquare();
	//	if (chi2<chi2min){
	//		chi2min=chi2;
	//		index = i;
	//	}
	//}
	//int a = index;
	//for (int j = 0; j<n_xy; j++){
	//	int pow3j = pow(3,n_xy-j-1);
	//	int b = a/pow3j;
	//	a-=pow3j*b;
	//	ar_x[j] = aar_x[b][j];
	//	ar_y[j] = aar_y[b][j];
	//}
	//if (g_x) delete g_x; g_x = new TGraph(n_xy,&(ar_x[0]),&(ar_y[0]));
	if (g_x) delete g_x; g_x = new TGraph(n_xy,&(aar_y[0][0]),&(aar_x[0][0]));
	g_x->Fit("f_x","qN0","");
	return 0;
	//g_x->Print();
	//f_x->Print();
}

int main(int argc, char** argv){

	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	if (argc>=3) nEventMax = (int)strtol(argv[2],NULL,10);
	std::string suffix = "";
	if (argc>=4){
		suffix  = argv[3];
		suffix=suffix+".";
	}

	//===================Get wire position============================
	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.v3.root");
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
	int map_wid[NCHT];
	int map_lid[NCHT];
	for ( int ch = 0; ch<NCHT; ch++ ){
		map_lid[ch] = -1;
		map_wid[ch] = -1;
	}
	for (int i = 0; i<NLAY; i++){
		wmax[i] = -1;
	}

	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<=7){
			map_xhv[wp_lid][wp_wid] = wp_xhv/10.;
			map_yhv[wp_lid][wp_wid] = wp_yhv/10.;
			map_xro[wp_lid][wp_wid] = wp_xro/10.;
			map_yro[wp_lid][wp_wid] = wp_yro/10.;
			map_wid[wp_ch+NCHS*wp_bid] = wp_wid;
			map_lid[wp_ch+NCHS*wp_bid] = wp_lid;
			if (wmax[wp_lid]<wp_wid) wmax[wp_lid] = wp_wid;
		}
	}
	for (int k = 1; k<NLAY-1; k++){
		printf("layer %d & %d\n",k,k+1);
		for (int i = 0; i<=wmax[k]; i++){
			for (int j = 0; j<=wmax[k+1]; j++){
				if ((-(map_xro[k+1][j]-map_xhv[k+1][j])+(map_xro[k][i]-map_xhv[k][i]))){
					map_k[k][i][j] = ((map_xro[k+1][j]+map_xhv[k+1][j])-(map_xro[k][i]+map_xhv[k][i]))/(-(map_xro[k+1][j]-map_xhv[k+1][j])+(map_xro[k][i]-map_xhv[k][i]));
				}
				else{
					map_k[k][i][j] = 10;
				}
				printf("  [%d,%d]: %lf\n",i,j,map_k[k][i][j]);
			}
		}
	}

	//===================Get input ROOT file============================
	TChain * c = new TChain("t","t");
	char fileName[128];
	sprintf(fileName,"../root/d_%d.root",runNo);
	c->Add(fileName);
	int triggerNumber;
	int d_i[NCHT];
	double d_dt[NCHT];
	int d_h[NCHT];
	double d_aa[NCHT];
	c->SetBranchAddress("triggerNumber",&triggerNumber);
	c->SetBranchAddress("i",d_i);
	c->SetBranchAddress("dt",d_dt);
	c->SetBranchAddress("h",d_h);
	c->SetBranchAddress("aa",d_aa);

	//===================Prepare output ROOT file============================
	std::vector<int> * O_wireID = 0;
	std::vector<int> * O_layerID = 0;
	std::vector<double> * O_driftD = 0;
	std::vector<double> * O_driftT = 0;
	std::vector<double> * O_sum = 0;
	std::vector<int> * O_peak = 0;
	double O_tx1;
	double O_tx2;
	double O_tz1;
	double O_tz2;
	int nHits;

	sprintf(fileName,"../root/h_%d.root",runNo);
	TFile * of = new TFile(fileName,"RECREATE");
	TTree * ot = new TTree("t","t");

	ot->Branch("nHits",&nHits);
	ot->Branch("triggerNumber",&triggerNumber);
	ot->Branch("tx1",&O_tx1);
	ot->Branch("tx2",&O_tx2);
	ot->Branch("tz1",&O_tz1);
	ot->Branch("tz2",&O_tz2);
	ot->Branch("driftT",&O_driftT);
	ot->Branch("driftD",&O_driftD);
	ot->Branch("wireID",&O_wireID);
	ot->Branch("layerID",&O_layerID);
	ot->Branch("sum",&O_sum);
	ot->Branch("peak",&O_peak);

	//===================Efficiency Counter============================
	int N_layerhit[NLAY];
	for (int i = 0; i<NLAY; i++){
		N_layerhit[i] = 0;
	}
	int N_allLayerhit = 0;
	int N_trigger = 0;
	int N_hit = 0;
	int N_found = 0;

	//===================Track Finding============================
	for ( int i = 0; i<c->GetEntries(); i++){
		if (i%1000==0) printf("%lf%...\n",(double)i/c->GetEntries()*100);
		// FIXME
		//printf("*****Event %d*******\n",i);
		c->GetEntry(i);
		if(O_driftT) delete O_driftT; O_driftT = new std::vector<double>;
		if(O_driftD) delete O_driftD; O_driftD = new std::vector<double>;
		if(O_wireID) delete O_wireID; O_wireID = new std::vector<int>;
		if(O_layerID) delete O_layerID; O_layerID = new std::vector<int>;
		if(O_sum) delete O_sum; O_sum = new std::vector<double>;
		if(O_peak) delete O_peak; O_peak = new std::vector<int>;

		// get basical cdc hit information
		nHits= 0;
		int nGoodHits = 0;
		for (int j = 0; j<NLAY; j++){
			c_layerhit[j] = 0;
		}
		for (int ch = 0; ch<NCHT; ch++){
			if (d_i[ch]<0) continue;
			int lid = map_lid[ch];
			int wid = map_wid[ch];
			double dt = d_dt[ch];
			//if (wid==0||wid==wmax[lid]) continue;
			if (lid<0) continue;
			if (dt<-0||dt>driftTMax) continue;
			i_layerhit[lid][c_layerhit[lid]]=O_wireID->size();
			c_layerhit[lid]++;
			O_wireID->push_back(wid);
			O_layerID->push_back(lid);
			O_driftD->push_back(t2x(dt));
			O_driftT->push_back(dt);
			O_peak->push_back(d_h[ch]);
			O_sum->push_back(d_aa[ch]);
			nHits++;
			if (lid!=4) nGoodHits++;
		}
		N_trigger++;

		int c_allLayerhit=0;
		for (int j = 1; j<NLAY; j++){
			//std::cout<<"  ["<<j<<"]: ";
			if (c_layerhit[j]){
				N_layerhit[j]++;
				if (j!=4)
					c_allLayerhit++;
				//std::cout<<"hit"<<std::endl;
			}
			else{
				//std::cout<<"no hit"<<std::endl;
			}
		}
		if (c_allLayerhit==6) N_allLayerhit++;
		//std::cout<<"c_allLayerhit = "<<c_allLayerhit<<std::endl;

		if (c_allLayerhit!=6) continue;
		if (nHits<6||nHits>7) continue;
		//if (nHits!=8) continue;
		//if (!c_layerhit[4]) continue;
		N_hit++;

		// fit y-z plane
		if (fityz(O_wireID)) continue;
		//printf("After y-z fit 1:\n");
		//g_k->Print();
		//f_k->Print();
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int wid = (*O_wireID)[index];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = ((1+k_temp)*xro+(1-k_temp)*xhv)/2.;
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*Lchamebr/2);
		//	}
		//}

		// fit x-y plane
		if (fitxy(O_wireID,O_driftD)) continue;
		//printf("After x-y fit 1:\n");
		//g_x->Print();
		//f_x->Print();
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int wid = (*O_wireID)[index];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = f_x->Eval(y_temp);
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*Lchamebr/2);
		//	}
		//}

		// correct z positions and fit y-z again
		//fityz(O_wireID,O_driftD);
		//printf("After y-z fit 2:\n");
		//g_k->Print();
		//f_k->Print();
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int wid = (*O_wireID)[index];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = f_x->Eval(y_temp);
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*Lchamebr/2);
		//	}
		//}

		// fit x-y again
		//fitxy(O_wireID,O_driftD);
		//printf("After x-y fit 2:\n");
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int wid = (*O_wireID)[index];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = f_x->Eval(y_temp);
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*Lchamebr/2);
		//	}
		//}

		N_found++;
		O_tx1 = f_x->Eval(ydown);
		O_tx2 = f_x->Eval(yup);
		O_tz1 = f_k->Eval(0)*Lchamebr/2;
		O_tz2 = f_k->Eval(6)*Lchamebr/2;
		ot->Fill();
	}// end of event loop

	ot->Write();
	of->Close();
	printf("Triggered Events: %d\n",N_trigger);
	printf("layer 1 hit Events: %d\n",N_layerhit[1]);
	printf("layer 2 hit Events: %d\n",N_layerhit[2]);
	printf("layer 3 hit Events: %d\n",N_layerhit[3]);
	printf("layer 4 hit Events: %d\n",N_layerhit[4]);
	printf("layer 5 hit Events: %d\n",N_layerhit[5]);
	printf("layer 6 hit Events: %d\n",N_layerhit[6]);
	printf("layer 7 hit Events: %d\n",N_layerhit[7]);
	printf("All layer hit Events: %d\n",N_allLayerhit);
	printf("7-9 hits Events: %d\n",N_hit);
	printf("Found Events: %d\n",N_found);
	return 0;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
