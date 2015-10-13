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

#define NLAY 8
#define NCEL 11
#define NCHT 96
#define NCHS 48

#define use_extra

//===================About fitting============================
int thelayer = 0;
TMinuit *gMinuit = 0;
Double_t arglist[10];
Int_t ierflg = 0;
Double_t amin,edm,errdef;
Int_t nvpar,nparx,icstat;
double slerX = 0.42*2;
double inerX = 2*2;

// ________About Cell___________
double U = 0.8; // cm
double zhv = -59.917/2;
double zro = 59.917/2;

//===================About tracking============================
std::vector<int> * O_wireID = 0;
std::vector<int> * O_layerID = 0;
std::vector<double> * O_driftD = 0;
std::vector<double> vec_i; // k VS i: z/Lhalf VS ilayer
std::vector<double> vec_k; // k VS i: z/Lhalf VS ilayer
double ar_x[50]; // x VS y
double ar_y[50]; // x VS y
double aar_x[3][50]; // x VS y: all 3 possibilities
double aar_y[3][50]; // x VS y: all 3 possibilities
int n_xy = 0; // number of paris of adjacent hits; to make x-y point
double yup = 62.397007;
double ydown = 52.760011;
double deltay = yup-ydown;
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

//TF1 * f_k = new TF1("f_k","pol1",-0.5,6.5); // k VS i: z/Lhalf VS ilayer
TF1 * f_k = new TF1("f_k","0",-0.5,6.5); // k VS i: z/Lhalf VS ilayer
TGraph * g_k = 0; // k VS i: z/Lhalf VS ilayer
TF1 * f_x = new TF1("f_x","pol1",50,64); // x VS y
TGraph * g_x = 0; // x VS y

//===================About geometry============================
double map_wid[NCHT];
double map_lid[NCHT];
double map_xhv[NLAY][NCEL];
double map_yhv[NLAY][NCEL];
double map_xro[NLAY][NCEL];
double map_yro[NLAY][NCEL];
double map_k[NLAY][NCEL][NCEL];
double errord[NLAY][NCEL];
int wmax[NLAY]; // maximum wireID

//===================About chamber============================
double Lchamebr = 59.917; // cm
// FIXME: currently only use costantant drift velocity. Will implement an interface to load xt cureves later.
double tres = 10;
double driftTMax = 0; // ns
std::vector<double> v_x_all_right;
std::vector<double> v_t_all_right;
std::vector<double> v_x_all_left;
std::vector<double> v_t_all_left;
TGraph * g_xt_all_right = 0;
TGraph * g_xt_all_left = 0;
TF1 * f_xt_all_right = 0;
TF1 * f_xt_all_left = 0;

//===================About hits============================
int i_layerhit[NLAY][NCEL];
int c_layerhit[NLAY];

//______________________________________________________________________________
Double_t get_dist(int lid, int wid, Double_t slx, Double_t inx)
{
	double xdown = inx-slx*deltay;
	vTrackU.SetXYZ(inx,yup,0);
	vTrackD.SetXYZ(xdown,ydown,0);
	vWireHV.SetXYZ(map_xhv[lid][wid],map_yhv[lid][wid],zhv);
	vWireRO.SetXYZ(map_xro[lid][wid],map_yro[lid][wid],zro);
	vTrack = vTrackD-vTrackU;
	vWire = vWireRO-vWireHV;
	vDist = vWireHV-vTrackU;
	vAxis = vWire.Cross(vTrack);
	double value = vDist*(vAxis.Unit());
	return value;
}

//______________________________________________________________________________
void getchi2(Double_t &f, Double_t slx, Double_t inx,bool all = true)
{
	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	Double_t dfit;
	int N = O_driftD->size();

	for (Int_t i=0;i<N; i++) {
		if ((*O_layerID)[i]==thelayer&&!all) continue;
		dfit = get_dist((*O_layerID)[i],(*O_wireID)[i],slx,inx);
		delta  = (fabs((*O_driftD)[i])-fabs(dfit))/errord[(*O_layerID)[i]][(*O_wireID)[i]];
		chisq += delta*delta;
	}
	f = chisq;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	getchi2(f,*par,*(par+1));
}

void print_usage(char* prog_name);
double getk(int lid, int wid, double dd){
	double xhv = map_xhv[lid][wid];
	double yhv = map_yhv[lid][wid];
	double xro = map_xro[lid][wid];
	double yro = map_yro[lid][wid];
	double y1 = yup;
	double y2 = ydown;
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
	double sina = (y2-y1)/dist;
	double cosa = (x2-x1)/dist;
	double sinb = (r1-r2)/dist;
	double cosb = sqrt(1-sinb*sinb);
	double rbar = (r1+r2)/2.;
	ox2 = ox0+rbar*(sina*cosb-cosa*sinb);
	oy2 = oy0-rbar*(cosa*cosb+sina*sinb);
	ox3 = ox0-rbar*(sina*cosb+cosa*sinb);
	oy3 = oy0+rbar*(cosa*cosb-sina*sinb);
//	printf("(%.3lf,%.3lf;%.3lf)->(%.3lf,%.3lf;%.3lf): (%.3lf,%.3lf)(%.3lf,%.3lf)(%.3lf,%.3lf)\n",x1,y1,r1,x2,y2,r2,ox1,oy1,ox2,oy2,ox3,oy3);
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
					if (fabs(k)>0.06) k = 0.06*k/fabs(k);
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

int fitxy2(){
	int wid1,wid2;
	n_xy = 0;
	for (int j = 0; j<6; j++){
		for (int hid1 = 0; hid1<c_layerhit[j]; hid1++){
			int index1 = i_layerhit[j][hid1];
			wid1 = (*O_wireID)[index1];
			double k1 = f_k->Eval(j);
			double x1 = ((1+k1)*map_xro[j][wid1]+(1-k1)*map_xhv[j][wid1])/2.;
			double y1 = ((1+k1)*map_yro[j][wid1]+(1-k1)*map_yhv[j][wid1])/2.;
			double r1 = (*O_driftD)[index1];
			for (int hid2 = 0; hid2<c_layerhit[j+1]; hid2++){
				int index2 = i_layerhit[j+1][hid2];
				wid2 = (*O_wireID)[index2];
				double k2 = f_k->Eval(j+1);
				double x2 = ((1+k2)*map_xro[j+1][wid2]+(1-k2)*map_xhv[j+1][wid2])/2.;
				double y2 = ((1+k2)*map_yro[j+1][wid2]+(1-k2)*map_yhv[j+1][wid2])/2.;
				double r2 = (*O_driftD)[index2];
				getxy(x1,y1,r1,x2,y2,r2,aar_x[0][n_xy],aar_x[1][n_xy],aar_x[2][n_xy],aar_y[0][n_xy],aar_y[1][n_xy],aar_y[2][n_xy]);
				n_xy++;
			}
		}
	}
	if (g_x) delete g_x; g_x = new TGraph(n_xy,&(aar_y[0][0]),&(aar_x[0][0]));
	g_x->Fit("f_x","qN0","");

	if(gMinuit) delete gMinuit;
	gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
	gMinuit->SetFCN(fcn);
	arglist[0] = 0;
	gMinuit->SetPrintLevel(-1); // no print
	gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	double sliX = f_x->GetParameter(1);
	double iniX = f_x->GetParameter(0);
	gMinuit->mnparm(0, "slopeX", sliX, slerX/1.e3, -slerX,slerX,ierflg);
	gMinuit->mnparm(1, "interceptX", iniX, inerX/1.e3, -inerX,inerX,ierflg);

	// Now ready for minimization step
	arglist[0] = 500.0;
	arglist[1] = 1.0;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Print results
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	printf("====Rrestul====\n");
	gMinuit->mnprin(3,amin);

	double temp;
	gMinuit->GetParameter(0, sliX, temp);
	gMinuit->GetParameter(1, iniX, temp);

	f_x->SetParameter(0,iniX);
	f_x->SetParameter(1,sliX);
}

int fitxy(std::vector<int> * wireID,std::vector<double> * driftD){
	int wid1,wid2;
	n_xy = 0;
	for (int j = 1-thelayer%2; j<=7; j+=2){
		for (int hid1 = 0; hid1<c_layerhit[j]; hid1++){
			int index1 = i_layerhit[j][hid1];
			wid1 = (*wireID)[index1];
			double k1 = f_k->Eval(j);
			double x1 = ((1+k1)*map_xro[j][wid1]+(1-k1)*map_xhv[j][wid1])/2.;
			double y1 = ((1+k1)*map_yro[j][wid1]+(1-k1)*map_yhv[j][wid1])/2.;
			double r1 = (*driftD)[index1];
			for (int hid2 = 0; hid2<c_layerhit[j+2]; hid2++){
				int index2 = i_layerhit[j+2][hid2];
				wid2 = (*wireID)[index2];
				double k2 = f_k->Eval(j+2);
				double x2 = ((1+k2)*map_xro[j+2][wid2]+(1-k2)*map_xhv[j+2][wid2])/2.;
				double y2 = ((1+k2)*map_yro[j+2][wid2]+(1-k2)*map_yhv[j+2][wid2])/2.;
				double r2 = (*driftD)[index2];
				getxy(x1,y1,r1,x2,y2,r2,aar_x[0][n_xy],aar_x[1][n_xy],aar_x[2][n_xy],aar_y[0][n_xy],aar_y[1][n_xy],aar_y[2][n_xy]);
				n_xy++;
			}
		}
		if (thelayer%2&&thelayer+2<=7&&j+1==thelayer+2) j--;
		if (thelayer%2&&thelayer+2<=7&&j==thelayer+2) j--;
		if (thelayer%2&&thelayer+2>7&&j+1==thelayer-2) j--;
		if (thelayer%2&&thelayer+2>7&&j==thelayer-2) j--;
	}
	if (n_xy<=1) return 1;

	double chi2min = 1e9;
	int index = -1;
	for (int i = 0; i<pow(3,n_xy); i++){
		int a = i;
		for (int j = 0; j<n_xy; j++){
			int pow3j = pow(3,n_xy-j-1);
			int b = a/pow3j;
			a-=pow3j*b;
			ar_x[j] = aar_x[b][j];
			ar_y[j] = aar_y[b][j];
		}
		if (g_x) delete g_x; g_x = new TGraph(n_xy,&(ar_y[0]),&(ar_x[0]));
		g_x->Fit("f_x","qN0","");
		double chi2 = f_x->GetChisquare();
		if (chi2<chi2min){
			chi2min=chi2;
			index = i;
		}
	}
	int a = index;
	for (int j = 0; j<n_xy; j++){
		int pow3j = pow(3,n_xy-j-1);
		int b = a/pow3j;
		a-=pow3j*b;
		ar_x[j] = aar_x[b][j];
		ar_y[j] = aar_y[b][j];
	}
	if (g_x) delete g_x; g_x = new TGraph(n_xy,&(ar_y[0]),&(ar_x[0]));
	//if (g_x) delete g_x; g_x = new TGraph(n_xy,&(aar_y[0][0]),&(aar_x[0][0]));
	g_x->Fit("f_x","qN0","");
	return 0;
	//g_x->Print();
	//f_x->Print();
}

double t2x(double time, double sign){
	TF1 * f = 0;
//	if (sign>0) f = f_xt_all_right;
//	else f = f_xt_all_left;
	f = f_xt_all_right;
	double dd = f->Eval(time);
	if (dd<0) dd = 0;
	return dd;
	// FIXME
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

	if (argc<3){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	thelayer = (int)strtol(argv[2],NULL,10);
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
			errord[wp_lid][wp_wid] = 0.02;
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
		// FIXME
		//i_xt_t*=(228.+t0shift)/223.;
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
	f_xt_all_right = new TF1("f_xt_all_right","pol9",0-tres,driftTMax+tres);
	f_xt_all_left = new TF1("f_xt_all_left","pol9",0-tres,driftTMax+tres);
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
	std::vector<double> * i_driftT = 0;
	std::vector<int> * i_layerID = 0;
	std::vector<int> * i_wireID = 0;
	std::vector<int> * i_type = 0;
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
	buf<<"../root/i_"<<runNo<<suffix<<".root";
	TFile * of = new TFile(buf.str().c_str(),"RECREATE"); 
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
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	for ( int i = 0; i<N; i++){
		if (i%1000==0) std::cout<<(double)i/N*100<<"%..."<<std::endl;
		// FIXME
		//printf("*****Event %d*******\n",i);
		c->GetEntry(i);
		if(O_driftT) delete O_driftT; O_driftT = new std::vector<double>;
		if(O_driftD) delete O_driftD; O_driftD = new std::vector<double>;
		if(O_wireID) delete O_wireID; O_wireID = new std::vector<int>;
		if(O_layerID) delete O_layerID; O_layerID = new std::vector<int>;

		// get basical cdc hit information
		nHits= 0;
		int nGoodHits = 0;
		for (int j = 0; j<NLAY; j++){
			c_layerhit[j] = 0;
		}
		for (int ihit = 0; ihit<i_type->size(); ihit++){
			if ((*i_type)[ihit]!=0) continue;
			double dt = (*i_driftT)[ihit]+t0shift;
			if (dt<0-tres||dt>driftTMax+tres) continue;
			if (dt>driftTMax) dt = driftTMax;
			int lid = (*i_layerID)[ihit];
			int wid = (*i_wireID)[ihit];
			i_layerhit[lid][c_layerhit[lid]]=O_wireID->size();
			c_layerhit[lid]++;
			O_wireID->push_back(wid);
			O_layerID->push_back(lid);
			O_driftT->push_back(dt);
			O_driftD->push_back(t2x(dt,1));
			nHits++;
			if (lid!=thelayer){
				nGoodHits++;
			}
		}
		N_trigger++;

		int c_allLayerhit=0;
		for (int j = 1; j<NLAY; j++){
			//std::cout<<"  ["<<j<<"]: ";
			if (c_layerhit[j]&&j!=thelayer){
				N_layerhit[j]++;
				c_allLayerhit++;
				//std::cout<<"hit"<<std::endl;
			}
			else{
				//std::cout<<"no hit"<<std::endl;
			}
		}
		if (c_allLayerhit!=6) continue;
		N_allLayerhit++;
		if (nGoodHits!=6) continue;
		N_hit++;

		// fit y-z plane
		//if (fityz(O_wireID)) continue;
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
		//fitxy2();
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

		for ( int ihit = 0; ihit<O_driftT->size(); ihit++ ){
			int wid = (*O_wireID)[ihit];
			int lid = (*O_layerID)[ihit];
			double dt = (*O_driftT)[ihit];
			double ytemp = (map_yhv[lid][wid]+map_yro[lid][wid])/2.; // FIXME: without consideration of z
			double xtemp = f_x->Eval(ytemp);
			(*O_driftD)[ihit]= t2x(dt,xtemp);
		}
		N_found++;
		O_tx1 = f_x->Eval(yup);
		O_tx2 = f_x->Eval(ydown);
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
