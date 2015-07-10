#include <vector>
#include <sstream>
#include <iostream>
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

#define use_extra

std::vector<double> vec_i;
std::vector<double> vec_k;
double aar_x[3][50];
double aar_y[3][50];
double ar_x[50];
double ar_y[50];
int n_xy = 0;
TF1 * f_k = new TF1("f_k","pol1",-0.5,6.5);
TGraph * g_k = 0;
TF1 * f_x = new TF1("f_x","pol1",50,64);
TGraph * g_x = 0;
int i_layerhit[7][11];
int c_layerhit[7];
double map_xhv[7][11];
double map_yhv[7][11];
double map_xro[7][11];
double map_yro[7][11];
int map_wid[7][240];
double map_k[6][11][11];
int wmax[7];

double getk(int lid, int cid, double dd){
	int wid = map_wid[lid][cid];
	double xhv = map_xhv[lid][wid];
	double yhv = map_yhv[lid][wid];
	double xro = map_xro[lid][wid];
	double yro = map_yro[lid][wid];
	double y1 = 53;
	double y2 = 62.6;
	double k1 = f_k->Eval(0);
	double k2 = f_k->Eval(6);
	double x1 = f_x->Eval(y1);
	double x2 = f_x->Eval(y2);
	double z1 = k1*59.917/2;
	double z2 = k2*59.917/2;
	TVector3 vt1t2(x2-x1,y2-y1,z2-z1);
	TVector3 vw1w2(xro-xhv,yro-yhv,59.917);
	vt1t2 = vt1t2.Unit();
	vw1w2 = vw1w2.Unit();
	TVector3 vt1w2(xro-x1,yro-y1,59.917/2-z1);
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
	double z = k*59.917/2;
	double deltaz = dd>d0?sqrt(dd*dd-d0*d0)/sina:0;
	double zini = f_k->Eval(lid)*59.917/2;
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

void fityz(std::vector<int> * cellID, std::vector<double> * driftD){
	vec_k.clear();
	vec_i.clear();
	for (int j = 0; j<7; j++){
		for (int hid = 0; hid<c_layerhit[j]; hid++){
			int index = i_layerhit[j][hid];
			vec_i.push_back(j);
			vec_k.push_back(getk(j,(*cellID)[index],(*driftD)[index]));
		}
	}
	if (g_k) delete g_k; g_k = new TGraph(vec_k.size(),&(vec_i[0]),&(vec_k[0]));
	g_k->Fit("f_k","qN0","");
}

int fityz(std::vector<int> * cellID){
	vec_k.clear();
	vec_i.clear();
	int cid1,cid2,wid1,wid2;
	double k;
	int nlayers = 0;
	for (int j = 0; j<6; j++){
		bool check = false;
		for (int hid1 = 0; hid1<c_layerhit[j]; hid1++){
			int index1 = i_layerhit[j][hid1];
			cid1 = (*cellID)[index1];
			wid1 = map_wid[j][cid1];
			for (int hid2 = 0; hid2<c_layerhit[j+1]; hid2++){
				int index2 = i_layerhit[j+1][hid2];
				cid2 = (*cellID)[index2];
				wid2 = map_wid[j+1][cid2];
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

int fitxy(std::vector<int> * cellID,std::vector<double> * driftD){
	int cid1,cid2,wid1,wid2;
	n_xy = 0;
	for (int j = 0; j<6; j++){
		for (int hid1 = 0; hid1<c_layerhit[j]; hid1++){
			int index1 = i_layerhit[j][hid1];
			cid1 = (*cellID)[index1];
			wid1 = map_wid[j][cid1];
			double k1 = f_k->Eval(j);
			double x1 = ((1+k1)*map_xro[j][wid1]+(1-k1)*map_xhv[j][wid1])/2.;
			double y1 = ((1+k1)*map_yro[j][wid1]+(1-k1)*map_yhv[j][wid1])/2.;
			double r1 = (*driftD)[index1];
			for (int hid2 = 0; hid2<c_layerhit[j+1]; hid2++){
				int index2 = i_layerhit[j+1][hid2];
				cid2 = (*cellID)[index2];
				wid2 = map_wid[j+1][cid2];
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

	// for input
	TH1D * h_xt = (TH1D*)(new TFile("/home/chen/MyWorkArea/Simulate/comet/data/xt.root"))->Get("xt");

	TFile * ifile = new TFile("/home/chen/MyWorkArea/Simulate/proto/output/run.002.root");
	TTree * it = (TTree*) ifile->Get("tree");

	int CdcCell_nHits = 0;
	std::vector<double> * CdcCell_t = 0;
	std::vector<double> * CdcCell_wx = 0;
	std::vector<double> * CdcCell_wy = 0;
	std::vector<double> * CdcCell_wz = 0;
	std::vector<double> * CdcCell_x = 0;
	std::vector<double> * CdcCell_y = 0;
	std::vector<double> * CdcCell_z = 0;
	std::vector<double> * CdcCell_px = 0;
	std::vector<double> * CdcCell_py = 0;
	std::vector<double> * CdcCell_pz = 0;
	std::vector<double> * CdcCell_driftD = 0;
	std::vector<double> * CdcCell_driftDtrue = 0;
	std::vector<double> * CdcCell_tstart = 0;
	std::vector<double> * CdcCell_tstop = 0;
	std::vector<int> * CdcCell_cellID = 0;
	std::vector<int> * CdcCell_layerID = 0;
	std::vector<double> * CdcCell_edep = 0;
	std::vector<int> * CdcCell_tid = 0;
	std::vector<int> * CdcCell_posflag = 0;

	int M_nHits = 0;
	std::vector<int> * M_tid = 0;
	std::vector<double> * M_x = 0;
	std::vector<double> * M_y = 0;
	std::vector<double> * M_z = 0;
	std::vector<double> * M_t = 0;
	std::vector<int> * M_id = 0;
	std::vector<std::string> * M_name = 0;
	int evt_num;
	int run_num;

	it->SetBranchAddress("evt_num",&evt_num);
	it->SetBranchAddress("run_num",&run_num);
	it->SetBranchAddress("CdcCell_nHits",&CdcCell_nHits);
	it->SetBranchAddress("CdcCell_t",&CdcCell_t);
	it->SetBranchAddress("CdcCell_wx",&CdcCell_wx);
	it->SetBranchAddress("CdcCell_wy",&CdcCell_wy);
	it->SetBranchAddress("CdcCell_wz",&CdcCell_wz);
	it->SetBranchAddress("CdcCell_x",&CdcCell_x);
	it->SetBranchAddress("CdcCell_y",&CdcCell_y);
	it->SetBranchAddress("CdcCell_z",&CdcCell_z);
	it->SetBranchAddress("CdcCell_px",&CdcCell_px);
	it->SetBranchAddress("CdcCell_py",&CdcCell_py);
	it->SetBranchAddress("CdcCell_pz",&CdcCell_pz);
	it->SetBranchAddress("CdcCell_driftD",&CdcCell_driftD);
	it->SetBranchAddress("CdcCell_driftDtrue",&CdcCell_driftDtrue);
	it->SetBranchAddress("CdcCell_tstart",&CdcCell_tstart);
	it->SetBranchAddress("CdcCell_tstop",&CdcCell_tstop);
	it->SetBranchAddress("CdcCell_cellID",&CdcCell_cellID);
	it->SetBranchAddress("CdcCell_layerID",&CdcCell_layerID);
	it->SetBranchAddress("CdcCell_edep",&CdcCell_edep);
	it->SetBranchAddress("CdcCell_tid",&CdcCell_tid);
	it->SetBranchAddress("CdcCell_posflag",&CdcCell_posflag);

	it->SetBranchAddress("M_nHits",&M_nHits);
	it->SetBranchAddress("M_volName",&M_name);
	it->SetBranchAddress("M_volID",&M_id);
	it->SetBranchAddress("M_tid",&M_tid);
	it->SetBranchAddress("M_x",&M_x);
	it->SetBranchAddress("M_y",&M_y);
	it->SetBranchAddress("M_z",&M_z);
	it->SetBranchAddress("M_t",&M_t);

	# ifdef use_extra
	int wire_nHits = 0;
	int T_nHits = 0;
	it->SetBranchAddress("T_nHits",&T_nHits);
	it->SetBranchAddress("wire_nHits",&wire_nHits);
	# endif

	// For wire position
	TFile * TFile_wirepos = new TFile("../info/wire-position.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");

	int wp_hid;
	int wp_wid;
	int wp_lid;
	double wp_xhv;
	double wp_yhv;
	double wp_xro;
	double wp_yro;

	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("h",&wp_hid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
	TTree_wirepos->SetBranchAddress("xro",&wp_xro);
	TTree_wirepos->SetBranchAddress("yro",&wp_yro);

	for (int i = 0; i<7; i++){
		for (int j = 0; j<240; j++){
			map_wid[i][j] = -1;
		}
		wmax[i] = -1;
	}

	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		if (wp_lid>=1&&wp_lid<=7){
			map_xhv[wp_lid-1][wp_wid] = wp_xhv/10.;
			map_yhv[wp_lid-1][wp_wid] = wp_yhv/10.;
			map_xro[wp_lid-1][wp_wid] = wp_xro/10.;
			map_yro[wp_lid-1][wp_wid] = wp_yro/10.;
			map_wid[wp_lid-1][wp_hid/2] = wp_wid;
			if (wmax[wp_lid-1]<wp_wid) wmax[wp_lid-1] = wp_wid;
//			printf("x[%d][%d]: (%lf,%lf)\n",wp_lid,wp_hid/2,wp_xhv/10.,wp_xro/10.);
		}
	}
	for (int k = 0; k<6; k++){
		printf("layer %d & %d\n",k+1,k+2);
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

	// for output
	double O_mx = 0;
	double O_my = 0;
	double O_mz = 0;
	double O_mt = 0;
	int O_nHits = 0;
	std::vector<int> * O_hittype = 0;
	std::vector<double> * O_t = 0;
	std::vector<double> * O_tof = 0;
	std::vector<double> * O_wx = 0;
	std::vector<double> * O_wy = 0;
	std::vector<double> * O_wz = 0;
	std::vector<double> * O_wux = 0;
	std::vector<double> * O_wuy = 0;
	std::vector<double> * O_wcx = 0;
	std::vector<double> * O_wcy = 0;
	std::vector<double> * O_wdx = 0;
	std::vector<double> * O_wdy = 0;
	std::vector<double> * O_wfx = 0;
	std::vector<double> * O_wfy = 0;
	std::vector<double> * O_wfz = 0;
	std::vector<double> * O_x = 0;
	std::vector<double> * O_y = 0;
	std::vector<double> * O_z = 0;
	std::vector<double> * O_px = 0;
	std::vector<double> * O_py = 0;
	std::vector<double> * O_pz = 0;
	std::vector<double> * O_driftD = 0;
	std::vector<double> * O_driftDtrue = 0;
	std::vector<double> * O_tstop = 0;
	std::vector<double> * O_tstart = 0;
	std::vector<int> * O_cellID = 0;
	std::vector<int> * O_wireID = 0;
	std::vector<int> * O_layerID = 0;
	std::vector<double> * O_edep = 0;
	std::vector<int> * O_tid = 0;
	std::vector<int> * O_posflag = 0;
	double O_tx1;
	double O_tx2;
	double O_tz1;
	double O_tz2;
	double O_txr1;
	double O_tyr1;
	double O_tzr1;
	double O_sliX;
	double O_sliZ;
	double O_slirX;
	double O_slirZ;

	TFile * of = new TFile("output.root","RECREATE");
	TTree * ot = new TTree("tree","tree");

	ot->Branch("evt_num",&evt_num);
	ot->Branch("run_num",&run_num);
	ot->Branch("CdcCell_nHits",&O_nHits);
	ot->Branch("CdcCell_t",&O_t);
	ot->Branch("CdcCell_tx1",&O_tx1);
	ot->Branch("CdcCell_tx2",&O_tx2);
	ot->Branch("CdcCell_tz1",&O_tz1);
	ot->Branch("CdcCell_tz2",&O_tz2);
	ot->Branch("CdcCell_txr1",&O_txr1);
	ot->Branch("CdcCell_tyr1",&O_tyr1);
	ot->Branch("CdcCell_tzr1",&O_tzr1);
	ot->Branch("CdcCell_sliX",&O_sliX);
	ot->Branch("CdcCell_slirX",&O_slirX);
	ot->Branch("CdcCell_sliZ",&O_sliZ);
	ot->Branch("CdcCell_slirZ",&O_slirZ);
	ot->Branch("CdcCell_wx",&O_wx);
	ot->Branch("CdcCell_wy",&O_wy);
	ot->Branch("CdcCell_wz",&O_wz);
	ot->Branch("CdcCell_wfx",&O_wfx);
	ot->Branch("CdcCell_wfy",&O_wfy);
	ot->Branch("CdcCell_wfz",&O_wfz);
	ot->Branch("CdcCell_wux",&O_wux);
	ot->Branch("CdcCell_wuy",&O_wuy);
	ot->Branch("CdcCell_wcx",&O_wcx);
	ot->Branch("CdcCell_wcy",&O_wcy);
	ot->Branch("CdcCell_wdx",&O_wdx);
	ot->Branch("CdcCell_wdy",&O_wdy);
	ot->Branch("CdcCell_x",&O_x);
	ot->Branch("CdcCell_y",&O_y);
	ot->Branch("CdcCell_z",&O_z);
	ot->Branch("CdcCell_px",&O_px);
	ot->Branch("CdcCell_py",&O_py);
	ot->Branch("CdcCell_pz",&O_pz);
	ot->Branch("CdcCell_driftD",&O_driftD);
	ot->Branch("CdcCell_driftDtrue",&O_driftDtrue);
	ot->Branch("CdcCell_tof",&O_tof);
	ot->Branch("CdcCell_tstart",&O_tstart);
	ot->Branch("CdcCell_tstop",&O_tstop);
	ot->Branch("CdcCell_cellID",&O_cellID);
	ot->Branch("CdcCell_wireID",&O_wireID);
	ot->Branch("CdcCell_layerID",&O_layerID);
	ot->Branch("CdcCell_edep",&O_edep);
	ot->Branch("CdcCell_posflag",&O_posflag);
	ot->Branch("CdcCell_hittype",&O_hittype);

	ot->Branch("CdcCell_mx",&O_mx);
	ot->Branch("CdcCell_my",&O_my);
	ot->Branch("CdcCell_mz",&O_mz);
	ot->Branch("CdcCell_mt",&O_mt);

	ot->Branch("M_nHits",&M_nHits);
	ot->Branch("M_volName",&M_name);
	ot->Branch("M_volID",&M_id);
	ot->Branch("M_tid",&M_tid);
	ot->Branch("M_x",&M_x);
	ot->Branch("M_y",&M_y);
	ot->Branch("M_z",&M_z);
	ot->Branch("M_t",&M_t);

	# ifdef use_extra
	int O_T_nHits = 0;
	int O_wire_nHits = 0;
	ot->Branch("T_nHits",&O_T_nHits);
	ot->Branch("wire_nHits",&O_wire_nHits);
	# endif

	int dict[18][306];

	int N_layerhit[7];
	for (int i = 0; i<7; i++){
		N_layerhit[i] = 0;
	}
	int N_allLayerhit = 0;
	int N_trigger = 0;

	bool triggerd = false;
	double firsthittime = 0;
	int nGoodHit = 0;
	for ( int i = 0; i<it->GetEntries(); i++){
		for(int j = 0; j<18; j++){
			for (int k = 0; k<306; k++){
				dict[j][k]=-1;
			}
		}
		if (i%1000==0) printf("%lf%...\n",(double)i/it->GetEntries()*100);
		// FIXME
		//printf("*****Event %d*******\n",i);
		it->GetEntry(i);
		# ifdef use_extra
		O_T_nHits = T_nHits-1;
		O_wire_nHits = wire_nHits;
		# endif
		if (CdcCell_nHits==0) continue;
		O_nHits = 0;
		if(O_hittype) delete O_hittype; O_hittype = new std::vector<int>;
		if(O_t) delete O_t; O_t = new std::vector<double>;
		if(O_tof) delete O_tof; O_tof = new std::vector<double>;
		if(O_wx) delete O_wx; O_wx = new std::vector<double>;
		if(O_wy) delete O_wy; O_wy = new std::vector<double>;
		if(O_wz) delete O_wz; O_wz = new std::vector<double>;
		if(O_wfx) delete O_wfx; O_wfx = new std::vector<double>;
		if(O_wfy) delete O_wfy; O_wfy = new std::vector<double>;
		if(O_wfz) delete O_wfz; O_wfz = new std::vector<double>;
		if(O_wux) delete O_wux; O_wux = new std::vector<double>;
		if(O_wuy) delete O_wuy; O_wuy = new std::vector<double>;
		if(O_wcx) delete O_wcx; O_wcx = new std::vector<double>;
		if(O_wcy) delete O_wcy; O_wcy = new std::vector<double>;
		if(O_wdx) delete O_wdx; O_wdx = new std::vector<double>;
		if(O_wdy) delete O_wdy; O_wdy = new std::vector<double>;
		if(O_x) delete O_x; O_x = new std::vector<double>;
		if(O_y) delete O_y; O_y = new std::vector<double>;
		if(O_z) delete O_z; O_z = new std::vector<double>;
		if(O_px) delete O_px; O_px = new std::vector<double>;
		if(O_py) delete O_py; O_py = new std::vector<double>;
		if(O_pz) delete O_pz; O_pz = new std::vector<double>;
		if(O_driftD) delete O_driftD; O_driftD = new std::vector<double>;
		if(O_driftDtrue) delete O_driftDtrue; O_driftDtrue = new std::vector<double>;
		if(O_tstart) delete O_tstart; O_tstart = new std::vector<double>;
		if(O_tstop) delete O_tstop; O_tstop = new std::vector<double>;
		if(O_cellID) delete O_cellID; O_cellID = new std::vector<int>;
		if(O_wireID) delete O_wireID; O_wireID = new std::vector<int>;
		if(O_layerID) delete O_layerID; O_layerID = new std::vector<int>;
		if(O_edep) delete O_edep; O_edep = new std::vector<double>;
		if(O_tid) delete O_tid; O_tid = new std::vector<int>;
		if(O_posflag) delete O_posflag; O_posflag = new std::vector<int>;

		// get basical cdc hit information
		triggerd = false;
		firsthittime = 0;
		nGoodHit = 0;
		for (int j = 0; j<7; j++){
			c_layerhit[j] = 0;
		}
		for (int j = 0; j<CdcCell_nHits; j++){
			int lid = (*CdcCell_layerID)[j];
			if ((*CdcCell_tid)[j]==1&&lid<7&&map_wid[lid][(*CdcCell_cellID)[j]]!=-1){
				nGoodHit++;
				if (!firsthittime){
					firsthittime = (*CdcCell_t)[j];
				}
				i_layerhit[lid][c_layerhit[lid]]=j;
				//printf("i_layerhit[%d][%d] = %d\n",lid,c_layerhit[lid],j);
				c_layerhit[lid]=c_layerhit[lid]+1;
			}
		}
		if (nGoodHit<1) continue;

		// fit y-z plane
		if (fityz(CdcCell_cellID)) continue;
		//printf("After y-z fit 1:\n");
		//g_k->Print();
		//f_k->Print();
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int cid = (*CdcCell_cellID)[index];
		//		int wid = map_wid[j][cid];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = ((1+k_temp)*xro+(1-k_temp)*xhv)/2.;
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*59.917/2);
		//	}
		//}

		// fit x-y plane
		fitxy(CdcCell_cellID,CdcCell_driftD);
		//printf("After x-y fit 1:\n");
		//g_x->Print();
		//f_x->Print();
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int cid = (*CdcCell_cellID)[index];
		//		int wid = map_wid[j][cid];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = f_x->Eval(y_temp);
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*59.917/2);
		//	}
		//}

		// correct z positions and fit y-z again
		//fityz(CdcCell_cellID,CdcCell_driftD);
		//printf("After y-z fit 2:\n");
		//g_k->Print();
		//f_k->Print();
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int cid = (*CdcCell_cellID)[index];
		//		int wid = map_wid[j][cid];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = f_x->Eval(y_temp);
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*59.917/2);
		//	}
		//}

		// fit x-y again
		//fitxy(CdcCell_cellID,CdcCell_driftD);
		//printf("After x-y fit 2:\n");
		//for (int j = 0; j<7; j++){
		//	for (int hid = 0; hid<c_layerhit[j]; hid++){
		//		int index = i_layerhit[j][hid];
		//		int cid = (*CdcCell_cellID)[index];
		//		int wid = map_wid[j][cid];
		//		double xhv = map_xhv[j][wid];
		//		double yhv = map_yhv[j][wid];
		//		double xro = map_xro[j][wid];
		//		double yro = map_yro[j][wid];
		//		double k_temp = f_k->Eval(j);
		//		double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
		//		double x_temp = f_x->Eval(y_temp);
		//		printf("  hit[%d,%d,%d]: @(%d,%d,%d) (%lf,\t%lf,\t%lf)\n",j,hid,index,j,cid,wid,x_temp,y_temp,k_temp*59.917/2);
		//	}
		//}

		// get trigger information
		// FIXME
		//M_nHits = 0;
		for (int j = 0; j<M_nHits; j++){
			if ((*M_tid)[j]==1){
				if ((*M_t)[j]<firsthittime){
					triggerd = true;
					O_mt = (*M_t)[j];
					O_mx = (*M_x)[j];
					O_my = (*M_y)[j];
					O_mz = (*M_z)[j];
					break;
				}
			}
		}
		//FIXME
		if (M_nHits<2||!triggerd) continue;
		N_trigger++;
		int c_allLayerhit=0;
		for (int j = 0; j<7; j++){
			if (c_layerhit[j]){
				N_layerhit[j]++;
				c_allLayerhit++;
			}
		}
		if (c_allLayerhit==7) N_allLayerhit++;

		// sort according to tof
		std::vector<int> hitindice;
		hitindice.resize(CdcCell_nHits);
		for (int j = 0; j<CdcCell_nHits; j++){
			hitindice[j] = j;
		}
		int temp;
		for (int j = 0; j<CdcCell_nHits; j++){
			for (int k = j+1; k<CdcCell_nHits; k++){
				if ((*CdcCell_t)[hitindice[j]]>(*CdcCell_t)[hitindice[k]]){
					temp = hitindice[j];
					hitindice[j] = hitindice[k];
					hitindice[k]=temp;
				}
			}
		}

		double hittime;
		double tof;
		double starttime;
		double stoptime;
		int hittype;
		// FIXME
		for (int j = 0; j<CdcCell_nHits; j++){
			if ((*CdcCell_tid)[hitindice[j]]==1){
				hittype = 0;
			}
			else{
				hittype = 1;
			}
			starttime = (*CdcCell_tstart)[hitindice[j]];
			stoptime = (*CdcCell_tstop)[hitindice[j]];
			hittime = (*CdcCell_t)[hitindice[j]];
			tof = (*CdcCell_t)[hitindice[j]];
			int layerID = (*CdcCell_layerID)[hitindice[j]];
			int cellID = (*CdcCell_cellID)[hitindice[j]];
			int wireID = map_wid[layerID][cellID];
			if (layerID>=7||wireID==-1){ 
				//printf("[%d,%d]->[%d,%d] out of range!\n",layerID,cellID,layerID+1,cellID*2+1);
				continue;
			}
			//else{
			//	printf("[%d,%d]->[%d,%d] Good!\n",layerID,cellID,layerID+1,cellID*2+1);
			//}
			if (dict[layerID][cellID]==-1){
				double px = (*CdcCell_px)[hitindice[j]];
				double py = (*CdcCell_py)[hitindice[j]];
				double pz = (*CdcCell_pz)[hitindice[j]];
				double pa = sqrt(px*px+py*py+pz*pz);
				//FIXME
				//if (tof>7||pa<0.103||hittype!=0) continue; // only first turn
				if (hittype!=0) continue; // only hits from the signal track

				dict[layerID][cellID]=O_nHits;
				O_nHits++;
				O_t->push_back(hittime);
				O_tof->push_back(tof);
				O_tstart->push_back(starttime);
				O_tstop->push_back(stoptime);
				O_wx->push_back((*CdcCell_wx)[hitindice[j]]);
				O_wy->push_back((*CdcCell_wy)[hitindice[j]]);
				O_wz->push_back((*CdcCell_wz)[hitindice[j]]);
				double xhv = map_xhv[layerID][wireID];
				double yhv = map_yhv[layerID][wireID];
				double xro = map_xro[layerID][wireID];
				double yro = map_yro[layerID][wireID];
				O_wux->push_back(xhv);
				O_wuy->push_back(yhv);
				O_wdx->push_back(xro);
				O_wdy->push_back(yro);
				O_wcx->push_back((xhv+xro)/2.);
				O_wcy->push_back((yhv+yro)/2.);
				double k_temp = f_k->Eval(layerID);
				double y_temp = ((1+k_temp)*yro+(1-k_temp)*yhv)/2.;
				double x_temp = ((1+k_temp)*xro+(1-k_temp)*xhv)/2.;
//				double x_temp = f_x->Eval(y_temp);
				O_wfx->push_back(x_temp);
				O_wfy->push_back(y_temp);
				O_wfz->push_back(k_temp*59.917/2);
				O_x->push_back((*CdcCell_x)[hitindice[j]]);
				O_y->push_back((*CdcCell_y)[hitindice[j]]);
				O_z->push_back((*CdcCell_z)[hitindice[j]]);
				O_px->push_back((*CdcCell_px)[hitindice[j]]);
				O_py->push_back((*CdcCell_py)[hitindice[j]]);
				O_pz->push_back((*CdcCell_pz)[hitindice[j]]);
				// FIXME: we can add smear
				//O_driftD->push_back((*CdcCell_driftD)[hitindice[j]]+gRandom->Gaus(0,0.02)); // no need to monify, by now...
				O_driftD->push_back((*CdcCell_driftD)[hitindice[j]]); // no need to monify, by now...
				O_driftDtrue->push_back((*CdcCell_driftDtrue)[hitindice[j]]);
				O_cellID->push_back(cellID);
				O_wireID->push_back(wireID);
				O_layerID->push_back(layerID);
				O_edep->push_back((*CdcCell_edep)[hitindice[j]]);
				O_posflag->push_back((*CdcCell_posflag)[hitindice[j]]);
				O_hittype->push_back(hittype);
				if (layerID==6){
					O_txr1 = (*CdcCell_x)[hitindice[j]];
					O_tyr1 = (*CdcCell_y)[hitindice[j]];
					O_tzr1 = (*CdcCell_z)[hitindice[j]];
				}
			}
			else{
				std::cout<<"WTH?"<<std::endl;
			}
		}
		O_tx1 = f_x->Eval(53);
		O_tx2 = f_x->Eval(62.6);
		O_tz1 = f_k->Eval(0)*59.917/2;
		O_tz2 = f_k->Eval(6)*59.917/2;
		O_sliX = (O_tx2-O_tx1)/9.6;
		O_sliZ = (O_tz2-O_tz1)/9.6;
		O_slirX = -(*CdcCell_px)[hitindice[0]]/(*CdcCell_py)[hitindice[0]];
		O_slirZ = -(*CdcCell_pz)[hitindice[0]]/(*CdcCell_py)[hitindice[0]];
		//FIXME
		if (O_nHits<=0) continue;

		ot->Fill();
		//std::cout<<"Filled! ot->GetEntries() = "<<ot->GetEntries()<<std::endl;
	}// end of event loop
	ot->Write();
	of->Close();
	printf("Triggered Events: %d\n",N_trigger);
	printf("layer 1 hit Events: %d\n",N_layerhit[0]);
	printf("layer 2 hit Events: %d\n",N_layerhit[1]);
	printf("layer 3 hit Events: %d\n",N_layerhit[2]);
	printf("layer 4 hit Events: %d\n",N_layerhit[3]);
	printf("layer 5 hit Events: %d\n",N_layerhit[4]);
	printf("layer 6 hit Events: %d\n",N_layerhit[5]);
	printf("layer 7 hit Events: %d\n",N_layerhit[6]);
	printf("All layer hit Events: %d\n",N_allLayerhit);
	return 0;
}
