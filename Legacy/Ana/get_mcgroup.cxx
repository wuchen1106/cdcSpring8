#include "TH1D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TF1.h"
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define NCELL 250
#define NLAY 9
#define NWIRE 11

int main(int argc, char** argv){

	//===================Get wire position============================
	// For wire position
	double Lchamebr = 599.17; // mm
	double zhv = -Lchamebr/2;
	double zro = Lchamebr/2;
	double yup=623.97007;
	double U = 8; // mm
	TFile * TFile_wirepos = new TFile("../info/wire-position.root");
	TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
	int wp_wid;
	int wp_lid;
	int wp_cid;
	double wp_xhv;
	double wp_yhv;
	double wp_xro;
	double wp_yro;
	TTree_wirepos->SetBranchAddress("l",&wp_lid);
	TTree_wirepos->SetBranchAddress("h",&wp_cid);
	TTree_wirepos->SetBranchAddress("w",&wp_wid);
	TTree_wirepos->SetBranchAddress("xro",&wp_xro);
	TTree_wirepos->SetBranchAddress("yro",&wp_yro);
	TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
	TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
	double map_wid[NLAY][NCELL];
	double map_xc[NLAY][NWIRE];
	double map_yc[NLAY][NWIRE];
	double map_xhv[NLAY][NWIRE];
	double map_yhv[NLAY][NWIRE];
	double map_xro[NLAY][NWIRE];
	double map_yro[NLAY][NWIRE];
	int map_check[NLAY][NCELL];
	for(int ilayer = 0; ilayer<NLAY; ilayer++){
		for (int icell = 0; icell<NCELL; icell++){
			map_check[ilayer][icell]=0;
		}
	}
	for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
		TTree_wirepos->GetEntry(i);
		wp_cid=wp_cid/2;
		if (wp_lid>=1&&wp_lid<NLAY){
			map_xhv[wp_lid][wp_wid] = wp_xhv;
			map_yhv[wp_lid][wp_wid] = wp_yhv;
			map_xro[wp_lid][wp_wid] = wp_xro;
			map_yro[wp_lid][wp_wid] = wp_yro;
			map_xc[wp_lid][wp_wid] = (wp_xhv+wp_xro)/2.;
			map_yc[wp_lid][wp_wid] = (wp_yhv+wp_yro)/2.;
			map_wid[wp_lid][wp_cid] = wp_wid;
			map_check[wp_lid][wp_cid] = 1;
		}
	}

	int nEventMax = 0;
	TChain* ichain = new TChain("tree");
	ichain->Add("/home/chen/MyWorkArea/Simulate/proto/output/signal.general.root");

	std::vector<int> * tritid = 0;
	std::vector<int> * tid = 0;
	std::vector<int> * layerID = 0;
	std::vector<int> * cellID = 0;
	std::vector<double> * driftD = 0;
	std::vector<double> * edep = 0;
	std::vector<double> * wx = 0;
	std::vector<double> * x = 0;
	std::vector<double> * y = 0;
	std::vector<double> * z = 0;
	std::vector<double> * px = 0;
	std::vector<double> * py = 0;
	std::vector<double> * pz = 0;

	ichain->SetBranchAddress("M_tid",&tritid);
	ichain->SetBranchAddress("CdcCell_wx",&wx);
	ichain->SetBranchAddress("CdcCell_x",&x);
	ichain->SetBranchAddress("CdcCell_y",&y);
	ichain->SetBranchAddress("CdcCell_z",&z);
	ichain->SetBranchAddress("CdcCell_px",&px);
	ichain->SetBranchAddress("CdcCell_py",&py);
	ichain->SetBranchAddress("CdcCell_pz",&pz);
	ichain->SetBranchAddress("CdcCell_edep",&edep);
	ichain->SetBranchAddress("CdcCell_tid",&tid);
	ichain->SetBranchAddress("CdcCell_layerID",&layerID);
	ichain->SetBranchAddress("CdcCell_cellID",&cellID);
	ichain->SetBranchAddress("CdcCell_driftD",&driftD);

	TFile *ofile = new TFile("output.root","RECREATE");
	TTree * ot= new TTree("t","t");
	int nGoodPairs;
	int nOKPairs;
	int nLayersWith1GoodHit;
	int nLayersWith1OKHit;
	int nLayersWith1Hit;
	int nLayersWithMoreHits;
	int nHitsGood;
	int nHitsOK;
	int nHits;
	int nHitsAndNoise;
	int nTracks;
	int nTriHits;
	double inx;
	double inz;
	double slx;
	double slz;
	std::vector<double> * zmc= 0;
	std::vector<double> * zcalc = 0;
	std::vector<double> * zy = 0;
	std::vector<double> * zdd1 = 0;
	std::vector<double> * zdd2 = 0;
	std::vector<int> * zndiff = 0;
	std::vector<int> * zwid1 = 0;
	std::vector<int> * zwid2 = 0;
	std::vector<int> * zlid1 = 0;
	std::vector<int> * zlid2 = 0;
	double chi2;
	double slzc;
	double inzc;
	int nz;
	ot->Branch("ntri",&nTriHits);
	ot->Branch("nt",&nTracks);
	ot->Branch("ngp",&nGoodPairs);
	ot->Branch("nop",&nOKPairs);
	ot->Branch("n1g",&nLayersWith1GoodHit);
	ot->Branch("n1o",&nLayersWith1OKHit);
	ot->Branch("n1",&nLayersWith1Hit);
	ot->Branch("nm",&nLayersWithMoreHits);
	ot->Branch("N",&nHits);
	ot->Branch("No",&nHitsOK);
	ot->Branch("Ng",&nHitsGood);
	ot->Branch("Na",&nHitsAndNoise);
	ot->Branch("inx",&inx);
	ot->Branch("inz",&inz);
	ot->Branch("slx",&slx);
	ot->Branch("slz",&slz);
	ot->Branch("chi2",&chi2);
	ot->Branch("inzc",&inzc);
	ot->Branch("slzc",&slzc);
	ot->Branch("zmc",&zmc);
	ot->Branch("zcalc",&zcalc);
	ot->Branch("zndiff",&zndiff);
	ot->Branch("zdd1",&zdd1);
	ot->Branch("zdd2",&zdd2);
	ot->Branch("zwid1",&zwid1);
	ot->Branch("zlid1",&zlid1);
	ot->Branch("zwid2",&zwid2);
	ot->Branch("zlid2",&zlid2);
	ot->Branch("nz",&nz);

	int isgood[NLAY];
	int isok[NLAY];
	int isall[NLAY];

	TGraph * g_z = 0;
	TF1 * f_z = new TF1("f_z","pol1",500,640);

	Long64_t N = ichain->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" events..."<<std::endl;
	for ( int i = 0 ; i<N; i++){
		if (i%1000==0) printf("%lf%...\n",(double)i/N*100);
		ichain->GetEntry(i);
		// is triggered?
		nTriHits = 0;
		for (int ihit = 0; ihit<tritid->size(); ihit++){
			if ((*tritid)[ihit]==1) nTriHits++;
		}

		nOKPairs=0;
		nGoodPairs=0;
		nLayersWith1GoodHit=0;
		nLayersWith1OKHit=0;
		nLayersWith1Hit=0;
		nLayersWithMoreHits=0;
		nHitsGood=0;
		nHitsOK=0;
		nHits=0;
		nHitsAndNoise=0;
		nTracks = 0;
		for (int ilayer = 0; ilayer<NLAY; ilayer++){
			isgood[ilayer] = 0;
			isok[ilayer] = 0;
			isall[ilayer] = 0;
		}
		int pretid = 0;
		inx=1000;
		inz=1000;
		slx=1000;
		slz=1000;
		for (int ihit = 0; ihit<layerID->size(); ihit++){
			if ((*edep)[ihit]<1.85-07) continue; // too small energy deposit
			int lid = (*layerID)[ihit]+1;
			if (!map_check[lid][(*cellID)[ihit]]) continue; // this wire is not read out
			if (lid==8){ // the outmost layer
				inx = (*x)[ihit]*10;
				inz = (*z)[ihit]*10;
				slx = (*px)[ihit]/(*py)[ihit];
				slz = (*pz)[ihit]/(*py)[ihit];
//				inx += (yup-(*y)[ihit]*10)*slx;
//				inz += (yup-(*y)[ihit]*10)*slz;
			}
			if ((*tid)[ihit]!=pretid){
				pretid=(*tid)[ihit];
				nTracks++;
			}
			if (lid==0||lid==5) continue;
			if (((*driftD)[ihit])>0.1&&((*driftD)[ihit])<0.77){
				isgood[lid]++;
				nHitsGood++;
			}
			if (((*driftD)[ihit])>0.038&&((*driftD)[ihit])<0.81){
				isok[lid]++;
				nHitsOK++;
			}
			if (((*driftD)[ihit])>0.038&&((*driftD)[ihit])<1.6){
				isall[lid]++;
				nHits++;
			}
			nHitsAndNoise++;
		}
		if (i==103)std::cout<<"inx = "<<inx<<std::endl;
		if (nTracks==1&&pretid!=1) nTracks = 0; // no signal track!
		for (int ilayer = 0; ilayer<NLAY; ilayer++){
			if (isgood[ilayer]==1) nLayersWith1GoodHit++;
			if (isall[ilayer]==1) nLayersWith1Hit++;
			if (isall[ilayer]>=3) nLayersWithMoreHits++;
			if (isok[ilayer]==1) nLayersWith1OKHit++;
			if (ilayer<NLAY-1){
				if (isgood[ilayer]==1&&isgood[ilayer+1]==1) nGoodPairs++;
				if (isok[ilayer]==1&&isok[ilayer+1]==1) nOKPairs++;
			}
		}
		if (zmc) delete zmc; zmc = new std::vector<double>;
		if (zcalc) delete zcalc; zcalc = new std::vector<double>;
		if (zy) delete zy; zy = new std::vector<double>;
		if (zndiff) delete zndiff; zndiff = new std::vector<int>;
		if (zdd1) delete zdd1; zdd1 = new std::vector<double>;
		if (zdd2) delete zdd2; zdd2 = new std::vector<double>;
		if (zwid1) delete zwid1; zwid1 = new std::vector<int>;
		if (zlid1) delete zlid1; zlid1 = new std::vector<int>;
		if (zwid2) delete zwid2; zwid2 = new std::vector<int>;
		if (zlid2) delete zlid2; zlid2 = new std::vector<int>;
		chi2=1e9;
		inzc=0;
		slzc=0;
		nz = 0;
		if (nHitsOK==nLayersWith1OKHit&&nHitsOK>=2){
//			for (int ihit = layerID->size()-1; ihit>0; ihit--){
//				if (d_down>8.1||d_down<0.38) continue;
//				v_index.push_back(ihit);
//			}
			for (int ihit = layerID->size()-1; ihit>0; ihit--){
				int lid_down = (*layerID)[ihit]+1;
				if (lid_down==5) continue;
				double d_down = (*driftD)[ihit]*10;
				if (d_down>8.1||d_down<0.38) continue;
				for (int jhit = ihit-1; jhit>0; jhit--){
					int lid_up = (*layerID)[jhit]+1;
					if (lid_up==lid_down) continue;
					if (lid_up<lid_down) continue;
					if (lid_up-lid_down>1) continue;
					if (lid_up==5) continue;
					double d_up = (*driftD)[jhit]*10;
					if (d_up>8.1||d_up<0.38) continue;
					int wid_down = map_wid[lid_down][(*cellID)[ihit]];
					int wid_up = map_wid[lid_up][(*cellID)[jhit]];
					if ((*x)[ihit]<(*wx)[ihit]) d_down*=-1;
					if ((*x)[jhit]<(*wx)[jhit]) d_up*=-1;
					double deltax = (d_down-d_up);
					//double deltax = (d_down-d_up)/cos(atan(slx))+(lid_up-lid_down)*16*sin(atan(slx));
					double xuro = map_xro[lid_up][wid_up];
					double xuhv = map_xhv[lid_up][wid_up];
					double xdro = map_xro[lid_down][wid_down];
					double xdhv = map_xhv[lid_down][wid_down];
					double k = 0;
					if (xuhv-xdhv==xuro-xdro) continue;
					else
						k = (2*deltax-(xuhv-xdhv)-(xuro-xdro))/(-(xuhv-xdhv)+(xuro-xdro));
					if (fabs(k)>1) continue;
					zcalc->push_back(k*Lchamebr/2.);
					double yhvu = map_yhv[lid_up][wid_up];
					double yhvd = map_yhv[lid_down][wid_down];
					double yrou = map_yro[lid_up][wid_up];
					double yrod = map_yro[lid_down][wid_down];
					zy->push_back((yhvu+yhvd)*(1-k)/4.+(yrou+yrod)*(1+k)/4.);
					zmc->push_back(((*z)[ihit]+(*z)[jhit])/2.*10);
					zndiff->push_back(lid_up-lid_down);
					zwid1->push_back(wid_up);
					zwid2->push_back(wid_down);
					zdd1->push_back(d_up);
					zdd2->push_back(d_down);
					zlid1->push_back(lid_up);
					zlid2->push_back(lid_down);
					//std::cout<<i<<": "<<k<<", "<<k*Lchamebr/2.+(d_up-d_down)*pow(-1,lid_down)*3.938<<std::endl;
					//std::cout<<"  "<<ihit<<": ["<<lid_down<<","<<wid_down<<"], "<<(*x)[ihit]<<", "<<(*wx)[ihit]<<", "<<(*z)[ihit]<<std::endl;
					//std::cout<<"  "<<jhit<<": ["<<lid_up<<","<<wid_up<<"], "<<(*x)[jhit]<<", "<<(*wx)[jhit]<<", "<<(*z)[jhit]<<std::endl;
					break;
				}
			}
			if (zcalc->size()>=2){
				if (g_z) delete g_z; g_z = new TGraph(zcalc->size(),&((*zy)[0]),&((*zcalc)[0]));
				g_z->Fit("f_z","qN0","");
				chi2=f_z->GetChisquare();
				inzc=f_z->Eval(yup);
				slzc=f_z->GetParameter(1);
			}
			else if (zcalc->size()==1){
				inzc=(*zcalc)[0];
				slzc=0;
			}
			else if (zcalc->size()==1){
				inzc=0;
				slzc=0;
			}
			nz=zcalc->size();
		}
		ot->Fill();
	}
	ot->Write();
	ofile->Close();
	return 0;
}
