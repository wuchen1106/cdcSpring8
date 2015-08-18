#include <iostream>
#include <iomanip>  
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TEllipse.h"
#include "TChain.h"
#include "TLine.h"
#include "TText.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TVector3.h"

//#define DEBUG
//#define MC

#define NSAM 32
#define NCHS 48
#define NCHT 96
#define MAXNHIT 40
#define NBD 2
#define MIN_ADC 180
#define MAX_ADC 700
#define NBINS  256
#define NLAY    8
#define NCEL    11

double PI = TMath::Pi();
int power2_15 = pow(2,15);

// ________About XT___________
TF1 * XTL_func[NLAY][NCEL];
TF1 * XTR_func[NLAY][NCEL];

// ________About Cell___________
Double_t U=0.8;
std::vector<std::vector<double> > Yhv,Xhv,errord;// High Voltage side
std::vector<std::vector<double> > Yro,Xro; // Read Out side
double zhv = -59.917/2;
double zro = 59.917/2;

// ________About Track___________
int thelayer = 5;
std::vector<double> * i_driftD = 0;
std::vector<double> * i_driftT = 0;
std::vector<int> * i_wireID = 0;
std::vector<int> * i_layerID = 0;
std::vector<int> * i_peak = 0;
std::vector<double> * i_sum = 0;
double yup = 62.397007;
double ydown = 52.760011;
double deltay = yup-ydown;
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

// ________About Fitting___________
TMinuit *gMinuit = 0;
Double_t arglist[10];
Int_t ierflg = 0;
Double_t amin,edm,errdef;
Int_t nvpar,nparx,icstat;
double slerZ = 0.42*2;
double slerX = 0.42*2;
double inerZ = 2*2;
double inerX = 2*2;

//______________________________________________________________________________
Double_t get_dist(int lid, int wid, Double_t slx, Double_t inx, Double_t slz, Double_t inz);
void getchi2(Double_t &f, Double_t slx, Double_t inx, Double_t slz, Double_t inz,bool all = false);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void do_fit(Double_t slix, Double_t inix,Double_t sliz, Double_t iniz);
void print_usage(char* prog_name);

//______________________________________________________________________________
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
	int nEventMax = 0;
	if (argc>=4) nEventMax = (int)strtol(argv[3],NULL,10);
	int iterationNo = 0;
	if (argc>=5) iterationNo = (int)strtol(argv[4],NULL,10);
	int N_MIN = 0;
	if (argc>=6) N_MIN = (int)strtol(argv[5],NULL,10);
	double chi2Xmax,chi2Ymax;
	chi2Xmax=1e9;
	if (argc>=7) chi2Xmax = (int)strtol(argv[6],NULL,10);
	chi2Ymax=chi2Xmax;

	double eff_p2_average = 0;
	double eff_p3_average = 0;

	//===================Get run info============================
//	TFile * if_run = new TFile("../run_summary/run.eff.root");
//	TTree * t_run = (TTree*) if_run->Get("t");
//	double i_runNo, HVTXY1,HVTXY2,HVP2,HVP3,THRTXY1,THRTXY2,THRP2,THRP3;
//	t_run->SetBranchAddress("runNo",&i_runNo);
//	t_run->SetBranchAddress("HVTXY1",&HVTXY1);
//	t_run->SetBranchAddress("HVTXY2",&HVTXY2);
//	t_run->SetBranchAddress("HVP2",&HVP2);
//	t_run->SetBranchAddress("HVP3",&HVP3);
//	t_run->SetBranchAddress("THRTXY1",&THRTXY1);
//	t_run->SetBranchAddress("THRTXY2",&THRTXY2);
//	t_run->SetBranchAddress("THRP2",&THRP2);
//	t_run->SetBranchAddress("THRP3",&THRP3);
//	//for(int i = 0; i<t_run->GetEntries(); i++){
//	for(int i = 0; i<10; i++){
//		t_run->GetEntry(i);
//		if (i_runNo == runNo) break;
//	}
//	std::cout<<"runNo#"<<runNo<<": "<<(int)HVTXY1<<" V"<<std::endl;
//	TString gastype = "He:CH_{4}(80:20)";
//	double npair = 17.96;
//	if (runNo>=227||runNo<=147){
//		gastype = "He:iC_{4}H_{10}(90:10)";
//		npair = 27.96;
//	}
//	else if (runNo<=194){
//		gastype = "He:C_{2}H_{4}(50:50)";
//		npair = 56.10;
//	}

	//===================Get a2c============================
	TF1 * f_a2c = new TF1("a2c","5.98739+2.4852*x+0.000573394*x*x-5.21769e-05*x*x*x+3.05897e-07*x*x*x*x-7.54057e-10*x*x*x*x*x+8.60252e-13*x*x*x*x*x*x-3.68603e-16*x*x*x*x*x*x*x",-10,800);

	//===================Set Geometry============================
	// The z values
	// get the inital values
	TFile * if_geom = new TFile("../info/wire-position.v3.root");
	TTree * t_geom = (TTree*) if_geom->Get("t");
	double x1_geom, y1_geom, z1_geom, x2_geom, y2_geom, z2_geom;
	int wid_geom,lid_geom;
	t_geom->SetBranchAddress("l",&lid_geom);
	t_geom->SetBranchAddress("w",&wid_geom);
	t_geom->SetBranchAddress("xhv",&x1_geom);
	t_geom->SetBranchAddress("yhv",&y1_geom);
	t_geom->SetBranchAddress("xro",&x2_geom);
	t_geom->SetBranchAddress("yro",&y2_geom);
	Float_t error = 0.02;
	for ( int i = 0; i<NLAY; i++){
		std::vector<double> xhv;
		std::vector<double> yhv;
		std::vector<double> xro;
		std::vector<double> yro;
		std::vector<double> err;
		for (int j = 0; j<NCEL; j++){
			yhv.push_back(0);
			xhv.push_back(0);
			yro.push_back(0);
			xro.push_back(0);
			err.push_back(0);
		}
		Yhv.push_back(yhv);
		Xhv.push_back(xhv);
		Yro.push_back(yro);
		Xro.push_back(xro);
		errord.push_back(err);
	}
	for(int i = 0; i<t_geom->GetEntries(); i++){
		t_geom->GetEntry(i);
		if (lid_geom<1||lid_geom>7) continue;
		Xhv[lid_geom][wid_geom] = x1_geom/10.;
		Yhv[lid_geom][wid_geom] = y1_geom/10.;
		Xro[lid_geom][wid_geom] = x2_geom/10.;
		Yro[lid_geom][wid_geom] = y2_geom/10.;
		errord[lid_geom][wid_geom]=error;
	}
	if_geom->Close();

	//===================Prepare Histograms============================
	//chi2
	TH1D * hist_chi2 = new TH1D("chi2","chi2",128,0,20);

	//===================Get ROOT File============================
	//TChain * c = new TChain("t","t");
	TChain * c = new TChain("t","t");
	std::stringstream buf;
	buf.str(""); buf.clear();
	buf<<"../root/i_"<<runNo<<suffix<<".root";
	c->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	int triggerNumber;
	double xup, xdown;
	double zup, zdown;
	c->SetBranchAddress("wireID",&i_wireID);
	c->SetBranchAddress("layerID",&i_layerID);
	c->SetBranchAddress("driftD",&i_driftD);
	c->SetBranchAddress("driftT",&i_driftT);
	c->SetBranchAddress("peak",&i_peak);
	c->SetBranchAddress("sum",&i_sum);
	c->SetBranchAddress("tx1",&xup);
	c->SetBranchAddress("tz1",&zup);
	c->SetBranchAddress("tx2",&xdown);
	c->SetBranchAddress("tz2",&zdown);
	c->SetBranchAddress("triggerNumber",&triggerNumber);
#ifdef MC
	double slirX, inirX, inirY;
	double slirZ, inirZ;
	c->SetBranchAddress("slirX",&slirX);
	c->SetBranchAddress("slirZ",&slirZ);
	c->SetBranchAddress("txr1",&inirX);
	c->SetBranchAddress("tyr1",&inirY);
	c->SetBranchAddress("tzr1",&inirZ);
#endif

	//===================Output file============================
	buf.str(""); buf.clear();
	//buf<<"../root/t_"<<runNo<<"."<<suffix<<iterationNo<<".root";
	buf<<"../root/t_"<<runNo<<suffix<<".root";
	TFile * f = new TFile(buf.str().c_str(),"RECREATE"); 
	TTree * t = new TTree("t","t");
	Double_t chi2;
	Double_t chi2i;
	Double_t slZ, inZ;
	Double_t slX, inX;
	Double_t sliZ, iniZ;
	Double_t sliX, iniX;
	std::vector<double> * o_fitD = 0;

	t->Branch("triggerNumber",&triggerNumber);

	t->Branch("driftD",&i_driftD);
	t->Branch("driftT",&i_driftT);
	t->Branch("wireID",&i_wireID);
	t->Branch("layerID",&i_layerID);
	t->Branch("peak",&i_peak);
	t->Branch("sum",&i_sum);

	t->Branch("fitD",&o_fitD);
	t->Branch("chi2",&chi2);
	t->Branch("chi2i",&chi2i);
	t->Branch("slZ",&slZ);
	t->Branch("slX",&slX);
	t->Branch("inX",&inX);
	t->Branch("inZ",&inZ);
	t->Branch("sliZ",&sliZ);
	t->Branch("sliX",&sliX);
	t->Branch("iniX",&iniX);
	t->Branch("iniZ",&iniZ);

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (Long64_t iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		c->GetEntry(iEvent);

		//===================Set Initial============================
		sliX = (xup-xdown)/deltay;
		sliZ = (zup-zdown)/deltay;
		iniX = xup;
		iniZ = zup;
		//getchi2(chi2i,sliX,iniX,sliZ,iniZ,true);
		getchi2(chi2i,sliX,iniX,sliZ,iniZ);

		//===================Do The Fitting============================
		do_fit(sliX,iniX,sliZ,iniZ);
		double temp;
		gMinuit->GetParameter(0, slX, temp);
		gMinuit->GetParameter(1, inX, temp);
		gMinuit->GetParameter(2, slZ, temp);
		gMinuit->GetParameter(3, inZ, temp);
		//getchi2(chi2,slX,inX,slZ,inZ,true);
		getchi2(chi2,slX,inX,slZ,inZ);

		//===================Get Info for output============================
		o_fitD = new std::vector<double>;
		for (int ihit = 0; ihit<i_driftD->size(); ihit++){
			o_fitD->push_back(get_dist((*i_layerID)[ihit],(*i_wireID)[ihit],slX,inX,slZ,inZ));
		}
		t->Fill();

#ifdef DEBUG
#ifdef MC
		Double_t chi2r;
		//getchi2(chi2r,slirX,inirX,slirZ,inirZ,true);
		getchi2(chi2r,slirX,inirX,slirZ,inirZ);
#endif
		std::cout<<"*********************************************************"<<iEvent<<std::endl;
		std::cout<<"                  iEvent = "<<iEvent<<", triggerNumber = "<<triggerNumber<<std::endl;
		std::cout<<"Initial:"<<std::endl;
		std::cout<<"  Track: ("<<xup<<","<<yup<<","<<zup<<") --> ("<<xdown<<","<<ydown<<","<<zdown<<")"<<std::endl;
		std::cout<<"  para: ("<<sliX<<","<<iniX<<","<<sliZ<<","<<iniZ<<")"<<std::endl;
		std::cout<<"  Chi2 = "<<chi2i<<std::endl;
		std::cout<<"After Fitting:"<<std::endl;
		std::cout<<"  para: ("<<slX<<","<<inX<<","<<slZ<<","<<inZ<<")"<<std::endl;
		std::cout<<"  Chi2 = "<<chi2<<std::endl;
#ifdef MC
		std::cout<<"McTruth:"<<std::endl;
		std::cout<<"  para: ("<<slirX<<","<<inirX<<","<<slirZ<<","<<inirZ<<")"<<std::endl;
		std::cout<<"  Chi2 = "<<chi2r<<std::endl;
#endif
		std::cout<<"Hits:"<<std::endl;
		for (int ihit = 0; ihit<i_driftD->size(); ihit++){
			std::cout<<"  ("<<(*i_layerID)[ihit]<<","<<(*i_wireID)[ihit]<<"): "<<(*i_driftD)[ihit]<<"-"<<(*o_fitD)[ihit]<<"="<<(*i_driftD)[ihit]-(*o_fitD)[ihit]<<std::endl;
		}
#endif
	}
	t->Write();
	f->Close();

	return 0;
}

//______________________________________________________________________________
Double_t get_dist(int lid, int wid, Double_t slx, Double_t inx, Double_t slz, Double_t inz)
{
	double xdown = inx-slx*deltay;
	double zdown = inz-slz*deltay;
	vTrackU.SetXYZ(inx,yup,inz);
	vTrackD.SetXYZ(xdown,ydown,zdown);
	vWireHV.SetXYZ(Xhv[lid][wid],Yhv[lid][wid],zhv);
	vWireRO.SetXYZ(Xro[lid][wid],Yro[lid][wid],zro);
	vTrack = vTrackD-vTrackU;
	vWire = vWireRO-vWireHV;
	vDist = vWireHV-vTrackU;
	vAxis = vWire.Cross(vTrack);
	double value = vDist*(vAxis.Unit());
#ifdef DEBUG
	std::cout<<"  >> get_dist("<<lid<<","<<wid<<","<<slx<<","<<inx<<","<<slz<<","<<inz<<")"<<std::endl;
	std::cout<<"    >>> ("<<vTrackU.x()<<","<<vTrackU.y()<<","<<vTrackU.z()
		<<") --> ("
		<<vTrackD.x()<<","<<vTrackD.y()<<","<<vTrackD.z()
		<<")  to  ("
		<<vWireHV.x()<<","<<vWireHV.y()<<","<<vWireHV.z()
		<<") --> ("
		<<vWireRO.x()<<","<<vWireRO.y()<<","<<vWireRO.z()
		<<")  is "<<value<<std::endl;
	std::cout<<"    >>> Yhv["<<lid<<"]["<<wid<<"] = "<<Yhv[lid][wid]<<std::endl;
#endif
	return value;
}

//______________________________________________________________________________
void getchi2(Double_t &f, Double_t slx, Double_t inx, Double_t slz, Double_t inz,bool all)
{
	//calculate chisquare
	Double_t chisq = 0;
	Double_t delta;
	Double_t dfit;
	int N = i_driftD->size();

	for (Int_t i=0;i<N; i++) {
		if ((*i_layerID)[i]==thelayer&&!all) continue;
		dfit = get_dist((*i_layerID)[i],(*i_wireID)[i],slx,inx,slz,inz);
		delta  = (fabs((*i_driftD)[i])-fabs(dfit))/errord[(*i_layerID)[i]][(*i_wireID)[i]];
		chisq += delta*delta;
	}
	f = chisq;
}

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	getchi2(f,*par,*(par+1),*(par+2),*(par+3));
}

//______________________________________________________________________________
void do_fit(Double_t sliX, Double_t iniX,Double_t sliZ, Double_t iniZ){
	if(gMinuit) delete gMinuit;
	gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);
	arglist[0] = 0;
	gMinuit->SetPrintLevel(-1); // no print
	gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
//	gMinuit->mnparm(0, "slopeX", sliX, slerX/1.e3, sliX-slerX,sliX+slerX,ierflg);
//	gMinuit->mnparm(1, "interceptX", iniX, inerX/1.e3, iniX-inerX,iniX+fabs(inerX),ierflg);
//	gMinuit->mnparm(2, "slopeZ", sliZ, slerZ/1.e3, sliZ-slerZ,sliZ+slerZ,ierflg);
//	gMinuit->mnparm(3, "interceptZ", iniZ, inerZ/1.e3, iniZ-inerZ,iniZ+inerZ,ierflg);
	gMinuit->mnparm(0, "slopeX", sliX, slerX/1.e4, -slerX,slerX,ierflg);
	gMinuit->mnparm(1, "interceptX", iniX, inerX/1.e4, -inerX,inerX,ierflg);
	gMinuit->mnparm(2, "slopeZ", sliZ, slerZ/1.e4, -slerZ,slerZ,ierflg);
	gMinuit->mnparm(3, "interceptZ", iniZ, inerZ/1.e4, -inerZ,inerZ,ierflg);
//	std::cout<<"gMinuit->mnparm(0, \"slopeX\","<< sliX<<", "<<slerX/1.e2<<", "<<sliX-slerX<<","<<sliX+slerX<<","<<ierflg<<");"<<std::endl;
//	std::cout<<"gMinuit->mnparm(1, \"interceptX\","<< iniX<<", "<<inerX/1.e2<<", "<<iniX-inerX<<","<<iniX+fabs(inerX)<<","<<ierflg<<");"<<std::endl;
//	std::cout<<"gMinuit->mnparm(2, \"slopeZ\","<< sliZ<<", "<<slerZ/1.e2<<", "<<sliZ-slerZ<<","<<sliZ+slerZ<<","<<ierflg<<");"<<std::endl;
//	std::cout<<"gMinuit->mnparm(3, \"interceptZ\","<< iniZ<<", "<<inerZ/1.e2<<", "<<iniZ-inerZ<<","<<iniZ+inerZ<<","<<ierflg<<");"<<std::endl;

	// Now ready for minimization step
	arglist[0] = 500.0;
	arglist[1] = 1.0;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Print results
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//printf("====Rrestul====\n");
	//gMinuit->mnprin(3,amin);
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[suffix] [nEventMax] [iterationNo] [N_Min]\n",prog_name);
}
