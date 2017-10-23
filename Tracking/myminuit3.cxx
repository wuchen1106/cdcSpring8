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

#define NLAY    9
#define NCEL    11

#define NCAND   4

// ________About Cell___________
double map_yhv[NLAY][NCEL]; // High Voltage side
double map_xhv[NLAY][NCEL];
double map_yro[NLAY][NCEL]; // Read Out side
double map_xro[NLAY][NCEL];
double errord[NLAY][NCEL];
double chamberHL= 599.17/2;

// ________About Track___________
int theCand = 0;
std::vector<int> * i_layerID = 0;
std::vector<int> * i_wireID = 0;
std::vector<double> * o_driftD[NCAND] = {0};
std::vector<int> * o_sel[NCAND] = {0};
double yup = 623.97007;
double ydown = 527.60011;
double deltay = yup-ydown;
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;
TF1 * funcErr;

// ________About Fitting___________
TMinuit *gMinuit = 0;
double arglist[10];
int ierflg = 0;
double amin,edm,errdef;
int nvpar,nparx,icstat;
double slZmax = 0.2;
double slXmax = 0.1;
double inZmax = 20;
double inXmax = 20;

//______________________________________________________________________________
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz);
void getchi2(double &f, double slx, double inx, double slz, double inz,bool all = false);
void fcn(int &npar, double *gin, double &f, double *par, int iflag);
void do_fit(double slix, double inix,double sliz, double iniz);
void print_usage(char* prog_name);

//______________________________________________________________________________
int main(int argc, char** argv){
	if (argc<3){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int thelayer = (int)strtol(argv[2],NULL,10);
	std::string suffix = "";
	if (argc>=4){
		suffix = argv[3];
	}
	int iterationNo = 0;
	if (argc>=5) iterationNo = (int)strtol(argv[4],NULL,10);
	int nEventMax = 0;
	if (argc>=6) nEventMax = (int)strtol(argv[5],NULL,10);

	double eff_p2_average = 0;
	double eff_p3_average = 0;

	//===================Get error============================
	funcErr = new TF1("funcErr","0.346904-0.221775*x+0.080226*x*x-0.0128037*x*x*x+0.000755738*x*x*x*x",0,9);

	//===================Set Geometry============================
	// The z values
	// get the inital values
	TFile * if_geom = new TFile("../info/wire-position.root");
	TTree * t_geom = (TTree*) if_geom->Get("t");
	double x1_geom, y1_geom, z1_geom, x2_geom, y2_geom, z2_geom;
	int wid_geom,lid_geom;
	t_geom->SetBranchAddress("l",&lid_geom);
	t_geom->SetBranchAddress("w",&wid_geom);
	t_geom->SetBranchAddress("xhv",&x1_geom);
	t_geom->SetBranchAddress("yhv",&y1_geom);
	t_geom->SetBranchAddress("xro",&x2_geom);
	t_geom->SetBranchAddress("yro",&y2_geom);
	Float_t error = 0.2;
	for(int i = 0; i<t_geom->GetEntries(); i++){
		t_geom->GetEntry(i);
		map_xhv[lid_geom][wid_geom] = x1_geom;
		map_yhv[lid_geom][wid_geom] = y1_geom;
		map_xro[lid_geom][wid_geom] = x2_geom;
		map_yro[lid_geom][wid_geom] = y2_geom;
		errord[lid_geom][wid_geom]=error;
	}
	if_geom->Close();

	//===================Get ROOT File============================
	//TChain * c = new TChain("t","t");
	TChain * c = new TChain("t","t");
	std::stringstream buf;
	buf.str(""); buf.clear();
	buf<<"../root/i_"<<runNo<<".layer"<<thelayer<<"."<<suffix<<".root";
	c->Add(buf.str().c_str());
	std::cout<<"Adding \""<<buf.str()<<"\""<<std::endl;
	int triggerNumber;
    int i_nHits;
    int i_nHitsgood;
    int i_nHitLayers;
    std::vector<double> * i_driftT = 0;
    std::vector<int> * i_type = 0;
    std::vector<int> * i_np = 0;
    std::vector<int> * i_ip = 0;
    std::vector<int> * i_clk = 0;
    std::vector<int> * i_width = 0;
    std::vector<int> * i_peak = 0;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_aa = 0;
    std::vector<double> * i_dxl = 0;
    std::vector<double> * i_dxr = 0;
    std::vector<double> * i_calD[NCAND] = {0};
    int i_icombi[NCAND];
    int i_iselec[NCAND];
    int i_npairs[NCAND];
    double i_slx[NCAND];
    double i_inx[NCAND];
    double i_slz[NCAND];
    double i_inz[NCAND];
    double i_chi2x[NCAND];
    double i_chi2z[NCAND];
    c->SetBranchAddress("triggerNumber",&triggerNumber);
    c->SetBranchAddress("nHits",&i_nHits);
	c->SetBranchAddress("nHitsgood",&i_nHitsgood);
	c->SetBranchAddress("nHitLayers",&i_nHitLayers);
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
    c->SetBranchAddress("dxl",&i_dxl);
    c->SetBranchAddress("dxr",&i_dxr);
    for(int iCand = 0; iCand<NCAND; iCand++){
        c->SetBranchAddress(Form("calD%d",iCand),&(i_calD[iCand]));
    }
    c->SetBranchAddress(Form("icom[%d]",NCAND),i_icombi);
    c->SetBranchAddress(Form("isel[%d]",NCAND),i_iselec);
    c->SetBranchAddress(Form("npairs[%d]",NCAND),i_npairs);
    c->SetBranchAddress(Form("slx[%d]",NCAND),i_slx);
    c->SetBranchAddress(Form("inx[%d]",NCAND),i_inx);
    c->SetBranchAddress(Form("slz[%d]",NCAND),i_slz);
    c->SetBranchAddress(Form("inz[%d]",NCAND),i_inz);
    c->SetBranchAddress(Form("chi2x[%d]",NCAND),i_chi2x);
    c->SetBranchAddress(Form("chi2z[%d]",NCAND),i_chi2z);

	//===================Output file============================
	buf.str(""); buf.clear();
	//buf<<"../root/t_"<<runNo<<"."<<suffix<<iterationNo<<".root";
	buf<<"../root/t_"<<runNo<<".layer"<<thelayer<<"."<<suffix<<".root";
	TFile * f = new TFile(buf.str().c_str(),"RECREATE"); 
	TTree * t = new TTree("t","t");
    std::vector<double> * o_fitD[NCAND] = {0};
    int o_nHitsSel[NCAND];
    double o_slx[NCAND];
    double o_inx[NCAND];
    double o_slz[NCAND];
    double o_inz[NCAND];
    double o_chi2i[NCAND];
    double o_chi2[NCAND];

	t->Branch("triggerNumber",&triggerNumber);
	t->Branch("nHits",&i_nHits);
	t->Branch("nHitsgood",&i_nHitsgood);
	t->Branch("nHitLayers",&i_nHitLayers);
	t->Branch("driftT",&i_driftT);
	t->Branch("wireID",&i_wireID);
	t->Branch("layerID",&i_layerID);
    t->Branch("type",&i_type); // -1: dummy layer; 1: guard layer; 2: left end; 3: right end; 0: center cell;
    t->Branch("np",&i_np);
    t->Branch("ip",&i_ip);
    t->Branch("clk",&i_clk);
    t->Branch("width",&i_width);
	t->Branch("peak",&i_peak);
	t->Branch("sum",&i_sum);
    t->Branch("aa",&i_aa);
    t->Branch("dxl",&i_dxl);
    t->Branch("dxr",&i_dxr);
    for(int iCand = 0; iCand<NCAND; iCand++){
        t->Branch(Form("i_icom%d",iCand),&(i_icombi[iCand]));
        t->Branch(Form("i_isel%d",iCand),&(i_iselec[iCand]));
        t->Branch(Form("i_npairs%d",iCand),&(i_npairs[iCand]));
        t->Branch(Form("i_slx%d",iCand),&(i_slx[iCand]));
        t->Branch(Form("i_inx%d",iCand),&(i_inx[iCand]));
        t->Branch(Form("i_slz%d",iCand),&(i_slz[iCand]));
        t->Branch(Form("i_inz%d",iCand),&(i_inz[iCand]));
        t->Branch(Form("i_chi2x%d",iCand),&(i_chi2x[iCand]));
        t->Branch(Form("i_chi2z%d",iCand),&(i_chi2z[iCand]));
        t->Branch(Form("i_calD%d",iCand),&(i_calD[iCand]));
        t->Branch(Form("nHitsSel%d",iCand),&(o_nHitsSel[iCand]));
        t->Branch(Form("slx%d", iCand),&(o_slx[iCand]));
        t->Branch(Form("inx%d", iCand),&(o_inx[iCand]));
        t->Branch(Form("slz%d", iCand),&(o_slz[iCand]));
        t->Branch(Form("inz%d", iCand),&(o_inz[iCand]));
        t->Branch(Form("chi2i%d",iCand),&(o_chi2i[iCand]));
        t->Branch(Form("chi2%d",iCand),&(o_chi2[iCand]));
        t->Branch(Form("sel%d",iCand),&(o_sel[iCand]));
        t->Branch(Form("driftD%d",iCand),&(o_driftD[iCand]));
        t->Branch(Form("fitD%d",iCand),&(o_fitD[iCand]));
    }
    for (int iCand = 0; iCand<NCAND; iCand++){
        o_fitD[iCand] = new std::vector<double>;
        o_driftD[iCand] = new std::vector<double>;
        o_sel[iCand] = new std::vector<int>;
    }

    // copies for chi2 sorting
    int    ti_icombi[NCAND];
	int    ti_iselec[NCAND];
	int    ti_npairs[NCAND];
	double ti_slx[NCAND];
	double ti_inx[NCAND];
	double ti_slz[NCAND];
	double ti_inz[NCAND];
	double ti_chi2x[NCAND];
	double ti_chi2z[NCAND];
	std::vector<double> * ti_calD[NCAND];
	for (int iCand = 0; iCand<NCAND; iCand++){
	    ti_calD[iCand] = new std::vector<double>;
    }
    int    to_nHitsSel[NCAND];
	double to_slx[NCAND];
	double to_inx[NCAND];
	double to_slz[NCAND];
	double to_inz[NCAND];
	double to_chi2[NCAND];
	int sortCand[NCAND];

	//===================Loop in Events============================
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	std::cout<<"Processing "<<N<<" Events ..."<<std::endl;
	for (Long64_t iEvent = 0;iEvent<N; iEvent++){
		if (iEvent%1000==0) std::cout<<(double)iEvent/N*100<<"..."<<std::endl;
		c->GetEntry(iEvent);

		// pre settings
        for(int iCand = 0; iCand<NCAND; iCand++){
            o_driftD[iCand]->clear();
            o_fitD[iCand]->clear();
            o_sel[iCand]->clear();
            o_nHitsSel[iCand] = 0;
            o_slx[iCand] = 0;
            o_inx[iCand] = 0;
            o_slz[iCand] = 0;
            o_inz[iCand] = 0;
            o_chi2i[iCand] = 0;
            o_chi2[iCand] = 0;
            sortCand[iCand] = iCand;
            ti_calD[iCand]->clear();
        }

        for (theCand = 0; theCand<NCAND; theCand++){
            // get hit list
            for (int ihit = 0; ihit<i_nHits; ihit++){
                double fitd = (*i_calD[theCand])[ihit];
                double dd;
                if (fitd>0) dd = (*i_dxr)[ihit];
                else dd = (*i_dxl)[ihit];
                int selected = 0;
                if (fabs(fitd-dd)<2&&thelayer!=(*i_layerID)[ihit]){
                    selected = 1;
                    o_nHitsSel[theCand]++;
                }
                o_driftD[theCand]->push_back(dd);
                o_sel[theCand]->push_back(selected);
            }
            if (o_nHitsSel[theCand]<5){
//                printf("Cand[%d]: nhits selected is %d!\n",theCand,o_nHitsSel[theCand]);
                continue;
            }
            // do the fitting
            do_fit(i_slx[theCand],i_inx[theCand],i_slz[theCand],i_inz[theCand]);
            double temp;
            gMinuit->GetParameter(0, o_slx[theCand], temp);
            gMinuit->GetParameter(1, o_inx[theCand], temp);
            gMinuit->GetParameter(2, o_slz[theCand], temp);
            gMinuit->GetParameter(3, o_inz[theCand], temp);
            getchi2(o_chi2[theCand],o_slx[theCand],o_inx[theCand],o_slz[theCand],o_inz[theCand]);
        }
        // copy candidates
        for (int iCand = 0; iCand<NCAND; iCand++){
            ti_icombi[iCand] = i_icombi[iCand];
            ti_iselec[iCand] = i_iselec[iCand];
            ti_npairs[iCand] = i_npairs[iCand];
            ti_slx[iCand] = i_slx[iCand];
            ti_inx[iCand] = i_inx[iCand];
            ti_slz[iCand] = i_slz[iCand];
            ti_inz[iCand] = i_inz[iCand];
            ti_chi2x[iCand] = i_chi2x[iCand];
            ti_chi2z[iCand] = i_chi2z[iCand];
            to_nHitsSel[iCand] = o_nHitsSel[iCand];
            to_slx[iCand] = o_slx[iCand];
            to_inx[iCand] = o_inx[iCand];
            to_slz[iCand] = o_slz[iCand];
            to_inz[iCand] = o_inz[iCand];
            to_chi2[iCand] = o_chi2[iCand];
            for (int ihit = 0; ihit<i_nHits; ihit++){
                ti_calD[iCand]->push_back((*i_calD[iCand])[ihit]);
            }
        }
        // sort candidates
        for (int iSort = 0; iSort<NCAND-1; iSort++){
            for (int jSort = iSort+1; jSort<NCAND; jSort++){
                if (o_nHitsSel[sortCand[iSort]]<o_nHitsSel[sortCand[jSort]]
                 ||(o_nHitsSel[sortCand[iSort]]==o_nHitsSel[sortCand[jSort]]&&o_chi2[sortCand[iSort]]>o_chi2[sortCand[jSort]])
                   ){
                    int tCand = sortCand[iSort]; sortCand[iSort] = sortCand[jSort]; sortCand[jSort] = tCand;
                }
            }
        }
//        printf("#%d rank: %d %d %d %d\n",iEvent,sortCand[0],sortCand[1],sortCand[2],sortCand[3]);
        for (int iSort = 0; iSort<NCAND; iSort++){
            int iCand = sortCand[iSort];
            i_icombi[iSort] = ti_icombi[iCand];
            i_iselec[iSort] = ti_iselec[iCand];
            i_npairs[iSort] = ti_npairs[iCand];
            i_slx[iSort] = ti_slx[iCand];
            i_inx[iSort] = ti_inx[iCand];
            i_slz[iSort] = ti_slz[iCand];
            i_inz[iSort] = ti_inz[iCand];
            i_chi2x[iSort] = ti_chi2x[iCand];
            i_chi2z[iSort] = ti_chi2z[iCand];
            o_nHitsSel[iSort] = to_nHitsSel[iCand];
            o_slx[iSort] = to_slx[iCand];
            o_inx[iSort] = to_inx[iCand];
            o_slz[iSort] = to_slz[iCand];
            o_inz[iSort] = to_inz[iCand];
            o_chi2[iSort] = to_chi2[iCand];
            for (int ihit = 0; ihit<i_nHits; ihit++){
                (*i_calD[iSort])[ihit] = (*ti_calD[iCand])[ihit];
            }
            // set output
            getchi2(o_chi2i[iSort],i_slx[iSort],i_inx[iSort],i_slz[iSort],i_inz[iSort]);
            for (int ihit = 0; ihit<i_nHits; ihit++){
                o_fitD[iSort]->push_back(get_dist((*i_layerID)[ihit],(*i_wireID)[ihit],o_slx[iSort],o_inx[iSort],o_slz[iSort],o_inz[iSort]));
            }
        }
		t->Fill();
	}
	t->Write();
	f->Close();

	return 0;
}

//______________________________________________________________________________
void do_fit(double sliX, double iniX,double sliZ, double iniZ){
	if(gMinuit) delete gMinuit;
	gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);
	arglist[0] = 0;
	gMinuit->SetPrintLevel(-1); // no print
	gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	// Set starting values and step sizes for parameters
	gMinuit->mnparm(0, "slopeX", sliX, slXmax/1.e4, -slXmax,slXmax,ierflg);
	gMinuit->mnparm(1, "interceptX", iniX, inXmax/1.e4, -inXmax,inXmax,ierflg);
	gMinuit->mnparm(2, "slopeZ", sliZ, slZmax/1.e4, -slZmax,slZmax,ierflg);
	gMinuit->mnparm(3, "interceptZ", iniZ, inZmax/1.e4, -inZmax,inZmax,ierflg);

	// Now ready for minimization step
	arglist[0] = 500.0;
	arglist[1] = 1.0;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	// Print results
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//printf("====Rrestul====\n");
	//gMinuit->mnprin(3,amin);
}

//______________________________________________________________________________
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
	getchi2(f,*par,*(par+1),*(par+2),*(par+3));
}

//______________________________________________________________________________
void getchi2(double &f, double slx, double inx, double slz, double inz,bool all)
{
	//calculate chisquare
	double chisq = 0;
	double delta;
	double dfit;
	int N = o_driftD[theCand]->size();

	for (int i=0;i<N; i++) {
		if ((*o_sel[theCand])[i]==0) continue;
		dfit = get_dist((*i_layerID)[i],(*i_wireID)[i],slx,inx,slz,inz);
		// FIXME: we should consider about the error
//		double error = errord[(*i_layerID)[i]][(*i_wireID)[i]]*(fabs(fabs(dfit)-4)+2/1.5)*1.5/4;
//		double error = errord[(*i_layerID)[i]][(*i_wireID)[i]];
//		double error = funcErr->Eval(fabs(dfit));
		double error = 0.2;
        delta  = ((*o_driftD[theCand])[i]-dfit)/error;
		chisq += delta*delta;
	}
	f = chisq;
}

//______________________________________________________________________________
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz)
{
	double xdown = inx-slx*deltay;
	double zdown = inz-slz*deltay;
	vTrackU.SetXYZ(inx,yup,inz);
	vTrackD.SetXYZ(xdown,ydown,zdown);
	vWireHV.SetXYZ(map_xhv[lid][wid],map_yhv[lid][wid],-chamberHL);
	vWireRO.SetXYZ(map_xro[lid][wid],map_yro[lid][wid],chamberHL);
	vTrack = vTrackD-vTrackU;
	vWire = vWireRO-vWireHV;
	vDist = vWireHV-vTrackU;
	vAxis = vWire.Cross(vTrack);
	double value = -vDist*(vAxis.Unit());
	return value;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] [thelayer] <[suffix] [nEventMax] [iterationNo]\n",prog_name);
}
