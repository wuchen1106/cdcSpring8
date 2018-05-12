#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom1.h"
#include "TF1.h"
#include "TGraph.h"

#include "header.h"

//===================Chamber Parameter============================
double chamberHL = 599.17/2; // mm
double chamberHH = 170.05/2; // mm
double chamberCY = 572; // mm
double sciYup = 0;
double sciYdown = 0;

// map for wire position
double  map_x[NLAY][NCEL][2];
double  map_y[NLAY][NCEL][2];
bool    map_exist[NLAY][NCEL];

// for get dist
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

// for smearing
int smearType = 0;
TRandom1 random1;

void printUsage(char* name);
int getwid(int layerID, int cellID);
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz);
double getResT(double x);
double getResX(double x);
double t2x(double t, bool isRight=true);
double x2t(double x);
void smearXT(double driftDmc,double & driftD, double & driftT);

TF1 * f_xt_l = 0;
TF1 * f_xt_r = 0;
TF1 * f_sigT = 0;
TGraph * gr_sigX = 0;

int main(int argc, char ** argv){

    if (argc<3){
        printUsage(argv[0]);
        return -1;
    }
    TString xt_filename = argv[2];
    int geoSetup = 0;
    TString outputFileName = "output.root";
    if (argc>=4) outputFileName = argv[3];
    if (argc>=5) geoSetup = atoi(argv[4]);
    double offset_sigma = 0.05;
    if (argc>=6) offset_sigma = atof(argv[5]);
    double offset_max = offset_sigma*2;
    if (argc>=7) offset_max = atof(argv[6]);
    if (offset_max<0){
        printf("offset_max = %.3e < 0! Will set positive\n",offset_max);
        offset_max *= -1;
    }
    int maxLayer = 0;
    if (argc>=8) maxLayer = atoi(argv[7]);
    if (maxLayer>NLAY-1){
    	maxLayer=NLAY-1;
    	printf("WARNING: maxLayer[%d] is exceeding the range! Totally support %d layers counting from 0\n",maxLayer,NLAY);
	}
	smearType = 0; // default: smear driftT
    if (argc>=9) smearType = atoi(argv[8]);
    int inputFileType = 0; // 0: using simulation file as input; 1: using tracking file (ana_XXX) as input
    if (argc>=10) inputFileType = atoi(argv[9]);
    int iEntryStart = 0;
    int iEntryStop = 0;
    if (argc>=12){
        iEntryStart = (int)strtol(argv[10],NULL,10);
        iEntryStop = (int)strtol(argv[11],NULL,10);
    }
    double chi2max = 2;
    if (argc>=13) chi2max = atof(argv[12]);;
    int nHitsSMin = 7;
    if (argc>=14) nHitsSMin = (int)strtol(argv[13],NULL,10);
    int nHitsMax = 0;
    if (argc>=15) nHitsMax = (int)strtol(argv[14],NULL,10);
    double slzMax = 0.1;
    if (argc>=16) slzMax = atof(argv[15]);
    printf("##############%s with %d Parameters##################\n",argv[0],argc);
	printf("geoSetup = %d\n",geoSetup);
	printf("offSig   = %.3e\n",offset_sigma);
	printf("offMax   = %.3e\n",offset_max);
	printf("maxLayer = %d\n",maxLayer);
	printf("smear on \"%s\"\n",smearType==0?"T":"X");
	printf("inputFile: %s, \"%s\"\n",inputFileType==0?"MC":"Data",argv[1]);
	printf("Entries: %d~%d (%s~%s)\n",iEntryStart,iEntryStop,argv[9],argv[10]);
	printf("max chi2:  %.3e\n",chi2max);
	printf("min nHitsS:%d\n",nHitsSMin);
	printf("max nHitsG:%d\n",nHitsMax);
	printf("max slz:   %.3e\n",slzMax);

    TString HOME=getenv("CDCS8WORKING_DIR");

	//===================Geometry Parameter============================
	if (geoSetup==0){
		// normal scintillator
		sciYup = chamberCY+chamberHH+180; // mm
        sciYdown = chamberCY-chamberHH-180; 
	}
	else{
		// finger scintillator
		sciYup = chamberCY+chamberHH+250; // mm
        sciYdown = chamberCY-chamberHH-195; 
	}

    //===================Prepare Maps============================
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            map_x[lid][wid][0] = 0;
            map_y[lid][wid][0] = 0;
            map_x[lid][wid][1] = 0;
            map_y[lid][wid][1] = 0;
            map_exist[lid][wid] = false;
        }
    }

    //===================Get XT and error=============================
	TFile * i_xt_file = new TFile(xt_filename);
	f_xt_l = (TF1*) i_xt_file->Get("fl_0");
	f_xt_r = (TF1*) i_xt_file->Get("fr_0");
	if (smearType==0){
        f_sigT = (TF1*) i_xt_file->Get("f_sigT_0");
    }
    else{
        gr_sigX = (TGraph*) i_xt_file->Get("gr_resIni");
        if (!gr_sigX){
            fprintf(stderr,"WARNING: Cannot find gr_resIni!\n");
            return 1;
        }
    }

    //===================Get Wire Position============================
    // prepare offset
    double deltaX[NLAY][NCEL];
    if (offset_sigma){
        for ( int lid = 0; lid<NLAY; lid++){
            for ( int wid = 0; wid<NCEL; wid++){
                if (lid<=8&&wid>=11) continue; // to be consistant with previous runs
                //deltaX[lid][wid] = random1.Gaus(0,0.1);//
                double dx = random1.Gaus(0,offset_sigma);
                while(offset_max&&fabs(dx)>offset_max){
                    dx = random1.Gaus(0,offset_sigma);
                }
                deltaX[lid][wid] = dx;
                printf("%d %d %.3e\n",lid,wid,deltaX[lid][wid]);
            }
        }
    }

    TFile * TFile_wirepos = new TFile(HOME+"/info/wire-position.root");
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
            map_x[wp_lid][wp_wid][0] = wp_xhv+deltaX[wp_lid][wp_wid];
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_x[wp_lid][wp_wid][1] = wp_xro+deltaX[wp_lid][wp_wid];
            map_y[wp_lid][wp_wid][1] = wp_yro;
            map_exist[wp_lid][wp_wid] = true;
        }
    }
    TFile_wirepos->Close();

    //==============================Get input file==============================
    TChain * ichain = 0;

    // trigger
    std::vector<double> * M_x = 0;
    std::vector<double> * M_y = 0;
    std::vector<double> * M_z = 0;
    std::vector<int>    * M_vid = 0;
    std::vector<int>    * M_tid = 0;
    // cdc
    std::vector<int>    * CdcCell_tid = 0;
    std::vector<int>    * CdcCell_layerID = 0;
    std::vector<int>    * CdcCell_cellID = 0;
    std::vector<double> * CdcCell_driftDtrue = 0;
    std::vector<double> * CdcCell_x = 0;
    std::vector<double> * CdcCell_y = 0;
    std::vector<double> * CdcCell_z = 0;
    // wire
    std::vector<double> * wire_x = 0;
    std::vector<double> * wire_y = 0;
    std::vector<double> * wire_z = 0;
    std::vector<int>    * wire_tid = 0;
    // truth
    std::vector<int>    * McTruth_tid = 0;
    std::vector<double> * McTruth_x = 0;
    std::vector<double> * McTruth_y = 0;
    std::vector<double> * McTruth_z = 0;
    std::vector<double> * McTruth_px = 0;
    std::vector<double> * McTruth_py = 0;
    std::vector<double> * McTruth_pz = 0;
    // Data input
    double i_chi2;
    int    i_nHitsS;
    int    i_nHitsG;
    double i_slx;
    double i_slz;
    double i_inx;
    double i_inz;
    int i_triggerNumber;

	if (inputFileType==0){
        ichain=new TChain("tree","tree");
        ichain->Add(argv[1]);
		ichain->SetBranchAddress("M_x",&M_x);
		ichain->SetBranchAddress("M_y",&M_y);
		ichain->SetBranchAddress("M_z",&M_z);
		ichain->SetBranchAddress("M_volID",&M_vid);
		ichain->SetBranchAddress("M_tid",&M_tid);
		ichain->SetBranchAddress("CdcCell_layerID",&CdcCell_layerID);
		ichain->SetBranchAddress("CdcCell_cellID",&CdcCell_cellID);
		ichain->SetBranchAddress("CdcCell_driftDtrue",&CdcCell_driftDtrue);
		ichain->SetBranchAddress("CdcCell_tid",&CdcCell_tid);
		ichain->SetBranchAddress("CdcCell_x",&CdcCell_x);
		ichain->SetBranchAddress("CdcCell_y",&CdcCell_y);
		ichain->SetBranchAddress("CdcCell_z",&CdcCell_z);
		ichain->SetBranchAddress("wire_tid",&wire_tid);
		ichain->SetBranchAddress("wire_x",&wire_x);
		ichain->SetBranchAddress("wire_y",&wire_y);
		ichain->SetBranchAddress("wire_z",&wire_z);
		ichain->SetBranchAddress("McTruth_tid",&McTruth_tid);
		ichain->SetBranchAddress("McTruth_x",&McTruth_x);
		ichain->SetBranchAddress("McTruth_y",&McTruth_y);
		ichain->SetBranchAddress("McTruth_z",&McTruth_z);
		ichain->SetBranchAddress("McTruth_px",&McTruth_px);
		ichain->SetBranchAddress("McTruth_py",&McTruth_py);
		ichain->SetBranchAddress("McTruth_pz",&McTruth_pz);
	}
	else{
        ichain=new TChain("t","t");
        ichain->Add(argv[1]);
		ichain->SetBranchAddress("chi2",&i_chi2);
		ichain->SetBranchAddress("nHitsS",&i_nHitsS);
		ichain->SetBranchAddress("nHitsG",&i_nHitsG);
		ichain->SetBranchAddress("slz",&i_slz);
		ichain->SetBranchAddress("slx",&i_slx);
		ichain->SetBranchAddress("inz",&i_inz);
		ichain->SetBranchAddress("inx",&i_inx);
		ichain->SetBranchAddress("triggerNumber",&i_triggerNumber);
	}

    //==============================Prepare output file==============================
    TFile * ofile = new TFile(outputFileName,"RECREATE");
    TTree * otree = new TTree("t","t");
	int triggerNumber;
	int o_nHits;
	int o_nLayers;
	std::vector<int> * o_layerID = 0;
	std::vector<int> * o_wireID = 0;
	std::vector<int> * o_type = 0;
	std::vector<int> * o_np = 0;
	std::vector<int> * o_ip = 0;
	std::vector<int> * o_clk = 0;
	std::vector<int> * o_width = 0;
	std::vector<int> * o_peak = 0;
	std::vector<int> * o_height = 0;
	std::vector<int> * o_mpn = 0;
	std::vector<int> * o_mpi = 0;
	std::vector<int> * o_rank = 0;
	std::vector<double> * o_ped = 0;
	std::vector<double> * o_sum = 0;
	std::vector<double> * o_aa = 0;
	std::vector<double> * o_driftT = 0;
	std::vector<double> * o_driftD = 0;
	std::vector<double> * o_driftDmc = 0;
	double inx,inz;
	double slx,slz;
	otree->Branch("triggerNumber",&triggerNumber);
	otree->Branch("nHits",&o_nHits);
	otree->Branch("nLayers",&o_nLayers);
	otree->Branch("layerID",&o_layerID);
	otree->Branch("wireID",&o_wireID);
	otree->Branch("type",&o_type);
	otree->Branch("np",&o_np);
	otree->Branch("ip",&o_ip);
	otree->Branch("clk",&o_clk);
	otree->Branch("width",&o_width);
	otree->Branch("peak",&o_peak);
	otree->Branch("height",&o_height);
	otree->Branch("mpn",&o_mpn);
	otree->Branch("mpi",&o_mpi);
	otree->Branch("rank",&o_rank);
	otree->Branch("ped",&o_ped);
	otree->Branch("sum",&o_sum);
	otree->Branch("aa",&o_aa);
	otree->Branch("driftT",&o_driftT);
	otree->Branch("driftDmc",&o_driftDmc);
	otree->Branch("driftD",&o_driftD);
	otree->Branch("slxmc",&slx);
	otree->Branch("slzmc",&slz);
	otree->Branch("inxmc",&inx);
	otree->Branch("inzmc",&inz);
	o_layerID = new std::vector<int>;
	o_wireID = new std::vector<int>;
	o_type = new std::vector<int>;
	o_ip = new std::vector<int>;
	o_np = new std::vector<int>;
	o_clk = new std::vector<int>;
	o_width = new std::vector<int>;
	o_peak = new std::vector<int>;
	o_height = new std::vector<int>;
	o_mpn = new std::vector<int>;
	o_mpi = new std::vector<int>;
	o_rank = new std::vector<int>;
	o_ped = new std::vector<double>;
	o_sum = new std::vector<double>;
	o_aa = new std::vector<double>;
	o_driftT = new std::vector<double>;
	o_driftD = new std::vector<double>;
	o_driftDmc = new std::vector<double>;

	double m1x,m1y,m1z;
	double m2x,m2y,m2z;
    
    // counters
    bool checkup,checkdown;
    int nHitLayer[NLAY];

    Long64_t nEntries = ichain->GetEntries();
    if (!iEntryStop&&!iEntryStart){iEntryStart = 0; iEntryStop=nEntries-1;}
    for (int iEntry = iEntryStart; iEntry<=iEntryStop; iEntry++){
        if (iEntry%10000==0) printf("iEntry %d\n",iEntry);
        ichain->GetEntry(iEntry);
        double slx,slz,inx,inz;
        if (inputFileType==0){ // MC: check trigger
			checkdown = false;
			checkup = false;
			for (int i = 0; i<M_x->size(); i++){
				if ((*M_tid)[i]!=2) continue;
				if ((*M_y)[i]>100) continue;
				if ((*M_vid)[i]==0){
					m1y = (*M_y)[i]*10+chamberCY;
					m1z = (*M_z)[i]*10;
					m1x = (*M_x)[i]*10;
					checkup = true;
				}
				else if ((*M_vid)[i]==1){
					m2y = (*M_y)[i]*10+chamberCY;
					m2z = (*M_z)[i]*10;
					m2x = (*M_x)[i]*10;
					checkdown = true;
				}
			}
			if (!checkdown||!checkup) continue;
			// get initial value
			slx = m2y-m1y==0?0:(m2x-m1x)/(m2y-m1y);
			slz = m2y-m1y==0?0:(m2z-m1z)/(m2y-m1y);
			inx = slx*(sciYup-m1y)+m1x;
			inz = slz*(sciYup-m1y)+m1z;
			triggerNumber = iEntry;
		}
		else{ // Data: check fitting quality
			if (i_chi2>chi2max) continue;
			if (i_nHitsS<nHitsSMin) continue;
			if (nHitsMax&&i_nHitsG>nHitsMax) continue;
			if (fabs(i_slz)>slzMax) continue;
			slx = i_slx;
			slz = i_slz;
			inx = i_inx;
			inz = i_inz;
			triggerNumber = i_triggerNumber;
		}

        // preset
        o_nHits = 0;
        o_nLayers = 0;
		o_layerID->clear();
		o_wireID->clear();
		o_type->clear();
		o_ip->clear();
		o_np->clear();
		o_clk->clear();
		o_width->clear();
		o_peak->clear();
		o_height->clear();
		o_mpn->clear();
		o_mpi->clear();
		o_rank->clear();
		o_ped->clear();
		o_sum->clear();
		o_aa->clear();
		o_driftT->clear();
		o_driftD->clear();
		o_driftDmc->clear();
        for (int i = 0; i<NLAY; i++){
            nHitLayer[i] = 0;
        }

        // get cdc hits
        if (!maxLayer&&inputFileType==0){
			for (int ihit = 0; ihit<CdcCell_tid->size(); ihit++){
				// is this the wanted track?
				if ((*CdcCell_tid)[ihit]!=2) continue; // looking for e- from gamma conversion, tid==2
				// get cell ID
				int lid = (*CdcCell_layerID)[ihit];
				int cellID = (*CdcCell_cellID)[ihit];
				int wid = getwid(lid,cellID);
				// get driftDmc
				// FIXME: now we want to ignore scattering and corner effect, so recalculating driftDmc by DOCA!
				//double driftDmc = (*CdcCell_driftDtrue)[ihit]*10;
				double driftDmc = get_dist(lid,wid,slx,inx,slz,inz);
				double driftT = 1e9;
				double driftD = 1e9;
				smearXT(driftDmc,driftD,driftT);

				// set output
				o_layerID->push_back(lid);
				o_wireID->push_back(wid);
				o_type->push_back(0);
				o_ip->push_back(0);
				o_np->push_back(1);
				o_clk->push_back(0);
				o_width->push_back(1);
				o_peak->push_back(300);
				o_height->push_back(300);
				o_mpn->push_back(1);
				o_mpi->push_back(0);
				o_rank->push_back(0);
				o_ped->push_back(200);
				o_sum->push_back(100);
				o_aa->push_back(100);
				o_driftT->push_back(driftT);
				o_driftD->push_back(driftD);
				o_driftDmc->push_back(driftDmc);

				// increment counters
				nHitLayer[lid]++;
				o_nHits++;
			}
		}
		else{
			for (int lid = 1; lid<=maxLayer; lid++){
				double driftDmc = 1e9;
				int wid = -1;
				for (int iw = 1; iw<NCEL; iw++){
					double doca = get_dist(lid,iw,slx,inx,slz,inz);
					if (fabs(doca)<10&&fabs(doca)<fabs(driftDmc)){
						driftDmc = doca;
						wid = iw;
					}
				}
				if (wid==-1) continue;
				// get driftT
				double driftT = 1e9;
				double driftD = 1e9;
				smearXT(driftDmc,driftD,driftT);

				// set output
				o_layerID->push_back(lid);
				o_wireID->push_back(wid);
				o_type->push_back(0);
				o_ip->push_back(0);
				o_np->push_back(1);
				o_clk->push_back(0);
				o_width->push_back(1);
				o_peak->push_back(300);
				o_height->push_back(300);
				o_mpn->push_back(1);
				o_mpi->push_back(0);
				o_rank->push_back(0);
				o_ped->push_back(200);
				o_sum->push_back(100);
				o_aa->push_back(100);
				o_driftT->push_back(driftT);
				o_driftD->push_back(driftD);
				o_driftDmc->push_back(driftDmc);

				// increment counters
				nHitLayer[lid]++;
				o_nHits++;
			}
			for (int lid = 0; lid < NLAY; lid++ ){
				if (nHitLayer[lid]) o_nLayers++;
			}
		}

        otree->Fill();
    }
    otree->Write();
    ofile->Close();
    return 0;
}

void printUsage(char* name){
    printf("%s [inputFile] [xtfile] <[outputfile] [geoSetup (0):general; 1: finger] [sigma_posX (0.05 mm)] [sigma_posXmax (2*sigma)] [maxLayer (0):same as MC; n] [SmearXorT (0):t; 1:x]  [inputFileType: (0): MC, 1: Data]>\n",name);
}

int getwid(int layerID, int cellID){
    int wid = 0;
    if (layerID==0)      wid = cellID-189;
    else if (layerID==1) wid = cellID>=194?cellID-194:cellID+4;
    else if (layerID==2) wid = cellID-194;
    else if (layerID==3) wid = cellID>=206?cellID-206:cellID+4;
    else if (layerID==4) wid = cellID-205;
    else if (layerID==5) wid = cellID>=217?cellID-217:cellID+5;
    else if (layerID==6) wid = cellID-217;
    else if (layerID==7) wid = cellID>=228?cellID-228:cellID+6;
    else if (layerID==8) wid = cellID-229;
    return wid;
}

double get_dist(int lid, int wid, double slx, double inx, double slz, double inz)
{
	if (!map_exist[lid][wid]) return 1e9;
	double xdown = inx-slx*(sciYup-sciYdown);
	double zdown = inz-slz*(sciYup-sciYdown);
	vTrackU.SetXYZ(inx,sciYup,inz);
	vTrackD.SetXYZ(xdown,sciYdown,zdown);
	vWireHV.SetXYZ(map_x[lid][wid][0],map_y[lid][wid][0],-chamberHL);
	vWireRO.SetXYZ(map_x[lid][wid][1],map_y[lid][wid][1],chamberHL);
	vTrack = vTrackD-vTrackU;
	vWire = vWireRO-vWireHV;
	vDist = vWireHV-vTrackU;
	vAxis = vWire.Cross(vTrack);
	double value = -vDist*(vAxis.Unit());
	return value;
}

double x2t(double x){
	TF1 * f = f_xt_r;
	if (x<0) f = f_xt_l;
	double t;
	if (fabs(x)<1)
		t = f->GetX(x,-25,100);
	else
		t = f->GetX(x,-25,800);
}

double t2x(double t, bool isRight){
	TF1 * f = f_xt_r;
	if (!isRight) f = f_xt_l;
	double x = f->Eval(t);
	return x;
}

double getResT(double x){
	return f_sigT->Eval(fabs(x));
}

double getResX(double x){
//    printf("x = %.3e\n",x);
    x = fabs(x);
	double residual = 0;

    bool foundPre = false;
    bool foundPost = false;
    double preX = 0;
    double postX = 0;
    double preRes = 0;
    double postRes = 0;
    int nPoints = gr_sigX->GetN();
    for (int i = 0; i<nPoints; i++){
        double res, doca;
        gr_sigX->GetPoint(i,doca,res);
        if (doca>x){
            postX = doca;
            postRes = res;
            foundPost = true;
            break;
        }
        preX = doca;
        preRes = res;
        foundPre = true;
    }
    if (!foundPre){
        if (foundPost){
            residual = postRes;
        }
        else{
            residual = 0;
        }
    }
    else if (!foundPost){
        residual = preRes;
    }
    else{
        residual = (postRes*(x-preX)+preRes*(postX-x))/(postX-preX);
    }
//    printf("found pre?%s: %.3e %.3e, found post?%s: %.3e %.3e, => %.3e\n",foundPre?"yes":"no",preX,preRes,foundPost?"yes":"no",postX,postRes,residual);

	double res = residual;
    return res;
}

void smearXT(double driftDmc,double & driftD, double & driftT){
    if (smearType==0){// smear driftT
        driftT = x2t(driftDmc);
        double resT = getResT(driftDmc);
        driftT += random1.Gaus((resT-4)*3,resT);
        driftD = t2x(driftT,driftDmc>0);
    }
    else{// smear driftD if needed
        double resX = getResX(driftDmc);
//        driftD = fabs(driftDmc)+random1.Gaus(0,resX);
//        if (driftD<0) driftD = 0; // FIXME: don't smear driftD across 0
//        if (driftDmc<0) driftD*=-1; // keep the left/right info
        driftD = driftDmc+random1.Gaus(0,resX);
        driftT = x2t(driftD);
    }
}
