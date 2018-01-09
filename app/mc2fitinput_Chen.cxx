#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom1.h"

#define NLAY 9
#define	NCEL 11

//===================Chamber Parameter============================
double chamberHL = 599.17/2; // mm
double chamberHH = 170.05/2; // mm
double chamberCY = 572; // mm
double sciYup = 0;
double sciYdown = 0;

// map for wire position
double  map_x[NLAY][NCEL][2];
double  map_y[NLAY][NCEL][2];

// for get dist
TVector3 vTrackU, vTrackD, vTrack;
TVector3 vWireHV, vWireRO, vWire;
TVector3 vDist;
TVector3 vAxis;

void printUsage(char* name);
int getwid(int layerID, int cellID);
double get_dist(int lid, int wid, double slx, double inx, double slz, double inz);

int main(int argc, char ** argv){

    if (argc<2){
        printUsage(argv[0]);
        return -1;
    }
    int geoSetup = 0;
    if (argc>=3) geoSetup = atoi(argv[2]);

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
        }
    }

    //===================Get Wire Position============================
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
            map_x[wp_lid][wp_wid][0] = wp_xhv;
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_x[wp_lid][wp_wid][1] = wp_xro;
            map_y[wp_lid][wp_wid][1] = wp_yro;
        }
    }
    TFile_wirepos->Close();

    //==============================Get input file==============================
    TChain * ichain = new TChain("tree","tree");
    ichain->Add(argv[1]);

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

    //==============================Prepare output file==============================
    TFile * ofile = new TFile("output.root","RECREATE");
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
	otree->Branch("driftDmc",&o_driftD);
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

	double m1x,m1y,m1z;
	double m2x,m2y,m2z;

	// Random for resolution
    TRandom1 random1;
    
    // counters
    bool checkup,checkdown;
    int nHitLayer[NLAY];

    int nEntries = ichain->GetEntries();
    for (int iEntry = 0; iEntry<nEntries; iEntry++){
        if (iEntry%1000==0) printf("iEntry %d\n",iEntry);
        ichain->GetEntry(iEntry);
        // check trigger
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

        // preset
        triggerNumber = iEntry;
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
        for (int i = 0; i<NLAY; i++){
            nHitLayer[i] = 0;
        }

        // get cdc hits
        for (int ihit = 0; ihit<CdcCell_tid->size(); ihit++){
        	// is this the wanted track?
            if ((*CdcCell_tid)[ihit]!=2) continue; // looking for e- from gamma conversion, tid==2
            // get cell ID
            int layerID = (*CdcCell_layerID)[ihit];
            int cellID = (*CdcCell_cellID)[ihit];
            int lid = layerID;
            int wid = getwid(layerID,cellID);
            // get driftD
            // FIXME: now we want to ignore scattering and corner effect, so recalculating driftD by DOCA!
            //double driftD = (*CdcCell_driftDtrue)[ihit]*10;
            double driftD = get_dist(lid,wid,slx,inx,slz,inz);
            // get driftT
            double driftT=fabs(driftD/0.023);
            // smear driftT
            driftT += random1.Gaus(0,7);

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

            // increment counters
            nHitLayer[lid]++;
            o_nHits++;
        }
		for (int lid = 0; lid < NLAY; lid++ ){
            if (nHitLayer[lid]) o_nLayers++;
        }

        otree->Fill();
    }
    otree->Write();
    ofile->Close();
    return 0;
}

void printUsage(char* name){
    printf("%s [inputFile]\n",name);
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
