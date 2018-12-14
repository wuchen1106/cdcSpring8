#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom1.h"

#include "header.hxx"

//===================Chamber Parameter============================
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

    double m1x, m1y, m1z;
    double m2x, m2y, m2z;
    double ox,oy,oz;
    double opx,opy,opz;
    int nHit;
    int triggerNumber;
    std::vector<double> * o_driftD = 0;
    std::vector<double> * o_driftT = 0;
    std::vector<int>    * o_layerID = 0;
    std::vector<int>    * o_cellID = 0;
    std::vector<int>    * o_nHitLayer = 0;
    std::vector<int>    * o_nPeaks = 0;
    std::vector<int>    * o_peakWidth = 0;
    std::vector<double> * o_adcPeak = 0;
    std::vector<double> * o_q = 0;
    std::vector<double> * o_t0 = 0;
    std::vector<double> * o_ped = 0;
    std::vector<bool> * o_goodHitFlag = 0;
    std::vector<double> * o_adccut = 0;
    std::vector<bool> * o_choosePeakFlag = 0;
    std::vector<std::vector<double> > * o_driftTAll = 0;
    otree->Branch("m1x",&m1x);
    otree->Branch("m1y",&m1y);
    otree->Branch("m1z",&m1z);
    otree->Branch("m2x",&m2x);
    otree->Branch("m2y",&m2y);
    otree->Branch("m2z",&m2z);
    otree->Branch("ox",&ox);
    otree->Branch("oy",&oy);
    otree->Branch("oz",&oz);
    otree->Branch("opx",&opx);
    otree->Branch("opy",&opy);
    otree->Branch("opz",&opz);
    otree->Branch("triggerNumber",&triggerNumber);
    otree->Branch("nHit",&nHit);
    otree->Branch("nHitLayer",&o_nHitLayer);
    otree->Branch("nPeaks",&o_nPeaks);
    otree->Branch("peakWidth",&o_peakWidth);
    otree->Branch("driftTime",&o_driftT);
    otree->Branch("driftTime2",&o_driftT);
    otree->Branch("driftTime3",&o_driftT);
    otree->Branch("driftD",&o_driftD);
    otree->Branch("layerID",&o_layerID);
    otree->Branch("cellID",&o_cellID);
    otree->Branch("adcPeak",&o_adcPeak);
    otree->Branch("q",&o_q);
    otree->Branch("choosePeakFlag",&o_choosePeakFlag);
    otree->Branch("goodHitFlag",&o_goodHitFlag);
    otree->Branch("ped",&o_ped);
    otree->Branch("adccut",&o_adccut);
    otree->Branch("t0offset",&o_t0);
    otree->Branch("t0",&o_t0);
    otree->Branch("driftTAll",&o_driftTAll);
    o_driftD    = new std::vector<double>;
    o_driftT    = new std::vector<double>;
    o_layerID   = new std::vector<int>;
    o_cellID    = new std::vector<int>;
    o_nHitLayer = new std::vector<int>;
    o_nPeaks    = new std::vector<int>;
    o_peakWidth = new std::vector<int>;
    o_adcPeak   = new std::vector<double>;
    o_q         = new std::vector<double>;
    o_t0        = new std::vector<double>;
    o_goodHitFlag = new std::vector<bool>;
    o_ped = new std::vector<double>;
    o_adccut = new std::vector<double>;
    o_choosePeakFlag = new std::vector<bool>;
    o_driftTAll = new std::vector<std::vector<double> >;

	std::vector<double> t_driftT;

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

        // get truth
        bool found = false;
        for (int i = 0; i<McTruth_tid->size(); i++){
            if ((*McTruth_tid)[i]==2){
                ox = (*McTruth_x)[i]*10;
                oy = (*McTruth_y)[i]*10;
                oz = (*McTruth_z)[i]*10;
                opx = (*McTruth_px)[i]*1000;
                opy = (*McTruth_py)[i]*1000;
                opz = (*McTruth_pz)[i]*1000;
                found = true;
                break;
            }
        }
        if (!found) continue;

        // get initial value
        double slx = m2y-m1y==0?0:(m2x-m1x)/(m2y-m1y);
        double slz = m2y-m1y==0?0:(m2z-m1z)/(m2y-m1y);
        double inx = slx*(sciYup-m1y)+m1x;
        double inz = slz*(sciYup-m1y)+m1z;

        // preset
        triggerNumber = iEntry;
        nHit = 0;
        o_driftD->clear();
        o_driftT->clear();
        o_layerID->clear();
        o_cellID->clear();
        o_nHitLayer->clear();
        o_nPeaks->clear();
        o_peakWidth->clear();
        o_adcPeak->clear();
        o_q->clear();
        o_t0->clear();
        o_goodHitFlag->clear();
        o_ped->clear();
        o_adccut->clear();
        o_choosePeakFlag->clear();
        o_driftTAll->clear();
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
            // ...
            //driftD+=random1.Gaus(0,0.2);

			// set output
			t_driftT.clear();
            t_driftT.push_back(driftT);
            o_driftT->push_back(driftT);
            o_driftD->push_back(driftD);
            o_layerID->push_back(lid);
            o_cellID->push_back(wid);
            o_nHitLayer->push_back(1);
            o_nPeaks->push_back(1);
            o_peakWidth->push_back(1);
            o_adcPeak->push_back(500);
            o_q->push_back(500);
            o_t0->push_back(0);
            o_goodHitFlag->push_back(1);
            o_ped->push_back(1);
            o_adccut->push_back(1);
            o_choosePeakFlag->push_back(1);
			o_driftTAll->push_back(t_driftT);
            nHitLayer[lid]++;
            nHit++;
        }
        for (int ihit = 0; ihit<o_nHitLayer->size(); ihit++){
            int lid = (*o_layerID)[ihit];
            (*o_nHitLayer)[ihit] = nHitLayer[lid];
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
