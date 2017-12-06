#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"

#define NLAY 9

void printUsage(char* name);
int getwid(int layerID, int cellID);

int main(int argc, char ** argv){

    if (argc<2){
        printUsage(argv[0]);
        return -1;
    }
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
                m1y = (*M_y)[i]*10;
                m1z = (*M_z)[i]*10;
                m1x = (*M_x)[i]*10;
                checkup = true;
            }
            else if ((*M_vid)[i]==1){
                m2y = (*M_y)[i]*10;
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
            if ((*CdcCell_tid)[ihit]!=2) continue;
            double driftD = (*CdcCell_driftDtrue)[ihit]*10;
            driftD+=random1.Gaus(0,0.2);
            double driftT=driftD/0.023;
            int layerID = (*CdcCell_layerID)[ihit];
            int cellID = (*CdcCell_cellID)[ihit];
            int lid = layerID;
            int wid = getwid(layerID,cellID);
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
