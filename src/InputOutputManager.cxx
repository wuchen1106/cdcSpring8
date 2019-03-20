#include <stdio.h>  /* printf, getenv */
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "Log.hxx"
#include "InputOutputManager.hxx"
#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"

InputOutputManager* InputOutputManager::fInputOutputManager = NULL;

InputOutputManager::InputOutputManager():
    fCurrentEntry(0),
    triggerNumber(0),
    nHits(0),
    LayerID(0),
    CellID(0),
    TDCClock(0),
    DriftT(0),
    DriftDmc(0),
    DOCA(0),
    Pedestal(0),
    ADCheight(0),
    ADCpeak(0),
    //rank(0),
    ADCsumPacket(0),
    ADCsumAll(0),
    PacketWidth(0),
    nPeaksInChannel(0),
    iPeakInChannel(0),
    nPeaksInPacket(0),
    iPeakInPacket(0),
    interceptXmc(0),
    interceptZmc(0),
    slopeXmc(0),
    slopeZmc(0),
    nHitsG(0),
    nCandidatesFound(0),
    fOutputTree(0),
    fOutputFile(0),
    fInputTChain(0)
{
}

InputOutputManager::~InputOutputManager(){
}

bool InputOutputManager::Initialize(){
    fCurrentEntry = 0;
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int runNo = RunInfoManager::Get().runNo;
    TString runName = RunInfoManager::Get().runName;
    int testLayer = RunInfoManager::Get().testLayer;
    InputType inputType = ParameterManager::Get().inputType;
    if (fOutputTree) delete fOutputTree; // FIXME: double delete?
    if (fOutputFile) fOutputFile->Close();
    if (fInputTChain) delete fInputTChain;
    fInputTChain = new TChain("t","t");
    fOutputFile = new TFile(Form("%s/root/tracks/t_%d.%s.layer%d.root",HOME.Data(),runNo,runName.Data(),testLayer),"RECREATE");
    fOutputTree = new TTree("t","t");

    if (inputType==kData)
        fInputTChain->Add(HOME+Form("/root/hits/h_%d.root",runNo));
    else
        fInputTChain->Add(HOME+Form("/root/hits/h_%d.MC.root",runNo));
    fInputTChain->SetBranchAddress("triggerNumber",&triggerNumber);
    fInputTChain->SetBranchAddress("nHits",&nHits);
    fInputTChain->SetBranchAddress("driftT",&DriftT);
    if (inputType==kMCDriftD||inputType==kMCDriftT) fInputTChain->SetBranchAddress("DOCA",&DOCA);
    if (inputType==kMCDriftD) fInputTChain->SetBranchAddress("driftD",&DriftDmc);
    fInputTChain->SetBranchAddress("layerID",&LayerID);
    fInputTChain->SetBranchAddress("wireID",&CellID);
    //fInputTChain->SetBranchAddress("type",&type); // 0 center, 1 left, 2 right, 3 guard, 4 dummy
    fInputTChain->SetBranchAddress("ped",&Pedestal);
    fInputTChain->SetBranchAddress("height",&ADCheight);
    fInputTChain->SetBranchAddress("peak",&ADCpeak);
    //fInputTChain->SetBranchAddress("rank",&rank);
    fInputTChain->SetBranchAddress("sum",&ADCsumPacket);
    fInputTChain->SetBranchAddress("aa",&ADCsumAll);
    fInputTChain->SetBranchAddress("width",&PacketWidth);
    fInputTChain->SetBranchAddress("np",&nPeaksInChannel);
    fInputTChain->SetBranchAddress("ip",&iPeakInChannel);
    fInputTChain->SetBranchAddress("mpn",&nPeaksInPacket);
    fInputTChain->SetBranchAddress("mpi",&iPeakInPacket);
    fInputTChain->SetBranchAddress("clk",&TDCClock);
    if (inputType==kMCDriftD||inputType==kMCDriftT){
		fInputTChain->SetBranchAddress("inxmc",&interceptXmc);
		fInputTChain->SetBranchAddress("inzmc",&interceptZmc);
		fInputTChain->SetBranchAddress("slxmc",&slopeXmc);
		fInputTChain->SetBranchAddress("slzmc",&slopeZmc);
    }

    //===================Prepare output ROOT file============================
    // from h_XXX
    fOutputTree->Branch("triggerNumber",&triggerNumber);
    // basic
    int nHitsG;
    fOutputTree->Branch("nHitsG",&nHitsG); // number of good hits in layers other than the test one: in t region and with good peak quality
    fOutputTree->Branch("nFind",&nCandidatesFound);
    fOutputTree->Branch("nPairs",nPairs,"nPairs[nFind]/I");
    fOutputTree->Branch("nHitsS",nHitsS,"nHitsS[nFind]/I"); // number of hits selected from finding and fed to fitting
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        fOutputTree->Branch(Form("hitIndexSelectedInLayer%d",iLayer),hitIndexSelected[iLayer],Form("hitIndexSelectedInLayer%d[nFind]/I",iLayer)); // number of hits selected from finding and fed to fitting
    }
    fOutputTree->Branch("t0offset",t0offset,"t0offset[nFind]/D");
    fOutputTree->Branch("interceptXInput",interceptXInput,"interceptXInput[nFind]/D");
    fOutputTree->Branch("interceptZInput",interceptZInput,"interceptZInput[nFind]/D");
    fOutputTree->Branch("slopeXInput",slopeXInput,"slopeXInput[nFind]/D");
    fOutputTree->Branch("slopeZInput",slopeZInput,"slopeZInput[nFind]/D");
    fOutputTree->Branch("chi2XInput",chi2XInput,"chi2XInput[nFind]/D");
    fOutputTree->Branch("chi2ZInput",chi2ZInput,"chi2ZInput[nFind]/D");
    fOutputTree->Branch("chi2Input",chi2Input,"chi2Input[nFind]/D");
    fOutputTree->Branch("chi2WithTestLayerInput",chi2WithTestLayerInput,"chi2WithTestLayerInput[nFind]/D");
    fOutputTree->Branch("pValueInput",pValueInput,"pValueInput[nFind]/D");
    fOutputTree->Branch("interceptX",interceptX,"interceptX[nFind]/D");
    fOutputTree->Branch("interceptZ",interceptZ,"interceptZ[nFind]/D");
    fOutputTree->Branch("slopeX",slopeX,"slopeX[nFind]/D");
    fOutputTree->Branch("slopeZ",slopeZ,"slopeZ[nFind]/D");
    fOutputTree->Branch("chi2",chi2,"chi2[nFind]/D");
    fOutputTree->Branch("chi2WithTestLayer",chi2WithTestLayer,"chi2WithTestLayer[nFind]/D");
    fOutputTree->Branch("pValue",pValue,"pValue[nFind]/D");
    if (inputType){
        fOutputTree->Branch("chi2mc",chi2mc,"chi2mc[nFind]/D");
        fOutputTree->Branch("chi2WithTestLayermc",chi2WithTestLayermc,"chi2WithTestLayermc[nFind]/D");
        fOutputTree->Branch("pValuemc",pValuemc,"pValuemc[nFind]/D");
		fOutputTree->Branch("inxmc",&interceptXmc);
		fOutputTree->Branch("inzmc",&interceptZmc);
		fOutputTree->Branch("slxmc",&slopeXmc);
		fOutputTree->Branch("slzmc",&slopeZmc);
	}

    return true;
}

void InputOutputManager::Reset(){ // called at the beginning of every event
    // prepare
    nHitsG = 0;
    nCandidatesFound = 0;
    for (int iCand = 0; iCand<NCAND; iCand++){
        nPairs[iCand] = 0;
        nHitsS[iCand] = 0;
        for (int iLayer = 0; iLayer<NLAY; iLayer++){
            hitIndexSelected[iLayer][iCand] = -1;
        }
        t0offset[iCand] = 0; // in case t0 is set free to adjustment
        interceptXInput[iCand] = 0;
        interceptZInput[iCand] = 0;
        slopeXInput[iCand] = 0;
        slopeZInput[iCand] = 0;
        chi2XInput[iCand] = -1;
        chi2ZInput[iCand] = -1;
        chi2Input[iCand] = -1;
        chi2WithTestLayerInput[iCand] = -1;
        pValueInput[iCand] = -1;
        interceptX[iCand] = 0;
        interceptZ[iCand] = 0;
        slopeX[iCand] = 0;
        slopeZ[iCand] = 0;
        chi2[iCand] = -1;
        chi2WithTestLayer[iCand] = -1;
        pValue[iCand] = -1;
        chi2mc[iCand] = -1;
        chi2WithTestLayermc[iCand] = -1;
        pValuemc[iCand] = -1;
    }
}

void InputOutputManager::Fill(){
    fOutputTree->Fill();
}

void InputOutputManager::Write(){
    fOutputTree->Write();
}

void InputOutputManager::Close(){
    fOutputFile->Close();
}

void InputOutputManager::GetEntry(Long64_t iEntry){
    if (fInputTChain) fInputTChain->GetEntry(iEntry); 
    fCurrentEntry = iEntry;
}

Long64_t InputOutputManager::GetEntries(){
    if (fInputTChain)
        return fInputTChain->GetEntries();
    else
        return 0;
}

void InputOutputManager::Print(TString opt){
    printf("Entry %d, triggerNumber %d:\n",fCurrentEntry,triggerNumber);
    printf("  Total hits: %d, good hits after cuts: %d\n",nHits,nHitsG);
    printf("  Total hits: %d, good hits after cuts: %d\n",LayerID->size(),nHitsG);
    if (opt.Contains("h")){
        int iPacket = 0;
        for (int iHit = 0; iHit<nHits; iHit++){
            int    lid    = LayerID->at(iHit);
            int    cid    = CellID->at(iHit);
            int    clock  = TDCClock->at(iHit);
            int    nPeakA = nPeaksInChannel->at(iHit);
            int    iPeakA = iPeakInChannel->at(iHit);
            int    nPeakL = nPeaksInPacket->at(iHit);
            int    iPeakL = iPeakInPacket->at(iHit);
            int    width  = PacketWidth->at(iHit);
            int    height = ADCheight->at(iHit);
            int    peak   = ADCpeak->at(iHit);
            double ped    = Pedestal->at(iHit);
            double sumA   = ADCsumAll->at(iHit);
            double sumL   = ADCsumPacket->at(iHit);
            double driftT = DriftT->at(iHit);
            if (iPeakA==0){
                printf("    layer %3d cell %3d:  %3d TDCs,  total ADC sum %.1f, pedestal %.1f\n",lid,cid,nPeakA,sumA,ped);
                iPacket = 0;
            }
            if (iPeakL==0){
                printf("        wave packet %3d: %3d peaks, local ADC sum %.1f, ADC peak %d, width %d ticks.\n",iPacket,nPeakL,sumL,peak,width);
                iPacket++;
            }
            printf("  %4d      peak %d @ clock %d: driftT = %.1f ns, ADC height at trigger point %d\n",iHit,iPeakL,clock,driftT,height);
        }
    }
}
