#include <stdio.h>  /* printf, getenv */
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "Log.hxx"
#include "InputOutputManager.hxx"
#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"
#include "XTManager.hxx"
#include "Track.hxx"

InputOutputManager* InputOutputManager::fInputOutputManager = NULL;

InputOutputManager::InputOutputManager():
    readRawFile(false),
    readPeakFile(false),
    readHitFile(false),
    readTrackFile(false),
    writeHitFile(false),
    writeTrackFile(false),
    writeAnaFile(false),
    hitFileIsMC(false),
    suffixHitFile(""),
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
    t0mc(0),
    nHitsG(0),
    nCandidatesFound(0),
    fOutputHitTree(0),
    fOutputHitFile(0),
    fOutputTrackTree(0),
    fOutputTrackFile(0),
    fInputRawChain(0),
    fInputPeakChain(0),
    fInputHitChain(0),
    fInputTrackChain(0)
{
}

InputOutputManager::~InputOutputManager(){
    if (LayerID) {delete LayerID;}
    if (CellID) {delete CellID;}
    if (TDCClock) {delete TDCClock;}
    if (DriftT) {delete DriftT;}
    if (DriftDmc) {delete DriftDmc;}
    if (DOCA) {delete DOCA;}
    if (Pedestal) {delete Pedestal;}
    if (ADCheight) {delete ADCheight;}
    if (ADCpeak) {delete ADCpeak;}
    if (ADCsumPacket) {delete ADCsumPacket;}
    if (ADCsumAll) {delete ADCsumAll;}
    if (PacketWidth) {delete PacketWidth;}
    if (nPeaksInChannel) {delete nPeaksInChannel;}
    if (iPeakInChannel) {delete iPeakInChannel;}
    if (nPeaksInPacket) {delete nPeaksInPacket;}
    if (iPeakInPacket) {delete iPeakInPacket;}
}

/// This is to set branches for input and output files
bool InputOutputManager::Initialize(bool withTrivialBranches){
    fCurrentEntry = 0;
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int runNo = RunInfoManager::Get().runNo;
    TString runName = RunInfoManager::Get().runName;
    int testLayer = RunInfoManager::Get().testLayer;

    if (!readHitFile){
        if (LayerID) {delete LayerID;} LayerID = new std::vector<int>;
        if (CellID) {delete CellID;} CellID = new std::vector<int>;
        if (TDCClock) {delete TDCClock;} TDCClock = new std::vector<int>;
        if (DriftT) {delete DriftT;} DriftT = new std::vector<double>;
        if (DriftDmc) {delete DriftDmc;} DriftDmc = new std::vector<double>; // For MC input
        if (DOCA) {delete DOCA;} DOCA = new std::vector<double>; // For MC input
        if (Pedestal) {delete Pedestal;} Pedestal = new std::vector<double>;
        if (ADCheight) {delete ADCheight;} ADCheight = new std::vector<int>;
        if (ADCpeak) {delete ADCpeak;} ADCpeak = new std::vector<int>;
        if (ADCsumPacket) {delete ADCsumPacket;} ADCsumPacket = new std::vector<double>;
        if (ADCsumAll) {delete ADCsumAll;} ADCsumAll = new std::vector<double>;
        if (PacketWidth) {delete PacketWidth;} PacketWidth = new std::vector<int>;
        if (nPeaksInChannel) {delete nPeaksInChannel;} nPeaksInChannel = new std::vector<int>;
        if (iPeakInChannel) {delete iPeakInChannel;} iPeakInChannel = new std::vector<int>;
        if (nPeaksInPacket) {delete nPeaksInPacket;} nPeaksInPacket = new std::vector<int>;
        if (iPeakInPacket) {delete iPeakInPacket;} iPeakInPacket = new std::vector<int>;
    }

    if (readRawFile){
        if (fInputRawChain) delete fInputRawChain;
        fInputRawChain = new TChain("tree","tree");
        fInputRawChain->Add(HOME+Form("/root/raw/run_%0.6d_built.root",runNo));
        fInputRawChain->SetBranchAddress("triggerNumber",&triggerNumber);
        fInputRawChain->SetBranchAddress("tdcNhit",tdcNhit);
        fInputRawChain->SetBranchAddress("clockNumberDriftTime",clockNumber);
        fInputRawChain->SetBranchAddress("driftTime",tdc);
        fInputRawChain->SetBranchAddress("adc",adc);
    }
    if (readHitFile){
        if (fInputHitChain) delete fInputHitChain;
        fInputHitChain = new TChain("t","t");
        fInputHitChain->Add(HOME+Form("/root/hits/h_%d%s.root",runNo,suffixHitFile.Data()));
        fInputHitChain->SetBranchAddress("triggerNumber",&triggerNumber);
        fInputHitChain->SetBranchAddress("nHits",&nHits);
        fInputHitChain->SetBranchAddress("driftT",&DriftT);
        if (hitFileIsMC){
            fInputHitChain->SetBranchAddress("DOCA",&DOCA);
            fInputHitChain->SetBranchAddress("driftD",&DriftDmc);
        }
        fInputHitChain->SetBranchAddress("layerID",&LayerID);
        fInputHitChain->SetBranchAddress("wireID",&CellID);
        //fInputHitChain->SetBranchAddress("type",&type); // 0 center, 1 left, 2 right, 3 guard, 4 dummy
        fInputHitChain->SetBranchAddress("ped",&Pedestal);
        fInputHitChain->SetBranchAddress("height",&ADCheight);
        fInputHitChain->SetBranchAddress("peak",&ADCpeak);
        //fInputHitChain->SetBranchAddress("rank",&rank);
        fInputHitChain->SetBranchAddress("sum",&ADCsumPacket);
        fInputHitChain->SetBranchAddress("aa",&ADCsumAll);
        fInputHitChain->SetBranchAddress("width",&PacketWidth);
        fInputHitChain->SetBranchAddress("np",&nPeaksInChannel);
        fInputHitChain->SetBranchAddress("ip",&iPeakInChannel);
        fInputHitChain->SetBranchAddress("mpn",&nPeaksInPacket);
        fInputHitChain->SetBranchAddress("mpi",&iPeakInPacket);
        fInputHitChain->SetBranchAddress("clk",&TDCClock);
        if (hitFileIsMC){
            fInputHitChain->SetBranchAddress("inxmc",&interceptXmc);
            fInputHitChain->SetBranchAddress("inzmc",&interceptZmc);
            fInputHitChain->SetBranchAddress("slxmc",&slopeXmc);
            fInputHitChain->SetBranchAddress("slzmc",&slopeZmc);
            fInputHitChain->SetBranchAddress("t0mc",&t0mc);
        }
    }
    if (readTrackFile){
        if (fInputTrackChain) delete fInputTrackChain; // FIXME: double delete?
        fInputTrackChain = new TChain("t","t");
        fInputTrackChain->Add(Form("%s/root/tracks/t_%d.%s.layer%d.root",HOME.Data(),runNo,runName.Data(),testLayer));
        // from h_XXX
        fInputTrackChain->SetBranchAddress("triggerNumber",&triggerNumber);
        // basic
        fInputTrackChain->SetBranchAddress("nHitsG",&nHitsG); // number of good hits in layers other than the test one: in t region and with good peak quality
        fInputTrackChain->SetBranchAddress("nFind",&nCandidatesFound);
        fInputTrackChain->SetBranchAddress("nPairs",nPairs);
        fInputTrackChain->SetBranchAddress("iSelection",iSelection);
        fInputTrackChain->SetBranchAddress("iCombination",iCombination);
        fInputTrackChain->SetBranchAddress("nHitsS",nHitsS); // number of hits selected from finding and fed to fitting
        for (int iLayer = 0; iLayer<NLAY; iLayer++){
            fInputTrackChain->SetBranchAddress(Form("hitIndexSelectedInLayer%d",iLayer),hitIndexSelected[iLayer]); // number of hits selected from finding and fed to fitting
        }
        fInputTrackChain->SetBranchAddress("t0Offset",t0Offset);
        fInputTrackChain->SetBranchAddress("interceptX",interceptX);
        fInputTrackChain->SetBranchAddress("interceptZ",interceptZ);
        fInputTrackChain->SetBranchAddress("slopeX",slopeX);
        fInputTrackChain->SetBranchAddress("slopeZ",slopeZ);
        fInputTrackChain->SetBranchAddress("chi2",chi2);
        fInputTrackChain->SetBranchAddress("chi2a",chi2a);
        fInputTrackChain->SetBranchAddress("chi2WithTestLayer",chi2WithTestLayer);
        fInputTrackChain->SetBranchAddress("pValue",pValue);
        if (withTrivialBranches){
            fInputTrackChain->SetBranchAddress("interceptXInput",interceptXInput);
            fInputTrackChain->SetBranchAddress("interceptZInput",interceptZInput);
            fInputTrackChain->SetBranchAddress("slopeXInput",slopeXInput);
            fInputTrackChain->SetBranchAddress("slopeZInput",slopeZInput);
            fInputTrackChain->SetBranchAddress("chi2XInput",chi2XInput);
            fInputTrackChain->SetBranchAddress("chi2ZInput",chi2ZInput);
            fInputTrackChain->SetBranchAddress("chi2Input",chi2Input);
            fInputTrackChain->SetBranchAddress("chi2aInput",chi2aInput);
            fInputTrackChain->SetBranchAddress("chi2WithTestLayerInput",chi2WithTestLayerInput);
            fInputTrackChain->SetBranchAddress("pValueInput",pValueInput);
        }
    }

    //===================Prepare output ROOT file============================
    if (writeHitFile){
        if (fOutputHitTree) delete fOutputHitTree;
        if (fOutputHitFile) fOutputHitFile->Close();
        fInputHitChain = new TChain("t","t");
        fOutputHitFile = new TFile(Form("%s/root/hits/h_%d%s.root",HOME.Data(),runNo,suffixHitFile.Data()),"RECREATE");
        fOutputHitTree = new TTree("t","t");
        fOutputHitTree->Branch("triggerNumber",&triggerNumber);
        fOutputHitTree->Branch("nHits",&nHits);
        fOutputHitTree->Branch("driftT",&DriftT);
        if (hitFileIsMC){
            fOutputHitTree->Branch("DOCA",&DOCA);
            fOutputHitTree->Branch("driftD",&DriftDmc);
        }
        fOutputHitTree->Branch("layerID",&LayerID);
        fOutputHitTree->Branch("wireID",&CellID);
        //fOutputHitTree->Branch("type",&type); // 0 center, 1 left, 2 right, 3 guard, 4 dummy
        fOutputHitTree->Branch("ped",&Pedestal);
        fOutputHitTree->Branch("height",&ADCheight);
        fOutputHitTree->Branch("peak",&ADCpeak);
        //fOutputHitTree->Branch("rank",&rank);
        fOutputHitTree->Branch("sum",&ADCsumPacket);
        fOutputHitTree->Branch("aa",&ADCsumAll);
        fOutputHitTree->Branch("width",&PacketWidth);
        fOutputHitTree->Branch("np",&nPeaksInChannel);
        fOutputHitTree->Branch("ip",&iPeakInChannel);
        fOutputHitTree->Branch("mpn",&nPeaksInPacket);
        fOutputHitTree->Branch("mpi",&iPeakInPacket);
        fOutputHitTree->Branch("clk",&TDCClock);
        if (hitFileIsMC){
            fOutputHitTree->Branch("inxmc",&interceptXmc);
            fOutputHitTree->Branch("inzmc",&interceptZmc);
            fOutputHitTree->Branch("slxmc",&slopeXmc);
            fOutputHitTree->Branch("slzmc",&slopeZmc);
            fOutputHitTree->Branch("t0mc",&t0mc);
        }
    }
    if (writeTrackFile){
        if (fOutputTrackTree) delete fOutputTrackTree; // FIXME: double delete?
        if (fOutputTrackFile) fOutputTrackFile->Close();
        fOutputTrackFile = new TFile(Form("%s/root/tracks/t_%d.%s.layer%d.root",HOME.Data(),runNo,runName.Data(),testLayer),"RECREATE");
        fOutputTrackTree = new TTree("t","t");
        // from h_XXX
        fOutputTrackTree->Branch("triggerNumber",&triggerNumber);
        // basic
        fOutputTrackTree->Branch("nHitsG",&nHitsG); // number of good hits in layers other than the test one: in t region and with good peak quality
        fOutputTrackTree->Branch("nFind",&nCandidatesFound);
        fOutputTrackTree->Branch("nPairs",nPairs,"nPairs[nFind]/I");
        fOutputTrackTree->Branch("iSelection",iSelection,"iSelection[nFind]/I");
        fOutputTrackTree->Branch("iCombination",iCombination,"iCombination[nFind]/I");
        fOutputTrackTree->Branch("nHitsS",nHitsS,"nHitsS[nFind]/I"); // number of hits selected from finding and fed to fitting
        for (int iLayer = 0; iLayer<NLAY; iLayer++){
            fOutputTrackTree->Branch(Form("hitIndexSelectedInLayer%d",iLayer),hitIndexSelected[iLayer],Form("hitIndexSelectedInLayer%d[nFind]/I",iLayer)); // number of hits selected from finding and fed to fitting
        }
        fOutputTrackTree->Branch("t0Offset",t0Offset,"t0Offset[nFind]/D");
        fOutputTrackTree->Branch("interceptX",interceptX,"interceptX[nFind]/D");
        fOutputTrackTree->Branch("interceptZ",interceptZ,"interceptZ[nFind]/D");
        fOutputTrackTree->Branch("slopeX",slopeX,"slopeX[nFind]/D");
        fOutputTrackTree->Branch("slopeZ",slopeZ,"slopeZ[nFind]/D");
        fOutputTrackTree->Branch("chi2",chi2,"chi2[nFind]/D");
        fOutputTrackTree->Branch("chi2a",chi2a,"chi2a[nFind]/D");
        fOutputTrackTree->Branch("chi2WithTestLayer",chi2WithTestLayer,"chi2WithTestLayer[nFind]/D");
        fOutputTrackTree->Branch("pValue",pValue,"pValue[nFind]/D");
        if (withTrivialBranches){
            fOutputTrackTree->Branch("interceptXInput",interceptXInput,"interceptXInput[nFind]/D");
            fOutputTrackTree->Branch("interceptZInput",interceptZInput,"interceptZInput[nFind]/D");
            fOutputTrackTree->Branch("slopeXInput",slopeXInput,"slopeXInput[nFind]/D");
            fOutputTrackTree->Branch("slopeZInput",slopeZInput,"slopeZInput[nFind]/D");
            fOutputTrackTree->Branch("chi2XInput",chi2XInput,"chi2XInput[nFind]/D");
            fOutputTrackTree->Branch("chi2ZInput",chi2ZInput,"chi2ZInput[nFind]/D");
            fOutputTrackTree->Branch("chi2Input",chi2Input,"chi2Input[nFind]/D");
            fOutputTrackTree->Branch("chi2aInput",chi2aInput,"chi2aInput[nFind]/D");
            fOutputTrackTree->Branch("chi2WithTestLayerInput",chi2WithTestLayerInput,"chi2WithTestLayerInput[nFind]/D");
            fOutputTrackTree->Branch("pValueInput",pValueInput,"pValueInput[nFind]/D");
        }
    }

    return true;
}

/// Will set initial values to input information so that if the GetEntry() function doesn't work, the user will get illegal values instead of misleading values from a previous entry
void InputOutputManager::Reset(){ // called at the beginning of every event
    // prepare
    nHitsG = 0;
    nCandidatesFound = 0;
    for (int iCand = 0; iCand<NCAND; iCand++){
        nPairs[iCand] = 0;
        iSelection[iCand] = -1;
        iCombination[iCand] = -1;
        nHitsS[iCand] = 0;
        for (int iLayer = 0; iLayer<NLAY; iLayer++){
            hitIndexSelected[iLayer][iCand] = -1;
        }
        t0Offset[iCand] = 0; // in case t0 is set free to adjustment
        interceptXInput[iCand] = 0;
        interceptZInput[iCand] = 0;
        slopeXInput[iCand] = 0;
        slopeZInput[iCand] = 0;
        chi2XInput[iCand] = -1;
        chi2ZInput[iCand] = -1;
        chi2Input[iCand] = -1;
        chi2aInput[iCand] = -1;
        chi2WithTestLayerInput[iCand] = -1;
        pValueInput[iCand] = -1;
        interceptX[iCand] = 0;
        interceptZ[iCand] = 0;
        slopeX[iCand] = 0;
        slopeZ[iCand] = 0;
        chi2[iCand] = -1;
        chi2a[iCand] = -1;
        chi2WithTestLayer[iCand] = -1;
        pValue[iCand] = -1;
        chi2mc[iCand] = -1;
        chi2amc[iCand] = -1;
        chi2WithTestLayermc[iCand] = -1;
        pValuemc[iCand] = -1;
    }
    if (!readHitFile&&writeHitFile){
        LayerID->clear();
        CellID->clear();
        TDCClock->clear();
        DriftT->clear();
        DriftDmc->clear();
        DOCA->clear();
        Pedestal->clear();
        ADCheight->clear();
        ADCpeak->clear();
        ADCsumPacket->clear();
        ADCsumAll->clear();
        PacketWidth->clear();
        nPeaksInChannel->clear();
        iPeakInChannel->clear();
        nPeaksInPacket->clear();
        iPeakInPacket->clear();
        nHits = 0;
    }
    interceptXmc = 0;
    interceptZmc = 0;
    slopeXmc = 0;
    slopeZmc = 0;
    t0mc = 0;
}

void InputOutputManager::Fill(){
    if (writeHitFile) fOutputHitTree->Fill();
    if (writeTrackFile) fOutputTrackTree->Fill();
}

void InputOutputManager::Write(){
    if (writeHitFile) fOutputHitTree->Write();
    if (writeTrackFile) fOutputTrackTree->Write();
}

void InputOutputManager::Close(){
    if (writeHitFile) fOutputHitFile->Close();
    if (writeTrackFile) fOutputTrackFile->Close();
}

void InputOutputManager::GetEntry(Long64_t iEntry){
    if (readRawFile&&fInputRawChain) fInputRawChain->GetEntry(iEntry); 
    if (readPeakFile&&fInputPeakChain) fInputPeakChain->GetEntry(iEntry); 
    if (readHitFile&&fInputHitChain) fInputHitChain->GetEntry(iEntry); 
    if (readTrackFile&&fInputTrackChain) fInputTrackChain->GetEntry(iEntry); 
    fCurrentEntry = iEntry;
}

bool InputOutputManager::IsRawFileReady(void){
    bool isReady = false;
    if (readRawFile&&fInputRawChain&&fInputRawChain->GetEntries()>0) isReady = true;
    return isReady;
}

bool InputOutputManager::IsPeakFileReady(void){
    bool isReady = false;
    if (readPeakFile&&fInputPeakChain&&fInputPeakChain->GetEntries()>0) isReady = true;
    return isReady;
}

bool InputOutputManager::IsHitFileReady(void){
    bool isReady = false;
    if (readHitFile&&fInputHitChain&&fInputHitChain->GetEntries()>0) isReady = true;
    return isReady;
}

bool InputOutputManager::IsTrackFileReady(void){
    bool isReady = false;
    if (readTrackFile&&fInputTrackChain&&fInputTrackChain->GetEntries()>0) isReady = true;
    return isReady;
}

Long64_t InputOutputManager::GetEntries(){
    Long64_t N = 0;
    TChain * chain = getChain();
    if (chain){
        N = chain->GetEntries();
    }
    return N;
}

Long64_t InputOutputManager::GetTriggerNumberMax(){
    Long64_t TriMax = 1;
    TChain * chain = getChain();
    if (chain&&chain->GetEntries()>0){
        chain->GetEntry(chain->GetEntries()-1);
        TriMax = triggerNumber;
    }
    return TriMax;
}

void InputOutputManager::Print(TString opt){
    printf("Entry %d, triggerNumber %d:\n",fCurrentEntry,triggerNumber);
    if (readHitFile||writeHitFile){
        printf("  Total hits: %d, good hits after cuts: %d\n",nHits,nHitsG);
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
                // TODO: add MC hit info
                if (iPeakA==0){
                    printf("[%2d,%2d] :  %3d TDCs,  total ADC sum %.1f, pedestal %.1f\n",lid,cid,nPeakA,sumA,ped);
                    iPacket = 0;
                }
                if (iPeakL==0){
                    printf("             packet with %2d peaks : ADC sum %.1f, peak %d, %d ticks wide.\n",nPeakL,sumL,peak,width);
                    iPacket++;
                }
                printf("%4d        %.1f ns, ADC %d @ clk %d\n",iHit,driftT,height,clock);
            }
        }
    }
}

bool InputOutputManager::SetTrack(int iFound, const Track3D * track3D, const Track2D * track2D){
    if (iFound>=NCAND){
        MyError("The fitting result index "<<iFound<<" exceeded the maximum capacity "<<NCAND);
        return false;
    }
    iSelection[iFound] = track3D->iSelection;
    iCombination[iFound] = track3D->iCombination;
    nHitsS[iFound] = track3D->hitIndexSelected.size();
    t0Offset[iFound] = track3D->t0Offset;
    interceptX[iFound] = track3D->interceptX;
    interceptZ[iFound] = track3D->interceptZ;
    slopeX[iFound] = track3D->slopeX;
    slopeZ[iFound] = track3D->slopeZ;
    chi2[iFound] = track3D->chi2;
    pValue[iFound] = track3D->pValue;
    chi2a[iFound] = track3D->chi2a;
    chi2WithTestLayer[iFound] = track3D->chi2WithTestLayer;
    for (int lid = 0; lid<NLAY; lid++){
        hitIndexSelected[lid][iFound] = -1;
    }
    for (size_t iHit = 0; iHit<track3D->hitIndexSelected.size(); iHit++){
        int lid = LayerID->at(track3D->hitIndexSelected[iHit]);
        if (lid>=0&&lid<NLAY){
            hitIndexSelected[lid][iFound] = track3D->hitIndexSelected[iHit];
        }
    }
    if (track2D){
        nPairs[iFound] = track2D->nPairs;
        interceptXInput[iFound] = track2D->interceptX;
        interceptZInput[iFound] = track2D->interceptZ;
        slopeXInput[iFound] = track2D->slopeX;
        slopeZInput[iFound] = track2D->slopeZ;
        chi2XInput[iFound] = track2D->chi2X;
        chi2ZInput[iFound] = track2D->chi2Z;
        chi2Input[iFound] = track2D->chi2;
        chi2aInput[iFound] = track2D->chi2a;
        chi2WithTestLayerInput[iFound] = track2D->chi2WithTestLayer;
        pValueInput[iFound] = track2D->pValue;
    }
    return true;
}

void InputOutputManager::PushHitMC(int lid, int wid, double driftT, double driftD, double doca){
    // TODO: currently this is for kMCDriftT; Add support for kMCDriftD
    int status;
    LayerID->push_back(lid);
    CellID->push_back(wid);
    TDCClock->push_back(0);
    DriftT->push_back(driftT);
    DriftDmc->push_back(driftD);
    DOCA->push_back(doca); // For MC input
    Pedestal->push_back(0);
    ADCheight->push_back(0);
    ADCpeak->push_back(0);
    ADCsumPacket->push_back(0);
    ADCsumAll->push_back(0);
    PacketWidth->push_back(0);
    nPeaksInChannel->push_back(0);
    iPeakInChannel->push_back(0);
    nPeaksInPacket->push_back(0);
    iPeakInPacket->push_back(0);
    nHits++;
}

TChain * InputOutputManager::getChain(){
    if (readTrackFile&&fInputTrackChain)
        return fInputTrackChain;
    else if (readHitFile&&fInputHitChain)
        return fInputHitChain;
    else if (readPeakFile&&fInputPeakChain)
        return fInputPeakChain;
    else if (readRawFile&&fInputRawChain)
        return fInputRawChain;
    else
        return 0;
}

