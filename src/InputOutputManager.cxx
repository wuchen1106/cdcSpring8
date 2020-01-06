#include <stdio.h>  /* printf, getenv */
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "Log.hxx"
#include "InputOutputManager.hxx"
#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"
#include "Track.hxx"

InputOutputManager* InputOutputManager::fInputOutputManager = NULL;

InputOutputManager::InputOutputManager():
    readRawFile(false),
    readPeakFile(false),
    readHitFile(false),
    readTrackFile(false),
    writeTrackFile(false),
    writeAnaFile(false),
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
    fOutputTrackTree(0),
    fOutputTrackFile(0),
    fInputHitChain(0),
    fInputTrackChain(0)
{
}

InputOutputManager::~InputOutputManager(){
}

/// This is to set branches for input and output files
bool InputOutputManager::Initialize(bool withTrivialBranches){
    fCurrentEntry = 0;
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int runNo = RunInfoManager::Get().runNo;
    TString runName = RunInfoManager::Get().runName;
    int testLayer = RunInfoManager::Get().testLayer;
    InputHitType inputHitType = ParameterManager::Get().inputHitType;

    // Need to load hit file?
    if (readHitFile){
        if (fInputHitChain) delete fInputHitChain;
        fInputHitChain = new TChain("t","t");
        if (inputHitType==kData)
            fInputHitChain->Add(HOME+Form("/root/hits/h_%d.root",runNo));
        else
            fInputHitChain->Add(HOME+Form("/root/hits/h_%d.MC.root",runNo));
        fInputHitChain->SetBranchAddress("triggerNumber",&triggerNumber);
        fInputHitChain->SetBranchAddress("nHits",&nHits);
        fInputHitChain->SetBranchAddress("driftT",&DriftT);
        if (inputHitType==kMCDriftD||inputHitType==kMCDriftT) fInputHitChain->SetBranchAddress("DOCA",&DOCA);
        if (inputHitType==kMCDriftD) fInputHitChain->SetBranchAddress("driftD",&DriftDmc);
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
        if (inputHitType==kMCDriftD||inputHitType==kMCDriftT){
            fInputHitChain->SetBranchAddress("inxmc",&interceptXmc);
            fInputHitChain->SetBranchAddress("inzmc",&interceptZmc);
            fInputHitChain->SetBranchAddress("slxmc",&slopeXmc);
            fInputHitChain->SetBranchAddress("slzmc",&slopeZmc);
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
        fInputTrackChain->SetBranchAddress("nPairsG",nGoodPairs);
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
        fInputTrackChain->SetBranchAddress("chi2WithTestLayer",chi2WithTestLayer);
        fInputTrackChain->SetBranchAddress("pValue",pValue);
        if (inputHitType!=kData){
            fInputTrackChain->SetBranchAddress("chi2mc",chi2mc);
            fInputTrackChain->SetBranchAddress("chi2WithTestLayermc",chi2WithTestLayermc);
            fInputTrackChain->SetBranchAddress("pValuemc",pValuemc);
            fInputTrackChain->SetBranchAddress("inxmc",&interceptXmc);
            fInputTrackChain->SetBranchAddress("inzmc",&interceptZmc);
            fInputTrackChain->SetBranchAddress("slxmc",&slopeXmc);
            fInputTrackChain->SetBranchAddress("slzmc",&slopeZmc);
        }
        if (withTrivialBranches){
            fInputTrackChain->SetBranchAddress("interceptXInput",interceptXInput);
            fInputTrackChain->SetBranchAddress("interceptZInput",interceptZInput);
            fInputTrackChain->SetBranchAddress("slopeXInput",slopeXInput);
            fInputTrackChain->SetBranchAddress("slopeZInput",slopeZInput);
            fInputTrackChain->SetBranchAddress("chi2XInput",chi2XInput);
            fInputTrackChain->SetBranchAddress("chi2ZInput",chi2ZInput);
            fInputTrackChain->SetBranchAddress("chi2Input",chi2Input);
            fInputTrackChain->SetBranchAddress("chi2WithTestLayerInput",chi2WithTestLayerInput);
            fInputTrackChain->SetBranchAddress("pValueInput",pValueInput);
        }
    }

    //===================Prepare output ROOT file============================
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
        fOutputTrackTree->Branch("nPairsG",nGoodPairs,"nPairsG[nFind]/I");
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
        fOutputTrackTree->Branch("chi2WithTestLayer",chi2WithTestLayer,"chi2WithTestLayer[nFind]/D");
        fOutputTrackTree->Branch("pValue",pValue,"pValue[nFind]/D");
        if (inputHitType!=kData){
            fOutputTrackTree->Branch("chi2mc",chi2mc,"chi2mc[nFind]/D");
            fOutputTrackTree->Branch("chi2WithTestLayermc",chi2WithTestLayermc,"chi2WithTestLayermc[nFind]/D");
            fOutputTrackTree->Branch("pValuemc",pValuemc,"pValuemc[nFind]/D");
            fOutputTrackTree->Branch("inxmc",&interceptXmc);
            fOutputTrackTree->Branch("inzmc",&interceptZmc);
            fOutputTrackTree->Branch("slxmc",&slopeXmc);
            fOutputTrackTree->Branch("slzmc",&slopeZmc);
        }
        if (withTrivialBranches){
            fOutputTrackTree->Branch("interceptXInput",interceptXInput,"interceptXInput[nFind]/D");
            fOutputTrackTree->Branch("interceptZInput",interceptZInput,"interceptZInput[nFind]/D");
            fOutputTrackTree->Branch("slopeXInput",slopeXInput,"slopeXInput[nFind]/D");
            fOutputTrackTree->Branch("slopeZInput",slopeZInput,"slopeZInput[nFind]/D");
            fOutputTrackTree->Branch("chi2XInput",chi2XInput,"chi2XInput[nFind]/D");
            fOutputTrackTree->Branch("chi2ZInput",chi2ZInput,"chi2ZInput[nFind]/D");
            fOutputTrackTree->Branch("chi2Input",chi2Input,"chi2Input[nFind]/D");
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
        nGoodPairs[iCand] = 0;
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
    if (writeTrackFile) fOutputTrackTree->Fill();
}

void InputOutputManager::Write(){
    if (writeTrackFile) fOutputTrackTree->Write();
}

void InputOutputManager::Close(){
    if (writeTrackFile) fOutputTrackFile->Close();
}

void InputOutputManager::GetEntry(Long64_t iEntry){
    if (readHitFile&&fInputHitChain) fInputHitChain->GetEntry(iEntry); 
    if (readTrackFile&&fInputTrackChain) fInputTrackChain->GetEntry(iEntry); 
    fCurrentEntry = iEntry;
}

Long64_t InputOutputManager::GetEntries(){
    if (readTrackFile&&fInputTrackChain)
        return fInputTrackChain->GetEntries();
    else if (readHitFile&&fInputHitChain)
        return fInputHitChain->GetEntries();
    else
        return 0;
}

void InputOutputManager::Print(TString opt){
    printf("Entry %d, triggerNumber %d:\n",fCurrentEntry,triggerNumber);
    if (readHitFile){
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

bool InputOutputManager::SetTrack(int iFound, const TrackResult * trackResult){
    if (iFound>=NCAND){
        MyError("The fitting result index "<<iFound<<" exceeded the maximum capacity "<<NCAND);
        return false;
    }
    nHitsS[iFound] = trackResult->NDF;
    t0Offset[iFound] = trackResult->t0Offset;
    iSelection[iFound] = trackResult->initialTrackCandidate.iSelection;
    iCombination[iFound] = trackResult->initialTrackCandidate.iCombination;
    nPairs[iFound] = trackResult->initialTrackCandidate.nPairs;
    nGoodPairs[iFound] = trackResult->initialTrackCandidate.nGoodPairs;
    interceptXInput[iFound] = trackResult->initialTrackCandidate.interceptX;
    interceptZInput[iFound] = trackResult->initialTrackCandidate.interceptZ;
    slopeXInput[iFound] = trackResult->initialTrackCandidate.slopeX;
    slopeZInput[iFound] = trackResult->initialTrackCandidate.slopeZ;
    chi2XInput[iFound] = trackResult->initialTrackCandidate.chi2X;
    chi2ZInput[iFound] = trackResult->initialTrackCandidate.chi2Z;
    chi2Input[iFound] = trackResult->initialTrackCandidate.chi2;
    pValueInput[iFound] = trackResult->initialTrackCandidate.pValue;
    interceptX[iFound] = trackResult->interceptX;
    interceptZ[iFound] = trackResult->interceptZ;
    slopeX[iFound] = trackResult->slopeX;
    slopeZ[iFound] = trackResult->slopeZ;
    chi2[iFound] = trackResult->chi2;
    pValue[iFound] = trackResult->pValue;
    chi2WithTestLayer[iFound] = trackResult->chi2WithTestLayer;
    for (int lid = 0; lid<NLAY; lid++){
        hitIndexSelected[lid][iFound] = -1;
    }
    for (size_t iHit = 0; iHit<trackResult->hitIndexSelected.size(); iHit++){
        int lid = LayerID->at(trackResult->hitIndexSelected[iHit]);
        if (lid>=0&&lid<NLAY){
            hitIndexSelected[lid][iFound] = trackResult->hitIndexSelected[iHit];
        }
    }
    return true;
}
