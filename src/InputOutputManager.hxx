#ifndef InputOutputManager_hxx_seen
#define InputOutputManager_hxx_seen

#define NCAND 10

#include <TString.h>

#include "GeometryManager.hxx"

class TFile;
class TTree;
class TChain;

/// This class is to manage input and output files
/// This is a singleton and you can use InputOutputManager::Get() to get its reference
/// enum type InputType and OutputType are defined here to indicate the format of input and output files
/// Initialize() function will set up input and output files according to types given by ParameterManager
/// In each event the user should follow the example bellow:
///\code{.cpp}
///InputOutputManager::Get().Reset();
///InputOutputManager::Get().GetEntry(iEntry);
///// User code to get values and set to InputOutputManager...
///InputOutputManager::Get().Fill();
///\endcode
/// Don't forget to call Write() and Close() at the end of your user code.
class InputOutputManager{
public:
    InputOutputManager();
    virtual ~InputOutputManager();
    
    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static InputOutputManager& Get(void) {
        if (!fInputOutputManager)
            fInputOutputManager = new InputOutputManager();
        return *fInputOutputManager;
    }

    enum InputType{
        kData,
        kMCDriftT,
        kMCDriftD
    };

    bool Initialize(); ///< Should be called at the beginning of the run
    void Reset(); ///< Should be called at the beginning of every event
    void Fill(); ///< Fill the branches in the output file
    void Write(); ///< Write the tree to the output file
    void Close(); ///< Close the output file
    void GetEntry(Long64_t iEntry); ///< Should be called at the beginning of every event AFTER Reset()
    Long64_t GetEntries();
    void Print(TString opt = "");

    int                   fCurrentEntry;

    // raw hits from input
    int                   triggerNumber;
    int                   nHits;
    std::vector<int> *    LayerID;
    std::vector<int> *    CellID;
    std::vector<int> *    TDCClock;
    std::vector<double> * DriftT;
    std::vector<double> * DriftDmc; // For MC input
    std::vector<double> * DOCA; // For MC input
    std::vector<double> * Pedestal;
    std::vector<int> *    ADCheight;
    std::vector<int> *    ADCpeak;
    //std::vector<int> *    rank;
    std::vector<double> * ADCsumPacket;
    std::vector<double> * ADCsumAll;
    std::vector<int> *    PacketWidth;
    std::vector<int> *    nPeaksInChannel;
    std::vector<int> *    iPeakInChannel;
    std::vector<int> *    nPeaksInPacket;
    std::vector<int> *    iPeakInPacket;
    double                interceptXmc; // For MC input
    double                interceptZmc; // For MC input
    double                slopeXmc; // For MC input
    double                slopeZmc; // For MC input

    // recon hits/tracks for output
    int                   nHitsG;
    int                   nCandidatesFound;
    int                   nPairs[NCAND];
    int                   nGoodPairs[NCAND];
    int                   iSelection[NCAND];
    int                   iCombination[NCAND];
    int                   nHitsS[NCAND];
    int                   hitIndexSelected[NLAY][NCAND];
    double                t0Offset[NCAND]; // in case t0 is set free to adjustment
    double                interceptXInput[NCAND];
    double                interceptZInput[NCAND];
    double                slopeXInput[NCAND];
    double                slopeZInput[NCAND];
    double                chi2XInput[NCAND];
    double                chi2ZInput[NCAND];
    double                chi2Input[NCAND];
    double                chi2WithTestLayerInput[NCAND];
    double                pValueInput[NCAND];
    double                interceptX[NCAND];
    double                interceptZ[NCAND];
    double                slopeX[NCAND];
    double                slopeZ[NCAND];
    double                chi2[NCAND];
    double                chi2WithTestLayer[NCAND];
    double                pValue[NCAND];
    double                chi2mc[NCAND];
    double                chi2WithTestLayermc[NCAND];
    double                pValuemc[NCAND];

private:
    /// The static pointer to the singleton instance.
    static InputOutputManager* fInputOutputManager;

    TFile  * fOutputFile;
    TTree  * fOutputTree;
    TChain * fInputTChain;
};

#endif
