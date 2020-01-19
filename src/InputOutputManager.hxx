#ifndef InputOutputManager_hxx_seen
#define InputOutputManager_hxx_seen

#define NCAND  10000
#define NCHMAX 5000
#define NSAM   32

#include <TString.h>

#include "GeometryManager.hxx"

class TFile;
class TTree;
class TChain;
class Track2D;
class Track3D;

/// This class is to manage input and output files
///
/// This is a singleton and you can use InputOutputManager::Get() to get its reference \n
/// Initialize() function will set up input and output files according to types given by ParameterManager \n
/// This manager is supposed to manage all kinds of input and output files including \n
/// -# raw data file as input
/// -# peaks file as input
/// -# hits file as input, including Data, MCDriftT, MCDriftD, defined as InputHitType
/// -# track file as input
/// -# track file as output
/// -# analysis file as output
///
/// Set flags directly to enable/disable inputs/outputs. By default none of them are enabled \n
/// In each event the user should follow the example bellow: \n
///\code{.cpp}
///    InputOutputManager::Get().Reset();
///    InputOutputManager::Get().GetEntry(iEntry);
///    // User code to get values and set to InputOutputManager...
///    InputOutputManager::Get().Fill();
///\endcode
///
/// Don't forget to call Write() and Close() at the end of your user code. \n
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

    enum InputHitType{
        kData,
        kMCDriftT,
        kMCDriftD
    };

    bool Initialize(bool withTrivialBranches = false); ///< Should be called at the beginning of the run
    void Reset(); ///< Should be called at the beginning of every event
    void Fill(); ///< Fill the branches in the output file
    void Write(); ///< Write the tree to the output file
    void Close(); ///< Close the output file
    void GetEntry(Long64_t iEntry); ///< Should be called at the beginning of every event AFTER Reset()
    Long64_t GetEntries();
    void Print(TString opt = "");
    bool SetTrack(int iFound, const Track3D* track3D, const Track2D* track2D = 0);

    /// flags about which to read and which to write
    bool                  readRawFile;
    bool                  readPeakFile;
    bool                  readHitFile;
    bool                  readTrackFile;
    bool                  writeTrackFile;
    bool                  writeAnaFile;

    /// the current entry; increment once GetEntry(iEntry) is called
    int                   fCurrentEntry;

    int                   triggerNumber;
    /// raw
    int                   tdcNhit[NCHMAX];
    int                   clockNumber[NCHMAX][NSAM];
    int                   adc[NCHMAX][NSAM];
    int                   tdc[NCHMAX][NSAM];

    /// hits
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

    /// pre-fitting information
    int                   nHitsG;
    int                   nCandidatesFound;
    // fitting results
    int                   iSelection[NCAND];
    int                   iCombination[NCAND];
    int                   nHitsS[NCAND];
    int                   hitIndexSelected[NLAY][NCAND];
    double                t0Offset[NCAND]; // in case t0 is set free to adjustment
    double                interceptX[NCAND];
    double                interceptZ[NCAND];
    double                slopeX[NCAND];
    double                slopeZ[NCAND];
    double                chi2[NCAND];
    double                chi2a[NCAND];
    double                chi2WithTestLayer[NCAND];
    double                pValue[NCAND];
    double                chi2mc[NCAND];
    double                chi2amc[NCAND];
    double                chi2WithTestLayermc[NCAND];
    double                pValuemc[NCAND];
    // trivial fitting informaiton
    int                   nPairs[NCAND];
    double                interceptXInput[NCAND];
    double                interceptZInput[NCAND];
    double                slopeXInput[NCAND];
    double                slopeZInput[NCAND];
    double                chi2XInput[NCAND];
    double                chi2ZInput[NCAND];
    double                chi2Input[NCAND];
    double                pValueInput[NCAND];
    double                chi2aInput[NCAND];
    double                chi2WithTestLayerInput[NCAND];

private:
    /// The static pointer to the singleton instance.
    static InputOutputManager* fInputOutputManager;

    TTree  * fOutputTrackTree;
    TFile  * fOutputTrackFile;
    TChain * fInputRawChain;
    TChain * fInputHitChain;
    TChain * fInputTrackChain;
};

#endif
