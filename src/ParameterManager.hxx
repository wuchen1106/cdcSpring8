#ifndef ParameterManager_hxx_seen
#define ParameterManager_hxx_seen

#include <TString.h>
#include "GeometryManager.hxx"
#include "BeamManager.hxx"
#include "InputOutputManager.hxx"

// parameter block for tracking
class TrackingPara{
public:
    enum PeakType{
        kFirstPeak,
        kHighestPeak,
        kAllPeaks
    };

    TrackingPara();
    virtual ~TrackingPara(){};
    void Print();
    int    nHitsMax;
    int    nHitsSMin;
    int    nPairsMin;
    int    t0shift0;
    int    t0shift1;
    int    tmin;
    int    tmax;
    double sumCut;
    double aaCut;
    PeakType peakType;
    int    workType;
    int    BlindLayer;
    double t0error;
    double inislx;
    double iniinz;
    int    lidStart;
    int    lidStop;
};

// parameter block for calibration
class CalibPara{
public:
    CalibPara();
    virtual ~CalibPara(){};
    void Print();
};

// parameter block for XTAnalyzer
class XTAnalylzerPara{
public:
    XTAnalylzerPara();
    virtual ~XTAnalylzerPara(){};
    void Print();
};

// parameter block for Ana
class AnaPara{
public:
    AnaPara();
    virtual ~AnaPara(){};
    void Print();
};


/// This is a class to load and manage parameters for tracking purpose
/// There are several blocks for detailed purposes:
///      Tracking, Calib, XTAnalyzer, Ana
class ParameterManager {
public:
    enum ParaBlock{
        kTracking,
        kCalib,
        kXTAnalyzer,
        kAna,
        kAll
    };

    ParameterManager();
    virtual ~ParameterManager(){};

    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static ParameterManager& Get(void) {
        if (!fParameterManager)
            fParameterManager = new ParameterManager();
        return *fParameterManager;
    }

    /// Reads parameters from input files.
    /// Function can be used to read in extra parameter files.
    void ReadInputFile(TString filename, TString dirName="", bool tryFile = false, bool fixParameters = false);  

    /// Load parameters from a given
    void LoadParameters(ParaBlock theParaBlock = kAll);

    void Print();

    GeometryManager::GeoSetup geoSetup;
    GeometryManager::ChamberType chamberType;
    GeometryManager::ConnectionType connectionType;
    BeamManager::BeamType beamType;
    InputOutputManager::InputType inputType;
    TrackingPara::PeakType peakType;

    TrackingPara    TrackingParameters;
    CalibPara       CalibParameters;
    XTAnalylzerPara XTAnalylzerParameters;
    AnaPara         AnaParameters;

private:
    /// The static pointer to the singleton instance.
    static ParameterManager* fParameterManager;
};

#endif
