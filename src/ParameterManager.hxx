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
class XTAnalyzerPara{
public:
    XTAnalyzerPara();
    virtual ~XTAnalyzerPara(){};
    void Print();
    int XTType;
    bool AsymXT;
    TString CandSelBy;
    bool RequireInTriggerCounter;
    bool RequireAllGoldenHits;
    bool FirstGoodPeak;
    bool UseGoodHit;
    bool AllGoodHitsUsed;
    int nHits_max;
    int nHitsS_min;
    double chi2_max;
    double pValue_min;
    double slz_min;
    double slz_max;
    double gold_t_min;
    double gold_t_max;
    // about minimum entries to apply fitting
    int    bin_n_min;
    // about binning
    double bin_t_min;
    double bin_t_max;
    int    bin_t_num;
    double bin_x_min;
    double bin_x_max;
    int    bin_x_num;
    // about projection
    int    bin_t_fit_num;
    double bin_t_tailTime;
    int    bin_t_fit_num_tail;
    int    bin_x_fit_num;
    // about the range for using Landau function
    double bin_t_landTmin;
    double bin_t_landTmax;
    double bin_x_landXmin;
    double bin_x_landXmax;
    // about fitting range
    double bin_t_ratio;
    double bin_x_ratio;
    // about forming graphs
    int    graph_n_min;
    double graph_chi2_max;
    double graph_sepX;
    // about XT function
    int    xt_center_nPol;
    double xt_center_tLeft;
    double xt_center_tRight;
    int    xt_middle_nPol;
    double xt_middle_tLeft;
    double xt_middle_tRight;
    int    xt_end_nPol;
    double xt_end_tLeft;
    double xt_end_tRight;
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
    InputOutputManager::InputHitType inputHitType;
    TrackingPara::PeakType peakType;

    TrackingPara    TrackingParameters;
    CalibPara       CalibParameters;
    XTAnalyzerPara XTAnalyzerParameters;
    AnaPara         AnaParameters;

private:
    /// The static pointer to the singleton instance.
    static ParameterManager* fParameterManager;
};

#endif
