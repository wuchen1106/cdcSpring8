#ifndef ParameterManager_hxx_seen
#define ParameterManager_hxx_seen

#include <TString.h>
#include "GeometryManager.hxx"
#include "BeamManager.hxx"
#include "InputOutputManager.hxx"

#define NRANGES 10

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
    int    nHitsGMax;
    int    nPairsMin;
    int    t0shift0;
    int    t0shift1;
    int    tmin;
    int    tmax;
    double sumCut;
    double aaCut;
    PeakType peakType;
    double t0error;
    double inislx;
    double iniinz;
    int    lidStart;
    int    lidStop;
    int    BlindLayer;
};

// parameter block for calibration
class CalibPara{
public:
    CalibPara();
    virtual ~CalibPara(){};
    void Print();
};

// parameter block for XTManager
class XTManagerPara{
public:
    XTManagerPara();
    virtual ~XTManagerPara(){};
    void Print();
    int  xtType;
    int  defaultLayer;
    int  evenLayer;
    int  oddLayer;
};

// parameter block for XTAnalyzer
class XTAnalyzerPara{
public:
    XTAnalyzerPara();
    virtual ~XTAnalyzerPara(){};
    void Print();
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
    // binning of the 2D histograms
    double bin_t_min;
    double bin_t_max;
    int    bin_t_num;
    double bin_x_min;
    double bin_x_max;
    int    bin_x_num;
    // regroup T bins to fit X
    int    fitX_nRanges;
    double fitX_tSep[NRANGES+1];
    int    fitX_minEntries[NRANGES];
    int    fitX_nBins[NRANGES];
    int    fitX_smooth[NRANGES];
    bool   fitX_fitBoth[NRANGES];
    bool   fitX_SetEmptyBins[NRANGES];
    // about forming graphs
    int    graph_n_min;
    double graph_chi2_max;
    double graph_prob_min;
    double graph_sepX;
    // about XT function
    int    xtfunc_nRanges;
    double xtfunc_tHighEdge;
    double xtfunc_tLeft[NRANGES];
    double xtfunc_tRight[NRANGES];
    double xtfunc_tLowEdge[NRANGES];
    int    xtfunc_nPol[NRANGES];
    // about drawing samples and fitting results
    double draw_tmin;
    double draw_tmax;
    double draw_xmin;
    double draw_xmax;
};

// parameter block for HistogramAnalyzer
class HistogramAnalyzerPara{
public:
    HistogramAnalyzerPara();
    virtual ~HistogramAnalyzerPara(){};
    void Print();
    // regroup T bins to fit X
    int    functionType[NRANGES];
    double peak_height_middle[NRANGES];
    double peak_height_left[NRANGES];
    double peak_height_right[NRANGES];
    double peak_sigma_middle[NRANGES];
    double peak_sigma_left[NRANGES];
    double peak_sigma_right[NRANGES];
    double peak_mean_range[NRANGES];
    double base_height_middle[NRANGES];
    double base_height_left[NRANGES];
    double base_height_right[NRANGES];
    double base_sigma_middle[NRANGES];
    double base_sigma_left[NRANGES];
    double base_sigma_right[NRANGES];
    double base_mean_range[NRANGES];
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
///      Tracking, Calib, XTManager, XTAnalyzer, HistogramAnalyzer, Ana
class ParameterManager {
public:
    enum ParaBlock{
        kTracking,
        kCalib,
        kXTManager,
        kXTAnalyzer,
        kHistogramAnalyzer,
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
    // TODO: either name them as hitFileType & trackFileType, or add output* values to separately control input and output data types
    InputOutputManager::DataType inputHitType;
    InputOutputManager::DataType inputTrackType;
    TrackingPara::PeakType peakType;

    TrackingPara    TrackingParameters;
    CalibPara       CalibParameters;
    XTManagerPara   XTManagerParameters;
    XTAnalyzerPara  XTAnalyzerParameters;
    HistogramAnalyzerPara HistogramAnalyzerParameters;
    AnaPara         AnaParameters;

private:
    int getFunctionType(TString name);
    int getXTType(TString name);

    /// The static pointer to the singleton instance.
    static ParameterManager* fParameterManager;
};

#endif
