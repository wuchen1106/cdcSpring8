#ifndef Tracker_hxx_seen
#define Tracker_hxx_seen

#include <vector>
#include <map>

#include "Track.hxx"
#include "GeometryManager.hxx"
#include "InputOutputManager.hxx"

class TF1;
class TGraphErrors;

class Tracker{
public:
    Tracker(InputOutputManager::InputHitType theInputHitType);
    virtual ~Tracker();

    void Reset(); /// prepare for tracking. To be called every time given a new event
    void DoTracking(); /// the main function to get tracks from the current event with given hit list
    void SetOutput(); /// set the fitted results to InputOutputManager

    static TrackResult currentTrackResult;
    static TrackResult trackResults[NCAND];

    static std::vector<int> * hitIndexInTestLayer; /// store a list of hit indice in the test layer
    static std::vector<std::vector<int>*> * hitLayerIndexMap; /// store a list of hit indice for each layer (except for the test layer)
    static std::map <int, int>    hitIndexLeftRight; /// store a map from hit index to left right. Set upon l/r combination selection
    static std::map <int, double> hitIndexDriftDLeftMap; /// store a map from hit index to driftD calculated by XT left side (negative value)
    static std::map <int, double> hitIndexDriftDRightMap; /// store a map from hit index to driftD calculated by XT right side (positive value)

    int nGoodTracks; /// number of good tracks found in one event
    std::vector<int> * pairableLayers; /// a list of layers that may form a pair, i.e. having a non-empty neighbour layer (excluding test layer)
    int nPairs; /// number of pairs found in current pick & l/r combination
    int nGoodPairs; /// number of good pairs that is roughly on track. Updated after Fit2D is called.
    int pickIndex[NLAY]; /// the current picked hit index sorted in the sequence of pairable layers
    int pickLR[NLAY]; /// left right of the picked hits
    double pickWireY[NLAY]; /// wire position upon the picked hits
    double pairX[NLAY]; /// pair position (formed by one hit from the given layer and another from its upper layer)
    double pairY[NLAY]; /// pair position (formed by one hit from the given layer and another from its upper layer)
    double pairZ[NLAY]; /// pair position (formed by one hit from the given layer and another from its upper layer)

private:
    void updateDriftD(); /// calculate drift distance for every hit in the given maps with both left and right assumptions. Will be stored in hitIndexDriftDLeftMap and hitIndexDriftDRightMap
    int tracking(int ipick,int & iselection); /// Loop in the given layer hits map (hitLayerIndexMap). Called recursively. In each iterative call, pick up one hit per layer in the pairable layers and perform tracking.
    int fitting(int iselection); /// Fit the track with given selection of hits
    void setLeftRight(int icombi); /// get left/right from the given combination index
    void Reset2DFunctions(double MoveRatioX = 0, double MoveRatioZ = 0); /// get the 2-D fitting functions reset to default values. If arguements are given, then set with them as offsets
    bool Fit2D(double safetyFactor, bool fitWithError, double & chi2X, double & chi2Z, bool & inScint, bool & fromSource); /// do the 2-D fitting. Firstly the pair positions will be recalculated according to the track parameters stored in the 2-D functions. Then the 2-D functions will be updated with new fitting.
    int getChi2XZ(double & chi2x, double & chi2z); /// get chi2 for 2-D fittings on Y-X and Y-Z planes
    int updateWirePositionOnHit(); /// update wire positions upon picked hits
    int updatePairPositions(); /// update pair positions
    int setPairPositionGraphs(bool noError); /// set graphs containing pair position information; errors are given by calculating the difference between pair position and the current track parameters (stored in the 2-D fitting functions)
    void pickUpHitsForFitting(double slx, double inx, double slz, double inz, double residualCut); /// pick up hits from the given layer hit map. In each layer only choose one hit that is closest to the track. Abandon some layers if the residual of the closest one is still too larger than the residual cut
    void doFitting(double sliX, double iniX,double sliZ, double iniZ); /// The core part of track fitting
    static void fcn(int &npar, double *gin, double &f, double *par, int iflag); /// The function to be used by TMinuit
    static void getchi2(double &f, double & cp, double & ca, double slx, double inx, double slz, double inz, double t0offset,bool all); /// get the chi2 with 3-D track
    bool checkResults(int nHitsSel, int icombi, int iselection); /// Compare the new tracking result (written in currentTrackResult) with previously stored candidates (stored in trackResults) and put it in correct place if needed.

    TF1 * func_pairYX; /// 2-D fitting function on Y-X plane
    TF1 * func_pairYZ; /// 2-D fitting function on Y-Z plane
    TGraphErrors * graph_pairX; /// a graph to store pair positions with error on Y-X plane;
    TGraphErrors * graph_pairZ; /// a graph to store pair positions with error on Y-X plane

    InputOutputManager::InputHitType inputHitType;

    /// a list of parameters used for track fitting with TMinuit
    double arglist[10];
    int    ierflg;
    double amin;
    double edm;
    double errdef;
    int    nvpar;
    int    nparx;
    int    icstat;
    double slxStep;
    double slxMin;
    double slxMax;
    double slzStep;
    double slzMin;
    double slzMax;
    double inxStep;
    double inxMin;
    double inxMax;
    double inzStep;
    double inzMin;
    double inzMax;
};

#endif
