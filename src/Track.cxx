#include "Track.hxx"

Track::Track():
    slopeX(0),
    slopeZ(0),
    interceptX(0),
    interceptZ(0)
{
}

TrackCandidate::TrackCandidate()
{
    Reset();
}

void TrackCandidate::Reset(){
    hitIndexSelected.clear();
    hitLeftRightSelected.clear();
    t0Offset = 0;
    nPairs = 0;
    nGoodPairs = 0;
    slopeX = 0;
    slopeZ = 0;
    interceptX = 0;
    interceptZ = 0;
    chi2X = 1e9;
    chi2Z = 1e9;
    chi2 = 1e9;
    chi2WithTestLayer = 1e9;
    pValue = 0;
    NDF = 0;
}

bool TrackCandidate::operator ==(const TrackCandidate & aCand){
    if (hitIndexSelected.size() != aCand.hitIndexSelected.size()) return false;
    for (int i = 0; i<hitIndexSelected.size(); i++){
        if (hitIndexSelected[i] != aCand.hitIndexSelected[i]) return false;
        if (hitLeftRightSelected[i] != aCand.hitLeftRightSelected[i]) return false;
    }
    return true;
}

void TrackCandidate::operator =(const TrackCandidate & aCand){
    hitIndexSelected = aCand.hitIndexSelected;
    hitLeftRightSelected = aCand.hitLeftRightSelected;
    t0Offset = aCand.t0Offset;
    nPairs = aCand.nPairs;
    nGoodPairs = aCand.nGoodPairs;
    chi2X = aCand.chi2X;
    chi2Z = aCand.chi2Z;
    slopeX = aCand.slopeX;
    slopeZ = aCand.slopeZ;
    interceptX = aCand.interceptX;
    interceptZ = aCand.interceptZ;
}

TrackResult::TrackResult()
{
    Reset();
}

void TrackResult::Reset(){
    hitIndexSelected.clear();
    hitLeftRightSelected.clear();
    slopeX = 0;
    slopeZ = 0;
    interceptX = 0;
    interceptZ = 0;
    chi2X = 1e9;
    chi2Z = 1e9;
    chi2 = 1e9;
    chi2WithTestLayer = 1e9;
    pValue = 0;
    NDF = 0;
    initialTrackCandidate.Reset();
}

bool TrackResult::operator ==(const TrackResult & aResult){
    if (hitIndexSelected.size() != aResult.hitIndexSelected.size()) return false;
    for (int i = 0; i<hitIndexSelected.size(); i++){
        if (hitIndexSelected[i] != aResult.hitIndexSelected[i]) return false;
        if (hitLeftRightSelected[i] != aResult.hitLeftRightSelected[i]) return false;
    }
    return true;
}

void TrackResult::operator =(const TrackResult & aResult){
    hitIndexSelected = aResult.hitIndexSelected;
    hitLeftRightSelected = aResult.hitLeftRightSelected;
    slopeX = aResult.slopeX;
    slopeZ = aResult.slopeZ;
    interceptX = aResult.interceptX;
    interceptZ = aResult.interceptZ;
    chi2 = aResult.chi2;
    chi2WithTestLayer = aResult.chi2WithTestLayer;
    pValue = aResult.pValue;
    NDF = aResult.NDF;
    initialTrackCandidate = aResult.initialTrackCandidate;
}
