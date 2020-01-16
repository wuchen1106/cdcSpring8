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
    iSelection = -1;
    iCombination = -1;
    slopeX = 0;
    slopeZ = 0;
    interceptX = 0;
    interceptZ = 0;
    chi2 = 1e9;
    chi2a = 1e9;
    pValue = 0;
}

bool TrackCandidate::operator ==(const TrackCandidate & aCand){
    if (hitIndexSelected.size() != aCand.hitIndexSelected.size()) return false;
    for (std::size_t i = 0; i<hitIndexSelected.size(); i++){
        if (hitIndexSelected[i] != aCand.hitIndexSelected[i]) return false;
        if (hitLeftRightSelected[i] != aCand.hitLeftRightSelected[i]) return false;
    }
    return true;
}

Track2D::Track2D()
{
    Reset();
}

void Track2D::Reset(){
    TrackCandidate::Reset();
    nPairs = 0;
    chi2X = 1e9;
    chi2Z = 1e9;
}

Track2D& Track2D::operator =(const TrackCandidate& aCand){
    TrackCandidate::operator =(aCand);
    return *this;
}

Track3D::Track3D()
{
    Reset();
}

void Track3D::Reset(){
    TrackCandidate::Reset();
}

Track3D& Track3D::operator =(const TrackCandidate& aCand){
    TrackCandidate::operator =(aCand);
    return *this;
}
