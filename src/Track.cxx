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
}

void TrackCandidate::Reset(){
    hitIndexSelected.clear();
    slopeX = 0;
    slopeZ = 0;
    interceptX = 0;
    interceptZ = 0;
}

TrackResult::TrackResult():
    chi2(0),
    NDF(0)
{
}

void TrackResult::Reset(){
}
