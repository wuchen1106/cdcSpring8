#ifndef Track_hxx_seen
#define Track_hxx_seen

#include <vector>

class Track{
public:
    Track();
    virtual ~Track(){};

    void Reset(){
        slopeX = 0;
        slopeZ = 0;
        interceptX = 0;
        interceptZ = 0;
    };

    double slopeX;
    double slopeZ;
    double interceptX;
    double interceptZ;
};

class TrackCandidate:public Track{
public:
    TrackCandidate();
    virtual ~TrackCandidate(){};

    void Reset();

    std::vector<int> hitIndexSelected;
};

class TrackResult:public TrackCandidate{
public:
    TrackResult();
    virtual ~TrackResult(){};

    void Reset();

    double chi2;
    double NDF;
    TrackResult * fInitialTrackResult;
};

#endif
