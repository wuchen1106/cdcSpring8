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

    int    nPars(void){return 4;};

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
    bool operator ==(const TrackCandidate &);
    void operator =(const TrackCandidate &);

    std::vector<int> hitIndexSelected;
    std::vector<int> hitLeftRightSelected;

    double t0Offset;
    int    nPairs;
    int    nGoodPairs;
    int    iSelection;
    int    iCombination;
    double chi2X;
    double chi2Z;
    double chi2;
    double chi2WithTestLayer;
    double pValue;
    double NDF;
};

class TrackResult:public TrackCandidate{
public:
    TrackResult();
    virtual ~TrackResult(){};

    void Reset();
    bool operator ==(const TrackResult &);
    void operator =(const TrackResult &);

    TrackCandidate initialTrackCandidate;
};

#endif
