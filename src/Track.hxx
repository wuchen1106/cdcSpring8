#ifndef Track_hxx_seen
#define Track_hxx_seen

#include <vector>

class Track{
public:
    Track();
    virtual ~Track(){};

    virtual void Reset(){
        slopeX = 0;
        slopeZ = 0;
        interceptX = 0;
        interceptZ = 0;
    };

    virtual int nPars(void) {return 4;};

    double slopeX;
    double slopeZ;
    double interceptX;
    double interceptZ;
};

class TrackCandidate:public Track{
public:
    TrackCandidate();
    virtual ~TrackCandidate(){};
    virtual void Reset();

    bool operator ==(const TrackCandidate &);

    std::vector<int> hitIndexSelected;
    std::vector<int> hitLeftRightSelected;
    int    iSelection;
    int    iCombination;
    double t0Offset;
    double chi2;
    double chi2a;
    double chi2WithTestLayer;
    double pValue;
    double NDF;
};

class Track2D:public TrackCandidate{
public:
    Track2D();
    virtual ~Track2D(){};
    virtual void Reset();

    Track2D& operator =(const TrackCandidate &);

    virtual int nPars(void) {return 2;};

    int    nPairs;
    double chi2X;
    double chi2Z;
};

class Track3D:public TrackCandidate{
public:
    Track3D();
    virtual ~Track3D(){};
    virtual void Reset();

    Track3D& operator =(const TrackCandidate &);

    virtual int nPars(void) {return 4;};
};

#endif
