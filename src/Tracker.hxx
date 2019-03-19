#ifndef Tracker_hxx_seen
#define Tracker_hxx_seen

#define NCAND 10

#include <vector>

#include "Track.hxx"

class Tracker{
public:
    Tracker();
    virtual ~Tracker();

    void Reset();
    void DoTracking();

    std::vector<std::vector<int>*> * hitLayerIndexMap;
    std::vector<int> * pairableLayers;

    TrackCandidate trackCandidates[NCAND];
};

#endif
