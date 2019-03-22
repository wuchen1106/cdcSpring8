#include "GeometryManager.hxx"
#include "Tracker.hxx"

Tracker::Tracker():
    hitLayerIndexMap(0),
    pairableLayers(0)
{
    hitLayerIndexMap = new std::vector<std::vector<int>*>;
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        hitLayerIndexMap->push_back(new std::vector<int>);
    }
    pairableLayers = new std::vector<int>;
}

Tracker::~Tracker(){
    delete hitLayerIndexMap;
    delete pairableLayers;
}

void Tracker::Reset(){
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        hitLayerIndexMap->at(iLayer)->clear();
    }
    pairableLayers->clear();
    for (int iCand = 0; iCand<NCAND; iCand++){
        trackCandidates[iCand].Reset();
    }
}

void Tracker::DoTracking(){
    // TODO: add tracking stuffs here
}
