#ifndef Hit_hxx_seen
#define Hit_hxx_seen

class RawHit{
public:
    RawHit();
    virtual ~RawHit(){};

    int    LayerID;
    int    CellID;
    double DriftT;
    double Pedestal;
    int    ADCheight;
    double ADCsumPacket;
    double ADCsumAll;
    int    PacketWidth;
    int    nPeaksInChannel;
    int    iPeakInChannel;
    int    nPeaksInPacket;
    int    iPeakInPacket;
};

class ReconHit{
public:
    ReconHit();
    virtual ~ReconHit(){};

    int RawHitIndex;
    double DOCA;
};

class HitPair{
public:
    HitPair(){};
    virtual ~HitPair(){};

    int    iHit[2];
    int    layerID[2];
    double wireX[2];
    double wireY[2];
    double wireZ[2];
};

RawHit::RawHit():
    LayerID(-1),
    CellID(-1),
    DriftT(0),
    Pedestal(0),
    ADCheight(0),
    ADCsumPacket(0),
    ADCsumAll(0),
    PacketWidth(0),
    nPeaksInChannel(0),
    iPeakInChannel(0),
    nPeaksInPacket(0),
    iPeakInPacket(0)
{
}

ReconHit::ReconHit():
    RawHitIndex(-1),
    DOCA(0)
{
}

#endif
