#ifndef BeamManager_hxx_seen
#define BeamManager_hxx_seen

class BeamManager{
public:
    enum BeamType{
        kSPring8,
        kSPring8Tilted,
        kCosmic
    };

    BeamManager();
    virtual ~BeamManager(){};

    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static BeamManager& Get(void) {
        if (!fBeamManager)
            fBeamManager = new BeamManager();
        return *fBeamManager;
    }

    bool Initialize(BeamType theBeamType);
    bool IsInBeam();
    void Print();

private:
    /// The static pointer to the singleton instance.
    static BeamManager* fBeamManager;

    BeamType beamType;
    double   beamSlz;
    double   beamSlx;
    double   beamInz;
    double   beamInx;
    double   beamSlzRange;
    double   beamSlxRange;
    double   beamInzRange;
    double   beamInxRange;
};

#endif
