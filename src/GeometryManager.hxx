#ifndef GeometryManager_hxx_seen
#define GeometryManager_hxx_seen

#include <TVector3.h>

// Single board parameters
#define NSAM 32
#define NCHS 48
#define MIN_ADCL -50  // used in getH
#define MAX_ADCL 1000 // used in getH
#define MIN_ADC 50   // used in getP
#define MAX_ADC 750  // used in getP
#define NBINS  256
#define MAX_WIDTH 10 // about merging. Depending on drift velocity: HV and Gas

// Choose number of layers
#define NLAY 9
//#define NLAY  11
//#define NLAY  19

// TODO: consider to make a large map for general cases
#if NLAY == 9
#define NCEL  11
#define NCELA 99
#define NBRD 2
#define NCHT 96
#define NZXP 8
#define NITERSMAX 100
#elif NLAY == 11
#define NCEL  12
#define NCELA 132
#define NBRD 3
#define NCHT 144
#define NZXP 10
#define NITERSMAX 100
#elif NLAY == 19 // CRT
#define NCEL  12
#define NCELA 228
#define NBRD 13
#define NCHT 624
#define NZXP 18
#define NITERSMAX 25
#endif

class TString;
class Scintillator;
class Chamber;

class GeometryManager{
public:
    enum ChamberType{
        kProto4,
        kCDC
    };

    enum ConnectionType{
        kSPring8,
        kCosmic
    };

    enum GeoSetup{
        kNormal,
        kTilted,
        kFinger
    };

    GeometryManager();
    virtual ~GeometryManager(){};

    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static GeometryManager& Get(void) {
        if (!fGeometryManager)
            fGeometryManager = new GeometryManager();
        return *fGeometryManager;
    }

    bool IsInScinti(double saftyFactor,double inx, double slx, double inz, double slz);
    bool Initialize(GeoSetup theGeoSetup = kNormal, ConnectionType theConnectionType = kSPring8, ChamberType theChamberType = kProto4);
    bool AdjustWirePosition(TString file);
    Scintillator * GetScintillator(){return fScintillator;};
    Chamber* GetChamber(){return fChamber;};
    double GetDOCA(int lid, int wid, double slx, double inx, double slz, double inz);
    void Print();

    // the refernece X-Z plane position
    // TODO: for real chamber case we need reference radius
    double ReferenceY;

    GeoSetup       fGeoSetup;
    Scintillator * fScintillator;
    Chamber      * fChamber;

private:
    /// The static pointer to the singleton instance.
    static GeometryManager* fGeometryManager;

    TVector3 vTrackU;
    TVector3 vTrackD;
    TVector3 vTrack;
    TVector3 vWireHV;
    TVector3 vWireRO;
    TVector3 vWire;
    TVector3 vDist;
    TVector3 vAxis;
};

class Scintillator{
public:
    Scintillator(GeometryManager::GeoSetup theGeoSetup = GeometryManager::kNormal);
    virtual ~Scintillator(){};

    void Print();
    void SetGeometry(GeometryManager::GeoSetup theGeoSetup);

    double Yup;
    double Ydown;
    double HalfLength;
    double HalfWidth;
};

class Chamber{
public:
    Chamber(GeometryManager::GeoSetup theGeoSetup = GeometryManager::kNormal, GeometryManager::ConnectionType theConnectionType = GeometryManager::kSPring8, GeometryManager::ChamberType theChamberType = GeometryManager::kProto4);
    virtual ~Chamber(){};

    void Print();
    void SetGeometry(GeometryManager::GeoSetup theGeoSetup, GeometryManager::ConnectionType theConnectionType, GeometryManager::ChamberType theChamberType);
    void Initialize();
    bool LoadWireMap(TString file);
    bool LoadCrossPoints(TString file);
    bool AdjustWirePosition(TString file);

    // Cell shape
    double cellHeight;
    double cellWidth;

    // chamber geometry
    double chamberLength;
    double chamberHeight;
    double chamberPositionX;
    double chamberPositionY;
    double chamberPositionZ;

    // map for wires
    double  wire_x[NLAY][NCEL][2];
    double  wire_y[NLAY][NCEL][2];
    double  wire_z[NLAY][NCEL][2];
    int     wire_ch[NLAY][NCEL];
    int     wire_bid[NLAY][NCEL];
    double  wire_theta[NLAY][NCEL];
    int     wire_lid[NBRD][NCHS];
    int     wire_wid[NBRD][NCHS];
    // map for wire position adjustment
    // TODO: consider about wire sag and wire rotation calibration in the future
    double  wire_adjustX[NLAY][NCEL];
    double  wire_adjustY[NLAY][NCEL];
    double  wire_adjustZ[NLAY][NCEL];
    // map for cross points of two wire with the two given cellIDs from the given layer and the layer above it.
    double wirecross_x[NLAY][NCEL][NCEL];
    double wirecross_y[NLAY][NCEL][NCEL];
    double wirecross_z[NLAY][NCEL][NCEL];
};

#endif
