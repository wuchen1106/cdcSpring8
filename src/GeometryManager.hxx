#ifndef GeometryManager_hxx_seen
#define GeometryManager_hxx_seen

#include "header.hxx"

class TString;
class Scintillator;
class Chamber;

class GeometryManager{
public:
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

    bool Initialize(GeoSetup theGeoSetup = kNormal);
    bool AdjustWirePosition(TString file);
    Scintillator * GetScintillator(){return fScintillator;};
    Chamber* GetChamber(){return fChamber;};
    void Print();

private:
    /// The static pointer to the singleton instance.
    static GeometryManager* fGeometryManager;

    GeoSetup       fGeoSetup;
    Scintillator * fScintillator;
    Chamber      * fChamber;
};

class Scintillator{
public:
    Scintillator(GeometryManager::GeoSetup theGeoSetup = GeometryManager::kNormal);
    virtual ~Scintillator(){};

    void Print();
    void SetGeometry(GeometryManager::GeoSetup theGeoSetup);
    bool IsInScinti(double saftyFactor,double inx, double slx, double inz, double slz, double y0 = YREFERENCE);

    double Yup;
    double Ydown;
    double HalfLength;
    double HalfWidth;
};

class Chamber{
public:
    Chamber(GeometryManager::GeoSetup theGeoSetup = GeometryManager::kNormal);
    virtual ~Chamber(){};

    void Print();
    void SetGeometry(GeometryManager::GeoSetup theGeoSetup);
    void Initialize();
    bool LoadWireMap(TString file);
    bool LoadCrossPoints(TString file);
    bool AdjustWirePosition(TString file);

    // map for wiremap
    double  map_x[NLAY][NCEL][2];
    double  map_y[NLAY][NCEL][2];
    double  map_z[NLAY][NCEL][2];
    int     map_ch[NLAY][NCEL];
    int     map_bid[NLAY][NCEL];
    double  map_theta[NLAY][NCEL];
    int     map_lid[NBRD][NCHS];
    int     map_wid[NBRD][NCHS];
    // map for wire position adjustment
    double  map_adjust[NLAY][NCEL];
    // map for cross points
    double mcp_xc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
    double mcp_zc[NZXP][NCEL][NCEL]; // z-x planes corresponding to the layerID of the lower layer counting from 1 
};

#endif
