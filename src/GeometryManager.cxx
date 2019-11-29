#include <TChain.h>
#include <TString.h>

#include "GeometryManager.hxx"
#include "Log.hxx"

Scintillator::Scintillator(GeometryManager::GeoSetup theGeoSetup):
    Yup(0),
    Ydown(0),
    HalfLength(0),
    HalfWidth(0)
{
    SetGeometry(theGeoSetup);
}

void Scintillator::SetGeometry(GeometryManager::GeoSetup theGeoSetup){
    double chamberHH = 170.05/2;
    double chamberCY = 572;
    if (theGeoSetup == GeometryManager::kNormal){
        // normal scintillator
        Yup = chamberCY+chamberHH+180; // mm
        Ydown = chamberCY-chamberHH-180; 
        HalfLength = 300/2.;
        HalfWidth = 90/2.;
    }
    else if (theGeoSetup == GeometryManager::kTilted){
        // normal scintillator for tilted run
        // TODO: position should be modified
        Yup = chamberCY+chamberHH+180; // mm
        Ydown = chamberCY-chamberHH-180; 
        HalfLength = 3000/2.;
        HalfWidth = 9000/2.;
    }
    else if (theGeoSetup == GeometryManager::kFinger){
        // finger scintillator
        Yup = chamberCY+chamberHH+250; // mm
        Ydown = chamberCY-chamberHH-195; 
        HalfLength = 33/2.;
        HalfWidth = 33/2.;
    }
}

void Scintillator::Print(){
    printf("  Scintillator geometry:\n");
    printf("    Yup        = %.3e\n",Yup);
    printf("    Ydown      = %.3e\n",Ydown);
    printf("    HalfLength = %.3e\n",HalfLength);
    printf("    HalfWidth  = %.3e\n",HalfWidth);
}

Chamber::Chamber(GeometryManager::GeoSetup theGeoSetup, GeometryManager::ConnectionType theConnectionType, GeometryManager::ChamberType theChamberType):
    cellHeight(0),
    cellWidth(0),
    chamberLength(0),
    chamberHeight(0),
    chamberPositionX(0),
    chamberPositionY(0),
    chamberPositionZ(0)
{
    SetGeometry(theGeoSetup, theConnectionType, theChamberType);
}

void Chamber::SetGeometry(GeometryManager::GeoSetup theGeoSetup, GeometryManager::ConnectionType theConnectionType, GeometryManager::ChamberType theChamberType){
    cellHeight = 16;
    cellWidth = 16.8;
    if (theChamberType == GeometryManager::kProto4){
        chamberLength = 599.17;
        chamberHeight = 170.05;
        chamberPositionX = 0;
        chamberPositionY = 572;
        chamberPositionZ = 0;
    }
}

void Chamber::Initialize(){
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            wire_x[lid][wid][0] = 0;
            wire_y[lid][wid][0] = 0;
            wire_z[lid][wid][0] = 0;
            wire_x[lid][wid][1] = 0;
            wire_y[lid][wid][1] = 0;
            wire_z[lid][wid][1] = 0;
            wire_ch[lid][wid] = -1;
            wire_bid[lid][wid] = -1;
            wire_adjustX[lid][wid] = 0;
            wire_adjustY[lid][wid] = 0;
            wire_adjustZ[lid][wid] = 0;
            for (int wjd = 0; wjd<NCEL; wjd++){
                wirecross_x[lid][wid][wjd] = 999;
                wirecross_z[lid][wid][wjd] = 999;
            }
        }
    }
}

bool Chamber::LoadWireMap(TString file){
    // TODO: this wire map is for proto type alone. Should consider the full wire map format later
    TChain * iChain = new TChain("t");
    iChain->Add(file);
    if (!iChain->GetEntries()){
        return false;
    }
    int     bid;
    int     ch;
    int     wid;
    int     lid;
    double  xro;
    double  yro;
    double  xhv;
    double  yhv;
    iChain->SetBranchAddress("b",&bid);
    iChain->SetBranchAddress("ch",&ch);
    iChain->SetBranchAddress("l",&lid);
    iChain->SetBranchAddress("w",&wid);
    iChain->SetBranchAddress("xhv",&xhv);
    iChain->SetBranchAddress("yhv",&yhv);
    iChain->SetBranchAddress("xro",&xro);
    iChain->SetBranchAddress("yro",&yro);
    for (Long64_t i = 0; i<iChain->GetEntries(); i++){
        iChain->GetEntry(i);
        if (lid>=0&&lid<NLAY&&wid>=0&&wid<NCEL){
            wire_x[lid][wid][0] = xhv+wire_adjustX[lid][wid];
            wire_y[lid][wid][0] = yhv+wire_adjustY[lid][wid];
            wire_z[lid][wid][0] = -chamberLength/2+wire_adjustZ[lid][wid];
            wire_x[lid][wid][1] = xro+wire_adjustX[lid][wid];
            wire_y[lid][wid][1] = yro+wire_adjustY[lid][wid];
            wire_z[lid][wid][1] = chamberLength/2+wire_adjustZ[lid][wid];
            wire_ch[lid][wid] = ch;
            wire_bid[lid][wid] = bid;
            wire_theta[lid][wid] = atan(-(xhv-xro)/chamberLength); // rotation angle w.r.t the dart plane: read out  plane; positive rotation angle point to -x direction
            MyNamedInfo("WireMap",Form("wire_theta[%d][%d] = atan(-(%.3e-%.3e)/%.3e) = %.3e",lid,wid,xhv,xro,chamberLength,wire_theta[lid][wid]));
        }
        if (bid>=0&&bid<NBRD&&ch>=0&&ch<NCHS){
            wire_lid[bid][ch] = lid;
            wire_wid[bid][ch] = wid;
        }
    }
    return true;
}

bool Chamber::LoadCrossPoints(TString file){
    TChain * iChain = new TChain("t");
    iChain->Add(file);
    if (!iChain->GetEntries()){
        return false;
    }
    int     l1;
    int     l2;
    int     w1;
    int     w2;
    double  zc;
    double  xc;
    iChain->SetBranchAddress("l1",&l1);
    iChain->SetBranchAddress("l2",&l2);
    iChain->SetBranchAddress("w1",&w1);
    iChain->SetBranchAddress("w2",&w2);
    iChain->SetBranchAddress("z",&zc);
    iChain->SetBranchAddress("x",&xc);
    int nEntries = iChain->GetEntries();
    for (Long64_t iEntry = 0; iEntry<nEntries; iEntry++){
        iChain->GetEntry(iEntry);
        if (l1>=0&&l1<NLAY&&w1>=0&&w1<NCEL&&l2>=0&&l2<NLAY&&w2>=0&&w2<NCEL){
            wirecross_x[l1][w1][w2] = xc;
            wirecross_z[l1][w1][w2] = zc;
        }
    }
    return true;
}

bool Chamber::AdjustWirePosition(TString file){
    //===================get old wire map adjustment============================
    TChain * iChain = new TChain("t");
    iChain->Add(file);
    if (iChain->GetEntries()){
        return false;
    }
    // TODO: add the support to more dimensions, and wire sag + wire rotation?
    double off_adjustment = 0;
    int off_lid = 0;
    int off_wid = 0;
    iChain->SetBranchAddress("adjust",&off_adjustment);
    iChain->SetBranchAddress("lid",&off_lid);
    iChain->SetBranchAddress("wid",&off_wid);
    for (Long64_t iEntry = 0; iEntry<iChain->GetEntries(); iEntry++){
        iChain->GetEntry(iEntry);
        wire_adjustX[off_lid][off_wid] = off_adjustment;
    }
    return true;
}

void Chamber::Print(){
    printf("  Chamber geometry:\n");
}

GeometryManager* GeometryManager::fGeometryManager = NULL;

GeometryManager::GeometryManager():
    ReferenceY(0),
    fScintillator(0),
    fChamber(0)
{
    fGeoSetup = kNormal;
    fScintillator = new Scintillator();
    fChamber = new Chamber();
}

void GeometryManager::Print(){
    printf("#########################GeometryManager###########################\n");
    printf("  geometry setup %d:\n",fGeoSetup);
    fScintillator->Print();
    fChamber->Print();
}

bool GeometryManager::Initialize(GeoSetup theGeoSetup, ConnectionType theConnectionType, ChamberType theChamberType){
    fGeoSetup = theGeoSetup;
    fScintillator->SetGeometry(theGeoSetup);
    fChamber->SetGeometry(theGeoSetup, theConnectionType,theChamberType);
    fChamber->Initialize();
    ReferenceY = fChamber->chamberPositionY+fChamber->chamberHeight/2.;
    TString HOME=getenv("CDCS8WORKING_DIR");
    // TODO: add support for other wire maps according to chamber type and connection type
    TString file = HOME+"/Input/wire-position.root";
    if (!fChamber->LoadWireMap(file)) return false;
    file = HOME+"/Input/crosspoint.root";
    if (!fChamber->LoadCrossPoints(file)) return false;

    return true;
}

bool GeometryManager::AdjustWirePosition(TString file){
    return fChamber->AdjustWirePosition(file);
}

bool GeometryManager::IsInScinti(double saftyFactor,double inx, double slx, double inz, double slz){
    double xtop = 1/saftyFactor*(inx+slx*(fScintillator->Yup-ReferenceY));
    double xbot = 1/saftyFactor*(inx+slx*(fScintillator->Ydown-ReferenceY));
    double ztop = 1/saftyFactor*(inz+slz*(fScintillator->Yup-ReferenceY));
    double zbot = 1/saftyFactor*(inz+slz*(fScintillator->Ydown-ReferenceY));
    double HalfWidth = fScintillator->HalfWidth;
    double HalfLength = fScintillator->HalfLength;
    if (xtop>HalfWidth||xtop<-HalfWidth||xbot>HalfWidth||xbot<-HalfWidth||ztop>HalfLength||ztop<-HalfLength||zbot>HalfLength||zbot<-HalfLength) return false;
    else return true;
}

double GeometryManager::GetDOCA(int lid, int wid, double slx, double inx, double slz, double inz)
{
    double ydown = fScintillator->Ydown;
    double xdown = inx-slx*(ReferenceY-ydown);
    double zdown = inz-slz*(ReferenceY-ydown);
    vTrackU.SetXYZ(inx,ReferenceY,inz);
    vTrackD.SetXYZ(xdown,ydown,zdown);
    vWireHV.SetXYZ(fChamber->wire_x[lid][wid][0],fChamber->wire_y[lid][wid][0],fChamber->wire_z[lid][wid][0]);
    vWireRO.SetXYZ(fChamber->wire_x[lid][wid][1],fChamber->wire_y[lid][wid][1],fChamber->wire_z[lid][wid][1]);
    vTrack = vTrackD-vTrackU;
    vWire = vWireRO-vWireHV;
    vDist = vWireHV-vTrackU;
    vAxis = vWire.Cross(vTrack);
    double value = -vDist*(vAxis.Unit());
    return value;
}
