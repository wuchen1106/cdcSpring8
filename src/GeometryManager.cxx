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

bool Scintillator::IsInScinti(double saftyFactor,double inx, double slx, double inz, double slz, double y0){
    double xtop = 1/saftyFactor*(inx+slx*(Yup-y0));
    double xbot = 1/saftyFactor*(inx+slx*(Ydown-y0));
    double ztop = 1/saftyFactor*(inz+slz*(Yup-y0));
    double zbot = 1/saftyFactor*(inz+slz*(Ydown-y0));
    if (xtop>HalfWidth||xtop<-HalfWidth||xbot>HalfWidth||xbot<-HalfWidth||ztop>HalfLength||ztop<-HalfLength||zbot>HalfLength||zbot<-HalfLength) return false;
    else return true;
}

void Scintillator::Print(){
    printf("  Scintillator geometry:\n");
    printf("    Yup        = %.3e\n",Yup);
    printf("    Ydown      = %.3e\n",Ydown);
    printf("    HalfLength = %.3e\n",HalfLength);
    printf("    HalfWidth  = %.3e\n",HalfWidth);
}

Chamber::Chamber(GeometryManager::GeoSetup theGeoSetup)
{
    SetGeometry(theGeoSetup);
}

void Chamber::SetGeometry(GeometryManager::GeoSetup theGeoSetup){
}

void Chamber::Initialize(){
    for(int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            map_x[lid][wid][0] = 0;
            map_y[lid][wid][0] = 0;
            map_z[lid][wid][0] = 0;
            map_x[lid][wid][1] = 0;
            map_y[lid][wid][1] = 0;
            map_z[lid][wid][1] = 0;
        	map_ch[lid][wid] = -1;
        	map_bid[lid][wid] = -1;
        	map_adjust[lid][wid] = 0;
            if (lid <NZXP){ // z-x planes corresponding to the layerID of the lower layer counting from 1 
                for (int wjd = 0; wjd<NCEL; wjd++){
                    mcp_xc[lid][wid][wjd] = 999;
                    mcp_zc[lid][wid][wjd] = 999;
                }
            }
        }
    }
}

bool Chamber::LoadWireMap(TString file){
    TChain * iChain = new TChain("t");
    iChain->Add(file);
    if (!iChain->GetEntries()){
        return false;
    }
    int     wp_bid;
    int     wp_ch;
    int     wp_wid;
    int     wp_lid;
    double  wp_xro;
    double  wp_yro;
    double  wp_xhv;
    double  wp_yhv;
    iChain->SetBranchAddress("b",&wp_bid);
    iChain->SetBranchAddress("ch",&wp_ch);
    iChain->SetBranchAddress("l",&wp_lid);
    iChain->SetBranchAddress("w",&wp_wid);
    iChain->SetBranchAddress("xhv",&wp_xhv);
    iChain->SetBranchAddress("yhv",&wp_yhv);
    iChain->SetBranchAddress("xro",&wp_xro);
    iChain->SetBranchAddress("yro",&wp_yro);
    for (int i = 0; i<iChain->GetEntries(); i++){
        iChain->GetEntry(i);
        if (wp_lid>=0&&wp_lid<NLAY&&wp_wid>=0&&wp_wid<NCEL){
            map_x[wp_lid][wp_wid][0] = wp_xhv+map_adjust[wp_lid][wp_wid];
            map_y[wp_lid][wp_wid][0] = wp_yhv;
            map_z[wp_lid][wp_wid][0] = -chamberHL;
            map_x[wp_lid][wp_wid][1] = wp_xro+map_adjust[wp_lid][wp_wid];
            map_y[wp_lid][wp_wid][1] = wp_yro;
            map_z[wp_lid][wp_wid][1] = chamberHL;
            map_ch[wp_lid][wp_wid] = wp_ch;
            map_bid[wp_lid][wp_wid] = wp_bid;
            map_theta[wp_lid][wp_wid] = atan(-(wp_xhv-wp_xro)/chamberHL/2); // rotation angle w.r.t the dart plane: read out  plane; positive rotation angle point to -x direction
            MyNamedInfo("WireMap",Form("map_theta[%d][%d] = atan(-(%.3e-%.3e)/%.3e/2) = %.3e",wp_lid,wp_wid,wp_xhv,wp_xro,chamberHL,map_theta[wp_lid][wp_wid]));
        }
        if (wp_bid>=0&&wp_bid<NBRD&&wp_ch>=0&&wp_ch<NCHS){
            map_lid[wp_bid][wp_ch] = wp_lid;
            map_wid[wp_bid][wp_ch] = wp_wid;
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
    int     cp_l1;
    int     cp_l2;
    int     cp_w1;
    int     cp_w2;
    double  cp_zc;
    double  cp_xc;
    iChain->SetBranchAddress("l1",&cp_l1);
    iChain->SetBranchAddress("l2",&cp_l2);
    iChain->SetBranchAddress("w1",&cp_w1);
    iChain->SetBranchAddress("w2",&cp_w2);
    iChain->SetBranchAddress("z",&cp_zc);
    iChain->SetBranchAddress("x",&cp_xc);
    int nEntries = iChain->GetEntries();
    for (int iEntry = 0; iEntry<nEntries; iEntry++){
        iChain->GetEntry(iEntry);
        if (cp_l1>=0&&cp_l1<NLAY&&cp_w1>=0&&cp_w1<NCEL&&cp_l2>=0&&cp_l2<NLAY&&cp_w2>=0&&cp_w2<NCEL){
            mcp_xc[cp_l1][cp_w1][cp_w2] = cp_xc;
            mcp_zc[cp_l1][cp_w1][cp_w2] = cp_zc;
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
    double off_adjustment = 0;
    int off_lid = 0;
    int off_wid = 0;
    iChain->SetBranchAddress("adjust",&off_adjustment);
    iChain->SetBranchAddress("lid",&off_lid);
    iChain->SetBranchAddress("wid",&off_wid);
    for (int iEntry = 0; iEntry<iChain->GetEntries(); iEntry++){
        iChain->GetEntry(iEntry);
        map_adjust[off_lid][off_wid] = off_adjustment;
    }
    return true;
}

void Chamber::Print(){
    printf("  Chamber geometry:\n");
}

GeometryManager* GeometryManager::fGeometryManager = NULL;

GeometryManager::GeometryManager():
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

bool GeometryManager::Initialize(GeometryManager::GeoSetup theGeoSetup){
    fGeoSetup = theGeoSetup;
    fScintillator->SetGeometry(theGeoSetup);
    fChamber->SetGeometry(theGeoSetup);
    fChamber->Initialize();
    TString HOME=getenv("CDCS8WORKING_DIR");
    TString file = HOME+"/Input/wire-position.root";
    if (!fChamber->LoadWireMap(file)) return false;
    file = HOME+"/Input/crosspoint.root";
    if (!fChamber->LoadCrossPoints(file)) return false;

    return true;
}

bool GeometryManager::AdjustWirePosition(TString file){
    return fChamber->AdjustWirePosition(file);
}
