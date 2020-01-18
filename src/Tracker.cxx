#include <TGraphErrors.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TBackCompFitter.h>

#include "GeometryManager.hxx"
#include "ParameterManager.hxx"
#include "Tracker.hxx"
#include "InputOutputManager.hxx"
#include "BeamManager.hxx"
#include "XTManager.hxx"
#include "Log.hxx"

std::vector<int>               * Tracker::hitIndexInTestLayer = NULL;
std::vector<std::vector<int>*> * Tracker::hitLayerIndexMap = NULL;
Track2D Tracker::currentTrack2D;
Track3D Tracker::currentTrack3D;
Track2D Tracker::track2Ds[NCAND];
Track3D Tracker::track3Ds[NCAND];
std::map <int, double> Tracker::hitIndexDriftDLeftMap;
std::map <int, double> Tracker::hitIndexDriftDRightMap;
TGraphErrors * Tracker::graph_pairX = NULL;
TGraphErrors * Tracker::graph_pairZ = NULL;

Tracker::Tracker(InputOutputManager::InputHitType theInputHitType):
    nGoodTracks(0),
    func_pairYX(0),
    func_pairYZ(0),
    inputHitType(theInputHitType),
    fMaxResults(NCAND),
    ierflg(0),
    amin(0),
    edm(0),
    errdef(0),
    nvpar(0),
    nparx(0),
    icstat(0),
    slxStep(1),
    slxMin(0),
    slxMax(0),
    slzStep(1),
    slzMin(0),
    slzMax(0),
    inxStep(1),
    inxMin(0),
    inxMax(0),
    inzStep(1),
    inzMin(0),
    inzMax(0),
    t0Offset(0),
    t0OffsetRange(0),
    skipLayerAllowed(false),
    minChi2Input(0),
    sortType(NDFchi2)
{
    hitIndexInTestLayer = new std::vector<int>;
    hitLayerIndexMap = new std::vector<std::vector<int>*>;
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        hitLayerIndexMap->push_back(new std::vector<int>);
    }
    func_pairYX = new TF1("func_pairYX","pol1",GeometryManager::Get().fScintillator->Ydown,GeometryManager::Get().fScintillator->Yup); // x VS y
    func_pairYZ = new TF1("func_pairYZ","pol1",GeometryManager::Get().fScintillator->Ydown,GeometryManager::Get().fScintillator->Yup); // x VS y
    graph_pairZ = new TGraphErrors(NLAY);
    graph_pairX = new TGraphErrors(NLAY);

    Reset();
}

Tracker::~Tracker(){
    delete hitIndexInTestLayer;
    delete hitLayerIndexMap;
    delete func_pairYX;
    delete func_pairYZ;
    delete graph_pairX;
    delete graph_pairZ;
}

void Tracker::Reset(){
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        hitLayerIndexMap->at(iLayer)->clear();
    }
    pairs.clear();
    hitIndexInTestLayer->clear();
    hitIndexDriftDLeftMap.clear();
    hitIndexDriftDRightMap.clear();

    nGoodTracks = 0;
    currentTrack2D.Reset();
    currentTrack3D.Reset();
    func_pairYX->SetParameter(0,0);
    func_pairYX->SetParameter(1,0);
    func_pairYZ->SetParameter(0,0);
    func_pairYZ->SetParameter(1,0);
    slxMin = BeamManager::Get().beamSlx-BeamManager::Get().beamSlxRange;
    slxMax = BeamManager::Get().beamSlx+BeamManager::Get().beamSlxRange;
    slzMin = BeamManager::Get().beamSlz-BeamManager::Get().beamSlzRange;
    slzMax = BeamManager::Get().beamSlz+BeamManager::Get().beamSlzRange;
    inxMin = BeamManager::Get().beamInx-BeamManager::Get().beamInxRange;
    inxMax = BeamManager::Get().beamInx+BeamManager::Get().beamInxRange;
    inzMin = BeamManager::Get().beamInz-BeamManager::Get().beamInzRange;
    inzMax = BeamManager::Get().beamInz+BeamManager::Get().beamInzRange;
    double stereoAngle = fabs(GeometryManager::Get().fChamber->wire_theta[4][4]); // TODO Later: this is just a random cell to get a stereo angle
    inxStep = 0.01; // 10 um as the torlerance
    inzStep = inxStep/tan(stereoAngle);
    slxStep = inxStep/GeometryManager::Get().fChamber->chamberHeight;
    slzStep = inzStep/GeometryManager::Get().fChamber->chamberHeight;
    MyNamedVerbose("Tracking","Reset fitting parameters (stereo angle is taken from cell [4,4] as "<<stereoAngle<<" rad");
    MyNamedVerbose("Tracking","  slopeX:     "<<slxStep<<", "<<slxMin<<" ~ "<<slxMax);
    MyNamedVerbose("Tracking","  slopeZ:     "<<slzStep<<", "<<slzMin<<" ~ "<<slzMax);
    MyNamedVerbose("Tracking","  interceptX: "<<inxStep*1000<<" um, "<<inxMin<<" ~ "<<inxMax<<" mm");
    MyNamedVerbose("Tracking","  interceptZ: "<<inzStep*1000<<" um, "<<inzMin<<" ~ "<<inzMax<<" mm");
}

void Tracker::SetSortType(SortType t){
    sortType = t;
}

bool Tracker::SetMaxResults(int n){
    if (n>NCAND){
        MyError("Maximum number of results required "<<n<<" is greater than its capacity "<<NCAND<<". Will ignore this setting and keep using "<<fMaxResults);
        return false;
    }
    fMaxResults = n;
    return true;
}

bool Tracker::GoodForTracking(){
    int nPairsMin = ParameterManager::Get().TrackingParameters.nPairsMin;
    int nHitsGMax = ParameterManager::Get().TrackingParameters.nHitsGMax;
    /// Loop in layers, see how many pairs we can get at most
    /// TODO Later: In this way we are picking one hit per layer to form pairs. Support to multiple tracks/turns/in-out-segments to be implemented.
    int nPairs = 0;
    int nHitsG = 0;
    int prev_lid = -1;
    for (int lid = 0; lid < NLAY; lid++){
        if (hitLayerIndexMap->at(lid)->size()==0) continue;
        nHitsG += hitLayerIndexMap->at(lid)->size();
        if (prev_lid>=0&&(lid+prev_lid)%2==1){ // a cross between even and odd layer
            nPairs++;
        }
        prev_lid = lid;
    }
    if (nPairs<nPairsMin||nHitsG>nHitsGMax){
        return false;
    }
    return true;
}

void Tracker::Print(TString opt){
    if (opt.Contains("h")){ // hits
        for (int lid = 0; lid<NLAY; lid++){
            printf(" => Layer %d: %d good hits\n",lid,(int)hitLayerIndexMap->at(lid)->size());
            for (size_t iter = 0; iter<hitLayerIndexMap->at(lid)->size(); iter++){
                int iHit = hitLayerIndexMap->at(lid)->at(iter);
                int wid = InputOutputManager::Get().CellID->at(iHit);
                printf("    %d [%d,%d]\n",iHit,lid,wid);
            }
        }
    }
    if (opt.Contains("p")){ // picked hits
        for (size_t iter = 0; iter<currentTrack3D.hitIndexSelected.size(); iter++){
            int iHit = currentTrack3D.hitIndexSelected[iter];
            int lid = InputOutputManager::Get().LayerID->at(iHit);
            int wid = InputOutputManager::Get().CellID->at(iHit);
            printf(" => Pick %d: %d [%d,%d]\n",(int)iter,iHit,lid,wid);
        }
    }
    if (opt.Contains("t")){ // tracking results
        printf("There is %d tracks reconstructed!\n",nGoodTracks);
        for (int i = 0; i<nGoodTracks; i++){
            printf(" => %d: [%d,%d] nHitsS %d chi2 %.1e chi2a %.1e chi2-wt %.1e slx %.2e slz %.2e inx %.2e inz %.2e\n",
                    i,track2Ds[i].iSelection,track2Ds[i].iCombination,
                    (int)track3Ds[i].hitIndexSelected.size(),track3Ds[i].chi2,track3Ds[i].chi2a,track3Ds[i].chi2WithTestLayer,
                    track3Ds[i].slopeX,track3Ds[i].slopeZ,track3Ds[i].interceptX,track3Ds[i].interceptZ);
        }
    }
}

void Tracker::DoTracking(){
    size_t nSelections = 0;
    for (int it0 = -t0OffsetRange; it0<=t0OffsetRange; it0++){
        t0Offset = it0;
        MyNamedInfo("Tracking","Setting t0 offset to "<<t0Offset<<" ns");
        updateDriftD();
        tracking(0,nSelections); // 0 means starting from layer 0; nSelections is the number of possible choices by selecting one hit per layer (will increment in the recursive function)
    }
}

void Tracker::updateDriftD(){
    MyNamedVerbose("Tracking","  Assigning driftD to good hits");
    for (size_t lid = 0; lid<hitLayerIndexMap->size(); lid++){
        size_t nHits = hitLayerIndexMap->at(lid)->size();
        for (size_t i = 0; i<nHits; i++){
            int hitIndex = hitLayerIndexMap->at(lid)->at(i);
            double driftT = InputOutputManager::Get().DriftT->at(hitIndex)+t0Offset;
            int wid = InputOutputManager::Get().CellID->at(hitIndex);
            int status;
            hitIndexDriftDLeftMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,-1,status);
            hitIndexDriftDRightMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,1,status);
            MyNamedVerbose("Tracking","   ["<<lid<<","<<wid<<"] "<<driftT<<" ns, L "<<hitIndexDriftDLeftMap[hitIndex]<<" mm, R "<<hitIndexDriftDRightMap[hitIndex]<<" mm");
        }
    }
    size_t nHits = hitIndexInTestLayer->size();
    MyNamedVerbose("Tracking","  Assigning driftD to "<<hitIndexInTestLayer->size()<<" hits in test layer");
    for (size_t i = 0; i<nHits; i++){
        int hitIndex = hitIndexInTestLayer->at(i);
        double driftT = InputOutputManager::Get().DriftT->at(hitIndex)+t0Offset;
        int lid = InputOutputManager::Get().LayerID->at(hitIndex);
        int wid = InputOutputManager::Get().CellID->at(hitIndex);
        int status;
        hitIndexDriftDLeftMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,-1,status);
        hitIndexDriftDRightMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,1,status);
        MyNamedVerbose("Tracking","   ["<<lid<<","<<wid<<"] "<<driftT<<" ns, L "<<hitIndexDriftDLeftMap[hitIndex]<<" mm, R "<<hitIndexDriftDRightMap[hitIndex]<<" mm");
    }
}

int Tracker::tracking(int iLayer,size_t & iselection){
    if (iLayer == NLAY){ // finished picking hits
        int nHitsSel = currentTrack3D.hitIndexSelected.size();
        MyNamedVerbose("Tracking",Form(" Finished picking selection %d: nHitsSel = %d",(int)iselection,nHitsSel));
        int NDF = nHitsSel-currentTrack3D.nPars();
        if (NDF>=1){
            fitting(iselection);
            iselection++;
        }
    }
    else{
        for (int i = -1; i<(int)hitLayerIndexMap->at(iLayer)->size(); i++){ // -1 means don't pick up any hit in this layer 
            if (i==-1){
                if (hitLayerIndexMap->at(iLayer)->size()==0||(skipLayerAllowed&&NLAY-iLayer-1+(int)currentTrack3D.hitIndexSelected.size()-currentTrack3D.nPars()>=1)){
                    MyNamedVerbose("Tracking",Form(" => skip layer %d",iLayer));
                    tracking(iLayer+1,iselection);
                }
            }
            else{
                int ihit = hitLayerIndexMap->at(iLayer)->at(i);
                int wid = InputOutputManager::Get().CellID->at(ihit);
                MyNamedVerbose("Tracking",Form(" => pick %d: %d [%d,%d]",(int)(currentTrack3D.hitIndexSelected.size()),ihit,iLayer,wid));
                currentTrack3D.hitIndexSelected.push_back(ihit);
                currentTrack3D.hitLeftRightSelected.push_back(0);
                tracking(iLayer+1,iselection);
                currentTrack3D.hitIndexSelected.pop_back();
                currentTrack3D.hitLeftRightSelected.pop_back();
            }
        }
    }
    return 0;
}

int Tracker::fitting(int iselection){
    int nPairsMin = ParameterManager::Get().TrackingParameters.nPairsMin;
    /// calculate number of left/right combinations and start a loop in them and tracking them individually
    int ncombi = pow(2,currentTrack3D.hitIndexSelected.size());
    if (Log::GetLogLevel("Tracking")>=Log::VerboseLevel) Print("p");
    MyNamedVerbose("Tracking",Form("  %d picked hits -> %d L/R combinations",(int)(currentTrack3D.hitIndexSelected.size()),ncombi));
    for (int icombi = 0; icombi<ncombi; icombi++){ // each combination corresponds to a unique left/right selection set
        currentTrack3D.iSelection = iselection;
        currentTrack3D.iCombination = icombi;
        currentTrack3D.t0Offset = t0Offset;
        bool fittingSucceeded = false;
        double chi2X = 0;
        double chi2Z = 0;

        /// 1. set left right for picked hits and reset 2-D functions
        setLeftRight(icombi);
        reset2DFunctions();

        /// 2.A Fit pairs on Y-X and Y-Z plains without error
        formPairs(); // using the default 2-D functions
        MyNamedVerbose("Tracking","     combi "<<icombi<<" nPairs (default) = "<<pairs.size());
        if (pairs.size()<nPairsMin){
            continue;
        }
        fittingSucceeded = fit2D(chi2X,chi2Z); // fit without error
        if (!fittingSucceeded) continue;
        double iinx = func_pairYX->Eval(GeometryManager::Get().ReferenceY);
        double iinz = func_pairYZ->Eval(GeometryManager::Get().ReferenceY);
        double islx = func_pairYX->GetParameter(1);
        double islz = func_pairYZ->GetParameter(1);
        double chi2i, chi2pi, chi2ai, chi2wi;
        getchi2(chi2i,chi2pi,chi2ai,chi2wi,islx,iinx,islz,iinz,true);
        MyNamedVerbose("Tracking",Form("       2D FITTING RESULT: x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e; chi2 %.3e chi2a %.3e chi2w",islx,GeometryManager::Get().ReferenceY,iinx,chi2X,islz,GeometryManager::Get().ReferenceY,iinz,chi2Z,chi2i,chi2ai,chi2wi));
        // if the quality is too bad, skip 3-D fitting to save time
        if (minChi2Input&&chi2i>minChi2Input) continue;


        /// 4. if 2-D fitting is successful, do 3-D fitting with TMinuit using 2-D fitting result as input
        doFitting(islx,iinx,islz,iinz);
        // check the result
        double temp;
        double slx,inx,slz,inz;
        gMinuit->GetParameter(0, slx, temp);
        gMinuit->GetParameter(1, inx, temp);
        gMinuit->GetParameter(2, slz, temp);
        gMinuit->GetParameter(3, inz, temp);
        double chi2, chi2p, chi2a, chi2w;
        getchi2(chi2,chi2p,chi2a,chi2w,slx,inx,slz,inz,true);
        int nHitsSel = currentTrack3D.hitIndexSelected.size();
        MyNamedVerbose("Tracking",Form("       TMinuit fitting RESULT: nHitsSel = %d, x=%.3e*(y-%.3e)+%.3e, z=%.3e*(y-%.3e)+%.3e, chi2i = %.3e chi2 = %.3e",nHitsSel,slx,GeometryManager::Get().ReferenceY,inx,slz,GeometryManager::Get().ReferenceY,inz,chi2i,chi2));

        /// 5. At last, if the newly fitted track meets our requirments, store the fitting result.
        currentTrack2D = (TrackCandidate)currentTrack3D;
        currentTrack2D.nPairs = pairs.size();
        currentTrack2D.interceptX = iinx;
        currentTrack2D.interceptZ = iinz;
        currentTrack2D.slopeX = islx;
        currentTrack2D.slopeZ = islz;
        currentTrack2D.chi2X = chi2X;
        currentTrack2D.chi2Z = chi2Z;
        currentTrack2D.chi2 = chi2i;
        currentTrack2D.chi2a = chi2ai;
        currentTrack2D.chi2WithTestLayer = chi2wi;
        currentTrack2D.pValue = chi2pi;
        currentTrack3D.slopeX = slx;
        currentTrack3D.slopeZ = slz;
        currentTrack3D.interceptX = inx;
        currentTrack3D.interceptZ = inz;
        currentTrack3D.chi2 = chi2;
        currentTrack3D.pValue = chi2p;
        currentTrack3D.chi2a = chi2a;
        currentTrack3D.chi2WithTestLayer = chi2w;
        if (fMaxResults){ // there is a limit on number of fitting results to save. Sort by chi2 and NDF.
            checkAndFitIn();
        }
        else{ // there is no limit (except for its capacity NCAND) we don't have to sort.
            if (nGoodTracks<NCAND){
                track3Ds[nGoodTracks] = currentTrack3D;
                nGoodTracks++;
            }
        }

        // TODO LATER: get mc input
        if (inputHitType == InputOutputManager::kMCDriftD || inputHitType == InputOutputManager::kMCDriftT){
            //getchi2(chi2mc,chi2pmc,chi2amc,chi2wmc,i_slxmc,i_inxmc,i_slzmc,i_inzmc,true);
        }
    }
    return 0;
}

void Tracker::setLeftRight(int icombi){
    for (size_t iter = 0; iter<currentTrack3D.hitIndexSelected.size(); iter++){
        int ilr = (icombi&(1<<iter))>>iter;
        if (ilr==0) ilr = -1; // 1 for right, -1 for left
        currentTrack3D.hitLeftRightSelected[iter] = ilr;
    }
}

void Tracker::reset2DFunctions(double MoveRatioX, double MoveRatioZ){
    func_pairYX->SetParameters(BeamManager::Get().beamInx-BeamManager::Get().beamSlx*GeometryManager::Get().ReferenceY+BeamManager::Get().beamInxRange*MoveRatioX,
                               BeamManager::Get().beamSlx);
    func_pairYZ->SetParameters(BeamManager::Get().beamInz-BeamManager::Get().beamSlz*GeometryManager::Get().ReferenceY+BeamManager::Get().beamInzRange*MoveRatioZ,
                               BeamManager::Get().beamSlz);
}

bool Tracker::fit2D(double & chi2X, double & chi2Z ){
    setPairPositionGraphs();

    // prepare to do the 1-D fitting
    // Do the fitting on Z pairs
    if (gMinuit) delete gMinuit;
    gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
    gMinuit->SetFCN(fcnZ);
    gMinuit->SetPrintLevel(-1); // no print
    arglist[0] = 0;
    gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    double temp; double inz = func_pairYZ->GetParameter(0); double slz = func_pairYZ->GetParameter(1);
    gMinuit->mnparm(0, "inz", inz, inzStep, inzMin,inzMax,ierflg);
    gMinuit->mnparm(1, "slz", slz, slzStep, slzMin,slzMax,ierflg);
    arglist[0] = 5000.0;
    arglist[1] = 1.0;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    gMinuit->GetParameter(0, inz, temp);
    gMinuit->GetParameter(1, slz, temp);
    func_pairYZ->SetParameter(0,inz); func_pairYZ->SetParameter(1,slz);

    // Do the fitting on X pairs
    if (gMinuit) delete gMinuit;
    gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
    gMinuit->SetFCN(fcnX);
    gMinuit->SetPrintLevel(-1); // no print
    arglist[0] = 0;
    gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    double inx = func_pairYX->GetParameter(0); double slx = func_pairYX->GetParameter(1);
    gMinuit->mnparm(0, "inx", inx, inxStep, inxMin,inxMax,ierflg);
    gMinuit->mnparm(1, "slx", slx, slxStep, slxMin,slxMax,ierflg);
    arglist[0] = 5000.0;
    arglist[1] = 1.0;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    gMinuit->GetParameter(0, inx, temp);
    gMinuit->GetParameter(1, slx, temp);
    func_pairYX->SetParameter(0,inx); func_pairYX->SetParameter(1,slx);

    double iinx = func_pairYX->Eval(GeometryManager::Get().ReferenceY);
    double iinz = func_pairYZ->Eval(GeometryManager::Get().ReferenceY);
    double islx = func_pairYX->GetParameter(1);
    double islz = func_pairYZ->GetParameter(1);
    getChi2XZ(chi2X,chi2Z);
    TString debugContent = Form("%.3e %.3e %.3e %.3e",islx,iinx,islz,iinz);
    for (int ipair = 0; ipair<pairs.size(); ipair++){
        debugContent += Form(" %.3e %.3e %.3e %.3e",pairs[ipair].pairX,func_pairYX->Eval(pairs[ipair].pairY),pairs[ipair].pairZ,func_pairYZ->Eval(pairs[ipair].pairY));
    }
    MyNamedDebug("check",debugContent);
    return true;
}

///////////////////////////////////////////////////////////////////////////////////////
void Tracker::formPairs(){
    pairs.clear();
    // NOTE: this is assuming hitIndexSelected is sorted in order of layer ID
    for (int iter = 0; iter<(int)currentTrack3D.hitIndexSelected.size()-1; iter++){
        int ihit = currentTrack3D.hitIndexSelected[iter];
        int jhit = currentTrack3D.hitIndexSelected[iter+1];
        int iLR = currentTrack3D.hitLeftRightSelected[iter];
        int jLR = currentTrack3D.hitLeftRightSelected[iter+1];
        Pair aPair(ihit,jhit,iLR,jLR);
        if(updatePairPosition(aPair)){
            pairs.push_back(aPair);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////
bool Tracker::updatePairPosition(Pair & aPair){
    int ihit = aPair.hitIndexL;
    int jhit = aPair.hitIndexH;
    int lid = InputOutputManager::Get().LayerID->at(ihit);
    int ljd = InputOutputManager::Get().LayerID->at(jhit);
    int wid = InputOutputManager::Get().CellID->at(ihit);
    int wjd = InputOutputManager::Get().CellID->at(jhit);
    double theta1 = GeometryManager::Get().fChamber->wire_thetaX[lid][wid];
    double theta2 = GeometryManager::Get().fChamber->wire_thetaX[ljd][wjd];
    if (theta1==theta2){
        return false; // cannot make a cross;
    }
    double y1 = getWireY(lid,wid);
    double y2 = getWireY(ljd,wjd);
    double deltaY = y2-y1;
    double dd1 = aPair.LRL>=0?hitIndexDriftDRightMap[ihit]:hitIndexDriftDLeftMap[ihit];
    double dd2 = aPair.LRH>=0?hitIndexDriftDRightMap[jhit]:hitIndexDriftDLeftMap[jhit];
    double sintheta12 = sin(theta1-theta2);
    double zc_fix_slx = deltaY*func_pairYX->GetParameter(1)/(tan(theta2)-tan(theta1));
    double xc = GeometryManager::Get().fChamber->wirecross_x[lid][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
    double yc = (y1+y2)/2;
    double zc = GeometryManager::Get().fChamber->wirecross_z[lid][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
    MyNamedVerbose("makePair",Form("        xc = %.3e+%.3e*sin(%.3e)/(-sin(%.3e-%.3e))+%.3e*sin(%.3e)/sin(%.3e-%.3e)",GeometryManager::Get().fChamber->wirecross_x[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2));
    MyNamedVerbose("makePair",Form("        zc = %.3e+%.3e*cos(%.3e)/(-sin(%.3e-%.3e))+%.3e*cos(%.3e)/sin(%.3e-%.3e)+%.3e",GeometryManager::Get().fChamber->wirecross_z[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2,zc_fix_slx));
    MyNamedVerbose("makePair",Form("       cp[%d,%d]: w(%d,%d) i(%d,%d) dd(%f,%f) xyz(%f,%f,%f)",lid,ljd,wid,wjd,ihit,jhit,dd1,dd2,xc,yc,zc));
    double chamberHL = GeometryManager::Get().fChamber->chamberLength/2;
    if (zc<-chamberHL||zc>chamberHL){
        MyNamedVerbose("makePair",Form("       BAD @ layer %d & layer %d!",lid,ljd));
        return false;
    }
    aPair.pairX = xc;
    aPair.pairY = yc;
    aPair.pairZ = zc;
    return true;
}

///////////////////////////////////////////////////////////////////////////////////////
double Tracker::getWireY(int lid, int wid){
    double wyro = GeometryManager::Get().fChamber->wire_y[lid][wid][1];
    double wzro = GeometryManager::Get().fChamber->wire_z[lid][wid][1];
    double wyhv = GeometryManager::Get().fChamber->wire_y[lid][wid][0];
    double wzhv = GeometryManager::Get().fChamber->wire_z[lid][wid][0];
    double wy0 = (wyro+wyhv)/2.;// assume wy
    double wz = func_pairYZ->Eval(wy0);// get wz by extrapolating the track to wy
    // correct wy according to wz
    double wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
    wz = func_pairYZ->Eval(wy);
    wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
    MyNamedVerbose("getWireY","    wire["<<lid<<"]["<<wid<<"]: "<<wy0<<" -> "<<wy);
    return wy;
}

///////////////////////////////////////////////////////////////////////////////////////
int Tracker::setPairPositionGraphs(){
    size_t nPairs = pairs.size();
    graph_pairX->Set(nPairs);
    graph_pairZ->Set(nPairs);
    for (int i = 0; i<nPairs; i++){
        graph_pairX->SetPoint(i,pairs[i].pairY,pairs[i].pairX);
        graph_pairZ->SetPoint(i,pairs[i].pairY,pairs[i].pairZ);
    }
    return 0;
}

void Tracker::getChi2XZ(double & chi2X, double & chi2Z){
    // calculate pairXyz
    chi2X = 0;
    chi2Z = 0;
    for (int ipair = 0; ipair<pairs.size(); ipair++){
        double tchi2Z = pairs[ipair].pairZ-func_pairYZ->Eval(pairs[ipair].pairY);
        double tchi2X = pairs[ipair].pairX-func_pairYX->Eval(pairs[ipair].pairY);
        MyNamedVerbose("getChi2XZ",Form("           getChi2XZ pair[%d]: error x = %.3e, error z = %.3e",ipair,tchi2X,tchi2Z));
        chi2Z += pow(tchi2Z,2);
        chi2X += pow(tchi2X,2);
    }
    if (pairs.size()){
        chi2Z/=pairs.size();
        chi2X/=pairs.size();
    }
}

void Tracker::doFitting(double sliX, double iniX,double sliZ, double iniZ){
    if (gMinuit) delete gMinuit;
    // TODO LATER: consider about different track representation
    gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
    gMinuit->SetFCN(fcn);
    gMinuit->SetPrintLevel(-1); // no print

    arglist[0] = 0;
    gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    // Set starting values and step sizes for parameters
    gMinuit->mnparm(0, "slopeX", sliX, slxStep, slxMin,slxMax,ierflg);
    gMinuit->mnparm(1, "interceptX", iniX, inxStep, inxMin,inxMax,ierflg);
    gMinuit->mnparm(2, "slopeZ", sliZ, slzStep, slzMin,slzMax,ierflg);
    gMinuit->mnparm(3, "interceptZ", iniZ, inzStep, inzMin,inzMax,ierflg);

    // Now ready for minimization step
    arglist[0] = 5000.0;
    arglist[1] = 1.0;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    // Print results
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    //printf("====Rrestul====\n");
    //gMinuit->mnprin(3,amin);
}

void Tracker::fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
    // TODO LATER: consider about different track representation
    double cp,ca,cw;
    if (npar<4){
        MyError("Not enough par! npar = "<<npar<<", but we need 4");
    }
    else{
        getchi2(f,cp,ca,cw,*par,*(par+1),*(par+2),*(par+3),false);
    }
}

void Tracker::getchi2(double &f, double & cp, double & ca, double & cw, double slx, double inx, double slz, double inz, bool all)
{
    //calculate chisquare
    double chisq = 0;
    double delta;
    double dfit;

    int N = -currentTrack3D.nPars();
    bool LayerUsed[NLAY] = {false};
    for (size_t i=0;i<currentTrack3D.hitIndexSelected.size(); i++) {
        int iHit = currentTrack3D.hitIndexSelected[i];
        int lr = currentTrack3D.hitLeftRightSelected[i];
        int lid = InputOutputManager::Get().LayerID->at(iHit);
        int wid = InputOutputManager::Get().CellID->at(iHit);
        dfit = GeometryManager::Get().GetDOCA(lid,wid,slx,inx,slz,inz);
        double dd = lr>=0?hitIndexDriftDRightMap[iHit]:hitIndexDriftDLeftMap[iHit];
        double error = XTManager::Get().GetError(fabs(dd));
        delta  = (dfit-dd)/error;
        chisq += delta*delta;
        LayerUsed[lid] = true;
        N++;
    }
    f = chisq;
    if (N>0) f/=N;

    if (all){ // should calculate chi2_pValue (cp) and chi2_all (ca)
        cp = TMath::Prob(f*N,N);
        // now find the closest hit in the test layer and add it into ca
        int nLayerAdded = 0;
        bool found = false;
        double minres = 1e9;
        double theDD = 0;
        for (int lid = 0; lid<NLAY; lid++){
            if (LayerUsed[lid]) continue; // already used this layer in chi2
            found = false;
            minres = 1e9;
            theDD = 0;
            for (size_t i=0;i<hitLayerIndexMap->at(lid)->size(); i++) {
                int iHit = hitLayerIndexMap->at(lid)->at(i);
                int wid = InputOutputManager::Get().CellID->at(iHit);
                dfit = GeometryManager::Get().GetDOCA(lid,wid,slx,inx,slz,inz);
                double dd = dfit>=0?hitIndexDriftDRightMap[iHit]:hitIndexDriftDLeftMap[iHit];
                if (fabs(minres)>fabs(dfit-dd)){
                    minres = dfit-dd;
                    theDD = dd;
                    found = true;
                }
            }
            if (found){
                double error = XTManager::Get().GetError(fabs(theDD));
                chisq += minres*minres/error/error;
                nLayerAdded++;
            }
        }
        ca = chisq;
        if (N+nLayerAdded>0) ca/=(N+nLayerAdded);
        // then calculate with test layer hit
        found = false;
        minres = 1e9;
        theDD = 0;
        for (size_t i=0;i<hitIndexInTestLayer->size(); i++) {
            int iHit = hitIndexInTestLayer->at(i);
            int lid = InputOutputManager::Get().LayerID->at(iHit);
            int wid = InputOutputManager::Get().CellID->at(iHit);
            dfit = GeometryManager::Get().GetDOCA(lid,wid,slx,inx,slz,inz);
            double dd = dfit>=0?hitIndexDriftDRightMap[iHit]:hitIndexDriftDLeftMap[iHit];
            if (fabs(minres)>fabs(dfit-dd)){
                minres = dfit-dd;
                theDD = dd;
                found = true;
            }
        }
        if (found){
            double error = XTManager::Get().GetError(fabs(theDD));
            chisq += minres*minres/error/error;
            cw = chisq/(N+nLayerAdded+1);
        }
        else{
            cw = ca;
        }
    }
}

void Tracker::fcnX(int &npar, double *gin, double &f, double *par, int iflag){
    f = getchi2Graph(graph_pairX,*par,*(par+1));
}

void Tracker::fcnZ(int &npar, double *gin, double &f, double *par, int iflag){
    f = getchi2Graph(graph_pairZ,*par,*(par+1));
}

double Tracker::getchi2Graph(TGraphErrors* graph, double v0, double sl){
    double chi2 = 0;
    for (size_t i = 0; i<graph->GetN(); i++){
        double x,y; double yerr;
        graph->GetPoint(i,x,y);
        yerr = graph->GetErrorY(i);
        double pred = v0+sl*x;
        chi2+=pow((pred-y)/(yerr?yerr:1),2);
    }
    return chi2;
}

bool Tracker::checkAndFitIn(){
    MyNamedVerbose("Tracking"," checking new result with "<<currentTrack3D.hitIndexSelected.size()<<" hits and chi2 = "<<currentTrack3D.chi2<<" chi2a = "<<currentTrack3D.chi2a<<" chi2w = "<<currentTrack3D.chi2WithTestLayer);
    int insertAt = -1;
    int takeOut = -1;
    // scan through current results rank list
    // mark insertAt as the rank of the new result (keep -1 if not good enough or the list is empty)
    // mark takeOut if the new result if identical to one on the list and the new result is better (keep -1 if not)
    for (int i = 0; i<nGoodTracks; i++){
        if (currentTrack3D.hitIndexSelected.size()<track3Ds[i].hitIndexSelected.size()) continue;
        // WARNING: now we rely on total chi2 including test layer hit, a slight bias
        if ((sortType==chi2a&&currentTrack3D.chi2a<track3Ds[i].chi2a)
          ||(sortType==chi2WithTestLayer&&currentTrack3D.chi2WithTestLayer<track3Ds[i].chi2WithTestLayer)
          ||(sortType==NDFchi2
             &&(currentTrack3D.hitIndexSelected.size()>track3Ds[i].hitIndexSelected.size()
                ||(currentTrack3D.hitIndexSelected.size()==track3Ds[i].hitIndexSelected.size()&&currentTrack3D.chi2<track3Ds[i].chi2)
               )
            )
           ){
            if (insertAt<0){ // modify the index to be insert at if it's not set yet. Keep on searching in case we may find a bad candidate later with the same hits.
                MyNamedVerbose("Tracking"," better than Cand#"<<i<<" with "<<track3Ds[i].hitIndexSelected.size()<<" hits and chi2 = "<<track3Ds[i].chi2<<" chi2a = "<<track3Ds[i].chi2a<<" chi2w = "<<track3Ds[i].chi2WithTestLayer);
                insertAt = i;
            }
            if (currentTrack3D == track3Ds[i]){ // they have used the same hits (with same left/right)
                MyNamedVerbose("Tracking"," same with Cand#"<<i);
                takeOut = i;
                break; // no need to continue
            }
        }
    }
    if (insertAt<0){ // no one on the list is worse than the new result
        if (nGoodTracks<fMaxResults){ // the list is not full yet
            MyNamedVerbose("Tracking"," ==> Add to the end, nGoodTracks "<<nGoodTracks<<" -> "<<(nGoodTracks+1));
            insertAt = nGoodTracks; // put the current result at the bottom
            takeOut = -1; // kick out nobody
            nGoodTracks++; // now we have a new track result
        }
        else{ // the list is full
            MyNamedVerbose("Tracking"," ==> Skip");
            return false;
        }
    }
    else if (takeOut<0){ // no one to be replaced by the new result
        if(nGoodTracks==fMaxResults){ // the list is full, then kick out the last one
            takeOut = nGoodTracks-1;
            MyNamedVerbose("Tracking"," ==> Insert at "<<insertAt<<", take out the last one at "<<takeOut);
        }
        else{ // the list is not full yet, make a room for it
            takeOut = nGoodTracks;
            MyNamedVerbose("Tracking"," ==> Insert at "<<insertAt<<", nGoodTracks "<<nGoodTracks<<" -> "<<(nGoodTracks+1));
            nGoodTracks++;
        }
    }
    else{
        MyNamedVerbose("Tracking"," ==> Insert at "<<insertAt<<", take out "<<takeOut);
    }
    for (int i = takeOut; i>insertAt; i--){ // move the candidates back by 1 and kick out the one to be replaced
        track3Ds[i] = track3Ds[i-1];
        track2Ds[i] = track2Ds[i-1];
    }
    track3Ds[insertAt] = currentTrack3D; // put the new result in position
    track2Ds[insertAt] = currentTrack2D; // put the new result in position
    return true;
}
