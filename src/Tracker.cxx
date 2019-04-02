#include <TGraphErrors.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>

#include "GeometryManager.hxx"
#include "Tracker.hxx"
#include "InputOutputManager.hxx"
#include "BeamManager.hxx"
#include "XTManager.hxx"
#include "Log.hxx"

std::vector<int>               * Tracker::hitIndexInTestLayer = NULL;
std::vector<std::vector<int>*> * Tracker::hitLayerIndexMap = NULL;
TrackResult Tracker::currentTrackResult;
TrackResult Tracker::trackResults[NCAND];
std::map <int, int>    Tracker::hitIndexLeftRight;
std::map <int, double> Tracker::hitIndexDriftDLeftMap;
std::map <int, double> Tracker::hitIndexDriftDRightMap;

Tracker::Tracker(InputOutputManager::InputType theInputType):
    pairableLayers(0),
    nPairs(0),
    nGoodPairs(0),
    func_pairYX(0),
    func_pairYZ(0),
    graph_pairX(0),
    graph_pairZ(0),
    inputType(theInputType),
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
    inzMax(0)
{
    hitIndexInTestLayer = new std::vector<int>;
    hitLayerIndexMap = new std::vector<std::vector<int>*>;
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        hitLayerIndexMap->push_back(new std::vector<int>);
    }
    pairableLayers = new std::vector<int>;
    func_pairYX = new TF1("func_pairYX","pol1",GeometryManager::Get().fScintillator->Ydown,GeometryManager::Get().fScintillator->Yup); // x VS y
    func_pairYZ = new TF1("func_pairYZ","pol1",GeometryManager::Get().fScintillator->Ydown,GeometryManager::Get().fScintillator->Yup); // x VS y
    graph_pairZ = new TGraphErrors(NLAY);
    graph_pairX = new TGraphErrors(NLAY);

    gMinuit = new TMinuit(5);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(fcn);

	arglist[0] = 0;
	gMinuit->SetPrintLevel(-1); // no print
	gMinuit->mnexcm("SET NOW", arglist ,1,ierflg); // no warning 

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    Reset();
}

Tracker::~Tracker(){
    delete hitIndexInTestLayer;
    delete hitLayerIndexMap;
    delete pairableLayers;
    delete func_pairYX;
    delete func_pairYZ;
    delete graph_pairX;
    delete graph_pairZ;
}

void Tracker::Reset(){
    for (int iLayer = 0; iLayer<NLAY; iLayer++){
        hitLayerIndexMap->at(iLayer)->clear();
        pickIndex[iLayer] = -1;
        pickLR[iLayer] = 0;
        pickWireY[iLayer] = 0;
        pairX[iLayer] = 0;
        pairY[iLayer] = 0;
        pairZ[iLayer] = 0;
    }
    hitIndexInTestLayer->clear();
    pairableLayers->clear();
    hitIndexLeftRight.clear();
    hitIndexDriftDLeftMap.clear();
    hitIndexDriftDRightMap.clear();

    currentTrackResult.Reset();
    nPairs = 0;
    nGoodPairs = 0;
    for (int iCand = 0; iCand<NCAND; iCand++){
        trackResults[iCand].Reset();
    }
    func_pairYX->SetParameter(0,0);
    func_pairYX->SetParameter(1,0);
    func_pairYZ->SetParameter(0,0);
    func_pairYZ->SetParameter(1,0);
    slxStep = BeamManager::Get().beamSlxRange/1e4; // FIXME: should leave the step as an option
    slxMin = BeamManager::Get().beamSlx-BeamManager::Get().beamSlxRange;
    slxMax = BeamManager::Get().beamSlx+BeamManager::Get().beamSlxRange;
    slzStep = BeamManager::Get().beamSlzRange/1e4;
    slzMin = BeamManager::Get().beamSlz-BeamManager::Get().beamSlzRange;
    slzMax = BeamManager::Get().beamSlz+BeamManager::Get().beamSlzRange;
    inxStep = BeamManager::Get().beamInxRange/1e4;
    inxMin = BeamManager::Get().beamInx-BeamManager::Get().beamInxRange;
    inxMax = BeamManager::Get().beamInx+BeamManager::Get().beamInxRange;
    inzStep = BeamManager::Get().beamInzRange/1e4;
    inzMin = BeamManager::Get().beamInz-BeamManager::Get().beamInzRange;
    inzMax = BeamManager::Get().beamInz+BeamManager::Get().beamInzRange;
}

void Tracker::DoTracking(){
    int nSelections = 0;
    updateDriftD();
    tracking(0,nSelections); // 0 means starting from the 1st pick; nSelections is the number of possible choices by selecting one hit per layer;
}

void Tracker::updateDriftD(){
    for (int lid = 0; lid<hitLayerIndexMap->size(); lid++){
        int nHits = hitLayerIndexMap->at(lid)->size();
        for (int i = 0; i<nHits; i++){
            int hitIndex = hitLayerIndexMap->at(lid)->at(i);
            double driftT = InputOutputManager::Get().DriftT->at(hitIndex); // FIXME: should add t0 offset consideration here
            int wid = InputOutputManager::Get().CellID->at(hitIndex);
            int status;
            hitIndexLeftRight[hitIndex] = 0;
            hitIndexDriftDLeftMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,-1,status);
            hitIndexDriftDRightMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,1,status);
        }
    }
    int nHits = hitIndexInTestLayer->size();
    for (int i = 0; i<nHits; i++){
        int hitIndex = hitIndexInTestLayer->at(i);
        double driftT = InputOutputManager::Get().DriftT->at(hitIndex); // FIXME: should add t0 offset consideration here
        int lid = InputOutputManager::Get().LayerID->at(hitIndex);
        int wid = InputOutputManager::Get().CellID->at(hitIndex);
        int status;
        hitIndexLeftRight[hitIndex] = 0;
        hitIndexDriftDLeftMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,-1,status);
        hitIndexDriftDRightMap[hitIndex] = XTManager::Get().t2x(driftT,lid,wid,1,status);
    }
}

int Tracker::tracking(int ipick,int & iselection){
    if (ipick == pairableLayers->size()){ // finished picking hits
        MyNamedVerbose("Tracking",Form(" Finished picking selection %d:",iselection));
        fitting(iselection);
        iselection++;
    }
    else{
        int lid = pairableLayers->at(ipick);
        for (int i = 0; i<hitLayerIndexMap->at(lid)->size(); i++){
            int ihit = hitLayerIndexMap->at(lid)->at(i);
            int wid = InputOutputManager::Get().CellID->at(ihit);
            MyNamedVerbose("Tracking",Form(" => pick # %d, layer %d, wire %d, hit[%d], ihit = %d",ipick,lid,wid,i,ihit));
            pickIndex[ipick] = hitLayerIndexMap->at(lid)->at(i);
            pickLR[ipick] = 0;
            pickWireY[ipick] = 0;
            pairX[ipick] = 0;
            pairY[ipick] = 0;
            pairZ[ipick] = 0;
            tracking(ipick+1,iselection);
        }
    }
    return 0;
}

int Tracker::fitting(int iselection){
    /// calculate number of left/right combinations and start a loop in them and tracking them individually
    int nPicks = pairableLayers->size();
    int ncombi = pow(2,nPicks);
    MyNamedVerbose("Tracking",Form("  %d picked layers -> %d combinations",nPicks,ncombi));
    for (int icombi = 0; icombi<ncombi; icombi++){ // each combination corresponds to a unique left/right selection set
        MyNamedVerbose("Tracking",Form("     combi %d",icombi));
        nPairs = 0;
        nGoodPairs = 0;
        bool inScint = false;
        bool fromSource = false;
        bool fittingSucceeded = false;
        double chi2X = 0;
        double chi2Z = 0;

        /// 1. Prepare drift distance according to left/right choices
        InputOutputManager::Get().nCandidatesFound++;
        setLeftRight(icombi); // set left right for picked hits

        /// 2. Get pair positions according to the initial track parameter. Fit pairs on Y-X and Y-Z plains
        Reset2DFunctions();
        fittingSucceeded = Fit2D(2.5,false,chi2X,chi2Z,inScint,fromSource); // fit without error
        if (!fittingSucceeded) continue;
        ///    If fitting result is not good, reset the initial values of these 2-D functions with 1/2 offset (on Z direction)
        ///    Do the 2-D fitting again based on the previous fitting resutl
        if (!fromSource||!inScint) Reset2DFunctions(0,0.5);
        fittingSucceeded = Fit2D(2.5,true,chi2X,chi2Z,inScint,fromSource); // fit with error
        if (!fittingSucceeded) continue;
        ///    If fitting result is not good, reset the initial values of these 2-D functions with -1/2 offset (on Z direction)
        ///    Do the 2-D fitting again based on the previous fitting resutl
        if (!fromSource||!inScint) Reset2DFunctions(0,-0.5);
        fittingSucceeded = Fit2D(2.5,true,chi2X,chi2Z,inScint,fromSource); // fit with error
        if (!fittingSucceeded) continue;
        
        /// 3. If the 2-D fitting is successfull, proceed to 3-D fitting
        if (inScint&&fromSource&&nGoodPairs>=3){ // good candidate
            ///    First, pick up closest hits (with 2 mm cut) based on the given track parameters from 2-D fitting
            pickUpHitsForFitting(2); // FIXME: should tune the error limit
            int nHitsSel = currentTrackResult.hitIndexSelected.size();
            MyNamedVerbose("Tracking",Form("       Selection finished: nHitsSel = %d",nHitsSel));
            ///    If the number of hits picked is greater than 5, then go to doFitting with current track parameters as initial fitting values
            if (nHitsSel>=5){ // FIXME: should consider change the cut
                // set the initial track candidate from 2-D fitting
                double iinx = func_pairYX->Eval(GeometryManager::Get().ReferenceY);
                double iinz = func_pairYZ->Eval(GeometryManager::Get().ReferenceY);
                double islx = func_pairYX->GetParameter(1);
                double islz = func_pairYZ->GetParameter(1);
                currentTrackResult.initialTrackCandidate.interceptX = iinx;
                currentTrackResult.initialTrackCandidate.interceptZ = iinz;
                currentTrackResult.initialTrackCandidate.slopeX = islx;
                currentTrackResult.initialTrackCandidate.slopeZ = islz;
                currentTrackResult.initialTrackCandidate.chi2X = chi2X;
                currentTrackResult.initialTrackCandidate.chi2Z = chi2Z;
                currentTrackResult.initialTrackCandidate.hitIndexSelected = currentTrackResult.hitIndexSelected;
                currentTrackResult.initialTrackCandidate.hitLeftRightSelected = currentTrackResult.hitLeftRightSelected;
                double chi2i, chi2pi, chi2ai;
                getchi2(chi2i,chi2pi,chi2ai,islx,iinx,islz,iinz,0,true);
                currentTrackResult.initialTrackCandidate.chi2 = chi2i;
                currentTrackResult.initialTrackCandidate.chi2WithTestLayer = chi2ai;
                currentTrackResult.initialTrackCandidate.pValue = chi2pi;
                // fitting with TMinuit
                doFitting(islx,iinx,islz,iinz);
                double temp;
                double slx,inx,slz,inz;
                gMinuit->GetParameter(0, slx, temp);
                gMinuit->GetParameter(1, inx, temp);
                gMinuit->GetParameter(2, slz, temp);
                gMinuit->GetParameter(3, inz, temp);
                // update fitD
                // reselect
                double chi2, chi2p, chi2a;
                getchi2(chi2,chi2p,chi2a,slx,inx,slz,inz,0,true);
                MyNamedVerbose("Tracking",Form("       1st fitting RESULT: nHitsSel = %d, x=%.3e*(y-%.3e)+%.3e, z=%.3e*(y-%.3e)+%.3e, chi2i = %.3e chi2 = %.3e",nHitsSel,slx,GeometryManager::Get().ReferenceY,inx,slz,GeometryManager::Get().ReferenceY,inz,chi2i,chi2));
                ///    Second, pick up hits again (with 1 mm cut) based on the updated track parameters
                pickUpHitsForFitting(1); // FIXME: should tune the error limit
                nHitsSel = currentTrackResult.hitIndexSelected.size();
                ///    If the number of hits picked is greater than 5, then go to doFitting with new track parameters
                if (nHitsSel>=5){ // FIXME: should consider change the cut
                    // fitting with TMinuit
                    doFitting(islx,iinx,islz,iinz);
                    double temp;
                    gMinuit->GetParameter(0, slx, temp);
                    gMinuit->GetParameter(1, inx, temp);
                    gMinuit->GetParameter(2, slz, temp);
                    gMinuit->GetParameter(3, inz, temp);
                    ///    At last, if the newly fitted track meets our requirments, store the fitting result.
                    inScint = GeometryManager::Get().IsInScinti(1.5,inx,slx,inz,slz); // FIXME: error limit should be tuned
                    fromSource = BeamManager::Get().IsInBeam(slx,slz);
                    if (inScint&&fromSource){
                        // update chi2
                        if (inputType == InputOutputManager::kMCDriftD || inputType == InputOutputManager::kMCDriftT){
                            // FIXME: get mc input
							//getchi2(chi2mc,chi2pmc,chi2amc,i_slxmc,i_inxmc,i_slzmc,i_inzmc,0,true);
                        }
                        getchi2(chi2,chi2p,chi2a,slx,inx,slz,inz,0,true);
                        // check chi2 and see where the result fits
                        MyNamedVerbose("Tracking",Form("       2nd fitting RESULT: nHitsSel = %d, x=%.3e*(y-%.3e)+%.3e, z=%.3e*(y-%.3e)+%.3e, chi2i = %.3e chi2 = %.3e",nHitsSel,slx,GeometryManager::Get().ReferenceY,inx,slz,GeometryManager::Get().ReferenceY,inz,chi2i,chi2));
                        currentTrackResult.slopeX = slx;
                        currentTrackResult.slopeZ = slz;
                        currentTrackResult.interceptX = inx;
                        currentTrackResult.interceptZ = inz;
                        currentTrackResult.chi2 = chi2;
                        currentTrackResult.pValue = chi2p;
                        currentTrackResult.chi2WithTestLayer = chi2a;
                        checkResults(nHitsSel,icombi,iselection);
                    }
                }
            }
        }
    }
    return 0;
}

void Tracker::setLeftRight(int icombi){
    int nPicks = pairableLayers->size();
    for (int iPick = 0; iPick<nPicks; iPick++){
        int ilr = (icombi&(1<<iPick))>>iPick;
        if (ilr==0) ilr = -1; // 1 for right, -1 for left
        pickLR[iPick] = ilr;
        hitIndexLeftRight[pickIndex[iPick]] = ilr;
    }
}

void Tracker::Reset2DFunctions(double MoveRatioZ, double MoveRatioX){
    func_pairYX->SetParameters(BeamManager::Get().beamInx+BeamManager::Get().beamInxRange*MoveRatioX,
                               BeamManager::Get().beamSlx);
    func_pairYZ->SetParameters(BeamManager::Get().beamInz+BeamManager::Get().beamInzRange*MoveRatioZ,
                               BeamManager::Get().beamSlz);
}

bool Tracker::Fit2D(double safetyFactor, bool fitWithError, double & chi2X, double & chi2Z, bool & inScint, bool & fromSource){
    int nPicks = pairableLayers->size();
    updateWirePositionOnHit(); // fix wy positions
    if (updatePairPositions()) return false; // cannot make pairs
    setPairPositionGraphs(fitWithError);
    graph_pairZ->Fit("func_pairYZ","QN0G","");
    graph_pairX->Fit("func_pairYX","QN0G","");
    double iinx = func_pairYX->Eval(GeometryManager::Get().ReferenceY);
    double iinz = func_pairYZ->Eval(GeometryManager::Get().ReferenceY);
    double islx = func_pairYX->GetParameter(1);
    double islz = func_pairYZ->GetParameter(1);
    inScint = GeometryManager::Get().IsInScinti(safetyFactor,iinx,islx,iinz,islz); // FIXME: need to tune
    fromSource = BeamManager::Get().IsInBeam(islx,islz);
    nGoodPairs = getChi2XZ(chi2X,chi2Z);
    MyNamedVerbose("Tracking",Form("       2D FITTING RESULT: nGoodPairs = %d, inScint? %s, fromSource? %s; x=%.3e*(y-%.3e)+%.3e, chi2 = %.3e, z=%.3e*(y-%.3e)+%.3e, chi2 = %.3e",nGoodPairs,inScint?"yes":"no",fromSource?"yes":"no",islx,GeometryManager::Get().ReferenceY,iinx,chi2X,islz,GeometryManager::Get().ReferenceY,iinz,chi2Z));
    TString debugContent = Form("%.3e %.3e %.3e %.3e",islx,iinx,islz,iinz);
    for (int ipair = 0; ipair<nPairs; ipair++){
        debugContent += Form(" %.3e %.3e %.3e %.3e",pairX[ipair],func_pairYX->Eval(pairY[ipair]),pairZ[ipair],func_pairYZ->Eval(pairY[ipair]));
    }
    MyNamedDebug("check",debugContent);
    return true;
}

int Tracker::updateWirePositionOnHit(){
    int nPicks = pairableLayers->size();
    // calculate pickWireY
    for (int ipick = 0; ipick<nPicks; ipick++){
        // Get hit information
        int ihit = pickIndex[ipick];
        int lid = InputOutputManager::Get().LayerID->at(ihit);
        int wid = InputOutputManager::Get().CellID->at(ihit);
        double wyro = GeometryManager::Get().fChamber->wire_y[lid][wid][1];
        double wzro = GeometryManager::Get().fChamber->wire_z[lid][wid][1];
        double wyhv = GeometryManager::Get().fChamber->wire_y[lid][wid][0];
        double wzhv = GeometryManager::Get().fChamber->wire_z[lid][wid][0];
        double wy = (wyro+wyhv)/2.;// assume wy
        double wz = func_pairYZ->Eval(wy);// get wz by extrapolating the track to wy
        // correct wy according to wz
        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
        wz = func_pairYZ->Eval(wy);
        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
        MyNamedVerbose("updateWirePositionOnHit","    pickWireY["<<ipick<<"]: "<<pickWireY[ipick]<<" -> "<<wy);
        pickWireY[ipick] = wy;
    }
    return 0;
}

int Tracker::updatePairPositions(){
    double chamberHL = GeometryManager::Get().fChamber->chamberLength/2;
    int nPicks = pairableLayers->size();
    // calculate pairXyz
    nPairs = 0;
    int ipick = 0;
    for (; ipick<nPicks-1; ipick++){
        int ihit = pickIndex[ipick];
        int jhit = pickIndex[ipick+1];
        int lid = InputOutputManager::Get().LayerID->at(ihit);
        int ljd = InputOutputManager::Get().LayerID->at(jhit);
        if (lid+1!=ljd) continue; // not adjacent
        double deltaY = pickWireY[ipick+1]-pickWireY[ipick];
        double dd1 = pickLR[ipick]>=0?hitIndexDriftDRightMap[ihit]:hitIndexDriftDLeftMap[ihit];
        double dd2 = pickLR[ipick+1]>=0?hitIndexDriftDRightMap[jhit]:hitIndexDriftDLeftMap[jhit];
        int wid = InputOutputManager::Get().CellID->at(ihit);
        int wjd = InputOutputManager::Get().CellID->at(jhit);
        double theta1 = GeometryManager::Get().fChamber->wire_theta[lid][wid];
        double theta2 = GeometryManager::Get().fChamber->wire_theta[ljd][wjd];
        double sintheta12 = sin(theta1-theta2);
        double zc_fix_slx = deltaY*func_pairYX->GetParameter(1)/(tan(theta2)-tan(theta1));
        double xc = GeometryManager::Get().fChamber->wirecross_x[lid][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
        double zc = GeometryManager::Get().fChamber->wirecross_z[lid][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
        pairX[nPairs] = xc;
        pairY[nPairs] = (pickWireY[ipick+1]+pickWireY[ipick])/2.;
        pairZ[nPairs] = zc;
        MyNamedVerbose("updatePairPositions",Form("        xc = %.3e+%.3e*sin(%.3e)/(-sin(%.3e-%.3e))+%.3e*sin(%.3e)/sin(%.3e-%.3e)",GeometryManager::Get().fChamber->wirecross_x[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2));
        MyNamedVerbose("updatePairPositions",Form("        zc = %.3e+%.3e*cos(%.3e)/(-sin(%.3e-%.3e))+%.3e*cos(%.3e)/sin(%.3e-%.3e)+%.3e",GeometryManager::Get().fChamber->wirecross_z[lid][wid][wjd],dd1,theta2,theta1,theta2,dd2,theta1,theta1,theta2,zc_fix_slx,zc_fix_slx));
        MyNamedVerbose("updatePairPositions",Form("       cp[%d,%d]: w(%d,%d) i(%d,%d) dd(%f,%f) xyz(%f,%f,%f)",lid,ljd,wid,wjd,ihit,jhit,dd1,dd2,xc,(pickWireY[ipick+1]+pickWireY[ipick])/2.,zc));
        if (zc<-chamberHL||zc>chamberHL){
            MyNamedVerbose("updatePairPositions",Form("       BAD combination!"));
            break;
        }
        nPairs++;
    }
    if (ipick==nPicks-1){
        MyNamedVerbose("updatePairPositions",Form("       GOOD combination!"));
        return 0;
    }
    else{
        MyNamedVerbose("updatePairPositions",Form("       BAD @ %d!",ipick));
        return 1;
    }
}

int Tracker::setPairPositionGraphs(bool withError){
    graph_pairX->Set(nPairs);
    graph_pairZ->Set(nPairs);
    for (int i = 0; i<nPairs; i++){
        graph_pairX->SetPoint(i,pairY[i],pairX[i]);
        graph_pairZ->SetPoint(i,pairY[i],pairZ[i]);
    }
    double errorzMax0 = 0;
    double errorzMax1 = 0;
    int errorzMax0_i = -1;
    int errorzMax1_i = -1;
    for (int ipair = 0; ipair<nPairs; ipair++){
        double errorz = 0;
        double errorx = 0;
        if (withError){
            errorz = fabs(func_pairYZ->Eval(pairY[ipair])-pairZ[ipair]);
            errorx = fabs(func_pairYX->Eval(pairY[ipair])-pairX[ipair]);
        }
        if (errorzMax0<errorz){
            errorzMax0 = errorz;
            errorzMax0_i = ipair;
        }
        else if (errorzMax1<errorz){
            errorzMax1 = errorz;
            errorzMax1_i = ipair;
        }
	    graph_pairZ->SetPointError(ipair,0,0.1);
	    graph_pairX->SetPointError(ipair,0,0.1);
    }
    if (errorzMax0>4) {
        graph_pairZ->SetPointError(errorzMax0_i,0,10);
        graph_pairX->SetPointError(errorzMax0_i,0,10);
        MyNamedVerbose("Tracking",Form("           setPairPositionGraphs pair[%d]: error x = ->10, error z = %.3e->10",errorzMax0_i,errorzMax0));
    }
    if (errorzMax1>4) {
        graph_pairZ->SetPointError(errorzMax1_i,0,10);
        graph_pairX->SetPointError(errorzMax1_i,0,10);
        MyNamedVerbose("Tracking",Form("           setPairPositionGraphs pair[%d]: error x = ->10, error z = %.3e->10",errorzMax1_i,errorzMax1));
    }
    return 0;
}

int Tracker::getChi2XZ(double & chi2X, double & chi2Z){
    // calculate pairXyz
    chi2X = 0;
    chi2Z = 0;
    int nCount = 0;
    for (int ipair = 0; ipair<nPairs; ipair++){
        double tchi2Z = pairZ[ipair]-func_pairYZ->Eval(pairY[ipair]);
        double tchi2X = pairX[ipair]-func_pairYX->Eval(pairY[ipair]);
        MyNamedVerbose("Tracking",Form("           getChi2XZ pair[%d]: error x = %.3e, error z = %.3e",ipair,tchi2X,tchi2Z));
        if (fabs(tchi2Z)<4&&fabs(tchi2X)<1){ // FIXME: error limit should be tuned
            chi2Z += pow(tchi2Z,2);
            chi2X += pow(tchi2X,2);
            nCount++;
        }
    }
    if (nCount){
        chi2Z/=nCount;
        chi2X/=nCount;
    }
    return nCount;
}

void Tracker::pickUpHitsForFitting(double residualCut){
    // get track parameters from the fitted 2-D functions
    double iinx = func_pairYX->Eval(GeometryManager::Get().ReferenceY);
    double iinz = func_pairYZ->Eval(GeometryManager::Get().ReferenceY);
    double islx = func_pairYX->GetParameter(1);
    double islz = func_pairYZ->GetParameter(1);
    currentTrackResult.hitIndexSelected.clear();
    // loop in all hits and pick up hits close to the track
    for (int lid = 0; lid<hitLayerIndexMap->size(); lid++){ // note that the hits here are already selected in main function, and hits from the test layer are not included in this map
        double resMin = 1e9;
        int    theHitIndex = -1;
        int    theLR = 0;
        for (int i = 0; i<hitLayerIndexMap->at(lid)->size(); i++){
            int iHit = hitLayerIndexMap->at(lid)->at(i);
            int wid = InputOutputManager::Get().CellID->at(iHit);
            double doca = GeometryManager::Get().GetDOCA(lid,wid,islx,iinx,islz,iinz);
            double driftD = doca>=0?hitIndexDriftDRightMap[iHit]:hitIndexDriftDLeftMap[iHit];
            double residual = doca-driftD;
            if (fabs(residual)<fabs(resMin)){
                resMin = residual;
                theHitIndex = iHit;
                if (hitIndexLeftRight[iHit])
                    theLR = hitIndexLeftRight[iHit];
                else
                    theLR = doca>=0?1:-1;
            }
        }
        if (fabs(resMin)<residualCut){
            // FIXME: URGENT, should consider about the l/r from pick list!
            currentTrackResult.hitIndexSelected.push_back(theHitIndex);
            currentTrackResult.hitLeftRightSelected.push_back(theLR);
        }
    }
}

void Tracker::doFitting(double sliX, double iniX,double sliZ, double iniZ){
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
	double cp,ca;
    getchi2(f,cp,ca,*par,*(par+1),*(par+2),*(par+3),0,false);
}

void Tracker::getchi2(double &f, double & cp, double & ca, double slx, double inx, double slz, double inz, double t0offset,bool all)
{
	//calculate chisquare
	double chisq = 0;
	double delta;
	double dfit;

    int N = -4;
	for (int i=0;i<currentTrackResult.hitIndexSelected.size(); i++) {
	    int iHit = currentTrackResult.hitIndexSelected[i];
	    int lr = currentTrackResult.hitLeftRightSelected[i];
	    int lid = InputOutputManager::Get().LayerID->at(iHit);
	    int wid = InputOutputManager::Get().CellID->at(iHit);
		dfit = GeometryManager::Get().GetDOCA(lid,wid,slx,inx,slz,inz);
		double dd = lr>=0?hitIndexDriftDRightMap[iHit]:hitIndexDriftDLeftMap[iHit];
		double error = XTManager::Get().GetError(fabs(dd));
        delta  = (dfit-dd)/error;
		chisq += delta*delta;
		N++;
	}
	if (N>0) chisq/=N;
	f = chisq;
	if (all){ // should calculate chi2_pValue (cp) and chi2_all (ca)
		cp = TMath::Prob(f*N,N);
		// now find the closest hit in the test layer and add it into ca
		double minres = 1e9;
		bool found = false;
		for (int i=0;i<hitIndexInTestLayer->size(); i++) {
            int iHit = hitIndexInTestLayer->at(i);
            int lid = InputOutputManager::Get().LayerID->at(iHit);
            int wid = InputOutputManager::Get().CellID->at(iHit);
			dfit = GeometryManager::Get().GetDOCA(lid,wid,slx,inx,slz,inz);
            double dd = dfit>=0?hitIndexDriftDRightMap[iHit]:hitIndexDriftDLeftMap[iHit];
			if (fabs(minres)>fabs(dfit-dd)){
				minres = dfit-dd;
				found = true;
			}
		}
		if (found){
			double error = 0.2;
			ca = f*N+minres*minres/error/error;
			ca/=(N+1);
		}
		else{
			ca = f;
		}
	}
}

bool Tracker::checkResults(int nHitsSel, int icombi, int iselection){
	bool covered = false;
	int insertAt = -1;
	int takeOut = -1;
    for (int i = 0; i<NCAND; i++){
        if (currentTrackResult.NDF<trackResults[i].NDF) continue;
        if (currentTrackResult == trackResults[i]){ // they have used the same hits (with same left/right)
            MyNamedVerbose("Tracking"," same with Cand#"<<i);
			// FIXME: WARNING, now we rely on total chi2 including test layer hit, a slight bias
			// TODO: make this an option
    		if (currentTrackResult.chi2WithTestLayer<trackResults[i].chi2WithTestLayer){
                MyNamedVerbose("Tracking","   better than Cand#"<<i);
                takeOut = i;
			}
			break;
    	}
    	if (insertAt<0){ // modify the index to be insert at if it's not set yet. Keep on searching in case we may find a bad candidate later with the same hits.
            if (currentTrackResult.NDF>trackResults[i].NDF
                    ||currentTrackResult.chi2WithTestLayer<trackResults[i].chi2WithTestLayer){
                MyNamedVerbose("Tracking"," better than Cand#"<<i);
                insertAt = i;
            }
        }
        takeOut = i; // by default if no candidate with same selection was found, the last one will be considered to be moved out (in case insertAt is valid)
        if (trackResults[i].pValue==0){ // FIXME: not a very reliable check here. Currently pValue = 0 means this result is empty
            break; // stop searching after an empty result is reached
        }
	}
	if (insertAt<0) return false;
	for (int i = takeOut; i>insertAt; i--){ // move the candidates back by 1 and kick out the one to be replaced
	    trackResults[i] = trackResults[i-1];
	}
	trackResults[insertAt] = currentTrackResult; // put the new result in position
	return true;
}
