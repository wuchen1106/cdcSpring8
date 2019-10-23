#include "ParameterManager.hxx"
#include "MyRuntimeParameters.hxx"

ParameterManager* ParameterManager::fParameterManager = NULL;

ParameterManager::ParameterManager()
{
    geoSetup = GeometryManager::kNormal;
    chamberType = GeometryManager::kProto4;
    connectionType = GeometryManager::kSPring8;
    beamType = BeamManager::kSPring8;
    inputHitType = InputOutputManager::kData;
    peakType = TrackingPara::kFirstPeak;
}

void ParameterManager::LoadParameters(ParaBlock theParaBlock){
    if (MyRuntimeParameters::Get().HasParameter("geoSetup")){
        TString name = MyRuntimeParameters::Get().GetParameterS("geoSetup");
        name.ToLower();
        if (name == "tilted") geoSetup = GeometryManager::kTilted;
        else if (name == "finger") geoSetup = GeometryManager::kFinger;
        else if (name == "normal") geoSetup = GeometryManager::kNormal;
    }
    if (MyRuntimeParameters::Get().HasParameter("chamberType")){
        TString name = MyRuntimeParameters::Get().GetParameterS("chamberType");
        name.ToLower();
        if (name == "cdc") chamberType = GeometryManager::kCDC;
        else if (name == "proto4") chamberType = GeometryManager::kProto4;
    }
    if (MyRuntimeParameters::Get().HasParameter("connectionType")){
        TString name = MyRuntimeParameters::Get().GetParameterS("connectionType");
        name.ToLower();
        if (name == "spring8") connectionType = GeometryManager::kSPring8;
        else if (name == "cosmic") connectionType = GeometryManager::kCosmic;
    }
    if (MyRuntimeParameters::Get().HasParameter("beam")){
        TString name = MyRuntimeParameters::Get().GetParameterS("beam");
        name.ToLower();
        if (name == "spring8") beamType = BeamManager::kSPring8;
        else if (name == "spring8tilted") beamType = BeamManager::kSPring8;
        else if (name == "cosmic") beamType = BeamManager::kCosmic;
    }
    if (MyRuntimeParameters::Get().HasParameter("inputHitType")){
        TString name = MyRuntimeParameters::Get().GetParameterS("inputHitType");
        name.ToLower();
        if (name == "data") inputHitType = InputOutputManager::kData;
        else if (name == "mc" || name == "mcdriftt") inputHitType = InputOutputManager::kMCDriftT;
        else if (name == "mcdriftd") inputHitType = InputOutputManager::kMCDriftD;
    }
    if (MyRuntimeParameters::Get().HasParameter("peakType")){
        TString name = MyRuntimeParameters::Get().GetParameterS("peakType");
        name.ToLower();
        if (name == "firstpeak") peakType = TrackingPara::kFirstPeak;
        else if (name == "highestpeak") peakType = TrackingPara::kHighestPeak;
        else if (name == "allpeaks") peakType = TrackingPara::kAllPeaks;
    }
    if (theParaBlock==kAll||theParaBlock==kTracking){
        if (MyRuntimeParameters::Get().HasParameter("tracking.nHitsMax")) TrackingParameters.nHitsMax = MyRuntimeParameters::Get().GetParameterI("tracking.nHitsMax");
        if (MyRuntimeParameters::Get().HasParameter("tracking.nHitsSMin")) TrackingParameters.nHitsSMin = MyRuntimeParameters::Get().GetParameterI("tracking.nHitsSMin");
        if (MyRuntimeParameters::Get().HasParameter("tracking.nPairsMin")) TrackingParameters.nPairsMin = MyRuntimeParameters::Get().GetParameterI("tracking.nPairsMin");
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0shift0")) TrackingParameters.t0shift0 = MyRuntimeParameters::Get().GetParameterI("tracking.t0shift0");
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0shift1")) TrackingParameters.t0shift1 = MyRuntimeParameters::Get().GetParameterI("tracking.t0shift1");
        if (MyRuntimeParameters::Get().HasParameter("tracking.tmin")) TrackingParameters.tmin = MyRuntimeParameters::Get().GetParameterI("tracking.tmin");
        if (MyRuntimeParameters::Get().HasParameter("tracking.tmax")) TrackingParameters.tmax = MyRuntimeParameters::Get().GetParameterI("tracking.tmax");
        if (MyRuntimeParameters::Get().HasParameter("tracking.sumCut")) TrackingParameters.sumCut = MyRuntimeParameters::Get().GetParameterD("tracking.sumCut");
        if (MyRuntimeParameters::Get().HasParameter("tracking.aaCut")) TrackingParameters.aaCut = MyRuntimeParameters::Get().GetParameterD("tracking.aaCut");
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0error")) TrackingParameters.t0error = MyRuntimeParameters::Get().GetParameterD("tracking.t0error");;
        if (MyRuntimeParameters::Get().HasParameter("tracking.lidStart")) TrackingParameters.lidStart = MyRuntimeParameters::Get().GetParameterI("tracking.lidStart");;
        if (MyRuntimeParameters::Get().HasParameter("tracking.lidStop")) TrackingParameters.lidStop = MyRuntimeParameters::Get().GetParameterI("tracking.lidStop");;
    }
    if (theParaBlock==kAll||theParaBlock==kXTAnalyzer){
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.XTType")) XTAnalylzerParameters.XTType = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.XTType");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.AsymXT")) XTAnalylzerParameters.AsymXT = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.AsymXT");
    }
}

void ParameterManager::ReadInputFile(TString filename, TString dirName, bool tryFile, bool fixParameters){
    MyRuntimeParameters::Get().ReadInputFile(filename,dirName,tryFile,fixParameters);
}

void ParameterManager::Print(){
    printf("#########################ParameterManager###########################\n");
    printf("  Chamber type:      %d\n",chamberType);
    printf("  connection type:   %d\n",connectionType);
    printf("  geometry setup:    %d\n",geoSetup);
    printf("  beam type:         %d\n",beamType);
    printf("  input type:        %d\n",inputHitType);
    printf("  peak type:         %d\n",peakType);
    TrackingParameters.Print();
    CalibParameters.Print();
    XTAnalylzerParameters.Print();
    AnaParameters.Print();
}

TrackingPara::TrackingPara(){
    nHitsMax = 30;
    nHitsSMin = 5;
    nPairsMin = 3;
    t0shift0 = 0;
    t0shift1 = 0;
    tmin = -20;
    tmax = 800;
    sumCut = -10;
    aaCut = 0;
    BlindLayer = 0; // Don't use this layer for tracking
    t0error = 0; // if t0 error is 0, then don't set it as a free parameter in fitting
    inislx = 0; // initial guess for slope x;
    iniinz = 100;
    lidStart = 1;
    lidStop = NLAY-1;
}

void TrackingPara::Print(){
    printf("Tracking Parameters:\n");
    printf("  nHitsMax    = %d\n",nHitsMax);
    printf("  nHitsSMin   = %d\n",nHitsSMin);
    printf("  nPairsMin   = %d\n",nPairsMin);
    printf("  t0shift b0  = %d\n",t0shift0);
    printf("  t0shift b1  = %d\n",t0shift1);
    printf("  tmin        = %d\n",tmin);
    printf("  tmax        = %d\n",tmax);
    printf("  sumCut      = %f\n",sumCut);
    printf("  aaCut       = %f\n",aaCut);
    printf("  BlindLayer  = %d\n",BlindLayer);
}

CalibPara::CalibPara(){
}

void CalibPara::Print(){
}

XTAnalylzerPara::XTAnalylzerPara(){
    XTType = 55; // use pol5 and pol5 to describe an XT function
    AsymXT = false; // don't use asymmetric XT
}

void XTAnalylzerPara::Print(){
    printf("XTAnalyzer Parameters:\n");
    printf("  XTType = %d\n",XTType);
    printf("  AsymXT = %d\n",AsymXT);
}

AnaPara::AnaPara(){
}

void AnaPara::Print(){
}
