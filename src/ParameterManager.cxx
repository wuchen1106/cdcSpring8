#include "ParameterManager.hxx"
#include "MyRuntimeParameters.hxx"

ParameterManager* ParameterManager::fParameterManager = NULL;

ParameterManager::ParameterManager()
{
    geoSetup = GeometryManager::kNormal;
    beamType = BeamManager::kSPring8;
    inputType = InputOutputManager::kData;
}

void ParameterManager::LoadParameters(ParaBlock theParaBlock){
    if (MyRuntimeParameters::Get().HasParameter("geoSetup")){
        TString name = MyRuntimeParameters::Get().GetParameterS("geoSetup");
        name.ToLower();
        if (name == "tilted") geoSetup = GeometryManager::kTilted;
        else if (name == "finger") geoSetup = GeometryManager::kFinger;
        else if (name == "normal") geoSetup = GeometryManager::kNormal;
    }
    if (MyRuntimeParameters::Get().HasParameter("beam")){
        TString name = MyRuntimeParameters::Get().GetParameterS("beam");
        name.ToLower();
        if (name == "spring8") beamType = BeamManager::kSPring8;
        else if (name == "spring8tilted") beamType = BeamManager::kSPring8;
        else if (name == "cosmic") beamType = BeamManager::kCosmic;
    }
    if (theParaBlock==kAll||theParaBlock==kTracking){
        if (MyRuntimeParameters::Get().HasParameter("tracking.inputType")){
            TString name = MyRuntimeParameters::Get().GetParameterS("tracking.inputType");
            name.ToLower();
            if (name == "data") inputType = InputOutputManager::kData;
            else if (name == "mc" || name == "mcdriftt") inputType = InputOutputManager::kMCDriftT;
            else if (name == "mcdriftd") inputType = InputOutputManager::kMCDriftD;
        }
        if (MyRuntimeParameters::Get().HasParameter("peakType")) TrackingParameters.peakType = MyRuntimeParameters::Get().GetParameterI("peakType");
        if (MyRuntimeParameters::Get().HasParameter("tracking.nHitsMax")) TrackingParameters.nHitsMax = MyRuntimeParameters::Get().GetParameterI("tracking.nHitsMax");
        if (MyRuntimeParameters::Get().HasParameter("tracking.nHitsSMin")) TrackingParameters.nHitsSMin = MyRuntimeParameters::Get().GetParameterI("tracking.nHitsSMin");
        if (MyRuntimeParameters::Get().HasParameter("tracking.nPairsMin")) TrackingParameters.nPairsMin = MyRuntimeParameters::Get().GetParameterI("tracking.nPairsMin");
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0shift0")) TrackingParameters.t0shift0 = MyRuntimeParameters::Get().GetParameterI("tracking.t0shift0");
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0shift1")) TrackingParameters.t0shift1 = MyRuntimeParameters::Get().GetParameterI("tracking.t0shift1");
        if (MyRuntimeParameters::Get().HasParameter("tracking.tmin")) TrackingParameters.tmin = MyRuntimeParameters::Get().GetParameterI("tracking.tmin");
        if (MyRuntimeParameters::Get().HasParameter("tracking.tmax")) TrackingParameters.tmax = MyRuntimeParameters::Get().GetParameterI("tracking.tmax");
        if (MyRuntimeParameters::Get().HasParameter("tracking.sumCut")) TrackingParameters.sumCut = MyRuntimeParameters::Get().GetParameterD("tracking.sumCut");
        if (MyRuntimeParameters::Get().HasParameter("tracking.aaCut")) TrackingParameters.aaCut = MyRuntimeParameters::Get().GetParameterD("tracking.aaCut");
        if (MyRuntimeParameters::Get().HasParameter("tracking.workType")) TrackingParameters.workType = MyRuntimeParameters::Get().GetParameterI("tracking.workType");;
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0error")) TrackingParameters.t0error = MyRuntimeParameters::Get().GetParameterD("tracking.t0error");;
        if (MyRuntimeParameters::Get().HasParameter("tracking.lidStart")) TrackingParameters.lidStart = MyRuntimeParameters::Get().GetParameterI("tracking.lidStart");;
        if (MyRuntimeParameters::Get().HasParameter("tracking.lidStop")) TrackingParameters.lidStop = MyRuntimeParameters::Get().GetParameterI("tracking.lidStop");;
    }
}

void ParameterManager::ReadInputFile(TString filename, TString dirName, bool tryFile, bool fixParameters){
    MyRuntimeParameters::Get().ReadInputFile(filename,dirName,tryFile,fixParameters);
}

void ParameterManager::Print(){
    printf("#########################ParameterManager###########################\n");
    printf("  geometry setup:    %d\n",geoSetup);
    printf("  beam type:         %d\n",beamType);
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
    inputType = 0; // 2 for MC without layer specified; 3 for MC with layer specified; 0 for data
    peakType = 0; // 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
    workType = 0; // fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
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
    printf("  workType    = %d, %s\n",workType,workType==0?"all as 0":(workType==1?"even/odd":(workType==-1?"even/odd reversed":"all layers")));
    printf("  inputType   = %d, %s\n",inputType,inputType==0?"Real Data":(inputType==2?"MC X":"MC T"));
    printf("  peakType    = %d, %s\n",peakType,peakType==0?"First peak over threshold":"All peaks over threshold");
}

CalibPara::CalibPara(){
}

void CalibPara::Print(){
}

XTAnalylzerPara::XTAnalylzerPara(){
}

void XTAnalylzerPara::Print(){
}

AnaPara::AnaPara(){
}

void AnaPara::Print(){
}
