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
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.CandSelBy")) XTAnalyzerParameters.CandSelBy = MyRuntimeParameters::Get().GetParameterS("XTAnalyzer.CandSelBy");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.XTType")) XTAnalyzerParameters.XTType = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.XTType");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.AsymXT")) XTAnalyzerParameters.AsymXT = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.AsymXT");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.RequireInTriggerCounter")) XTAnalyzerParameters.RequireInTriggerCounter = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.RequireInTriggerCounter");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.RequireAllGoldenHits")) XTAnalyzerParameters.RequireAllGoldenHits = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.RequireAllGoldenHits");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.ClosestPeak")) XTAnalyzerParameters.ClosestPeak = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.ClosestPeak");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.UseGoodHit")) XTAnalyzerParameters.UseGoodHit = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.UseGoodHit");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.AllGoodHitsUsed")) XTAnalyzerParameters.AllGoodHitsUsed = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.AllGoodHitsUsed");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.nHits_max")) XTAnalyzerParameters.nHits_max = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.nHits_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.nHitsS_min")) XTAnalyzerParameters.nHitsS_min = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.nHitsS_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.chi2_max")) XTAnalyzerParameters.chi2_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.chi2_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.pValue_min")) XTAnalyzerParameters.pValue_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.pValue_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.slz_min")) XTAnalyzerParameters.slz_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.slz_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.slz_max")) XTAnalyzerParameters.slz_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.slz_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.gold_t_min")) XTAnalyzerParameters.gold_t_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.gold_t_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.gold_t_max")) XTAnalyzerParameters.gold_t_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.gold_t_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_n_min")) XTAnalyzerParameters.bin_n_min = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_n_min");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_min")) XTAnalyzerParameters.bin_t_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_min");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_max")) XTAnalyzerParameters.bin_t_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_max");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_num")) XTAnalyzerParameters.bin_t_num = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_t_num");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_min")) XTAnalyzerParameters.bin_x_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_min");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_max")) XTAnalyzerParameters.bin_x_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_max");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_num")) XTAnalyzerParameters.bin_x_num = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_x_num");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_landTmin")) XTAnalyzerParameters.bin_t_landTmin = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_landTmin");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_landTmax")) XTAnalyzerParameters.bin_t_landTmax = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_landTmax");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_landXmin")) XTAnalyzerParameters.bin_x_landXmin = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_landXmin");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_landXmax")) XTAnalyzerParameters.bin_x_landXmax = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_landXmax");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_fit_num")) XTAnalyzerParameters.bin_x_fit_num = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_x_fit_num");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_fit_num")) XTAnalyzerParameters.bin_t_fit_num = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_t_fit_num");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_fit_num_tail")) XTAnalyzerParameters.bin_t_fit_num_tail = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_t_fit_num_tail");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_tailTime")) XTAnalyzerParameters.bin_t_tailTime = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_tailTime");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_ratio")) XTAnalyzerParameters.bin_t_ratio = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_ratio");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_ratio")) XTAnalyzerParameters.bin_x_ratio = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_ratio");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_n_min")) XTAnalyzerParameters.graph_n_min = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.graph_n_min");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_chi2_max")) XTAnalyzerParameters.graph_chi2_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.graph_chi2_max");;
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_sepX")) XTAnalyzerParameters.graph_sepX = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.graph_sepX");;
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
    XTAnalyzerParameters.Print();
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

XTAnalyzerPara::XTAnalyzerPara(){
    XTType = 55; // use pol5 and pol5 to describe an XT function
    AsymXT = false; // don't use asymmetric XT
    CandSelBy = "Original"; // Just use the first candidate
    RequireInTriggerCounter = true;
    RequireAllGoldenHits = false;
    ClosestPeak = false;
    UseGoodHit = false;
    AllGoodHitsUsed = false;
    nHits_max = 30;
    nHitsS_min = 7;
    chi2_max = 0;
    pValue_min = 0.05;
    slz_min = -0.1;
    slz_max = 0.1;
    gold_t_min = 50;
    gold_t_max = 300;
    bin_n_min = 0; // minimum number of entries in one slice to apply fitting function; Otherwise use mean value & RMS instead.
    bin_t_min = -25-1/0.96/2; // t range for one x bin
    bin_t_max = 800+1/0.96/2;
    bin_t_num = 792+1;
    bin_t_fit_num = 3; //  number of T bins to project
    bin_t_tailTime = 350; // starting tail mode from here: different number of bins as a T slice to fit X
    bin_t_fit_num_tail = 3; //  number of T bins to project
    bin_x_min = -0.02;
    bin_x_max = 10.02; // x range for one t bin
    bin_x_num = 501;
    bin_t_landTmin = 330;
    bin_t_landTmax = 450;
    bin_x_landXmin = 7;
    bin_x_landXmax = 8.5;
    bin_x_fit_num = 3; //  number of X bins to project
    bin_t_ratio = 0.5;
    bin_x_ratio = 0.5;
    graph_n_min = 50;
    graph_chi2_max = 100;
    graph_sepX = 1;
}

void XTAnalyzerPara::Print(){
    printf("XTAnalyzer Parameters:\n");
    printf(" About event selection:\n");
    printf("  Select candidate by = %s\n",CandSelBy.Data());
    printf("  RequireInTriggerCounter = %s\n",RequireInTriggerCounter?"true":"false");
    printf("  RequireAllGoldenHits = %s\n",RequireAllGoldenHits?"true":"false");
    printf("  UseGoodHit = %s\n",UseGoodHit?"true":"false");
    printf("  AllGoodHitsUsed = %s\n",AllGoodHitsUsed?"true":"false");
    printf("  nHits_max = %d\n",nHits_max);
    printf("  nHitsS_min = %d\n",nHitsS_min);
    printf("  chi2_max = %f\n",chi2_max);
    printf("  pValue_min = %f\n",pValue_min);
    printf("  slz_min = %f\n",slz_min);
    printf("  slz_max = %f\n",slz_max);
    printf(" About hit selection:\n");
    printf("  ClosestPeak = %s\n",ClosestPeak?"true":"false");
    printf("  golden hit t_min = %f\n",gold_t_min);
    printf("  golden hit t_max = %f\n",gold_t_max);
    printf(" About bin-by-bin analysis :\n");
    printf("  X axis: number of bins = %d\n",bin_x_num);
    printf("  X axis: minimum = %f\n",bin_x_min);
    printf("  X axis: maximum = %f\n",bin_x_max);
    printf("  T axis: number of bins = %d\n",bin_t_num);
    printf("  T axis: minimum = %f\n",bin_t_min);
    printf("  T axis: maximum = %f\n",bin_t_max);
    printf("  X fitting: Landau function range lower edge = %f\n",bin_t_landTmin);
    printf("  X fitting: Landau function range lower edge = %f\n",bin_t_landTmax);
    printf("  X fitting: tail starting from t = %.1f ns\n",bin_t_tailTime);
    printf("  X fitting: number of T bins to project before tail: %d\n",bin_t_fit_num);
    printf("  X fitting: number of T bins to project in tail : %d\n",bin_t_fit_num_tail);
    printf("  X fitting: fit in X range according to cut ratio: %.1f\n",bin_t_ratio);
    printf("  T fitting: number of X bins to project: %d\n",bin_x_fit_num);
    printf("  T fitting: fit in T range according to cut ratio: %.1f\n",bin_x_ratio);
    printf("  T fitting: Landau function range lower edge = %f\n",bin_x_landXmin);
    printf("  T fitting: Landau function range lower edge = %f\n",bin_x_landXmax);
    printf(" About XT graphs:\n");
    printf("  Minimum number of entries of the sample point to be included in graph: %d\n",graph_n_min);
    printf("  Maximum chi2 of the sample point to be included in graph: %.3e\n",graph_chi2_max);
    printf("  Separation X to combine space samples and time samples: %.3e\n",graph_sepX);
    printf(" About XT function:\n");
    printf("  XTType = %d\n",XTType);
    printf("  AsymXT = %s\n",AsymXT?"true":"false");
}

AnaPara::AnaPara(){
}

void AnaPara::Print(){
}
