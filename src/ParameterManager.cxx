#include "ParameterManager.hxx"
#include "XTAnalyzer.hxx"
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
        if (MyRuntimeParameters::Get().HasParameter("tracking.t0error")) TrackingParameters.t0error = MyRuntimeParameters::Get().GetParameterD("tracking.t0error");
        if (MyRuntimeParameters::Get().HasParameter("tracking.lidStart")) TrackingParameters.lidStart = MyRuntimeParameters::Get().GetParameterI("tracking.lidStart");
        if (MyRuntimeParameters::Get().HasParameter("tracking.lidStop")) TrackingParameters.lidStop = MyRuntimeParameters::Get().GetParameterI("tracking.lidStop");
    }
    if (theParaBlock==kAll||theParaBlock==kXTAnalyzer){
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.CandSelBy")) XTAnalyzerParameters.CandSelBy = MyRuntimeParameters::Get().GetParameterS("XTAnalyzer.CandSelBy");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.XTType")) XTAnalyzerParameters.XTType = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.XTType");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.AsymXT")) XTAnalyzerParameters.AsymXT = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.AsymXT");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.RequireInTriggerCounter")) XTAnalyzerParameters.RequireInTriggerCounter = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.RequireInTriggerCounter");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.RequireAllGoldenHits")) XTAnalyzerParameters.RequireAllGoldenHits = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.RequireAllGoldenHits");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.FirstGoodPeak")) XTAnalyzerParameters.FirstGoodPeak = MyRuntimeParameters::Get().GetParameterB("XTAnalyzer.FirstGoodPeak");
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
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_min")) XTAnalyzerParameters.bin_t_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_max")) XTAnalyzerParameters.bin_t_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_t_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_t_num")) XTAnalyzerParameters.bin_t_num = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_t_num");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_min")) XTAnalyzerParameters.bin_x_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_max")) XTAnalyzerParameters.bin_x_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.bin_x_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.bin_x_num")) XTAnalyzerParameters.bin_x_num = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.bin_x_num");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.fitX_nRanges")){
            XTAnalyzerParameters.fitX_nRanges = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.fitX_nRanges");
            if (XTAnalyzerParameters.fitX_nRanges > NRANGES){
                MyWarn("XTAnalyzerParameters.fitX_nRanges is set to "<<XTAnalyzerParameters.fitX_nRanges<<" but cannot be larger than "<<NRANGES);
                XTAnalyzerParameters.fitX_nRanges = NRANGES;
            }
        }
        for (int iRange = 0; iRange<XTAnalyzerParameters.fitX_nRanges; iRange++){
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_nBins",iRange))) {XTAnalyzerParameters.fitX_nBins[iRange] = MyRuntimeParameters::Get().GetParameterI(Form("XTAnalyzer.fitX_%d_nBins",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_nBins[jRange] = XTAnalyzerParameters.fitX_nBins[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_smooth",iRange))) {XTAnalyzerParameters.fitX_smooth[iRange] = MyRuntimeParameters::Get().GetParameterI(Form("XTAnalyzer.fitX_%d_smooth",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_smooth[jRange] = XTAnalyzerParameters.fitX_smooth[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_minEntries",iRange))) {XTAnalyzerParameters.fitX_minEntries[iRange] = MyRuntimeParameters::Get().GetParameterI(Form("XTAnalyzer.fitX_%d_minEntries",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_minEntries[jRange] = XTAnalyzerParameters.fitX_minEntries[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_tSep",iRange))) {
                XTAnalyzerParameters.fitX_tSep[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_tSep",iRange));
                if (XTAnalyzerParameters.fitX_tSep[iRange]<XTAnalyzerParameters.bin_t_min){
                    MyWarn("XTAnalyzerParameters.fitX_tSep is set to "<<XTAnalyzerParameters.fitX_tSep[iRange]<<" but cannot be smaller than XTAnalyzerParameters.bin_t_min = "<<XTAnalyzerParameters.bin_t_min);
                    XTAnalyzerParameters.fitX_tSep[iRange]=XTAnalyzerParameters.bin_t_min;
                }
                if (XTAnalyzerParameters.fitX_tSep[iRange]>XTAnalyzerParameters.bin_t_max){
                    MyWarn("XTAnalyzerParameters.fitX_tSep is set to "<<XTAnalyzerParameters.fitX_tSep[iRange]<<" but cannot be larger than XTAnalyzerParameters.bin_t_max = "<<XTAnalyzerParameters.bin_t_max);
                    XTAnalyzerParameters.fitX_tSep[iRange]=XTAnalyzerParameters.bin_t_max;
                }
                for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_tSep[jRange] = XTAnalyzerParameters.fitX_tSep[iRange];}
            }
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_fitBoth",iRange))) {XTAnalyzerParameters.fitX_fitBoth[iRange] = MyRuntimeParameters::Get().GetParameterB(Form("XTAnalyzer.fitX_%d_fitBoth",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_fitBoth[jRange] = XTAnalyzerParameters.fitX_fitBoth[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_SetEmptyBins",iRange))) {XTAnalyzerParameters.fitX_SetEmptyBins[iRange] = MyRuntimeParameters::Get().GetParameterB(Form("XTAnalyzer.fitX_%d_SetEmptyBins",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_SetEmptyBins[jRange] = XTAnalyzerParameters.fitX_SetEmptyBins[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_functionType",iRange))) {XTAnalyzerParameters.fitX_functionType[iRange] = getFunctionType(MyRuntimeParameters::Get().GetParameterS(Form("XTAnalyzer.fitX_%d_functionType",iRange))); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_functionType[jRange] = XTAnalyzerParameters.fitX_functionType[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_height_middle",iRange))) {XTAnalyzerParameters.fitX_peak_height_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_height_middle",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_height_middle[jRange] = XTAnalyzerParameters.fitX_peak_height_middle[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_height_left",iRange))) {XTAnalyzerParameters.fitX_peak_height_left[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_height_left",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_height_left[jRange] = XTAnalyzerParameters.fitX_peak_height_left[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_height_right",iRange))) {XTAnalyzerParameters.fitX_peak_height_right[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_height_right",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_height_right[jRange] = XTAnalyzerParameters.fitX_peak_height_right[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_sigma_middle",iRange))) {XTAnalyzerParameters.fitX_peak_sigma_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_sigma_middle",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_sigma_middle[jRange] = XTAnalyzerParameters.fitX_peak_sigma_middle[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_sigma_left",iRange))) {XTAnalyzerParameters.fitX_peak_sigma_left[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_sigma_left",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_sigma_left[jRange] = XTAnalyzerParameters.fitX_peak_sigma_left[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_sigma_right",iRange))) {XTAnalyzerParameters.fitX_peak_sigma_right[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_sigma_right",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_sigma_right[jRange] = XTAnalyzerParameters.fitX_peak_sigma_right[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_peak_mean_range",iRange))) {XTAnalyzerParameters.fitX_peak_mean_range[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_peak_mean_range",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_mean_range[jRange] = XTAnalyzerParameters.fitX_peak_mean_range[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_height_middle",iRange))) {XTAnalyzerParameters.fitX_base_height_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_height_middle",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_height_middle[jRange] = XTAnalyzerParameters.fitX_base_height_middle[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_height_left",iRange))) {XTAnalyzerParameters.fitX_base_height_left[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_height_left",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_height_left[jRange] = XTAnalyzerParameters.fitX_base_height_left[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_height_right",iRange))) {XTAnalyzerParameters.fitX_base_height_right[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_height_right",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_height_right[jRange] = XTAnalyzerParameters.fitX_base_height_right[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_sigma_middle",iRange))) {XTAnalyzerParameters.fitX_base_sigma_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_sigma_middle",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_sigma_middle[jRange] = XTAnalyzerParameters.fitX_base_sigma_middle[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_sigma_left",iRange))) {XTAnalyzerParameters.fitX_base_sigma_left[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_sigma_left",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_sigma_left[jRange] = XTAnalyzerParameters.fitX_base_sigma_left[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_sigma_right",iRange))) {XTAnalyzerParameters.fitX_base_sigma_right[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_sigma_right",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_sigma_right[jRange] = XTAnalyzerParameters.fitX_base_sigma_right[iRange];}}
            if (MyRuntimeParameters::Get().HasParameter(Form("XTAnalyzer.fitX_%d_base_mean_range",iRange))) {XTAnalyzerParameters.fitX_base_mean_range[iRange] = MyRuntimeParameters::Get().GetParameterD(Form("XTAnalyzer.fitX_%d_base_mean_range",iRange)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_mean_range[jRange] = XTAnalyzerParameters.fitX_base_mean_range[iRange];}}
        }
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_n_min")) XTAnalyzerParameters.graph_n_min = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.graph_n_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_chi2_max")) XTAnalyzerParameters.graph_chi2_max = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.graph_chi2_max");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_prob_min")) XTAnalyzerParameters.graph_prob_min = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.graph_prob_min");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.graph_sepX")) XTAnalyzerParameters.graph_sepX = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.graph_sepX");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_center_nPol")) XTAnalyzerParameters.xt_center_nPol = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.xt_center_nPol");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_center_tLeft")) XTAnalyzerParameters.xt_center_tLeft = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.xt_center_tLeft");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_center_tRight")) XTAnalyzerParameters.xt_center_tRight = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.xt_center_tRight");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_middle_nPol")) XTAnalyzerParameters.xt_middle_nPol = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.xt_middle_nPol");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_middle_tLeft")) XTAnalyzerParameters.xt_middle_tLeft = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.xt_middle_tLeft");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_middle_tRight")) XTAnalyzerParameters.xt_middle_tRight = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.xt_middle_tRight");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_end_nPol")) XTAnalyzerParameters.xt_end_nPol = MyRuntimeParameters::Get().GetParameterI("XTAnalyzer.xt_end_nPol");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_end_tLeft")) XTAnalyzerParameters.xt_end_tLeft = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.xt_end_tLeft");
        if (MyRuntimeParameters::Get().HasParameter("XTAnalyzer.xt_end_tRight")) XTAnalyzerParameters.xt_end_tRight = MyRuntimeParameters::Get().GetParameterD("XTAnalyzer.xt_end_tRight");
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
    FirstGoodPeak = false;
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
    bin_t_min = -25-1/0.96/2; // t range for one x bin
    bin_t_max = 800+1/0.96/2;
    bin_t_num = 792+1;
    bin_x_min = -0.02;
    bin_x_max = 10.02; // x range for one t bin
    bin_x_num = 501;

    fitX_nRanges = 1;
    for (int iRange = 0; iRange<NRANGES; iRange++){
        fitX_minEntries[iRange] = 0; // minimum number of entries in one slice to apply fitting function; Otherwise use mean value & RMS instead.
        fitX_nBins[iRange] = 1;
        fitX_smooth[iRange] = 0; // by default no smoothing
        fitX_tSep[iRange] = 0;
        fitX_fitBoth[iRange] = false;
        fitX_SetEmptyBins[iRange] = false;
        fitX_functionType[iRange] = XTAnalyzer::kGaussian;
        fitX_peak_height_middle[iRange] = 0.9;
        fitX_peak_height_left[iRange] = 0;
        fitX_peak_height_right[iRange] = 2;
        fitX_peak_sigma_middle[iRange] = 0;
        fitX_peak_sigma_left[iRange] = 0;
        fitX_peak_sigma_right[iRange] = 2;
        fitX_peak_mean_range[iRange] = 1;
        fitX_base_height_middle[iRange] = 0.1;
        fitX_base_height_left[iRange] = 0;
        fitX_base_height_right[iRange] = 1;
        fitX_base_sigma_middle[iRange] = 0;
        fitX_base_sigma_left[iRange] = 0.5;
        fitX_base_sigma_right[iRange] = 5;
        fitX_base_mean_range[iRange] = 3;
    }

    graph_n_min = 50;
    graph_chi2_max = 0; // no cut
    graph_prob_min = 0.1;
    graph_sepX = 1;
    xt_center_nPol = 3;
    xt_center_tLeft = 0;
    xt_center_tRight = 100;
    xt_middle_nPol = 3;
    xt_middle_tLeft = 50;
    xt_middle_tRight = 340;
    xt_end_nPol = 3;
    xt_end_tLeft = 330;
    xt_end_tRight = 800;
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
    printf("  FirstGoodPeak = %s\n",FirstGoodPeak?"true":"false");
    printf("  golden hit t_min = %f\n",gold_t_min);
    printf("  golden hit t_max = %f\n",gold_t_max);
    printf(" About 2-D histogram binning :\n");
    printf("  X axis: number of bins = %d\n",bin_x_num);
    printf("  X axis: minimum = %f\n",bin_x_min);
    printf("  X axis: maximum = %f\n",bin_x_max);
    printf("  T axis: number of bins = %d\n",bin_t_num);
    printf("  T axis: minimum = %f\n",bin_t_min);
    printf("  T axis: maximum = %f\n",bin_t_max);
    printf(" About X fitting on T bins :\n");
    printf("  Number of ranges %d\n",fitX_nRanges);
    std::string functionName("");
    for (int iRange = 0; iRange<fitX_nRanges; iRange++){
        if (fitX_functionType[iRange]==XTAnalyzer::kGaussian) functionName = "Gaussian";
        else if (fitX_functionType[iRange]==XTAnalyzer::kLandau) functionName = "Landau";
        else if (fitX_functionType[iRange]==XTAnalyzer::kDoubleGaussian) functionName = "Gaussian+Gaussian";
        else if (fitX_functionType[iRange]==XTAnalyzer::kDoubleLandau) functionName = "Landau+Landau";
        else if (fitX_functionType[iRange]==XTAnalyzer::kGaussianPlusLandau) functionName = "Gaussian+Landau";
        else if (fitX_functionType[iRange]==XTAnalyzer::kLandauPlusGaussian) functionName = "Landau+Gaussian";
        printf("    %d: %.1f ~ %.1f ns, minEntries = %d, nBins = %d, smooth %d, fit both sides together? %s, set empty bins? %s, fit function %s\n",iRange,iRange==0?bin_t_min:fitX_tSep[iRange-1],fitX_tSep[iRange],fitX_minEntries[iRange],fitX_nBins[iRange],fitX_smooth[iRange],fitX_fitBoth[iRange]?"yes":"no",fitX_SetEmptyBins[iRange]?"yes":"no",functionName.c_str());
        printf("        peak part: height (rel to hist) %.2f ~ %.2f ~ %.2f, sigma (mm) %.2f ~ %.2f ~ %.2f, x offset (to hist) range (mm) %.2f\n",fitX_peak_height_left[iRange],fitX_peak_height_middle[iRange],fitX_peak_height_right[iRange],fitX_peak_sigma_left[iRange],fitX_peak_sigma_middle[iRange],fitX_peak_sigma_right[iRange],fitX_peak_mean_range[iRange]);
        if (fitX_functionType[iRange]!=XTAnalyzer::kGaussian&&fitX_functionType[iRange]!=XTAnalyzer::kLandau){
            printf("        base part: height (rel to peak) %.2f ~ %.2f ~ %.2f, rel-sigma (peak sigma) %.2f ~ %.2f ~ %.2f, x offset (to peak) rel-range (peak sigma) %.2f\n",fitX_base_height_left[iRange],fitX_base_height_middle[iRange],fitX_base_height_right[iRange],fitX_base_sigma_left[iRange],fitX_base_sigma_middle[iRange],fitX_base_sigma_right[iRange],fitX_base_mean_range[iRange]);
        }
    }
    printf(" About XT graphs:\n");
    printf("  Minimum number of entries of the sample point to be included in graph: %d\n",graph_n_min);
    printf("  Maximum chi2 of the sample point to be included in graph: %.3e\n",graph_chi2_max);
    printf("  Minimum p-value of the sample point to be included in graph: %.3e\n",graph_prob_min);
    printf("  Separation X to combine space samples and time samples: %.3e\n",graph_sepX);
    printf(" About XT function:\n");
    printf("  XTType = %d\n",XTType);
    printf("  AsymXT = %s\n",AsymXT?"true":"false");
    printf("  center: Pol%d, T %.0f ~ %.0f ns",xt_center_nPol,xt_center_tLeft,xt_center_tRight);
    printf("  middle: Pol%d, T %.0f ~ %.0f ns",xt_middle_nPol,xt_middle_tLeft,xt_middle_tRight);
    printf("  end: Pol%d, T %.0f ~ %.0f ns",xt_end_nPol,xt_end_tLeft,xt_end_tRight);
}

AnaPara::AnaPara(){
}

void AnaPara::Print(){
}

int ParameterManager::getFunctionType(TString name){
    int functionType = 0;
    if (name.Contains("gaus",TString::kIgnoreCase)
      &&name.Contains("land",TString::kIgnoreCase)){
        if (name.Index("gaus",0,TString::kIgnoreCase)<name.Index("land",0,TString::kIgnoreCase)){
            functionType = XTAnalyzer::kGaussianPlusLandau;
        }
        else{
            functionType = XTAnalyzer::kLandauPlusGaussian;
        }
    }
    else if (name.Contains("gaus",TString::kIgnoreCase)){
        if (name.Contains("double",TString::kIgnoreCase)){
            functionType = XTAnalyzer::kDoubleGaussian;
        }
        else{
            functionType = XTAnalyzer::kGaussian;
        }
    }
    else if (name.Contains("land",TString::kIgnoreCase)){
        if (name.Contains("double",TString::kIgnoreCase)){
            functionType = XTAnalyzer::kDoubleLandau;
        }
        else{
            functionType = XTAnalyzer::kLandau;
        }
    }
    return functionType;
}
