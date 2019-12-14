#include "ParameterManager.hxx"
#include "XTManager.hxx"
#include "XTAnalyzer.hxx"
#include "MyRuntimeParameters.hxx"
#include "HEPUnits.hxx"

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
    std::string parName("");
    parName="geoSetup";
    if (MyRuntimeParameters::Get().HasParameter(parName)){
        TString name = MyRuntimeParameters::Get().GetParameterS(parName);
        name.ToLower();
        if (name == "tilted") geoSetup = GeometryManager::kTilted;
        else if (name == "finger") geoSetup = GeometryManager::kFinger;
        else if (name == "normal") geoSetup = GeometryManager::kNormal;
    }
    parName="chamberType";
    if (MyRuntimeParameters::Get().HasParameter(parName)){
        TString name = MyRuntimeParameters::Get().GetParameterS(parName);
        name.ToLower();
        if (name == "cdc") chamberType = GeometryManager::kCDC;
        else if (name == "proto4") chamberType = GeometryManager::kProto4;
    }
    parName="connectionType";
    if (MyRuntimeParameters::Get().HasParameter(parName)){
        TString name = MyRuntimeParameters::Get().GetParameterS(parName);
        name.ToLower();
        if (name == "spring8") connectionType = GeometryManager::kSPring8;
        else if (name == "cosmic") connectionType = GeometryManager::kCosmic;
    }
    parName="beam";
    if (MyRuntimeParameters::Get().HasParameter(parName)){
        TString name = MyRuntimeParameters::Get().GetParameterS(parName);
        name.ToLower();
        if (name == "spring8") beamType = BeamManager::kSPring8;
        else if (name == "spring8tilted") beamType = BeamManager::kSPring8;
        else if (name == "cosmic") beamType = BeamManager::kCosmic;
    }
    parName="inputHitType";
    if (MyRuntimeParameters::Get().HasParameter(parName)){
        TString name = MyRuntimeParameters::Get().GetParameterS(parName);
        name.ToLower();
        if (name == "data") inputHitType = InputOutputManager::kData;
        else if (name == "mc" || name == "mcdriftt") inputHitType = InputOutputManager::kMCDriftT;
        else if (name == "mcdriftd") inputHitType = InputOutputManager::kMCDriftD;
    }
    parName="peakType";
    if (MyRuntimeParameters::Get().HasParameter(parName)){
        TString name = MyRuntimeParameters::Get().GetParameterS(parName);
        name.ToLower();
        if (name == "firstpeak") peakType = TrackingPara::kFirstPeak;
        else if (name == "highestpeak") peakType = TrackingPara::kHighestPeak;
        else if (name == "allpeaks") peakType = TrackingPara::kAllPeaks;
    }
    if (theParaBlock==kAll||theParaBlock==kTracking){
        parName="tracking.nHitsGMax";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.nHitsGMax = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.nHitsSMin";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.nHitsSMin = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.nPairsMin";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.nPairsMin = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.t0shift0";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.t0shift0 = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.t0shift1";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.t0shift1 = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.tmin";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.tmin = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.tmax";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.tmax = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.sumCut";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.sumCut = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="tracking.aaCut";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.aaCut = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="tracking.t0error";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.t0error = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="tracking.lidStart";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.lidStart = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="tracking.lidStop";if (MyRuntimeParameters::Get().HasParameter(parName)) TrackingParameters.lidStop = MyRuntimeParameters::Get().GetParameterI(parName);
    }
    if (theParaBlock==kAll||theParaBlock==kXTManager){
        parName="XTManager.xtType";if (MyRuntimeParameters::Get().HasParameter(parName)) XTManagerParameters.xtType = getXTType(MyRuntimeParameters::Get().GetParameterS(parName));
        parName="XTManager.defaultLayer";if (MyRuntimeParameters::Get().HasParameter(parName)) XTManagerParameters.defaultLayer = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTManager.evenLayer";if (MyRuntimeParameters::Get().HasParameter(parName)) XTManagerParameters.evenLayer = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTManager.oddLayer";if (MyRuntimeParameters::Get().HasParameter(parName)) XTManagerParameters.oddLayer = MyRuntimeParameters::Get().GetParameterI(parName);
    }
    if (theParaBlock==kAll||theParaBlock==kXTAnalyzer){
        parName="XTAnalyzer.CandSelBy";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.CandSelBy = MyRuntimeParameters::Get().GetParameterS(parName);
        parName="XTAnalyzer.RequireInTriggerCounter";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.RequireInTriggerCounter = MyRuntimeParameters::Get().GetParameterB(parName);
        parName="XTAnalyzer.RequireAllGoldenHits";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.RequireAllGoldenHits = MyRuntimeParameters::Get().GetParameterB(parName);
        parName="XTAnalyzer.FirstGoodPeak";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.FirstGoodPeak = MyRuntimeParameters::Get().GetParameterB(parName);
        parName="XTAnalyzer.UseGoodHit";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.UseGoodHit = MyRuntimeParameters::Get().GetParameterB(parName);
        parName="XTAnalyzer.AllGoodHitsUsed";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.AllGoodHitsUsed = MyRuntimeParameters::Get().GetParameterB(parName);
        parName="XTAnalyzer.nHits_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.nHits_max = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTAnalyzer.nHitsS_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.nHitsS_min = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTAnalyzer.chi2_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.chi2_max = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.pValue_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.pValue_min = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.slz_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.slz_min = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.slz_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.slz_max = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.gold_t_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.gold_t_min = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.gold_t_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.gold_t_max = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.bin_t_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.bin_t_min = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.bin_t_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.bin_t_max = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.bin_t_num";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.bin_t_num = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTAnalyzer.bin_x_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.bin_x_min = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.bin_x_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.bin_x_max = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.bin_x_num";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.bin_x_num = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTAnalyzer.fitX_nRanges";
        if (MyRuntimeParameters::Get().HasParameter(parName)){
            XTAnalyzerParameters.fitX_nRanges = MyRuntimeParameters::Get().GetParameterI(parName);
            if (XTAnalyzerParameters.fitX_nRanges > NRANGES){
                MyWarn("XTAnalyzerParameters.fitX_nRanges is set to "<<XTAnalyzerParameters.fitX_nRanges<<" but cannot be larger than "<<NRANGES);
                XTAnalyzerParameters.fitX_nRanges = NRANGES;
            }
        }
        for (int iRange = 0; iRange<XTAnalyzerParameters.fitX_nRanges; iRange++){
            parName=Form("XTAnalyzer.fitX_%d_nBins",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_nBins[iRange] = MyRuntimeParameters::Get().GetParameterI(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_nBins[jRange] = XTAnalyzerParameters.fitX_nBins[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_smooth",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_smooth[iRange] = MyRuntimeParameters::Get().GetParameterI(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_smooth[jRange] = XTAnalyzerParameters.fitX_smooth[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_minEntries",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_minEntries[iRange] = MyRuntimeParameters::Get().GetParameterI(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_minEntries[jRange] = XTAnalyzerParameters.fitX_minEntries[iRange];}}
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
            parName=Form("XTAnalyzer.fitX_%d_fitBoth",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_fitBoth[iRange] = MyRuntimeParameters::Get().GetParameterB(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_fitBoth[jRange] = XTAnalyzerParameters.fitX_fitBoth[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_SetEmptyBins",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_SetEmptyBins[iRange] = MyRuntimeParameters::Get().GetParameterB(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_SetEmptyBins[jRange] = XTAnalyzerParameters.fitX_SetEmptyBins[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_functionType",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_functionType[iRange] = getFunctionType(MyRuntimeParameters::Get().GetParameterS(parName)); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_functionType[jRange] = XTAnalyzerParameters.fitX_functionType[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_height_middle",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_height_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_height_middle[jRange] = XTAnalyzerParameters.fitX_peak_height_middle[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_height_left",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_height_left[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_height_left[jRange] = XTAnalyzerParameters.fitX_peak_height_left[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_height_right",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_height_right[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_height_right[jRange] = XTAnalyzerParameters.fitX_peak_height_right[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_sigma_middle",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_sigma_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_sigma_middle[jRange] = XTAnalyzerParameters.fitX_peak_sigma_middle[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_sigma_left",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_sigma_left[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_sigma_left[jRange] = XTAnalyzerParameters.fitX_peak_sigma_left[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_sigma_right",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_sigma_right[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_sigma_right[jRange] = XTAnalyzerParameters.fitX_peak_sigma_right[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_peak_mean_range",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_peak_mean_range[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_peak_mean_range[jRange] = XTAnalyzerParameters.fitX_peak_mean_range[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_height_middle",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_height_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_height_middle[jRange] = XTAnalyzerParameters.fitX_base_height_middle[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_height_left",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_height_left[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_height_left[jRange] = XTAnalyzerParameters.fitX_base_height_left[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_height_right",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_height_right[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_height_right[jRange] = XTAnalyzerParameters.fitX_base_height_right[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_sigma_middle",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_sigma_middle[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_sigma_middle[jRange] = XTAnalyzerParameters.fitX_base_sigma_middle[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_sigma_left",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_sigma_left[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_sigma_left[jRange] = XTAnalyzerParameters.fitX_base_sigma_left[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_sigma_right",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_sigma_right[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_sigma_right[jRange] = XTAnalyzerParameters.fitX_base_sigma_right[iRange];}}
            parName=Form("XTAnalyzer.fitX_%d_base_mean_range",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.fitX_base_mean_range[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.fitX_nRanges; jRange++){XTAnalyzerParameters.fitX_base_mean_range[jRange] = XTAnalyzerParameters.fitX_base_mean_range[iRange];}}
        }
        // get the last tSep
        int iRange = XTAnalyzerParameters.fitX_nRanges;
        parName=Form("XTAnalyzer.fitX_%d_tSep",iRange);
        if (MyRuntimeParameters::Get().HasParameter(parName)) {
            XTAnalyzerParameters.fitX_tSep[iRange] = MyRuntimeParameters::Get().GetParameterD(parName);
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
        parName="XTAnalyzer.graph_n_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.graph_n_min = MyRuntimeParameters::Get().GetParameterI(parName);
        parName="XTAnalyzer.graph_chi2_max";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.graph_chi2_max = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.graph_prob_min";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.graph_prob_min = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.graph_sepX";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.graph_sepX = MyRuntimeParameters::Get().GetParameterD(parName);

        parName="XTAnalyzer.xtfunc_nRanges";
        if (MyRuntimeParameters::Get().HasParameter(parName)){
            XTAnalyzerParameters.xtfunc_nRanges = MyRuntimeParameters::Get().GetParameterI(parName);
            if (XTAnalyzerParameters.xtfunc_nRanges > NRANGES){
                MyWarn("XTAnalyzerParameters.xtfunc_nRanges is set to "<<XTAnalyzerParameters.xtfunc_nRanges<<" but cannot be larger than "<<NRANGES);
                XTAnalyzerParameters.xtfunc_nRanges = NRANGES;
            }
        }
        parName="XTAnalyzer.xtfunc_tHighEdge";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.xtfunc_tHighEdge = MyRuntimeParameters::Get().GetParameterD(parName);
        for (int iRange = 0; iRange<XTAnalyzerParameters.xtfunc_nRanges; iRange++){
            parName=Form("XTAnalyzer.xtfunc_%d_nPol",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.xtfunc_nPol[iRange] = MyRuntimeParameters::Get().GetParameterI(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.xtfunc_nRanges; jRange++){XTAnalyzerParameters.xtfunc_nPol[jRange] = XTAnalyzerParameters.xtfunc_nPol[iRange];}}
            parName=Form("XTAnalyzer.xtfunc_%d_tLeft",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.xtfunc_tLeft[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.xtfunc_nRanges; jRange++){XTAnalyzerParameters.xtfunc_tLeft[jRange] = XTAnalyzerParameters.xtfunc_tLeft[iRange];}}
            parName=Form("XTAnalyzer.xtfunc_%d_tRight",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.xtfunc_tRight[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.xtfunc_nRanges; jRange++){XTAnalyzerParameters.xtfunc_tRight[jRange] = XTAnalyzerParameters.xtfunc_tRight[iRange];}}
            parName=Form("XTAnalyzer.xtfunc_%d_tLowEdge",iRange);if (MyRuntimeParameters::Get().HasParameter(parName)) {XTAnalyzerParameters.xtfunc_tLowEdge[iRange] = MyRuntimeParameters::Get().GetParameterD(parName); for (int jRange = iRange+1; jRange<XTAnalyzerParameters.xtfunc_nRanges; jRange++){XTAnalyzerParameters.xtfunc_tLowEdge[jRange] = XTAnalyzerParameters.xtfunc_tLowEdge[iRange];}}
        }
        parName="XTAnalyzer.draw_tmin";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.draw_tmin = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.draw_tmax";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.draw_tmax = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.draw_xmin";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.draw_xmin = MyRuntimeParameters::Get().GetParameterD(parName);
        parName="XTAnalyzer.draw_xmax";if (MyRuntimeParameters::Get().HasParameter(parName)) XTAnalyzerParameters.draw_xmax = MyRuntimeParameters::Get().GetParameterD(parName);
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
    XTManagerParameters.Print();
    XTAnalyzerParameters.Print();
    AnaParameters.Print();
}

TrackingPara::TrackingPara(){
    nHitsGMax = 30;
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
    printf("  nHitsGMax   = %d\n",nHitsGMax);
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

XTManagerPara::XTManagerPara(){
    xtType = XTManager::kSingleFolded;
    defaultLayer = 4;
    evenLayer = 4;
    oddLayer = 5;
}

void XTManagerPara::Print(){
    printf("XTManager Parameters:\n");
    TString xttypename("");
    if (xtType==XTManager::kSingleFolded) xttypename = "SingleFolded";
    else if (xtType==XTManager::kSingleLeftRight) xttypename = "SingleLeftRight";
    else if (xtType==XTManager::kAllFolded) xttypename = "AllFolded";
    else if (xtType==XTManager::kAllLeftRight) xttypename = "AllLeftRight";
    else if (xtType==XTManager::kEvenOddFolded) xttypename = "EvenOddFolded";
    else if (xtType==XTManager::kEvenOddLeftRight) xttypename = "EvenOddLeftRight";
    else xttypename = "UNRECOGONISED";
    printf(" xtType = %d (%s)\n",xtType,xttypename.Data());
    printf(" Default Layer: %d\n",defaultLayer);
    printf(" Even Layer: %d\n",evenLayer);
    printf(" Odd Layer: %d\n",oddLayer);
}

XTAnalyzerPara::XTAnalyzerPara(){
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

    xtfunc_nRanges = 1;
    xtfunc_tHighEdge = 0;
    for (int iRange = 0; iRange<NRANGES; iRange++){
        xtfunc_nPol[iRange] = 1;
        xtfunc_tLeft[iRange] = 0;
        xtfunc_tRight[iRange] = 0;
        xtfunc_tLowEdge[iRange] = 0;
    }

    draw_tmin = -25*unit::ns;
    draw_tmax = 800*unit::ns;
    draw_xmin = -10*unit::mm;
    draw_xmax = 10*unit::mm;
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
        else if (fitX_functionType[iRange]==XTAnalyzer::kOptimal) functionName = "Optimal: Gaussian+Landau or Landau+Gaussian";
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
    printf(" About XT functions fitting: valid range %.1f ~ %.1f ns\n",xtfunc_tLowEdge[0],xtfunc_tHighEdge);
    printf("  Number of ranges %d\n",xtfunc_nRanges);
    for (int iRange = 0; iRange<xtfunc_nRanges; iRange++){
        printf("    %d: %.1f ~ %.1f ns, nPol = %d, connect to previous at %.1f ns\n",iRange,xtfunc_tLeft[iRange],xtfunc_tRight[iRange],xtfunc_nPol[iRange],xtfunc_tLowEdge[iRange]);
    }
    printf(" About drawing samples and fitting results: t %.1f ~ %.1f ns, x %.1f ~ %.1f mm\n",draw_tmin/unit::ns,draw_tmax/unit::ns,draw_xmin/unit::mm,draw_xmax/unit::mm);
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
    else if (name.Contains("optimal",TString::kIgnoreCase)){
        functionType = XTAnalyzer::kOptimal;
    }
    return functionType;
}

int ParameterManager::getXTType(TString name){
    int xtType = 0;
    if (name.Contains("single",TString::kIgnoreCase)){
        if (name.Contains("fold",TString::kIgnoreCase)){
            xtType = XTManager::kSingleFolded;
        }
        else{
            xtType = XTManager::kSingleLeftRight;
        }
    }
    else if (name.Contains("all",TString::kIgnoreCase)){
        if (name.Contains("fold",TString::kIgnoreCase)){
            xtType = XTManager::kAllFolded;
        }
        else{
            xtType = XTManager::kAllLeftRight;
        }
    }
    else if (name.Contains("even",TString::kIgnoreCase)&&name.Contains("odd",TString::kIgnoreCase)){
        if (name.Contains("fold",TString::kIgnoreCase)){
            xtType = XTManager::kEvenOddFolded;
        }
        else{
            xtType = XTManager::kEvenOddLeftRight;
        }
    }
    return xtType;
}
