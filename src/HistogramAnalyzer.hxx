#ifndef HISTOGRAMANALYZER_H
#define HISTOGRAMANALYZER_H

#include "TString.h"

class TH1D;
class TF1;
class TCanvas;

class HistogramAnalyzer{
public:
    enum FunctionType{
        kNone = -3, // keep this at first
        kOptimal = -2,
        kBasic = -1, // keep this before basic functions and keep it as -1
        kGaussian,
        kLandau,
        kMultiple, // keep this between basic functions and multiple functions
        kDoubleGaussian,
        kDoubleLandau,
        kGaussianPlusLandau,
        kLandauPlusGaussian,
        kBothSides, // keep this between single side and double sides functions
        // keep the same order as single side functions
        kGaussianBothSides,
        kLandauBothSides,
        kMultipleBothSides, // keep this between basic functions and multiple functions
        kDoubleGaussianBothSides,
        kDoubleLandauBothSides,
        kGaussianPlusLandauBothSides,
        kLandauPlusGaussianBothSides
    };

    HistogramAnalyzer(void);
    virtual ~HistogramAnalyzer(void);

    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static HistogramAnalyzer& Get(void) {
        if (!fHistogramAnalyzer)
            fHistogramAnalyzer = new HistogramAnalyzer();
        return *fHistogramAnalyzer;
    }

    void SetFittingParameters(int iRange);
    int  FitSlice(TH1D * h, double & chi2, double & prob, bool fitBoth, bool isLeft);
    void DrawFitting(TH1D* h);

    // Getters
    int    get_result(void){return fit_result;};
    double get_chi2(void){return fit_chi2;};
    double get_pValue(void){return fit_pValue;};
    double get_x(void){return fit_x;};
    double get_sig(void){return fit_sig;};
    double get_xerr(void){return fit_xerr;};
    double get_xL(void){return fit_xL;};
    double get_sigL(void){return fit_sigL;};
    double get_xerrL(void){return fit_xerrL;};
    double get_xR(void){return fit_xR;};
    double get_sigR(void){return fit_sigR;};
    double get_xerrR(void){return fit_xerrR;};
    int    get_functionType(void){return cur_functionType;};

private:
    bool isFittingGood(TF1 * f);
    void getMeanRMS(TF1 * f, const TH1D * hist, double & mean, double & sigma);
    bool chooseCurrentFunctions(int functionType, bool fitBoth, bool isLeft);
    int  doFitSlice(TH1D * hist);
    void setFittingResults(TH1D * hist,bool isLeft);

private:
    /// The static pointer to the singleton instance.
    static HistogramAnalyzer * fHistogramAnalyzer;

    // parameters (instructions from the user)
    double par_peak_height_middle;
    double par_peak_height_left;
    double par_peak_height_right;
    double par_peak_sigma_middle;
    double par_peak_sigma_left;
    double par_peak_sigma_right;
    double par_peak_mean_range;
    double par_base_height_middle;
    double par_base_height_left;
    double par_base_height_right;
    double par_base_sigma_middle;
    double par_base_sigma_left;
    double par_base_sigma_right;
    double par_base_mean_range;
    int    par_functionType; // can be kOptimal

    // functions to fit x/t slices
    TF1 * f_gausL; // simple or peak component
    TF1 * f_gausR; // simple or peak component
    TF1 * f_gaus2L; // base component
    TF1 * f_gaus2R; // base component
    TF1 * f_gausBoth; // simple
    TF1 * f_landL; // simple or peak component
    TF1 * f_landR; // simple or peak component
    TF1 * f_land2L; // base component
    TF1 * f_land2R; // base component
    TF1 * f_landBoth; // simple
    TF1 * f_combDoubleGausL; // multiple
    TF1 * f_combDoubleGausR; // multiple
    TF1 * f_combDoubleGausBoth; // multiple
    TF1 * f_combDoubleLandL; // multiple
    TF1 * f_combDoubleLandR; // multiple
    TF1 * f_combDoubleLandBoth; // multiple
    TF1 * f_combGausLandL; // multiple
    TF1 * f_combGausLandR; // multiple
    TF1 * f_combGausLandBoth; // multiple
    TF1 * f_combLandGausL; // multiple
    TF1 * f_combLandGausR; // multiple
    TF1 * f_combLandGausBoth; // multiple

    // to record the current status
    TF1 * cur_function;
    TF1 * cur_peak_function;
    TF1 * cur_base_function;
    TF1 * cur_functionL;
    TF1 * cur_peak_functionL;
    TF1 * cur_base_functionL;
    TF1 * cur_functionR;
    TF1 * cur_peak_functionR;
    TF1 * cur_base_functionR;
    int   cur_functionType;

    // current fitting results
    int    fit_result;
    double fit_chi2;
    double fit_pValue;
    double fit_x;
    double fit_sig;
    double fit_xerr;
    double fit_xL;
    double fit_sigL;
    double fit_xerrL;
    double fit_xR;
    double fit_sigR;
    double fit_xerrR;
};

#endif
