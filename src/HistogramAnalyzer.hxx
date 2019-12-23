#ifndef HISTOGRAMANALYZER_H
#define HISTOGRAMANALYZER_H

#include "TString.h"

class TH1D;
class TF1;
class TCanvas;

class HistogramAnalyzer{
public:
    enum FunctionType{
        kNone,
        kOptimal,
        kGaussian,
        kLandau,
        kDoubleGaussian,
        kDoubleLandau,
        kGaussianPlusLandau,
        kGaussianBothSides,
        kLandauBothSides,
        kDoubleGaussianBothSides,
        kDoubleLandauBothSides,
        kGaussianPlusLandauBothSides,
        kLandauPlusGaussian,
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

    TF1 * FitSliceBothSides(TH1D * h, double & x1,double & xerr1,double & sig1,double & x2,double & xerr2,double & sig2,double & chi2,double & prob, int & result, int & functionType, int iRange);
    TF1 * FitSliceSingleSide(TH1D * h, double & x,double & xerr,double & sig,double & chi2,double & prob, int & result, int & functionType, int iRange, bool isLeft);
    void DrawFitting(TH1D* h,TCanvas * c,TString title, TString filename, int function, double center1, double center2, bool isLeft = false);

private:
    bool isFittingGood(TF1 * f);
    void getMeanRMS(TF1 * f, const TH1D * hist, double & mean, double & sigma);

private:
    /// The static pointer to the singleton instance.
    static HistogramAnalyzer * fHistogramAnalyzer;

    // functions to fit x/t slices
    TF1 * f_gausL;
    TF1 * f_gausR;
    TF1 * f_gaus2L;
    TF1 * f_gaus2R;
    TF1 * f_gausBoth;
    TF1 * f_landL;
    TF1 * f_landR;
    TF1 * f_land2L;
    TF1 * f_land2R;
    TF1 * f_landBoth;
    TF1 * f_combDoubleGausL;
    TF1 * f_combDoubleGausR;
    TF1 * f_combDoubleGausBoth;
    TF1 * f_combDoubleLandL;
    TF1 * f_combDoubleLandR;
    TF1 * f_combDoubleLandBoth;
    TF1 * f_combGausLandL;
    TF1 * f_combGausLandR;
    TF1 * f_combGausLandBoth;
    TF1 * f_combLandGausL;
    TF1 * f_combLandGausR;
    TF1 * f_combLandGausBoth;
};

#endif
