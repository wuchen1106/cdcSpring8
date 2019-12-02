#ifndef XTANALYZER_H
#define XTANALYZER_H

#define NRANGE 10
#define NPOL   10

#include <vector>

#include "TString.h"

class TFile;
class TF1;
class TGraph;
class TLine;
class TCanvas;
class TTree;
class TH1D;
class TH2D;
class TDirectory;
class TGraphErrors;

/// This is a class to collect XT relation samples and to analyze them bin by bin
///
/// To use this class, prepare a writable ROOT file as output file first, construct and tell this class what whether to save and draw histograms. \n
/// Then everytime you want to collect a new set of XT relation samples, just call Initialize(layerID) first (where 2D histograms got created), and then call Fill(t,x) consequtively. \n
/// These 2-D histograms include an unfolded one and a foled one. \n
/// After collection, call Process() to analyze the collected relation samples bin-by-bin. \n
/// This analysis including anlaysis on T slices and T slices. \n
/// Eventually analysis results will be saved in a TTree object in the output file named "XTBins" \n
/// Notice: the output tree is created in the constructor. After that it only keeps pushing in new bins but never reset again. At last it will be saved into the output file regardless of the existence of a duplicated tree already saved there. \n
/// TODO: add support to cell specific XT relation analysis
class XTAnalyzer{
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

        XTAnalyzer(TString runname, TFile * outfile, bool drawDetails = false);
        virtual ~XTAnalyzer(void);

        void SetDrawDetails(bool draw){mDrawDetails = draw;}
        void SetSuffix(TString suf){m_suffix = suf;}

        void Initialize(void);
        int  Prepare2DHists(bool reLoad=false); ///< If reLoad is set to true, then the objects will be loaded instead of being created
        int  PrepareTree(bool reLoad=false); ///< If reLoad is set to true, then the objects will be loaded instead of being created
        int  PrepareXTFunctions();
        void Fill(double t, double x);
        void BinAnalysis(void);
        void FitXT(void);
        void Write(void);

    private:
        TF1 * fitSliceBothSides(TH1D * h, double & x1,double & xerr1,double & sig1,double & x2,double & xerr2,double & sig2,double & chi2,double & prob, int & result, int & functionType, int iRange);
        TF1 * fitSliceSingleSide(TH1D * h, double & x,double & xerr,double & sig,double & chi2,double & prob, int & result, int & functionType, int iRange, bool isLeft);
        void getMaximum(TH1D * h, double & position, double & maximum, double left, double right);
        void plusGraph(TGraphErrors * gn, const TGraphErrors * gl, const TGraphErrors * gr, double sl = 1, double sr = 1);
        double interpolate(const TGraphErrors * graph, double theX);
        TF1 * myNewTF1(TString name, TString form, double left, double right);
        void drawFitting(TH1D* h,TCanvas * c,TString title, TString filename, int function, double center1, double center2, bool isLeft = false);
        void drawSamples(void);
        void formXTGraphs(void); //< scan through the tree and select points to form up the graphs.
        void getMeanRMS(TF1 * f, const TH1D * h, double & mean, double & sigma);

    private:
        // options
        TString mRunName;
        bool mDrawDetails;

        // for output xt tree
        TFile * mOutFile;
        TTree * mOutTree;
        TString m_suffix;
        double mX;
        double mXerr;
        double mT;
        double mTerr;
        double mSig;
        double mChi2;
        double mProb;
        int mEntries;
        int mFunction;

        // Histograms for data points
        TH2D * h2_xt;
        TGraphErrors * gr_left;
        TGraphErrors * gr_right;
        TGraphErrors * gr_rightMinusLeft;

        // XT range for drawing
        double mDrawTmin;
        double mDrawTmax;
        double mDrawXmin;
        double mDrawXmax;

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

        // functions to fit XT relation
        TF1 * f_basicXT[NRANGE][NPOL];
        TF1 * f_left;
        TF1 * f_right;
};

#endif
