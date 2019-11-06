#ifndef XTANALYZER_H
#define XTANALYZER_H

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
class XTBinAnalyzer{
    public:
        enum SampleType{
            kTimeSliceUnfolded,
            kSpaceSliceUnfolded,
            kTimeSliceFolded,
            kSpaceSliceFolded
        };
        enum FunctionType{
            kNone,
            kGaussian,
            kLandau,
            kLandauFlipped,
            kGaussianPlusLandau,
            kGaussianPlusLandauFlipped,
            kDoubleGaussian
        };

        XTBinAnalyzer(TString runname, TFile * outfile, bool savehists = false, bool drawDetails = false);
        virtual ~XTBinAnalyzer(void);

        void SetSaveHists(bool save){mSaveHists = save;}
        void SetDrawDetails(bool draw){mDrawDetails = draw;}
        void SetTestLayer(int lid){mLayerID = lid;}

        int  Initialize(void);
        void Process(void);

        void Fill(double t, double x);

    private:
        void fitSlice(TH1D * h,TCanvas * canv, TString drawTitle, double &mean, double &error, double &sigma, double &chi2, double &left, double &right, int &function, int &previousFunction, double ratio = 0.3, bool tryLandau = false, bool flippedLandau = false);
        void getRange(TH1D * h, int & bmax, int & binl, int & binr, double ratio = 0.3);
        void fitSliceFloat(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right, double ratio = 0.3);
        TF1 * fitSliceGaus(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right);
        TF1 * fitSliceLand(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right);
        TF1 * fitSliceLandF(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right);
        TF1 * fitSlice2Gaus(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right);
        TF1 * myNewTF1(TString name, TString form, double left, double right);
        void drawFitting(TH1D* h,TF1 * f,TCanvas * c,TString title, TString filename,double left, double center, double right);

    private:
        // options
        TString mRunName;
        bool mSaveHists;
        bool mDrawDetails;

        // for output xt tree
        TFile * mOutFile;
        TTree * mOutTree;
        int mLayerID;
        int mCellID;
        double mX;
        double mXerr;
        double mT;
        double mTerr;
        double mSig;
        double mChi2;
        int mEntries;
        int mType;
        int mFunction;

        int    minEntries;
        // about binning
        double m_bin_t_min;
        double m_bin_t_max;
        int    m_bin_t_num;
        double m_bin_x_min;
        double m_bin_x_max;
        int    m_bin_x_num;
        // about projection
        int    m_bin_t_fit_num;
        int    m_bin_t_fit_num_tail;
        double m_bin_t_tailTime;
        int    m_bin_x_fit_num;
        // about the range for using Landau function
        //   time slice to fit space
        double m_bin_t_landTmin;
        double m_bin_t_landTmax;
        //   space slice to fit time
        double m_bin_x_landXmin;
        double m_bin_x_landXmax;
        // about fitting range
        double m_bin_t_ratio; // for the given T slice, get a range on X using this ratio*maximum_height as a threshold, do the fitting on X within this range
        double m_bin_x_ratio; // for the given X slice, get a range on T using this ratio*maximum_height as a threshold, do the fitting on T within this range

        // Histograms for data points
        TH2D * h2_xt;
        TH2D * h2_xtn; // for neutral driftD

        // functions to fit x/t slices
        TF1 * f_gaus;
        TF1 * f_land;
        TF1 * f_landF;
        TF1 * f_2gaus;

};

#endif
