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
    XTAnalyzer(TString runname, TFile * outfile, bool drawDetails = false);
    virtual ~XTAnalyzer(void);

    void SetDrawDetails(bool draw){mDrawDetails = draw;}
    void SetSuffix(TString suf){m_suffix = suf;}

    void Initialize(void);
    int  Prepare2DHists(bool reLoad=false); ///< If reLoad is set to true, then the objects will be loaded instead of being created
    int  PrepareTree(bool reLoad=false); ///< If reLoad is set to true, then the objects will be loaded instead of being created
    int  PrepareXTFunctions();
    void Fill(double t, double x);
    void Write(void);
    void BinAnalysis(void);
    void FitXT(void);

private:
    TF1 * myNewTF1(TString name, TString form, double left, double right);
    void drawSample2D(bool withFunction = false);
    void drawSampleCombined2D(void);
    void drawSampleAtt(void);
    void formXTGraphs(void); //< scan through the tree and select points to form up the graphs.
    void doFitXT(TGraph * gr, TF1 * f, const char * suf);
    void combineLeftAndRight(double offset);

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
    TH2D * h2_xt_combined;
    TGraphErrors * gr_left;
    TGraphErrors * gr_right;
    TGraphErrors * gr_rightMinusLeft;
    TGraphErrors * gr_combined;

    // XT range for drawing
    double mDrawTmin;
    double mDrawTmax;
    double mDrawXmin;
    double mDrawXmax;

    // functions to fit XT relation
    TF1 * f_basicXT[NRANGE][NPOL];
    TF1 * f_left;
    TF1 * f_right;
    TF1 * f_combinedLeft;
    TF1 * f_combinedRight;
};

#endif
