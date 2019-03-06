#ifndef XTANALYZER_H
#define XTANALYZER_H

#include <vector>

#include "TString.h"

#define NSLICEX 401 // 401 bins from -8.02 mm to 8.02 mm, binning every 40 um
//#define NSLICEX 201 // 201 bins from -8.04 mm to 8.04 mm, binning every 80 um
#define NSLICET 513 // 257 bins from -1.5625 ns to 801.5625 ns including left and right side, bining every 3.125 ns (3 TDC)
//#define NSLICET 309 // 155 bins from -2.60416 ns to 802.60416 ns including left and right side, bining every 5.2083 ns (5 TDC)
//#define NSLICET 221 // 110 bins from -3.64583 ns to 803.64583 ns including left and right side, bining every 7.2916 ns (7 TDC)
//#define NSLICET 171 // 86 bins from -4.6875 ns to 804.6875 ns including left and right side, bining every 9.375 ns (9 TDC)

class TFile;
class TF1;
class TGraph;
class TLine;
class TCanvas;
class TTree;
class TH1D;
class TH2D;

class XTAnalyzer{
	public:
		XTAnalyzer(int gasID, int debug = 0);
		virtual ~XTAnalyzer(void);

		void SetXTType(int type);
		void SetSaveHists(int save);

		int  Initialize(TString runname, int lid, TFile * infile, TFile * outfile, TTree * otree, int xttype, bool sym, int savehists, bool drawDetails, bool saveXT0 = false,int saveOddEven = 0, bool updateXT = true);
		void Process(void);

		void Push(double t, double x);

		int t2d(double t, double & d, bool isRight);

        void SetBinning(int nt, double tmin, double tmax, int nx, double xmin, double xmax);

	private:

		int t2i(double t, bool positive);
		int x2i(double x);
		void i2x(int i,double & fdmin, double & fdmid, double & fdmax);
		void i2t(int i,double & dtmin, double & dtmid, double & dtmax);
		double getMean(TH1D * h, TH1D * hy, double left, double right);
		void fitSliceHistFloat(TH1D * h, double ratio, double & mean, double & sigma, double & chi2, double & left, double & right);
		TF1 * fitSliceGaus(TH1D * h, double & mean, double & sigma, double & chi2, double & left, double & right);
		TF1 * fitSliceLand(TH1D * h, double & mean, double & sigma, double & chi2, double & left, double & right);
		TF1 * fitSlice2Gaus(TH1D * h, double & mean, double & sigma, double & chi2, double & left, double & right);
		double findFirstZero(TF1 * f, double xmin, double xmax, double delta);
		void getTT(int & iT, std::vector<double> & vx, std::vector<double> & vt, std::vector<double> & vsig, std::vector<double> & vn, bool negtive = false);
		TF1 * myNewTF1(TString name, TString form, double left, double right);
		TGraph * myNewTGraph(TString name, int n, const double * x, const double * y, TString title, TString Xtitle, TString Ytitle, int mType, double mSize, int mColor, double lSize, int lColor);
		TF1 * scalePolN(TF1 * f, double factor, TString name = "f_test");
		TF1 * minusPolN(TString name, TF1 * f1, TF1 * f2, double xmin, double xmax);
		TF1 * combinePolN(TString name, TF1 * f1, TF1 * f2, double x0, double x1, double x2, double xmin, double xmax);
		TF1 * combinePolN(TString name, TF1 * f1, TF1 * f2, TF1 * f3, double x0, double x1, double x2, double x3, double xmin, double xmax);
		TString formPolN(int start,int n);
		void createGraphs();
		void drawFitting(TH1D* h,TF1 * f,TCanvas * c,TString title, TString filename,double left, double center, double right);
		void drawSamplingsLR();
		void drawSamplingsB();
		void drawLRB();
		void drawIteration();
		void writeObjects();
		void sigmaXReset();
		void sigmaXIncrement(double t, double sig, double n, std::vector<double> & vt, std::vector<double> & vs);
		void sigmaXFinalcheck(std::vector<double> & vt, std::vector<double> & vs);

	private:
		// options
		int mGasID;
		int mDebugLevel;
		int mSaveHists;
		bool mDrawDetails;
		bool mSaveXT0;
		int mSaveXTEO;
		bool mUpdateXT;
		bool mSymmetric;
		int mXTType;
		int mCentPolN;
		int mMidPolN;
		int mEndPolN;
		TString mRunName;
		TString mEOsuffix;

		int mEntriesMin;
		double mSigXmax;
		double mSigTmax;

		// input file
		TFile * mInFile;

		// for output xt tree
		TFile * mOutFile;
		TTree * mOutTree;
		double mX;
		double mT;
		int mLayerID;
		double mSig;
		double mChi2;
		int mEntries;
		int mType;

		// about slicing
		double mBWT;
		double mBWX;
		double mXLEFT;
		double mXRIGHT;
		double mTLEFT;
		double mTRIGHT;

		// about binning
		double mTmin;
		double mTmax;
		double mNbint;

		double mXmin;
		double mXmax;
		double mNbinx;

		// Histograms for data points
		TH2D * h2_xt;
		TH2D * h2_xtn; // for neutral driftD

		// Histograms for bin analysis
		TH1D * h_t[NSLICEX];
		TH1D * h_x[NSLICET];
		TH1D * h_mx[NSLICET]; // to record pure negative x for landau fitting;
		TH1D * h_tn[NSLICEX]; // for neutral driftD
		TH1D * h_xn[NSLICET]; // for neutral driftD
		TH1D * h_mxn[NSLICET]; // to record pure negative x for landau fitting;
		// to record bin centers
		TH1D * h_t_xsum[NSLICEX];
		TH1D * h_x_tsum[NSLICET];
		TH1D * h_tn_xsum[NSLICEX]; // for neutral driftD
		TH1D * h_xn_tsum[NSLICET]; // for neutral driftD

		// functions to fit x/t slices
		TF1 * f_gaus;
		TF1 * f_land;
		TF1 * f_2gaus;

		// about drawing the histograms
		TLine * l_left;
		TLine * l_center;
		TLine * l_right;

		// vectors to cantain sample points
		std::vector<double> v_x_slicex;
		std::vector<double> v_t_slicex;
		std::vector<double> v_sig_slicex;
		std::vector<double> v_n_slicex;
		std::vector<double> v_chi2_slicex;
		std::vector<double> v_x_slicet;
		std::vector<double> v_t_slicet;
		std::vector<double> v_sig_slicet;
		std::vector<double> v_n_slicet;
		std::vector<double> v_chi2_slicet;
		// for neutral driftD
		std::vector<double> v_x_slicexn;
		std::vector<double> v_t_slicexn;
		std::vector<double> v_sig_slicexn;
		std::vector<double> v_n_slicexn;
		std::vector<double> v_chi2_slicexn;
		std::vector<double> v_x_slicetn;
		std::vector<double> v_t_slicetn;
		std::vector<double> v_sig_slicetn;
		std::vector<double> v_n_slicetn;
		std::vector<double> v_chi2_slicetn;

		std::vector<double> v_t_slicetls;
		std::vector<double> v_t_slicetrs;
		std::vector<double> v_t_slicetns;
		std::vector<double> v_sig_slicetls;
		std::vector<double> v_sig_slicetrs;
		std::vector<double> v_sig_slicetns;

		double m_sig_sel;
		double m_t_sel;
		double m_n_sel;
		int m_i_sel;

		// graphs for slice analysis
		TGraph * gr_xt_slicet;
		TGraph * gr_xt_slicex;
		TGraph * gr_xn_slicex;
		TGraph * gr_xsig_slicex;
		TGraph * gr_xchi2_slicex;
		TGraph * gr_nt_slicetl;
		TGraph * gr_sigt_slicetl;
		TGraph * gr_chi2t_slicetl;
		TGraph * gr_nt_slicetr;
		TGraph * gr_sigt_slicetr;
		TGraph * gr_chi2t_slicetr;
		TGraph * gr_xt_slicetn;
		TGraph * gr_xt_slicexn;
		TGraph * gr_xn_slicexn;
		TGraph * gr_xsig_slicexn;
		TGraph * gr_xchi2_slicexn;
		TGraph * gr_nt_slicetn;
		TGraph * gr_sigt_slicetn;
		TGraph * gr_chi2t_slicetn;

		TGraph * gr_sigts_slicetl;
		TGraph * gr_sigts_slicetr;
		TGraph * gr_sigts_slicetl0;
		TGraph * gr_sigts_slicetr0;
		TGraph * gr_sigts_slicetlEO;
		TGraph * gr_sigts_slicetrEO;

		// vectors selected for XT graphs
		std::vector<double> v_left_cen_x;
		std::vector<double> v_left_cen_t;
		std::vector<double> v_right_cen_x;
		std::vector<double> v_right_cen_t;
		std::vector<double> v_both_cen_x;
		std::vector<double> v_both_cen_t;
		std::vector<double> v_bothL_cen_x;
		std::vector<double> v_left_mid_x;
		std::vector<double> v_left_mid_t;
		std::vector<double> v_right_mid_x;
		std::vector<double> v_right_mid_t;
		std::vector<double> v_both_mid_x;
		std::vector<double> v_both_mid_t;
		std::vector<double> v_bothL_mid_x;
		std::vector<double> v_left_end_x;
		std::vector<double> v_left_end_t;
		std::vector<double> v_right_end_x;
		std::vector<double> v_right_end_t;
		std::vector<double> v_both_end_x;
		std::vector<double> v_both_end_t;
		std::vector<double> v_bothL_end_x;

		// for comparing left/right/both-side
		std::vector<double> v_LmB_func_dx;
		std::vector<double> v_LmB_func_t;
		std::vector<double> v_RmB_func_dx;
		std::vector<double> v_RmB_func_t;

		double m_RLmB_dx_max;

		// for comparing this/last xt
		std::vector<double> v_TmL_left_dx;
		std::vector<double> v_TmL_left_t;
		std::vector<double> v_TmL_right_dx;
		std::vector<double> v_TmL_right_t;
		std::vector<double> v_TmL_both_dx;
		std::vector<double> v_TmL_both_t;

		double m_TmL_LR_max;
		double m_TmL_B_max;

		// for comparing sample points with functions
		std::vector<double> v_SmF_left_t;
		std::vector<double> v_SmF_left_dx;
		std::vector<double> v_SmF_right_t;
		std::vector<double> v_SmF_right_dx;
		std::vector<double> v_SmF_both_t;
		std::vector<double> v_SmF_both_dx;

		// Graphs for XT fitting
		TGraph * gr_right_cen;
		TGraph * gr_left_cen;
		TGraph * gr_both_cen; // for neutral driftD
		TGraph * gr_bothL_cen; // for neutral driftD
		TGraph * gr_right_mid;
		TGraph * gr_left_mid;
		TGraph * gr_both_mid; // for neutral driftD
		TGraph * gr_bothL_mid; // for neutral driftD
		TGraph * gr_right_end;
		TGraph * gr_left_end;
		TGraph * gr_both_end; // for neutral driftD
		TGraph * gr_bothL_end; // for neutral driftD

		// for comparing left/right/both-side
		TGraph * gr_RmB_func;
		TGraph * gr_LmB_func;

		// for comparing this/last xt
		TGraph * gr_TmL_left;
		TGraph * gr_TmL_right;
		TGraph * gr_TmL_both;

		// for comparing sample points with functions
		TGraph * gr_SmF_left;
		TGraph * gr_SmF_right;
		TGraph * gr_SmF_both;

		// Functions for XT fitting
		TF1 * f_right_cen;
		TF1 * f_left_cen;
		TF1 * f_both_cen; // for neutral driftD
		TF1 * f_right_mid;
		TF1 * f_left_mid;
		TF1 * f_both_mid; // for neutral driftD
		TF1 * f_right_end;
		TF1 * f_left_end;
		TF1 * f_both_end; // for neutral driftD

		TF1 * f_left_deltac;
		TF1 * f_right_deltac;
		TF1 * f_both_deltac;
		TF1 * f_left_delta;
		TF1 * f_right_delta;
		TF1 * f_both_delta;

		TF1 * f_left_com;
		TF1 * f_right_com;
		TF1 * f_both_com;
		TF1 * f_bothL_com;

		TF1 * fo_right; // previous XT curve
		TF1 * fo_left; // previous XT curve
		TF1 * fo_both; // previous XT curve
		TF1 * f_right;
		TF1 * f_left;
		TF1 * f_right0;
		TF1 * f_left0;
		TF1 * f_rightEO;
		TF1 * f_leftEO;
};

#endif
