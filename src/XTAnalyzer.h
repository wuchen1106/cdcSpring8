#ifndef XTANALYZER_H
#define XTANALYZER_H

#include <vector>

#include "TString.h"

#define NSLICEX 161 // 161 bins from -8.05 mm to 8.05 mm, binning every 100 um
#define NSLICET 528 // 264 bins from -25 ns to 800 ns including left and right side, bining every 3.125 ns (3 TDC)

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
		XTAnalyzer(int debug = 0);
		virtual ~XTAnalyzer(void);

		void SetXTType(int type);
		void SetSaveHists(int save);

		int  Initialize(TString runname, int lid, TFile * infile, TFile * outfile, TTree * otree, int xttype, bool savehists,bool saveXT0 = false);
		void Process(void);

		void Push(double t, double x);

		int t2d(double t, double & d, bool isRight);

	private:

		int t2i(double t, bool positive);
		int x2i(double x);
		void i2x(int i,double & fdmin, double & fdmid, double & fdmax);
		void i2t(int i,double & dtmin, double & dtmid, double & dtmax);
		void getOneThirdLines(TH1D* h, double & left, double & right);
		double findFirstZero(TF1 * f, double xmin, double xmax, double delta);
		void getT8(double & t8left, double & t8right, double & t8both);
		TF1 * myNewTF1(TString name, TString form, double left, double right);
		TGraph * myNewTGraph(TString name, int n, const double * x, const double * y, TString title, TString Xtitle, TString Ytitle, int mType, double mSize, int mColor, double lSize, int lColor);
		TF1 * scalePolN(TF1 * f, double factor, TString name = "f_test");
		TF1 * minusPolN(TString name, TF1 * f1, TF1 * f2, double xmin, double xmax);
		TF1 * combinePolN(TString name, TF1 * f1, TF1 * f2, double x0, double x1, double x2, double xmin, double xmax);
		TString formPolN(TF1 * f);
		void createGraphs();
		void drawFitting(TH1D* h,TF1 * f,TCanvas * c,TString title, TString filename,double left, double right);
		void drawSamplingsLR();
		void drawSamplingsB();
		void drawLRB();
		void drawIteration();
		void writeObjects();

	private:
		// options
		int mDebugLevel;
		bool mSaveHists;
		bool mSaveXT0;
		int mXTType;
		TString mRunName;

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

		double mXmax;
		double mNbinx;

		// Histograms for data points
		TH2D * h2_xt;
		TH2D * h2_xtn; // for neutral driftD

		// Histograms for bin analysis
		TH1D * h_t[NSLICEX];
		TH1D * h_x[NSLICET];
		TH1D * h_tn[NSLICEX]; // for neutral driftD
		TH1D * h_xn[NSLICET]; // for neutral driftD

		// functions to fit x/t slices
		TF1 * f_x;
		TF1 * f_t;

		// about drawing the histograms
		TLine * l_left;
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

		// vectors selected for XT graphs
		std::vector<double> v_left_mid_x;
		std::vector<double> v_leftR_mid_x;
		std::vector<double> v_left_mid_t;
		std::vector<double> v_right_mid_x;
		std::vector<double> v_right_mid_t;
		std::vector<double> v_left_end_x;
		std::vector<double> v_leftR_end_x;
		std::vector<double> v_left_end_t;
		std::vector<double> v_right_end_x;
		std::vector<double> v_right_end_t;
		std::vector<double> v_both_mid_x;
		std::vector<double> v_both_mid_t;
		std::vector<double> v_both_end_x;
		std::vector<double> v_both_end_t;

		// for comparing left/right/both-side
		std::vector<double> v_LmB_mid_dt;
		std::vector<double> v_LmB_mid_x;
		std::vector<double> v_RmB_mid_dt;
		std::vector<double> v_RmB_mid_x;
		std::vector<double> v_LmB_end_dx;
		std::vector<double> v_LmB_end_t;
		std::vector<double> v_RmB_end_dx;
		std::vector<double> v_RmB_end_t;

		std::vector<double> v_LmB_func_mid_dt;
		std::vector<double> v_LmB_func_mid_x;
		std::vector<double> v_RmB_func_mid_dt;
		std::vector<double> v_RmB_func_mid_x;
		std::vector<double> v_LmB_func_end_dx;
		std::vector<double> v_LmB_func_end_t;
		std::vector<double> v_RmB_func_end_dx;
		std::vector<double> v_RmB_func_end_t;

		double m_RLmB_dx_max;
		double m_RLmB_dt_max;

		// for comparing this/last xt
		std::vector<double> v_TmL_left_dx;
		std::vector<double> v_TmL_left_t;
		std::vector<double> v_TmL_right_dx;
		std::vector<double> v_TmL_right_t;
		std::vector<double> v_TmL_both_dx;
		std::vector<double> v_TmL_both_t;

		double m_TmL_LR_max;
		double m_TmL_B_max;

		// Graphs for XT fitting
		TGraph * gr_right_mid;
		TGraph * gr_left_mid;
		TGraph * gr_leftR_mid;
		TGraph * gr_right_end;
		TGraph * gr_left_end;
		TGraph * gr_leftR_end;
		TGraph * gr_both_mid; // for neutral driftD
		TGraph * gr_both_end; // for neutral driftD

		// for comparing left/right/both-side
		TGraph * gr_RmB_mid;
		TGraph * gr_LmB_mid;
		TGraph * gr_RmB_end;
		TGraph * gr_LmB_end;

		TGraph * gr_RmB_func_mid;
		TGraph * gr_LmB_func_mid;
		TGraph * gr_RmB_func_end;
		TGraph * gr_LmB_func_end;

		// for comparing this/last xt
		TGraph * gr_TmL_left;
		TGraph * gr_TmL_right;
		TGraph * gr_TmL_both;

		// Functions for XT fitting
		TF1 * f_right_mid;
		TF1 * f_left_mid;
		TF1 * f_right_end;
		TF1 * f_left_end;
		TF1 * f_both_mid; // for neutral driftD
		TF1 * f_both_end; // for neutral driftD

		TF1 * f_left_delta;
		TF1 * f_right_delta;
		TF1 * f_both_delta;
		TF1 * f_left_com;
		TF1 * f_leftR_com;
		TF1 * f_right_com;
		TF1 * f_both_com;

		TF1 * f_RmB;
		TF1 * f_LmB;

		TF1 * fo_right; // previous XT curve
		TF1 * fo_left; // previous XT curve
		TF1 * fo_both; // previous XT curve
		TF1 * f_right;
		TF1 * f_left;
		TF1 * f_right0;
		TF1 * f_left0;
};

#endif
