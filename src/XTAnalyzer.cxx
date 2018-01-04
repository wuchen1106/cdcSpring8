#include "XTAnalyzer.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TPad.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TStyle.h"

XTAnalyzer::XTAnalyzer(int debug)
	:mDebugLevel(debug)
{
	v_x_slicex.resize(NSLICEX);
	v_t_slicex.resize(NSLICEX);
	v_sig_slicex.resize(NSLICEX);
	v_n_slicex.resize(NSLICEX);
	v_chi2_slicex.resize(NSLICEX);
	v_x_slicet.resize(NSLICET);
	v_t_slicet.resize(NSLICET);
	v_sig_slicet.resize(NSLICET);
	v_n_slicet.resize(NSLICET);
	v_chi2_slicet.resize(NSLICET);
	v_x_slicexn.resize(NSLICEX);
	v_t_slicexn.resize(NSLICEX);
	v_sig_slicexn.resize(NSLICEX);
	v_n_slicexn.resize(NSLICEX);
	v_chi2_slicexn.resize(NSLICEX);
	v_x_slicetn.resize(NSLICET);
	v_t_slicetn.resize(NSLICET);
	v_sig_slicetn.resize(NSLICET);
	v_n_slicetn.resize(NSLICET);
	v_chi2_slicetn.resize(NSLICET);
}

XTAnalyzer::~XTAnalyzer(void){
}

void XTAnalyzer::SetXTType(int type){
	mXTType = type;
}

void XTAnalyzer::SetSaveHists(int save){
	mSaveHists = save;
}

int XTAnalyzer::Initialize(TString runname, int lid, TFile * infile, TFile * outfile, TTree * otree, int xttype, int savehists, bool saveXT0){
	// Set options
	mRunName = runname;
	mLayerID = lid;
	mInFile = infile;
	mOutFile = outfile;
	mOutTree = otree;
	mXTType = xttype;
	mSaveHists = savehists;
	mSaveXT0 = saveXT0;

	mEntriesMin = 30;
	mSigTmax = 20;
	mSigXmax = 0.5;

	// load previous xt curves
	if (!mInFile) {
		fprintf(stderr,"WARNING: input XT file is not valid\n");
		return 1;
	}
	fo_left = (TF1*)mInFile->Get(Form("flc_%d",mLayerID));
	fo_right = (TF1*)mInFile->Get(Form("frc_%d",mLayerID));
	fo_both = (TF1*)mInFile->Get(Form("fbc_%d",mLayerID));
	if (!fo_left||!fo_right||!fo_both){
        fo_left = (TF1*)mInFile->Get(Form("fl_%d",0));
        fo_right = (TF1*)mInFile->Get(Form("fr_%d",0));
        fo_both = (TF1*)mInFile->Get(Form("fr_%d",0));
    }
	if (!fo_left||!fo_right||!fo_both){
		fprintf(stderr,"WARNING: cannot find flc_%d/frc_%d/fbc_%d in input XT file\n",mLayerID,mLayerID,mLayerID);
		return 2;
	}
	fo_left->SetName(Form("fl_old_%d",mLayerID));
	fo_right->SetName(Form("fr_old_%d",mLayerID));
	fo_both->SetName(Form("fb_old_%d",mLayerID));

	// set branches for the output tree.
	if (!mOutTree){
		fprintf(stderr,"WARNING: output XT tree is not initialized!\n");
		return 3;
	}
	mOutTree->SetBranchAddress("x",&mX);
	mOutTree->SetBranchAddress("t",&mT);
	mOutTree->SetBranchAddress("lid",&mLayerID);
	mOutTree->SetBranchAddress("sig",&mSig);
	mOutTree->SetBranchAddress("chi2",&mChi2);
	mOutTree->SetBranchAddress("n",&mEntries);
	mOutTree->SetBranchAddress("type",&mType);

	// set for slicing
//	mBWX = 0.08;
	mBWX = 0.04;
	mXLEFT = -mBWX*NSLICEX/2;
	mXRIGHT = mBWX*NSLICEX/2;
	mBWT = 3/0.96;
//	mBWT = 5/0.96;
//	mBWT = 7/0.96;
//	mBWT = 9/0.96;
	mTLEFT = -mBWT*NSLICET/2;
	mTRIGHT = mBWT*NSLICET/2;
	// set for binning
	mTmin = -25-1/0.96/2; // t range for one x bin
	mTmax = 800+1/0.96/2;
	mNbint = 792+1;
	mXmin = -0.02;
	mXmax = 10.02; // x range for one t bin
	mNbinx = 501;

	// prepare 2D histograms
	h2_xt = new TH2D(Form("h2_xt_%d",mLayerID),"XT Relation",mNbint,mTmin,mTmax,mNbinx,-mXmax,mXmax);
	h2_xt->GetXaxis()->SetTitle("T [ns]");
	h2_xt->GetYaxis()->SetTitle("X [mm]");
	h2_xtn = new TH2D(Form("h2_xtn_%d",mLayerID),"XT Relation",mNbint,mTmin,mTmax,mNbinx/2,-mXmin,mXmax);
	h2_xtn->GetXaxis()->SetTitle("T [ns]");
	h2_xtn->GetYaxis()->SetTitle("X [mm]");

	// prepare histograms for slice analysis
	for (int ix = 0; ix<NSLICEX; ix++){
        double fdmin,fdmid,fdmax;
        i2x(ix,fdmin,fdmid,fdmax);
		h_t[ix] = new TH1D(Form("h_t_%d_%d",mLayerID,ix),Form("T_{drift} distribution with DOCA [%.2f~%.2f] mm",fdmin,fdmax),mNbint,mTmin,mTmax);
		h_t_xsum[ix] = new TH1D(Form("h_t_xsum_%d_%d",mLayerID,ix),Form("Average DOCA in slice [%.2f~%.2f] mm",fdmin,fdmax),mNbint,mTmin,mTmax);
		h_tn[ix] = new TH1D(Form("h_tn_%d_%d",mLayerID,ix),Form("T_{drift} distribution with DOCA [%.2f~%.2f] mm",fdmin,fdmax),mNbint,mTmin,mTmax);
		h_tn_xsum[ix] = new TH1D(Form("h_tn_xsum_%d_%d",mLayerID,ix),Form("Average DOCA in slice [%.2f~%.2f] mm",fdmin,fdmax),mNbint,mTmin,mTmax);
	}
	for (int it = 0; it<NSLICET; it++){
        double dtmin,dtmid,dtmax;
        i2t(it,dtmin,dtmid,dtmax);
		h_x[it] = new TH1D(Form("h_x_%d_%d",mLayerID,it),Form("DOCA (left) distribution with T_{drift} [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
		h_mx[it] = new TH1D(Form("h_mx_%d_%d",mLayerID,it),Form("-DOCA (left) distribution with T_{drift} [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
		h_x_tsum[it] = new TH1D(Form("h_x_tsum_%d_%d",mLayerID,it),Form("Average T_{drift} in slice [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
		h_xn[it] = new TH1D(Form("h_xn_%d_%d",mLayerID,it),Form("DOCA (left) distribution with T_{drift} [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
		h_mxn[it] = new TH1D(Form("h_mxn_%d_%d",mLayerID,it),Form("-DOCA (left) distribution with T_{drift} [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
		h_xn_tsum[it] = new TH1D(Form("h_xn_tsum_%d_%d",mLayerID,it),Form("Average T_{drift} in slice [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
	}

	// prepare functions for slice analysis
	f_gaus = myNewTF1("fgaus","gaus",-mTmax,mTmax);
	f_land = myNewTF1("fland","landau",-mTmax,mTmax);
	f_2gaus = myNewTF1("f2gaus","(x>[1])*[0]*exp(-0.5*((x-[1])/[2])**2)+(x<=[1])*[0]*exp(-0.5*((x-[1])/[3])**2)",-mTmax,mTmax);
	l_left = new TLine(0,0,1,1);
	l_right = new TLine(0,0,1,1);
	l_center = new TLine(0,0,1,1);
	l_center->SetLineColor(kRed);
	l_center->SetLineStyle(2);

	// clear vectors for chosen XT sample points
	v_left_cen_x.clear();
	v_left_cen_t.clear();
	v_right_cen_x.clear();
	v_right_cen_t.clear();
	v_both_cen_x.clear();
	v_both_cen_t.clear();
	v_bothL_cen_x.clear();
	v_left_mid_x.clear();
	v_left_mid_t.clear();
	v_right_mid_x.clear();
	v_right_mid_t.clear();
	v_both_mid_x.clear();
	v_both_mid_t.clear();
	v_bothL_mid_x.clear();
	v_left_end_x.clear();
	v_left_end_t.clear();
	v_right_end_x.clear();
	v_right_end_t.clear();
	v_both_end_x.clear();
	v_both_end_t.clear();
	v_bothL_end_x.clear();
	v_LmB_func_dx.clear();
	v_LmB_func_t.clear();
	v_RmB_func_dx.clear();
	v_RmB_func_t.clear();
	v_TmL_left_dx.clear();
	v_TmL_left_t.clear();
	v_TmL_right_dx.clear();
	v_TmL_right_t.clear();
	v_TmL_both_dx.clear();
	v_TmL_both_t.clear();

	v_t_slicetls.clear();
	v_t_slicetrs.clear();
	v_t_slicetns.clear();
	v_sig_slicetls.clear();
	v_sig_slicetrs.clear();
	v_sig_slicetns.clear();

	v_SmF_left_t.clear();
	v_SmF_left_dx.clear();
	v_SmF_right_t.clear();
	v_SmF_right_dx.clear();
	v_SmF_both_t.clear();
	v_SmF_both_dx.clear();

	m_RLmB_dx_max = 0;

	m_TmL_LR_max = 0;
	m_TmL_B_max = 0;

	// prepare new XT functions
	f_left_cen = myNewTF1(Form("flce_%d",mLayerID),"pol9",mTmin,mTmax);
	f_right_cen = myNewTF1(Form("frce_%d",mLayerID),"pol9",mTmin,mTmax);
	f_both_cen = myNewTF1(Form("fbce_%d",mLayerID),"pol9",mTmin,mTmax);
	f_left_mid = myNewTF1(Form("flm_%d",mLayerID),"pol9",mTmin,mTmax);
	f_right_mid = myNewTF1(Form("frm_%d",mLayerID),"pol9",mTmin,mTmax);
	f_both_mid = myNewTF1(Form("fbm_%d",mLayerID),"pol9",mTmin,mTmax);
	f_left_end = myNewTF1(Form("fle_%d",mLayerID),"pol5",mTmin,mTmax);
	f_right_end = myNewTF1(Form("fre_%d",mLayerID),"pol5",mTmin,mTmax);
	f_both_end = myNewTF1(Form("fbe_%d",mLayerID),"pol5",mTmin,mTmax);

	if (mDebugLevel>=1) {printf("XTAnalyzer successfully initialized!\n");fflush(stdout);}
	return 0;
}

void XTAnalyzer::Push(double t, double x){
	if (mDebugLevel>=10) printf("XTAnalyzer::push(%.2e,%.2e)\n",t,x);
	double absx = fabs(x);
	h2_xt->Fill(t,x);
	h2_xtn->Fill(t,absx);
	int ix = x2i(x);
	for (int i = (ix>0?ix-1:0); i<=(ix+1>=NSLICEX?NSLICEX-1:ix+1); i++){
		h_t[i]->Fill(t);
		int ibin = h_t_xsum[i]->FindBin(t);
		h_t_xsum[i]->AddBinContent(ibin,x);
	}
	int it = t2i(t,x>0);
	for (int i = (it>0?it-1:0); i<=(it+1>=NSLICET?NSLICET-1:it+1); i++){
		h_x[i]->Fill(x);
		h_mx[i]->Fill(-absx);
		int ibin = h_x_tsum[i]->FindBin(x);
		h_x_tsum[i]->AddBinContent(ibin,t);
	}
	ix = x2i(absx);
	for (int i = (ix>NSLICEX/2?ix-1:NSLICEX/2); i<=(ix+1>=NSLICEX?NSLICEX-1:ix+1); i++){
		h_tn[i]->Fill(t);
		int ibin = h_tn_xsum[i]->FindBin(t);
		h_tn_xsum[i]->AddBinContent(ibin,absx);
	}
	it = t2i(t,true);
	for (int i = (it>NSLICET/2?it-1:NSLICET/2); i<=(it+1>=NSLICET?NSLICET-1:it+1); i++){
		h_xn[i]->Fill(absx);
		h_mxn[i]->Fill(-absx);
		int ibin = h_xn_tsum[i]->FindBin(absx);
		h_xn_tsum[i]->AddBinContent(ibin,t);
	}
	if (mDebugLevel>=10) printf("            pushed\n",t,x);
}

void XTAnalyzer::Process(void){
	if (mDebugLevel>0) {printf("In XTAnalyzer::Process\n");fflush(stdout);}
	//==========================Taking Samples from Slices==============================
	// fit x histograms, and push to vectors & tree
    TCanvas * canv_fitting = new TCanvas("cfit","cfit",1024,768);
	int minEntries = 100;
	for (int i = 0; i<NSLICET; i++){
		if (mDebugLevel>0) {printf("=>h_x[%d\]\n",i);fflush(stdout);}
		mType = 0;
		double left, right;
		double divleft,divright;
		i2t(i,divleft,mT,divright);
		mEntries = h_x[i]->Integral();
		fitSliceHistFloat(h_x[i],0.5,mX,mSig,mChi2,left,right);
		if (mDebugLevel>0) {printf("h_x[%d] (%d) after fitSliceHistFloat: x=%.2f, sig=%.2f, chi2=%.2f, left = %.2f, right = %.2f\n",i,(int)mEntries,mX,mSig,mChi2,left,right);fflush(stdout);}
		if (mEntries>minEntries){
			TF1 * f = 0;
			bool flipped = false;
			if (fabs(mX)>7&&fabs(mX)<7.8){ // FIXME: boundary up to tuning
				if (mX>0){
					flipped = true;
					fitSliceHistFloat(h_mx[i],0.3,mX,mSig,mChi2,left,right);
					f = fitSliceLand(h_mx[i],mX,mSig,mChi2,left,right);
					double temp = left;
					left = -right;
					right = -temp;
					mX*=-1;
				}
				else{
					fitSliceHistFloat(h_x[i],0.3,mX,mSig,mChi2,left,right);
					f = fitSliceLand(h_x[i],mX,mSig,mChi2,left,right);
				}
			}
			else{
				f = fitSliceGaus(h_x[i],mX,mSig,mChi2,left,right);
			}
			if (mDebugLevel>0) {printf("h_x[%d] after fitSlice: x=%.2f, sig=%.2f, chi2=%.2f, left = %.2f, right = %.2f\n",i,mX,mSig,mChi2,left,right);fflush(stdout);}
			mT = getMean(h_x[i],h_x_tsum[i],left,right);
			if (flipped){
				h_mx[i]->GetXaxis()->SetRangeUser(-mX-mSigXmax*3,-mX+mSigXmax*3);
				if (mSaveHists>=2) drawFitting(h_mx[i],f,canv_fitting,Form("%.1f-%.1f-%.1f ns, N=%d, x=%.2f mm, #sigma=%.0f um, #chi^{2}=%.1f",divleft,mT,divright,mEntries,mX,mSig*1000,mChi2),Form("h_x%d_%s.png",i,mRunName.Data()),-right,-mX,-left);
			}
			else{
				h_x[i]->GetXaxis()->SetRangeUser(mX-mSigXmax*3,mX+mSigXmax*3);
				if (mSaveHists>=2) drawFitting(h_x[i],f,canv_fitting,Form("%.1f-%.1f-%.1f ns, N=%d, x=%.2f mm, #sigma=%.0f um, #chi^{2}=%.1f",divleft,mT,divright,mEntries,mX,mSig*1000,mChi2),Form("h_x%d_%s.png",i,mRunName.Data()),left,mX,right);
			}
		}
		else mT = getMean(h_x[i],h_x_tsum[i],left,right);
		v_n_slicet[i] = mEntries;
		v_sig_slicet[i] = mSig;
		v_chi2_slicet[i] = mChi2;
		v_x_slicet[i] = mX;
		v_t_slicet[i] = mT;
		mOutTree->Fill();
	}
	// fit t histograms, and push to vectors & tree
	for (int i = 0; i<NSLICEX; i++){
		mType = 1;
		double left, right;
		double divleft,divright;
		i2x(i,divleft,mX,divright);
		mEntries = h_t[i]->Integral();
		fitSliceHistFloat(h_t[i],0.5,mT,mSig,mChi2,left,right);
		if (mEntries>minEntries){
			TF1 * f = 0;
			if (fabs(mX)<1){
				fitSliceHistFloat(h_t[i],0.3,mT,mSig,mChi2,left,right);
				f = fitSliceLand(h_t[i],mT,mSig,mChi2,left,right);
			}
			else{
				f = fitSliceGaus(h_t[i],mT,mSig,mChi2,left,right);
			}
			mX = getMean(h_t[i],h_t_xsum[i],left,right);
			h_t[i]->GetXaxis()->SetRangeUser(mT-mSigTmax*3,mT+mSigTmax*3);
			if (mSaveHists>=2) drawFitting(h_t[i],f,canv_fitting,Form("%.2f-%.2f-%.2f mm, N=%d, t=%.1f ns, #sigma=%.1f ns, #chi^{2}=%.1f",divleft,mX,divright,mEntries,mT,mSig,mChi2),Form("h_t%d_%s.png",i,mRunName.Data()),left,mT,right);
		}
		else mX = getMean(h_t[i],h_t_xsum[i],left,right);
		v_n_slicex[i] = mEntries;
		v_sig_slicex[i] = mSig;
		v_chi2_slicex[i] = mChi2;
		v_x_slicex[i] = mX;
		v_t_slicex[i] = mT;
		mOutTree->Fill();
	}
	// fit x histograms for both-side case, and push to vectors & tree
	for (int i = NSLICET/2; i<NSLICET; i++){
		if (mDebugLevel>0) {printf("=>h_xn[%d\]\n",i);fflush(stdout);}
		mType = 2;
		double left, right;
		double divleft,divright;
		i2t(i,divleft,mT,divright);
		mEntries = h_xn[i]->Integral();
		fitSliceHistFloat(h_xn[i],0.5,mX,mSig,mChi2,left,right);
		if (mDebugLevel>0) {printf("h_xn[%d] (%d) after fitSliceHistFloat: x=%.2f, sig=%.2f, chi2=%.2f, left = %.2f, right = %.2f\n",i,(int)mEntries,mX,mSig,mChi2,left,right);fflush(stdout);}
		if (mEntries>minEntries){
			TF1 * f = 0;
			if (fabs(mX)>7.15&&fabs(mX)<7.8){ // FIXME: boundary up to tuning
				fitSliceHistFloat(h_mxn[i],0.3,mX,mSig,mChi2,left,right);
				f = fitSliceLand(h_mxn[i],mX,mSig,mChi2,left,right);
				mX*=-1;
				double temp = left;
				left = -right;
				right = -temp;
				mT = getMean(h_xn[i],h_xn_tsum[i],left,right);
				if (mDebugLevel>0) {printf("h_xn[%d] after fitSlice: x=%.2f, sig=%.2f, chi2=%.2f, left = %.2f, right = %.2f\n",i,mX,mSig,mChi2,left,right);fflush(stdout);}
				h_mxn[i]->GetXaxis()->SetRangeUser(-mX-mSigXmax*3,-mX+mSigXmax*3);
				if (mSaveHists>=1) drawFitting(h_mxn[i],f,canv_fitting,Form("%.1f-%.1f-%.1f ns, N=%d, x=%.2f mm, #sigma=%.0f um, #chi^{2}=%.1f",divleft,mT,divright,mEntries,mX,mSig*1000,mChi2),Form("h_xn%d_%s.png",i,mRunName.Data()),-right,-mX,-left);
			}
			else{
				f = fitSliceGaus(h_xn[i],mX,mSig,mChi2,left,right);
				mT = getMean(h_xn[i],h_xn_tsum[i],left,right);
				if (mDebugLevel>0) {printf("h_xn[%d] after fitSlice: x=%.2f, sig=%.2f, chi2=%.2f, left = %.2f, right = %.2f\n",i,mX,mSig,mChi2,left,right);fflush(stdout);}
				h_xn[i]->GetXaxis()->SetRangeUser(mX-mSigXmax*3,mX+mSigXmax*3);
				if (mSaveHists>=1) drawFitting(h_xn[i],f,canv_fitting,Form("%.1f-%.1f-%.1f ns, N=%d, x=%.2f mm, #sigma=%.0f um, #chi^{2}=%.1f",divleft,mT,divright,mEntries,mX,mSig*1000,mChi2),Form("h_xn%d_%s.png",i,mRunName.Data()),left,mX,right);
			}
		}
		else mT = getMean(h_xn[i],h_xn_tsum[i],left,right);
		v_n_slicetn[i] = mEntries;
		v_sig_slicetn[i] = mSig;
		v_chi2_slicetn[i] = mChi2;
		v_x_slicetn[i] = mX;
		v_t_slicetn[i] = mT;
		mOutTree->Fill();
	}
	// fit t histograms for both-side case, and push to vectors & tree
	for (int i = NSLICEX/2; i<NSLICEX; i++){
		mType = 3;
		double left, right;
		double divleft,divright;
		i2x(i,divleft,mX,divright);
		mEntries = h_tn[i]->Integral();
		fitSliceHistFloat(h_tn[i],0.5,mT,mSig,mChi2,left,right);
		if (mEntries>minEntries){
			TF1 * f = 0;
			if (fabs(mX)<1){
				fitSliceHistFloat(h_tn[i],0.3,mT,mSig,mChi2,left,right);
				f = fitSliceLand(h_tn[i],mT,mSig,mChi2,left,right);
			}
			else{
				f = fitSliceGaus(h_tn[i],mT,mSig,mChi2,left,right);
			}
			mX = getMean(h_tn[i],h_tn_xsum[i],left,right);
			h_tn[i]->GetXaxis()->SetRangeUser(mT-mSigTmax*3,mT+mSigTmax*3);
			if (mSaveHists>=1) drawFitting(h_tn[i],f,canv_fitting,Form("%.2f-%.2f-%.2f mm, N=%d, t=%.1f ns, #sigma=%.1f ns, #chi^{2}=%.1f",divleft,mX,divright,mEntries,mT,mSig,mChi2),Form("h_tn%d_%s.png",i,mRunName.Data()),left,mT,right);
		}
		else mX = getMean(h_tn[i],h_tn_xsum[i],left,right);
		v_n_slicexn[i] = mEntries;
		v_sig_slicexn[i] = mSig;
		v_chi2_slicexn[i] = mChi2;
		v_x_slicexn[i] = mX;
		v_t_slicexn[i] = mT;
		mOutTree->Fill();
	}

	//==========================Select Samples==============================
	// select sample points and make graphs
	// FIXME: Currently seperate the mid/end graphs by 8 mm line, and search for the real one from 7 mm line.
	// FIXME: Currently seperate the cen/mid graphs by 1 mm line, and search for the real one from mTmin.
	double xCenter2Mid = 1.5;
	double xStart2Turn = 7;
	double xMargin = 0.2;
	double t8Left = v_t_slicex[0];
	double t8Right = v_t_slicex[NSLICEX-1];
	double t8Both = v_t_slicexn[NSLICEX-1];
	int iTLeft, iTRight, iTBoth;
	getTT(iTLeft,v_x_slicet,v_t_slicet,v_sig_slicet,v_n_slicet,true);
	getTT(iTRight,v_x_slicet,v_t_slicet,v_sig_slicet,v_n_slicet,false);
	getTT(iTBoth,v_x_slicetn,v_t_slicetn,v_sig_slicetn,v_n_slicetn,false);
	double tMargin = 10;
	double t7Left = v_t_slicex[1./mBWX];
	double t7Right = v_t_slicex[NSLICEX-1./mBWX];
	double t7Both = v_t_slicexn[NSLICEX-1./mBWX];
	if (mDebugLevel>=1){
		printf("Before selecting samples:\n");
		printf(" tTl: %.1f @%d, tTr: %.1f @%d, tTb: %.1f @%d \n",v_t_slicet[iTLeft],iTLeft,v_t_slicet[iTRight],iTRight,v_t_slicetn[iTBoth],iTBoth);
		printf(" t7l:%.1f, t7r:%.1f, t7b:%.1f\n",t7Left,t7Right,t7Both);
	}
	sigmaXReset();
	for (int i = 0; i<NSLICET/2; i++){ // x samples in t slices, left
		if (mDebugLevel>=2) printf("  LR T slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.2f\n",i,v_x_slicet[i],v_t_slicet[i],v_n_slicet[i],v_sig_slicet[i]);
		if (v_n_slicet[i]<mEntriesMin||v_sig_slicet[i]>mSigXmax||v_sig_slicet[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (i<=iTLeft+1){ // left end
			if (mDebugLevel>=2) printf("                  t>%.1f, push to left_end!\n",t8Left);
			v_left_end_x.push_back(v_x_slicet[i]);
			v_left_end_t.push_back(v_t_slicet[i]);
		}
		if (i>=iTLeft-1&&v_x_slicet[i]<=-xStart2Turn){ // turning part
			if (mDebugLevel>=2) printf("                  t<%.1f x<%.2f, push to left_mid!\n",t8Left+tMargin,-xStart2Turn);
			v_left_mid_x.push_back(v_x_slicet[i]);
			v_left_mid_t.push_back(v_t_slicet[i]);
		}
		sigmaXIncrement(v_t_slicet[i],v_sig_slicet[i],v_n_slicet[i],v_t_slicetls,v_sig_slicetls);
	}
	sigmaXFinalcheck(v_t_slicetls,v_sig_slicetls);
	// t samples in x slices, find out the minimum t first
	double mint_t_slicex_l = mTmax;
	double mint_x_slicex_l = 0;
	double mint_t_slicex_r = mTmax;
	double mint_x_slicex_r = 0;
	for (int i = 0; i<NSLICEX; i++){
		double t = v_t_slicex[i];
		double x = v_x_slicex[i];
		if (i<=NSLICEX/2){ // left
			if (mint_t_slicex_l>t){
				mint_t_slicex_l = t;
				mint_x_slicex_l = x;
			}
		}
		if (i>=NSLICEX/2){ // right
			if (mint_t_slicex_r>t){
				mint_t_slicex_r = t;
				mint_x_slicex_r = x;
			}
		}
	}
	for (int i = 0; i<NSLICEX; i++){ // t samples in x slices
		double t = v_t_slicex[i];
		double x = v_x_slicex[i];
		if (mDebugLevel>=2) printf("  LR X slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.1f\n",i,x,t,v_n_slicex[i],v_sig_slicex[i]);
		if (v_n_slicex[i]<mEntriesMin||v_sig_slicex[i]>mSigTmax||v_sig_slicex[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (i<=NSLICEX/2){ // left
			if (x>mint_x_slicex_l) t = mint_t_slicex_l;
			if (x>-xStart2Turn){ // middle part
				if (x>-xCenter2Mid){
					if (mDebugLevel>=2) printf("                  %.2f<x, push to left_cen!\n",-xCenter2Mid);
					v_left_cen_x.push_back(x);
					v_left_cen_t.push_back(t);
				}
				if (x<-xCenter2Mid+xMargin){
					if (mDebugLevel>=2) printf("                  %.2f<x<%.2f, push to left_mid!\n",-xStart2Turn,-xCenter2Mid+xMargin);
					v_left_mid_x.push_back(x);
					v_left_mid_t.push_back(t);
				}
			}
		}
		if (i>=NSLICEX/2){ // right
			if (x<mint_x_slicex_r) t = mint_t_slicex_r;
			if (x<xStart2Turn){ // middle part
				if (x<xCenter2Mid){
					if (mDebugLevel>=2) printf("                  x<%.2f, push to right_cen!\n",xCenter2Mid);
					v_right_cen_x.push_back(x);
					v_right_cen_t.push_back(t);
				}
				if (x>xCenter2Mid-xMargin){
					if (mDebugLevel>=2) printf("                  %.2f<x<%.2f, push to right_mid!\n",xCenter2Mid-xMargin,xStart2Turn);
					v_right_mid_x.push_back(x);
					v_right_mid_t.push_back(t);
				}
			}
		}
	}
	sigmaXReset();
	for (int i = NSLICET/2; i<NSLICET; i++){ // x samples in t slices, right 
		if (mDebugLevel>=2) printf("  LR T slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.2f\n",i,v_x_slicet[i],v_t_slicet[i],v_n_slicet[i],v_sig_slicet[i]);
		if (v_n_slicet[i]<mEntriesMin||v_sig_slicet[i]>mSigXmax||v_sig_slicet[i]<=0) continue;
		if (i>=iTRight-1){ // right end
			if (mDebugLevel>=2) printf("                  t>%.1f, push to right_end!\n",t8Right);
			v_right_end_x.push_back(v_x_slicet[i]);
			v_right_end_t.push_back(v_t_slicet[i]);
		}
		if (i<=iTRight+1&&v_x_slicet[i]>=xStart2Turn){ // turning part
			if (mDebugLevel>=2) printf("                  t<%.1f x>%.2f, push to right_mid!\n",t8Right+tMargin,xStart2Turn);
			v_right_mid_x.push_back(v_x_slicet[i]);
			v_right_mid_t.push_back(v_t_slicet[i]);
		}
		sigmaXIncrement(v_t_slicet[i],v_sig_slicet[i],v_n_slicet[i],v_t_slicetrs,v_sig_slicetrs);
	}
	sigmaXFinalcheck(v_t_slicetrs,v_sig_slicetrs);
	// t samples in x slices, both, find out the minimum t first
	double mint_t_slicex_b = mTmax;
	double mint_x_slicex_b = 0;
	for (int i = NSLICEX/2; i<NSLICEX; i++){
		double t = v_t_slicexn[i];
		double x = v_x_slicexn[i];
		if (mint_t_slicex_b>t){
			mint_t_slicex_b = t;
			mint_x_slicex_b = x;
		}
	}
	for (int i = NSLICEX/2; i<NSLICEX; i++){ // t samples in x slices, both-side
		double t = v_t_slicexn[i];
		double x = v_x_slicexn[i];
		if (x<mint_x_slicex_b) t = mint_t_slicex_b;
		if (mDebugLevel>=2) printf("  BS X slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.1f\n",i,x,t,v_n_slicexn[i],v_sig_slicexn[i]);
		if (v_n_slicexn[i]<mEntriesMin||v_sig_slicexn[i]>mSigTmax||v_sig_slicexn[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (x<xStart2Turn){ // middle part
			if (x<xCenter2Mid){
				if (mDebugLevel>=2) printf("                  %.2f<x, push to both_cen!\n",xCenter2Mid);
				v_both_cen_x.push_back(x);
				v_bothL_cen_x.push_back(-x);
				v_both_cen_t.push_back(t);
			}
			if (x>xCenter2Mid-xMargin){
				if (mDebugLevel>=2) printf("                  %.2f<x<%.2f, push to both_mid!\n",xCenter2Mid-xMargin,xStart2Turn);
				v_both_mid_x.push_back(x);
				v_bothL_mid_x.push_back(-x);
				v_both_mid_t.push_back(t);
			}
		}
	}
	sigmaXReset();
	for (int i = NSLICET/2; i<NSLICET; i++){ // x samples in t slices, both-side
		if (mDebugLevel>=2) printf("  BS T slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.2f\n",i,v_x_slicetn[i],v_t_slicetn[i],v_n_slicetn[i],v_sig_slicetn[i]);
		if (v_n_slicetn[i]<mEntriesMin||v_sig_slicetn[i]>mSigXmax||v_sig_slicetn[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (i>=iTBoth-1){ // both-side end
			if (mDebugLevel>=2) printf("                  t>%.1f, push to both_end!\n",t8Both);
			v_both_end_x.push_back(v_x_slicetn[i]);
			v_bothL_end_x.push_back(-v_x_slicetn[i]);
			v_both_end_t.push_back(v_t_slicetn[i]);
		}
		if (i<=iTBoth+1&&v_x_slicetn[i]>=xStart2Turn){ // turning part
			if (mDebugLevel>=2) printf("                  t<%.1f, x>=%.2f, push to both_mid!\n",t8Both+tMargin,xStart2Turn);
			v_both_mid_x.push_back(v_x_slicetn[i]);
			v_bothL_mid_x.push_back(-v_x_slicetn[i]);
			v_both_mid_t.push_back(v_t_slicetn[i]);
		}
		sigmaXIncrement(v_t_slicetn[i],v_sig_slicetn[i],v_n_slicetn[i],v_t_slicetns,v_sig_slicetns);
	}
	sigmaXFinalcheck(v_t_slicetns,v_sig_slicetns);

	//==========================Prepare graphs==============================
	createGraphs();

	//==========================Get XT Functions==============================
	// fit xt functions
	// FIXME: should think about fixing some parameters to let left/right sides XTs go through the same 0 point
	//        now we have poor shape of XT near 0 point, and we'd better rely on both-side XT
	if (!gr_left_cen->GetN()) fprintf(stderr,"WARNING: gr_left_cen is empty!\n"); else gr_left_cen->Fit(Form("flce_%d",mLayerID),"qN0","");;
	if (!gr_right_cen->GetN()) fprintf(stderr,"WARNING: gr_right_cen is empty!\n"); else gr_right_cen->Fit(Form("frce_%d",mLayerID),"qN0","");;
	if (!gr_both_cen->GetN()) fprintf(stderr,"WARNING: gr_both_cen is empty!\n"); else gr_both_cen->Fit(Form("fbce_%d",mLayerID),"qN0","");;
	if (!gr_left_mid->GetN()) fprintf(stderr,"WARNING: gr_left_mid is empty!\n"); else gr_left_mid->Fit(Form("flm_%d",mLayerID),"qN0","");;
	if (!gr_right_mid->GetN()) fprintf(stderr,"WARNING: gr_right_mid is empty!\n"); else gr_right_mid->Fit(Form("frm_%d",mLayerID),"qN0","");;
	if (!gr_both_mid->GetN()) fprintf(stderr,"WARNING: gr_both_mid is empty!\n"); else gr_both_mid->Fit(Form("fbm_%d",mLayerID),"qN0","");;
	if (!gr_left_end->GetN()) fprintf(stderr,"WARNING: gr_left_end is empty!\n"); else gr_left_end->Fit(Form("fle_%d",mLayerID),"qN0","");;
	if (!gr_right_end->GetN()) fprintf(stderr,"WARNING: gr_right_end is empty!\n"); else gr_right_end->Fit(Form("fre_%d",mLayerID),"qN0","");;
	if (!gr_both_end->GetN()) fprintf(stderr,"WARNING: gr_both_end is empty!\n"); else gr_both_end->Fit(Form("fbe_%d",mLayerID),"qN0","");;

	// get ranges for XT functions
	f_left_deltac = minusPolN(Form("fldc_%d",mLayerID),f_left_cen,f_left_mid,mTmin,mTmax);
	f_right_deltac = minusPolN(Form("frdc_%d",mLayerID),f_right_cen,f_right_mid,mTmin,mTmax);
	f_both_deltac = minusPolN(Form("fbdc_%d",mLayerID),f_both_cen,f_both_mid,mTmin,mTmax);
	f_left_delta = minusPolN(Form("fld_%d",mLayerID),f_left_mid,f_left_end,0,mTmax);
	f_right_delta = minusPolN(Form("frd_%d",mLayerID),f_right_mid,f_right_end,0,mTmax);
	f_both_delta = minusPolN(Form("fbd_%d",mLayerID),f_both_mid,f_both_end,0,mTmax);
	double tZeroLeft = findFirstZero(f_left_cen,mTmin,mTmax,1);
	double tZeroRight = findFirstZero(f_right_cen,mTmin,mTmax,1);
	double tZeroBoth = findFirstZero(f_both_cen,mTmin,mTmax,1);
	double tCentLeft = findFirstZero(f_left_deltac,5,mTmax,1); // FIXME: better to use driftT at 0.5 mm to start searching
	double tCentRight = findFirstZero(f_right_deltac,5,mTmax,1);
	double tCentBoth = findFirstZero(f_both_deltac,5,mTmax,1);
	double tTurnLeft = findFirstZero(f_left_delta,t7Left,mTmax,1);
	double tTurnRight = findFirstZero(f_right_delta,t7Right,mTmax,1);
	double tTurnBoth = findFirstZero(f_both_delta,t7Both,mTmax,1);
	double tEndLeft = v_left_end_t.size()>0?v_left_end_t[0]:0;
	double tEndRight = v_right_end_t.size()>0?v_right_end_t[v_right_end_t.size()-1]:0;
	double tEndBoth = v_both_end_t.size()>0?v_both_end_t[v_both_end_t.size()-1]:0;
	if (mDebugLevel>=1){
		printf("After selecting samples:\n");
		printf(" t0l:%.1f, t0r:%.1f, t0b:%.1f\n",tZeroLeft,tZeroRight,tZeroBoth);
		printf(" tcl:%.1f, tcr:%.1f, tcb:%.1f\n",tCentLeft,tCentRight,tCentBoth);
		printf(" ttl:%.1f, ttr:%.1f, ttb:%.1f\n",tTurnLeft,tTurnRight,tTurnBoth);
		printf(" tel:%.1f, ter:%.1f, teb:%.1f\n",tEndLeft,tEndRight,tEndBoth);
	}

	// combine functions
	f_left_com = combinePolN(Form("flc_%d",mLayerID),f_left_cen,f_left_mid,f_left_end,tZeroLeft,tCentLeft,tTurnLeft,tEndLeft,mTmin,tEndLeft);
	f_right_com = combinePolN(Form("frc_%d",mLayerID),f_right_cen,f_right_mid,f_right_end,tZeroRight,tCentRight,tTurnRight,tEndRight,mTmin,tEndRight);
	f_both_com = combinePolN(Form("fbc_%d",mLayerID),f_both_cen,f_both_mid,f_both_end,tZeroBoth,tCentBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
	f_bothL_com = combinePolN(Form("fblc_%d",mLayerID),scalePolN(f_both_cen,-1),scalePolN(f_both_mid,-1),scalePolN(f_both_end,-1),tZeroBoth,tCentBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
	f_left_com->SetLineColor(kMagenta);
	f_both_com->SetLineColor(kBlack);
	f_bothL_com->SetLineColor(kBlack);

	// get the final xt functions
	if (mXTType==0){ // use Left/Right case
		f_left = combinePolN(Form("fl_%d",mLayerID),f_left_cen,f_left_mid,f_left_end,tZeroLeft,tCentLeft,tTurnLeft,tEndLeft,mTmin,tEndLeft);
		f_right = combinePolN(Form("fr_%d",mLayerID),f_right_cen,f_right_mid,f_right_end,tZeroRight,tCentRight,tTurnRight,tEndRight,mTmin,tEndRight);
		if (mSaveXT0){
			f_left0 = combinePolN(Form("fl_%d",0),f_left_cen,f_left_mid,f_left_end,tZeroLeft,tCentLeft,tTurnLeft,tEndLeft,mTmin,tEndLeft);
			f_right0 = combinePolN(Form("fr_%d",0),f_right_cen,f_right_mid,f_right_end,tZeroRight,tCentRight,tTurnRight,tEndRight,mTmin,tEndRight);
		}
	}
	else{ // use Both-Side case
		f_left = combinePolN(Form("fl_%d",mLayerID),scalePolN(f_both_cen,-1),scalePolN(f_both_mid,-1),scalePolN(f_both_end,-1),tZeroBoth,tCentBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
		f_right = combinePolN(Form("fr_%d",mLayerID),f_both_cen,f_both_mid,f_both_end,tZeroBoth,tCentBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
		if (mSaveXT0){
			f_left0 = combinePolN(Form("fl_%d",0),scalePolN(f_both_cen,-1),scalePolN(f_both_mid,-1),scalePolN(f_both_end,-1),tZeroBoth,tCentBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
			f_right0 = combinePolN(Form("fr_%d",0),f_both_cen,f_both_mid,f_both_end,tZeroBoth,tCentBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
		}
	}

	//==========================Prepare more graphs==============================
	// get Left/Right/Both-sides differences by functions
	for (int i = 0; i<1024; i++){
		double t = mTmax*(i+0.5)/1024;
		double xl = f_left_com->Eval(t);
		double xr = f_right_com->Eval(t);
		double xb = f_both_com->Eval(t);
		if (t<=tEndLeft&&t<=tEndBoth){
			v_LmB_func_t.push_back(t);
			v_LmB_func_dx.push_back((-xl-xb)*1000);
		}
		if (t<=tEndRight&&t<=tEndBoth){
			v_RmB_func_t.push_back(t);
			v_RmB_func_dx.push_back((xr-xb)*1000);
		}
	}
	// get this/last xt difference by functions
	double tEndLeftPre = fo_left->GetXmax();
	double tEndRightPre = fo_right->GetXmax();
	double tEndBothPre = fo_both->GetXmax();
	for (int i = 0; i<1024; i++){
		double t = mTmax*(i+0.5)/1024;
		if (t<=tEndLeft&&t<=tEndLeftPre){
			double x = f_left_com->Eval(t);
			double xo = fo_left->Eval(t);
			double dx = (x-xo)*1000;
			v_TmL_left_t.push_back(t);
			v_TmL_left_dx.push_back(dx);
			if (fabs(dx)>m_TmL_LR_max) m_TmL_LR_max = fabs(dx);
		}
		if (t<=tEndRight&&t<=tEndRightPre){
			double x = f_right_com->Eval(t);
			double xo = fo_right->Eval(t);
			double dx = (x-xo)*1000;
			v_TmL_right_t.push_back(t);
			v_TmL_right_dx.push_back(dx);
			if (fabs(dx)>m_TmL_LR_max) m_TmL_LR_max = fabs(dx);
		}
		if (t<=tEndBoth&&t<=tEndBothPre){
			double x = f_both_com->Eval(t);
			double xo = fo_both->Eval(t);
			double dx = (x-xo)*1000;
			v_TmL_both_t.push_back(t);
			v_TmL_both_dx.push_back(dx);
			if (fabs(dx)>m_TmL_B_max) m_TmL_B_max = fabs(dx);
		}
	}
	// get Sample/Function differences
	for (int i = 0; i<v_left_end_t.size(); i++){
		double t = v_left_end_t[i];
		double x = v_left_end_x[i];
		double xf = f_left_com->Eval(t);
		v_SmF_left_t.push_back(t);
		v_SmF_left_dx.push_back((x-xf)*1000);
	}
	for (int i = 3; i<v_left_mid_t.size(); i++){ // 3 points overlaping with end
		double t = v_left_mid_t[i];
		double x = v_left_mid_x[i];
		if (x>-xCenter2Mid) continue; // marginal region overlaping with center
		double xf = f_left_com->Eval(t);
		v_SmF_left_t.push_back(t);
		v_SmF_left_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_left_cen_t.size(); i++){
		double t = v_left_cen_t[i];
		double x = v_left_cen_x[i];
		double xf = f_left_com->Eval(t);
		v_SmF_left_t.push_back(t);
		v_SmF_left_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_right_cen_t.size(); i++){
		double t = v_right_cen_t[i];
		double x = v_right_cen_x[i];
		double xf = f_right_com->Eval(t);
		v_SmF_right_t.push_back(t);
		v_SmF_right_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_right_mid_t.size()-3; i++){ // 3 points overlaping with end
		double t = v_right_mid_t[i];
		double x = v_right_mid_x[i];
		if (x<xCenter2Mid) continue; // marginal region overlaping with center
		double xf = f_right_com->Eval(t);
		v_SmF_right_t.push_back(t);
		v_SmF_right_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_right_end_t.size(); i++){
		double t = v_right_end_t[i];
		double x = v_right_end_x[i];
		double xf = f_right_com->Eval(t);
		v_SmF_right_t.push_back(t);
		v_SmF_right_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_both_cen_t.size(); i++){
		double t = v_both_cen_t[i];
		double x = v_both_cen_x[i];
		double xf = f_both_com->Eval(t);
		v_SmF_both_t.push_back(t);
		v_SmF_both_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_both_mid_t.size()-3; i++){ // 3 points overlaping with end
		double t = v_both_mid_t[i];
		double x = v_both_mid_x[i];
		if (x<xCenter2Mid) continue; // marginal region overlaping with center
		double xf = f_both_com->Eval(t);
		v_SmF_both_t.push_back(t);
		v_SmF_both_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_both_end_t.size(); i++){
		double t = v_both_end_t[i];
		double x = v_both_end_x[i];
		double xf = f_both_com->Eval(t);
		v_SmF_both_t.push_back(t);
		v_SmF_both_dx.push_back((x-xf)*1000);
	}

	gr_LmB_func= myNewTGraph(Form("gr_LmB_func_%d",mLayerID),v_LmB_func_t.size(),&(v_LmB_func_t[0]),&(v_LmB_func_dx[0]),
			"XT Differences with Both-side Combined Case","T [ns]","#Delta_{X} [um]",20,0.5,kMagenta,0.5,kMagenta);
	gr_RmB_func= myNewTGraph(Form("gr_RmB_func_%d",mLayerID),v_RmB_func_t.size(),&(v_RmB_func_t[0]),&(v_RmB_func_dx[0]),
			"XT Differences with Both-side Combined Case","T [ns]","#Delta_{X} [um]",20,0.5,kRed,0.5,kRed);
	gr_TmL_left = myNewTGraph(Form("gr_TmL_left_%d",mLayerID),v_TmL_left_t.size(),&(v_TmL_left_t[0]),&(v_TmL_left_dx[0]),
			"XT Differences Comparing with Last Time","T [ns]","#Delta_X [um]",20,0.5,kMagenta,0.5,kMagenta);
	gr_TmL_right = myNewTGraph(Form("gr_TmL_right_%d",mLayerID),v_TmL_right_t.size(),&(v_TmL_right_t[0]),&(v_TmL_right_dx[0]),
			"XT Differences Comparing with Last Time","T [ns]","#Delta_X [um]",20,0.5,kRed,0.5,kRed);
	gr_TmL_both = myNewTGraph(Form("gr_TmL_both_%d",mLayerID),v_TmL_both_t.size(),&(v_TmL_both_t[0]),&(v_TmL_both_dx[0]),
			"XT Differences Comparing with Last Time","T [ns]","#Delta_X [um]",20,0.5,kBlack,0.5,kBlack);
	gr_SmF_left = myNewTGraph(Form("gr_SmF_left_%d",mLayerID),v_SmF_left_t.size(),&(v_SmF_left_t[0]),&(v_SmF_left_dx[0]),
			"XT Sampling Points Subtracted by fitted function","T [ns]","#Delta_X [um]",20,0.5,kMagenta,0.5,kMagenta);
	gr_SmF_right = myNewTGraph(Form("gr_SmF_right_%d",mLayerID),v_SmF_right_t.size(),&(v_SmF_right_t[0]),&(v_SmF_right_dx[0]),
			"XT Sampling Points Subtracted by fitted function","T [ns]","#Delta_X [um]",20,0.5,kRed,0.5,kRed);
	gr_SmF_both = myNewTGraph(Form("gr_SmF_both_%d",mLayerID),v_SmF_both_t.size(),&(v_SmF_both_t[0]),&(v_SmF_both_dx[0]),
			"XT Sampling Points Subtracted by fitted function","T [ns]","#Delta_X [um]",20,0.5,kBlack,0.5,kBlack);

	//==========================Draw Plots==============================
	// draw plots to compare left/right/both-side
	drawLRB();
	// Draw the XT histogram and plots
	drawSamplingsLR();
	// Draw the XT histogram and plots for both-side case
	drawSamplingsB();
	// draw pltos to check iterations
	drawIteration();

	//==========================Save==============================
	writeObjects();
}

int XTAnalyzer::t2d(double t, double & d, bool isRight){
	int status = 0;
	TF1 * f = 0;
	if (isRight) f = f_right;
	else f = f_left;
	double tmax = f->GetXmax();
	double tmin = f->GetXmin();
	if (t<tmin){
		d = 0;
		status = -1;
	}
	else if (t>tmax){
		d = f->Eval(tmax);
		status = 1;
	}
	else {
		d = f->Eval(t);
		status = 0;
	}
	return status;
}

void XTAnalyzer::i2x(int i, double & fdmin, double & fdmid, double & fdmax){
    fdmin = i*mBWX+mXLEFT;
    fdmid = fdmin+mBWX/2.;
    fdmax = fdmid+mBWX/2.;
    fdmin-=mBWX; // take adjacent two slices into account
    fdmax+=mBWX;
}

void XTAnalyzer::i2t(int i, double & dtmin, double & dtmid, double & dtmax){
    dtmin = i*mBWT+mTLEFT;
    if (dtmin<-mBWT/2){ // 
    	dtmin = fabs((i+1)*mBWT+mTLEFT);
    }
	else{
		dtmin = dtmin;
	}
    dtmid = dtmin+mBWT/2.;
    dtmax = dtmid+mBWT/2.;
    dtmin-=mBWT; // take adjacent two slices into account
    dtmax+=mBWT;
}

int XTAnalyzer::x2i(double fitD){
    int i=(fitD-mXLEFT)/mBWX;
    if (i>=NSLICEX) i=-999;
    return i;
}

int XTAnalyzer::t2i(double driftT,bool positive){
	if (!positive) driftT = -driftT;
    int i=(driftT-mTLEFT)/mBWT;
    if (fabs(driftT)<mBWT/2) i = NSLICET/2;
    if (i>=NSLICET||i<0) i=-999;
    return i;
}

double XTAnalyzer::getMean(TH1D * h, TH1D * hy, double left, double right){
	int bleft = h->FindBin(left);
	int bright = h->FindBin(right);
	double ysum = hy->Integral(bleft,bright);
	int entries = h->Integral(bleft,bright);
	double y = ysum;
	if (entries) y/=entries;
	return y;
}

void XTAnalyzer::fitSliceHistFloat(TH1D * h, double ratio, double & mean, double & sigma, double & chi2, double & left, double & right){
	int bmax = h->GetMaximumBin();
    double max = h->GetBinContent(bmax)*ratio;
    int binl = bmax-1;
    for (;binl>=3; binl--){
    	double height3bins = h->GetBinContent(binl);
    	height3bins+=h->GetBinContent(binl-1);
    	height3bins+=h->GetBinContent(binl-2);
        if (height3bins/3<max) break;
    }
    int binr = bmax+1;
    int nbins = h->GetNbinsX();
    for (;binr<=nbins-2; binr++){
    	double height3bins = h->GetBinContent(binr);
    	height3bins+=h->GetBinContent(binr+1);
    	height3bins+=h->GetBinContent(binr+2);
        if (height3bins/3<max) break;
    }
    left = h->GetBinCenter(binl);
    right = h->GetBinCenter(binr);
	h->GetXaxis()->SetRangeUser(left,right);
	mean = h->GetMean();
	sigma = h->GetRMS();
	chi2 = 0;
}

TF1 * XTAnalyzer::fitSliceGaus(TH1D * h, double & mean, double & sigma, double & chi2, double & left, double & right){
	if (mDebugLevel>0) printf("in fitSliceGaus: \"%s\" has %d(%d) entries\n",h->GetName(),h->GetEntries(),h->Integral());
	TF1 * f = f_gaus;
	int bmax = h->GetMaximumBin();
	mean = h->GetBinCenter(bmax);
    double lrange = mean-left;
    double rrange = right-mean;
	if (mDebugLevel>0) printf("  %.2f - %.2f - %.2f, lrange = %.2f, rrange = %.2f\n",left,mean,right,lrange,rrange);
	h->Fit("fgaus","qN0","",left,right);
	for (int i = 0; i<10; i++){
		mean = f->GetParameter(1);
		sigma = fabs(f->GetParameter(2));
		if (mDebugLevel>0) printf("      sigma: %.2f mean: %.2f\n",i,sigma,mean);
		if (sigma<rrange*5&&sigma>rrange/5&&mean>left&&mean<right) break;
		left-=lrange/10;
		right+=rrange/10;
		if (mDebugLevel>0) printf("   -> %d: %.2f -- %.2f\n",i,left,right);
		h->Fit("fgaus","qN0","",left,right);
		mean = f->GetParameter(1);
		sigma = fabs(f->GetParameter(2));
	}
	if (mean>right||mean<left) sigma = 1e9;
	chi2 = f->GetChisquare();
	return f;
}

TF1 * XTAnalyzer::fitSliceLand(TH1D * h, double & mean, double & sigma, double & chi2, double & left, double & right){
	if (mDebugLevel>0) printf("in fitSliceLand: \"%s\" has %d(%d) entries\n",h->GetName(),h->GetEntries(),h->Integral());
	TF1 * f = f_land;
	int bmax = h->GetMaximumBin();
	mean = h->GetBinCenter(bmax);
    double lrange = mean-left;
    double rrange = right-mean;
	if (mDebugLevel>0) printf("  %.2f - %.2f - %.2f, lrange = %.2f, rrange = %.2f\n",left,mean,right,lrange,rrange);
	h->Fit("fland","qN0","",left,right);
	for (int i = 0; i<10; i++){
		mean = f->GetParameter(1);
		sigma = fabs(f->GetParameter(2));
		if (mDebugLevel>0) printf("      sigma: %.2f mean: %.2f\n",i,sigma,mean);
		if (sigma<rrange*5&&sigma>rrange/5&&mean>left&&mean<right) break;
		left-=lrange/10;
		right+=rrange/10;
		if (mDebugLevel>0) printf("   -> %d: %.2f -- %.2f\n",i,left,right);
		h->Fit("fland","qN0","",left,right);
		mean = f->GetParameter(1);
		sigma = fabs(f->GetParameter(2));
	}
	if (mean>right||mean<left) sigma = 1e9;
	chi2 = f->GetChisquare();
	return f;
}

TF1 * XTAnalyzer::fitSlice2Gaus(TH1D * h, double & mean, double & sigma, double & chi2, double & left, double & right){
	TF1* f = f_2gaus;
	int bmax = h->GetMaximumBin();
	mean = h->GetBinCenter(bmax);
	double height = h->GetBinContent(bmax);
    double lrange = mean-left;
    double rrange = right-mean;
	sigma = h->GetRMS();
	f->SetParameters(height,mean,sigma,sigma);
	h->Fit("f2gaus","qN0","",left,right);
	double sigma1;
	double sigma2;
	for (int i = 0; i<10; i++){
		mean = f->GetParameter(1);
		sigma1 = fabs(f->GetParameter(2));
		sigma2 = fabs(f->GetParameter(3));
		if (sigma1<rrange*5&&sigma2<lrange*5&&sigma1>rrange/5&&sigma2>lrange/5&&mean>left&&mean<right) break;
		left-=lrange/10;
		right+=rrange/10;
		h->Fit("f2gaus","qN0","",left,right);
		mean = f->GetParameter(1);
		sigma1 = fabs(f->GetParameter(2));
		sigma2 = fabs(f->GetParameter(3));
	}
	// set results
	sigma = sqrt((sigma1*sigma1+sigma2*sigma2)/2);
	if (mean>right||mean<left) sigma = 1e9;
	chi2 = f->GetChisquare();
	return f;
}

double XTAnalyzer::findFirstZero(TF1 * f, double xmin, double xmax, double delta){
	double theX = 0;
	for (double x = xmin+delta; x<xmax; x+=delta){ // At least two solutions. Scan to find the smallest one
		theX = f->GetX(0,xmin,x);
		if (fabs(theX-x)>delta/10.&&fabs(theX-xmin)>delta/10.){
			break;
		}
	}
	if (mDebugLevel>=3) printf("findFirstZero: theX = %.2f from [%.2f,%.2f]\n",theX,xmin,xmax);
	return theX;
}

void XTAnalyzer::getTT(int & iT, std::vector<double> & vx, std::vector<double> & vt, std::vector<double> & vsig, std::vector<double> & vn, bool negtive){
	double d7 = 7;
	double t7 = 250;
	for (int i = NSLICET/2; i<NSLICET; i++){
		if (vsig[i]>mSigXmax||vn[i]<mEntriesMin) continue;
		if (vx[i]>7){
			d7 = vx[i];
			t7 = vt[i];
			break;
		}
	}
	double velocity = d7/t7;
	iT = NSLICET/2;
	int direction = 1; if (negtive) direction = -1;
	for (; iT>=0&&iT<NSLICET; iT+=direction){
		if (vsig[iT]>mSigXmax||vn[iT]<mEntriesMin) continue;
		if (fabs(vx[iT])>6) break;
	}
	for (; iT>=0&&iT<NSLICET; iT+=direction){
		if (vsig[iT]>mSigXmax||vn[iT]<mEntriesMin) continue;
		int j = iT+direction;
		for (;j>=0&&j<NSLICET; j+=direction){
			if (vsig[j]<=mSigTmax&&vn[j]>=mEntriesMin) break;
		}
		double vel = 0;
		if (vt[iT]-vt[j]) vel = (vx[iT]-vx[j])/(vt[iT]-vt[j]);
		else{
			fprintf(stderr,"WARNING: in XTAnalyzer::getTT, vt @ %d and %d are the same! Ignoring these two hits\n",iT,j);
			vel = velocity;
		}
		if (negtive) vel*=-1;
		if (vel<velocity/3) break;
	}
}

void XTAnalyzer::sigmaXReset(){
	m_sig_sel = 0;
	m_t_sel = 0;
	m_n_sel = 0;
	m_i_sel = 0;
}

void XTAnalyzer::sigmaXIncrement(double t, double sig, double n, std::vector<double> & vt, std::vector<double> & vs){
	m_sig_sel += sig*n;
	m_t_sel += t*n;
	m_n_sel += n;
	if (m_i_sel%4==3){
		if (m_n_sel){
			m_t_sel/=m_n_sel;
			m_sig_sel/=m_n_sel;
		}
		vt.push_back(m_t_sel);
		vs.push_back(m_sig_sel);
		m_t_sel = 0; m_sig_sel = 0; m_n_sel = 0;
	}
	m_i_sel++;
}

void XTAnalyzer::sigmaXFinalcheck(std::vector<double> & vt, std::vector<double> & vs){
	if (m_i_sel%4!=0){
		if (m_n_sel){
			m_t_sel/=m_n_sel;
			m_sig_sel/=m_n_sel;
		}
		vt.push_back(m_t_sel);
		vs.push_back(m_sig_sel);
	}
}

TF1 * XTAnalyzer::myNewTF1(TString name, TString form, double left, double right){
	TF1 * f = new TF1(name,form,left,right);
	f->SetNpx(1024);
	f->SetNumberFitPoints(1024);
	f->SetLineWidth(0.3);
	return f;
}

TGraph * XTAnalyzer::myNewTGraph(TString name, int n, const double * x, const double * y, TString title, TString Xtitle, TString Ytitle, int mType, double mSize, int mColor, double lSize, int lColor){
	TGraph * gr = new TGraph(n,x,y);
	gr->SetName(name);
	gr->SetTitle(title);
	gr->GetXaxis()->SetTitle(Xtitle);
	gr->GetYaxis()->SetTitle(Ytitle);
	gr->SetMarkerStyle(mType);
	gr->SetMarkerSize(mSize);
	gr->SetMarkerColor(mColor);
	gr->SetLineWidth(lSize);
	gr->SetLineColor(lColor);
	return gr;
}

TF1 * XTAnalyzer::scalePolN(TF1 * f, double factor, TString name){
	int n = f->GetNpar();
	double xmin = f->GetXmin();
	double xmax = f->GetXmax();
	TF1 * ff = new TF1(name,Form("pol%d",n),xmin,xmax);
	for (int i = 0; i<n; i++){
		double p=f->GetParameter(i);
		ff->SetParameter(i,p*factor);
	}
	return ff;
}

TF1 * XTAnalyzer::minusPolN(TString name, TF1 * f1, TF1 * f2, double xmin, double xmax){
	int n1 = f1->GetNpar();
	int n2 = f2->GetNpar();
	int n = n1>n2?n1:n2;
	TF1 * f = new TF1(name,Form("pol%d",n-1),xmin,xmax);
	for (int i = 0; i<n; i++){
		double p = 0;
		if (i<n1) p+=f1->GetParameter(i);
		if (i<n2) p-=f2->GetParameter(i);
		f->SetParameter(i,p);
	}
	f->SetNpx(1024);
	f->SetNumberFitPoints(1024);
	f->SetLineWidth(0.3);
	return f;
}

TF1 * XTAnalyzer::combinePolN(TString name, TF1 * f1, TF1 * f2, TF1 * f3, double x0, double x1, double x2, double x3, double xmin, double xmax){
	TString form1 = formPolN(f1);
	TString form2 = formPolN(f2);
	TString form3 = formPolN(f3);
	double ymax = f3->Eval(x3);
	TF1 * f = new TF1(name,Form("(x>=%.14f&&x<%.14f)*(%s)+(x>=%.14f&&x<%.14f)*(%s)+(x>=%.14f&&x<%.14f)*(%s)+(x>=%.14f)*%.14e",x0,x1,form1.Data(),x1,x2,form2.Data(),x2,x3,form3.Data(),x3,ymax),xmin,xmax);
	f->SetNpx(1024);
	f->SetNumberFitPoints(1024);
	f->SetLineWidth(0.3);
	return f;
}

TString XTAnalyzer::formPolN(TF1 * f){
	int n = f->GetNpar();
	double p0 = f->GetParameter(0);
	TString form = Form("%.14e",p0);
	TString x = "x";
	for (int i = 1; i<n; i++){
		double p = f->GetParameter(i);
		form+=Form("+(%.14e*%s)",p,x.Data());
		x+="*x";
	}
	return form;
}

void XTAnalyzer::createGraphs(){
	// L/R slices
	gr_xt_slicet = myNewTGraph(Form("gr_xt_slicet_%d",mLayerID),NSLICET,&(v_t_slicet[0]),&(v_x_slicet[0]),
			"XT Relation (X fitting in T slice)","X [mm]","T [ns]",
			4,0.3,kRed,0.1,kRed);
	gr_xt_slicex = myNewTGraph(Form("gr_xt_slicex_%d",mLayerID),NSLICEX,&(v_t_slicex[0]),&(v_x_slicex[0]),
			"XT Relation (T fitting in X slice)","X [mm]","T [ns]",
			4,0.3,kBlack,0.1,kBlack);
	gr_xn_slicex = myNewTGraph(Form("gr_xn_slicex_%d",mLayerID),NSLICEX,&(v_x_slicex[0]),&(v_n_slicex[0]),
			"Number of Entries in each X slice","X [mm]","Entries",
			20,0.5,kBlack,0.5,kBlack);
	gr_xsig_slicex = myNewTGraph(Form("gr_xsig_slicex_%d",mLayerID),NSLICEX,&(v_x_slicex[0]),&(v_sig_slicex[0]),
			"Sigma of T in each X slice","X [mm]","#sigma_{T} [ns]",
			20,0.5,kBlack,0.5,kBlack);
	gr_xchi2_slicex = myNewTGraph(Form("gr_xchi2_slicex_%d",mLayerID),NSLICEX,&(v_x_slicex[0]),&(v_chi2_slicex[0]),
			"#chi^{2} of T in each X slice","X [mm]","#chi^{2}_{T}",
			20,0.5,kBlack,0.5,kBlack);
	gr_nt_slicetl = myNewTGraph(Form("gr_nt_slicetl_%d",mLayerID),NSLICET/2,&(v_t_slicet[0]),&(v_n_slicet[0]),
			"Number of Entries in each T slice","T [ns]","Entries",
			20,0.5,kMagenta,0.5,kMagenta);
	gr_sigt_slicetl = myNewTGraph(Form("gr_sigt_slicetl_%d",mLayerID),NSLICET/2,&(v_t_slicet[0]),&(v_sig_slicet[0]),
			"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
			20,0.5,kMagenta,0.5,kMagenta);
	gr_chi2t_slicetl = myNewTGraph(Form("gr_chi2t_slicetl_%d",mLayerID),NSLICET/2,&(v_t_slicet[0]),&(v_chi2_slicet[0]),
			"#chi^{2} of X in each T slice","T [ns]","#chi^{2}_{X}",
			20,0.5,kMagenta,0.5,kMagenta);
	gr_nt_slicetr = myNewTGraph(Form("gr_nt_slicetr_%d",mLayerID),NSLICET/2,&(v_t_slicet[NSLICET/2]),&(v_n_slicet[NSLICET/2]),
			"Number of Entries in each X slice","T [ns]","Entries",
			20,0.5,kRed,0.5,kRed);
	gr_sigt_slicetr = myNewTGraph(Form("gr_sigt_slicetr_%d",mLayerID),NSLICET/2,&(v_t_slicet[NSLICET/2]),&(v_sig_slicet[NSLICET/2]),
			"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
			20,0.5,kRed,0.5,kRed);
	gr_chi2t_slicetr = myNewTGraph(Form("gr_chi2t_slicetr_%d",mLayerID),NSLICET/2,&(v_t_slicet[NSLICET/2]),&(v_chi2_slicet[NSLICET/2]),
			"#chi^{2} of X in each T slice","T [ns]","#chi^{2}_{X}",
			20,0.5,kRed,0.5,kRed);
	// Both-sides slices
	gr_xt_slicetn = myNewTGraph(Form("gr_xt_slicetn_%d",mLayerID),NSLICET/2,&(v_t_slicetn[NSLICET/2]),&(v_x_slicetn[NSLICET/2]),
			"XT Relation (X fitting in T slice)","X [mm]","T [ns]",
			4,0.3,kRed,0.1,kRed);
	gr_xt_slicexn = myNewTGraph(Form("gr_xt_slicexn_%d",mLayerID),NSLICEX/2+1,&(v_t_slicexn[NSLICEX/2]),&(v_x_slicexn[NSLICEX/2]), // NSLICEX should be an odd number!
			"XT Relation (T fitting in X slice)","X [mm]","T [ns]",
			4,0.3,kBlack,0.1,kBlack);
	gr_xn_slicexn = myNewTGraph(Form("gr_xn_slicexn_%d",mLayerID),NSLICEX/2+1,&(v_x_slicexn[NSLICEX/2]),&(v_n_slicexn[NSLICEX/2]),
			"Number of Entries in each X slice","X [mm]","Entries",
			20,0.5,kBlack,0.5,kBlack);
	gr_xsig_slicexn = myNewTGraph(Form("gr_xsig_slicexn_%d",mLayerID),NSLICEX/2+1,&(v_x_slicexn[NSLICEX/2]),&(v_sig_slicexn[NSLICEX/2]),
			"Sigma of T in each X slice","X [mm]","#sigma_{T} [ns]",
			20,0.5,kBlack,0.5,kBlack);
	gr_xchi2_slicexn = myNewTGraph(Form("gr_xchi2_slicexn_%d",mLayerID),NSLICEX/2+1,&(v_x_slicexn[NSLICEX/2]),&(v_chi2_slicexn[NSLICEX/2]),
			"#chi^{2} of T in each X slice","X [mm]","#chi^{2}_{T}",
			20,0.5,kBlack,0.5,kBlack);
	gr_nt_slicetn = myNewTGraph(Form("gr_nt_slicetn_%d",mLayerID),NSLICET/2,&(v_t_slicetn[NSLICET/2]),&(v_n_slicetn[NSLICET/2]),
			"Number of Entries in each T slice","T [ns]","Entries",
			20,0.5,kBlack,0.5,kBlack);
	gr_sigt_slicetn = myNewTGraph(Form("gr_sigt_slicetn_%d",mLayerID),NSLICET/2,&(v_t_slicetn[NSLICET/2]),&(v_sig_slicetn[NSLICET/2]),
			"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
			20,0.5,kBlack,0.5,kBlack);
	gr_chi2t_slicetn = myNewTGraph(Form("gr_chi2t_slicetn_%d",mLayerID),NSLICET/2,&(v_t_slicetn[NSLICET/2]),&(v_chi2_slicetn[NSLICET/2]),
			"#chi^{2} of X in each T slice","T [ns]","#chi^{2}_{X}",
			20,0.5,kBlack,0.5,kBlack);
	// selected graphs for XT fitting (and showing)
	gr_left_cen = myNewTGraph(Form("gr_xt_lce_%d",mLayerID),v_left_cen_x.size(),&(v_left_cen_t[0]),&(v_left_cen_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kMagenta,0.3,kMagenta);
	gr_right_cen = myNewTGraph(Form("gr_xt_rce_%d",mLayerID),v_right_cen_x.size(),&(v_right_cen_t[0]),&(v_right_cen_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kRed,0.3,kRed);
	gr_both_cen = myNewTGraph(Form("gr_xt_bce_%d",mLayerID),v_both_cen_x.size(),&(v_both_cen_t[0]),&(v_both_cen_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_bothL_cen = myNewTGraph(Form("gr_xt_blm_%d",mLayerID),v_bothL_cen_x.size(),&(v_both_cen_t[0]),&(v_bothL_cen_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_left_mid = myNewTGraph(Form("gr_xt_lm_%d",mLayerID),v_left_mid_x.size(),&(v_left_mid_t[0]),&(v_left_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kMagenta,0.3,kMagenta);
	gr_right_mid = myNewTGraph(Form("gr_xt_rm_%d",mLayerID),v_right_mid_x.size(),&(v_right_mid_t[0]),&(v_right_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kRed,0.3,kRed);
	gr_both_mid = myNewTGraph(Form("gr_xt_bm_%d",mLayerID),v_both_mid_x.size(),&(v_both_mid_t[0]),&(v_both_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_bothL_mid = myNewTGraph(Form("gr_xt_blm_%d",mLayerID),v_bothL_mid_x.size(),&(v_both_mid_t[0]),&(v_bothL_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_left_end = myNewTGraph(Form("gr_xt_le_%d",mLayerID),v_left_end_x.size(),&(v_left_end_t[0]),&(v_left_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kMagenta,0.3,kMagenta);
	gr_right_end = myNewTGraph(Form("gr_xt_re_%d",mLayerID),v_right_end_x.size(),&(v_right_end_t[0]),&(v_right_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kRed,0.3,kRed);
	gr_both_end = myNewTGraph(Form("gr_xt_be_%d",mLayerID),v_both_end_x.size(),&(v_both_end_t[0]),&(v_both_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_bothL_end = myNewTGraph(Form("gr_xt_ble_%d",mLayerID),v_bothL_end_x.size(),&(v_both_end_t[0]),&(v_bothL_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	// selected graphs for error reference
	if (mXTType==0){
		gr_sigts_slicetl = myNewTGraph(Form("gr_sigts_slicetl_%d",mLayerID),v_t_slicetls.size(),&(v_t_slicetls[0]),&(v_sig_slicetls[0]),
				"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
				20,0.5,kMagenta,0.5,kMagenta);
		gr_sigts_slicetr = myNewTGraph(Form("gr_sigts_slicetr_%d",mLayerID),v_t_slicetrs.size(),&(v_t_slicetrs[0]),&(v_sig_slicetrs[0]),
				"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
				20,0.5,kRed,0.5,kRed);
		if (mSaveXT0){
			gr_sigts_slicetl0 = myNewTGraph(Form("gr_sigts_slicetl_%d",0),v_t_slicetls.size(),&(v_t_slicetls[0]),&(v_sig_slicetls[0]),
					"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
					20,0.5,kMagenta,0.5,kMagenta);
			gr_sigts_slicetr0 = myNewTGraph(Form("gr_sigts_slicetr_%d",0),v_t_slicetrs.size(),&(v_t_slicetrs[0]),&(v_sig_slicetrs[0]),
					"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
					20,0.5,kRed,0.5,kRed);
		}
	}
	else{
		gr_sigts_slicetl = myNewTGraph(Form("gr_sigts_slicetl_%d",mLayerID),v_t_slicetns.size(),&(v_t_slicetns[0]),&(v_sig_slicetns[0]),
				"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
				20,0.5,kMagenta,0.5,kMagenta);
		gr_sigts_slicetr = myNewTGraph(Form("gr_sigts_slicetr_%d",mLayerID),v_t_slicetns.size(),&(v_t_slicetns[0]),&(v_sig_slicetns[0]),
				"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
				20,0.5,kRed,0.5,kRed);
		if (mSaveXT0){
			gr_sigts_slicetl0 = myNewTGraph(Form("gr_sigts_slicetl_%d",0),v_t_slicetns.size(),&(v_t_slicetns[0]),&(v_sig_slicetns[0]),
					"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
					20,0.5,kMagenta,0.5,kMagenta);
			gr_sigts_slicetr0 = myNewTGraph(Form("gr_sigts_slicetr_%d",0),v_t_slicetns.size(),&(v_t_slicetns[0]),&(v_sig_slicetns[0]),
					"Sigma of X in each T slice","T [ns]","#sigma_{X} [mm]",
					20,0.5,kRed,0.5,kRed);
		}
	}
}

void XTAnalyzer::drawFitting(TH1D* h,TF1 * f, TCanvas * c,TString title, TString filename,double left, double center, double right){
	if (!h) fprintf(stderr,"ERROR: in drawFitting, input histogram does not exist!\n");
	if (!f) fprintf(stderr,"ERROR: in drawFitting, input function does not exist!\n");
	if (!c) fprintf(stderr,"ERROR: in drawFitting, input canvas does not exist!\n");
	if (!h||!f||!c) return;
	if (mDebugLevel>1) printf("drawFitting %s",h->GetName());
	c->cd();
	h->SetTitle(title);
	h->Draw();
	double max = h->GetMaximum();
	l_left->SetX1(left);
	l_left->SetX2(left);
	l_left->SetY1(0);
	l_left->SetY2(max);
	l_center->SetX1(center);
	l_center->SetX2(center);
	l_center->SetY1(0);
	l_center->SetY2(max);
	l_right->SetX1(right);
	l_right->SetX2(right);
	l_right->SetY1(0);
	l_right->SetY2(max);
	l_left->Draw();
	l_center->Draw();
	l_right->Draw();
	f->Draw("SAME");
	c->SaveAs(filename);
}

void XTAnalyzer::drawSamplingsLR(){
	TCanvas * canv_xtsamples = new TCanvas("canv_xtsamples","canv_xtsamples",800,1000);
	TPad * pad_xtsamples[2];
	for (int il = 0; il<1; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*1+il;
			pad_xtsamples[index] = new TPad(Form("pad_xtsamples%d",index),"pad",il/1.,(2-ir)/2.,(il+1)/1.,(1-ir)/2.);
			pad_xtsamples[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_xtsamples[index]->SetGridx(1);
			pad_xtsamples[index]->SetGridy(1);
		}
	}
	pad_xtsamples[0]->cd();
	h2_xt->Draw("COLZ");
	gr_xt_slicet->Draw("PSAME");
	gr_xt_slicex->Draw("PSAME");
	gr_left_cen->Draw("PSAME");
	gr_right_cen->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
	gr_right_mid->Draw("PSAME");
	gr_left_end->Draw("PSAME");
	gr_right_end->Draw("PSAME");
	f_left_com->Draw("SAME");
	f_right_com->Draw("SAME");
	pad_xtsamples[1]->cd();
	TH2D * h2_samplings = new TH2D("h2_samplings",gr_SmF_left->GetTitle(),mNbint,mTmin,mTmax,512,-200,200);
	h2_samplings->Draw();
	gr_SmF_left->Draw("PLSAME");
	gr_SmF_right->Draw("PLSAME");
	canv_xtsamples->SaveAs("xtsamples_"+mRunName+".png");
	canv_xtsamples->SaveAs("xtsamples_"+mRunName+".pdf");
	// Add zoom-in plots
	TCanvas * canv_xtsamplesz = new TCanvas("canv_xtsamplesz","canv_xtsamplesz",800,1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	gr_xt_slicet->SetMarkerSize(1);
	gr_xt_slicex->SetMarkerSize(1);
	gr_left_end->SetMarkerSize(1);
	gr_left_cen->SetMarkerSize(1);
	gr_right_cen->SetMarkerSize(1);
	gr_left_mid->SetMarkerSize(1);
	gr_right_mid->SetMarkerSize(1);
	gr_right_end->SetMarkerSize(1);
	gr_xt_slicet->Draw("PSAME");
	gr_xt_slicex->Draw("PSAME");
	gr_left_cen->Draw("PSAME");
	gr_right_cen->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
	gr_right_mid->Draw("PSAME");
	gr_right_end->Draw("PSAME");
	h2_xt->GetXaxis()->SetRangeUser(-10,60);
	h2_xt->GetYaxis()->SetRangeUser(-2,2);
	h2_xt->Draw("COLZ");
	gr_xt_slicet->Draw("PSAME");
	gr_xt_slicex->Draw("PSAME");
	gr_left_cen->Draw("PSAME");
	gr_right_cen->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
	gr_right_mid->Draw("PSAME");
	f_left_com->Draw("SAME");
	f_right_com->Draw("SAME");
	canv_xtsamplesz->SaveAs("xtsamples_center_"+mRunName+".png");
	canv_xtsamplesz->SaveAs("xtsamples_center_"+mRunName+".pdf");
	h2_xt->GetXaxis()->SetRangeUser(250,450);
	h2_xt->GetYaxis()->SetRangeUser(6.5,8);
	h2_xt->Draw("COLZ");
	gr_xt_slicet->Draw("PSAME");
	gr_xt_slicex->Draw("PSAME");
	gr_right_mid->Draw("PSAME");
	gr_right_end->Draw("PSAME");
	f_right_com->Draw("SAME");
	canv_xtsamplesz->SaveAs("xtsamples_endR_"+mRunName+".png");
	canv_xtsamplesz->SaveAs("xtsamples_endR_"+mRunName+".pdf");
	h2_xt->GetXaxis()->SetRangeUser(250,450);
	h2_xt->GetYaxis()->SetRangeUser(-8,-6.5);
	h2_xt->Draw("COLZ");
	gr_xt_slicet->Draw("PSAME");
	gr_xt_slicex->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
	gr_left_end->Draw("PSAME");
	f_left_com->Draw("SAME");
	canv_xtsamplesz->SaveAs("xtsamples_endL_"+mRunName+".png");
	canv_xtsamplesz->SaveAs("xtsamples_endL_"+mRunName+".pdf");
	h2_xt->GetXaxis()->UnZoom();
	h2_xt->GetYaxis()->UnZoom();
	// fitting status of slices 
	TCanvas * canv_xtslices = new TCanvas("canv_xtslices","canv_xtslices",1364,768);
	TPad * pad_xtslices[6];
	for (int il = 0; il<3; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*3+il;
			pad_xtslices[index] = new TPad(Form("pad_xtslices%d",index),"pad",il/3.,(2-ir)/2.,(il+1)/3.,(1-ir)/2.);
			pad_xtslices[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_xtslices[index]->SetGridx(1);
			pad_xtslices[index]->SetGridy(1);
		}
	}
	pad_xtslices[0]->cd();
	gr_xn_slicex->Draw("AP");
	pad_xtslices[1]->cd();
	gr_xsig_slicex->GetYaxis()->SetRangeUser(0,mSigTmax);
	gr_xsig_slicex->Draw("AP");
	pad_xtslices[2]->cd();
	gr_xchi2_slicex->Draw("AP");
	pad_xtslices[3]->cd();
	pad_xtslices[3]->SetLogy(1);
	gr_nt_slicetl->Draw("AP");
	gr_nt_slicetr->Draw("PSAME");
	pad_xtslices[4]->cd();
	gr_sigt_slicetl->GetYaxis()->SetRangeUser(0,mSigXmax);
	gr_sigt_slicetl->Draw("AP");
	gr_sigt_slicetr->Draw("PSAME");
	pad_xtslices[5]->cd();
	gr_chi2t_slicetl->Draw("AP");
	gr_chi2t_slicetr->Draw("PSAME");
	canv_xtslices->SaveAs("xtslices_"+mRunName+".png");
	canv_xtslices->SaveAs("xtslices_"+mRunName+".pdf");
}

void XTAnalyzer::drawSamplingsB(){
	TCanvas * canv_xtsamplesn = new TCanvas("canv_xtsamplesn","canv_xtsamplesn",800,1000);
	TPad * pad_xtsamplesn[2];
	for (int il = 0; il<1; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*1+il;
			pad_xtsamplesn[index] = new TPad(Form("pad_xtsamplesn%d",index),"pad",il/1.,(2-ir)/2.,(il+1)/1.,(1-ir)/2.);
			pad_xtsamplesn[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_xtsamplesn[index]->SetGridx(1);
			pad_xtsamplesn[index]->SetGridy(1);
		}
	}
	pad_xtsamplesn[0]->cd();
	h2_xtn->Draw("COLZ");
	gr_xt_slicetn->Draw("PSAME");
	gr_xt_slicexn->Draw("PSAME");
	gr_both_cen->Draw("PSAME");
	gr_both_mid->Draw("PSAME");
	gr_both_end->Draw("PSAME");
	f_both_com->Draw("SAME");
	pad_xtsamplesn[1]->cd();
	TH2D * h2_samplingsn = new TH2D("h2_samplingsn",gr_SmF_both->GetTitle(),mNbint,mTmin,mTmax,512,-200,200);
	h2_samplingsn->Draw();
	gr_SmF_both->Draw("PLSAME");
	canv_xtsamplesn->SaveAs("xtsamplesn_"+mRunName+".png");
	canv_xtsamplesn->SaveAs("xtsamplesn_"+mRunName+".pdf");
	// Add zoom-in plots
	TCanvas * canv_xtsamplesnz = new TCanvas("canv_xtsamplesnz","canv_xtsamplesnz",800,1000);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	gr_xt_slicetn->SetMarkerSize(1);
	gr_xt_slicexn->SetMarkerSize(1);
	gr_both_cen->SetMarkerSize(1);
	gr_both_mid->SetMarkerSize(1);
	gr_both_end->SetMarkerSize(1);
	h2_xtn->GetXaxis()->SetRangeUser(-10,60);
	h2_xtn->GetYaxis()->SetRangeUser(0,2);
	h2_xtn->Draw("COLZ");
	gr_xt_slicetn->Draw("PSAME");
	gr_xt_slicexn->Draw("PSAME");
	gr_both_cen->Draw("PSAME");
	gr_both_mid->Draw("PSAME");
	f_both_com->Draw("SAME");
	canv_xtsamplesnz->SaveAs("xtsamplesn_center_"+mRunName+".png");
	canv_xtsamplesnz->SaveAs("xtsamplesn_center_"+mRunName+".pdf");
	h2_xtn->GetXaxis()->SetRangeUser(250,450);
	h2_xtn->GetYaxis()->SetRangeUser(6.5,8);
	h2_xtn->Draw("COLZ");
	gr_xt_slicetn->Draw("PSAME");
	gr_xt_slicexn->Draw("PSAME");
	gr_both_mid->Draw("PSAME");
	gr_both_end->Draw("PSAME");
	f_both_com->Draw("SAME");
	canv_xtsamplesnz->SaveAs("xtsamplesn_end_"+mRunName+".png");
	canv_xtsamplesnz->SaveAs("xtsamplesn_end_"+mRunName+".pdf");
	h2_xtn->GetXaxis()->UnZoom();
	h2_xtn->GetYaxis()->UnZoom();
	// fitting status of slices 
	TCanvas * canv_xtslicesn = new TCanvas("canv_xtslicesn","canv_xtslicesn",1364,768);
	TPad * pad_xtslicesn[6];
	for (int il = 0; il<3; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*3+il;
			pad_xtslicesn[index] = new TPad(Form("pad_xtslicesn%d",index),"pad",il/3.,(2-ir)/2.,(il+1)/3.,(1-ir)/2.);
			pad_xtslicesn[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_xtslicesn[index]->SetGridx(1);
			pad_xtslicesn[index]->SetGridy(1);
		}
	}
	pad_xtslicesn[0]->cd();
	gr_xn_slicexn->Draw("AP");
	pad_xtslicesn[1]->cd();
	gr_xsig_slicexn->GetYaxis()->SetRangeUser(0,mSigTmax);
	gr_xsig_slicexn->Draw("AP");
	pad_xtslicesn[2]->cd();
	gr_xchi2_slicexn->Draw("AP");
	pad_xtslicesn[3]->cd();
	pad_xtslicesn[3]->SetLogy(1);
	gr_nt_slicetn->Draw("AP");
	pad_xtslicesn[4]->cd();
	gr_sigt_slicetn->GetYaxis()->SetRangeUser(0,mSigXmax);
	gr_sigt_slicetn->Draw("AP");
	pad_xtslicesn[5]->cd();
	gr_chi2t_slicetn->Draw("AP");
	canv_xtslicesn->SaveAs("xtslicesn_"+mRunName+".png");
	canv_xtslicesn->SaveAs("xtslicesn_"+mRunName+".pdf");
}

void XTAnalyzer::drawLRB(){
	TCanvas * canv_LRB = new TCanvas("canv_LRB","canv_LRB",800,1000);
	TPad * pad_LRB[2];
	for (int il = 0; il<1; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*1+il;
			pad_LRB[index] = new TPad(Form("pad_LRB%d",index),"pad",il/1.,(2-ir)/2.,(il+1)/1.,(1-ir)/2.);
			pad_LRB[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_LRB[index]->SetGridx(1);
			pad_LRB[index]->SetGridy(1);
		}
	}
	pad_LRB[0]->cd();
	h2_xt->Draw("COLZ");
	gr_both_cen->Draw("PSAME");
	gr_bothL_cen->Draw("PSAME");
	gr_both_mid->Draw("PSAME");
	gr_bothL_mid->Draw("PSAME");
	gr_both_end->Draw("PSAME");
	gr_bothL_end->Draw("PSAME");
	gr_left_end->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
	gr_left_cen->Draw("PSAME");
	gr_right_cen->Draw("PSAME");
	gr_right_mid->Draw("PSAME");
	gr_right_end->Draw("PSAME");
	f_left_com->Draw("SAME");
	f_right_com->Draw("SAME");
	f_both_com->Draw("SAME");
	f_bothL_com->Draw("SAME");
	pad_LRB[1]->cd();
	if (m_RLmB_dx_max<200) m_RLmB_dx_max = 200; // set constant range as +-200 um for later iteration phases.
	TH2D * h2_LRmB_dx = new TH2D("h2_LRmB_dx","XT Differences with Both-side Combined Case",mNbint,mTmin,mTmax,512,-m_RLmB_dx_max*1.05,m_RLmB_dx_max*1.05);
	h2_LRmB_dx->GetXaxis()->SetTitle("T [ns]");
	h2_LRmB_dx->GetYaxis()->SetTitle("#Delta_{X} [um]");
	h2_LRmB_dx->Draw();
	gr_RmB_func->Draw("LSAME");
	gr_LmB_func->Draw("LSAME");
	canv_LRB->SaveAs("LRB_"+mRunName+".png");
	canv_LRB->SaveAs("LRB_"+mRunName+".pdf");
}

void XTAnalyzer::drawIteration(){
	TCanvas * canv_xtiterationn = new TCanvas("canv_xtiterationn","canv_xtiterationn",800,1000);
	TPad * pad_IterN[2];
	for (int il = 0; il<1; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*1+il;
			pad_IterN[index] = new TPad(Form("pad_IterN%d",index),"pad",il/1.,(2-ir)/2.,(il+1)/1.,(1-ir)/2.);
			pad_IterN[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_IterN[index]->SetGridx(1);
			pad_IterN[index]->SetGridy(1);
		}
	}
	pad_IterN[0]->cd();
	h2_xtn->Draw("COLZ");
	f_both_com->Draw("SAME");
	fo_both->SetLineColor(kRed);
	fo_both->SetLineStyle(2);
	fo_both->SetLineWidth(0.3);
	fo_both->Draw("SAME");
	pad_IterN[1]->cd();
	if (m_TmL_B_max<100) m_TmL_B_max = 100; // set constant range as +-100 um for later iteration phases.
	TH2D * h2_IterN = new TH2D("h2_IterN","XT Differences Comparing with Last Time",mNbint,mTmin,mTmax,1024,-m_TmL_B_max*1.05,m_TmL_B_max*1.05);
	h2_IterN->GetXaxis()->SetTitle("T [ns]");
	h2_IterN->GetYaxis()->SetTitle("#Delta_{X} [um]");
	h2_IterN->Draw();
	gr_TmL_both->Draw("LSAME");
	canv_xtiterationn->SaveAs("IterN_"+mRunName+".png");
	canv_xtiterationn->SaveAs("IterN_"+mRunName+".pdf");
	TCanvas * canv_xtiteration = new TCanvas("canv_xtiteration","canv_xtiteration",800,1000);
	TPad * pad_Iter[2];
	for (int il = 0; il<1; il++){
		for (int ir = 0; ir<2; ir++){
			int index = ir*1+il;
			pad_Iter[index] = new TPad(Form("pad_Iter%d",index),"pad",il/1.,(2-ir)/2.,(il+1)/1.,(1-ir)/2.);
			pad_Iter[index]->Draw();
			gStyle->SetPalette(1);
			gStyle->SetOptStat(0);
			gStyle->SetPadTickX(1);
			gStyle->SetPadTickY(1);
			pad_Iter[index]->SetGridx(1);
			pad_Iter[index]->SetGridy(1);
		}
	}
	pad_Iter[0]->cd();
	h2_xt->Draw("COLZ");
	f_left_com->Draw("SAME");
	f_right_com->Draw("SAME");
	fo_left->SetLineColor(kBlack);
	fo_right->SetLineColor(kBlack);
	fo_left->SetLineStyle(2);
	fo_right->SetLineStyle(2);
	fo_left->SetLineWidth(0.3);
	fo_right->SetLineWidth(0.3);
	fo_left->Draw("SAME");
	fo_right->Draw("SAME");
	pad_Iter[1]->cd();
	if (m_TmL_LR_max<100) m_TmL_LR_max = 100; // set constant range as +-100 um for later iteration phases.
	TH2D * h2_Iter = new TH2D("h2_Iter","XT Differences Comparing with Last Time",mNbint,mTmin,mTmax,1024,-m_TmL_LR_max*1.05,m_TmL_LR_max*1.05);
	h2_Iter->GetXaxis()->SetTitle("T [ns]");
	h2_Iter->GetYaxis()->SetTitle("#Delta_{X} [um]");
	h2_Iter->Draw();
	gr_TmL_left->Draw("LSAME");
	gr_TmL_right->Draw("LSAME");
	canv_xtiteration->SaveAs("Iter_"+mRunName+".png");
	canv_xtiteration->SaveAs("Iter_"+mRunName+".pdf");
}

void XTAnalyzer::writeObjects(){
	mOutFile->cd();
	for (int ix = 0; ix<NSLICEX; ix++){
		h_t[ix]->Write();
		h_t_xsum[ix]->Write();
		h_tn[ix]->Write();
		h_tn_xsum[ix]->Write();
	}
	for (int it = 0; it<NSLICET; it++){
		h_x[it]->Write();
		h_x_tsum[it]->Write();
		h_xn[it]->Write();
		h_xn_tsum[it]->Write();
	}
	h2_xt->Write();
	h2_xtn->Write();
	gr_xt_slicex->Write();
	gr_xt_slicet->Write();
	gr_xt_slicexn->Write();
	gr_xt_slicetn->Write();
	gr_sigt_slicetl->Write();
	gr_sigts_slicetl->Write();
	gr_nt_slicetl->Write();
	gr_chi2t_slicetl->Write();
	gr_sigt_slicetr->Write();
	gr_sigts_slicetr->Write();
	gr_nt_slicetr->Write();
	gr_chi2t_slicetr->Write();
	gr_sigt_slicetn->Write();
	gr_nt_slicetn->Write();
	gr_chi2t_slicetn->Write();
	gr_xsig_slicex->Write();
	gr_xn_slicex->Write();
	gr_xchi2_slicex->Write();
	gr_xsig_slicexn->Write();
	gr_xn_slicexn->Write();
	gr_xchi2_slicexn->Write();
	gr_left_end->Write();
	gr_left_mid->Write();
	gr_left_cen->Write();
	gr_right_cen->Write();
	gr_right_mid->Write();
	gr_right_end->Write();
	gr_both_cen->Write();
	gr_both_mid->Write();
	gr_both_end->Write();
	gr_bothL_cen->Write();
	gr_bothL_mid->Write();
	gr_bothL_end->Write();
	gr_LmB_func->Write();
	gr_RmB_func->Write();
	gr_TmL_left->Write();
	gr_TmL_right->Write();
	gr_TmL_both->Write();
	f_left_cen->Write();
	f_right_cen->Write();
	f_both_cen->Write();
	f_left_mid->Write();
	f_right_mid->Write();
	f_both_mid->Write();
	f_left_end->Write();
	f_right_end->Write();
	f_both_end->Write();
	f_left_com->Write();
	f_right_com->Write();
	f_both_com->Write();
	f_bothL_com->Write();
	f_left_deltac->Write();
	f_right_deltac->Write();
	f_both_deltac->Write();
	f_left_delta->Write();
	f_right_delta->Write();
	f_both_delta->Write();
	f_right->Write();
	f_left->Write();
	if (mSaveXT0){
		f_right0->Write();
		f_left0->Write();
		gr_sigts_slicetl0->Write();
		gr_sigts_slicetr0->Write();
	}
}
