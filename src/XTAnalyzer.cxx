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
	if (save) mSaveHists = true;
	else mSaveHists = false;
}

int XTAnalyzer::Initialize(TString runname, int lid, TFile * infile, TFile * outfile, TTree * otree, int xttype, bool savehists, bool saveXT0){
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
	mSigTmax = 15;
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
	mBWX = 0.1;
	mXLEFT = -mBWX*NSLICEX/2;
	mXRIGHT = mBWX*NSLICEX/2;
	mBWT = 3/0.96;
	mTLEFT = -mBWT*NSLICET/2;
	mTRIGHT = mBWT*NSLICET/2;
	// set for binning
	mTmin = -24.5; // t range for one x bin
	mTmax = 800.5;
	mNbint = 264;
	mXmax = 10; // x range for one t bin
	mNbinx = 256;

	// prepare 2D histograms
	h2_xt = new TH2D(Form("h2_xt_%d",mLayerID),"XT Relation",mNbint,mTmin,mTmax,mNbinx*2,-mXmax,mXmax);
	h2_xt->GetXaxis()->SetTitle("T [ns]");
	h2_xt->GetYaxis()->SetTitle("X [mm]");
	h2_xtn = new TH2D(Form("h2_xtn_%d",mLayerID),"XT Relation",mNbint,mTmin,mTmax,mNbinx,0,mXmax);
	h2_xtn->GetXaxis()->SetTitle("T [ns]");
	h2_xtn->GetYaxis()->SetTitle("X [mm]");

	// prepare histograms for slice analysis
	for (int ix = 0; ix<NSLICEX; ix++){
        double fdmin,fdmid,fdmax;
        i2x(ix,fdmin,fdmid,fdmax);
		h_t[ix] = new TH1D(Form("h_t_%d_%d",mLayerID,ix),Form("T_{drift} distribution with DOCA [%.2f~%.2f] mm",fdmin,fdmax),mNbint,mTmin,mTmax);
		h_tn[ix] = new TH1D(Form("h_tn_%d_%d",mLayerID,ix),Form("T_{drift} distribution with DOCA [%.2f~%.2f] mm",fdmin,fdmax),mNbint,mTmin,mTmax);
	}
	for (int it = 0; it<NSLICET; it++){
        double dtmin,dtmid,dtmax;
        i2t(it,dtmin,dtmid,dtmax);
		h_x[it] = new TH1D(Form("h_x_%d_%d",mLayerID,it),Form("DOCA (left) distribution with T_{drift} [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
		h_xn[it] = new TH1D(Form("h_xn_%d_%d",mLayerID,it),Form("DOCA (left) distribution with T_{drift} [%.0f~%.0f] ns",dtmin,dtmax),mNbinx,-mXmax,mXmax);
	}

	// prepare functions for slice analysis
	f_x = myNewTF1("f_x","gaus",-mXmax,mXmax);
	f_t = myNewTF1("f_t","gaus",mTmin,mTmax);
	l_right = new TLine(0,0,1,1);
	l_left = new TLine(0,0,1,1);

	// clear vectors for chosen XT sample points
	v_left_mid_x.clear();
	v_left_mid_t.clear();
	v_right_mid_x.clear();
	v_right_mid_t.clear();
	v_left_end_x.clear();
	v_left_end_t.clear();
	v_right_end_x.clear();
	v_right_end_t.clear();
	v_both_mid_x.clear();
	v_both_mid_t.clear();
	v_both_end_x.clear();
	v_both_end_t.clear();
	v_bothL_mid_x.clear();
	v_bothL_end_x.clear();
	v_LmB_dx.clear();
	v_LmB_t.clear();
	v_RmB_dx.clear();
	v_RmB_t.clear();
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
	f_left_mid = myNewTF1(Form("flm_%d",mLayerID),"pol9",mTmin,mTmax);
	f_right_mid = myNewTF1(Form("frm_%d",mLayerID),"pol9",mTmin,mTmax);
	f_both_mid = myNewTF1(Form("fbm_%d",mLayerID),"pol9",mTmin,mTmax);
	f_left_end = myNewTF1(Form("fle_%d",mLayerID),"pol5",mTmin,mTmax);
	f_right_end = myNewTF1(Form("fre_%d",mLayerID),"pol5",mTmin,mTmax);
	f_both_end = myNewTF1(Form("fbe_%d",mLayerID),"pol5",mTmin,mTmax);

	if (mDebugLevel>=1) printf("XTAnalyzer successfully initialized!\n");
	return 0;
}

void XTAnalyzer::Push(double t, double x){
	if (mDebugLevel>=10) printf("XTAnalyzer::push(%.2e,%.2e)\n",t,x);
	double absx = fabs(x);
	h2_xt->Fill(t,x);
	h2_xtn->Fill(t,absx);
	int ix = x2i(x);
	if (ix>=0&&ix<NSLICEX) h_t[ix]->Fill(t);
	int it = t2i(t,x>0);
	if (it>=0&&it<NSLICET) h_x[it]->Fill(x);
	ix = x2i(absx);
	if (ix>=0&&ix<NSLICEX) h_tn[ix]->Fill(t);
	it = t2i(t,true);
	if (it>=0&&it<NSLICET) h_xn[it]->Fill(absx);
	if (mDebugLevel>=10) printf("            pushed\n",t,x);
}

void XTAnalyzer::Process(void){
	//==========================Taking Samples from Slices==============================
	// fit x histograms, and push to vectors & tree
    TCanvas * canv_fitting = new TCanvas("cfit","cfit",1024,768);
	int midEntries = 100;
	for (int i = 0; i<NSLICET; i++){
		mType = 0;
		double left, right; // find the range: 1/3 of the highest bin
		double divleft,divright;
		i2t(i,divleft,mT,divright);
		mEntries = h_x[i]->Integral();
		mX = 0;
		mSig = 0;
		mChi2 = 0;
		if (mEntries>0){
			getOneThirdLines(h_x[i],left,right);
			h_x[i]->Fit("f_x","qN0","",left,right);
			mX = f_x->GetParameter(1);
			mSig = f_x->GetParameter(2);
			mChi2 = f_x->GetChisquare();
			if (mSaveHists) drawFitting(h_x[i],f_x,canv_fitting,Form("%.0f-%.0f ns, N=%d, x=%.2f mm, #sigma=%.0f um, #chi^{2}=%.1f",divleft,divright,mEntries,mX,mSig*1000,mChi2),Form("h_x%d_%s.png",i,mRunName.Data()),left,right);
		}
		if (mEntries<midEntries){
			mX = h_x[i]->GetMean();
			mSig = h_x[i]->GetRMS();
			mChi2 = 0;
		}
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
		double left, right; // find the range: 1/3 of the highest bin
		double divleft,divright;
		i2x(i,divleft,mX,divright);
		mEntries = h_t[i]->Integral();
		mT = 0;
		mSig = 0;
		mChi2 = 0;
		if (mEntries>0){
			getOneThirdLines(h_t[i],left,right);
			h_t[i]->Fit("f_t","qN0","",left,right);
			mT = f_t->GetParameter(1);
			mSig = f_t->GetParameter(2);
			mChi2 = f_t->GetChisquare();
			if (mSaveHists) drawFitting(h_t[i],f_t,canv_fitting,Form("%.2f-%.2f mm, N=%d, t=%.1f ns, #sigma=%.1f ns, #chi^{2}=%.1f",divleft,divright,mEntries,mT,mSig,mChi2),Form("h_t%d_%s.png",i,mRunName.Data()),left,right);
		}
		if (mEntries<midEntries){
			mT = h_t[i]->GetMean();
			mSig = h_t[i]->GetRMS();
			mChi2 = 0;
		}
		v_n_slicex[i] = mEntries;
		v_sig_slicex[i] = mSig;
		v_chi2_slicex[i] = mChi2;
		v_x_slicex[i] = mX;
		v_t_slicex[i] = mT;
		mOutTree->Fill();
	}
	// fit x histograms for both-side case, and push to vectors & tree
	for (int i = NSLICET/2; i<NSLICET; i++){
		mType = 2;
		double left, right; // find the range: 1/3 of the highest bin
		double divleft,divright;
		i2t(i,divleft,mT,divright);
		mEntries = h_xn[i]->Integral();
		mX = 0;
		mSig = 0;
		mChi2 = 0;
		if (mEntries>0){
			getOneThirdLines(h_xn[i],left,right);
			h_xn[i]->Fit("f_x","qN0","",left,right);
			mX = f_x->GetParameter(1);
			mSig = f_x->GetParameter(2);
			mChi2 = f_x->GetChisquare();
			if (mSaveHists) drawFitting(h_xn[i],f_x,canv_fitting,Form("%.0f-%.0f ns, N=%d, x=%.2f mm, #sigma=%.0f um, #chi^{2}=%.1f",divleft,divright,mEntries,mX,mSig*1000,mChi2),Form("h_xn%d_%s.png",i,mRunName.Data()),left,right);
		}
		if (mEntries<midEntries){
			mX = h_xn[i]->GetMean();
			mSig = h_xn[i]->GetRMS();
			mChi2 = 0;
		}
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
		double left, right; // find the range: 1/3 of the highest bin
		double divleft,divright;
		i2x(i,divleft,mX,divright);
		mEntries = h_tn[i]->Integral();
		mT = 0;
		mSig = 0;
		mChi2 = 0;
		if (mEntries>0){
			getOneThirdLines(h_tn[i],left,right);
			h_tn[i]->Fit("f_t","qN0","",left,right);
			mT = f_t->GetParameter(1);
			mSig = f_t->GetParameter(2);
			mChi2 = f_t->GetChisquare();
			if (mSaveHists) drawFitting(h_tn[i],f_t,canv_fitting,Form("%.2f-%.2f mm, N=%d, t=%.1f ns, #sigma=%.1f ns, #chi^{2}=%.1f",divleft,divright,mEntries,mT,mSig,mChi2),Form("h_tn%d_%s.png",i,mRunName.Data()),left,right);
		}
		if (mEntries<midEntries){
			mT = h_tn[i]->GetMean();
			mSig = h_tn[i]->GetRMS();
			mChi2 = 0;
		}
		v_n_slicexn[i] = mEntries;
		v_sig_slicexn[i] = mSig;
		v_chi2_slicexn[i] = mChi2;
		v_x_slicexn[i] = mX;
		v_t_slicexn[i] = mT;
		mOutTree->Fill();
	}

	//==========================Select Samples==============================
	// select sample points and make graphs
	// FIXME: Currently seperate the mid/end graphs by 7.8 mm line, and search for the real one from 7.5 line.
	double xStart2Turn = 7;
	double t8Left = v_t_slicex[0];
	double t8Right = v_t_slicex[NSLICEX-1];
	double t8Both = v_t_slicexn[NSLICEX-1];
	getT8(t8Left,t8Right,t8Both);
	double t7Left = v_t_slicex[1./mBWX];
	double t7Right = v_t_slicex[NSLICEX-1./mBWX];
	double t7Both = v_t_slicexn[NSLICEX-1./mBWX];
	if (mDebugLevel>=1){
		printf("Before selecting samples:\n");
		printf(" t8l:%.1f, t8r:%.1f, t8b:%.1f\n",t8Left,t8Right,t8Both);
		printf(" t7l:%.1f, t7r:%.1f, t7b:%.1f\n",t7Left,t7Right,t7Both);
	}
	sigmaXReset();
	for (int i = 0; i<NSLICET/2; i++){ // x samples in t slices, left
		if (mDebugLevel>=2) printf("  LR T slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.2f\n",i,v_x_slicet[i],v_t_slicet[i],v_n_slicet[i],v_sig_slicet[i]);
		if (v_n_slicet[i]<mEntriesMin||v_sig_slicet[i]>mSigXmax||v_sig_slicet[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (v_t_slicet[i]>t8Left){ // left end
			if (mDebugLevel>=2) printf("                  t>=%.1f, push to left_end!\n",t8Left);
			v_left_end_x.push_back(v_x_slicet[i]);
			v_left_end_t.push_back(v_t_slicet[i]);
		}
		else if (v_x_slicet[i]<=-xStart2Turn){ // turning part
			if (mDebugLevel>=2) printf("                  x<=%.2f, push to left_mid!\n",-xStart2Turn);
			v_left_mid_x.push_back(v_x_slicet[i]);
			v_left_mid_t.push_back(v_t_slicet[i]);
		}
		sigmaXIncrement(v_t_slicet[i],v_sig_slicet[i],v_n_slicet[i],v_t_slicetls,v_sig_slicetls);
	}
	sigmaXFinalcheck(v_t_slicetls,v_sig_slicetls);
	for (int i = 0; i<NSLICEX; i++){ // t samples in x slices
		if (mDebugLevel>=2) printf("  LR X slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.1f\n",i,v_x_slicex[i],v_t_slicex[i],v_n_slicex[i],v_sig_slicex[i]);
		if (v_n_slicex[i]<mEntriesMin||v_sig_slicex[i]>mSigTmax||v_sig_slicex[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (i<=NSLICEX/2){ // left
			if (v_x_slicex[i]>-xStart2Turn){ // middle part
				if (mDebugLevel>=2) printf("                  x>%.2f, push to left_mid!\n",-xStart2Turn);
				v_left_mid_x.push_back(v_x_slicex[i]);
				v_left_mid_t.push_back(v_t_slicex[i]);
			}
		}
		if (i>=NSLICEX/2){ // right
			if (v_x_slicex[i]<xStart2Turn){ // middle part
				if (mDebugLevel>=2) printf("                  x<%.2f, push to right_mid!\n",xStart2Turn);
				v_right_mid_x.push_back(v_x_slicex[i]);
				v_right_mid_t.push_back(v_t_slicex[i]);
			}
		}
	}
	sigmaXReset();
	for (int i = NSLICET/2; i<NSLICET; i++){ // x samples in t slices, right 
		if (mDebugLevel>=2) printf("  LR T slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.2f\n",i,v_x_slicet[i],v_t_slicet[i],v_n_slicet[i],v_sig_slicet[i]);
		if (v_n_slicet[i]<mEntriesMin||v_sig_slicet[i]>mSigXmax||v_sig_slicet[i]<=0) continue;
		if (v_t_slicet[i]>t8Right){ // right end
			if (mDebugLevel>=2) printf("                  t>=%.1f, push to right_end!\n",t8Right);
			v_right_end_x.push_back(v_x_slicet[i]);
			v_right_end_t.push_back(v_t_slicet[i]);
		}
		else if (v_x_slicet[i]>=xStart2Turn){ // turning part
			if (mDebugLevel>=2) printf("                  x>=%.2f, push to right_mid!\n",xStart2Turn);
			v_right_mid_x.push_back(v_x_slicet[i]);
			v_right_mid_t.push_back(v_t_slicet[i]);
		}
		sigmaXIncrement(v_t_slicet[i],v_sig_slicet[i],v_n_slicet[i],v_t_slicetrs,v_sig_slicetrs);
	}
	sigmaXFinalcheck(v_t_slicetrs,v_sig_slicetrs);
	for (int i = NSLICEX/2; i<NSLICEX; i++){ // t samples in x slices, both-side
		if (mDebugLevel>=2) printf("  BS X slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.1f\n",i,v_x_slicexn[i],v_t_slicexn[i],v_n_slicexn[i],v_sig_slicexn[i]);
		if (v_n_slicexn[i]<mEntriesMin||v_sig_slicexn[i]>mSigTmax||v_sig_slicexn[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (v_x_slicexn[i]<xStart2Turn){ // middle part
			if (mDebugLevel>=2) printf("                  x<%.2f, push to both_mid!\n",xStart2Turn);
			v_both_mid_x.push_back(v_x_slicexn[i]);
			v_bothL_mid_x.push_back(-v_x_slicexn[i]);
			v_both_mid_t.push_back(v_t_slicexn[i]);
		}
	}
	sigmaXReset();
	for (int i = NSLICET/2; i<NSLICET; i++){ // x samples in t slices, both-side
		if (mDebugLevel>=2) printf("  BS T slice[%d]: x=%.2f, t=%.1f, n=%.0f, sig=%.2f\n",i,v_x_slicetn[i],v_t_slicetn[i],v_n_slicetn[i],v_sig_slicetn[i]);
		if (v_n_slicetn[i]<mEntriesMin||v_sig_slicetn[i]>mSigXmax||v_sig_slicetn[i]<=0) continue;
		if (mDebugLevel>=2) printf("                  Passed!\n");
		if (v_t_slicetn[i]>t8Both){ // both-side end
			if (mDebugLevel>=2) printf("                  t>=%.1f, push to both_end!\n",t8Both);
			v_both_end_x.push_back(v_x_slicetn[i]);
			v_bothL_end_x.push_back(-v_x_slicetn[i]);
			v_both_end_t.push_back(v_t_slicetn[i]);
		}
		else if (v_x_slicetn[i]>=xStart2Turn){ // turning part
			if (mDebugLevel>=2) printf("                  x>=%.2f, push to both_mid!\n",xStart2Turn);
			v_both_mid_x.push_back(v_x_slicetn[i]);
			v_bothL_mid_x.push_back(-v_x_slicetn[i]);
			v_both_mid_t.push_back(v_t_slicetn[i]);
		}
		sigmaXIncrement(v_t_slicetn[i],v_sig_slicetn[i],v_n_slicetn[i],v_t_slicetns,v_sig_slicetns);
	}
	sigmaXFinalcheck(v_t_slicetns,v_sig_slicetns);

	// get Left/Right/Both-sides differences by samples
	for (int i = 0; i<v_left_mid_x.size(); i++){
		for (int j = 0; j<v_both_mid_x.size(); j++){
			if (fabs(-v_left_mid_x[i]-v_both_mid_x[j])<1e-7){
				double dt = v_left_mid_t[i]-v_both_mid_t[j];
				int n = 0;
				double vel = 0;
				if (j!=0){
					vel += (v_both_mid_x[j]-v_both_mid_x[j-1])/(v_both_mid_t[j]-v_both_mid_t[j-1]);
					n++;
				}
				if (j!=v_both_mid_x.size()-1){
					vel += (v_both_mid_x[j]-v_both_mid_x[j+1])/(v_both_mid_t[j]-v_both_mid_t[j+1]);
					n++;
				}
				if (n) vel/=n;
				double dx = -dt*vel*1000;
				v_LmB_dx.push_back(dx);
				v_LmB_t.push_back((v_left_mid_t[i]+v_both_mid_t[j])/2);
				if (m_RLmB_dx_max<fabs(dx)) m_RLmB_dx_max = fabs(dx);
				break;
			}
			else if (fabs(v_left_mid_t[i]-v_both_mid_t[j])<1e-7){
				v_LmB_dx.push_back((-v_left_mid_x[i]-v_both_mid_x[j])*1000); // turn to use um
				v_LmB_t.push_back(v_left_mid_t[i]);
				if (m_RLmB_dx_max<fabs(-v_left_mid_x[i]-v_both_mid_x[j])*1000) m_RLmB_dx_max = fabs(-v_left_mid_x[i]-v_both_mid_x[j])*1000;
				break;
			}
		}
	}
	for (int i = 0; i<v_left_end_t.size(); i++){
		for (int j = 0; j<v_both_end_t.size(); j++){
			if (fabs(v_left_end_t[i]-v_both_end_t[j])<1e-7){
				v_LmB_dx.push_back((-v_left_end_x[i]-v_both_end_x[j])*1000); // turn to use um
				v_LmB_t.push_back(v_left_end_t[i]);
				if (m_RLmB_dx_max<fabs(-v_left_end_x[i]-v_both_end_x[j])*1000) m_RLmB_dx_max = fabs(-v_left_end_x[i]-v_both_end_x[j])*1000;
				break;
			}
		}
	}
	for (int i = 0; i<v_right_mid_x.size(); i++){
		for (int j = 0; j<v_both_mid_x.size(); j++){
			if (fabs(v_right_mid_x[i]-v_both_mid_x[j])<1e-7){
				double dt = v_right_mid_t[i]-v_both_mid_t[j];
				int n = 0;
				double vel = 0;
				if (j!=0){
					vel += (v_both_mid_x[j]-v_both_mid_x[j-1])/(v_both_mid_t[j]-v_both_mid_t[j-1]);
					n++;
				}
				if (j!=v_both_mid_x.size()-1){
					vel += (v_both_mid_x[j]-v_both_mid_x[j+1])/(v_both_mid_t[j]-v_both_mid_t[j+1]);
					n++;
				}
				if (n) vel/=n;
				double dx = -dt*vel*1000;
				v_RmB_dx.push_back(dx);
				v_RmB_t.push_back((v_right_mid_t[i]+v_both_mid_t[j])/2);
				if (m_RLmB_dx_max<fabs(dx)) m_RLmB_dx_max = fabs(dx);
				break;
			}
			else if (fabs(v_right_mid_t[i]-v_both_mid_t[j])<1e-7){
				v_RmB_dx.push_back((v_right_mid_x[i]-v_both_mid_x[j])*1000); // turn to use um
				v_RmB_t.push_back(v_right_mid_t[i]);
				if (m_RLmB_dx_max<fabs(v_right_mid_x[i]-v_both_mid_x[j])*1000) m_RLmB_dx_max = fabs(v_right_mid_x[i]-v_both_mid_x[j])*1000;
				break;
			}
		}
	}
	for (int i = 0; i<v_right_end_t.size(); i++){
		for (int j = 0; j<v_both_end_t.size(); j++){
			if (fabs(v_right_end_t[i]-v_both_end_t[j])<1e-7){
				v_RmB_dx.push_back((v_right_end_x[i]-v_both_end_x[j])*1000); // turn to use um
				v_RmB_t.push_back(v_right_end_t[i]);
				if (m_RLmB_dx_max<fabs(v_right_end_x[i]-v_both_end_x[j])*1000) m_RLmB_dx_max = fabs(v_right_end_x[i]-v_both_end_x[j])*1000;
				break;
			}
		}
	}

	//==========================Prepare graphs==============================
	createGraphs();

	//==========================Get XT Functions==============================
	// fit xt functions
	// FIXME: should think about fixing some parameters to let left/right sides XTs go through the same 0 point
	//        now we have poor shape of XT near 0 point, and we'd better rely on both-side XT
	if (!gr_left_end->GetN()) fprintf(stderr,"WARNING: gr_left_end is empty!\n"); else gr_left_end->Fit(Form("fle_%d",mLayerID),"qN0","");;
	if (!gr_left_mid->GetN()) fprintf(stderr,"WARNING: gr_left_mid is empty!\n"); else gr_left_mid->Fit(Form("flm_%d",mLayerID),"qN0","");;
	if (!gr_right_mid->GetN()) fprintf(stderr,"WARNING: gr_right_mid is empty!\n"); else gr_right_mid->Fit(Form("frm_%d",mLayerID),"qN0","");;
	if (!gr_right_end->GetN()) fprintf(stderr,"WARNING: gr_right_end is empty!\n"); else gr_right_end->Fit(Form("fre_%d",mLayerID),"qN0","");;
	if (!gr_both_mid->GetN()) fprintf(stderr,"WARNING: gr_both_mid is empty!\n"); else gr_both_mid->Fit(Form("fbm_%d",mLayerID),"qN0","");;
	if (!gr_both_end->GetN()) fprintf(stderr,"WARNING: gr_both_end is empty!\n"); else gr_both_end->Fit(Form("fbe_%d",mLayerID),"qN0","");;

	// get ranges for XT functions
	f_left_delta = minusPolN(Form("fld_%d",mLayerID),f_left_mid,f_left_end,0,mTmax);
	f_right_delta = minusPolN(Form("frd_%d",mLayerID),f_right_mid,f_right_end,0,mTmax);
	f_both_delta= minusPolN(Form("fbd_%d",mLayerID),f_both_mid,f_both_end,0,mTmax);
	double tTurnLeft = findFirstZero(f_left_delta,t7Left,mTmax,1);
	double tTurnRight = findFirstZero(f_right_delta,t7Right,mTmax,1);
	double tTurnBoth = findFirstZero(f_both_delta,t7Both,mTmax,1);
	double tZeroLeft = findFirstZero(f_left_mid,mTmin,mTmax,1);
	double tZeroRight = findFirstZero(f_right_mid,mTmin,mTmax,1);
	double tZeroBoth = findFirstZero(f_both_mid,mTmin,mTmax,1);
	double tEndLeft = v_left_end_t.size()>0?v_left_end_t[0]:0;
	double tEndRight = v_right_end_t.size()>0?v_right_end_t[v_right_end_t.size()-1]:0;
	double tEndBoth = v_both_end_t.size()>0?v_both_end_t[v_both_end_t.size()-1]:0;
	if (mDebugLevel>=1){
		printf("After selecting samples:\n");
		printf(" t0l:%.1f, t0r:%.1f, t0b:%.1f\n",tZeroLeft,tZeroRight,tZeroBoth);
		printf(" ttl:%.1f, ttr:%.1f, ttb:%.1f\n",tTurnLeft,tTurnRight,tTurnBoth);
		printf(" tel:%.1f, ter:%.1f, teb:%.1f\n",tEndLeft,tEndRight,tEndBoth);
	}

	// combine functions
	f_left_com = combinePolN(Form("flc_%d",mLayerID),f_left_mid,f_left_end,tZeroLeft,tTurnLeft,tEndLeft,mTmin,tEndLeft);
	f_right_com = combinePolN(Form("frc_%d",mLayerID),f_right_mid,f_right_end,tZeroRight,tTurnRight,tEndRight,mTmin,tEndRight);
	f_both_com = combinePolN(Form("fbc_%d",mLayerID),f_both_mid,f_both_end,tZeroBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
	f_bothL_com = combinePolN(Form("fblc_%d",mLayerID),scalePolN(f_both_mid,-1),scalePolN(f_both_end,-1),tZeroBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
	f_left_com->SetLineColor(kMagenta);
	f_both_com->SetLineColor(kBlack);
	f_bothL_com->SetLineColor(kBlack);

	// get the final xt functions
	if (mXTType==0){ // use Left/Right case
		f_left = combinePolN(Form("fl_%d",mLayerID),f_left_mid,f_left_end,tZeroLeft,tTurnLeft,tEndLeft,mTmin,tEndLeft);
		f_right = combinePolN(Form("fr_%d",mLayerID),f_right_mid,f_right_end,tZeroRight,tTurnRight,tEndRight,mTmin,tEndRight);
		if (mSaveXT0){
			f_left0 = combinePolN(Form("fl_%d",0),f_left_mid,f_left_end,tZeroLeft,tTurnLeft,tEndLeft,mTmin,tEndLeft);
			f_right0 = combinePolN(Form("fr_%d",0),f_right_mid,f_right_end,tZeroRight,tTurnRight,tEndRight,mTmin,tEndRight);
		}
	}
	else{ // use Both-Side case
		f_left = combinePolN(Form("fl_%d",mLayerID),scalePolN(f_both_mid,-1),scalePolN(f_both_end,-1),tZeroBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
		f_right = combinePolN(Form("fr_%d",mLayerID),f_both_mid,f_both_end,tZeroBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
		if (mSaveXT0){
			f_left0 = combinePolN(Form("fl_%d",0),scalePolN(f_both_mid,-1),scalePolN(f_both_end,-1),tZeroBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
			f_right0 = combinePolN(Form("fr_%d",0),f_both_mid,f_both_end,tZeroBoth,tTurnBoth,tEndBoth,mTmin,tEndBoth);
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
	for (int i = 0; i<v_left_mid_t.size(); i++){
		double t = v_left_mid_t[i];
		double x = v_left_mid_x[i];
		double xf = f_left_com->Eval(t);
		v_SmF_left_t.push_back(t);
		v_SmF_left_dx.push_back((x-xf)*1000);
	}
	for (int i = 0; i<v_right_mid_t.size(); i++){
		double t = v_right_mid_t[i];
		double x = v_right_mid_x[i];
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
	for (int i = 0; i<v_both_mid_t.size(); i++){
		double t = v_both_mid_t[i];
		double x = v_both_mid_x[i];
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
	// Draw the XT histogram and plots
	drawSamplingsLR();
	// Draw the XT histogram and plots for both-side case
	drawSamplingsB();
	// draw plots to compare left/right/both-side
	drawLRB();
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
}

void XTAnalyzer::i2t(int i, double & dtmin, double & dtmid, double & dtmax){
    dtmin = i*mBWT+mTLEFT;
    if (dtmin<0){
    	dtmin = fabs((i+1)*mBWT+mTLEFT);
    }
	else{
		dtmin = dtmin;
	}
    dtmid = dtmin+mBWT/2.;
    dtmax = dtmid+mBWT/2.;
}

int XTAnalyzer::x2i(double fitD){
    int i=(fitD-mXLEFT)/mBWX;
    if (i>=NSLICEX) i=-999;
    return i;
}

int XTAnalyzer::t2i(double driftT,bool positive){
	if (driftT<0) return -999; // ignore negative time
	if (!positive) driftT = -driftT;
    int i=(driftT-mTLEFT)/mBWT;
    if (i>=NSLICET) i=-999;
    return i;
}

void XTAnalyzer::getOneThirdLines(TH1D* h, double & left, double & right){
    int bmax = h->GetMaximumBin();
    double max = h->GetBinContent(bmax);
    int binl = bmax;
    for (;binl>0; binl--){
        if (h->GetBinContent(binl)<max/3) break;
    }
    int binr = bmax;
    int nbins = h->GetNbinsX();
    for (;binr<=nbins; binr++){
        if (h->GetBinContent(binr)<max/3) break;
    }
    left = h->GetBinLowEdge(binl);
    right = h->GetBinLowEdge(binr+1);
}

double XTAnalyzer::findFirstZero(TF1 * f, double xmin, double xmax, double delta){
	double theX = 0;
	for (double x = xmin; x<xmax; x+=delta){ // At least two solutions. Scan to find the smallest one
		theX = f->GetX(0,xmin,x);
		if (abs(theX-x)>0.5){
			break;
		}
	}
	if (mDebugLevel>=3) printf("findFirstZero: theX = %.2f from [%.2f,%.2f]\n",theX,xmin,xmax);
	return theX;
}

void XTAnalyzer::getT8(double & t8left, double & t8right, double & t8both){
	double t1,t2,x1,x2;
	bool find1 = false;
	bool find2 = false;
	for (int i = 0; i<NSLICEX; i++){
		if (v_sig_slicex[i]<mSigTmax&&v_n_slicex[i]>mEntriesMin){
			if (!find1){
				t1 = v_t_slicex[i];
				x1 = v_x_slicex[i];
				find1 = true;
			}
			else{
				t2 = v_t_slicex[i];
				x2 = v_x_slicex[i];
				find2 = true;
				break;
			}
		}
	}
	if (find1&&fabs(x1+8)<1e-4) t8left = t1;
	else if (find1&&find2){
		t8left = t1+(t2-t1)*(-8-x1)/(x2-x1);
	}
	else{
		t8left = t1;
		fprintf(stderr,"WARNING: cannot find the valid point at -8 mm or two valid points in x slices!\n");
	}
	find1 = false;
	find2 = false;
	for (int i = NSLICEX-1; i>=0; i--){
		if (v_sig_slicex[i]<mSigTmax&&v_n_slicex[i]>mEntriesMin){
			if (!find1){
				t1 = v_t_slicex[i];
				x1 = v_x_slicex[i];
				find1 = true;
			}
			else{
				t2 = v_t_slicex[i];
				x2 = v_x_slicex[i];
				find2 = true;
				break;
			}
		}
	}
	if (find1&&fabs(x1-8)<1e-4) t8right = t1;
	else if (find1&&find2){
		t8right = t1+(t2-t1)*(8-x1)/(x2-x1);
	}
	else{
		t8right = t1;
		fprintf(stderr,"WARNING: cannot find the valid point at 8 mm or two valid points in x slices!\n");
	}
	find1 = false;
	find2 = false;
	for (int i = NSLICEX-1; i>=0; i--){
		if (v_sig_slicexn[i]<mSigTmax&&v_n_slicexn[i]>mEntriesMin){
			if (!find1){
				t1 = v_t_slicexn[i];
				x1 = v_x_slicexn[i];
				find1 = true;
			}
			else{
				t2 = v_t_slicexn[i];
				x2 = v_x_slicexn[i];
				find2 = true;
				break;
			}
		}
	}
	if (find1&&fabs(x1-8)<1e-4) t8both = t1;
	else if (find1&&find2){
		t8both = t1+(t2-t1)*(8-x1)/(x2-x1);
	}
	else{
		t8both = t1;
		fprintf(stderr,"WARNING: cannot find the valid point at 8 mm or two valid points in x slices (with both sides)!\n");
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
	TF1 * f = new TF1(name,Form("pol%d",n),xmin,xmax);
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

TF1 * XTAnalyzer::combinePolN(TString name, TF1 * f1, TF1 * f2, double x0, double x1, double x2, double xmin, double xmax){
	TString form1 = formPolN(f1);
	TString form2 = formPolN(f2);
	double ymax = f2->Eval(x2);
	TF1 * f = new TF1(name,Form("(x>=%.14f&&x<%.14f)*(%s)+(x>=%.14f&&x<%.14f)*(%s)+(x>=%.14f)*%.14e",x0,x1,form1.Data(),x1,x2,form2.Data(),x2,ymax),xmin,xmax);
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
	gr_left_end = myNewTGraph(Form("gr_xt_le_%d",mLayerID),v_left_end_x.size(),&(v_left_end_t[0]),&(v_left_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kMagenta,0.3,kMagenta);
	gr_left_mid = myNewTGraph(Form("gr_xt_lm_%d",mLayerID),v_left_mid_x.size(),&(v_left_mid_t[0]),&(v_left_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kMagenta,0.3,kMagenta);
	gr_right_mid = myNewTGraph(Form("gr_xt_rm_%d",mLayerID),v_right_mid_x.size(),&(v_right_mid_t[0]),&(v_right_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kRed,0.3,kRed);
	gr_right_end = myNewTGraph(Form("gr_xt_re_%d",mLayerID),v_right_end_x.size(),&(v_right_end_t[0]),&(v_right_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kRed,0.3,kRed);
	gr_both_mid = myNewTGraph(Form("gr_xt_bm_%d",mLayerID),v_both_mid_x.size(),&(v_both_mid_t[0]),&(v_both_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_both_end = myNewTGraph(Form("gr_xt_be_%d",mLayerID),v_both_end_x.size(),&(v_both_end_t[0]),&(v_both_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_bothL_mid = myNewTGraph(Form("gr_xt_blm_%d",mLayerID),v_bothL_mid_x.size(),&(v_both_mid_t[0]),&(v_bothL_mid_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_bothL_end = myNewTGraph(Form("gr_xt_ble_%d",mLayerID),v_bothL_end_x.size(),&(v_both_end_t[0]),&(v_bothL_end_x[0]),
			"XT Relation","T [ns]","X [mm]",20,0.3,kBlack,0.3,kBlack);
	gr_LmB = myNewTGraph(Form("gr_LmB_%d",mLayerID),v_LmB_t.size(),&(v_LmB_t[0]),&(v_LmB_dx[0]),
			"XT Differences with Both-side Combined Case","T [ns]","#Delta_{X} [um]",20,0.5,kMagenta,0.5,kMagenta);
	gr_RmB = myNewTGraph(Form("gr_RmB_%d",mLayerID),v_RmB_t.size(),&(v_RmB_t[0]),&(v_RmB_dx[0]),
			"XT Differences with Both-side Combined Case","T [ns]","#Delta_{X} [um]",20,0.5,kRed,0.5,kRed);
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

void XTAnalyzer::drawFitting(TH1D* h,TF1 * f,TCanvas * c,TString title, TString filename,double left, double right){
	c->cd();
	h->SetTitle(title);
	h->Draw();
	double max = h->GetMaximum();
	l_left->SetX1(left);
	l_left->SetX2(left);
	l_left->SetY1(0);
	l_left->SetY2(max);
	l_right->SetX1(right);
	l_right->SetX2(right);
	l_right->SetY1(0);
	l_right->SetY2(max);
	l_left->Draw();
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
	gr_left_end->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
	gr_right_mid->Draw("PSAME");
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
	gr_xsig_slicex->GetYaxis()->SetRangeUser(0,20); // fix at 20 ns for iteration comparisons
	gr_xsig_slicex->Draw("AP");
	pad_xtslices[2]->cd();
	gr_xchi2_slicex->Draw("AP");
	pad_xtslices[3]->cd();
	pad_xtslices[3]->SetLogy(1);
	gr_nt_slicetl->Draw("AP");
	gr_nt_slicetr->Draw("PSAME");
	pad_xtslices[4]->cd();
	gr_sigt_slicetl->GetYaxis()->SetRangeUser(0,0.5); // fix at 0.5 mm for iteration comparisons
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
	gr_both_mid->Draw("PSAME");
	gr_both_end->Draw("PSAME");
	f_both_com->Draw("SAME");
	pad_xtsamplesn[1]->cd();
	TH2D * h2_samplingsn = new TH2D("h2_samplingsn",gr_SmF_both->GetTitle(),mNbint,mTmin,mTmax,512,-200,200);
	h2_samplingsn->Draw();
	gr_SmF_both->Draw("PLSAME");
	canv_xtsamplesn->SaveAs("xtsamplesn_"+mRunName+".png");
	canv_xtsamplesn->SaveAs("xtsamplesn_"+mRunName+".pdf");
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
	gr_xsig_slicexn->GetYaxis()->SetRangeUser(0,20); // fix at 20 ns for iteration comparisons
	gr_xsig_slicexn->Draw("AP");
	pad_xtslicesn[2]->cd();
	gr_xchi2_slicexn->Draw("AP");
	pad_xtslicesn[3]->cd();
	pad_xtslicesn[3]->SetLogy(1);
	gr_nt_slicetn->Draw("AP");
	pad_xtslicesn[4]->cd();
	gr_sigt_slicetn->GetYaxis()->SetRangeUser(0,0.5); // fix at 0.5 mm for iteration comparisons
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
	gr_both_mid->Draw("PSAME");
	gr_both_end->Draw("PSAME");
	gr_bothL_mid->Draw("PSAME");
	gr_bothL_end->Draw("PSAME");
	gr_left_end->Draw("PSAME");
	gr_left_mid->Draw("PSAME");
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
	gr_RmB->Draw("PSAME");
	gr_LmB->Draw("PSAME");
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
		h_tn[ix]->Write();
	}
	for (int it = 0; it<NSLICET; it++){
		h_x[it]->Write();
		h_xn[it]->Write();
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
	gr_right_mid->Write();
	gr_right_end->Write();
	gr_both_mid->Write();
	gr_both_end->Write();
	gr_bothL_mid->Write();
	gr_bothL_end->Write();
	gr_LmB->Write();
	gr_RmB->Write();
	gr_LmB_func->Write();
	gr_RmB_func->Write();
	gr_TmL_left->Write();
	gr_TmL_right->Write();
	gr_TmL_both->Write();
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
