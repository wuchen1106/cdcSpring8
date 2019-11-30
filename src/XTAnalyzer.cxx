#include "XTAnalyzer.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGraphErrors.h"

#include <math.h>

#include "ParameterManager.hxx"
#include "XTManager.hxx"
#include "Log.hxx"

XTAnalyzer::XTAnalyzer(TString runname, TFile * outfile, bool drawDetails)
        :mRunName(runname), mDrawDetails(drawDetails),
        mOutFile(outfile), mOutTree(NULL),
        mX(0), mXerr(0), mT(0), mTerr(0), mSig(0), mChi2(0), mProb(0), mEntries(0), mFunction(0), 
        h2_xt(NULL), gr_left(NULL), gr_right(NULL), gr_rightMinusLeft(NULL),
        f_gausL(NULL), f_gausR(NULL), f_gaus2L(NULL), f_gaus2R(NULL), f_gausBoth(NULL), f_landL(NULL), f_landR(NULL), f_land2L(NULL), f_land2R(NULL), f_landBoth(NULL),
        f_combDoubleGausL(NULL), f_combDoubleGausR(NULL), f_combDoubleGausBoth(NULL),
        f_combDoubleLandL(NULL), f_combDoubleLandR(NULL), f_combDoubleLandBoth(NULL),
        f_combGausLandL(NULL), f_combGausLandR(NULL), f_combGausLandBoth(NULL),
        f_combLandGausL(NULL), f_combLandGausR(NULL), f_combLandGausBoth(NULL),
        f_left_cen(NULL), f_left_mid(NULL), f_left_end(NULL),
        f_right_cen(NULL), f_right_mid(NULL), f_right_end(NULL),
        f_folded_cen(NULL), f_folded_mid(NULL), f_folded_end(NULL),
        f_left(NULL), f_right(NULL), f_folded(NULL)
{
    TString formula("");
    // prepare functions for slice analysis
    // single functions: [0]: height, [1]: mean, [2]: sigma
    f_gausL = myNewTF1("f_gausL","gaus",-1000,1000); f_gausL->SetLineStyle(2); f_gausL->SetLineColor(kRed);
    f_gausR = myNewTF1("f_gausR","gaus",-1000,1000); f_gausR->SetLineStyle(2); f_gausR->SetLineColor(kRed);
    f_gaus2L = myNewTF1("f_gaus2L","gaus",-1000,1000); f_gaus2L->SetLineStyle(2); f_gaus2L->SetLineColor(kBlue);
    f_gaus2R = myNewTF1("f_gaus2R","gaus",-1000,1000); f_gaus2R->SetLineStyle(2); f_gaus2R->SetLineColor(kBlue);
    f_landL = myNewTF1("f_landL","[0]/0.18065564*TMath::Landau(x,[1],[2],false)",-1000,1000); f_landL->SetLineStyle(2); f_landL->SetLineColor(kRed);
    f_landR = myNewTF1("f_landR","[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)",-1000,1000); f_landR->SetLineStyle(2); f_landR->SetLineColor(kRed);
    f_land2L = myNewTF1("f_land2L","[0]/0.18065564*TMath::Landau(x,[1],[2],false)",-1000,1000); f_land2L->SetLineStyle(2); f_land2L->SetLineColor(kBlue);
    f_land2R = myNewTF1("f_land2R","[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)",-1000,1000); f_land2R->SetLineStyle(2); f_land2R->SetLineColor(kBlue);
    // combined functions: [0]: peak height, [1]: peak mean, [2]: peak sigma, [3]: base height (rel), [4]: base mean (rel) (positive means innerward), [5]: base sigma (rel)
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]/2*[0]*exp(-0.5*((x-[1]-[4]*[2])/[5]/[2])*((x-[1]-[4]*[2])/[5]/[2]))"; // base gaussian
    f_combDoubleGausL = myNewTF1("f_combDoubleGausL",formula,-1000,1000); f_combDoubleGausL->SetLineStyle(1); f_combDoubleGausL->SetLineColor(kMagenta);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]/2*[0]*exp(-0.5*((x-[1]+[4]*[2])/[5]/[2])*((x-[1]+[4]*[2])/[5]/[2]))"; // base gaussian
    f_combDoubleGausR = myNewTF1("f_combDoubleGausR",formula,-1000,1000); f_combDoubleGausR->SetLineStyle(1); f_combDoubleGausR->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combDoubleLandL = myNewTF1("f_combDoubleLandL",formula,-1000,1000); f_combDoubleLandL->SetLineStyle(1); f_combDoubleLandL->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)"; // peak landau
    formula += "+[3]*[0]/0.18065564*TMath::Landau(-x,-[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combDoubleLandR = myNewTF1("f_combDoubleLandR",formula,-1000,1000); f_combDoubleLandR->SetLineStyle(1); f_combDoubleLandR->SetLineColor(kMagenta);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combGausLandL = myNewTF1("f_combGausLandL",formula,-1000,1000); f_combGausLandL->SetLineStyle(1); f_combGausLandL->SetLineColor(kMagenta);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]/0.18065564*TMath::Landau(-x,-[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combGausLandR = myNewTF1("f_combGausLandR",formula,-1000,1000); f_combGausLandR->SetLineStyle(1); f_combGausLandR->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]*exp(-0.5*((x-[1]-[4]*[2])/[2]/[5])*((x-[1]-[4]*[2])/[2]/[5]))"; // base gaussian
    f_combLandGausL = myNewTF1("f_combLandGausL",formula,-1000,1000); f_combLandGausL->SetLineStyle(1); f_combLandGausL->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)"; // peak landau
    formula += "+[3]*[0]*exp(-0.5*((x-[1]+[4]*[2])/[2]/[5])*((x-[1]+[4]*[2])/[2]/[5]))"; // base gaussian
    f_combLandGausR = myNewTF1("f_combLandGausR",formula,-1000,1000); f_combLandGausR->SetLineStyle(1); f_combLandGausR->SetLineColor(kMagenta);
    // extend to both sides: add a mirrored part with new position and height
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*exp(-0.5*((x-[4])/[2])*((x-[4])/[2]))"; // peak gaussian
    f_gausBoth = myNewTF1("f_gausBoth",formula,-1000,1000); f_gausBoth->SetLineStyle(1); f_gausBoth->SetLineColor(kBlack);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]/0.18065564*TMath::Landau(-x,-[4],[2],false)"; // peak landau
    f_landBoth = myNewTF1("f_landBoth",formula,-1000,1000); f_landBoth->SetLineStyle(1); f_landBoth->SetLineColor(kBlack);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]*exp(-0.5*((x-[1]-[4]*[2])/[5]/[2])*((x-[1]-[4]*[2])/[5]/[2]))"; // base gaussian
    formula += "+[6]*exp(-0.5*((x-[7])/[2])*((x-[7])/[2]))"; // peak gaussian on the other side
    formula += "+[3]*[6]*exp(-0.5*((x-[7]+[4]*[2])/[5]/[2])*((x-[7]+[4]*[2])/[5]/[2]))"; // base gaussian on the other side
    f_combDoubleGausBoth = myNewTF1("f_combDoubleGausBoth",formula,-1000,1000); f_combDoubleGausBoth->SetLineStyle(1); f_combDoubleGausBoth->SetLineColor(kBlack);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    formula += "+[6]/0.18065564*TMath::Landau(-x,-[7],[2],false)"; // peak landau on the other side
    formula += "+[3]*[6]/0.18065564*TMath::Landau(-x,-[7]+[4]*[2],[2]*[5],false)"; // base landau on the other side
    f_combDoubleLandBoth = myNewTF1("f_combDoubleLandBoth",formula,-1000,1000); f_combDoubleLandBoth->SetLineStyle(1); f_combDoubleLandBoth->SetLineColor(kBlack);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    formula += "+[6]*exp(-0.5*((x-[7])/[2])*((x-[7])/[2]))"; // peak gaussian on the other side
    formula += "+[3]*[6]/0.18065564*TMath::Landau(-x,-[7]+[4]*[2],[2]*[5],false)"; // base landau on the other side
    f_combGausLandBoth = myNewTF1("f_combGausLandBoth",formula,-1000,1000); f_combGausLandBoth->SetLineStyle(1); f_combGausLandBoth->SetLineColor(kBlack);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]*exp(-0.5*((x-[1]-[4]*[2])/[2]/[5])*((x-[1]-[4]*[2])/[2]/[5]))"; // base gaussian
    formula += "+[6]/0.18065564*TMath::Landau(-x,-[7],[2],false)"; // peak landau on the other side
    formula += "+[3]*[6]*exp(-0.5*((x-[7]+[4]*[2])/[2]/[5])*((x-[7]+[4]*[2])/[2]/[5]))"; // base gaussian on the other side
    f_combLandGausBoth = myNewTF1("f_combLandGausBoth",formula,-1000,1000); f_combLandGausBoth->SetLineStyle(1); f_combLandGausBoth->SetLineColor(kBlack);
}

XTAnalyzer::~XTAnalyzer(void){
    if (f_gausL) delete f_gausL;
    if (f_gausR) delete f_gausR;
    if (f_gaus2L) delete f_gaus2L;
    if (f_gaus2R) delete f_gaus2R;
    if (f_gausBoth) delete f_gausBoth;
    if (f_landL) delete f_landL;
    if (f_landR) delete f_landR;
    if (f_land2L) delete f_land2L;
    if (f_land2R) delete f_land2R;
    if (f_landBoth) delete f_landBoth;
    if (f_combDoubleGausL) delete f_combDoubleGausL;
    if (f_combDoubleGausR) delete f_combDoubleGausR;
    if (f_combDoubleGausBoth) delete f_combDoubleGausBoth;
    if (f_combDoubleLandL) delete f_combDoubleLandL;
    if (f_combDoubleLandR) delete f_combDoubleLandR;
    if (f_combDoubleLandBoth) delete f_combDoubleLandBoth;
    if (f_combGausLandL) delete f_combGausLandL;
    if (f_combGausLandR) delete f_combGausLandR;
    if (f_combGausLandBoth) delete f_combGausLandBoth;
    if (f_combLandGausL) delete f_combLandGausL;
    if (f_combLandGausR) delete f_combLandGausR;
    if (f_combLandGausBoth) delete f_combLandGausBoth;
}

void XTAnalyzer::Initialize(void){
    if(mOutTree) {delete mOutTree; mOutTree = NULL;}
    if(h2_xt) {delete h2_xt; h2_xt = NULL;}
    if(gr_left) {delete gr_left; gr_left = NULL;}
    if(gr_right) {delete gr_right ; gr_right = NULL;}
    if(gr_rightMinusLeft) {delete gr_rightMinusLeft; gr_rightMinusLeft = NULL;}
    if(f_left_cen) {delete f_left_cen; f_left_cen = NULL;}
    if(f_right_cen) {delete f_right_cen; f_right_cen = NULL;}
    if(f_folded_cen) {delete f_folded_cen; f_folded_cen = NULL;}
    if(f_left_mid) {delete f_left_mid; f_left_mid = NULL;}
    if(f_right_mid) {delete f_right_mid; f_right_mid = NULL;}
    if(f_folded_mid) {delete f_folded_mid; f_folded_mid = NULL;}
    if(f_left_end) {delete f_left_end; f_left_end = NULL;}
    if(f_right_end) {delete f_right_end; f_right_end = NULL;}
    if(f_folded_end) {delete f_folded_end; f_folded_end = NULL;}
}

int XTAnalyzer::Prepare2DHists(bool reLoad){
    if (reLoad){
        if (!mOutFile){
            MyError("Cannot load XTBins from the file \""<<mOutFile->GetPath()<<"\"");
            return 1;
        }
        h2_xt = (TH2D*)mOutFile->Get("h2_xt"+m_suffix);
        if (!h2_xt){
            MyError("Cannot load \""<<"h2_xt"<<m_suffix<<"\" from the file \""<<mOutFile->GetPath()<<"\"");
            return 1;
        }
    }
    else{
        // about binning
        double bin_t_min = ParameterManager::Get().XTAnalyzerParameters.bin_t_min; // t range for one x bin
        double bin_t_max = ParameterManager::Get().XTAnalyzerParameters.bin_t_max;
        int    bin_t_num = ParameterManager::Get().XTAnalyzerParameters.bin_t_num;
        double bin_x_max = ParameterManager::Get().XTAnalyzerParameters.bin_x_max; // x range for one t bin
        int    bin_x_num = ParameterManager::Get().XTAnalyzerParameters.bin_x_num;
        mOutFile->cd();
        h2_xt = new TH2D("h2_xt"+m_suffix,"XT Relation",bin_t_num,bin_t_min,bin_t_max,bin_x_num,-bin_x_max,bin_x_max);
        h2_xt->GetXaxis()->SetTitle("T [ns]");
        h2_xt->GetYaxis()->SetTitle("X [mm]");
    }
    return 0;
}

int XTAnalyzer::PrepareTree(bool reLoad){
    if (reLoad){
        mOutTree = (TTree*) mOutFile->Get("XTBins"+m_suffix);
        if (!mOutFile){
            MyError("Cannot load \""<<"XTBins"<<m_suffix<<"\" from the file \""<<mOutFile->GetPath()<<"\"");
            return 1;
        }
        mOutTree->SetBranchAddress("x",&mX);
        mOutTree->SetBranchAddress("t",&mT);
        mOutTree->SetBranchAddress("xerr",&mXerr);
        mOutTree->SetBranchAddress("terr",&mTerr);
        mOutTree->SetBranchAddress("sig",&mSig);
        mOutTree->SetBranchAddress("chi2",&mChi2);
        mOutTree->SetBranchAddress("prob",&mProb);
        mOutTree->SetBranchAddress("n",&mEntries);
        mOutTree->SetBranchAddress("func",&mFunction);
        gr_left = (TGraphErrors*) mOutFile->Get(Form("gr%s_l",m_suffix.Data()));
        gr_right = (TGraphErrors*) mOutFile->Get(Form("gr%s_r",m_suffix.Data()));
        gr_rightMinusLeft = (TGraphErrors*) mOutFile->Get(Form("gr%s_rml",m_suffix.Data()));
        if (!gr_left||gr_right){
            MyError("Cannot load XT graphs gr"<<m_suffix<<"_l/r/rml from the file \""<<mOutFile->GetPath()<<"\"");
        }
    }
    else{
        mOutFile->cd();
        mOutTree = new TTree("XTBins"+m_suffix,"XTBins"+m_suffix);
        mOutTree->Branch("x",&mX);
        mOutTree->Branch("t",&mT);
        mOutTree->Branch("xerr",&mXerr);
        mOutTree->Branch("terr",&mTerr);
        mOutTree->Branch("sig",&mSig);
        mOutTree->Branch("chi2",&mChi2);
        mOutTree->Branch("prob",&mProb);
        mOutTree->Branch("n",&mEntries);
        mOutTree->Branch("func",&mFunction);
        double markerSize = 0.5;
        gr_left = new TGraphErrors(); gr_left->SetName(Form("gr%s_l",m_suffix.Data()));
        gr_left->SetMarkerStyle(24);gr_left->SetMarkerSize(markerSize);
        gr_left->SetMarkerColor(kBlue);gr_left->SetLineColor(kBlue);
        gr_right = new TGraphErrors(); gr_right->SetName(Form("gr%s_r",m_suffix.Data()));
        gr_right->SetMarkerStyle(24);gr_right->SetMarkerSize(markerSize);
        gr_right->SetMarkerColor(kRed);gr_right->SetLineColor(kRed);
        gr_rightMinusLeft = new TGraphErrors(); gr_rightMinusLeft->SetName(Form("gr%s_rml",m_suffix.Data()));
        gr_rightMinusLeft->SetMarkerStyle(20);gr_rightMinusLeft->SetMarkerSize(markerSize);
        gr_rightMinusLeft->SetMarkerColor(kBlue);gr_rightMinusLeft->SetLineColor(kBlue);
    }
    return 0;
}

int XTAnalyzer::PrepareXTFunctions(){
    double bin_t_min = ParameterManager::Get().XTAnalyzerParameters.bin_t_min; // t range for one x bin
    double bin_t_max = ParameterManager::Get().XTAnalyzerParameters.bin_t_max;
    int mCentPolN = ParameterManager::Get().XTAnalyzerParameters.xt_center_nPol;
    int mMidPolN = ParameterManager::Get().XTAnalyzerParameters.xt_middle_nPol;
    int mEndPolN = ParameterManager::Get().XTAnalyzerParameters.xt_end_nPol;
    f_left_cen = myNewTF1("flce"+m_suffix,Form("pol%d",mCentPolN),bin_t_min,bin_t_max);
    f_right_cen = myNewTF1("frce"+m_suffix,Form("pol%d",mCentPolN),bin_t_min,bin_t_max);
    f_folded_cen = myNewTF1("ffce"+m_suffix,Form("pol%d",mCentPolN),bin_t_min,bin_t_max);
    f_left_mid = myNewTF1("flm"+m_suffix,Form("pol%d",mMidPolN),bin_t_min,bin_t_max);
    f_right_mid = myNewTF1("frm"+m_suffix,Form("pol%d",mMidPolN),bin_t_min,bin_t_max);
    f_folded_mid = myNewTF1("ffm"+m_suffix,Form("pol%d",mMidPolN),bin_t_min,bin_t_max);
    f_left_end = myNewTF1("fle"+m_suffix,Form("pol%d",mEndPolN),bin_t_min,bin_t_max);
    f_right_end = myNewTF1("fre"+m_suffix,Form("pol%d",mEndPolN),bin_t_min,bin_t_max);
    f_folded_end = myNewTF1("ffe"+m_suffix,Form("pol%d",mEndPolN),bin_t_min,bin_t_max);
    return 0;
}

void XTAnalyzer::Fill(double t, double x){
    if (!h2_xt){
        MyError("h2_xt is not prepared yet!");
        return;
    }
    h2_xt->Fill(t,x);
}

void XTAnalyzer::BinAnalysis(void){
    if (!h2_xt){MyError("h2_xt is not prepared yet!"); return;}
    if (!mOutFile){MyError("output tree is not prepared yet!"); return;}
    MyNamedLog("XTAnalyzer","In XTAnalyzer::BinAnalysis, using h2_xt"<<m_suffix<<" ("<<h2_xt->GetName()<<") and filling to "<<mOutTree->GetName());
    //==========================Preparation work==============================
    // number of ranges to group bins and analyze
    int fitX_nRanges = ParameterManager::Get().XTAnalyzerParameters.fitX_nRanges;
    int bin_x_num = ParameterManager::Get().XTAnalyzerParameters.bin_x_num;
    // prepare a canvas to draw fittings if needed
    TCanvas * canv = new TCanvas("cfit","cfit",1024,768);canv->SetGridx(1);canv->SetGridy(1);
    // record the current style
    int oldStyle = gStyle->GetOptStat();
    gStyle->SetOptStat(1); // we want to draw the stat box for bin-by-bin histograms

    //==========================Fit X along T slices (unfolded)==============================
    int current_nbins_to_fit = ParameterManager::Get().XTAnalyzerParameters.fitX_nBins[0];
    int current_iRange = 0;
    mDrawTmin = 1e14; mDrawXmin = 1e14;
    mDrawTmax = -1e14; mDrawXmax = -1e14;
    for (int i = 1; i<=h2_xt->GetXaxis()->GetNbins(); i+=current_nbins_to_fit){
        bool fitBoth = ParameterManager::Get().XTAnalyzerParameters.fitX_fitBoth[current_iRange];
        // fist, get the T division
        double divleft = h2_xt->GetXaxis()->GetBinLowEdge(i);
        double divmiddle = current_nbins_to_fit%2==0?h2_xt->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xt->GetXaxis()->GetBinCenter(i+current_nbins_to_fit/2);
        double divright = h2_xt->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit);
        if (divleft<ParameterManager::Get().XTAnalyzerParameters.fitX_tSep[current_iRange]){ continue; } // not started yet
        h2_xt->GetXaxis()->SetRangeUser(divleft,divright);
        int nLR = 2; if (fitBoth) nLR=1; // if we fit both sides, then don't loop L/R
        for (int iLR = 0; iLR<nLR; iLR++){
            // second, get the mean time and the 1-D histogram
            double t1 = 0; double terr1 = 0; double t2 = 0; double terr2 = 0; double n1 = 0; double n2 = 0; double nmin = 0;
            TH1D * hist = NULL;
            if (fitBoth){ // fit both sides together
                h2_xt->GetYaxis()->SetRangeUser(-100,0);
                t1 = h2_xt->GetMean(1); terr1 = h2_xt->GetRMS(1); n1 = h2_xt->Integral();
                h2_xt->GetYaxis()->SetRangeUser(0,100);
                t2 = h2_xt->GetMean(1); terr2 = h2_xt->GetRMS(1); n2 = h2_xt->Integral();
                h2_xt->GetYaxis()->UnZoom();
                MyNamedInfo("XTAnalyzer",Form("=> Bin %d~%d: %.1f (%.1f~%.1f) ns, left %.0f entries, right %.0f entries",i,i+current_nbins_to_fit-1,divmiddle,divleft,divright,n1,n2));
                hist = h2_xt->ProjectionY(Form("h_t_B%s_%d",m_suffix.Data(),i));
                nmin = n1>n2?n2:n1;
            }
            else{
                if (iLR==0){
                    h2_xt->GetYaxis()->SetRangeUser(-100,0);
                    hist = h2_xt->ProjectionY(Form("h_t_L%s_%d",m_suffix.Data(),i));
                }
                else{
                    h2_xt->GetYaxis()->SetRangeUser(0,100);
                    hist = h2_xt->ProjectionY(Form("h_t_R%s_%d",m_suffix.Data(),i));
                }
                t1 = h2_xt->GetMean(1); double terr = h2_xt->GetRMS(1); n1 = h2_xt->Integral();
                h2_xt->GetYaxis()->UnZoom();
                MyNamedInfo("XTAnalyzer",Form("=> Bin %d~%d: %.1f (%.1f~%.1f) ns, %s %.0f entries",i,i+current_nbins_to_fit-1,divmiddle,divleft,divright,iLR==0?"left":"right",n1));
                nmin = n1;
            }
            // third, prepare the 1-D hist
            int smooth_times = ParameterManager::Get().XTAnalyzerParameters.fitX_smooth[current_iRange];
            if (smooth_times>0) hist->Smooth(smooth_times);
            if (ParameterManager::Get().XTAnalyzerParameters.fitX_SetEmptyBins[current_iRange]){
                for (int i = 1; i<=hist->GetXaxis()->GetNbins(); i++){
                    double content = hist->GetBinContent(i);
                    if (content == 0) content = 1e-7;
                    hist->SetBinContent(i,content);
                }
            }
            // at last , fit this histogram if needed
            int minEntries = ParameterManager::Get().XTAnalyzerParameters.fitX_minEntries[current_iRange]; // minimum number of entries in one slice to apply fitting function; Otherwise use mean value & RMS instead.
            if (!minEntries||nmin>minEntries){
                double x1,xerr1,sig1,x2,xerr2,sig2; int result;
                if (fitBoth){
                    fitSliceBothSides(hist,x1,xerr1,sig1,x2,xerr2,sig2,mChi2,mProb,result,mFunction,current_iRange);
                    if (mDrawDetails){
                        hist->GetXaxis()->SetRangeUser((x1>0?0:x1)-0.6,(x2<0?0:x2)+0.6);
                        TString x1string = x1<-1?Form("%.2f mm",x1):Form("%.0f um",x1*1000);
                        TString x2string = x2>1?Form("%.2f mm",x2):Form("%.0f um",x2*1000);
                        drawFitting(hist,canv
                                ,Form("t = %.1f (%.1f~%.1f) ns, x_{1}=%s (%.0f um), #sigma_{1}=%.0f um, x_{2}=%s (%.0f um), #sigma_{2}=%.0f um, #chi^{2}=%.0f, prob=%.2f, stat %d",divmiddle,divleft,divright,x1string.Data(),xerr1*1000,sig1*1000,x2string.Data(),xerr2*1000,sig2*1000,mChi2,mProb,result)
                                ,Form("%s.%s%s.png",hist->GetName(),mRunName.Data(),m_suffix.Data())
                                ,mFunction,x1,x2);
                    }
                }
                else{
                    fitSliceSingleSide(hist,x1,xerr1,sig1,mChi2,mProb,result,mFunction,current_iRange,iLR==0);
                    if (mDrawDetails){
                        hist->GetXaxis()->SetRangeUser(x1-0.6,x1+0.6);
                        TString xstring = x1<-1?Form("%.2f mm",x1):Form("%.0f um",x1*1000);
                        drawFitting(hist,canv
                                ,Form("t = %.1f (%.1f~%.1f) ns, x=%s (%.0f um), #sigma=%.0f um, #chi^{2}=%.0f, prob=%.2f, stat %d",t1,divleft,divright,xstring.Data(),xerr1*1000,sig1*1000,mChi2,mProb,result)
                                ,Form("%s.%s%s.png",hist->GetName(),mRunName.Data(),m_suffix.Data())
                                ,mFunction,x1,0,iLR==0);
                    }
                }
                // update the canvas range for drawing samples
                if (mDrawTmin>t1) mDrawTmin = t1;
                if (mDrawTmax<t1) mDrawTmax = t1;
                if (mDrawXmin>x1) mDrawXmin = x1;
                if (mDrawXmax<x1) mDrawXmax = x1;
                mT = t1; mTerr = terr1; mEntries = n1; mX = x1; mXerr = xerr1; mSig = sig1;
                mOutTree->Fill();
                if (fitBoth){
                    if (mDrawTmin>t2) mDrawTmin = t2;
                    if (mDrawTmax<t2) mDrawTmax = t2;
                    if (mDrawXmin>x2) mDrawXmin = x2;
                    if (mDrawXmax<x2) mDrawXmax = x2;
                    mT = t2; mTerr = terr2; mEntries = n2; mX = x2; mXerr = xerr2; mSig = sig2;
                    mOutTree->Fill();
                }
            }
        }
        h2_xt->GetXaxis()->UnZoom();
        if (divright>=ParameterManager::Get().XTAnalyzerParameters.fitX_tSep[current_iRange+1]){ // switch to the next range
            current_iRange++;
            if (current_iRange>=fitX_nRanges){
                MyNamedInfo("XTAnalyzer",Form("Division right side %.1f ns assed the edge of the last range[%d]: t = %.1f ns. Stop sampling.",divright,fitX_nRanges,ParameterManager::Get().XTAnalyzerParameters.fitX_tSep[current_iRange+1]));
                break;
            }
            current_nbins_to_fit = ParameterManager::Get().XTAnalyzerParameters.fitX_nBins[current_iRange];
            MyNamedInfo("XTAnalyzer",Form("Entering range %d, nbins to combine changed to %d",current_iRange,current_nbins_to_fit));
        }
    }
    mDrawTmin-=10;
    mDrawTmax+=10;
    mDrawXmin-=0.5;
    mDrawXmax+=0.5;

    //==========================Draw the samples==============================
    formXTGraphs();
    gStyle->SetOptStat(0);
    drawSamples();
    // set the style back
    gStyle->SetOptStat(oldStyle);
}

void XTAnalyzer::FitXT(){
}

void XTAnalyzer::Write(){
    mOutFile->cd();
    if(mOutTree) mOutTree->Write();
    if(h2_xt) h2_xt->Write();
    if(gr_left) gr_left->Write();
    if(gr_right) gr_right->Write();
    if(gr_rightMinusLeft) gr_rightMinusLeft->Write();
    if(f_left_cen) f_left_cen->Write();
    if(f_right_cen) f_right_cen->Write();
    if(f_folded_cen) f_folded_cen->Write();
    if(f_left_mid) f_left_mid->Write();
    if(f_right_mid) f_right_mid->Write();
    if(f_folded_mid) f_folded_mid->Write();
    if(f_left_end) f_left_end->Write();
    if(f_right_end) f_right_end->Write();
    if(f_folded_end) f_folded_end->Write();
}

double XTAnalyzer::interpolate(const TGraphErrors * graph, double theX){
    double theY = 0;
    // check the first point
    double x,y;
    graph->GetPoint(0,x,y);
    bool isLeft = x>theX; // if the asked position is on the left side of the first point.
    double firstX = x; double firstY = y;
    double prevX = x; double prevY = y;
    bool found = false;
    for (int i = 1; i<graph->GetN(); i++){
        graph->GetPoint(i,x,y);
        if ((isLeft&&x<theX)||(!isLeft&&x>theX)){ // moved to the other side of the asked position, then interpolate
            theY = (y*(x-theX)+prevY*(theX-prevX))/(x-prevX);
            found = true;
            break;
        }
        prevX = x;
        prevY = y;
    }
    if (!found){ // didn't cross the asked position
        theY = fabs(theX-prevX)>fabs(theX-firstX)?firstY:prevY;
    }
    return theY;
}

void XTAnalyzer::plusGraph(TGraphErrors * gn, const TGraphErrors * gl, const TGraphErrors * gr, double sl, double sr){
    gn->Set(gl->GetN());
    for (int iPoint = 0; iPoint<gl->GetN(); iPoint++){
        double xl,yl; double xerr,yerr;
        gl->GetPoint(iPoint,xl,yl);
        yl*=sl;
        xerr = gl->GetErrorX(iPoint);
        yerr = gl->GetErrorY(iPoint);
        double yr = sr*interpolate(gr,xl);
        gn->SetPoint(iPoint,xl,yl+yr);
        gn->SetPointError(iPoint,xerr,yerr);
    }
}

void XTAnalyzer::drawSamples(){
    double markerSize = 0.5;
    // draw the unfolded histogram first
    TCanvas * canv= new TCanvas("cgraph","cgraph",1024,768);canv->SetGridx(1);canv->SetGridy(1);
    gPad->SetGridx(1);gPad->SetGridy(1);
    h2_xt->GetXaxis()->SetRangeUser(mDrawTmin,mDrawTmax);
    h2_xt->GetYaxis()->SetRangeUser(mDrawXmin,mDrawXmax);
    h2_xt->Draw("COLZ");
    gr_left->Draw("PLSAME");
    gr_right->Draw("PLSAME");
    canv->SaveAs(Form("result/sample2D_%s%s.png",mRunName.Data(),m_suffix.Data()));
    canv->SaveAs(Form("result/sample2D_%s%s.pdf",mRunName.Data(),m_suffix.Data()));

    // show the difference among left and right
    TH2D * h2_bkg_rightMinusLeft = new TH2D("h2_bkg_rightMinusLeft","Mean of the two side of XT relation",512,mDrawTmin,mDrawTmax,1024,-1,1);
    h2_bkg_rightMinusLeft->Draw();
    gr_rightMinusLeft->Draw("PLSAME");
    canv->SaveAs(Form("result/sampleDiff_%s%s.png",mRunName.Data(),m_suffix.Data()));

    mOutTree->SetMarkerStyle(20);mOutTree->SetMarkerSize(markerSize);
    canv= new TCanvas("csample","csample",1024,768);
    // firstly, get the maximum number of entries in one sample
    int mEntriesMax = 0;
    for (Long64_t iEntry = 0; iEntry<mOutTree->GetEntries(); iEntry++){
        mOutTree->GetEntry(iEntry);
        if (mEntries>mEntriesMax) mEntriesMax = mEntries;
    }
    canv->Divide(2,2);
    canv->cd(1);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_entriesT = new TH2D("h2_bkg_entriesT","Number of entries in each T slice",512,mDrawTmin,mDrawTmax,1024,0,mEntriesMax*1.1);
    h2_bkg_entriesT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("n:t","x>=0&&n>0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("n:t","x<0&&n>0","PSAME");
    canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_sigT = new TH2D("h2_bkg_sigT","#sigma of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,0.8);
    h2_bkg_sigT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("sig:t","x>=0&&n>0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("sig:t","x<0&&n>0","PSAME");
    canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_chi2T = new TH2D("h2_bkg_chi2T","#chi^{2} of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,150);
    h2_bkg_chi2T->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("chi2:t","x>=0&&n>0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("chi2:t","x<0&&n>0","PSAME");
    canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_probT = new TH2D("h2_bkg_probT","p-value of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,1);
    h2_bkg_probT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("prob:t","x>=0&&n>0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("prob:t","x<0&&n>0","PSAME");
    canv->SaveAs(Form("result/sampleAtt_%s%s.png",mRunName.Data(),m_suffix.Data()));
}

void XTAnalyzer::formXTGraphs(){
    // set parameters
    int    graph_n_min = ParameterManager::Get().XTAnalyzerParameters.graph_n_min;
    double graph_chi2_max = ParameterManager::Get().XTAnalyzerParameters.graph_chi2_max;
    double graph_prob_min = ParameterManager::Get().XTAnalyzerParameters.graph_prob_min;

    // make graphs from different samplings: left/right/folded TIMES time/space
    gr_left->Set(mOutTree->GetEntries(Form("func!=0&&n>=%d&&(%d||chi2<=%.7e)&&prob>=%.7e&&x<%.7e",graph_n_min,graph_chi2_max?0:1,graph_chi2_max,graph_prob_min,0.)));
    gr_right->Set(mOutTree->GetEntries(Form("func!=0&&n>=%d&&(%d||chi2<=%.7e)&&prob>=%.7e&&x>%.7e",graph_n_min,graph_chi2_max?0:1,graph_chi2_max,graph_prob_min,0.)));
    int count_left = 0;
    int count_right = 0;
    MyNamedVerbose("XTAnalyzer","Looping in the output tree again, "<<mOutTree->GetEntries()<<" entries");
    for (Long64_t iEntry = 0; iEntry<mOutTree->GetEntries(); iEntry++){
        mOutTree->GetEntry(iEntry);
        if (mFunction!=0&&mEntries>=graph_n_min&&(!graph_chi2_max||mChi2<=graph_chi2_max)&&mProb>=graph_prob_min){
            if (mX<0){
                gr_left->SetPoint(count_left,mT,mX);
                gr_left->SetPointError(count_left,mTerr,mXerr);
                count_left++;
            }
            if (mX>0){
                gr_right->SetPoint(count_right,mT,mX);
                gr_right->SetPointError(count_right,mTerr,mXerr);
                count_right++;
            }
        }
    }
    plusGraph(gr_rightMinusLeft,gr_right,gr_left);
}

TF1 * XTAnalyzer::myNewTF1(TString name, TString form, double left, double right){
    TF1 * f = new TF1(name,form,left,right);
    f->SetNpx(1024);
    f->SetNumberFitPoints(1024);
    f->SetLineWidth(1);
    f->SetLineColor(kRed);
    return f;
}

TF1 * XTAnalyzer::fitSliceSingleSide(TH1D * hist, double & x1,double & xerr1,double & sig1,double & chi2,double & prob,int & result, int & functionType, int iRange, bool isLeft){
    MyNamedVerbose("XTAnalyzer",Form("fitSliceSingleSide %s side of histogram %s",hist->GetName(),isLeft?"left":"right"));
    // get the shape of the input histogram
    double h1 = 0;
    getMaximum(hist,x1,h1,-11,11);
    sig1 = hist->GetRMS(); xerr1 = sig1;
    // prepare the function
    TF1 * f = NULL; TF1 * fp = NULL; TF1 * fb = NULL;
    double p_h_m = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_height_middle[iRange];
    double p_h_l = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_height_left[iRange];
    double p_h_r = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_height_right[iRange];
    double p_s_m = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_sigma_middle[iRange];
    double p_s_l = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_sigma_left[iRange];
    double p_s_r = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_sigma_right[iRange];
    double p_x_d = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_mean_range[iRange];
    double b_h_m = ParameterManager::Get().XTAnalyzerParameters.fitX_base_height_middle[iRange];
    double b_h_l = ParameterManager::Get().XTAnalyzerParameters.fitX_base_height_left[iRange];
    double b_h_r = ParameterManager::Get().XTAnalyzerParameters.fitX_base_height_right[iRange];
    double b_s_m = ParameterManager::Get().XTAnalyzerParameters.fitX_base_sigma_middle[iRange];
    double b_s_l = ParameterManager::Get().XTAnalyzerParameters.fitX_base_sigma_left[iRange];
    double b_s_r = ParameterManager::Get().XTAnalyzerParameters.fitX_base_sigma_right[iRange];
    double b_x_d = ParameterManager::Get().XTAnalyzerParameters.fitX_base_mean_range[iRange];
    functionType = ParameterManager::Get().XTAnalyzerParameters.fitX_functionType[iRange];
    int nPars = 0;
    if (functionType==kGaussian||functionType==kLandau){ // single function
        if (functionType==kGaussian){ // fit with gaussian
            if (isLeft) {f = f_gausL;} else {f = f_gausR;}
        }
        else if (functionType==kLandau){ // fit with landau
            if (isLeft) {f = f_landL;} else {f = f_landR;}
        }
        f->SetParameters(h1*p_h_m,x1,p_s_m?p_s_m:sig1);
        // peak on the left
        f->SetParLimits(0,h1*p_h_l,h1*p_h_r);
        f->SetParLimits(1,x1-p_x_d,x1+p_x_d);
        f->SetParLimits(2,p_s_l,p_s_r);
        nPars = 3;
    }
    else{ // composite function
        // Note that if fitFunctionType is kOptimal, do fitting with kGaussianPlusLandau first, then test kLandauPlusGaussian
        if (functionType==kGaussianPlusLandau||functionType==kOptimal){ // fit with gaussian plus landau
            if (isLeft){
                f = f_combGausLandL;
                fp = f_gausL;
                fb = f_land2L;
            }
            else{
                f = f_combGausLandR;
                fp = f_gausR;
                fb = f_land2R;
            }
        }
        else if (functionType==kLandauPlusGaussian){ // fit with gaussian plus landau
            if (isLeft){
                f = f_combLandGausL;
                fp = f_landL;
                fb = f_gaus2L;
            }
            else{
                f = f_combLandGausR;
                fp = f_landR;
                fb = f_gaus2R;
            }
        }
        else if (functionType==kDoubleGaussian){ // fit with double gaussian
            if (isLeft){
                f = f_combDoubleGausL;
                fp = f_gausL;
                fb = f_gaus2L;
            }
            else{
                f = f_combDoubleGausR;
                fp = f_gausR;
                fb = f_gaus2R;
            }
        }
        else if (functionType==kDoubleLandau){ // fit with double gaussian
            if (isLeft){
                f = f_combDoubleLandL;
                fp = f_landL;
                fb = f_land2L;
            }
            else{
                f = f_combDoubleLandR;
                fp = f_landR;
                fb = f_land2R;
            }
        }
        f->SetParameters(h1*p_h_m,x1,p_s_m?p_s_m:sig1, // peak on the left
                         b_h_m,0.,b_s_m?b_s_m:0.1); // base part (relative to peak)
        // peak part
        f->SetParLimits(0,h1*p_h_l,h1*p_h_r);
        f->SetParLimits(1,x1-p_x_d,x1+p_x_d);
        f->SetParLimits(2,p_s_l,p_s_r);
        // base part (relative to peak)
        f->SetParLimits(3,b_h_l,b_h_r);
        f->SetParLimits(4,-b_x_d,b_x_d);
        f->SetParLimits(5,b_s_l,b_s_r);
        nPars = 6;
    }
    prob = 0;
    for (int iTry = 0; iTry<10; iTry++){
        if (prob>1e-20) break;
        MyNamedVerbose("XTAnalyzer","Before fitting. Set parameters:");
        for (int i = 0; i<nPars; i++){
            double parmin,parmax;
            f->GetParLimits(i,parmin,parmax);
            if (iTry!=0){
                double ratio = gRandom->Uniform(0.4,0.6);
                f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
            }
            MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
        }
        result = hist->Fit(f,"qN0");
        chi2 = f->GetChisquare();
        prob = f->GetProb();
        MyNamedVerbose("XTAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2,prob,result));
        for (int i = 0; i<nPars; i++){
            MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
        }
    }

    if (functionType==kOptimal){ // try LandauPlusGaussian in this case
        MyNamedVerbose("XTAnalyzer",Form("fitFunctionType is kOptimal. Try LandauPlusGaussian"));
        if (isLeft){
            f = f_combLandGausL;
            fp = f_landL;
            fb = f_gaus2L;
        }
        else{
            f = f_combLandGausR;
            fp = f_landR;
            fb = f_gaus2R;
        }
        f->SetParameters(h1*p_h_m,x1,p_s_m?p_s_m:sig1, // peak on the left
                         b_h_m,0.,b_s_m?b_s_m:0.1); // base part (relative to peak)
        // peak part
        f->SetParLimits(0,h1*p_h_l,h1*p_h_r);
        f->SetParLimits(1,x1-p_x_d,x1+p_x_d);
        f->SetParLimits(2,p_s_l,p_s_r);
        // base part (relative to peak)
        f->SetParLimits(3,b_h_l,b_h_r);
        f->SetParLimits(4,-b_x_d,b_x_d);
        f->SetParLimits(5,b_s_l,b_s_r);
        int resultTemp = 0;
        double probTemp = 0;
        double chi2Temp = 0;
        for (int iTry = 0; iTry<10; iTry++){
            if (probTemp>1e-20) break;
            MyNamedVerbose("XTAnalyzer","Before fitting. Set parameters:");
            for (int i = 0; i<nPars; i++){
                double parmin,parmax;
                f->GetParLimits(i,parmin,parmax);
                if (iTry!=0){
                    double ratio = gRandom->Uniform(0.4,0.6);
                    f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
                }
                MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
            }
            resultTemp = hist->Fit(f,"qN0");
            chi2Temp = f->GetChisquare();
            probTemp = f->GetProb();
            MyNamedVerbose("XTAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2Temp,probTemp,resultTemp));
            for (int i = 0; i<nPars; i++){
                MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
            }
        }
        if (chi2Temp<chi2){
            functionType = kLandauPlusGaussian;
            chi2 = chi2Temp;
            prob = probTemp;
            result = resultTemp;
        }
        else{ // change the function back
            functionType = kGaussianPlusLandau;
            if (isLeft){
                f = f_combGausLandL;
                fp = f_gausL;
                fb = f_land2L;
            }
            else{
                f = f_combGausLandR;
                fp = f_gausR;
                fb = f_land2R;
            }
        }
    }

    if (functionType==kGaussian||functionType==kLandau){
        double height = f->GetParameter(0);
        x1 = f->GetParameter(1); xerr1 = f->GetParError(1);
        sig1 = f->GetParameter(2);
        MyNamedVerbose("XTAnalyzer",Form(" peak: height %.3e mean %.3e sig %.3e",height,x1,sig1));
    }
    else{
        double fit_peak_height = f->GetParameter(0);
        double fit_peak_mean = f->GetParameter(1); double fit_peak_meanErr = f->GetParError(1);
        double fit_peak_sigma = f->GetParameter(2);
        double fit_base_heightRel = f->GetParameter(3);
        double fit_base_meanRel = f->GetParameter(4); double fit_base_meanRelErr = f->GetParError(4);
        double fit_base_sigmaRel = f->GetParameter(5);
        double fit_base_height = fit_peak_height*fit_base_heightRel;
        double fit_base_mean = fit_peak_mean+(isLeft?1:-1)*fit_base_meanRel*fit_peak_sigma;
        double fit_base_sigma = fit_base_sigmaRel*fit_peak_sigma;
        // get mean error and sigma
        xerr1 = fit_peak_meanErr;
        getMeanRMS(f,hist,x1,sig1);
        // set the single functions
        fp->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma);
        fb->SetParameters(fit_base_height,fit_base_mean,fit_base_sigma);
        MyNamedVerbose("XTAnalyzer",Form(" peak: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
        MyNamedVerbose("XTAnalyzer",Form(" base: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));
    }

    MyNamedVerbose("XTAnalyzer",Form(" Result:  mean %.3e +- %.3e sig %.3e",x1,xerr1,sig1));

    return f;
}

TF1 * XTAnalyzer::fitSliceBothSides(TH1D * hist, double & x1,double & xerr1,double & sig1,double & x2,double & xerr2,double & sig2,double & chi2,double & prob, int & result, int & functionType, int iRange){
    MyNamedVerbose("XTAnalyzer",Form("fitSliceBothSides of histogram %s",hist->GetName()));
    // get the shape of the input histogram
    double h1 = 0;
    getMaximum(hist,x1,h1,-11,0);
    hist->GetXaxis()->SetRangeUser(-11,0);
    sig1 = hist->GetRMS(); xerr1 = sig1;
    double h2 = 0;
    getMaximum(hist,x2,h2,0,11);
    hist->GetXaxis()->SetRangeUser(0,11);
    sig2 = hist->GetRMS(); xerr2 = sig2;
    hist->GetXaxis()->UnZoom();
    // prepare the function
    TF1 * f = NULL; TF1 * fpl = NULL; TF1 * fpr = NULL; TF1 * fbl = NULL; TF1 * fbr = NULL; TF1 * fl = NULL; TF1 * fr = NULL;
    double p_h_m = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_height_middle[iRange];
    double p_h_l = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_height_left[iRange];
    double p_h_r = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_height_right[iRange];
    double p_s_m = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_sigma_middle[iRange];
    double p_s_l = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_sigma_left[iRange];
    double p_s_r = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_sigma_right[iRange];
    double p_x_d = ParameterManager::Get().XTAnalyzerParameters.fitX_peak_mean_range[iRange];
    double b_h_m = ParameterManager::Get().XTAnalyzerParameters.fitX_base_height_middle[iRange];
    double b_h_l = ParameterManager::Get().XTAnalyzerParameters.fitX_base_height_left[iRange];
    double b_h_r = ParameterManager::Get().XTAnalyzerParameters.fitX_base_height_right[iRange];
    double b_s_m = ParameterManager::Get().XTAnalyzerParameters.fitX_base_sigma_middle[iRange];
    double b_s_l = ParameterManager::Get().XTAnalyzerParameters.fitX_base_sigma_left[iRange];
    double b_s_r = ParameterManager::Get().XTAnalyzerParameters.fitX_base_sigma_right[iRange];
    double b_x_d = ParameterManager::Get().XTAnalyzerParameters.fitX_base_mean_range[iRange];
    int fitFunctionType = ParameterManager::Get().XTAnalyzerParameters.fitX_functionType[iRange];
    int nPars = 0;
    if (fitFunctionType==kGaussian||fitFunctionType==kLandau){ // single function
        if (fitFunctionType==kGaussian){ // fit with gaussian
            functionType = kGaussianBothSides;
            f = f_gausBoth;
            fl = f_gausL;
            fr = f_gausR;
        }
        else if (fitFunctionType==kLandau){ // fit with landau
            functionType = kLandauBothSides;;
            f = f_landBoth;
            fl = f_landL;
            fr = f_landR;
        }
        f->SetParameters(h1*p_h_m,x1,p_s_m?p_s_m:sig1, // peak on the left
                         h2*p_h_m,x2); // peak on the right
        // peak on the left
        f->SetParLimits(0,h1*p_h_l,h1*p_h_r);
        f->SetParLimits(1,x1-p_x_d,x1+p_x_d);
        f->SetParLimits(2,p_s_l,p_s_r);
        // peak on the right
        f->SetParLimits(3,h2*p_h_l,h2*p_h_r);
        f->SetParLimits(4,x2-p_x_d,x2+p_x_d);
        nPars = 5;
    }
    else{ // composite function
        // Note that if fitFunctionType is kOptimal, do fitting with kGaussianPlusLandau first, then test kLandauPlusGaussian
        if (fitFunctionType==kGaussianPlusLandau||fitFunctionType==kOptimal){ // fit with gaussian plus landau
            functionType = kGaussianPlusLandauBothSides;
            f = f_combGausLandBoth;
            fpl = f_gausL;
            fpr = f_gausR;
            fbl = f_land2L;
            fbr = f_land2R;
            fl = f_combGausLandL;
            fr = f_combGausLandR;
        }
        else if (fitFunctionType==kLandauPlusGaussian){ // fit with gaussian plus landau
            functionType = kLandauPlusGaussianBothSides;
            f = f_combLandGausBoth;
            fpl = f_landL;
            fpr = f_landR;
            fbl = f_gaus2L;
            fbr = f_gaus2R;
            fl = f_combLandGausL;
            fr = f_combLandGausR;
        }
        else if (fitFunctionType==kDoubleGaussian){ // fit with double gaussian
            functionType = kDoubleGaussianBothSides;;
            f = f_combDoubleGausBoth;
            fpl = f_gausL;
            fpr = f_gausR;
            fbl = f_gaus2L;
            fbr = f_gaus2R;
            fl = f_combDoubleGausL;
            fr = f_combDoubleGausR;
        }
        else if (fitFunctionType==kDoubleLandau){ // fit with double gaussian
            functionType = kDoubleLandauBothSides;;
            f = f_combDoubleLandBoth;
            fpl = f_landL;
            fpr = f_landR;
            fbl = f_land2L;
            fbr = f_land2R;
            fl = f_combDoubleLandL;
            fr = f_combDoubleLandR;
        }
        f->SetParameters(h1*p_h_m,x1,p_s_m?p_s_m:sig1, // peak on the left
                         b_h_m,0.,b_s_m?b_s_m:0.1, // base part (relative to peak)
                         h2*p_h_m,x2); // peak on the right
        // peak on the left
        f->SetParLimits(0,h1*p_h_l,h1*p_h_r);
        f->SetParLimits(1,x1-p_x_d,x1+p_x_d);
        f->SetParLimits(2,p_s_l,p_s_r);
        // base part (relative to peak)
        f->SetParLimits(3,b_h_l,b_h_r);
        f->SetParLimits(4,-b_x_d,b_x_d);
        f->SetParLimits(5,b_s_l,b_s_r);
        // peak on the right
        f->SetParLimits(6,h2*p_h_l,h2*p_h_r);
        f->SetParLimits(7,x2-p_x_d,x2+p_x_d);
        nPars = 8;
    }
    prob = 0;
    for (int iTry = 0; iTry<10; iTry++){
        if (prob>1e-20) break;
        MyNamedVerbose("XTAnalyzer","Before fitting. Set parameters:");
        for (int i = 0; i<nPars; i++){
            double parmin,parmax;
            f->GetParLimits(i,parmin,parmax);
            if (iTry!=0){
                double ratio = gRandom->Uniform(0.4,0.6);
                f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
            }
            MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
        }
        result = hist->Fit(f,"qN0");
        chi2 = f->GetChisquare();
        prob = f->GetProb();
        MyNamedVerbose("XTAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2,prob,result));
        for (int i = 0; i<nPars; i++){
            MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
        }
    }

    if (fitFunctionType==kOptimal){ // try LandauPlusGaussian in this case
        MyNamedVerbose("XTAnalyzer",Form("fitFunctionType is kOptimal. Try LandauPlusGaussian"));
        functionType = kLandauPlusGaussianBothSides;
        f = f_combLandGausBoth;
        fpl = f_landL;
        fpr = f_landR;
        fbl = f_gaus2L;
        fbr = f_gaus2R;
        fl = f_combLandGausL;
        fr = f_combLandGausR;
        f->SetParameters(h1*p_h_m,x1,p_s_m?p_s_m:sig1, // peak on the left
                b_h_m,0.,b_s_m?b_s_m:0.1, // base part (relative to peak)
                h2*p_h_m,x2); // peak on the right
        // peak on the left
        f->SetParLimits(0,h1*p_h_l,h1*p_h_r);
        f->SetParLimits(1,x1-p_x_d,x1+p_x_d);
        f->SetParLimits(2,p_s_l,p_s_r);
        // base part (relative to peak)
        f->SetParLimits(3,b_h_l,b_h_r);
        f->SetParLimits(4,-b_x_d,b_x_d);
        f->SetParLimits(5,b_s_l,b_s_r);
        // peak on the right
        f->SetParLimits(6,h2*p_h_l,h2*p_h_r);
        f->SetParLimits(7,x2-p_x_d,x2+p_x_d);
        int resultTemp = 0;
        double probTemp = 0;
        double chi2Temp = 0;
        for (int iTry = 0; iTry<10; iTry++){
            if (probTemp>1e-20) break;
            MyNamedVerbose("XTAnalyzer","Before fitting. Set parameters:");
            for (int i = 0; i<nPars; i++){
                double parmin,parmax;
                f->GetParLimits(i,parmin,parmax);
                if (iTry!=0){
                    double ratio = gRandom->Uniform(0.4,0.6);
                    f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
                }
                MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
            }
            result = hist->Fit(f,"qN0");
            chi2Temp = f->GetChisquare();
            probTemp = f->GetProb();
            MyNamedVerbose("XTAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2Temp,probTemp,resultTemp));
            for (int i = 0; i<nPars; i++){
                MyNamedVerbose("XTAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
            }
        }
        if (chi2Temp<chi2){
            chi2 = chi2Temp;
            prob = probTemp;
            result = resultTemp;
        }
        else{ // change the function back
            functionType = kGaussianPlusLandauBothSides;
            f = f_combGausLandBoth;
            fpl = f_gausL;
            fpr = f_gausR;
            fbl = f_land2L;
            fbr = f_land2R;
            fl = f_combGausLandL;
            fr = f_combGausLandR;
        }
    }

    if (fitFunctionType==kGaussian||fitFunctionType==kLandau){
        // now let's get the left peak
        double height = f->GetParameter(0);
        x1 = f->GetParameter(1); xerr1 = f->GetParError(1);
        sig1 = f->GetParameter(2);
        fl->SetParameters(height,x1,sig1);
        MyNamedVerbose("XTAnalyzer",Form(" peakL: height %.3e mean %.3e sig %.3e",height,x1,sig1));
        height = f->GetParameter(3);
        x2 = f->GetParameter(4); xerr2 = f->GetParError(1);
        sig2 = sig1;
        fr->SetParameters(height,x2,sig2);
        MyNamedVerbose("XTAnalyzer",Form(" peakR: height %.3e mean %.3e sig %.3e",height,x2,sig2));
    }
    else{
        // first, get the left peak
        double fit_peak_height = f->GetParameter(0);
        double fit_peak_mean = f->GetParameter(1); double fit_peak_meanErr = f->GetParError(1);
        double fit_peak_sigma = f->GetParameter(2);
        double fit_base_heightRel = f->GetParameter(3);
        double fit_base_meanRel = f->GetParameter(4); double fit_base_meanRelErr = f->GetParError(4);
        double fit_base_sigmaRel = f->GetParameter(5);
        double fit_base_height = fit_peak_height*fit_base_heightRel;
        double fit_base_mean = fit_peak_mean+fit_base_meanRel*fit_peak_sigma;
        double fit_base_sigma = fit_base_sigmaRel*fit_peak_sigma;
        // get mean error and sigma
        xerr1 = fit_peak_meanErr;
        fl->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma,fit_base_heightRel,fit_base_meanRel,fit_base_sigmaRel);
        getMeanRMS(fl,hist,x1,sig1);
        // set the single functions
        fpl->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma);
        fbl->SetParameters(fit_base_height,fit_base_mean,fit_base_sigma);
        MyNamedVerbose("XTAnalyzer",Form(" peakL: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
        MyNamedVerbose("XTAnalyzer",Form(" baseL: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));

        // second, get the left peak
        fit_peak_height = f->GetParameter(6);
        fit_peak_mean = f->GetParameter(7); fit_peak_meanErr = f->GetParError(7);
        fit_base_height = fit_peak_height*fit_base_heightRel;
        fit_base_mean = fit_peak_mean-fit_base_meanRel*fit_peak_sigma; // flipped
        // get mean error and sigma
        xerr2 = fit_peak_meanErr;
        fr->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma,fit_base_heightRel,fit_base_meanRel,fit_base_sigmaRel);
        getMeanRMS(fr,hist,x2,sig2);
        // set the single functions
        fpr->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma);
        fbr->SetParameters(fit_base_height,fit_base_mean,fit_base_sigma);
        MyNamedVerbose("XTAnalyzer",Form(" peakR: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
        MyNamedVerbose("XTAnalyzer",Form(" baseR: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));
    }

    MyNamedVerbose("XTAnalyzer",Form(" Left Part:  mean %.3e +- %.3e sig %.3e",x1,xerr1,sig1));
    MyNamedVerbose("XTAnalyzer",Form(" Right Part: mean %.3e +- %.3e sig %.3e",x2,xerr2,sig2));

    return f;
}

void XTAnalyzer::drawFitting(TH1D* hist, TCanvas * c,TString title, TString filename, int function, double center1, double center2, bool isLeft){
    if (!hist) fprintf(stderr,"ERROR: in drawFitting, input histogram does not exist!\n");
    if (!c) fprintf(stderr,"ERROR: in drawFitting, input canvas does not exist!\n");
    if (!hist||!c) return;
    c->cd();
    int oldStyle = gStyle->GetOptStat();
    gStyle->SetOptStat(1);
    hist->SetTitle(title);
    hist->Draw();
    double max = hist->GetMaximum();
    TLine * l_center1= new TLine(center1,0,center1,max);
    l_center1->SetLineStyle(2); l_center1->SetLineColor(kRed);
    l_center1->Draw();
    TLine * l_center2= new TLine(center2,0,center2,max);
    l_center2->SetLineStyle(2); l_center2->SetLineColor(kRed);
    if (function>=kGaussianBothSides){ // fit both sides
        l_center2->Draw();
    }
    if (function==kGaussianPlusLandauBothSides){
        f_combGausLandBoth->Draw("SAME");
        f_gausL->Draw("SAME");
        f_land2L->Draw("SAME");
        f_combGausLandL->Draw("SAME");
        f_gausR->Draw("SAME");
        f_land2R->Draw("SAME");
        f_combGausLandR->Draw("SAME");
    }
    else if (function==kLandauPlusGaussianBothSides){
        f_combLandGausBoth->Draw("SAME");
        f_landL->Draw("SAME");
        f_gaus2L->Draw("SAME");
        f_combLandGausL->Draw("SAME");
        f_landR->Draw("SAME");
        f_gaus2R->Draw("SAME");
        f_combLandGausR->Draw("SAME");
    }
    else if (function==kDoubleGaussianBothSides){
        f_combDoubleGausBoth->Draw("SAME");
        f_gausL->Draw("SAME");
        f_gaus2L->Draw("SAME");
        f_combDoubleGausL->Draw("SAME");
        f_gausR->Draw("SAME");
        f_gaus2R->Draw("SAME");
        f_combDoubleGausR->Draw("SAME");
    }
    else if (function==kDoubleLandauBothSides){
        f_combDoubleLandBoth->Draw("SAME");
        f_landL->Draw("SAME");
        f_land2L->Draw("SAME");
        f_combDoubleLandL->Draw("SAME");
        f_landR->Draw("SAME");
        f_land2R->Draw("SAME");
        f_combDoubleLandR->Draw("SAME");
    }
    else if (function==kGaussianBothSides){
        f_gausBoth->Draw("SAME");
        f_gausL->Draw("SAME");
        f_gausR->Draw("SAME");
    }
    else if (function==kLandauBothSides){
        f_landBoth->Draw("SAME");
        f_landL->Draw("SAME");
        f_landR->Draw("SAME");
    }
    else if (function==kGaussianPlusLandau){
        if (isLeft){
            f_combGausLandL->Draw("SAME");
            f_gausL->Draw("SAME");
            f_land2L->Draw("SAME");
        }
        else{
            f_combGausLandR->Draw("SAME");
            f_gausR->Draw("SAME");
            f_land2R->Draw("SAME");
        }
    }
    else if (function==kLandauPlusGaussian){
        if (isLeft){
            f_combLandGausL->Draw("SAME");
            f_landL->Draw("SAME");
            f_gaus2L->Draw("SAME");
        }
        else{
            f_combLandGausR->Draw("SAME");
            f_landR->Draw("SAME");
            f_gaus2R->Draw("SAME");
        }
    }
    else if (function==kDoubleGaussian){
        if (isLeft){
            f_combDoubleGausL->Draw("SAME");
            f_gausL->Draw("SAME");
            f_gaus2L->Draw("SAME");
        }
        else{
            f_combDoubleGausR->Draw("SAME");
            f_gausR->Draw("SAME");
            f_gaus2R->Draw("SAME");
        }
    }
    else if (function==kDoubleLandau){
        if (isLeft){
            f_combDoubleLandL->Draw("SAME");
            f_landL->Draw("SAME");
            f_land2L->Draw("SAME");
        }
        else{
            f_combDoubleLandR->Draw("SAME");
            f_landR->Draw("SAME");
            f_land2R->Draw("SAME");
        }
    }
    else if (function==kGaussian){
        if (isLeft){
            f_gausL->Draw("SAME");
        }
        else{
            f_gausR->Draw("SAME");
        }
    }
    else if (function==kLandau){
        if (isLeft){
            f_landL->Draw("SAME");
        }
        else{
            f_landR->Draw("SAME");
        }
    }
    c->SaveAs(filename);
    gStyle->SetOptStat(oldStyle);
    delete l_center1;
    delete l_center2;
}

void XTAnalyzer::getMaximum(TH1D * hist, double & position, double & maximum, double left, double right){
    // find the highest bin in the given region
    // NOTICE: always using 2 adjacent bins to help smoothing
    position = 0;
    maximum = -1e14;
    int nBins = hist->GetNbinsX();
    for (int i = 1+1;i<=nBins-1; i++){
        double x = hist->GetBinCenter(i);
        if (x<left||x>right) continue;
        double value = (hist->GetBinContent(i-1)+hist->GetBinContent(i)+hist->GetBinContent(i+1))/3;
        if (maximum<value){
            position = x;
            maximum = value;
        }
    }
}

void XTAnalyzer::getMeanRMS(TF1 * f, const TH1D * hist, double & mean, double & sigma){
    mean = f->GetMaximumX(-10,10);
    TH1D * histNew = new TH1D(*hist);
    histNew->Reset();
    for (int i = 0; i<1e4; i++){
        double x = f->GetRandom();
        double content = hist->GetBinContent(histNew->FindBin(x));
        if (content) histNew->Fill(x);
    }
    sigma = histNew->GetRMS();
    histNew->GetXaxis()->SetRangeUser(mean-sigma*3,mean+sigma*3);
    sigma = histNew->GetRMS();
    delete histNew;
}
