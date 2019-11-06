#include "XTBinAnalyzer.hxx"

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

#include <math.h>

#include "ParameterManager.hxx"
#include "XTManager.hxx"
#include "Log.hxx"

XTBinAnalyzer::XTBinAnalyzer(TString runname, TFile * outfile, bool savehists, bool drawDetails)
        :mRunName(runname), mSaveHists(savehists), mDrawDetails(drawDetails),
        mOutFile(outfile), mOutTree(NULL), 
        mLayerID(0), mCellID(0), mX(0), mXerr(0), mT(0), mTerr(0), mSig(0), mChi2(0), mEntries(0), mType(0), mFunction(0), 
        minEntries(0),
        m_bin_t_min(0), m_bin_t_max(0), m_bin_t_num(0), m_bin_x_min(0), m_bin_x_max(0), m_bin_x_num(0), 
        m_bin_t_fit_num(1),m_bin_t_fit_num_tail(1),m_bin_t_tailTime(800),m_bin_x_fit_num(1),
        m_bin_t_landTmin(0),m_bin_t_landTmax(800),
        m_bin_t_ratio(0.5),m_bin_x_ratio(0.5),
        m_bin_x_landXmin(0),m_bin_x_landXmax(10),
        h2_xt(NULL), h2_xtn(NULL),
        f_gaus(0), 
        f_land(0), 
        f_landF(0), 
        f_2gaus(0)
{
    mOutFile->cd();
    mOutTree = new TTree("XTBins","XTBins");
    mOutTree->Branch("x",&mX);
    mOutTree->Branch("t",&mT);
    mOutTree->Branch("xerr",&mXerr);
    mOutTree->Branch("terr",&mTerr);
    mOutTree->Branch("lid",&mLayerID);
    mOutTree->Branch("wid",&mCellID);
    mOutTree->Branch("sig",&mSig);
    mOutTree->Branch("chi2",&mChi2);
    mOutTree->Branch("n",&mEntries);
    mOutTree->Branch("type",&mType);
    mOutTree->Branch("func",&mFunction);
}

XTBinAnalyzer::~XTBinAnalyzer(void){
    if (h2_xt) delete h2_xt;
    if (h2_xtn) delete h2_xtn;
}

int XTBinAnalyzer::Initialize(void){
    minEntries = ParameterManager::Get().XTAnalyzerParameters.bin_n_min; // minimum number of entries in one slice to apply fitting function; Otherwise use mean value & RMS instead.
    // about binning
    m_bin_t_min = ParameterManager::Get().XTAnalyzerParameters.bin_t_min; // t range for one x bin
    m_bin_t_max = ParameterManager::Get().XTAnalyzerParameters.bin_t_max;
    m_bin_t_num = ParameterManager::Get().XTAnalyzerParameters.bin_t_num;
    m_bin_x_min = ParameterManager::Get().XTAnalyzerParameters.bin_x_min;
    m_bin_x_max = ParameterManager::Get().XTAnalyzerParameters.bin_x_max; // x range for one t bin
    m_bin_x_num = ParameterManager::Get().XTAnalyzerParameters.bin_x_num;
    // about projection
    m_bin_t_fit_num = ParameterManager::Get().XTAnalyzerParameters.bin_t_fit_num;
    m_bin_t_tailTime = ParameterManager::Get().XTAnalyzerParameters.bin_t_tailTime;
    m_bin_t_fit_num_tail = ParameterManager::Get().XTAnalyzerParameters.bin_t_fit_num_tail;
    m_bin_x_fit_num = ParameterManager::Get().XTAnalyzerParameters.bin_x_fit_num;
    // about the range for using Landau function
    //   time slice to fit space
    m_bin_t_landTmin = ParameterManager::Get().XTAnalyzerParameters.bin_t_landTmin;
    m_bin_t_landTmax = ParameterManager::Get().XTAnalyzerParameters.bin_t_landTmax;
    //   space slice to fit time
    m_bin_x_landXmin = ParameterManager::Get().XTAnalyzerParameters.bin_x_landXmin;
    m_bin_x_landXmax = ParameterManager::Get().XTAnalyzerParameters.bin_x_landXmax;
    // about fitting range
    m_bin_t_ratio = ParameterManager::Get().XTAnalyzerParameters.bin_t_ratio;
    m_bin_x_ratio = ParameterManager::Get().XTAnalyzerParameters.bin_x_ratio;

    // prepare 2D histograms
    //if (h2_xt) delete h2_xt;
    //if (h2_xtn) delete h2_xtn;
    mOutFile->cd();
    h2_xt = new TH2D(Form("h2_xt_%d",mLayerID),"XT Relation",m_bin_t_num,m_bin_t_min,m_bin_t_max,m_bin_x_num,-m_bin_x_max,m_bin_x_max);
    h2_xt->GetXaxis()->SetTitle("T [ns]");
    h2_xt->GetYaxis()->SetTitle("X [mm]");
    h2_xtn = new TH2D(Form("h2_xtn_%d",mLayerID),"XT Relation",m_bin_t_num,m_bin_t_min,m_bin_t_max,m_bin_x_num/2+1,m_bin_x_min,m_bin_x_max);
    h2_xtn->GetXaxis()->SetTitle("T [ns]");
    h2_xtn->GetYaxis()->SetTitle("X [mm]");

    // prepare functions for slice analysis
    f_gaus = myNewTF1("fgaus","gaus",-m_bin_t_max,m_bin_t_max);
    f_land = myNewTF1("fland","landau",-m_bin_t_max,m_bin_t_max);
    f_landF = myNewTF1("flandF","landau(-x)",-m_bin_t_max,m_bin_t_max);
    f_2gaus = myNewTF1("f2gaus","(x>[1])*[0]*exp(-0.5*((x-[1])/[2])**2)+(x<=[1])*[0]*exp(-0.5*((x-[1])/[3])**2)",-m_bin_t_max,m_bin_t_max);

    MyNamedLog("XTBinAnalyzer","XTBinAnalyzer successfully initialized!");
    return 0;
}

void XTBinAnalyzer::Fill(double t, double x){
    double absx = fabs(x);
    h2_xt->Fill(t,x);
    h2_xtn->Fill(t,absx);
}

void XTBinAnalyzer::Process(void){
    MyNamedLog("XTBinAnalyzer","In XTBinAnalyzer::Process");
    // prepare a canvas to draw fittings if needed
    TCanvas * canv_fitting = new TCanvas("cfit","cfit",1024,768);canv_fitting->SetGridx(1);canv_fitting->SetGridy(1);
    TString drawTitle;
    int oldStyle = gStyle->GetOptStat();
    gStyle->SetOptStat(1);
    // make a directory to hold the slices
    TDirectory * oDir = NULL;mOutFile->cd();oDir = mOutFile->mkdir(Form("slices_%d",mLayerID));oDir->cd();
    int current_nbins_to_fit = 1;
    int previousFunction = kGaussian;
    //==========================Fit X along T slices (unfolded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kTimeSliceUnfolded; // set the type of this round of fittings
    for (int iLR = 0; iLR<2; iLR++){
        if (iLR == 0){ // left side
            h2_xt->GetYaxis()->SetRangeUser(-m_bin_x_max,-m_bin_x_min);
        }
        else{ // right side
            h2_xt->GetYaxis()->SetRangeUser(m_bin_x_min,m_bin_x_max);
        }
        previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
        current_nbins_to_fit = m_bin_t_fit_num;
        for (int i = 1; i<=h2_xt->GetXaxis()->GetNbins(); i+=current_nbins_to_fit){
            // fist, get the x histogram for this time slice
            double divleft = h2_xt->GetXaxis()->GetBinLowEdge(i);
            double divmiddle = current_nbins_to_fit%2==0?h2_xt->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xt->GetXaxis()->GetBinCenter(i+current_nbins_to_fit/2);
            double divright = h2_xt->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit);
            if (divright>=m_bin_t_tailTime) current_nbins_to_fit = m_bin_t_fit_num_tail; // switch to tail mode;
            TH1D * h = h2_xt->ProjectionY(Form("h_x_%d_%s_%d",mLayerID,iLR==0?"L":"R",i),i,i+current_nbins_to_fit-1);
            mEntries = h->Integral();
            // second, get the average T value within this T division and the X fit range
            h2_xt->GetXaxis()->SetRangeUser(divleft,divright);
            mT = h2_xt->GetMean(1);
            mTerr = h2_xt->GetRMS(1);
            h2_xt->GetXaxis()->UnZoom();
            drawTitle = Form("time slice (%.1f,%.1f) ns, mean %.1f ns",divleft,divright,mT);
            MyNamedInfo("XTBinAnalyzer","=>"<<h->GetName()<<": "<<drawTitle<<", "<<mEntries<<" entries");
            // then, fit this histogram
            double left,right;
            fitSlice(h,canv_fitting,drawTitle,mX,mXerr,mSig,mChi2,left,right,mFunction,previousFunction,m_bin_t_ratio,mT>=m_bin_t_landTmin&&mT<=m_bin_t_landTmax,mX>0);
            // at last, fill the tree
            mOutTree->Fill();
        }
    }
    // set the 2-D histogram back after this analysis
    h2_xt->GetYaxis()->SetRangeUser(-m_bin_x_max,m_bin_x_max);

    //==========================Fit X along T slices (folded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kTimeSliceFolded; // set the type of this round of fittings
    previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
    current_nbins_to_fit = m_bin_t_fit_num;
    for (int i = 1; i<=h2_xtn->GetXaxis()->GetNbins(); i+=current_nbins_to_fit){
        // fist, get the x histogram for this time slice
        double divleft = h2_xtn->GetXaxis()->GetBinLowEdge(i);
        double divmiddle = current_nbins_to_fit%2==0?h2_xtn->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xtn->GetXaxis()->GetBinCenter(i+current_nbins_to_fit/2);
        double divright = h2_xtn->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit);
        if (divright>=m_bin_t_tailTime) current_nbins_to_fit = m_bin_t_fit_num_tail; // switch to tail mode;
        TH1D * h = h2_xtn->ProjectionY(Form("h_x_%d_F_%d",mLayerID,i),i,i+current_nbins_to_fit-1);
        mEntries = h->Integral();
        // second, get the average T value within this T division and the X fit range
        h2_xtn->GetXaxis()->SetRangeUser(divleft,divright);
        mT = h2_xtn->GetMean(1);
        mTerr = h2_xtn->GetRMS(1);
        h2_xtn->GetXaxis()->UnZoom();
        drawTitle = Form("time slice (%.1f,%.1f) ns, mean %.1f ns, ",divleft,divright,mT);
        MyNamedInfo("XTBinAnalyzer","=>"<<h->GetName()<<": "<<drawTitle<<", "<<mEntries<<" entries");
        // then, fit this histogram
        double left,right;
        fitSlice(h,canv_fitting,drawTitle,mX,mXerr,mSig,mChi2,left,right,mFunction,previousFunction,m_bin_t_ratio,mT>=m_bin_t_landTmin&&mT<=m_bin_t_landTmax,true);
        // at last, fill the tree
        mOutTree->Fill();
    }

    //==========================Fit T along X slices (unfolded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kSpaceSliceUnfolded; // set the type of this round of fittings
    previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
    current_nbins_to_fit = m_bin_x_fit_num;
    for (int i = 1; i<=h2_xt->GetYaxis()->GetNbins(); i+=current_nbins_to_fit){
        // fist, get the x histogram for this time slice
        double divleft = h2_xt->GetYaxis()->GetBinLowEdge(i);
        double divmiddle = current_nbins_to_fit%2==0?h2_xt->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xt->GetYaxis()->GetBinCenter(i+current_nbins_to_fit/2);
        double divright = h2_xt->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit);
        TH1D * h = h2_xt->ProjectionX(Form("h_t_%d_U_%d",mLayerID,i),i,i+current_nbins_to_fit-1);
        mEntries = h->Integral();
        // second, get the average T value within this T division and the X fit range
        h2_xt->GetYaxis()->SetRangeUser(divleft,divright);
        mX = h2_xt->GetMean(2);
        mXerr = h2_xt->GetRMS(2);
        h2_xt->GetYaxis()->UnZoom();
        drawTitle = Form("x slice (%.3f,%.3f) mm, mean %.3f mm, ",divleft,divright,mX);
        MyNamedInfo("XTBinAnalyzer","=>"<<h->GetName()<<": "<<drawTitle<<", "<<mEntries<<" entries");
        // then, fit this histogram
        double left,right;
        fitSlice(h,canv_fitting,drawTitle,mT,mTerr,mSig,mChi2,left,right,mFunction,previousFunction,m_bin_t_ratio,fabs(mX)>=m_bin_x_landXmin&&fabs(mX)<=m_bin_x_landXmax,false);
        // at last, fill the tree
        mOutTree->Fill();
    }

    //==========================Fit T along X slices (folded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kSpaceSliceFolded; // set the type of this round of fittings
    previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
    current_nbins_to_fit = m_bin_x_fit_num;
    for (int i = 1; i<=h2_xtn->GetYaxis()->GetNbins(); i+=current_nbins_to_fit){
        // fist, get the x histogram for this time slice
        double divleft = h2_xtn->GetYaxis()->GetBinLowEdge(i);
        double divmiddle = current_nbins_to_fit%2==0?h2_xtn->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xtn->GetYaxis()->GetBinCenter(i+current_nbins_to_fit/2);
        double divright = h2_xtn->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit);
        TH1D * h = h2_xtn->ProjectionX(Form("h_t_%d_F_%d",mLayerID,i),i,i+current_nbins_to_fit-1);
        mEntries = h->Integral();
        // second, get the average T value within this T division and the X fit range
        h2_xtn->GetYaxis()->SetRangeUser(divleft,divright);
        mX = h2_xtn->GetMean(2);
        mXerr = h2_xtn->GetRMS(2);
        h2_xtn->GetYaxis()->UnZoom();
        drawTitle = Form("x slice (%.3f,%.3f) mm, mean %.3f mm, ",divleft,divright,mX);
        MyNamedInfo("XTBinAnalyzer","=>"<<h->GetName()<<": "<<drawTitle<<", "<<mEntries<<" entries");
        // then, fit this histogram
        double left,right;
        fitSlice(h,canv_fitting,drawTitle,mT,mTerr,mSig,mChi2,left,right,mFunction,previousFunction,m_bin_t_ratio,fabs(mX)>=m_bin_x_landXmin&&fabs(mX)<=m_bin_x_landXmax,false);
        // at last, fill the tree
        mOutTree->Fill();
    }
    gStyle->SetOptStat(oldStyle);
}

TF1 * XTBinAnalyzer::myNewTF1(TString name, TString form, double left, double right){
    TF1 * f = new TF1(name,form,left,right);
    f->SetNpx(1024);
    f->SetNumberFitPoints(1024);
    f->SetLineWidth(1);
    f->SetLineColor(kRed);
    return f;
}

void XTBinAnalyzer::fitSlice(TH1D * h,TCanvas * canv, TString drawTitle, double &mean, double &error, double &sigma, double &chi2, double &left, double &right, int &function, int &previousFunction, double ratio, bool tryLandau, bool flippedLandau){
    //    1. get mean position and left right boundary at the given ratio
    int entries = h->Integral();
    double RMS = h->GetRMS();
    fitSliceFloat(h,mean,error,sigma,chi2,left,right,ratio);
    MyNamedInfo("XTBinAnalyzer",Form("      after fitSliceFloat: mean=%.2e+-%.2e, half height range %.2e~%.2e, sig=%.2e, chi2=%.2e",mean,error,left,right,sigma,chi2));
    mFunction = kNone;
    //    2. if the number of entries within the fitting range is enough, go to fitting
    if (!minEntries||entries>minEntries){
        TF1 * f = fitSliceGaus(h,mean,error,sigma,chi2,left,right);// Try with Gaussian function first
        h->GetYaxis()->SetRangeUser(0,h->GetBinContent(h->GetMaximumBin())*1.1);h->GetXaxis()->SetRangeUser(mean-RMS*5,mean+RMS*5); // for drawing
        if (mSaveHists) drawFitting(h,f,canv,Form("%s, fit in (%.2e,%.2e), mean=%.2e+-%.2e, #sigma=%.2e, #chi^{2}=%.2e",drawTitle.Data(),left,right,mean,error,sigma,chi2),Form("%s_Gauss_%s.layer%d.png",h->GetName(),mRunName.Data(),mLayerID),left,mean,right);
        function = kGaussian;
        if (tryLandau){ // If it's within the Landau function range, then try with Landau function
            double tempMean = mean; double tempError = error; double tempSig = sigma; double tempChi2 = chi2; double tempLeft = left; double tempRight = right;
            int tempFunction = kNone;
            if (flippedLandau){ // use flipped Landau function
                f = fitSliceLandF(h,tempMean,tempError,tempSig,tempChi2,tempLeft,tempRight);
                tempFunction = kLandauFlipped;
            }
            else{ // use Landau function
                f = fitSliceLand(h,tempMean,tempError,tempSig,tempChi2,tempLeft,tempRight);
                tempFunction = kLandau;
            }
            h->GetYaxis()->SetRangeUser(0,h->GetBinContent(h->GetMaximumBin())*1.1);h->GetXaxis()->SetRangeUser(mean-RMS*5,mean+RMS*5); // for drawing
            if (mSaveHists) drawFitting(h,f,canv,Form("%s, fit in (%.2e,%.2e), mean=%.2e+-%.2e, #sigma=%.2e, #chi^{2}=%.2e",drawTitle.Data(),tempLeft,tempRight,tempMean,tempError,tempSig,tempChi2),Form("%s_LandauFlipped_%s.layer%d.png",h->GetName(),mRunName.Data(),mLayerID),tempLeft,tempMean,tempRight);
            if (tempChi2<chi2){// The Landau function is better
                MyNamedInfo("XTBinAnalyzer",Form("        prefer Landau over Gaussian: chi2 %.2e -> %.2e, mean %.2e+-%.2e -> %.2e+-%.2e",chi2,tempChi2,mean,error,tempMean,tempError));
                if (previousFunction == kGaussian){ // previous function is Gaussian, then mix Gaussian with Landau
                    MyNamedInfo("XTBinAnalyzer",Form("        Gaussian to Landau: chi2 (%.2e+%.2e)/2=%.2e, mean (%.2e+%.2e)/2=%.2e",chi2,tempChi2,(chi2+tempChi2)/2,mean,tempMean,(mean+tempMean)/2));
                    mean = (mean+tempMean)/2; error = (error+tempError)/2; sigma = (sigma+tempSig)/2; chi2 = (chi2+tempChi2)/2;
                    if (tempFunction==kLandauFlipped) mFunction = kGaussianPlusLandauFlipped;
                    else if (tempFunction==kLandau) mFunction = kGaussianPlusLandau;
                }
                else{
                    mean = tempMean; error = tempError; sigma = tempSig;
                    function = tempFunction;
                }
                chi2 = tempChi2; left = tempLeft; right = tempRight;
            }
            else{ // The Landau function is worse
                if (previousFunction == kLandau || previousFunction == kLandauFlipped){ // previous function is Landau, then mix Gaussian with Landau
                    MyNamedInfo("XTBinAnalyzer",Form("        Landau to Gaussian: chi2 (%.2e+%.2e)/2=%.2e, mean (%.2e+%.2e)/2=%.2e",chi2,tempChi2,(chi2+tempChi2)/2,mean,tempMean,(mean+tempMean)/2));
                    mean = (mean+tempMean)/2; error = (error+tempError)/2; sigma = (sigma+tempSig)/2; chi2 = (chi2+tempChi2)/2;
                    if (tempFunction==kLandauFlipped) mFunction = kGaussianPlusLandauFlipped;
                    else if (tempFunction==kLandau) mFunction = kGaussianPlusLandau;
                }
            }
        }
        previousFunction = mFunction;
    }
    MyNamedInfo("XTBinAnalyzer",Form("      after fitSlice: mean=%.2e+-%.2e, fit range %.2e~%.2e, sig=%.2e, chi2=%.3e",mean,error,left,right,sigma,chi2));
}

void XTBinAnalyzer::drawFitting(TH1D* h,TF1 * f, TCanvas * c,TString title, TString filename,double left, double center, double right){
    if (!h) fprintf(stderr,"ERROR: in drawFitting, input histogram does not exist!\n");
    if (!f) fprintf(stderr,"ERROR: in drawFitting, input function does not exist!\n");
    if (!c) fprintf(stderr,"ERROR: in drawFitting, input canvas does not exist!\n");
    if (!h||!f||!c) return;
    c->cd();
    int oldStyle = gStyle->GetOptStat();
    gStyle->SetOptStat(1);
    h->SetTitle(title);
    h->Draw();
    double max = h->GetMaximum();
    TLine * l_left = new TLine(left,0,left,max);
    TLine * l_center= new TLine(center,0,center,max);
    TLine * l_right = new TLine(right,0,right,max);
    l_left->Draw();
    l_center->Draw();
    l_right->Draw();
    f->Draw("SAME");
    c->SaveAs(filename);
    gStyle->SetOptStat(oldStyle);
    delete l_left;
    delete l_center;
    delete l_right;
}

void XTBinAnalyzer::getRange(TH1D * h, int & binm, int & binl, int & binr, double ratio){
    // find the highest bin in the given region
    // NOTICE: by default ratio is 0.3
    // NOTICE: always using 2 adjacent bins to help smoothing
    binm = 1;
    double max = -1e14;
    int nBins = h->GetNbinsX();
    for (int i = 1+1;i<=nBins-1; i++){
        double value = (h->GetBinContent(i-1)+h->GetBinContent(i)+h->GetBinContent(i+1))/3;
        if (max<value){
            max = value;
            binm = i;
        }
    }
    // look for the left edge
    binl = binm-1;
    for (;binl>=2; binl--){
        double height3bins = h->GetBinContent(binl);
        height3bins+=h->GetBinContent(binl-1);
        height3bins+=h->GetBinContent(binl+1);
        if (height3bins/3<max*ratio) break;
    }
    // look for the right edge
    binr = binm+1;
    for (;binr<=nBins-1; binr++){
        double height3bins = h->GetBinContent(binr);
        height3bins+=h->GetBinContent(binr+1);
        height3bins+=h->GetBinContent(binr-1);
        if (height3bins/3<max*ratio) break;
    }
    MyNamedVerbose("XTBinAnalyzer",Form("        in getRange: peak height %.3f peak bin %d left bin %d right bin %d",max,binm,binl,binr));
}

void XTBinAnalyzer::fitSliceFloat(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right, double ratio){
    // get the new range according to the given threshold
    int binm, binl, binr;
    getRange(h,binm,binl,binr,ratio);
    // set left and right if the new values are stricter
    double l = h->GetBinCenter(binl);
    double r = h->GetBinCenter(binr);
    double m = h->GetBinCenter(binm);
    MyNamedVerbose("XTBinAnalyzer",Form("        in fitSliceFloat: old range %.3f ~ %.3f, new range w.r.t ratio=%.3f %.3f ~ %.3f, mean position %.3f",left,right,ratio,l,r,m));
    // set the new range and get the according statistics
    left = l;
    right = r;
    h->GetXaxis()->SetRangeUser(left,right);
    // TODO: should we use the average value or peak value?
    //mean = m;
    mean = h->GetMean();
    // add sigma, error and chi2
    sigma = h->GetRMS();
    meanErr = sigma;
    chi2 = 0; // 0 chi2 means this slice is analyzed statistically, not from fitting
}

TF1 * XTBinAnalyzer::fitSliceGaus(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right){
    // set the range according to the initial range
    h->GetXaxis()->SetRangeUser(left,right);
    // prepare the function
    TF1 * f = f_gaus;
    // TODO: should we update the mean value as initial value or use the mean value given by the user?
    mean = h->GetMean();
    f->SetParameter(0,h->GetMaximum());
    f->SetParameter(1,mean);
    f->SetParameter(2,h->GetRMS());
    MyNamedVerbose("XTBinAnalyzer",Form("        in fitSliceGaus: range %.3f ~ %.3f, mean position %.3f, initial height %.3f, initial width %.3f",left,right,mean,h->GetMaximum(),h->GetRMS()));
    // fit within the new range
    double lrange = mean-left;
    double rrange = right-mean;
    for (int i = 0; i<10; i++){
        h->Fit("fgaus","qN0","",left,right);
        mean = f->GetParameter(1);
        sigma = fabs(f->GetParameter(2));
        chi2 = f->GetChisquare();
        MyNamedVerbose("XTBinAnalyzer",Form("          fit[%d]: height %.3f sigma %.3f mean %.3f chi2 %.3f",i,f->GetParameter(0),sigma,mean,chi2));
        if (sigma<rrange*5&&sigma>rrange/5&&mean>left&&mean<right) break;
        left-=lrange/10;
        right+=rrange/10;
        h->GetXaxis()->SetRangeUser(left,right);
        MyNamedVerbose("XTBinAnalyzer",Form("          range -> %.3f -- %.3f",left,right));
    }
    meanErr = f->GetParError(1);
    if (mean>right||mean<left){
        sigma = 1e9;
        meanErr = sigma;
    }
    return f;
}

TF1 * XTBinAnalyzer::fitSliceLand(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right){
    h->GetXaxis()->SetRangeUser(left,right);
    // prepare the function
    TF1 * f = f_land;
    mean = h->GetMean();
    f->SetParameter(0,h->GetMaximum()/0.18065564);
    f->SetParameter(1,mean);
    f->SetParameter(2,h->GetRMS());
    MyNamedVerbose("XTBinAnalyzer",Form("        in fitSliceLand: range %.3f ~ %.3f, mean position %.3f, initial height %.3f, initial width %.3f",left,right,mean,h->GetMaximum(),h->GetRMS()));
    // fit within the new range
    double lrange = mean-left;
    double rrange = right-mean;
    for (int i = 0; i<10; i++){
        h->Fit("fland","qN0","",left,right);
        mean = f->GetParameter(1);
        sigma = fabs(f->GetParameter(2));
        chi2 = f->GetChisquare();
        MyNamedVerbose("XTBinAnalyzer",Form("          fit[%d]: height %.3f sigma %.3f mean %.3f chi2 %.3f",i,f->GetParameter(0),sigma,mean,chi2));
        if (sigma<rrange*5&&sigma>rrange/5&&mean>left&&mean<right) break;
        left-=lrange/10;
        right+=rrange/10;
        h->GetXaxis()->SetRangeUser(left,right);
        MyNamedVerbose("XTBinAnalyzer",Form("          range -> %.3f -- %.3f",left,right));
    }
    meanErr = f->GetParError(1);
    if (mean>right||mean<left){
        sigma = 1e9;
        meanErr = sigma;
    }
    return f;
}

TF1 * XTBinAnalyzer::fitSliceLandF(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right){
    h->GetXaxis()->SetRangeUser(left,right);
    // prepare the function
    TF1 * f = f_landF;
    mean = h->GetMean();
    f->SetParameter(0,h->GetMaximum()/0.18065564);
    f->SetParameter(1,-mean);
    f->SetParameter(2,h->GetRMS());
    MyNamedVerbose("XTBinAnalyzer",Form("        in fitSliceLandF: range %.3f ~ %.3f, mean position %.3f, initial height %.3f, initial width %.3f",left,right,mean,h->GetMaximum(),h->GetRMS()));
    // fit within the new range
    double lrange = mean-left;
    double rrange = right-mean;
    for (int i = 0; i<10; i++){
        h->Fit("flandF","qN0","",left,right);
        mean = -f->GetParameter(1); // flipped
        sigma = fabs(f->GetParameter(2));
        chi2 = f->GetChisquare();
        MyNamedVerbose("XTBinAnalyzer",Form("          fit[%d]: height %.3f sigma %.3f mean %.3f chi2 %.3f",i,f->GetParameter(0),sigma,mean,chi2));
        if (sigma<rrange*5&&sigma>rrange/5&&mean>left&&mean<right) break;
        left-=lrange/10;
        right+=rrange/10;
        h->GetXaxis()->SetRangeUser(left,right);
        MyNamedVerbose("XTBinAnalyzer",Form("          range -> %.3f -- %.3f",left,right));
    }
    meanErr = f->GetParError(1);
    if (mean>right||mean<left){
        sigma = 1e9;
        meanErr = sigma;
    }
    return f;
}

TF1 * XTBinAnalyzer::fitSlice2Gaus(TH1D * h, double & mean, double & meanErr, double & sigma, double & chi2, double & left, double & right){
    h->GetXaxis()->SetRangeUser(left,right);
    // prepare the function
    TF1 * f = f_2gaus;
    mean = h->GetMean();
    f->SetParameter(0,h->GetMaximum());
    f->SetParameter(1,mean);
    f->SetParameter(2,h->GetRMS()/2);
    f->SetParameter(3,h->GetRMS()*2);
    MyNamedVerbose("XTBinAnalyzer",Form("        in fitSlice2Gaus: range %.3f ~ %.3f, mean position %.3f, initial height %.3f, initial width %.3f %.3f",left,right,mean,h->GetMaximum(),h->GetRMS()/2,h->GetRMS()*2));
    // fit within the new range
    double lrange = mean-left;
    double rrange = right-mean;
    double sigma1;
    double sigma2;
    for (int i = 0; i<10; i++){
        h->Fit("f2gaus","qN0","",left,right);
        mean = f->GetParameter(1);
        sigma1 = fabs(f->GetParameter(2));
        sigma2 = fabs(f->GetParameter(3));
        chi2 = f->GetChisquare();
        MyNamedVerbose("XTBinAnalyzer",Form("          fit[%d]: sigma %.3f mean %.3f chi2 %.3f",i,sigma,mean,chi2));
        if (sigma1<rrange*5&&sigma2<lrange*5&&sigma1>rrange/5&&sigma2>lrange/5&&mean>left&&mean<right) break;
        left-=lrange/10;
        right+=rrange/10;
        h->GetXaxis()->SetRangeUser(left,right);
        MyNamedVerbose("XTBinAnalyzer",Form("          range -> %.3f -- %.3f",left,right));
    }
    // set results
    sigma = sqrt((sigma1*sigma1+sigma2*sigma2)/2);
    meanErr = f->GetParError(1);
    if (mean>right||mean<left){
        sigma = 1e9;
        meanErr = sigma;
    }
    return f;
}
