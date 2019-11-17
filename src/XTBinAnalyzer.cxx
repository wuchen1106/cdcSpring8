#include "XTBinAnalyzer.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TStyle.h"
#include "TGraphErrors.h"

#include <math.h>

#include "ParameterManager.hxx"
#include "XTManager.hxx"
#include "Log.hxx"

XTBinAnalyzer::XTBinAnalyzer(TString runname, TFile * outfile, bool drawDetails)
        :mRunName(runname), mDrawDetails(drawDetails),
        mOutFile(outfile), mOutTree(NULL), mSliceDir(NULL),
        mTestLayerID(0),
        mLayerID(0), mCellID(0), mX(0), mXerr(0), mT(0), mTerr(0), mSig(0), mChi2(0), mEntries(0), mType(0), mFunction(0), 
        h2_xt(NULL), h2_xtn(NULL),
        f_gaus(NULL), f_land(NULL), f_landF(NULL), f_2gaus(NULL),
        f_left_cen(NULL), f_left_mid(NULL), f_left_end(NULL),
        f_right_cen(NULL), f_right_mid(NULL), f_right_end(NULL),
        f_folded_cen(NULL), f_folded_mid(NULL), f_folded_end(NULL),
        f_left(NULL), f_right(NULL), f_folded(NULL)
{
    // prepare functions for slice analysis
    f_gaus = myNewTF1("fgaus","gaus",-1000,1000);
    f_land = myNewTF1("fland","landau",-1000,1000);
    f_landF = myNewTF1("flandF","landau(-x)",-1000,1000);
    f_2gaus = myNewTF1("f2gaus","(x>[1])*[0]*exp(-0.5*((x-[1])/[2])**2)+(x<=[1])*[0]*exp(-0.5*((x-[1])/[3])**2)",-1000,1000);
}

XTBinAnalyzer::~XTBinAnalyzer(void){
    if (f_gaus) delete f_gaus;
    if (f_land) delete f_land;
    if (f_landF) delete f_landF;
    if (f_2gaus) delete f_2gaus;
}

int XTBinAnalyzer::Initialize(bool reLoad){
    // about binning
    double bin_t_min = ParameterManager::Get().XTAnalyzerParameters.bin_t_min; // t range for one x bin
    double bin_t_max = ParameterManager::Get().XTAnalyzerParameters.bin_t_max;
    int    bin_t_num = ParameterManager::Get().XTAnalyzerParameters.bin_t_num;
    double bin_x_min = ParameterManager::Get().XTAnalyzerParameters.bin_x_min;
    double bin_x_max = ParameterManager::Get().XTAnalyzerParameters.bin_x_max; // x range for one t bin
    int    bin_x_num = ParameterManager::Get().XTAnalyzerParameters.bin_x_num;
    if (reLoad){
        if (!mOutTree){
            mOutTree = (TTree*) mOutFile->Get("XTBins");
            mOutTree->SetBranchAddress("x",&mX);
            mOutTree->SetBranchAddress("t",&mT);
            mOutTree->SetBranchAddress("xerr",&mXerr);
            mOutTree->SetBranchAddress("terr",&mTerr);
            mOutTree->SetBranchAddress("lid",&mLayerID);
            mOutTree->SetBranchAddress("wid",&mCellID);
            mOutTree->SetBranchAddress("sig",&mSig);
            mOutTree->SetBranchAddress("chi2",&mChi2);
            mOutTree->SetBranchAddress("n",&mEntries);
            mOutTree->SetBranchAddress("type",&mType);
            mOutTree->SetBranchAddress("func",&mFunction);
        }
        if (!mOutFile){
            MyError("Cannot load XTBins from the file \""<<mOutFile->GetPath()<<"\"");
            return 1;
        }
        mSliceDir = (TDirectory*) mOutFile->Get(Form("slices_%d",mTestLayerID));
        if (!mSliceDir){
            MyError("Cannot load directory \""<<Form("slices_%d",mTestLayerID)<<"\" from the file \""<<mOutFile->GetPath()<<"\"");
            return 1;
        }
        h2_xt = (TH2D*)mOutFile->Get(Form("h2_xt_%d",mTestLayerID));
        h2_xtn = (TH2D*)mOutFile->Get(Form("h2_xtn_%d",mTestLayerID));
        if (!h2_xt||!h2_xtn){
            MyError("Cannot load \""<<Form("h2_xt_%d",mTestLayerID)<<"\" or \""<<Form("h2_xtn_%d",mTestLayerID)<<"\" from the file \""<<mOutFile->GetPath()<<"\"");
            return 1;
        }
    }
    else{
        if (!mOutTree){
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
        mOutFile->cd();
        mSliceDir = mOutFile->mkdir(Form("slices_%d",mTestLayerID));

        // prepare 2D histograms
        //if (h2_xt) delete h2_xt;
        //if (h2_xtn) delete h2_xtn;
        mOutFile->cd();
        h2_xt = new TH2D(Form("h2_xt_%d",mTestLayerID),"XT Relation",bin_t_num,bin_t_min,bin_t_max,bin_x_num,-bin_x_max,bin_x_max);
        h2_xt->GetXaxis()->SetTitle("T [ns]");
        h2_xt->GetYaxis()->SetTitle("X [mm]");
        h2_xtn = new TH2D(Form("h2_xtn_%d",mTestLayerID),"XT Relation",bin_t_num,bin_t_min,bin_t_max,bin_x_num/2+1,bin_x_min,bin_x_max);
        h2_xtn->GetXaxis()->SetTitle("T [ns]");
        h2_xtn->GetYaxis()->SetTitle("X [mm]");
    }

    // WARNING: this Initialize will only create these XT functions ONCE in XTBinAnalyzer's lifetime!!!
    int mCentPolN = ParameterManager::Get().XTAnalyzerParameters.xt_center_nPol;
    int mMidPolN = ParameterManager::Get().XTAnalyzerParameters.xt_middle_nPol;
    int mEndPolN = ParameterManager::Get().XTAnalyzerParameters.xt_end_nPol;
    if (!f_left_cen) f_left_cen = myNewTF1(Form("flce_%d",mLayerID),Form("pol%d",mCentPolN),bin_t_min,bin_t_max);
    if (!f_right_cen) f_right_cen = myNewTF1(Form("frce_%d",mLayerID),Form("pol%d",mCentPolN),bin_t_min,bin_t_max);
    if (!f_folded_cen) f_folded_cen = myNewTF1(Form("ffce_%d",mLayerID),Form("pol%d",mCentPolN),bin_t_min,bin_t_max);
    if (!f_left_mid) f_left_mid = myNewTF1(Form("flm_%d",mLayerID),Form("pol%d",mMidPolN),bin_t_min,bin_t_max);
    if (!f_right_mid) f_right_mid = myNewTF1(Form("frm_%d",mLayerID),Form("pol%d",mMidPolN),bin_t_min,bin_t_max);
    if (!f_folded_mid) f_folded_mid = myNewTF1(Form("ffm_%d",mLayerID),Form("pol%d",mMidPolN),bin_t_min,bin_t_max);
    if (!f_left_end) f_left_end = myNewTF1(Form("fle_%d",mLayerID),Form("pol%d",mEndPolN),bin_t_min,bin_t_max);
    if (!f_right_end) f_right_end = myNewTF1(Form("fre_%d",mLayerID),Form("pol%d",mEndPolN),bin_t_min,bin_t_max);
    if (!f_folded_end) f_folded_end = myNewTF1(Form("ffe_%d",mLayerID),Form("pol%d",mEndPolN),bin_t_min,bin_t_max);

    MyNamedLog("XTBinAnalyzer","XTBinAnalyzer successfully initialized!");
    return 0;
}

void XTBinAnalyzer::Fill(double t, double x){
    double absx = fabs(x);
    h2_xt->Fill(t,x);
    h2_xtn->Fill(t,absx);
}

void XTBinAnalyzer::BinAnalysis(void){
    MyNamedLog("XTBinAnalyzer","In XTBinAnalyzer::Process");
    //==========================Load the parameters first==============================
    // about projection
    int    bin_t_fit_num = ParameterManager::Get().XTAnalyzerParameters.bin_t_fit_num;
    double bin_t_tailTime = ParameterManager::Get().XTAnalyzerParameters.bin_t_tailTime;
    int    bin_t_fit_num_tail = ParameterManager::Get().XTAnalyzerParameters.bin_t_fit_num_tail;
    int    bin_x_fit_num = ParameterManager::Get().XTAnalyzerParameters.bin_x_fit_num;
    // about the range for using Landau function
    //   time slice to fit space
    double bin_t_landTmin = ParameterManager::Get().XTAnalyzerParameters.bin_t_landTmin;
    double bin_t_landTmax = ParameterManager::Get().XTAnalyzerParameters.bin_t_landTmax;
    //   space slice to fit time
    double bin_x_landXmin = ParameterManager::Get().XTAnalyzerParameters.bin_x_landXmin;
    double bin_x_landXmax = ParameterManager::Get().XTAnalyzerParameters.bin_x_landXmax;
    // about fitting range
    double bin_t_ratio = ParameterManager::Get().XTAnalyzerParameters.bin_t_ratio;
    double bin_x_ratio = ParameterManager::Get().XTAnalyzerParameters.bin_x_ratio;

    //==========================Preparation work==============================
    // prepare a canvas to draw fittings if needed
    TCanvas * canv = new TCanvas("cfit","cfit",1024,768);canv->SetGridx(1);canv->SetGridy(1);
    // record the current style
    int oldStyle = gStyle->GetOptStat();
    gStyle->SetOptStat(1); // we want to draw the stat box for bin-by-bin histograms
    // go to the output file and save the sampled histograms
    mOutFile->cd();
    h2_xt->Write();
    h2_xtn->Write();
    // enter the directory to hold the slices
    mSliceDir->cd();
    int current_nbins_to_fit = 1;
    int previousFunction = kGaussian;
    // set the layer ID for this process
    mLayerID = mTestLayerID;
    TString drawTitle;

    //==========================Fit T along X slices (unfolded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kSpaceSliceUnfolded; // set the type of this round of fittings
    previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
    current_nbins_to_fit = bin_x_fit_num;
    for (int iLR = 0; iLR<2; iLR++){
        for (int i = (iLR==0?h2_xt->GetYaxis()->FindBin(0.):h2_xt->GetYaxis()->FindBin(0.)+1); i>=1&&(iLR==0?true:i<=h2_xt->GetYaxis()->GetNbins()); i+=(iLR==0?-current_nbins_to_fit:current_nbins_to_fit)){
            // fist, get the time histogram for this x slice
            double divleft = h2_xt->GetYaxis()->GetBinLowEdge(i);
            double divmiddle = current_nbins_to_fit%2==0?h2_xt->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xt->GetYaxis()->GetBinCenter(i+current_nbins_to_fit/2);
            double divright = h2_xt->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit);
            TH1D * h = h2_xt->ProjectionX(Form("h_x_%d_%s_%d",mTestLayerID,iLR==0?"L":"R",i),i,i+current_nbins_to_fit-1);
            h->Write();
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
            fitSlice(h,canv,drawTitle,mT,mTerr,mSig,mChi2,left,right,mFunction,previousFunction,bin_t_ratio,fabs(mX)>=bin_x_landXmin&&fabs(mX)<=bin_x_landXmax,false);
            // at last, fill the tree
            mOutTree->Fill();
        }
    }

    //==========================Fit T along X slices (folded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kSpaceSliceFolded; // set the type of this round of fittings
    previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
    current_nbins_to_fit = bin_x_fit_num;
    for (int i = 1; i<=h2_xtn->GetYaxis()->GetNbins(); i+=current_nbins_to_fit){
        // fist, get the time histogram for this x slice
        double divleft = h2_xtn->GetYaxis()->GetBinLowEdge(i);
        double divmiddle = current_nbins_to_fit%2==0?h2_xtn->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xtn->GetYaxis()->GetBinCenter(i+current_nbins_to_fit/2);
        double divright = h2_xtn->GetYaxis()->GetBinLowEdge(i+current_nbins_to_fit);
        TH1D * h = h2_xtn->ProjectionX(Form("h_x_%d_F_%d",mTestLayerID,i),i,i+current_nbins_to_fit-1);
        h->Write();
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
        fitSlice(h,canv,drawTitle,mT,mTerr,mSig,mChi2,left,right,mFunction,previousFunction,bin_t_ratio,fabs(mX)>=bin_x_landXmin&&fabs(mX)<=bin_x_landXmax,false);
        // at last, fill the tree
        mOutTree->Fill();
    }

    //==========================Fit X along T slices (unfolded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kTimeSliceUnfolded; // set the type of this round of fittings
    double x_max = h2_xt->GetYaxis()->GetXmax();
    double x_min = h2_xt->GetYaxis()->GetXmin();
    for (int iLR = 0; iLR<2; iLR++){
        if (iLR == 0){ // left side
            h2_xt->GetYaxis()->SetRangeUser(x_min,0);
        }
        else{ // right side
            h2_xt->GetYaxis()->SetRangeUser(0,x_max);
        }
        previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
        current_nbins_to_fit = bin_t_fit_num;
        for (int i = 1; i<=h2_xt->GetXaxis()->GetNbins(); i+=current_nbins_to_fit){
            // fist, get the x histogram for this time slice
            double divleft = h2_xt->GetXaxis()->GetBinLowEdge(i);
            double divmiddle = current_nbins_to_fit%2==0?h2_xt->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xt->GetXaxis()->GetBinCenter(i+current_nbins_to_fit/2);
            double divright = h2_xt->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit);
            if (divright>=bin_t_tailTime) current_nbins_to_fit = bin_t_fit_num_tail; // switch to tail mode;
            TH1D * h = h2_xt->ProjectionY(Form("h_t_%d_%s_%d",mTestLayerID,iLR==0?"L":"R",i),i,i+current_nbins_to_fit-1);
            h->Write();
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
            fitSlice(h,canv,drawTitle,mX,mXerr,mSig,mChi2,left,right,mFunction,previousFunction,bin_t_ratio,mT>=bin_t_landTmin&&mT<=bin_t_landTmax,mX>0);
            // correct the X value if it passes the range
            if (iLR==0&&mX>0) mX = 0;
            else if (iLR==1&&mX<0) mX = 0;
            // at last, fill the tree
            mOutTree->Fill();
        }
    }
    // set the 2-D histogram back after this analysis
    h2_xt->GetYaxis()->UnZoom();

    //==========================Fit X along T slices (folded)==============================
    MyNamedInfo("XTBinAnalyzer","Fit time slice (unfolded)");
    mType = kTimeSliceFolded; // set the type of this round of fittings
    previousFunction = kGaussian;// record the previous used function to mark the function switching point so that the turning point can be smoothed out.
    current_nbins_to_fit = bin_t_fit_num;
    for (int i = 1; i<=h2_xtn->GetXaxis()->GetNbins(); i+=current_nbins_to_fit){
        // fist, get the x histogram for this time slice
        double divleft = h2_xtn->GetXaxis()->GetBinLowEdge(i);
        double divmiddle = current_nbins_to_fit%2==0?h2_xtn->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit/2):h2_xtn->GetXaxis()->GetBinCenter(i+current_nbins_to_fit/2);
        double divright = h2_xtn->GetXaxis()->GetBinLowEdge(i+current_nbins_to_fit);
        if (divright>=bin_t_tailTime) current_nbins_to_fit = bin_t_fit_num_tail; // switch to tail mode;
        TH1D * h = h2_xtn->ProjectionY(Form("h_t_%d_F_%d",mTestLayerID,i),i,i+current_nbins_to_fit-1);
        h->Write();
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
        fitSlice(h,canv,drawTitle,mX,mXerr,mSig,mChi2,left,right,mFunction,previousFunction,bin_t_ratio,mT>=bin_t_landTmin&&mT<=bin_t_landTmax,true);
        // at last, fill the tree
        mOutTree->Fill();
    }

    //==========================Form the graphs==============================
    gStyle->SetOptStat(0);
    drawSamples();
    formXTGraphs();

    //==========================Fit the XT relation functions==============================
    FitXT();
    f_left_cen->Write();
    f_right_cen->Write();
    f_folded_cen->Write();
    f_left_mid->Write();
    f_right_mid->Write();
    f_folded_mid->Write();
    f_left_end->Write();
    f_right_end->Write();
    f_folded_end->Write();

    // set the style back
    gStyle->SetOptStat(oldStyle);
}

void XTBinAnalyzer::FitXT(){
}

void XTBinAnalyzer::Write(){
    mOutFile->cd();
    mOutTree->Write();
}

double XTBinAnalyzer::interpolate(const TGraphErrors * graph, double theX){
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
        if (isLeft&&x<theX||!isLeft&&x>theX){ // moved to the other side of the asked position, then interpolate
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

TGraphErrors * XTBinAnalyzer::minusGraph(const TGraphErrors * gr_left, const TGraphErrors * gr_right, TString name, bool isNegative){
    TGraphErrors * gr_new = new TGraphErrors(); gr_new->SetName(name);
    gr_new->Set(gr_left->GetN());
    for (int iPoint = 0; iPoint<gr_left->GetN(); iPoint++){
        double x,y; double xerr,yerr;
        gr_left->GetPoint(iPoint,x,y);
        if (isNegative) y*=-1;
        xerr = gr_left->GetErrorX(iPoint);
        yerr = gr_left->GetErrorY(iPoint);
        double y0 = interpolate(gr_right,x);
        gr_new->SetPoint(iPoint,x,y-y0);
        gr_new->SetPointError(iPoint,xerr,yerr);
    }
    return gr_new;
}

void XTBinAnalyzer::drawSamples(){
    mOutTree->SetMarkerStyle(20);mOutTree->SetMarkerSize(0.5);
    TCanvas * canv= new TCanvas("csample","csample",1024,768);
    // firstly, get the maximum number of entries in one sample
    int mEntriesXMax = 0;
    int mEntriesTMax = 0;
    for (Long64_t iEntry = 0; iEntry<mOutTree->GetEntries(); iEntry++){
        mOutTree->GetEntry(iEntry);
        if (mLayerID!=mTestLayerID) continue;
        if ((mType==kSpaceSliceUnfolded||mType==kSpaceSliceFolded)&&mEntries>mEntriesXMax) mEntriesXMax = mEntries;
        if ((mType==kTimeSliceUnfolded||mType==kTimeSliceFolded)&&mEntries>mEntriesTMax) mEntriesTMax = mEntries;
    }
    canv->Divide(3,2);
    canv->cd(1);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_entriesX = new TH2D("h2_bkg_entriesX","Number of entries in each X slice",512,-10,10,1024,0,mEntriesXMax*1.1);
    h2_bkg_entriesX->Draw();
    mOutTree->SetMarkerColor(kBlack); mOutTree->Draw("n:x",Form("type==%d",kSpaceSliceFolded),"PSAME");
    mOutTree->SetMarkerColor(kMagenta); mOutTree->Draw("n:x",Form("type==%d",kSpaceSliceUnfolded),"PSAME");
    canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_sigX = new TH2D("h2_bkg_sigX","#sigma of T fitting in each X slice",512,-10,10,512,0,20);
    h2_bkg_sigX->Draw();
    mOutTree->SetMarkerColor(kBlack); mOutTree->Draw("sig:x",Form("type==%d",kSpaceSliceFolded),"PSAME");
    mOutTree->SetMarkerColor(kMagenta); mOutTree->Draw("sig:x",Form("type==%d",kSpaceSliceUnfolded),"PSAME");
    canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_chi2X = new TH2D("h2_bkg_chi2X","#chi^{2} of T fitting in each X slice",512,-10,10,512,0,50);
    h2_bkg_chi2X->Draw();
    mOutTree->SetMarkerColor(kBlack); mOutTree->Draw("chi2:x",Form("type==%d",kSpaceSliceFolded),"PSAME");
    mOutTree->SetMarkerColor(kMagenta); mOutTree->Draw("chi2:x",Form("type==%d",kSpaceSliceUnfolded),"PSAME");
    canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_entriesT = new TH2D("h2_bkg_entriesT","Number of entries in each T slice",512,-25,800,1024,0,mEntriesTMax*1.1);
    h2_bkg_entriesT->Draw();
    mOutTree->SetMarkerColor(kBlack); mOutTree->Draw("n:t",Form("type==%d",kTimeSliceFolded),"PSAME");
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("n:t",Form("type==%d&&x>=0",kTimeSliceUnfolded),"PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("n:t",Form("type==%d&&x<0",kTimeSliceUnfolded),"PSAME");
    canv->cd(5);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_sigT = new TH2D("h2_bkg_sigT","#sigma of X fitting in each T slice",512,-25,800,512,0,0.5);
    h2_bkg_sigT->Draw();
    mOutTree->SetMarkerColor(kBlack); mOutTree->Draw("sig:t",Form("type==%d",kTimeSliceFolded),"PSAME");
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("sig:t",Form("type==%d&&x>=0",kTimeSliceUnfolded),"PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("sig:t",Form("type==%d&&x<0",kTimeSliceUnfolded),"PSAME");
    canv->cd(6);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_chi2T = new TH2D("h2_bkg_chi2T","#chi^{2} of X fitting in each T slice",512,-25,800,512,0,100);
    h2_bkg_chi2T->Draw();
    mOutTree->SetMarkerColor(kBlack); mOutTree->Draw("chi2:t",Form("type==%d",kTimeSliceFolded),"PSAME");
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("chi2:t",Form("type==%d&&x>=0",kTimeSliceUnfolded),"PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("chi2:t",Form("type==%d&&x<0",kTimeSliceUnfolded),"PSAME");
    canv->SaveAs(Form("result/samples_%s.layer%d.png",mRunName.Data(),mTestLayerID));
}

void XTBinAnalyzer::formXTGraphs(){
    // set parameters
    double bin_x_min = ParameterManager::Get().XTAnalyzerParameters.bin_x_min;
    int    graph_n_min = ParameterManager::Get().XTAnalyzerParameters.graph_n_min;
    double graph_chi2_max = ParameterManager::Get().XTAnalyzerParameters.graph_chi2_max;
    double graph_sepX = ParameterManager::Get().XTAnalyzerParameters.graph_sepX;
    double markerSize = 0.5;

    // make graphs from different samplings: left/right/folded TIMES time/space
    TGraphErrors * gr_timeSlices_unfolded_left = new TGraphErrors(); gr_timeSlices_unfolded_left->SetName(Form("gr_%d_t_l",mTestLayerID));
    gr_timeSlices_unfolded_left->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x<%.7e",kTimeSliceUnfolded,graph_n_min,graph_chi2_max,0.)));
    gr_timeSlices_unfolded_left->SetMarkerStyle(24);gr_timeSlices_unfolded_left->SetMarkerSize(markerSize);
    gr_timeSlices_unfolded_left->SetMarkerColor(kRed);gr_timeSlices_unfolded_left->SetLineColor(kRed);
    TGraphErrors * gr_timeSlices_unfolded_right = new TGraphErrors(); gr_timeSlices_unfolded_right->SetName(Form("gr_%d_t_r",mTestLayerID));
    gr_timeSlices_unfolded_right->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x>%.7e",kTimeSliceUnfolded,graph_n_min,graph_chi2_max,0.)));
    gr_timeSlices_unfolded_right->SetMarkerStyle(24);gr_timeSlices_unfolded_right->SetMarkerSize(markerSize);
    gr_timeSlices_unfolded_right->SetMarkerColor(kRed);gr_timeSlices_unfolded_right->SetLineColor(kRed);
    TGraphErrors * gr_spaceSlices_unfolded_left = new TGraphErrors(); gr_spaceSlices_unfolded_left->SetName(Form("gr_%d_x_l",mTestLayerID));
    gr_spaceSlices_unfolded_left->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x<=%.7e",kSpaceSliceUnfolded,graph_n_min,graph_chi2_max,-bin_x_min)));
    gr_spaceSlices_unfolded_left->SetMarkerStyle(20);gr_spaceSlices_unfolded_left->SetMarkerSize(markerSize);
    gr_spaceSlices_unfolded_left->SetMarkerColor(kBlack);gr_spaceSlices_unfolded_left->SetLineColor(kBlack);
    TGraphErrors * gr_spaceSlices_unfolded_right = new TGraphErrors(); gr_spaceSlices_unfolded_right->SetName(Form("gr_%d_x_r",mTestLayerID));
    gr_spaceSlices_unfolded_right->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x>=%.7e",kSpaceSliceUnfolded,graph_n_min,graph_chi2_max,bin_x_min)));
    gr_spaceSlices_unfolded_right->SetMarkerStyle(20);gr_spaceSlices_unfolded_right->SetMarkerSize(markerSize);
    gr_spaceSlices_unfolded_right->SetMarkerColor(kBlack);gr_spaceSlices_unfolded_right->SetLineColor(kBlack);
    TGraphErrors * gr_timeSlices_folded = new TGraphErrors(); gr_timeSlices_folded->SetName(Form("gr_%d_t_f",mTestLayerID));
    gr_timeSlices_folded->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0",kTimeSliceFolded,graph_n_min,graph_chi2_max)));
    gr_timeSlices_folded->SetMarkerStyle(24);gr_timeSlices_folded->SetMarkerSize(markerSize);
    gr_timeSlices_folded->SetMarkerColor(kRed);gr_timeSlices_folded->SetLineColor(kRed);
    TGraphErrors * gr_spaceSlices_folded = new TGraphErrors(); gr_spaceSlices_folded->SetName(Form("gr_%d_x_f",mTestLayerID));
    gr_spaceSlices_folded->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0",kSpaceSliceFolded,graph_n_min,graph_chi2_max)));
    gr_spaceSlices_folded->SetMarkerStyle(20);gr_spaceSlices_folded->SetMarkerSize(markerSize);
    gr_spaceSlices_folded->SetMarkerColor(kBlack);gr_spaceSlices_folded->SetLineColor(kBlack);
    // combined graphs
    TGraphErrors * gr_combined_unfolded_left = new TGraphErrors(); gr_combined_unfolded_left->SetName(Form("gr_%d_c_l",mTestLayerID));
    gr_combined_unfolded_left->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x<%.7e",kTimeSliceUnfolded,graph_n_min,graph_chi2_max,-graph_sepX))
                                    +mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x<0&&x>%.7e",kSpaceSliceUnfolded,graph_n_min,graph_chi2_max,-graph_sepX)));
    TGraphErrors * gr_combined_unfolded_right = new TGraphErrors(); gr_combined_unfolded_right->SetName(Form("gr_%d_c_r",mTestLayerID));
    gr_combined_unfolded_right->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x>%.7e",kTimeSliceUnfolded,graph_n_min,graph_chi2_max,graph_sepX))
                                    +mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x>0&&x<%.7e",kSpaceSliceUnfolded,graph_n_min,graph_chi2_max,graph_sepX)));
    TGraphErrors * gr_combined_folded = new TGraphErrors(); gr_combined_folded->SetName(Form("gr_%d_c_f",mTestLayerID));
    gr_combined_folded->Set(mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x>%.7e",kTimeSliceFolded,graph_n_min,graph_chi2_max,graph_sepX))
                                    +mOutTree->GetEntries(Form("type==%d&&n>=%d&&chi2<=%.7e&&func!=0&&x<%.7e",kSpaceSliceFolded,graph_n_min,graph_chi2_max,graph_sepX)));
    int count_timeSlices_unfolded_left = 0;
    int count_timeSlices_unfolded_right = 0;
    int count_timeSlices_folded = 0;
    int count_spaceSlices_unfolded_left = 0;
    int count_spaceSlices_unfolded_right = 0;
    int count_spaceSlices_folded = 0;
    int count_combined_unfolded_left = 0;
    int count_combined_unfolded_right= 0;
    int count_combined_folded = 0;
    MyNamedVerbose("XTBinAnalyzer","Looping in the output tree again, "<<mOutTree->GetEntries()<<" entries");
    for (Long64_t iEntry = 0; iEntry<mOutTree->GetEntries(); iEntry++){
        mOutTree->GetEntry(iEntry);
        if (mLayerID!=mTestLayerID);
        if (mType==kTimeSliceUnfolded&&mEntries>=graph_n_min&&mChi2<=graph_chi2_max&&mFunction!=0){
            if (mX<0){
                printf("X %.14e < 0\n",mX);
                gr_timeSlices_unfolded_left->SetPoint(count_timeSlices_unfolded_left,mT,mX);
                gr_timeSlices_unfolded_left->SetPointError(count_timeSlices_unfolded_left,mTerr,mXerr);
                count_timeSlices_unfolded_left++;
            }
            if (mX>0){
                gr_timeSlices_unfolded_right->SetPoint(count_timeSlices_unfolded_right,mT,mX);
                gr_timeSlices_unfolded_right->SetPointError(count_timeSlices_unfolded_right,mTerr,mXerr);
                count_timeSlices_unfolded_right++;
            }
            if (mX>graph_sepX){
                gr_combined_unfolded_right->SetPoint(count_combined_unfolded_right,mT,mX);
                gr_combined_unfolded_right->SetPointError(count_combined_unfolded_right,mTerr,mXerr);
                count_combined_unfolded_right++;
            }
            if (mX<-graph_sepX){
                gr_combined_unfolded_left->SetPoint(count_combined_unfolded_left,mT,mX);
                gr_combined_unfolded_left->SetPointError(count_combined_unfolded_left,mTerr,mXerr);
                count_combined_unfolded_left++;
            }
        }
        else if (mType==kTimeSliceFolded&&mEntries>=graph_n_min&&mChi2<=graph_chi2_max&&mFunction!=0){
            gr_timeSlices_folded->SetPoint(count_timeSlices_folded,mT,mX);
            gr_timeSlices_folded->SetPointError(count_timeSlices_folded,mTerr,mXerr);
            count_timeSlices_folded++;
            if (mX>graph_sepX){
                gr_combined_folded->SetPoint(count_combined_folded,mT,mX);
                gr_combined_folded->SetPointError(count_combined_folded,mTerr,mXerr);
                count_combined_folded++;
            }
        }
        else if (mType==kSpaceSliceUnfolded&&mEntries>=graph_n_min&&mChi2<=graph_chi2_max&&mFunction!=0){
            if (mX<=-bin_x_min){
                gr_spaceSlices_unfolded_left->SetPoint(count_spaceSlices_unfolded_left,mT,mX);
                gr_spaceSlices_unfolded_left->SetPointError(count_spaceSlices_unfolded_left,mTerr,mXerr);
                count_spaceSlices_unfolded_left++;
                if (mX>-graph_sepX){
                    gr_combined_unfolded_left->SetPoint(count_combined_unfolded_left,mT,mX);
                    gr_combined_unfolded_left->SetPointError(count_combined_unfolded_left,mTerr,mXerr);
                    count_combined_unfolded_left++;
                }
            }
            if (mX>=bin_x_min){
                gr_spaceSlices_unfolded_right->SetPoint(count_spaceSlices_unfolded_right,mT,mX);
                gr_spaceSlices_unfolded_right->SetPointError(count_spaceSlices_unfolded_right,mTerr,mXerr);
                count_spaceSlices_unfolded_right++;
                if (mX<graph_sepX){
                    gr_combined_unfolded_right->SetPoint(count_combined_unfolded_right,mT,mX);
                    gr_combined_unfolded_right->SetPointError(count_combined_unfolded_right,mTerr,mXerr);
                    count_combined_unfolded_right++;
                }
            }
        }
        else if (mType==kSpaceSliceFolded&&mEntries>=graph_n_min&&mChi2<=graph_chi2_max&&mFunction!=0){
            gr_spaceSlices_folded->SetPoint(count_spaceSlices_folded,mT,mX);
            gr_spaceSlices_folded->SetPointError(count_spaceSlices_folded,mTerr,mXerr);
            count_spaceSlices_folded++;
            if (mX<graph_sepX){
                gr_combined_folded->SetPoint(count_combined_folded,mT,mX);
                gr_combined_folded->SetPointError(count_combined_folded,mTerr,mXerr);
                count_combined_folded++;
            }
        }
    }
    // go to the output file and save the graphs
    mOutFile->cd();
    gr_timeSlices_unfolded_left->Write();
    gr_timeSlices_unfolded_right->Write();
    gr_timeSlices_folded->Write();
    gr_spaceSlices_unfolded_left->Write();
    gr_spaceSlices_unfolded_right->Write();
    gr_spaceSlices_folded->Write();
    gr_combined_unfolded_left->Write();
    gr_combined_unfolded_right->Write();
    gr_combined_folded->Write();

    // draw the unfolded histogram first
    TCanvas * canv= new TCanvas("cgraph","cgraph",1024,768);canv->SetGridx(1);canv->SetGridy(1);
    h2_xt->Draw("COLZ");
    gr_timeSlices_unfolded_left->Draw("PLSAME");
    gr_timeSlices_unfolded_right->Draw("PLSAME");
    gr_spaceSlices_unfolded_left->Draw("PLSAME");
    gr_spaceSlices_unfolded_right->Draw("PLSAME");
    canv->SaveAs(Form("result/sampleUnfolded_%s.layer%d.png",mRunName.Data(),mTestLayerID));
    canv->SaveAs(Form("result/sampleUnfolded_%s.layer%d.pdf",mRunName.Data(),mTestLayerID));
    // draw the folded histogram first
    h2_xtn->Draw("COLZ");
    gr_timeSlices_folded->Draw("PLSAME");
    gr_spaceSlices_folded->Draw("PLSAME");
    canv->SaveAs(Form("result/sampleFolded_%s.layer%d.png",mRunName.Data(),mTestLayerID));
    canv->SaveAs(Form("result/sampleFolded_%s.layer%d.pdf",mRunName.Data(),mTestLayerID));

    // show the difference among left right and folded: left-folded and right-folded for time/space samplings
    TGraphErrors * gr_spaceSlices_leftMinusFolded = minusGraph(gr_spaceSlices_unfolded_left,gr_spaceSlices_folded,Form("gr_%d_x_lmf",mTestLayerID),true);
    gr_spaceSlices_leftMinusFolded->Print();
    gr_spaceSlices_leftMinusFolded->SetMarkerStyle(20);gr_spaceSlices_leftMinusFolded->SetMarkerSize(markerSize);
    gr_spaceSlices_leftMinusFolded->SetMarkerColor(kBlue);gr_spaceSlices_leftMinusFolded->SetLineColor(kBlue);
    TGraphErrors * gr_spaceSlices_rightMinusFolded = minusGraph(gr_spaceSlices_unfolded_right,gr_spaceSlices_folded,Form("gr_%d_x_rmf",mTestLayerID));
    gr_spaceSlices_rightMinusFolded->SetMarkerStyle(20);gr_spaceSlices_rightMinusFolded->SetMarkerSize(markerSize);
    gr_spaceSlices_rightMinusFolded->SetMarkerColor(kRed);gr_spaceSlices_rightMinusFolded->SetLineColor(kRed);
    TGraphErrors * gr_timeSlices_leftMinusFolded = minusGraph(gr_timeSlices_unfolded_left,gr_timeSlices_folded,Form("gr_%d_t_lmf",mTestLayerID),true);
    gr_timeSlices_leftMinusFolded->SetMarkerStyle(24);gr_timeSlices_leftMinusFolded->SetMarkerSize(markerSize);
    gr_timeSlices_leftMinusFolded->SetMarkerColor(kBlue);gr_timeSlices_leftMinusFolded->SetLineColor(kBlue);
    TGraphErrors * gr_timeSlices_rightMinusFolded = minusGraph(gr_timeSlices_unfolded_right,gr_timeSlices_folded,Form("gr_%d_t_rmf",mTestLayerID));
    gr_timeSlices_rightMinusFolded->SetMarkerStyle(24);gr_timeSlices_rightMinusFolded->SetMarkerSize(markerSize);
    gr_timeSlices_rightMinusFolded->SetMarkerColor(kRed);gr_timeSlices_rightMinusFolded->SetLineColor(kRed);
    // show left right difference
    canv->Clear();
    canv->Divide(2,2);
    TH2D * h2_bkg = new TH2D("h2_bkg","All samples",1024,-50,800,512,-10,10);
    h2_bkg->GetXaxis()->SetTitle("T [ns]");
    h2_bkg->GetYaxis()->SetTitle("X [mm]");
    canv->cd(1);gPad->SetGridx(1);gPad->SetGridy(1);
    h2_bkg->Draw();
    gr_spaceSlices_unfolded_left->Draw("PLSAME");
    gr_spaceSlices_unfolded_right->Draw("PLSAME");
    gr_spaceSlices_folded->Draw("PLSAME");
    gr_timeSlices_unfolded_left->Draw("PLSAME");
    gr_timeSlices_unfolded_right->Draw("PLSAME");
    gr_timeSlices_folded->Draw("PLSAME");
    TH2D * h2_bkg_diff = new TH2D("h2_bkg_diff","Difference between left/right side samplings and folded sampling",1024,-50,800,512,-0.5,0.5);
    h2_bkg_diff->GetXaxis()->SetTitle("T [ns]");
    h2_bkg_diff->GetYaxis()->SetTitle("#Delta_{X} [mm]");
    canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
    h2_bkg_diff->Draw();
    gr_spaceSlices_leftMinusFolded->Draw("PLSAME");
    gr_spaceSlices_rightMinusFolded->Draw("PLSAME");
    canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
    h2_bkg_diff->Draw();
    gr_timeSlices_leftMinusFolded->Draw("PLSAME");
    gr_timeSlices_rightMinusFolded->Draw("PLSAME");
    // show the difference between time and space samples: left right and folded 
    TGraphErrors * gr_timeMinusSpace_left = minusGraph(gr_timeSlices_unfolded_left,gr_spaceSlices_unfolded_left,Form("gr_%d_tmx_l",mTestLayerID));
    gr_timeMinusSpace_left->SetMarkerStyle(20);gr_timeMinusSpace_left->SetMarkerSize(markerSize);
    gr_timeMinusSpace_left->SetMarkerColor(kBlue);gr_timeMinusSpace_left->SetLineColor(kBlue);
    TGraphErrors * gr_timeMinusSpace_right = minusGraph(gr_timeSlices_unfolded_right,gr_spaceSlices_unfolded_right,Form("gr_%d_tmx_l",mTestLayerID));
    gr_timeMinusSpace_right->SetMarkerStyle(20);gr_timeMinusSpace_right->SetMarkerSize(markerSize);
    gr_timeMinusSpace_right->SetMarkerColor(kRed);gr_timeMinusSpace_right->SetLineColor(kRed);
    TGraphErrors * gr_timeMinusSpace_folded = minusGraph(gr_timeSlices_folded,gr_spaceSlices_folded,Form("gr_%d_tmx_l",mTestLayerID));
    gr_timeMinusSpace_folded->SetMarkerStyle(20);gr_timeMinusSpace_folded->SetMarkerSize(markerSize);
    gr_timeMinusSpace_folded->SetMarkerColor(kBlack);gr_timeMinusSpace_folded->SetLineColor(kBlack);
    canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1);
    h2_bkg_diff->Draw();
    gr_timeMinusSpace_left->Draw("PLSAME");
    gr_timeMinusSpace_right->Draw("PLSAME");
    gr_timeMinusSpace_folded->Draw("PLSAME");
    canv->SaveAs(Form("result/graph_%s.layer%d.png",mRunName.Data(),mTestLayerID));
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
    int minEntries = ParameterManager::Get().XTAnalyzerParameters.bin_n_min; // minimum number of entries in one slice to apply fitting function; Otherwise use mean value & RMS instead.
    if (!minEntries||entries>minEntries){
        TF1 * f = fitSliceGaus(h,mean,error,sigma,chi2,left,right);// Try with Gaussian function first
        h->GetYaxis()->SetRangeUser(0,h->GetBinContent(h->GetMaximumBin())*1.1);h->GetXaxis()->SetRangeUser(mean-RMS*5,mean+RMS*5); // for drawing
        if (mDrawDetails) drawFitting(h,f,canv,Form("%s, fit in (%.2e,%.2e), mean=%.2e+-%.2e, #sigma=%.2e, #chi^{2}=%.2e",drawTitle.Data(),left,right,mean,error,sigma,chi2),Form("%s_Gauss_%s.layer%d.png",h->GetName(),mRunName.Data(),mTestLayerID),left,mean,right);
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
            if (mDrawDetails) drawFitting(h,f,canv,Form("%s, fit in (%.2e,%.2e), mean=%.2e+-%.2e, #sigma=%.2e, #chi^{2}=%.2e",drawTitle.Data(),tempLeft,tempRight,tempMean,tempError,tempSig,tempChi2),Form("%s_LandauFlipped_%s.layer%d.png",h->GetName(),mRunName.Data(),mTestLayerID),tempLeft,tempMean,tempRight);
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
