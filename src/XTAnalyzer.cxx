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
#include "TLatex.h"
#include "TLegend.h"
#include "TGaxis.h"

#include <math.h>

#include "CommonTools.hxx"
#include "HistogramAnalyzer.hxx"
#include "ParameterManager.hxx"
#include "GeometryManager.hxx"
#include "Log.hxx"

XTAnalyzer::XTAnalyzer(TString runname, TFile * outfile, bool drawDetails)
        :mRunName(runname), mDrawDetails(drawDetails),
        mOutFile(outfile), mOutTree(NULL),
        mX(0), mXerr(0), mT(0), mTerr(0), mSig(0), mChi2(0), mProb(0), mEntries(0), mFunction(0), 
        h2_xt(NULL), h2_xt_combined(NULL), gr_left(NULL), gr_right(NULL), gr_rightMinusLeft(NULL), gr_combined(NULL),
        f_left(NULL), f_right(NULL), f_combinedLeft(NULL), f_combinedRight(NULL)
{
    // prepare XT basic functions
    for (int iRange = 0; iRange<NRANGE; iRange++){
        for (int iPol = 0; iPol<NPOL; iPol++){
            f_basicXT[iRange][iPol] = CommonTools::TF1New(Form("f_basicXT_%d_%d",iRange,iPol),Form("pol%d",iPol),-1000,1000);
            f_basicXT[iRange][iPol]->SetLineColor(kGray);
        }
    }
}

XTAnalyzer::~XTAnalyzer(void){
    for (int iLR = 0; iLR<2; iLR++){
        for (int iRange = 0; iRange<NRANGE; iRange++){
            for (int iPol = 0; iPol<NPOL; iPol++){
                if(f_basicXT[iRange][iPol]) delete f_basicXT[iRange][iPol];
            }
        }
    }
}

void XTAnalyzer::Initialize(void){
    if(mOutTree) {delete mOutTree; mOutTree = NULL;}
    if(h2_xt) {delete h2_xt; h2_xt = NULL;}
    if(h2_xt_combined) {delete h2_xt_combined; h2_xt_combined = NULL;}
    if(gr_left) {delete gr_left; gr_left = NULL;}
    if(gr_right) {delete gr_right ; gr_right = NULL;}
    if(gr_rightMinusLeft) {delete gr_rightMinusLeft; gr_rightMinusLeft = NULL;}
    if(gr_combined) {delete gr_combined; gr_combined = NULL;}
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
    mDrawTmin = ParameterManager::Get().XTAnalyzerParameters.draw_tmin;
    mDrawTmax = ParameterManager::Get().XTAnalyzerParameters.draw_tmax;
    mDrawXmin = ParameterManager::Get().XTAnalyzerParameters.draw_xmin;
    mDrawXmax = ParameterManager::Get().XTAnalyzerParameters.draw_xmax;
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
        gr_rightMinusLeft->SetMarkerColor(kMagenta);gr_rightMinusLeft->SetLineColor(kMagenta);
    }
    return 0;
}

int XTAnalyzer::PrepareXTFunctions(){
    int nRanges = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nRanges;
    TString formula("");
    int counter = 0; int previous_range_counter = 0;
    for (int iRange = 0; iRange<nRanges; iRange++){
        int nPol = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nPol[iRange];
        if (iRange!=0) formula+="+";
        formula+=Form("pol%d(%d)",nPol,counter);
        counter+=nPol+1;
        if (iRange==nRanges-1){
            formula+=Form("*(x>[%d])",previous_range_counter);
        }
        else{
            if (iRange==0){
                formula+=Form("*(x<=[%d])",counter);
            }
            else{
                formula+=Form("*(x>[%d]&&x<=[%d])",previous_range_counter,counter);
            }
            previous_range_counter = counter;
            counter+=1;
        }
    }
    formula+=Form("+[%d]",counter); // the last parameter is serving as an offset to cancel wire position offset if needed
    if (f_left) delete f_left;
    if (f_right) delete f_right;
    f_left = CommonTools::TF1New("fl"+m_suffix,formula,-25,800);
    f_right = CommonTools::TF1New("fr"+m_suffix,formula,-25,800);
    if (ParameterManager::Get().XTAnalyzerParameters.CombineLeftAndRight){
        f_combinedLeft = CommonTools::TF1New("flc"+m_suffix,formula,-25,800);
        f_combinedRight = CommonTools::TF1New("frc"+m_suffix,formula,-25,800);
    }
    return 0;
}

void XTAnalyzer::Fill(double t, double x){
    if (!h2_xt){
        MyError("h2_xt is not prepared yet!");
        return;
    }
    h2_xt->Fill(t,x);
}

void XTAnalyzer::Write(){
    h2_xt->Write();
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
                HistogramAnalyzer::Get().SetFittingParameters(current_iRange);
                double x1,xerr1,sig1,x2,xerr2,sig2; int result;
                result = HistogramAnalyzer::Get().FitSlice(hist,mChi2,mProb,fitBoth,iLR==0);
                mFunction = HistogramAnalyzer::Get().get_functionType();
                if (fitBoth){
                    x1 = HistogramAnalyzer::Get().get_xL();
                    xerr1 = HistogramAnalyzer::Get().get_xerrL();
                    sig1 = HistogramAnalyzer::Get().get_sigL();
                    x2 = HistogramAnalyzer::Get().get_xR();
                    xerr2 = HistogramAnalyzer::Get().get_xerrR();
                    sig2 = HistogramAnalyzer::Get().get_sigR();
                }
                else{
                    x1 = HistogramAnalyzer::Get().get_x();
                    xerr1 = HistogramAnalyzer::Get().get_xerr();
                    sig1 = HistogramAnalyzer::Get().get_sig();
                }
                if (mDrawDetails){
                    if (fitBoth){
                        hist->GetXaxis()->SetRangeUser((x1>0?0:x1)-0.6,(x2<0?0:x2)+0.6);
                        TString x1string = x1<-1?Form("%.2f mm",x1):Form("%.0f #mum",x1*1000);
                        TString x2string = x2>1?Form("%.2f mm",x2):Form("%.0f #mum",x2*1000);
                        hist->SetTitle(Form("t = %.1f (%.1f~%.1f) ns, x_{1}=%s (%.0f #mum), #sigma_{1}=%.0f #mum, x_{2}=%s (%.0f #mum), #sigma_{2}=%.0f #mum, #chi^{2}=%.0f, prob=%.2f, stat %d",divmiddle,divleft,divright,x1string.Data(),xerr1*1000,sig1*1000,x2string.Data(),xerr2*1000,sig2*1000,mChi2,mProb,result));
                    }
                    else{
                        hist->GetXaxis()->SetRangeUser(x1-0.6,x1+0.6);
                        TString xstring = x1<-1?Form("%.2f mm",x1):Form("%.0f #mum",x1*1000);
                        hist->SetTitle(Form("t = %.1f (%.1f~%.1f) ns, x=%s (%.0f #mum), #sigma=%.0f #mum, #chi^{2}=%.0f, prob=%.2f, stat %d",t1,divleft,divright,xstring.Data(),xerr1*1000,sig1*1000,mChi2,mProb,result));
                    }
                    canv->cd();
                    HistogramAnalyzer::Get().DrawFitting(hist);
                    canv->SaveAs(Form("%s.%s%s.png",hist->GetName(),mRunName.Data(),m_suffix.Data()));
                }
                // update the canvas range for drawing samples
                mT = t1; mTerr = terr1; mEntries = n1; mX = x1; mXerr = xerr1; mSig = sig1;
                mOutTree->Fill();
                if (fitBoth){
                    mT = t2; mTerr = terr2; mEntries = n2; mX = x2; mXerr = xerr2; mSig = sig2;
                    mOutTree->Fill();
                }
            }
            else{
                mFunction = HistogramAnalyzer::kNone;
                mT = t1; mTerr = terr1; mEntries = n1; mX = 0; mXerr = 0; mSig = 0;
                mOutTree->Fill();
                if (fitBoth){
                    mT = t2; mTerr = terr2; mEntries = n2; mX = 0; mXerr = 0; mSig = 0;
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

    //==========================Draw the samples==============================
    formXTGraphs();
    gStyle->SetOptStat(0);
    drawSampleAtt();
    drawSample2D(false);
    // set the style back
    gStyle->SetOptStat(oldStyle);

    //==========================Draw the samples==============================
    mOutFile->cd();
    mOutTree->Write();
}

void XTAnalyzer::FitXT(){
    int nRanges = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nRanges;
    if (nRanges>NRANGE){
        MyError("Cannot handle "<<nRanges<<" ranges. Maximum: "<<NRANGE);
        return;
    }

    doFitXT(gr_left,f_left,"L");
    doFitXT(gr_right,f_right,"R");

    drawSample2D(true);

    // move the XT to center: we don't want to leave the wire position offset information in this XT relation
    double tmin = ParameterManager::Get().XTAnalyzerParameters.gold_t_min;
    double tmax = ParameterManager::Get().XTAnalyzerParameters.gold_t_max;
    double offset = (f_left->Integral(tmin,tmax)+f_right->Integral(tmin,tmax))/2./(tmax-tmin);
    int nPars = f_left->GetNpar();
    f_left->SetParameter(nPars-1,-offset);
    nPars = f_right->GetNpar();
    f_right->SetParameter(nPars-1,-offset);

    // save the newly created objects
    mOutFile->cd();
    f_left->Write();
    f_right->Write();

    // Combine the two side together if needed
    if (ParameterManager::Get().XTAnalyzerParameters.CombineLeftAndRight){
        combineLeftAndRight(offset); // the offset should be considered while making gr_combined
        doFitXT(gr_combined,f_combinedRight,"C");
        int counter = 0;
        for (int iRange = 0; iRange<nRanges; iRange++){
            if (iRange>0){
                double par = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tLowEdge[iRange];
                f_combinedLeft->SetParameter(counter++,par);
            }
            int nPol = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nPol[iRange];
            for (int iPar = 0; iPar<=nPol; iPar++){
                double par = -f_basicXT[iRange][nPol]->GetParameter(iPar); // flip it
                f_combinedLeft->SetParameter(counter++,par);
            }
        }
        f_combinedLeft->SetParameter(counter++,0); // the offset is 0

        drawSampleCombined2D();

        h2_xt_combined->Write();
        gr_combined->Write();
        f_combinedRight->Write();
    }
}

void XTAnalyzer::doFitXT(TGraph * gr, TF1 * f, const char * suf){
    int nRanges = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nRanges;
    // to draw the basic function fittings
    TCanvas * canv_basics[NRANGE];
    TPad* pads[NRANGE];
    TLatex * tex[NRANGE];
    for (int iRange = 0; iRange<nRanges; iRange++){
        canv_basics[iRange] = new TCanvas(Form("canv_xtfunc_basic_%d",iRange),"",1024,1024);
        pads[iRange] = new TPad(Form("pad%d",iRange),"pad",0,0,1,0.95);
        pads[iRange]->Draw();
        pads[iRange]->Divide(3,3);
        tex[iRange]  = new TLatex(0.1,0.96,"");
        tex[iRange]->SetTextSize(0.025);
    }

    // loop in ranges to perform basic fittings
    for (int iRange = 0; iRange<nRanges; iRange++){
        // get the T region for this range
        double tLeft = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tLeft[iRange];
        double tRight = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tRight[iRange];
        MyNamedInfo("XTAnalyzer","==> Fitting in T range "<<tLeft<<" ~ "<<tRight<<" ns");
        int N=0; // number of points within this T region
        int iStart=-1; // the first point index in this region
        double minX = 1e9; double maxX=-1e9; // to decide the range of the histogram to be used as background to draw fitting result in this range
        double minT = 1e9; double maxT=-1e0;
        for (size_t iPoint = 0; iPoint<gr->GetN(); iPoint++){
            double x,t;
            gr->GetPoint(iPoint,t,x);
            if (t>=tLeft&&t<=tRight){
                if (iStart==-1) iStart=iPoint;
                N++;
                if (minX>x) minX = x;
                if (maxX<x) maxX = x;
                if (minT>t) minT = t;
                if (maxT<t) maxT = t;
            }
        }
        // update the title and draw to the basic fitting canvas
        tex[iRange]->SetText(0.1,0.96,Form("%s %s, Fit in T %.0f ~ %.0f ns",mRunName.Data(),gr->GetName(),tLeft,tRight));
        canv_basics[iRange]->cd(); tex[iRange]->Draw();
        // loop in basic polynomial functions with different orders
        for (int iPol = 1; iPol<=9; iPol++){
            // do the fitting
            int iTry = 0;
            int fitResult = 0;
            for (; iTry<10; iTry++){
                if (iPol>1){
                    for (int iPar = 0; iPar<iPol; iPar++){
                        double par = f_basicXT[iRange][iPol-1]->GetParameter(iPar);
                        if (fabs(par)>1e6) par=0;
                        //if (iTry==0) par*=(1+gRandom->Uniform(-0.01,0.01));
                        //if (iTry>0) par*=(1+gRandom->Uniform(-0.1,0.1));
                        par*=(1+gRandom->Uniform(-0.1,0.1));
                        f_basicXT[iRange][iPol]->SetParameter(iPar,par);
                    }
                    f_basicXT[iRange][iPol]->SetParameter(iPol,0);
                }
                else{
                    f_basicXT[iRange][iPol]->SetParameters(0,0);
                }
                if (iRange<=1)
                    fitResult = gr->Fit(f_basicXT[iRange][iPol],"qN0","",tLeft,tRight);
                else
                    fitResult = gr->Fit(f_basicXT[iRange][iPol],"qN0","",tLeft-10,tRight+10);
                double maxPar = 0;
                for (int iPar = 0; iPar<=iPol; iPar++){
                    double par = f_basicXT[iRange][iPol]->GetParameter(iPar);
                    if (fabs(par)>fabs(maxPar)) maxPar = par;
                }
                if (!fitResult&&fabs(maxPar)<1e6&&f_basicXT[iRange][iPol]->GetChisquare()<1e4){
                    break;
                }
            }
            MyNamedInfo("XTAnalyzer","After "<<iTry<<" fitting, Pol"<<iPol<<":"<<" chi2 = "<<f_basicXT[iRange][iPol]->GetChisquare()<<" result "<<fitResult);
            if (iTry==10) {
                MyNamedWarn("XTAnalyzer","Failed to fit Pol"<<iPol);
            }
            // now make a new graph with error to record the difference between the fitted function and the original graph points.
            TGraphErrors * gr_diff = new TGraphErrors(); gr_diff->SetName(Form("gr_diff%s.%s_%d_%.0f_%.0f",m_suffix.Data(),suf,iPol,tLeft,tRight)); gr_diff->SetTitle(Form("Fit in T in %.0f~%.0f ns",tLeft,tRight));
            gr_diff->SetLineColor(kMagenta); gr_diff->SetMarkerColor(kMagenta); gr_diff->SetMarkerStyle(24);
            gr_diff->Set(N);
            double delaXRange = 0.2;
            double scale = (maxX-minX+0.5)/delaXRange; double offset = (minX+maxX)/2;
            double deltaXMax = 0;
            for (int iPoint = 0; iPoint<N; iPoint++){
                double x,t;
                gr->GetPoint(iPoint+iStart,t,x);
                double terr=gr->GetErrorX(iPoint+iStart);
                double xerr=gr->GetErrorY(iPoint+iStart);
                double deltaX = f_basicXT[iRange][iPol]->Eval(t)-x;
                if (fabs(deltaXMax)<fabs(deltaX)) deltaXMax = deltaX;
                gr_diff->SetPoint(iPoint,t,deltaX*scale+offset);
                gr_diff->SetPointError(iPoint,terr,xerr*scale);
            }
            // draw the fitting result with basic functions
            pads[iRange]->cd(iPol);gPad->SetGridx(1);gPad->SetGridy(1);
            TH2D * hbkg = new TH2D(Form("hbkg_xtbasic%s.%s_%d_%d",m_suffix.Data(),suf,iRange,iPol),Form("Pol%d, #Delta_{X} RMS %.0f (%.0f) #mum",iPol,gr_diff->GetRMS(2)/scale*1000,deltaXMax*1000),1024,minT-5,maxT+5,512,minX-0.25,maxX+0.25);
            hbkg->GetXaxis()->SetTitle("T [ns]");
            hbkg->GetYaxis()->SetTitle("X [mm]");
            hbkg->GetYaxis()->SetTitleOffset(1.1);
            hbkg->Draw();
            TGaxis * axis = new TGaxis(maxT+5,minX-0.25,maxT+5,maxX+0.25,-delaXRange/2*1000,delaXRange/2*1000,510,"+L");
            axis->SetTitle("#Delta_{X} #mum");
            axis->SetTitleOffset(1.1);
            axis->Draw();
            axis->SetTitleFont(42); axis->SetLabelFont(42);
            gr->Draw("PLSAME");
            f_basicXT[iRange][iPol]->Draw("SAME");
            gr_diff->Draw("PLSAME");
        }
        canv_basics[iRange]->SaveAs(Form("result/fitXTBasics_%s%s_%s_%d.png",mRunName.Data(),m_suffix.Data(),suf,iRange));
    } // end of ranges loop

    // at last, let's pick up one basic polynomial function in each range according to the specific given orders and combine them into one wholestic XT
    int counter = 0;
    for (int iRange = 0; iRange<nRanges; iRange++){
        if (iRange>0){
            double par = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tLowEdge[iRange];
            f->SetParameter(counter++,par);
        }
        int nPol = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nPol[iRange];
        for (int iPar = 0; iPar<=nPol; iPar++){
            double par = f_basicXT[iRange][nPol]->GetParameter(iPar);
            f->SetParameter(counter++,par);
        }
    }
    f->SetParameter(counter++,0); // the offset is by default 0; will be updated after two functions are formed.
    // decide the range of the function
    double tmin = 1e14; double tmax =-1e14;
    for (size_t iPoint = 0; iPoint<gr->GetN(); iPoint++){
        double t,x;
        gr->GetPoint(iPoint,t,x);
        if (tmax<t) tmax = t;
        if (tmin>t) tmin = t;
    }
    double range_tl = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tLowEdge[0];
    double range_th = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tHighEdge;
    if (tmax<range_th) range_th = tmax;
    if (tmin>range_tl) range_tl = tmin;
    f->SetRange(range_tl,range_th);

    // make a new canvas to show the conjunction part
    TCanvas * canv2 = new TCanvas("cgraph2","cgraph",1024,768);
    canv2->Divide(nRanges-1);
    for (int iRange = 1; iRange<nRanges; iRange++){
        canv2->cd(iRange);gPad->SetGridx(1);gPad->SetGridy(1);
        double tmin = ParameterManager::Get().XTAnalyzerParameters.xtfunc_tLowEdge[iRange]-20;
        double tmax = tmin+40;
        double xmin,xmax;
        xmin = f->Eval(tmax)-0.5;
        xmax = f->Eval(tmin)+0.5;
        if (xmin>xmax){double temp = xmin; xmin = xmax; xmax = temp;}
        TH2D * hbkg = new TH2D(Form("hbkg_conjonction%s.%s.%d",m_suffix.Data(),suf,iRange),Form("%s: Range %d to Range %d",suf,iRange-1,iRange),512,tmin,tmax,512,xmin,xmax);
        hbkg->GetXaxis()->SetTitle("T [ns]");
        hbkg->GetYaxis()->SetTitle("X [mm]");
        hbkg->Draw();
        gr->Draw("PLSAME");
        int nPol = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nPol[iRange];
        f_basicXT[iRange][nPol]->Draw("SAME");
        int nPolPre = ParameterManager::Get().XTAnalyzerParameters.xtfunc_nPol[iRange-1];
        f_basicXT[iRange-1][nPolPre]->Draw("SAME");
        f->Draw("SAME");
    }
    canv2->SaveAs(Form("result/sampleConj_%s%s_%s.png",mRunName.Data(),m_suffix.Data(),suf));

    // remove the canvases
    for (int iRange = 0; iRange<nRanges; iRange++){
        if (canv_basics[iRange]) delete canv_basics[iRange];
    }
}

void XTAnalyzer::combineLeftAndRight(double offset){
    double combineAtDOCA = ParameterManager::Get().XTAnalyzerParameters.CombineAtDOCA;
    // get the 2D histogram first
    // about binning
    h2_xt_combined = new TH2D(*h2_xt);
    h2_xt_combined->SetName("h2_xtc"+m_suffix);
    for (int iBinY = 1; iBinY <= h2_xt->GetYaxis()->GetNbins(); iBinY++){
        double Y = h2_xt->GetYaxis()->GetBinCenter(iBinY);
        if ((Y<combineAtDOCA&&Y>=0)||(Y>=-combineAtDOCA&&Y<0)){ // TODO: to consider more combination patterns in the future (quite rarely used feature)
            int iBinY2Copy = h2_xt->GetYaxis()->FindBin(-Y);
            for (int iBinX = 1; iBinX <= h2_xt->GetXaxis()->GetNbins(); iBinX++){
                double c = h2_xt->GetBinContent(iBinX,iBinY);
                h2_xt_combined->SetBinContent(iBinX,iBinY,c);
                h2_xt_combined->SetBinContent(iBinX,iBinY2Copy,c);
            }
        }
    }

    // then get the graph
    gr_combined = new TGraphErrors(); gr_combined->SetName(Form("gr%s_c",m_suffix.Data()));
    gr_combined->SetMarkerStyle(24);gr_combined->SetMarkerSize(0.5);
    gr_combined->SetMarkerColor(kBlack);gr_combined->SetLineColor(kBlack);
    gr_combined->Set(gr_left->GetN()+gr_right->GetN());
    int count = 0;
    for (int iPoint = 0; iPoint<gr_left->GetN(); iPoint++){
        double T,X; gr_left->GetPoint(iPoint,T,X);
        double Terr = gr_left->GetErrorX(iPoint);
        double Xerr = gr_left->GetErrorY(iPoint);
        if (X>=-combineAtDOCA){
            gr_combined->SetPoint(count,T,-X);
            gr_combined->SetPointError(count,Terr,Xerr);
            count++;
        }
    }
    for (int iPoint = 0; iPoint<gr_right->GetN(); iPoint++){
        double T,X; gr_right->GetPoint(iPoint,T,X);
        double Terr = gr_right->GetErrorX(iPoint);
        double Xerr = gr_right->GetErrorY(iPoint);
        if (X>combineAtDOCA){
            gr_combined->SetPoint(count,T,X);
            gr_combined->SetPointError(count,Terr,Xerr);
            count++;
        }
    }
    gr_combined->Set(count);
}

void XTAnalyzer::drawSampleCombined2D(){
    double combineAtDOCA = ParameterManager::Get().XTAnalyzerParameters.CombineAtDOCA;
    // draw the unfolded histogram first
    TCanvas * canv = new TCanvas("cgraph","cgraph",1024,1024);
    canv->Divide(2,2);

    canv->cd(1);gPad->SetGridx(1);gPad->SetGridy(1);
    h2_xt_combined->GetXaxis()->SetRangeUser(mDrawTmin,mDrawTmax);
    h2_xt_combined->GetYaxis()->SetRangeUser(0,mDrawXmax);
    h2_xt_combined->Draw("COLZ");
    gr_combined->SetMarkerSize(0.2);
    gr_combined->Draw("PLSAME");

    canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_sigT = new TH2D("h2_bkg_sigT2","#sigma of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,500);
    h2_bkg_sigT->GetXaxis()->SetTitle("Drift Time [ns]");
    h2_bkg_sigT->GetYaxis()->SetTitle("#sigma [#mum]");
    h2_bkg_sigT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("sig*1000:t",Form("x>=%.7e&&n>0&&func>=0",combineAtDOCA),"PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("sig*1000:t",Form("x<0&&x>=%.7e&&n>0&&func>=0",-combineAtDOCA),"PSAME");

    canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * hbkg = new TH2D(Form("hbkg_diff%s",m_suffix.Data()),"#Delta_{X}: Function - Sample",1024,mDrawTmin,mDrawTmax,512,-500,500);
    hbkg->GetXaxis()->SetTitle("T [ns]");
    hbkg->GetYaxis()->SetTitle("#Delta_{X} [#mum]");
    hbkg->Draw();
    // now make a new graph with error to record the difference between the fitted function and the original graph points.
    TGraphErrors * gr_diff = new TGraphErrors(); gr_diff->SetName(Form("gr_diff%s.%s",m_suffix.Data(),"C"));
    gr_diff->SetLineColor(gr_combined->GetLineColor()); gr_diff->SetMarkerColor(gr_combined->GetMarkerColor()); gr_diff->SetMarkerStyle(24);
    int N = gr_combined->GetN();
    gr_diff->Set(N);
    double deltaXMax = 0;
    for (int iPoint = 0; iPoint<N; iPoint++){
        double x,t;
        gr_combined->GetPoint(iPoint,t,x);
        double terr=gr_combined->GetErrorX(iPoint);
        double xerr=gr_combined->GetErrorY(iPoint);
        double deltaX = f_combinedRight->Eval(t)-x;
        if (fabs(deltaXMax)<fabs(deltaX)) deltaXMax = deltaX;
        gr_diff->SetPoint(iPoint,t,deltaX*1000);
        gr_diff->SetPointError(iPoint,terr,xerr*1000);
    }
    gr_diff->Draw("PLSAME");

    canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1);
    int mEntriesMax = 0;
    for (Long64_t iEntry = 0; iEntry<mOutTree->GetEntries(); iEntry++){
        mOutTree->GetEntry(iEntry);
        if (mEntries>mEntriesMax) mEntriesMax = mEntries;
    }
    TH2D * h2_bkg_entriesT = new TH2D("h2_bkg_entriesT","Number of entries in each T slice",512,mDrawTmin,mDrawTmax,1024,0,mEntriesMax*1.1);
    h2_bkg_entriesT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("n:t",Form("x>=%.7e&&n>0&&func>=0",combineAtDOCA),"PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("n:t",Form("x<0&&x>=%.7e&&n>0&&func>=0",-combineAtDOCA),"PSAME");
    mOutTree->SetMarkerColor(kGray); mOutTree->Draw("n:t",Form("((x<0&&x>=%.7e)||x>=%.7e)&&func<0",-combineAtDOCA,combineAtDOCA),"PSAME");

    canv->SaveAs(Form("result/sample2DCombined_%s%s.png",mRunName.Data(),m_suffix.Data()));
    canv->SaveAs(Form("result/sample2DCombined_%s%s.pdf",mRunName.Data(),m_suffix.Data()));
    h2_xt_combined->GetXaxis()->UnZoom();
    h2_xt_combined->GetYaxis()->UnZoom();
}

void XTAnalyzer::drawSample2D(bool withFunction){
    // draw the unfolded histogram first
    TCanvas * canv = new TCanvas("cgraph","cgraph",1024,1024);
    canv->Divide(2,2);

    canv->cd(1);gPad->SetGridx(1);gPad->SetGridy(1);
    h2_xt->GetXaxis()->SetRangeUser(mDrawTmin,mDrawTmax);
    h2_xt->GetYaxis()->SetRangeUser(mDrawXmin,mDrawXmax);
    h2_xt->Draw("COLZ");
    gr_left->SetMarkerSize(0.2);
    gr_right->SetMarkerSize(0.2);
    gr_left->Draw("PLSAME");
    gr_right->Draw("PLSAME");

    canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_sigT = new TH2D("h2_bkg_sigT2","#sigma of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,500);
    h2_bkg_sigT->GetXaxis()->SetTitle("Drift Time [ns]");
    h2_bkg_sigT->GetYaxis()->SetTitle("#sigma [#mum]");
    h2_bkg_sigT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("sig*1000:t","x>=0&&n>0&&func>=0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("sig*1000:t","x<0&&n>0&&func>=0","PSAME");

    canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_rightMinusLeft = new TH2D("h2_bkg_rightMinusLeft","Mean of the two side of XT relation",512,mDrawTmin,mDrawTmax,1024,-1,1);
    h2_bkg_rightMinusLeft->Draw();
    gr_rightMinusLeft->Draw("PLSAME");

    if (withFunction){
        double xtrange_tmin = f_right->GetXmin()<f_left->GetXmin()?f_left->GetXmin():f_right->GetXmin();
        double xtrange_tmax = f_right->GetXmax()>f_left->GetXmax()?f_left->GetXmax():f_right->GetXmax();
        TF1 * f_diff = CommonTools::TF1New("f_diff",Form("(fl%s+fr%s)/2",m_suffix.Data(),m_suffix.Data()),xtrange_tmin,xtrange_tmax);
        int nParLeft=f_left->GetNpar();
        for (int iPar = 0; iPar<nParLeft; iPar++){
            f_diff->SetParameter(iPar,f_left->GetParameter(iPar));
        }
        int nParRight=f_right->GetNpar();
        for (int iPar = 0; iPar<nParRight; iPar++){
            f_diff->SetParameter(iPar+nParLeft,f_right->GetParameter(iPar));
        }
        double maxdiff = f_diff->GetMinimum(xtrange_tmin,xtrange_tmax);
        if (fabs(maxdiff)<f_diff->GetMaximum(xtrange_tmin,xtrange_tmax)) maxdiff=f_diff->GetMaximum(xtrange_tmin,xtrange_tmax);
        h2_bkg_rightMinusLeft->SetTitle(Form("Right + Left: mean %.0f #mum, max %.0f #mum",f_diff->Integral(xtrange_tmin,xtrange_tmax)/(xtrange_tmax-xtrange_tmin)*1000,maxdiff*1000));
        f_diff->SetLineColor(kMagenta);
        canv->cd(4);
        f_diff->Draw("SAME");

        canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
        TH2D * hbkg = new TH2D(Form("hbkg_diff%s",m_suffix.Data()),"#Delta_{X}: Function - Sample",1024,mDrawTmin,mDrawTmax,512,-500,500);
        hbkg->GetXaxis()->SetTitle("T [ns]");
        hbkg->GetYaxis()->SetTitle("#Delta_{X} [#mum]");
        hbkg->Draw();
        TLegend *legend_diff = new TLegend(0.6,0.7,0.9,0.9);
        legend_diff->Draw("SAME");
        for (int iLR = 0; iLR<2; iLR++){
            TF1 * f = NULL; TGraphErrors * gr = NULL;
            if (iLR==0){
                f = f_left;
                gr = gr_left;
            }
            else{
                f = f_right;
                gr = gr_right;
            }
            canv->cd(1);
            f->Draw("SAME");

            // now make a new graph with error to record the difference between the fitted function and the original graph points.
            TGraphErrors * gr_diff = new TGraphErrors(); gr_diff->SetName(Form("gr_diff%s.%s",m_suffix.Data(),iLR==0?"L":"R"));
            gr_diff->SetLineColor(gr->GetLineColor()); gr_diff->SetMarkerColor(gr->GetMarkerColor()); gr_diff->SetMarkerStyle(24);
            int N = gr->GetN();
            gr_diff->Set(N);
            double deltaXMax = 0;
            for (int iPoint = 0; iPoint<N; iPoint++){
                double x,t;
                gr->GetPoint(iPoint,t,x);
                double terr=gr->GetErrorX(iPoint);
                double xerr=gr->GetErrorY(iPoint);
                double deltaX = f->Eval(t)-x;
                if (fabs(deltaXMax)<fabs(deltaX)) deltaXMax = deltaX;
                gr_diff->SetPoint(iPoint,t,deltaX*1000);
                gr_diff->SetPointError(iPoint,terr,xerr*1000);
            }
            canv->cd(3);
            gr_diff->Draw("PLSAME");
            legend_diff->AddEntry(gr_diff,Form("%s: RMS %.0f, Max %.0f #mum",iLR==0?"Left":"Right",gr_diff->GetRMS(2),deltaXMax*1000),"PL");
        }
    }
    canv->SaveAs(Form("result/sample2D_%s%s.png",mRunName.Data(),m_suffix.Data()));
    canv->SaveAs(Form("result/sample2D_%s%s.pdf",mRunName.Data(),m_suffix.Data()));
    h2_xt->GetXaxis()->UnZoom();
    h2_xt->GetYaxis()->UnZoom();
}

void XTAnalyzer::drawSampleAtt(){
    mOutTree->SetMarkerStyle(20);mOutTree->SetMarkerSize(0.5);
    TCanvas * canv= new TCanvas("csample","csample",1024,768);
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
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("n:t","x>=0&&n>0&&func>=0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("n:t","x<0&&n>0&&func>=0","PSAME");
    mOutTree->SetMarkerColor(kGray); mOutTree->Draw("n:t","func<0","PSAME");
    canv->cd(2);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_sigT = new TH2D("h2_bkg_sigT","#sigma of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,0.8);
    h2_bkg_sigT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("sig:t","x>=0&&n>0&&func>=0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("sig:t","x<0&&n>0&&func>=0","PSAME");
    canv->cd(3);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_chi2T = new TH2D("h2_bkg_chi2T","#chi^{2} of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,150);
    h2_bkg_chi2T->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("chi2:t","x>=0&&n>0&&func>=0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("chi2:t","x<0&&n>0&&func>=0","PSAME");
    canv->cd(4);gPad->SetGridx(1);gPad->SetGridy(1);
    TH2D * h2_bkg_probT = new TH2D("h2_bkg_probT","p-value of X fitting in each T slice",512,mDrawTmin,mDrawTmax,512,0,1);
    h2_bkg_probT->Draw();
    mOutTree->SetMarkerColor(kRed); mOutTree->Draw("prob:t","x>=0&&n>0&&func>=0","PSAME");
    mOutTree->SetMarkerColor(kBlue); mOutTree->Draw("prob:t","x<0&&n>0&&func>=0","PSAME");
    canv->SaveAs(Form("result/sampleAtt_%s%s.png",mRunName.Data(),m_suffix.Data()));
}

void XTAnalyzer::formXTGraphs(){
    // set parameters
    int    graph_n_min = ParameterManager::Get().XTAnalyzerParameters.graph_n_min;
    double graph_chi2_max = ParameterManager::Get().XTAnalyzerParameters.graph_chi2_max;
    double graph_prob_min = ParameterManager::Get().XTAnalyzerParameters.graph_prob_min;

    // make graphs from different samplings: left/right/folded TIMES time/space
    gr_left->Set(mOutTree->GetEntries(Form("func>=0&&n>=%d&&(%d||chi2<=%.7e)&&prob>=%.7e&&x<%.7e",graph_n_min,graph_chi2_max?0:1,graph_chi2_max,graph_prob_min,0.)));
    gr_right->Set(mOutTree->GetEntries(Form("func>=0&&n>=%d&&(%d||chi2<=%.7e)&&prob>=%.7e&&x>%.7e",graph_n_min,graph_chi2_max?0:1,graph_chi2_max,graph_prob_min,0.)));
    int count_left = 0;
    int count_right = 0;
    MyNamedVerbose("XTAnalyzer","Looping in the output tree again, "<<mOutTree->GetEntries()<<" entries");
    for (Long64_t iEntry = 0; iEntry<mOutTree->GetEntries(); iEntry++){
        mOutTree->GetEntry(iEntry);
        if (mFunction>=0&&mEntries>=graph_n_min&&(!graph_chi2_max||mChi2<=graph_chi2_max)&&mProb>=graph_prob_min){
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
    CommonTools::TGraphErrorsPlus(gr_rightMinusLeft,gr_right,gr_left,0.5,0.5); // mean of left and right
    mOutFile->cd();
    gr_left->Write();
    gr_right->Write();
    gr_rightMinusLeft->Write();
}
