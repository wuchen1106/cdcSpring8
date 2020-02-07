#include <cstddef>

#include "TF1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TStyle.h"

#include "HistogramAnalyzer.hxx"
#include "Log.hxx"
#include "ParameterManager.hxx"
#include "CommonTools.hxx"

HistogramAnalyzer* HistogramAnalyzer::fHistogramAnalyzer= NULL;

HistogramAnalyzer::HistogramAnalyzer():
    par_peak_height_middle(0),par_peak_height_left(0),par_peak_height_right(0),par_peak_sigma_middle(0),par_peak_sigma_left(0),par_peak_sigma_right(0),par_peak_mean_range(0),
    par_base_height_middle(0),par_base_height_left(0),par_base_height_right(0),par_base_sigma_middle(0),par_base_sigma_left(0),par_base_sigma_right(0),par_base_mean_range(0),
    par_functionType(0),
    f_gausL(NULL), f_gausR(NULL), f_gaus2L(NULL), f_gaus2R(NULL), f_gausBoth(NULL), f_landL(NULL), f_landR(NULL), f_land2L(NULL), f_land2R(NULL), f_landBoth(NULL),
    f_combDoubleGausL(NULL), f_combDoubleGausR(NULL), f_combDoubleGausBoth(NULL),
    f_combDoubleLandL(NULL), f_combDoubleLandR(NULL), f_combDoubleLandBoth(NULL),
    f_combGausLandL(NULL), f_combGausLandR(NULL), f_combGausLandBoth(NULL),
    f_combLandGausL(NULL), f_combLandGausR(NULL), f_combLandGausBoth(NULL),
    cur_function(NULL),cur_peak_function(NULL),cur_base_function(NULL),cur_functionL(NULL),cur_peak_functionL(NULL),cur_base_functionL(NULL),cur_functionR(NULL),cur_peak_functionR(NULL),cur_base_functionR(NULL),cur_functionType(0),
    fit_result(0),fit_chi2(0),fit_pValue(0),fit_x(0),fit_sig(0),fit_xerr(0),fit_xL(0),fit_sigL(0),fit_xerrL(0),fit_xR(0),fit_sigR(0),fit_xerrR(0)
{
    // single functions: [0]: height, [1]: mean, [2]: sigma
    f_gausL = CommonTools::TF1New("f_gausL","gaus",-1000,1000); f_gausL->SetLineStyle(2); f_gausL->SetLineColor(kRed);
    f_gausR = CommonTools::TF1New("f_gausR","gaus",-1000,1000); f_gausR->SetLineStyle(2); f_gausR->SetLineColor(kRed);
    f_gaus2L = CommonTools::TF1New("f_gaus2L","gaus",-1000,1000); f_gaus2L->SetLineStyle(2); f_gaus2L->SetLineColor(kBlue);
    f_gaus2R = CommonTools::TF1New("f_gaus2R","gaus",-1000,1000); f_gaus2R->SetLineStyle(2); f_gaus2R->SetLineColor(kBlue);
    f_landL = CommonTools::TF1New("f_landL","[0]/0.18065564*TMath::Landau(x,[1],[2],false)",-1000,1000); f_landL->SetLineStyle(2); f_landL->SetLineColor(kRed);
    f_landR = CommonTools::TF1New("f_landR","[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)",-1000,1000); f_landR->SetLineStyle(2); f_landR->SetLineColor(kRed);
    f_land2L = CommonTools::TF1New("f_land2L","[0]/0.18065564*TMath::Landau(x,[1],[2],false)",-1000,1000); f_land2L->SetLineStyle(2); f_land2L->SetLineColor(kBlue);
    f_land2R = CommonTools::TF1New("f_land2R","[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)",-1000,1000); f_land2R->SetLineStyle(2); f_land2R->SetLineColor(kBlue);
    TString formula("");
    // combined functions: [0]: peak height, [1]: peak mean, [2]: peak sigma, [3]: base height (rel), [4]: base mean (rel) (positive means innerward), [5]: base sigma (rel)
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]/2*[0]*exp(-0.5*((x-[1]-[4]*[2])/[5]/[2])*((x-[1]-[4]*[2])/[5]/[2]))"; // base gaussian
    f_combDoubleGausL = CommonTools::TF1New("f_combDoubleGausL",formula,-1000,1000); f_combDoubleGausL->SetLineStyle(1); f_combDoubleGausL->SetLineColor(kMagenta);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]/2*[0]*exp(-0.5*((x-[1]+[4]*[2])/[5]/[2])*((x-[1]+[4]*[2])/[5]/[2]))"; // base gaussian
    f_combDoubleGausR = CommonTools::TF1New("f_combDoubleGausR",formula,-1000,1000); f_combDoubleGausR->SetLineStyle(1); f_combDoubleGausR->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combDoubleLandL = CommonTools::TF1New("f_combDoubleLandL",formula,-1000,1000); f_combDoubleLandL->SetLineStyle(1); f_combDoubleLandL->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)"; // peak landau
    formula += "+[3]*[0]/0.18065564*TMath::Landau(-x,-[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combDoubleLandR = CommonTools::TF1New("f_combDoubleLandR",formula,-1000,1000); f_combDoubleLandR->SetLineStyle(1); f_combDoubleLandR->SetLineColor(kMagenta);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combGausLandL = CommonTools::TF1New("f_combGausLandL",formula,-1000,1000); f_combGausLandL->SetLineStyle(1); f_combGausLandL->SetLineColor(kMagenta);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]/0.18065564*TMath::Landau(-x,-[1]+[4]*[2],[2]*[5],false)"; // base landau
    f_combGausLandR = CommonTools::TF1New("f_combGausLandR",formula,-1000,1000); f_combGausLandR->SetLineStyle(1); f_combGausLandR->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]*exp(-0.5*((x-[1]-[4]*[2])/[2]/[5])*((x-[1]-[4]*[2])/[2]/[5]))"; // base gaussian
    f_combLandGausL = CommonTools::TF1New("f_combLandGausL",formula,-1000,1000); f_combLandGausL->SetLineStyle(1); f_combLandGausL->SetLineColor(kMagenta);
    formula = "[0]/0.18065564*TMath::Landau(-x,-[1],[2],false)"; // peak landau
    formula += "+[3]*[0]*exp(-0.5*((x-[1]+[4]*[2])/[2]/[5])*((x-[1]+[4]*[2])/[2]/[5]))"; // base gaussian
    f_combLandGausR = CommonTools::TF1New("f_combLandGausR",formula,-1000,1000); f_combLandGausR->SetLineStyle(1); f_combLandGausR->SetLineColor(kMagenta);
    // extend to both sides: add a mirrored part with new position and height
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*exp(-0.5*((x-[4])/[2])*((x-[4])/[2]))"; // peak gaussian
    f_gausBoth = CommonTools::TF1New("f_gausBoth",formula,-1000,1000); f_gausBoth->SetLineStyle(1); f_gausBoth->SetLineColor(kBlack);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]/0.18065564*TMath::Landau(-x,-[4],[2],false)"; // peak landau
    f_landBoth = CommonTools::TF1New("f_landBoth",formula,-1000,1000); f_landBoth->SetLineStyle(1); f_landBoth->SetLineColor(kBlack);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]*exp(-0.5*((x-[1]-[4]*[2])/[5]/[2])*((x-[1]-[4]*[2])/[5]/[2]))"; // base gaussian
    formula += "+[6]*exp(-0.5*((x-[7])/[2])*((x-[7])/[2]))"; // peak gaussian on the other side
    formula += "+[3]*[6]*exp(-0.5*((x-[7]+[4]*[2])/[5]/[2])*((x-[7]+[4]*[2])/[5]/[2]))"; // base gaussian on the other side
    f_combDoubleGausBoth = CommonTools::TF1New("f_combDoubleGausBoth",formula,-1000,1000); f_combDoubleGausBoth->SetLineStyle(1); f_combDoubleGausBoth->SetLineColor(kBlack);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    formula += "+[6]/0.18065564*TMath::Landau(-x,-[7],[2],false)"; // peak landau on the other side
    formula += "+[3]*[6]/0.18065564*TMath::Landau(-x,-[7]+[4]*[2],[2]*[5],false)"; // base landau on the other side
    f_combDoubleLandBoth = CommonTools::TF1New("f_combDoubleLandBoth",formula,-1000,1000); f_combDoubleLandBoth->SetLineStyle(1); f_combDoubleLandBoth->SetLineColor(kBlack);
    formula = "[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))"; // peak gaussian
    formula += "+[3]*[0]/0.18065564*TMath::Landau(x,[1]+[4]*[2],[2]*[5],false)"; // base landau
    formula += "+[6]*exp(-0.5*((x-[7])/[2])*((x-[7])/[2]))"; // peak gaussian on the other side
    formula += "+[3]*[6]/0.18065564*TMath::Landau(-x,-[7]+[4]*[2],[2]*[5],false)"; // base landau on the other side
    f_combGausLandBoth = CommonTools::TF1New("f_combGausLandBoth",formula,-1000,1000); f_combGausLandBoth->SetLineStyle(1); f_combGausLandBoth->SetLineColor(kBlack);
    formula = "[0]/0.18065564*TMath::Landau(x,[1],[2],false)"; // peak landau
    formula += "+[3]*[0]*exp(-0.5*((x-[1]-[4]*[2])/[2]/[5])*((x-[1]-[4]*[2])/[2]/[5]))"; // base gaussian
    formula += "+[6]/0.18065564*TMath::Landau(-x,-[7],[2],false)"; // peak landau on the other side
    formula += "+[3]*[6]*exp(-0.5*((x-[7]+[4]*[2])/[2]/[5])*((x-[7]+[4]*[2])/[2]/[5]))"; // base gaussian on the other side
    f_combLandGausBoth = CommonTools::TF1New("f_combLandGausBoth",formula,-1000,1000); f_combLandGausBoth->SetLineStyle(1); f_combLandGausBoth->SetLineColor(kBlack);
}

HistogramAnalyzer::~HistogramAnalyzer(void){
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


void HistogramAnalyzer::SetFittingParameters(int iRange){
    if (iRange<0||iRange>NRANGES){
        MyError("iRange "<<iRange<<" should be within 0~"<<NRANGES<<"! Will use range 0 instead!");
        iRange = 0;
    }
    par_peak_height_middle = ParameterManager::Get().HistogramAnalyzerParameters.peak_height_middle[iRange];
    par_peak_height_left = ParameterManager::Get().HistogramAnalyzerParameters.peak_height_left[iRange];
    par_peak_height_right = ParameterManager::Get().HistogramAnalyzerParameters.peak_height_right[iRange];
    par_peak_sigma_middle = ParameterManager::Get().HistogramAnalyzerParameters.peak_sigma_middle[iRange];
    par_peak_sigma_left = ParameterManager::Get().HistogramAnalyzerParameters.peak_sigma_left[iRange];
    par_peak_sigma_right = ParameterManager::Get().HistogramAnalyzerParameters.peak_sigma_right[iRange];
    par_peak_mean_range = ParameterManager::Get().HistogramAnalyzerParameters.peak_mean_range[iRange];
    par_base_height_middle = ParameterManager::Get().HistogramAnalyzerParameters.base_height_middle[iRange];
    par_base_height_left = ParameterManager::Get().HistogramAnalyzerParameters.base_height_left[iRange];
    par_base_height_right = ParameterManager::Get().HistogramAnalyzerParameters.base_height_right[iRange];
    par_base_sigma_middle = ParameterManager::Get().HistogramAnalyzerParameters.base_sigma_middle[iRange];
    par_base_sigma_left = ParameterManager::Get().HistogramAnalyzerParameters.base_sigma_left[iRange];
    par_base_sigma_right = ParameterManager::Get().HistogramAnalyzerParameters.base_sigma_right[iRange];
    par_base_mean_range = ParameterManager::Get().HistogramAnalyzerParameters.base_mean_range[iRange];
    par_functionType = ParameterManager::Get().HistogramAnalyzerParameters.functionType[iRange];
}

int HistogramAnalyzer::FitSlice(TH1D * hist, double & chi2,double & pValue, bool fitBoth, bool isLeft){
    MyNamedVerbose("HistogramAnalyzer",Form("fitSliceSingleSide %s side of histogram %s",hist->GetName(),isLeft?"left":"right"));
    // do the fitting
    int functionType = par_functionType;
    if (par_functionType==kOptimal){// loop in all functions and choose the best one
        MyNamedVerbose("HistogramAnalyzer",Form("par_functionType is kOptimal. Try all functions"));
        double pValueMax = 0;
        for (int tryFunctionType = kBasic+1; tryFunctionType<kBothSides; tryFunctionType++){
            if (tryFunctionType == kMultiple) continue;
            MyNamedVerbose("HistogramAnalyzer","  tryFunctionType = "<<tryFunctionType);
            int result = -1;
            if (chooseCurrentFunctions(tryFunctionType,fitBoth,isLeft)) result = doFitSlice(hist);
            if (result==0&&pValueMax<fit_pValue){
                MyNamedVerbose("HistogramAnalyzer","    Good result: pValue "<<pValueMax<<" -> "<<fit_pValue);
                pValueMax = fit_pValue;
                functionType = tryFunctionType;
                fit_result = result;
            }
        }
        if (pValueMax==0){// didn't find any good fitting, so just use gaussian
            fit_result = -1;
            functionType = kGaussian;
        }
        MyNamedVerbose("HistogramAnalyzer","Choose functionType = "<<functionType);
        chooseCurrentFunctions(functionType,fitBoth,isLeft);
    }
    else{
        MyNamedVerbose("HistogramAnalyzer","Choose functionType = "<<functionType);
        fit_result = -1;
        if (chooseCurrentFunctions(functionType,fitBoth,isLeft)) fit_result = doFitSlice(hist);
    }
    // set the fit_* (except fit_result which is set above) according to cur_* status
    // also set up function components to be used in drawing
    setFittingResults(hist,isLeft);
    chi2 = fit_chi2;
    pValue = fit_pValue;

    return fit_result;
}

void HistogramAnalyzer::DrawFitting(TH1D* hist){
    if (!hist){
        fprintf(stderr,"ERROR: in drawFitting, input histogram does not exist!\n");
        return;
    }
    int oldStyle = gStyle->GetOptStat();
    gStyle->SetOptStat(1);
    hist->Draw();
    double max = hist->GetMaximum();
    TLine * l_center1= new TLine(fit_x,0,fit_x,max);
    l_center1->SetLineStyle(2); l_center1->SetLineColor(kRed);
    l_center1->Draw();

    if (cur_functionType<kBothSides){ // fit single side
        cur_function->Draw("SAME");
        if (cur_functionType>kMultiple){
            cur_peak_function->Draw("SAME");
            cur_base_function->Draw("SAME");
        }
    }
    else{
        cur_functionL->Draw("SAME");
        cur_functionR->Draw("SAME");
        if (cur_functionType>kMultipleBothSides){
            cur_peak_functionL->Draw("SAME");
            cur_base_functionL->Draw("SAME");
            cur_peak_functionR->Draw("SAME");
            cur_base_functionR->Draw("SAME");
        }
        TLine * l_center1= new TLine(fit_xL,0,fit_xL,max);
        l_center1->SetLineStyle(2); l_center1->SetLineColor(kRed);
        l_center1->Draw();
        TLine * l_center2= new TLine(fit_xR,0,fit_xR,max);
        l_center2->SetLineStyle(2); l_center2->SetLineColor(kRed);
        l_center2->Draw();
    }
    gStyle->SetOptStat(oldStyle);
}

int HistogramAnalyzer::doFitSlice(TH1D * hist){
    if (!hist){
        MyError("Invalid input histogram");
        return -1;
    }
    if (!cur_function){
        MyError("Current function not set yet");
        return -2;
    }
    // get the shape of the input histogram
    double h1 = 0;double sig1 = 0; double x1 = 0; double xerr1 = 0;
    double h2 = 0;double sig2 = 0; double x2 = 0; double xerr2 = 0;
    if (cur_functionType<kBothSides){
        CommonTools::TH1GetMaximum(hist,x1,h1,-11,11);
        sig1 = hist->GetRMS(); xerr1 = sig1;
    }
    else{
        CommonTools::TH1GetMaximum(hist,x1,h1,-11,0);
        hist->GetXaxis()->SetRangeUser(-11,0);
        sig1 = hist->GetRMS(); xerr1 = sig1;
        CommonTools::TH1GetMaximum(hist,x2,h2,0,11);
        hist->GetXaxis()->SetRangeUser(0,11);
        sig2 = hist->GetRMS(); xerr2 = sig2;
        hist->GetXaxis()->UnZoom();
    }
    // set range and set parameters of the function
    cur_function->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    cur_function->SetParameters(h1*par_peak_height_middle,x1,par_peak_sigma_middle?par_peak_sigma_middle:sig1);
    cur_function->SetParLimits(0,h1*par_peak_height_left,h1*par_peak_height_right);
    cur_function->SetParLimits(1,x1-par_peak_mean_range*sig1,x1+par_peak_mean_range*sig1);
    cur_function->SetParLimits(2,par_peak_sigma_left*sig1,par_peak_sigma_right*sig1);
    if (cur_functionType>kMultiple&&cur_functionType<kBothSides){
        cur_function->SetParameters(h1*par_peak_height_middle,x1,par_peak_sigma_middle?par_peak_sigma_middle:sig1, // peak on the left
                par_base_height_middle,0.,par_base_sigma_middle?par_base_sigma_middle:1); // base part (relative to peak)
        cur_function->SetParLimits(3,par_base_height_left,par_base_height_right);
        cur_function->SetParLimits(4,-par_base_mean_range,par_base_mean_range);
        cur_function->SetParLimits(5,par_base_sigma_left,par_base_sigma_right);
    }
    else if (cur_functionType>kBothSides&&cur_functionType<kMultipleBothSides){
        cur_function->SetParameters(h1*par_peak_height_middle,x1,par_peak_sigma_middle?par_peak_sigma_middle:sig1, // peak on the left
                h2*par_peak_height_middle,x2); // peak on the right
        cur_function->SetParLimits(3,h2*par_peak_height_left,h2*par_peak_height_right);
        cur_function->SetParLimits(4,x2-par_peak_mean_range*sig1,x2+par_peak_mean_range*sig1);
    }
    else if (cur_functionType>kMultipleBothSides){
        cur_function->SetParameters(h1*par_peak_height_middle,x1,par_peak_sigma_middle?par_peak_sigma_middle:sig1, // peak on the left
                par_base_height_middle,0.,par_base_sigma_middle?par_base_sigma_middle:1, // base part (relative to peak)
                h2*par_peak_height_middle,x2); // peak on the right
        cur_function->SetParLimits(3,par_base_height_left,par_base_height_right);
        cur_function->SetParLimits(4,-par_base_mean_range,par_base_mean_range);
        cur_function->SetParLimits(5,par_base_sigma_left,par_base_sigma_right);
        cur_function->SetParLimits(6,h2*par_peak_height_left,h2*par_peak_height_right);
        cur_function->SetParLimits(7,x2-par_peak_mean_range*sig1,x2+par_peak_mean_range*sig1);
    }
    if (cur_functionL) cur_functionL->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_functionR) cur_functionR->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_peak_function) cur_peak_function->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_base_function) cur_base_function->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_peak_functionL) cur_peak_functionL->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_base_functionL) cur_base_functionL->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_peak_functionR) cur_peak_functionR->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    if (cur_base_functionR) cur_base_functionR->SetRange(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

    int result = 0;
    for (int iTry = 0; iTry<10; iTry++){
        if (iTry&&isFittingGood(cur_function)) break;
        MyNamedVerbose("HistogramAnalyzer","Before fitting. Set parameters:");
        int nPars = cur_function->GetNpar();
        for (int i = 0; i<nPars; i++){
            double parmin,parmax;
            cur_function->GetParLimits(i,parmin,parmax);
            if (iTry!=0){
                double ratio = gRandom->Uniform(0.4,0.6);
                cur_function->SetParameter(i,parmin*ratio+parmax*(1-ratio));
            }
            MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,cur_function->GetParameter(i),parmin,parmax));
        }
        result = hist->Fit(cur_function,"qN0");
        fit_chi2 = cur_function->GetChisquare();
        fit_pValue = cur_function->GetProb();
        MyNamedVerbose("HistogramAnalyzer",Form("After fitting, chi2 = %.3e, pValue = %.3e, result %d",fit_chi2,fit_pValue,result));
        for (int i = 0; i<nPars; i++){
            MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e +- %.3e",i,cur_function->GetParameter(i),cur_function->GetParError(i)));
        }
    }

    return result;
}

void HistogramAnalyzer::setFittingResults(TH1D * hist, bool isLeft){
    // update these two values (yet again) just to make sure
    fit_chi2 = cur_function->GetChisquare();
    fit_pValue = cur_function->GetProb();
    // update other values related to function shape
    if (cur_functionType<kBothSides){ // single side
        if (cur_functionType<kMultiple){
            double height = cur_function->GetParameter(0);
            fit_x = cur_function->GetParameter(1);
            fit_xerr = cur_function->GetParError(1);
            fit_sig = cur_function->GetParameter(2);
            MyNamedVerbose("HistogramAnalyzer",Form(" Result:  mean %.3e +- %.3e sig %.3e height %.3e",fit_xL,fit_xerrL,fit_sigL,height));
        }
        else{
            double peak_height = cur_function->GetParameter(0);
            double peak_mean = cur_function->GetParameter(1); double peak_meanErr = cur_function->GetParError(1);
            double peak_sigma = cur_function->GetParameter(2);
            double base_heightRel = cur_function->GetParameter(3);
            double base_meanRel = cur_function->GetParameter(4); double base_meanRelErr = cur_function->GetParError(4);
            double base_sigmaRel = cur_function->GetParameter(5);
            double base_height = peak_height*base_heightRel;
            double base_mean = peak_mean+(isLeft?1:-1)*base_meanRel*peak_sigma;
            double base_sigma = base_sigmaRel*peak_sigma;
            // get mean error and sigma
            fit_xerr = peak_meanErr;
            getMeanRMS(cur_function,hist,fit_x,fit_sig);
            // set the single functions
            cur_peak_function->SetParameters(peak_height,peak_mean,peak_sigma);
            cur_base_function->SetParameters(base_height,base_mean,base_sigma);
            MyNamedVerbose("HistogramAnalyzer",Form(" Result:  mean %.3e +- %.3e sig %.3e",fit_xL,fit_xerrL,fit_sigL));
            MyNamedVerbose("HistogramAnalyzer",Form(" peak: height %.3e mean %.3e +- %.3e sig %.3e",peak_height,peak_mean,peak_meanErr,peak_sigma));
            MyNamedVerbose("HistogramAnalyzer",Form(" base: height %.3e mean %.3e +- %.3e sig %.3e",base_height,base_mean,base_meanRelErr*peak_sigma,base_sigma));
        }
    }
    else{ // both sides
        if (cur_functionType<kMultipleBothSides){
            // now let's get the left peak
            double height = cur_function->GetParameter(0);
            fit_xL = cur_function->GetParameter(1); fit_xerrL = cur_function->GetParError(1);
            fit_sigL = cur_function->GetParameter(2);
            cur_functionL->SetParameters(height,fit_xL,fit_sigL);
            MyNamedVerbose("HistogramAnalyzer",Form(" Left Part:  mean %.3e +- %.3e sig %.3e height %.3e",fit_xL,fit_xerrL,fit_sigL,height));
            height = cur_function->GetParameter(3);
            fit_xR = cur_function->GetParameter(4); fit_xerrR = cur_function->GetParError(1);
            fit_sigR = fit_sigL;
            cur_functionR->SetParameters(height,fit_xR,fit_sigR);
            MyNamedVerbose("HistogramAnalyzer",Form(" Right Part: mean %.3e +- %.3e sig %.3e height %.3e",fit_xR,fit_xerrR,fit_sigR,height));
        }
        else {
            // first, get the left peak
            double fit_peak_height = cur_function->GetParameter(0);
            double fit_peak_mean = cur_function->GetParameter(1); double fit_peak_meanErr = cur_function->GetParError(1);
            double fit_peak_sigma = cur_function->GetParameter(2);
            double fit_base_heightRel = cur_function->GetParameter(3);
            double fit_base_meanRel = cur_function->GetParameter(4); double fit_base_meanRelErr = cur_function->GetParError(4);
            double fit_base_sigmaRel = cur_function->GetParameter(5);
            double fit_base_height = fit_peak_height*fit_base_heightRel;
            double fit_base_mean = fit_peak_mean+fit_base_meanRel*fit_peak_sigma;
            double fit_base_sigma = fit_base_sigmaRel*fit_peak_sigma;
            // get mean error and sigma
            fit_xerrL = fit_peak_meanErr;
            cur_functionL->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma,fit_base_heightRel,fit_base_meanRel,fit_base_sigmaRel);
            getMeanRMS(cur_functionL,hist,fit_xL,fit_sigL);
            // set the single functions
            cur_peak_functionL->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma);
            cur_base_functionL->SetParameters(fit_base_height,fit_base_mean,fit_base_sigma);
            MyNamedVerbose("HistogramAnalyzer",Form(" Left Part:  mean %.3e +- %.3e sig %.3e",fit_xL,fit_xerrL,fit_sigL));
            MyNamedVerbose("HistogramAnalyzer",Form(" peakL: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
            MyNamedVerbose("HistogramAnalyzer",Form(" baseL: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));

            // second, get the left peak
            fit_peak_height = cur_function->GetParameter(6);
            fit_peak_mean = cur_function->GetParameter(7); fit_peak_meanErr = cur_function->GetParError(7);
            fit_base_height = fit_peak_height*fit_base_heightRel;
            fit_base_mean = fit_peak_mean-fit_base_meanRel*fit_peak_sigma; // flipped
            // get mean error and sigma
            fit_xerrR = fit_peak_meanErr;
            cur_functionR->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma,fit_base_heightRel,fit_base_meanRel,fit_base_sigmaRel);
            getMeanRMS(cur_functionR,hist,fit_xR,fit_sigR);
            // set the single functions
            cur_peak_functionR->SetParameters(fit_peak_height,fit_peak_mean,fit_peak_sigma);
            cur_base_functionR->SetParameters(fit_base_height,fit_base_mean,fit_base_sigma);
            MyNamedVerbose("HistogramAnalyzer",Form(" Right Part: mean %.3e +- %.3e sig %.3e",fit_xR,fit_xerrR,fit_sigR));
            MyNamedVerbose("HistogramAnalyzer",Form(" peakR: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
            MyNamedVerbose("HistogramAnalyzer",Form(" baseR: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));
        }
    }
}

bool HistogramAnalyzer::isFittingGood(TF1 * f){
    bool isGood = true;
    int nPars = f->GetNpar();
    for (int i = 0; i<nPars; i++){
        double parmin,parmax,par;
        f->GetParLimits(i,parmin,parmax);
        par = f->GetParameter(i);
        if (fabs(par-parmin)<(parmax-parmin)*1e-7||fabs(par-parmax)<(parmax-parmin)*1e-7){
            isGood = false;
        }
    }
    return isGood;
}

void HistogramAnalyzer::getMeanRMS(TF1 * f, const TH1D * hist, double & mean, double & sigma){
    mean = f->GetMaximumX();
    TH1D * histNew = new TH1D(*hist);
    histNew->Reset();
    for (int i = 0; i<1e6; i++){
        double x = f->GetRandom();
        double content = hist->GetBinContent(histNew->FindBin(x));
        if (content) histNew->Fill(x);
    }
    sigma = histNew->GetRMS();
    histNew->GetXaxis()->SetRangeUser(mean-sigma*5,mean+sigma*5);
    sigma = histNew->GetRMS();
    delete histNew;
}

bool HistogramAnalyzer::chooseCurrentFunctions(int functionType, bool fitBoth, bool isLeft){
    cur_functionType = functionType;
    if (fitBoth) cur_functionType += kBothSides+1; // this is assuming function types are arranged in the same order in single-side block and both-sides block
    MyNamedVerbose("HistogramAnalyzer","cur_functionType set to "<<cur_functionType);
    if (functionType==kGaussian||functionType==kLandau){ // single function
        if (functionType==kGaussian){ // fit with gaussian
            if (fitBoth){
                cur_function = f_gausBoth;
                cur_functionL = f_gausL;
                cur_functionR = f_gausR;
            }
            else{
                if (isLeft) {cur_function = f_gausL;} else {cur_function = f_gausR;}
            }
        }
        else if (functionType==kLandau){ // fit with landau
            if (fitBoth){
                cur_function = f_landBoth;
                cur_functionL = f_landL;
                cur_functionR = f_landR;
            }
            else{
            }
            if (isLeft) {cur_function = f_landL;} else {cur_function = f_landR;}
        }
    }
    else{ // composite function
        if (functionType==kGaussianPlusLandau){ // fit with gaussian plus landau
            if (fitBoth){
                cur_function = f_combGausLandBoth;
                cur_functionL = f_combGausLandL;
                cur_functionR = f_combGausLandR;
                cur_peak_functionL = f_gausL;
                cur_base_functionL = f_land2L;
                cur_peak_functionR = f_gausR;
                cur_base_functionR = f_land2R;
            }
            else{
                if (isLeft){
                    cur_function = f_combGausLandL;
                    cur_peak_function = f_gausL;
                    cur_base_function = f_land2L;
                }
                else{
                    cur_function = f_combGausLandR;
                    cur_peak_function = f_gausR;
                    cur_base_function = f_land2R;
                }
            }
        }
        else if (functionType==kLandauPlusGaussian){ // fit with gaussian plus landau
            if (fitBoth){
                cur_function = f_combLandGausBoth;
                cur_functionL = f_combLandGausL;
                cur_functionR = f_combLandGausR;
                cur_peak_functionL = f_landL;
                cur_base_functionL = f_gaus2L;
                cur_peak_functionR = f_landR;
                cur_base_functionR = f_gaus2R;
            }
            else{
                if (isLeft){
                    cur_function = f_combLandGausL;
                    cur_peak_function = f_landL;
                    cur_base_function = f_gaus2L;
                }
                else{
                    cur_function = f_combLandGausR;
                    cur_peak_function = f_landR;
                    cur_base_function = f_gaus2R;
                }
            }
        }
        else if (functionType==kDoubleGaussian){ // fit with double gaussian
            if (fitBoth){
                cur_function = f_combDoubleGausBoth;
                cur_functionL = f_combDoubleGausL;
                cur_functionR = f_combDoubleGausR;
                cur_peak_functionL = f_gausL;
                cur_base_functionL = f_gaus2L;
                cur_peak_functionR = f_gausR;
                cur_base_functionR = f_gaus2R;
            }
            else{
                if (isLeft){
                    cur_function = f_combDoubleGausL;
                    cur_peak_function = f_gausL;
                    cur_base_function = f_gaus2L;
                }
                else{
                    cur_function = f_combDoubleGausR;
                    cur_peak_function = f_gausR;
                    cur_base_function = f_gaus2R;
                }
            }
        }
        else if (functionType==kDoubleLandau){ // fit with double gaussian
            if (fitBoth){
                cur_function = f_combDoubleLandBoth;
                cur_functionL = f_combDoubleLandL;
                cur_functionR = f_combDoubleLandR;
                cur_peak_functionL = f_landL;
                cur_base_functionL = f_land2L;
                cur_peak_functionR = f_landR;
                cur_base_functionR = f_land2R;
            }
            else{
                if (isLeft){
                    cur_function = f_combDoubleLandL;
                    cur_peak_function = f_landL;
                    cur_base_function = f_land2L;
                }
                else{
                    cur_function = f_combDoubleLandR;
                    cur_peak_function = f_landR;
                    cur_base_function = f_land2R;
                }
            }
        }
        else{
            MyError("Unknown functionType "<<functionType);
            return false;
        }
    }
    return true;
}
