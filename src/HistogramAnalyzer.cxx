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
    f_gausL(NULL), f_gausR(NULL), f_gaus2L(NULL), f_gaus2R(NULL), f_gausBoth(NULL), f_landL(NULL), f_landR(NULL), f_land2L(NULL), f_land2R(NULL), f_landBoth(NULL),
    f_combDoubleGausL(NULL), f_combDoubleGausR(NULL), f_combDoubleGausBoth(NULL),
    f_combDoubleLandL(NULL), f_combDoubleLandR(NULL), f_combDoubleLandBoth(NULL),
    f_combGausLandL(NULL), f_combGausLandR(NULL), f_combGausLandBoth(NULL),
    f_combLandGausL(NULL), f_combLandGausR(NULL), f_combLandGausBoth(NULL)
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

TF1 * HistogramAnalyzer::FitSliceSingleSide(TH1D * hist, double & x1,double & xerr1,double & sig1,double & chi2,double & prob,int & result, int & functionType, int iRange, bool isLeft){
    MyNamedVerbose("HistogramAnalyzer",Form("fitSliceSingleSide %s side of histogram %s",hist->GetName(),isLeft?"left":"right"));
    // get the shape of the input histogram
    double h1 = 0;
    CommonTools::TH1GetMaximum(hist,x1,h1,-11,11);
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
    if (functionType==HistogramAnalyzer::kGaussian||functionType==HistogramAnalyzer::kLandau){ // single function
        if (functionType==HistogramAnalyzer::kGaussian){ // fit with gaussian
            if (isLeft) {f = f_gausL;} else {f = f_gausR;}
        }
        else if (functionType==HistogramAnalyzer::kLandau){ // fit with landau
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
        if (functionType==HistogramAnalyzer::kGaussianPlusLandau||functionType==HistogramAnalyzer::kOptimal){ // fit with gaussian plus landau
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
        else if (functionType==HistogramAnalyzer::kLandauPlusGaussian){ // fit with gaussian plus landau
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
        else if (functionType==HistogramAnalyzer::kDoubleGaussian){ // fit with double gaussian
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
        else if (functionType==HistogramAnalyzer::kDoubleLandau){ // fit with double gaussian
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
        if (iTry&&isFittingGood(f)) break;
        MyNamedVerbose("HistogramAnalyzer","Before fitting. Set parameters:");
        for (int i = 0; i<nPars; i++){
            double parmin,parmax;
            f->GetParLimits(i,parmin,parmax);
            if (iTry!=0){
                double ratio = gRandom->Uniform(0.4,0.6);
                f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
            }
            MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
        }
        result = hist->Fit(f,"qN0");
        chi2 = f->GetChisquare();
        prob = f->GetProb();
        MyNamedVerbose("HistogramAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2,prob,result));
        for (int i = 0; i<nPars; i++){
            MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
        }
    }

    if (functionType==HistogramAnalyzer::kOptimal){ // try LandauPlusGaussian in this case
        MyNamedVerbose("HistogramAnalyzer",Form("fitFunctionType is kOptimal. Try LandauPlusGaussian"));
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
            if (iTry&&isFittingGood(f)) break;
            MyNamedVerbose("HistogramAnalyzer","Before fitting. Set parameters:");
            for (int i = 0; i<nPars; i++){
                double parmin,parmax;
                f->GetParLimits(i,parmin,parmax);
                if (iTry!=0){
                    double ratio = gRandom->Uniform(0.4,0.6);
                    f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
                }
                MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
            }
            resultTemp = hist->Fit(f,"qN0");
            chi2Temp = f->GetChisquare();
            probTemp = f->GetProb();
            MyNamedVerbose("HistogramAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2Temp,probTemp,resultTemp));
            for (int i = 0; i<nPars; i++){
                MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
            }
        }
        if (chi2Temp<chi2){
            functionType = HistogramAnalyzer::kLandauPlusGaussian;
            chi2 = chi2Temp;
            prob = probTemp;
            result = resultTemp;
        }
        else{ // change the function back
            functionType = HistogramAnalyzer::kGaussianPlusLandau;
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

    if (functionType==HistogramAnalyzer::kGaussian||functionType==HistogramAnalyzer::kLandau){
        double height = f->GetParameter(0);
        x1 = f->GetParameter(1); xerr1 = f->GetParError(1);
        sig1 = f->GetParameter(2);
        MyNamedVerbose("HistogramAnalyzer",Form(" peak: height %.3e mean %.3e sig %.3e",height,x1,sig1));
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
        MyNamedVerbose("HistogramAnalyzer",Form(" peak: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
        MyNamedVerbose("HistogramAnalyzer",Form(" base: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));
    }

    MyNamedVerbose("HistogramAnalyzer",Form(" Result:  mean %.3e +- %.3e sig %.3e",x1,xerr1,sig1));

    return f;
}

TF1 * HistogramAnalyzer::FitSliceBothSides(TH1D * hist, double & x1,double & xerr1,double & sig1,double & x2,double & xerr2,double & sig2,double & chi2,double & prob, int & result, int & functionType, int iRange){
    MyNamedVerbose("HistogramAnalyzer",Form("fitSliceBothSides of histogram %s",hist->GetName()));
    // get the shape of the input histogram
    double h1 = 0;
    CommonTools::TH1GetMaximum(hist,x1,h1,-11,0);
    hist->GetXaxis()->SetRangeUser(-11,0);
    sig1 = hist->GetRMS(); xerr1 = sig1;
    double h2 = 0;
    CommonTools::TH1GetMaximum(hist,x2,h2,0,11);
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
    if (fitFunctionType==HistogramAnalyzer::kGaussian||fitFunctionType==HistogramAnalyzer::kLandau){ // single function
        if (fitFunctionType==HistogramAnalyzer::kGaussian){ // fit with gaussian
            functionType = HistogramAnalyzer::kGaussianBothSides;
            f = f_gausBoth;
            fl = f_gausL;
            fr = f_gausR;
        }
        else if (fitFunctionType==HistogramAnalyzer::kLandau){ // fit with landau
            functionType = HistogramAnalyzer::kLandauBothSides;;
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
        if (fitFunctionType==HistogramAnalyzer::kGaussianPlusLandau||fitFunctionType==HistogramAnalyzer::kOptimal){ // fit with gaussian plus landau
            functionType = HistogramAnalyzer::kGaussianPlusLandauBothSides;
            f = f_combGausLandBoth;
            fpl = f_gausL;
            fpr = f_gausR;
            fbl = f_land2L;
            fbr = f_land2R;
            fl = f_combGausLandL;
            fr = f_combGausLandR;
        }
        else if (fitFunctionType==HistogramAnalyzer::kLandauPlusGaussian){ // fit with gaussian plus landau
            functionType = HistogramAnalyzer::kLandauPlusGaussianBothSides;
            f = f_combLandGausBoth;
            fpl = f_landL;
            fpr = f_landR;
            fbl = f_gaus2L;
            fbr = f_gaus2R;
            fl = f_combLandGausL;
            fr = f_combLandGausR;
        }
        else if (fitFunctionType==HistogramAnalyzer::kDoubleGaussian){ // fit with double gaussian
            functionType = HistogramAnalyzer::kDoubleGaussianBothSides;;
            f = f_combDoubleGausBoth;
            fpl = f_gausL;
            fpr = f_gausR;
            fbl = f_gaus2L;
            fbr = f_gaus2R;
            fl = f_combDoubleGausL;
            fr = f_combDoubleGausR;
        }
        else if (fitFunctionType==HistogramAnalyzer::kDoubleLandau){ // fit with double gaussian
            functionType = HistogramAnalyzer::kDoubleLandauBothSides;;
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
        if (iTry&&isFittingGood(f)) break;
        MyNamedVerbose("HistogramAnalyzer","Before fitting. Set parameters:");
        for (int i = 0; i<nPars; i++){
            double parmin,parmax;
            f->GetParLimits(i,parmin,parmax);
            if (iTry!=0){
                double ratio = gRandom->Uniform(0.4,0.6);
                f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
            }
            MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
        }
        result = hist->Fit(f,"qN0");
        chi2 = f->GetChisquare();
        prob = f->GetProb();
        MyNamedVerbose("HistogramAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2,prob,result));
        for (int i = 0; i<nPars; i++){
            MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
        }
    }

    if (fitFunctionType==HistogramAnalyzer::kOptimal){ // try LandauPlusGaussian in this case
        MyNamedVerbose("HistogramAnalyzer",Form("fitFunctionType is kOptimal. Try LandauPlusGaussian"));
        functionType = HistogramAnalyzer::kLandauPlusGaussianBothSides;
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
            if (iTry&&isFittingGood(f)) break;
            MyNamedVerbose("HistogramAnalyzer","Before fitting. Set parameters:");
            for (int i = 0; i<nPars; i++){
                double parmin,parmax;
                f->GetParLimits(i,parmin,parmax);
                if (iTry!=0){
                    double ratio = gRandom->Uniform(0.4,0.6);
                    f->SetParameter(i,parmin*ratio+parmax*(1-ratio));
                }
                MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e (%.3e ~ %.3e)",i,f->GetParameter(i),parmin,parmax));
            }
            result = hist->Fit(f,"qN0");
            chi2Temp = f->GetChisquare();
            probTemp = f->GetProb();
            MyNamedVerbose("HistogramAnalyzer",Form("After fitting, chi2 = %.3e, prob = %.3e, result %d",chi2Temp,probTemp,resultTemp));
            for (int i = 0; i<nPars; i++){
                MyNamedVerbose("HistogramAnalyzer",Form("  %d: %.3e +- %.3e",i,f->GetParameter(i),f->GetParError(i)));
            }
        }
        if (chi2Temp<chi2){
            chi2 = chi2Temp;
            prob = probTemp;
            result = resultTemp;
        }
        else{ // change the function back
            functionType = HistogramAnalyzer::kGaussianPlusLandauBothSides;
            f = f_combGausLandBoth;
            fpl = f_gausL;
            fpr = f_gausR;
            fbl = f_land2L;
            fbr = f_land2R;
            fl = f_combGausLandL;
            fr = f_combGausLandR;
        }
    }

    if (fitFunctionType==HistogramAnalyzer::kGaussian||fitFunctionType==HistogramAnalyzer::kLandau){
        // now let's get the left peak
        double height = f->GetParameter(0);
        x1 = f->GetParameter(1); xerr1 = f->GetParError(1);
        sig1 = f->GetParameter(2);
        fl->SetParameters(height,x1,sig1);
        MyNamedVerbose("HistogramAnalyzer",Form(" peakL: height %.3e mean %.3e sig %.3e",height,x1,sig1));
        height = f->GetParameter(3);
        x2 = f->GetParameter(4); xerr2 = f->GetParError(1);
        sig2 = sig1;
        fr->SetParameters(height,x2,sig2);
        MyNamedVerbose("HistogramAnalyzer",Form(" peakR: height %.3e mean %.3e sig %.3e",height,x2,sig2));
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
        MyNamedVerbose("HistogramAnalyzer",Form(" peakL: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
        MyNamedVerbose("HistogramAnalyzer",Form(" baseL: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));

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
        MyNamedVerbose("HistogramAnalyzer",Form(" peakR: height %.3e mean %.3e +- %.3e sig %.3e",fit_peak_height,fit_peak_mean,fit_peak_meanErr,fit_peak_sigma));
        MyNamedVerbose("HistogramAnalyzer",Form(" baseR: height %.3e mean %.3e +- %.3e sig %.3e",fit_base_height,fit_base_mean,fit_base_meanRelErr*fit_peak_sigma,fit_base_sigma));
    }

    MyNamedVerbose("HistogramAnalyzer",Form(" Left Part:  mean %.3e +- %.3e sig %.3e",x1,xerr1,sig1));
    MyNamedVerbose("HistogramAnalyzer",Form(" Right Part: mean %.3e +- %.3e sig %.3e",x2,xerr2,sig2));

    return f;
}

void HistogramAnalyzer::DrawFitting(TH1D* hist, TCanvas * c,TString title, TString filename, int function, double center1, double center2, bool isLeft){
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
    if (function>=HistogramAnalyzer::kGaussianBothSides){ // fit both sides
        l_center2->Draw();
    }
    if (function==HistogramAnalyzer::kGaussianPlusLandauBothSides){
        f_combGausLandBoth->Draw("SAME");
        f_gausL->Draw("SAME");
        f_land2L->Draw("SAME");
        f_combGausLandL->Draw("SAME");
        f_gausR->Draw("SAME");
        f_land2R->Draw("SAME");
        f_combGausLandR->Draw("SAME");
    }
    else if (function==HistogramAnalyzer::kLandauPlusGaussianBothSides){
        f_combLandGausBoth->Draw("SAME");
        f_landL->Draw("SAME");
        f_gaus2L->Draw("SAME");
        f_combLandGausL->Draw("SAME");
        f_landR->Draw("SAME");
        f_gaus2R->Draw("SAME");
        f_combLandGausR->Draw("SAME");
    }
    else if (function==HistogramAnalyzer::kDoubleGaussianBothSides){
        f_combDoubleGausBoth->Draw("SAME");
        f_gausL->Draw("SAME");
        f_gaus2L->Draw("SAME");
        f_combDoubleGausL->Draw("SAME");
        f_gausR->Draw("SAME");
        f_gaus2R->Draw("SAME");
        f_combDoubleGausR->Draw("SAME");
    }
    else if (function==HistogramAnalyzer::kDoubleLandauBothSides){
        f_combDoubleLandBoth->Draw("SAME");
        f_landL->Draw("SAME");
        f_land2L->Draw("SAME");
        f_combDoubleLandL->Draw("SAME");
        f_landR->Draw("SAME");
        f_land2R->Draw("SAME");
        f_combDoubleLandR->Draw("SAME");
    }
    else if (function==HistogramAnalyzer::kGaussianBothSides){
        f_gausBoth->Draw("SAME");
        f_gausL->Draw("SAME");
        f_gausR->Draw("SAME");
    }
    else if (function==HistogramAnalyzer::kLandauBothSides){
        f_landBoth->Draw("SAME");
        f_landL->Draw("SAME");
        f_landR->Draw("SAME");
    }
    else if (function==HistogramAnalyzer::kGaussianPlusLandau){
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
    else if (function==HistogramAnalyzer::kLandauPlusGaussian){
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
    else if (function==HistogramAnalyzer::kDoubleGaussian){
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
    else if (function==HistogramAnalyzer::kDoubleLandau){
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
    else if (function==HistogramAnalyzer::kGaussian){
        if (isLeft){
            f_gausL->Draw("SAME");
        }
        else{
            f_gausR->Draw("SAME");
        }
    }
    else if (function==HistogramAnalyzer::kLandau){
        if (isLeft){
            f_landL->Draw("SAME");
        }
        else{
            f_landR->Draw("SAME");
        }
    }
    if ( filename != "" ) c->SaveAs(filename);
    gStyle->SetOptStat(oldStyle);
    delete l_center1;
    delete l_center2;
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

