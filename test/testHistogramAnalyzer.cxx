#include <stdio.h>

#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "HistogramAnalyzer.hxx"

int main(int argc, char **argv){

    TH1D * hist = 0;
    TF1  * f = 0;
    double x,xerr,sig,chi2,prob;
    int    result;
    bool   isLeft = false;
    // set the landau ratio
    double ratio[10] = {0,0.01,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45};

    HistogramAnalyzer::Get().SetFittingParameters(0);
    for (int i = 1; i<=9; i++){
        f = new TF1(Form("fgausland%d",i),"gaus(0)+landau(3)",-10,10);
        f->SetParameters(1,0,1,ratio[i],0,2);
        f->SetLineColor(kBlack); f->SetLineStyle(2);
        hist = new TH1D(Form("h%d",i),"",1024,-10,10);
        for (int i = 0; i<1e6; i++){
            hist->Fill(f->GetRandom());
        }
        isLeft = true;
        result = HistogramAnalyzer::Get().FitSlice(hist,chi2,prob,false,isLeft);
        x = HistogramAnalyzer::Get().get_x();
        xerr = HistogramAnalyzer::Get().get_xerrR();
        sig = HistogramAnalyzer::Get().get_sig();
        std::cout<<"Test "<<i<<": gaus+land("<<f->GetParameter(0)<<","<<f->GetParameter(1)<<","<<f->GetParameter(2)<<","<<f->GetParameter(3)<<","<<f->GetParameter(4)<<","<<f->GetParameter(5)<<")"<<std::endl;
        std::cout<<"  After fitting:"<<std::endl;
        std::cout<<"    x = "<<x<<std::endl;
        std::cout<<"    xerr = "<<xerr<<std::endl;
        std::cout<<"    sig = "<<sig<<std::endl;
        std::cout<<"    chi2 = "<<chi2<<std::endl;
        std::cout<<"    prob = "<<prob<<std::endl;
        std::cout<<"    result = "<<result<<std::endl;
        TCanvas * canv = new TCanvas(Form("canv%d",i),"",1024,1024);
        hist->Draw();
        HistogramAnalyzer::Get().DrawFitting(hist);
        canv->SaveAs(Form("test%d.png",i));
    }

    return 0;
}
