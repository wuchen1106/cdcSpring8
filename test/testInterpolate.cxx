#include "TCanvas.h"
#include "TGraph.h"

#include "CommonTools.hxx"

int main(){
    TGraph * grin = new TGraph(5);
    for (int i = 0; i<5; i++){
        grin->SetPoint(i,i,i);
    }
    TGraph * grout = new TGraph(50);
    for (int i = 0; i<50; i++){
        double x = i/10.;
        double y;
        CommonTools::TGraphInterpolate(grin,x,y);
        grout->SetPoint(i,x,y);
    }
    TCanvas * canv = new TCanvas();
    grin->SetMarkerStyle(20);
    grin->Draw("APL");
    grout->SetMarkerStyle(21);
    grout->SetMarkerColor(kRed);
    grout->SetLineColor(kRed);
    grout->Draw("PLSAME");
    canv->SaveAs("testInterpolate.png");
}
