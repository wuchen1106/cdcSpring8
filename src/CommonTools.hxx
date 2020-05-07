
#ifndef CommonTools_hxx_seen
#define CommonTools_hxx_seen

#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"

class CommonTools{
public:
    CommonTools(){};
    virtual ~CommonTools(){};

    struct HistogramData{
        double Mean;
        double MaxHeight;
        double RMS;
        double MPV;
        double FWHM;
        double left;
        double right;
        double cutRatio;
    };

    static HistogramData TH1Analyze(TH1 * h, bool ignoreFirstBin = false, double cutRatio = 0.5){
        HistogramData data;
        data.Mean = h->GetMean();
        data.RMS = h->GetRMS();
        data.cutRatio = cutRatio;
        if (ignoreFirstBin){
            h->GetXaxis()->SetRangeUser(h->GetBinLowEdge(2),h->GetBinLowEdge(h->GetNbinsX()+1));
        }
        data.MaxHeight = h->GetBinContent(h->GetMaximumBin());
        data.MPV = h->GetBinCenter(h->GetMaximumBin());
        bool passedMPV = false;
        bool foundLeft = false;
        bool foundRight = false;
        for (int i = (ignoreFirstBin?2:1); i<=h->GetNbinsX(); i++){
            double c = h->GetBinContent(i);
            double x = h->GetBinCenter(i);
            if (!foundLeft&&!passedMPV){
                if (c>data.MaxHeight*cutRatio){
                    data.left = x;
                    foundLeft = true;
                }
            }
            if (!passedMPV&&x>=data.MPV) passedMPV = true;
            if (passedMPV&&!foundRight){
                if (c<data.MaxHeight*cutRatio){
                    data.right = x;
                    foundRight = true;
                    break;
                }
            }
        }
        if (!foundLeft) data.left = h->GetBinLowEdge(h->GetNbinsX()+1);
        if (!foundRight) data.right = h->GetBinLowEdge(1);
        data.FWHM = data.right-data.left;
        return data;
    };

    static TF1 * TF1New(TString name, TString form, double left, double right){
        TF1 * f = new TF1(name,form,left,right);
        f->SetNpx(1024);
        f->SetNumberFitPoints(1024);
        f->SetLineWidth(1);
        f->SetLineColor(kRed);
        return f;
    };

    static void TH1GetMaximum(TH1D * hist, double & position, double & maximum, double left, double right){
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
    };

    static void TGraphErrorsPlus(TGraphErrors * gn, const TGraphErrors * gl, const TGraphErrors * gr, double sl, double sr){
        gn->Set(gl->GetN());
        int count = 0;
        for (int iPoint = 0; iPoint<gl->GetN(); iPoint++){
            double xl,yl; double xerr,yerr;
            gl->GetPoint(iPoint,xl,yl);
            yl*=sl;
            xerr = gl->GetErrorX(iPoint);
            yerr = gl->GetErrorY(iPoint);
            double yr = 0;
            if(!TGraphInterpolate(gr,xl,yr)){
                continue;
            }
            yr*=sr;
            gn->SetPoint(count,xl,yl+yr);
            gn->SetPointError(count,xerr,yerr);
            count++;
        }
        gn->Set(count);
    };

    static bool TGraphInterpolate(const TGraph* graph, double theX, double & theY){
        theY = 0;
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
                theY = (prevY*(x-theX)+y*(theX-prevX))/(x-prevX);
                found = true;
                break;
            }
            prevX = x;
            prevY = y;
        }
        if (!found){ // didn't cross the asked position
            // DON'T EXTRAPOLATE
            theY = fabs(theX-prevX)>fabs(theX-firstX)?firstY:prevY;
            return false;
        }
        return true;
    };

    static bool TGraphGetX(const TGraph* graph, double theY, double & theX){
        // check the first point
        double x,y;
        graph->GetPoint(0,x,y);
        bool isLeft = y>theY; // if the asked position is on the left side of the first point.
        double firstX = x; double firstY = y;
        double prevX = x; double prevY = y;
        bool found = false;
        for (int i = 1; i<graph->GetN(); i++){
            graph->GetPoint(i,x,y);
            if ((isLeft&&y<theY)||(!isLeft&&y>theY)){ // moved to the other side of the asked position, then interpolate
                theX = (prevX*(y-theY)+x*(theY-prevY))/(y-prevY);
                //std::cout<<i<<": "<<"("<<prevX<<"*("<<y<<"-"<<theY<<")+"<<x<<"*("<<theY<<"-"<<prevY<<"))/("<<y<<"-"<<prevY<<") = "<<theX<<std::endl;
                found = true;
                break;
            }
            prevX = x;
            prevY = y;
        }
        if (!found){ // didn't cross the asked position
            // DON'T EXTRAPOLATE
            theX = fabs(theY-prevY)>fabs(theY-firstY)?firstX:prevX;
            return false;
        }
        return true;
    };

    static void TGraphErrorsSortByX(TGraphErrors * graph){
        for (int i = 0; i<graph->GetN(); i++){
            for (int j = i; j<graph->GetN(); j++){
                double xi,yi,xerri,yerri;
                graph->GetPoint(i,xi,yi);
                xerri = graph->GetErrorX(i);
                yerri = graph->GetErrorY(i);
                double xj,yj,xerrj,yerrj;
                graph->GetPoint(j,xj,yj);
                xerrj = graph->GetErrorX(j);
                yerrj = graph->GetErrorY(j);
                if (xi>xj){
                    graph->SetPoint(i,xj,yj);
                    graph->SetPointError(i,xerrj,yerrj);
                    graph->SetPoint(j,xi,yi);
                    graph->SetPointError(j,xerri,yerri);
                }
            }
        }
    };

    static void TGraphGetPol1(const TGraph * graph, double & p0, double & p1, double left = 0, double right = 0){
        double x1 = 1e9;
        double y1 = 0;
        double x2 = -1e9;
        double y2 = 0;
        for (int i = 0; i<graph->GetN(); i++){
            double x,y; graph->GetPoint(i,x,y);
            if (left==right||x<x1){
                x1 = x;
                y1 = y;
            }
            if (left==right||x>x2){
                x2 = x;
                y2 = y;
            }
        }
        if (x1==x2){
            p0 = 0;
            p1 = 0;
        }
        else{
            p1 = (y2-y1)/(x2-x1);
            p0 = y1-p1*x1;
        }
    };

    static double TGraphDerivative(const TGraph * graph, double theX){
        double deriv = 0;
        // check the first point
        double x,y;
        graph->GetPoint(0,x,y);
        bool isLeft = x>theX; // if the asked position is on the left side of the first point.
        double firstX = x; double firstY = y;
        int N = graph->GetN();
        for (int i = 0; i<N; i++){
            graph->GetPoint(i,x,y);
            if ((isLeft&&x<theX)||(!isLeft&&x>theX)){ // moved to the other side of the asked position, then interpolate
                double x2,y2;
                if (i==N-1) graph->GetPoint(i-1,x2,y2);
                else graph->GetPoint(i+1,x2,y2);
                if (y2!=y) deriv = (y2-y)/(x2-x);
                break;
            }
        }
        return deriv;
    };

    static TString f2s(double value){
        int nDecimals = 0;
        for (int i = 0; i<7; i++){
            double residual = value - std::round(value*pow(10,i))/pow(10,i);
            if (fabs(residual) < 1e-7){
                nDecimals = i;
                break;
            }
        }
        TString reg("%.");
        reg+=Form("%df",nDecimals);
        TString valueStr = Form(reg.Data(),value);
        return valueStr;
    }

};

#endif
