
#ifndef CommonTools_hxx_seen
#define CommonTools_hxx_seen

#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"

class CommonTools{
public:
    CommonTools(){};
    virtual ~CommonTools(){};

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
    }

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
    }

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
            // theY = fabs(theX-prevX)>fabs(theX-firstX)?firstY:prevY;
            return false;
        }
        return true;
    }
};

#endif
