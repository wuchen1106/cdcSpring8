#include <stdio.h>  /* printf, getenv */
#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TGraph.h>

#include "Log.hxx"
#include "XTManager.hxx"
#include "RunInfoManager.hxx"
#include "ParameterManager.hxx"

XTManager* XTManager::fXTManager = NULL;

XTManager::XTManager():
    useSideXT(false),
    fInputFileXT(NULL),
    fInputFileRes(NULL),
    fXTLeftEven(NULL),
    fXTRightEven(NULL),
    fXTLeftOdd(NULL),
    fXTRightOdd(NULL),
    fXTLeftDefault(NULL),
    fXTRightDefault(NULL),
    fXTHistEven(NULL),
    fXTHistOdd(NULL),
    fXTHistDefault(NULL),
    fResIntrinsic(NULL),
    fResIntrinsicFunction(NULL),
    xtType(kSingleFolded)
{
    for (int i = 0; i<NLAY; i++){
        fXTLeft[i] = NULL; 
        fXTRight[i] = NULL;
        fXTHist[i] = NULL;
    }
}

XTManager::~XTManager(){
    if (fInputFileXT) fInputFileXT->Close();
    if (fInputFileRes) fInputFileRes->Close();
    // FIXME: will there be double delete?
    if (fXTLeftEven) delete fXTLeftEven;
    if (fXTRightEven) delete fXTRightEven;
    if (fXTLeftOdd) delete fXTLeftOdd;
    if (fXTRightOdd) delete fXTRightOdd;
    if (fXTLeftDefault) delete fXTLeftDefault;
    if (fXTRightDefault) delete fXTRightDefault;
    if (fXTHistEven) delete fXTHistEven;
    if (fXTHistOdd) delete fXTHistOdd;
    if (fXTHistDefault) delete fXTHistDefault;
    for (int i = 0; i<NLAY; i++){
        if (fXTLeft[i]) delete fXTLeft[i]; 
        if (fXTRight[i]) delete fXTRight[i];
        if (fXTHist[i]) delete fXTHist[i];
    }
    if (fResIntrinsic) delete fResIntrinsic;
    if (fResIntrinsicFunction) delete fResIntrinsicFunction;
}

bool XTManager::SetInputFileXT(TString file){
    fInputFileXT = new TFile(file);
    if (!fInputFileXT||fInputFileXT->IsZombie()){
        fInputFileXT = NULL;
        MyWarn("The given XT file \""<<file<<"\" is not available! will disard this setting and seek for default XT.");
        return false;
    }
    return true;
}

bool XTManager::Initialize(){
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int runNo = RunInfoManager::Get().runNo;
    TString gasTypeShort = RunInfoManager::Get().gasTypeShort;
    int HV = RunInfoManager::Get().HV;
    TString runName = RunInfoManager::Get().runName;

    xtType = ParameterManager::Get().XTManagerParameters.xtType;

    // Prepare XT functions
    if (!fInputFileXT||fInputFileXT->IsZombie()){
        if(!SetInputFileXT(HOME+Form("/info/xt.%d.%s.root",runNo,runName.Data()))){
            if(!SetInputFileXT(HOME+Form("/info/xt.%s.%d.root",gasTypeShort.Data(),HV))){
                MyError("Cannot find the default XT neither! XTManager failed to initialize.");
                return false;
            }
        }
    }
    for (int i = 1; i<NLAY; i++){ // first layer (0) is dummy
        if (useSideXT){
            fXTLeft[i] = (TF1*) fInputFileXT->Get(Form("fls_%d",i));
            fXTRight[i] = (TF1*) fInputFileXT->Get(Form("frs_%d",i));
        }
        else{
            fXTLeft[i] = (TF1*) fInputFileXT->Get(Form("fl_%d",i));
            fXTRight[i] = (TF1*) fInputFileXT->Get(Form("fr_%d",i));
        }
        fXTHist[i] = (TH2D*) fInputFileXT->Get(Form("h2_dt_%d",i));
        if (!fXTHist[i]) fXTHist[i] = (TH2D*) fInputFileXT->Get(Form("h2_xt_%d",i));
        if (fXTLeft[i]) fXTLeftDefault = fXTLeft[i];
        if (fXTRight[i]) fXTRightDefault = fXTRight[i];
        if (fXTHist[i]) fXTHistDefault = fXTHist[i];
    }
    int lid = ParameterManager::Get().XTManagerParameters.defaultLayer;
    if (lid>=NLAY||lid<0) {MyError("Invalid default layer "<<lid); return false;}
    if (fXTLeft[lid]) fXTLeftDefault = fXTLeft[lid];
    if (fXTRight[lid]) fXTRightDefault = fXTRight[lid];
    if (fXTHist[lid]) fXTHistDefault = fXTHist[lid];
    if (useSideXT){
        if (!fXTLeftDefault) fXTLeftDefault = (TF1*) fInputFileXT->Get("fls_0");
        if (!fXTRightDefault) fXTRightDefault = (TF1*) fInputFileXT->Get("frs_0");
    }
    else{
        if (!fXTLeftDefault) fXTLeftDefault = (TF1*) fInputFileXT->Get("fl_0");
        if (!fXTRightDefault) fXTRightDefault = (TF1*) fInputFileXT->Get("fr_0");
    }
    if (!fXTHistDefault) fXTHistDefault = (TH2D*) fInputFileXT->Get("h2_xt_0");
    if (!fXTLeftDefault||!fXTRightDefault){
        MyError("Cannot find the default xt functions in the given xt file!");
        return false;
    }
    if (!fXTHistDefault){
        MyWarn("Cannot find the default xt 2D histogram. XTManager will continue to serve for t2x function but will fail at RandomDriftT function!!!");
    }
    lid = ParameterManager::Get().XTManagerParameters.evenLayer;
    if (lid>=NLAY||lid<0) {MyError("Invalid even layer "<<lid); return false;}
    fXTLeftEven = fXTLeft[lid];
    fXTRightEven = fXTRight[lid];
    fXTHistEven = fXTHist[lid];
    lid = ParameterManager::Get().XTManagerParameters.oddLayer;
    if (lid>=NLAY||lid<0) {MyError("Invalid odd layer "<<lid); return false;}
    fXTLeftOdd = fXTLeft[lid];
    fXTRightOdd = fXTRight[lid];
    fXTHistOdd = fXTHist[lid];

    // Prepare error function
    fInputFileRes = new TFile(HOME+Form("/info/reso.%d.",runNo)+runName+".root");
    if (!fInputFileRes||fInputFileRes->IsZombie()){
        MyWarn("Cannot find reso file according to the given run name. Will use default reso instead.");
        fInputFileRes = new TFile(HOME+Form("/info/reso.%s.%d.root",gasTypeShort.Data(),HV));
        if (!fInputFileRes||fInputFileRes->IsZombie()){
            MyError("Cannot find the default reso: "<<HOME+Form("/info/reso.%s.%d.root",gasTypeShort.Data(),HV));
            return false;
        }
    }
    fResIntrinsic = (TGraph*)fInputFileRes->Get("gr_resInt");
    if(!fResIntrinsic) fResIntrinsic = (TGraph*)fInputFileRes->Get("gr_resIni"); // DPRECATED. For old compatebility
    fResIntrinsicFunction = (TF1*)fInputFileRes->Get("f_resInt");
    if (!fResIntrinsicFunction) fResIntrinsicFunction = (TF1*)fInputFileRes->Get("f_resIni"); // DPRECATED. For old compatebility

    return true;
}

// TODO: better to make status as return value and the to-be-calculated number (x) as an input reference. Same task goes with x2t and RandomDriftT
double XTManager::t2x(double time, int lid, int wid, double lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
    TF1* f = choseFunction(lid,wid,lr);
    if (!f){
        MyError("Cannot get XT curve for layer "<<lid<<"!\n");
        status = -2;
        return 0;
    }
    double tmax = f->GetXmax();
    double tmin = f->GetXmin();
    status = 0;
    double dd = 0;
    if (time>tmax){
        status = 1;
        dd = f->Eval(tmax);
    }
    else if (time<tmin){
        status = -1;
        dd = 0;
    }
    else {
        status = 0;
        dd = f->Eval(time);
    }
    return dd;
}

double XTManager::x2t(double doca, int lid, int wid){
    TF1* f = choseFunction(lid,wid,doca);
    if (!f){
        MyError("Cannot get XT curve for layer "<<lid<<"!\n");
        return 0;
    }
    // to avoid out of range error
    double docaMin = f->GetMinimum();
    double docaMax = f->GetMaximum();
    double dtMin = f->GetMinimumX();
    double dtMax = f->GetMaximumX();
    double dt;
    if (doca>docaMax){
        dt = dtMax;
    }
    else if (doca<docaMin){
        dt = dtMin;
    }
    else{
        dt = f->GetX(doca);
    }
    return dt;
}

double XTManager::RandomDriftT(double doca, int lid, int wid){
    TH2D* h=0;
    // FIXME: consider to add customizable parameters
    if (fXTHist[lid]) h = fXTHist[lid];
    else if (fXTHistDefault) h = fXTHistDefault;
    if (!h){
        MyError("Cannot get XT relation 2D histogram for layer "<<lid<<"!\n");
        return 0;
    }
    TH1D * hpy = (TH1D*)gFile->Get("RandomDriftT_py");
    if (!hpy) hpy = h->ProjectionY("RandomDriftT_py");
    int ibin = hpy->FindBin(doca);
    TH1D * hpx = (TH1D*)gFile->Get(Form("RandomDriftT_px%d_%d",lid,ibin));
    if (!hpx) hpx = h->ProjectionX(Form("RandomDriftT_px%d_%d",lid,ibin),ibin,ibin);
    return hpx->GetRandom();
}

double XTManager::GetError(double dd){
    dd = fabs(dd); // TODO: either consider about left/right asymmetry or add other parameter control
    double error = 0.2; // default value 200 um
    if (fResIntrinsicFunction){
        error = fResIntrinsicFunction->Eval(dd>7.5?7.5:dd);
    }
    else{
        int N = fResIntrinsic->GetN();
        for (int i = 0; i<N-1; i++){
            double d1,sig1;
            double d2,sig2;
            fResIntrinsic->GetPoint(i,d1,sig1);
            fResIntrinsic->GetPoint(i+1,d2,sig2);
            if (d2>7){
                error = sig1;
                break;
            }
            else if (d1<dd&&d2>=dd){
                error = (sig1*(d2-dd)+sig2*(dd-d1))/(d2-d1);
                break;
            }
        }
    }
    return error;
}

void XTManager::Print(){
    printf("XTManager: \n"); 
    for (int i = 1; i<NLAY; i++){ // first layer (0) is dummy
        printf("  XT in layer[%d]: ",i);
        PrintXTfunc(fXTLeft[i],fXTRight[i]);
    }
    printf("  XT in even layers: ");
    PrintXTfunc(fXTLeftEven,fXTRightEven);
    printf("  XT in odd layers: ");
    PrintXTfunc(fXTLeftOdd,fXTRightOdd);
    printf("  XT in default layers: ");
    PrintXTfunc(fXTLeftDefault,fXTRightDefault);

    // TODO print resolution
}

bool XTManager::PrintXTfunc(const TF1 * fl, const TF1 * fr){
    if (!fl||!fr){
        MyNamedLog("XTManager","Empty XT function!");
        return false;
    }
    double tmaxl = 0;
    double tmaxr = 0;
    double tminl = 0;
    double tminr = 0;
    double xmaxl = 0;
    double xmaxr = 0;
    double xminl = 0;
    double xminr = 0;
    tmaxl = fl->GetXmax();
    tminl = fl->GetXmin();
    xmaxl = fl->Eval(tmaxl);
    xminl = fl->Eval(tminl);
    tmaxr = fr->GetXmax();
    tminr = fr->GetXmin();
    xmaxr = fr->Eval(tmaxr);
    xminr = fr->Eval(tminr);
    printf("Left %s (%.1f ns, %.2f mm)~(%.1f ns, %.2f mm), Right %s (%.1f ns, %.2f mm)~(%.1f ns, %.2f mm)\n",fl->GetName(),tmaxl,xmaxl,tminl,xminl,fr->GetName(),tminr,xminr,tmaxr,xmaxr);
    return true;
}

TF1 * XTManager::choseFunction(int lid, int wid, double lr){
    TF1 * f = NULL;
    // FIXME: consider to use left/right case and folded case
    if (xtType == kSingleFolded){
        if (lr>=0) f = fXTRightDefault;
        else       f = fXTLeftDefault;
    }
    else if (xtType == kSingleLeftRight){
        if (lr>=0) f = fXTRightDefault;
        else       f = fXTLeftDefault;
    }
    else if (xtType == kEvenOddFolded){
        if (lid%2==0){
            if (lr>=0) f = fXTRightEven;
            else       f = fXTLeftEven;
        }
        else{
            if (lr>=0) f = fXTRightOdd;
            else       f = fXTLeftOdd;
        }
    }
    else if (xtType == kEvenOddLeftRight){
        if (lid%2==0){
            if (lr>=0) f = fXTRightEven;
            else       f = fXTLeftEven;
        }
        else{
            if (lr>=0) f = fXTRightOdd;
            else       f = fXTLeftOdd;
        }
    }
    else if (xtType == kAllFolded){
        if (lr>=0){
            f = fXTRight[lid];
        }
        else {
            f = fXTLeft[lid];
        }
    }
    else if (xtType == kAllLeftRight){
        if (lr>=0){
            f = fXTRight[lid];
        }
        else {
            f = fXTLeft[lid];
        }
    }
    return f;
}
