#include <stdio.h>  /* printf, getenv */
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>

#include "Log.hxx"
#include "XTManager.hxx"
#include "RunInfoManager.hxx"

XTManager* XTManager::fXTManager = NULL;

XTManager::XTManager():
    fInputFileXT(0),
    fInputFileRes(0),
	fXTLeftEven(0),
	fXTRightEven(0),
	fXTLeftOdd(0),
	fXTRightOdd(0),
	fXTLeftDefault(0),
	fXTRightDefault(0),
    fResIntrinsic(0)
{
    xtType = kSingleFolded;
    for (int i = 0; i<NLAY; i++){
        fXTLeft[i] = 0; 
        fXTRight[i] = 0;
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
    for (int i = 0; i<NLAY; i++){
        if (fXTLeft[i]) delete fXTLeft[i]; 
        if (fXTRight[i]) delete fXTRight[i];
    }
    if (fResIntrinsic) delete fResIntrinsic;
}

bool XTManager::Initialize(){
    TString HOME=getenv("CDCS8WORKING_DIR");;
    int runNo = RunInfoManager::Get().runNo;
    TString gasTypeShort = RunInfoManager::Get().gasTypeShort;
    int HV = RunInfoManager::Get().HV;
    TString runName = RunInfoManager::Get().preRunName;

    // Prepare XT functions
    fInputFileXT = new TFile(HOME+Form("/info/xt.%d.",runNo)+runName+".root");
    if (!fInputFileXT||fInputFileXT->IsZombie()){
        MyWarn("Cannot find xt file according to the given run name. Will use garfield xt instead.");
        fInputFileXT = new TFile(HOME+Form("/info/xt.%s.%d.root",gasTypeShort.Data(),HV));
        if (!fInputFileXT||fInputFileXT->IsZombie()){
            MyError("Cannot find the default garfield xt: "<<HOME+Form("/info/xt.%s.%d.root",gasTypeShort.Data(),HV));
            return false;
        }
    }
    for (int i = 1; i<NLAY; i++){ // first layer (0) is dummy
        fXTLeft[i] = (TF1*) fInputFileXT->Get(Form("fl_%d",i));
        fXTRight[i] = (TF1*) fInputFileXT->Get(Form("fr_%d",i));
    }
	fXTLeftEven = (TF1*) fInputFileXT->Get("fl_even");
	fXTRightEven = (TF1*) fInputFileXT->Get("fr_even");
	fXTLeftOdd = (TF1*) fInputFileXT->Get("fl_odd");
	fXTRightOdd = (TF1*) fInputFileXT->Get("fr_odd");
	fXTLeftDefault = (TF1*) fInputFileXT->Get("fl_0");
	fXTRightDefault = (TF1*) fInputFileXT->Get("fr_0");

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
	fResIntrinsic = (TGraph*)fInputFileRes->Get("gr_resIni");

    return true;
}

double XTManager::t2x(double time, int lid, int wid, int lr, int & status){ // 1: right; 2: right end; -1: left; -2: left end; 0 out of range
    TF1* f=0;
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

double XTManager::GetError(double dd){
	double error = 0.2; // default value 200 um
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
    if (fl||fr){
        printf("Empty!\n");
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
    printf("Left (%.1f ns, %.2f mm)~(%.1f ns, %.2f mm), Right (%.1f ns, %.2f mm)~(%.1f ns, %.2f mm)\n",tmaxl,xmaxl,tminl,xminl,tminr,xminr,tmaxr,xmaxr);
    return true;
}
