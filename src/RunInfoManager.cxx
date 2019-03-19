#include <stdio.h>  /* printf, getenv */
#include <TString.h>
#include <TChain.h>

#include "RunInfoManager.hxx"

RunInfoManager* RunInfoManager::fRunInfoManager = NULL;

RunInfoManager::RunInfoManager():
    runName(""),
    preRunName(""),
    testLayer(0),
    runNo(0),
    gasType(""),
    gasTypeShort(""),
    gasID(0),
    W(0),
    npair_per_cm(0),
    duration(""),
    durationTime(0),
    runGr(0),
    HV(0),
    THR(0),
    t00(0),
    t01(0),
    aacut(0),
    sumcut(0)
{
}

bool RunInfoManager::Initialize(int theRunNo, TString thePreRunName, TString theRunName, int theTestLayer){
    preRunName = thePreRunName;
    runName = theRunName;
    testLayer = theTestLayer;
    TString HOME=getenv("CDCS8WORKING_DIR");;
    TString file;
    file = HOME+"/Input/run-info.root";
    // get run info
    char runDu[128];
    TChain * iChain = new TChain("t");
    iChain->Add(file);
    if (!iChain->GetEntries()){
        return false;
    }
    iChain->SetBranchAddress("run_number",&runNo);
    iChain->SetBranchAddress("gas_mixture_id",&gasID);
    iChain->SetBranchAddress("hv_ch0",&HV);
    iChain->SetBranchAddress("recbe_th_input_bd0",&THR);
    iChain->SetBranchAddress("duration",&runDu);
    iChain->SetBranchAddress("run_grade",&runGr);
    iChain->SetBranchAddress("t00",&t00);
    iChain->SetBranchAddress("t01",&t01);
    iChain->SetBranchAddress("aa",&aacut);
    iChain->SetBranchAddress("sum",&sumcut);
    for(int i = 0; i<iChain->GetEntries(); i++){
        iChain->GetEntry(i);
        if (runNo == theRunNo) break;
    }

    npair_per_cm = 60;
    gasType = "He:C_{2}H_{6}(50:50)";
    gasTypeShort = "C2H6";
    W = 32; // eV
    if (gasID==1){
        gasType = "He:iC_{4}H_{10}(90:10)";
        gasTypeShort = "C4H10";
        npair_per_cm = 29;
        W = 39; // eV
    }
    else if (gasID==2){
        gasType = "He:CH_{4}(80:20)";
        gasTypeShort = "CH4";
        npair_per_cm = 17;
        W = 39; // eV
    }
    duration = runDu;
    durationTime = 0;
    const char *sep = ":";
    char * durationSep = strtok(runDu,sep);
    double timeunit = 3600;
    while(durationSep){
        durationTime += timeunit*strtol(durationSep,NULL,10);
        timeunit/=60;
        durationSep = strtok(NULL,sep);
    }
    return true;
}

void RunInfoManager::Print(){
    printf("runNo#%d: %s, %d, %s, %d V, %d mV, %.0f sec\n",runNo,gasType.Data(),runGr,duration.Data(),HV,THR,durationTime);
    printf("Test layer %d, previous run: \"%s\", this run: \"%s\"\n",testLayer,preRunName.Data(),runName.Data());
}
