#ifndef RunInfoManager_hxx_seen
#define RunInfoManager_hxx_seen

#include <TString.h>

class RunInfoManager{
public:
    RunInfoManager();
    virtual ~RunInfoManager(){};
    
    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static RunInfoManager& Get(void) {
        if (!fRunInfoManager)
            fRunInfoManager = new RunInfoManager();
        return *fRunInfoManager;
    }

    bool Initialize(int theRunNo, TString thePreRunName, TString theRunName, int theTestLayer);
    void Print();

    TString preRunName;
    TString runName;
    int     testLayer;
    int     runNo;
    TString gasType;
    TString gasTypeShort;
    int     gasID;
    TString duration;
    double  durationTime;
    int     runGr;
    double  W;
    double  npair_per_cm;
    int     HV;
    int     THR;
    double  t00;
    double  t01;
    double  aacut;
    double  sumcut;

private:
    /// The static pointer to the singleton instance.
    static RunInfoManager* fRunInfoManager;
};

#endif
