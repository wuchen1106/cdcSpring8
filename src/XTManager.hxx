#ifndef XTManager_hxx_seen
#define XTManager_hxx_seen

#include <TString.h>

#include "GeometryManager.hxx"

class TFile;
class TF1;
class TH2D;
class TGraph;

class XTManager{
public:
    enum XTType{
        kSingleFolded,
        kSingleLeftRight,
        kAllFolded,
        kAllLeftRight,
        kEvenOddFolded,
        kEvenOddLeftRight
    };

    XTManager();
    virtual ~XTManager();
    
    ///  Get a reference to the singleton instance of dummy 
    ///  database information.  If this is first attempt at reference
    ///  then singleton is instantiated and parameters are 
    ///  read from text files.
    static XTManager& Get(void) {
        if (!fXTManager)
            fXTManager = new XTManager();
        return *fXTManager;
    }

    bool Initialize();
    double t2x(double time, int lid, int wid, double lr, int & status);
    double x2t(double doca, int lid, int wid);
    double RandomDriftT(double doca, int lid, int wid);
    double GetError(double dd);
    bool SetInputFileXT(TString file);
    void Print();
    bool PrintXTfunc(const TF1 * fl, const TF1 * fr);
    TH2D * GetXTHistDefault(){return fXTHistDefault;};
    void UseSideXT(){useSideXT = true;};

    int xtType;

private:
    TF1 * choseFunction(int lid, int wid, double lr);

    /// The static pointer to the singleton instance.
    static XTManager* fXTManager;

    bool     useSideXT;

    TFile  * fInputFileXT;
    TFile  * fInputFileRes;
    TF1    * fXTLeft[NLAY]; // TODO: support multiple XT relations in the future
    TF1    * fXTRight[NLAY];
    TF1    * fXTLeftEven;
    TF1    * fXTRightEven;
    TF1    * fXTLeftOdd;
    TF1    * fXTRightOdd;
    TF1    * fXTLeftDefault;
    TF1    * fXTRightDefault;
    TGraph * fResIntrinsic; // TODO: support multiple resolution graphs in the future
    TF1    * fResIntrinsicFunction;
    TH2D   * fXTHist[NLAY];
    TH2D   * fXTHistOdd;
    TH2D   * fXTHistEven;
    TH2D   * fXTHistDefault;
};

#endif
