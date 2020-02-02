#include "TText.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TF1.h"

#include <math.h>
#include <stdlib.h>

#include "MyProcessManager.hxx"
#include "Log.hxx"

#include "ParameterManager.hxx"
#include "RunInfoManager.hxx"
#include "BeamManager.hxx"
#include "GeometryManager.hxx"
#include "XTManager.hxx"
#include "InputOutputManager.hxx"

//#define PRINT_CROSSPOINTS

#define MULTI 5 // maximum peak multiplicity of one cell to be considered in drawing

#define NLAY4DR 8
#define NWIRE4DR 8
#define NCELL4DR (NLAY4DR*NWIRE4DR)
#define IFWIRE4DR 1
#define IFLAY4DR 1

// for test 180404
//#define NLAY4DR 2
//#define NWIRE4DR 2
//#define NCELL4DR (NLAY4DR*NWIRE4DR)
//#define IFWIRE4DR 1
//#define IFLAY4DR 5

//===================About xt============================
double tdc2t(int tdc);
void cid4dr2lidwid(int cid, int & lid, int & wid);
bool lidwid2cid4dr(int lid, int wid, int & cid);

bool isGood(int iHit);
void getRunTimeParameters(TString configureFile);
void print_usage(char* prog_name);

MyProcessManager * pMyProcessManager;

int main(int argc, char** argv){
    TString HOME=getenv("CDCS8WORKING_DIR");
    int m_runNo = 0;
    TString m_runName = "";
    int m_iEntryStart = -1;
    int m_iEntryStop = -1;
    int m_nEntries = 0;
    int m_modulo = 100;
    bool m_memdebug = false;
    int m_testLayer = 4;
    double m_slxini = 0;
    TString m_suffixHitFile = "";
    TString m_wireAdjustmentFile = "";
    TString m_xtFile = "";
    bool    m_drawZX = false;

    // Load options
    int    opt_result;
    while((opt_result=getopt(argc,argv,"A:B:C:D:E:H:L:MN:P:R:V:X:Z"))!=-1){
        switch(opt_result){
            case 'M':
                m_memdebug = true;
                printf("Turning on memory debug\n");
                break;
            case 'P':
                m_modulo = atoi(optarg);
                printf("Printing modulo set to %d\n",m_modulo);
                break;
            case 'R':
                m_runNo = atoi(optarg);
                break;
            case 'L':
                m_testLayer = atoi(optarg);
                break;
            case 'B':
                m_iEntryStart = atoi(optarg);
                break;
            case 'E':
                m_iEntryStop = atoi(optarg);
                break;
            case 'N':
                m_nEntries = atoi(optarg);
                break;
            case 'C':
                getRunTimeParameters(optarg);
                printf("Using configure file \"%s\"\n",optarg);
                break;
            case 'H':
                m_suffixHitFile = optarg;
                printf("Added suffix \"%s\" to the output file\n",optarg);
                break;
            case 'A':
                m_wireAdjustmentFile = optarg;
                break;
            case 'X':
                m_xtFile = optarg;
                printf("Load xt curves from %s\n",m_xtFile.Data());
                break;
            case 'Z':
                m_drawZX = true;
                printf("Will draw ZX planes\n");
                break;
            case 'D':
                if (!Log::ConfigureD(optarg)) print_usage(argv[0]);
                break;
            case 'V':
                if (!Log::ConfigureV(optarg)) print_usage(argv[0]);
                break;
            case '?':
                printf("Wrong option! optopt=%c, optarg=%s\n", optopt, optarg);
            case 'h':
            default:
                print_usage(argv[0]);
                return 1;
        }
    }
    if (m_nEntries>0){
        if (m_iEntryStart<0) m_iEntryStart = 0;
        m_iEntryStop = m_iEntryStart+m_nEntries-1;
    }

    if (argc-optind<1){
        print_usage(argv[0]);
        return -1;
    }
    m_runName= argv[optind++];
    m_slxini = -18.4*M_PI/180;

    printf("##############%s##################\n",argv[0]);
    printf("runNo               = %d\n",m_runNo);
    printf("input XT File       = %s\n",m_xtFile==""?"self":m_xtFile.Data());
    printf("runName             = \"%s\"\n",m_runName.Data());
    printf("Test layer ID       = %d\n",m_testLayer);
    printf("Start Entry         = %d\n",m_iEntryStart);
    printf("Stop Entry          = %d\n",m_iEntryStop);
    printf("Using wire adjustment file \"%s\"\n",m_wireAdjustmentFile.Data());

    if (m_memdebug){
        pMyProcessManager = MyProcessManager::GetMyProcessManager();
        if (!pMyProcessManager) return -1;
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());

    // Prepare managers
    bool success = false;
    success = RunInfoManager::Get().Initialize(m_runNo,m_runName,m_testLayer);RunInfoManager::Get().Print();
    if (!success) {MyError("Cannot initialize RunInfoManager"); return 1;}
    success = BeamManager::Get().Initialize(ParameterManager::Get().beamType);BeamManager::Get().Print();
    if (!success) {MyError("Cannot initialize BeamManager"); return 1;}
    success = GeometryManager::Get().Initialize(ParameterManager::Get().geoSetup,ParameterManager::Get().connectionType,ParameterManager::Get().chamberType); GeometryManager::Get().Print();
    if (!success) {MyError("Cannot initialize GeometryManager"); return 1;}
    if (m_wireAdjustmentFile!=""){
        success = GeometryManager::Get().AdjustWirePosition(m_wireAdjustmentFile);
        if (!success) MyWarn("Cannot load offset file for wire adjustment. Will ignore this step.");
    }
    if (m_xtFile!=""){
        success = XTManager::Get().SetInputFileXT(m_xtFile);
        if (!success){MyError("Invalid input XT file"); return 1;}
    }
    success = XTManager::Get().Initialize();XTManager::Get().Print();
    if (!success) {MyError("Cannot initialize XTManager"); return 1;}
    XTManager::Get().Print();
    InputOutputManager::Get().readRawFile = true;
    InputOutputManager::Get().readPeakFile = true;
    InputOutputManager::Get().readHitFile = true;
    InputOutputManager::Get().readTrackFile = true;
    InputOutputManager::Get().SetHitFileSuffix(m_suffixHitFile); // the output file name will be like h_100SUFFIX.root
    success = InputOutputManager::Get().Initialize(true); // with trivial branches
    if (!success) {MyError("Cannot initialize InputOutputManager"); return 1;}
    bool m_foundRawFile = InputOutputManager::Get().IsRawFileReady();
    bool m_foundPeakFile = InputOutputManager::Get().IsPeakFileReady();
    bool m_foundHitFile = InputOutputManager::Get().IsHitFileReady();
    bool m_foundTrackFile = InputOutputManager::Get().IsTrackFileReady();
    printf("Available input files:\n");
    printf("        raw file?    %s\n",m_foundRawFile?"yes":"no");
    printf("        hit file?   %s\n",m_foundHitFile?"yes":"no");
    printf("        peak file?  %s\n",m_foundPeakFile?"yes":"no");
    printf("        track file? %s\n",m_foundTrackFile?"yes":"no");

    double XMAX = 130; // range of x-z plane
    double ZMAX = 350; // range of x-z plane

    //==================Prepare canvas for drawing==========================
    // run summary
    printf("Preparing canvas...\n");
    TLatex * text_runsum = new TLatex();
    text_runsum->SetTextSize(0.02);
    Long64_t nEntries = InputOutputManager::Get().GetEntries();
    Long64_t triggerNumberMax = InputOutputManager::Get().GetTriggerNumberMax();
    text_runsum->SetText(0.05,0.98,Form("run#%d ",m_runNo)+RunInfoManager::Get().gasType+Form(", %d V,%d mV, Grade#%d",RunInfoManager::Get().HV,RunInfoManager::Get().THR,RunInfoManager::Get().runGr)+", "+RunInfoManager::Get().duration+Form(", %lld events, Eff_{daq} = %2.2lf%%, Rate_{tri} = %1.1lfkHz",nEntries,((double)nEntries)/(triggerNumberMax+1)*100,(triggerNumberMax+1)/RunInfoManager::Get().durationTime/1000));

    TLatex * text_title = new TLatex(0,0,"");
    text_title->SetTextSize(0.02);
    //Prepare the Canvas for waveforms (by board)
    TCanvas * ca_WF[NBRD];
    TPad * pad_WF[NBRD][NCHS];
    for (int bid = 0; bid<NBRD; bid++){
        ca_WF[bid] = new TCanvas(Form("ca_WF_%d",bid),"ca_WF",896,896);
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        for (int i = 0; i<8; i++){
            for (int j = 0; j<6; j++){
                int index = j*8+i;
                pad_WF[bid][index] = new TPad(Form("pad_%d_%d_%d",bid,i,j),Form("pad_%d_%d_%d",bid,i,j),1./8*i,0.95/6*(5-j),1./8*(i+1),0.95/6*(6-j));
                pad_WF[bid][index]->Draw();
                pad_WF[bid][index]->SetGridx(1);
                pad_WF[bid][index]->SetGridy(1);
            }
        }
    }
    //Prepare the Canvas for waveforms (by layer)
    TCanvas * ca_WFL;
    TPad * pad_WFL[NCELL4DR];
    ca_WFL = new TCanvas("ca_WFL","ca_WFL",1024,1024.*NLAY4DR/NWIRE4DR);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    for (int i = 0; i<NWIRE4DR; i++){
        for (int j = 0; j<NLAY4DR; j++){
            int index = j*NWIRE4DR+i;
            pad_WFL[index] = new TPad(Form("padL_%d_%d",i,j),Form("padL_%d_%d",i,j),1./NWIRE4DR*i,0.95/NLAY4DR*(NLAY4DR-1-j),1./NWIRE4DR*(i+1),0.95/NLAY4DR*(NLAY4DR-j));
            pad_WFL[index]->Draw();
            //pad_WFL[index]->SetGridx(1);
            //pad_WFL[index]->SetGridy(1);
        }
    }
    //Prepare the Canvas for x-y plane and target chanel ADC
    TCanvas * ca_xyADC = new TCanvas("ca_xyADC","ca_xyADC",896,1024);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    TPad * pad_xyADC = new TPad("pad_xyADC","pad_xyADC",0,0,1,0.96);
    pad_xyADC->Draw();
    pad_xyADC->SetGridx(1);
    pad_xyADC->SetGridy(1);
    //Prepare the Canvas for z-x planes
    TCanvas* ca_zx_all = new TCanvas("ca_zx_all","ca_zx_all",1024,768);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    TCanvas* ca_zx[NZXP]; // z-x planes corresponding to the layerID of the lower layer counting from 1
    if (m_drawZX){
        for (int izx = 1; izx<NZXP; izx++){
            ca_zx[izx] = new TCanvas(Form("ca_zx_%d",izx),"ca_zx",1024,768);
            gStyle->SetPalette(1);
            gStyle->SetOptStat(0);
            gStyle->SetPadTickX(1);
            gStyle->SetPadTickY(1);
            gPad->SetGridx(1);
            gPad->SetGridy(1);
        }
    }
    // Prepare colors to be used for each wire
    int color[NCEL];
    int icolor = 0;
    color[icolor++] = kBlack;
    color[icolor++] = kMagenta+2;
    color[icolor++] = kMagenta;
    color[icolor++] = kBlue;
    color[icolor++] = kBlue-2;
    color[icolor++] = kCyan;
    color[icolor++] = kGreen;
    color[icolor++] = kYellow;
    color[icolor++] = kOrange;
    color[icolor++] = kRed;
    color[icolor++] = kRed+2;

    //==================Prepare drawing objects==========================
    //Prepare for ADC
    TGraph * gr_waveForm[NBRD][NCHS] = {0};
    int vSample[NSAM];
    for (int i=0; i<NSAM; i++){
        vSample[i] = i;
    }
    TLatex *textWF[NBRD][NCHS];
    TLatex *textTDC[NBRD][NCHS][NSAM];
    TMarker *markerTDC[NBRD][NCHS][NSAM];
    for (int bid = 0; bid<NBRD; bid++){
        for (int ch = 0; ch<NCHS; ch++){
            textWF[bid][ch] = new TLatex(0,0,"");
            textWF[bid][ch]->SetTextSize(0.04);
            for (int sid=0; sid<NSAM; sid++) {
                textTDC[bid][ch][sid] = new TLatex(0,0,"");
                textTDC[bid][ch][sid]->SetTextSize(0.04);
                markerTDC[bid][ch][sid] = new TMarker(0,0,20);
                markerTDC[bid][ch][sid]->SetMarkerSize(0.3);
            }
        } 
    }
    //Prepare for hit circles on x-y plane
    TLatex * text_xyhit[NLAY][NCEL];
    TEllipse * circle_driftD[NLAY][NCEL][MULTI];  // value from track finding/fitting;
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            text_xyhit[lid][wid] = new TLatex(0,0,"");
            text_xyhit[lid][wid]->SetTextSize(0.02);
            for (int ip = 0; ip<MULTI; ip++){
                circle_driftD[lid][wid][ip] = new TEllipse(0,0,1,1);
                circle_driftD[lid][wid][ip]->SetFillStyle(0);
            }
        }
    }
    //Prepare two tracks on x-y plane: l_track: track after fitting; l_itrack: track before fitting, after track finding.
    double yup = GeometryManager::Get().fChamber->chamberPositionY+GeometryManager::Get().fChamber->chamberHeight/2;
    double ydown = GeometryManager::Get().fChamber->chamberPositionY-GeometryManager::Get().fChamber->chamberHeight/2;
    TLine * l_track = new TLine();
    l_track->SetLineColor(kRed);
    l_track->SetY1(yup-5);
    l_track->SetY2(ydown+5);
    TLine * l_itrack = new TLine();
    l_itrack->SetLineColor(kBlue);
    l_itrack->SetY1(yup);
    l_itrack->SetY2(ydown);
    // Background graph for x-y plane
    TGraph * gr_wireCenter = new TGraph();
    gr_wireCenter->Set(NLAY*NCEL);
    int count_nwires = 0;
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            if (GeometryManager::Get().fChamber->wire_ch[lid][wid]!=-1){
                double xc = (GeometryManager::Get().fChamber->wire_x[lid][wid][0]+GeometryManager::Get().fChamber->wire_x[lid][wid][1])/2;
                double yc = (GeometryManager::Get().fChamber->wire_y[lid][wid][0]+GeometryManager::Get().fChamber->wire_y[lid][wid][1])/2;
                gr_wireCenter->SetPoint(count_nwires,xc,yc);
                count_nwires++;
            }
        }
    }
    gr_wireCenter->Set(count_nwires);
    gr_wireCenter->SetMarkerStyle(20);
    gr_wireCenter->SetMarkerSize(0.55);
    gr_wireCenter->GetXaxis()->SetTitle("x [mm]");
    gr_wireCenter->GetYaxis()->SetTitle("y [mm]");
    // Prepare driftT lines on z-x planes
    TLine * l_zx[NLAY][NCEL][MULTI][2];
    if (m_drawZX&&m_foundTrackFile){
        for (int lid = 0; lid<NLAY; lid++){
            for (int wid = 0; wid<NCEL; wid++){
                for (int ip = 0; ip<MULTI; ip++){
                    for (int ilr = 0; ilr<2; ilr++){
                        l_zx[lid][wid][ip][ilr] = new TLine();
                        l_zx[lid][wid][ip][ilr]->SetLineColor(color[wid]);
                        l_zx[lid][wid][ip][ilr]->SetLineStyle(3);
                        l_zx[lid][wid][ip][ilr]->SetX1(-GeometryManager::Get().fChamber->chamberLength/2);
                        l_zx[lid][wid][ip][ilr]->SetX2(GeometryManager::Get().fChamber->chamberLength/2);
                    }
                }
            }
        }
    }
    // Prepare cross points of driftT lines on z-x planes
    TMarker * point_cross_zx[NZXP][NCEL][NCEL][MULTI][MULTI][4];
    TText * text_cross_zx[NZXP][NCEL][NCEL][MULTI][MULTI][4];
    if (m_drawZX&&m_foundTrackFile){
        for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
            for (int wid = 0; wid<NCEL; wid++){
                for (int wjd = 0; wjd<NCEL; wjd++){
                    for (int ip = 0; ip<MULTI; ip++){
                        for (int jp = 0; jp<MULTI; jp++){
                            for (int icombi = 0; icombi<4; icombi++){
                                if (!m_foundTrackFile||icombi==3) // reverse izx when icombi is 0 or 1, reverse izx+1 when icombi is 0 or 2. Keep them unchanged only when icombi is 3
                                    point_cross_zx[izx][wid][wjd][ip][jp][icombi] = new TMarker(0,0,20);
                                else
                                    point_cross_zx[izx][wid][wjd][ip][jp][icombi] = new TMarker(0,0,4);
                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(kBlack);
                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerSize(0.55);
                                text_cross_zx[izx][wid][wjd][ip][jp][icombi] = new TText(0,0,Form("%d",izx));
                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextColor(color[izx]);
                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextSize(0.01);
                            }
                        }
                    }
                }
            }
        }
    }
    // Prepare track points on z-x planes
    TMarker * point_track_zx[NZXP];
    TMarker * point_itrack_zx[NZXP];
    if (m_drawZX&&m_foundTrackFile){
        for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
            point_track_zx[izx] = new TMarker(0,0,20);
            point_track_zx[izx]->SetMarkerColor(kRed);
            point_track_zx[izx]->SetMarkerSize(0.5);
            point_itrack_zx[izx] = new TMarker(0,0,4);
            point_itrack_zx[izx]->SetMarkerColor(kBlue);
            point_itrack_zx[izx]->SetMarkerSize(0.6);
        }
    }
    // Prepare wires in each layer on z-x planes
    TGraph * gr_wire[NLAY][NCEL];
    if (m_foundTrackFile){
        for (int lid = 1; lid<NLAY; lid++){
            for (int wid = 0; wid<NCEL; wid++){
                if (GeometryManager::Get().fChamber->wire_ch[lid][wid]!=-1){
                    gr_wire[lid][wid] = new TGraph(2);
                    gr_wire[lid][wid]->SetPoint(0,GeometryManager::Get().fChamber->wire_z[lid][wid][0],GeometryManager::Get().fChamber->wire_x[lid][wid][0]);
                    gr_wire[lid][wid]->SetPoint(1,GeometryManager::Get().fChamber->wire_z[lid][wid][1],GeometryManager::Get().fChamber->wire_x[lid][wid][1]);
                    gr_wire[lid][wid]->SetLineColor(color[wid]);
                    gr_wire[lid][wid]->SetMarkerColor(color[wid]);
                }
            }
        }
    }
    // Prepare texts and markers for wire on each z-x planes
    TMarker * point_cross_wire[NZXP][NCEL][NCEL];
    TLatex * text_cross_wire[NZXP][NCEL][NCEL];
    TLatex * text_cr_1[NZXP];
    TLatex * text_cr_1l[NZXP][NCEL];
    TLatex * text_cr_2[NZXP];
    TLatex * text_cr_2r[NZXP][NCEL];
    TGraph * gr_all[NZXP];
    double ar_cr_x[2] = {-ZMAX,ZMAX};
    double ar_cr_y[2] = {-XMAX,XMAX};
    gr_all[0] = new TGraph(2,ar_cr_x,ar_cr_y);
    gr_all[0]->SetTitle(Form("Cross points of all layers"));
    gr_all[0]->SetMarkerColor(kWhite);
    gr_all[0]->SetMarkerSize(0.1);
    gr_all[0]->GetXaxis()->SetRangeUser(-ZMAX,ZMAX);
    gr_all[0]->GetYaxis()->SetRangeUser(-XMAX,XMAX);
    gr_all[0]->GetXaxis()->SetTitle("z [mm]");
    gr_all[0]->GetYaxis()->SetTitle("x [mm]");
    if (m_drawZX&&m_foundTrackFile){
        for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
            // Background graph for z-x planes
            gr_all[izx] = new TGraph(2,ar_cr_x,ar_cr_y);
            gr_all[izx]->SetTitle(Form("layer #%d and layer #%d",izx,izx+1));
            gr_all[izx]->SetMarkerColor(kWhite);
            gr_all[izx]->SetMarkerSize(0.1);
            gr_all[izx]->GetXaxis()->SetRangeUser(-ZMAX,ZMAX);
            gr_all[izx]->GetYaxis()->SetRangeUser(-XMAX,XMAX);
            gr_all[izx]->GetXaxis()->SetTitle("z [mm]");
            gr_all[izx]->GetYaxis()->SetTitle("x [mm]");
            text_cr_1[izx] = new TLatex(GeometryManager::Get().fChamber->wire_z[izx][0][0]-30,XMAX-20,Form("Layer %d",izx));
            text_cr_1[izx]->SetTextSize(0.025);
            text_cr_2[izx] = new TLatex(GeometryManager::Get().fChamber->wire_z[izx+1][0][1]-30,XMAX-20,Form("Layer %d",izx+1));
            text_cr_2[izx]->SetTextSize(0.025);
            for (int wid = 0; wid<NCEL; wid++){
                if (GeometryManager::Get().fChamber->wire_ch[izx][wid]!=-1){
                    text_cr_1l[izx][wid] = new TLatex(GeometryManager::Get().fChamber->wire_z[izx][wid][0]-10,GeometryManager::Get().fChamber->wire_x[izx][wid][0],Form("%d",wid));
                    text_cr_1l[izx][wid]->SetTextColor(color[wid]);
                    text_cr_1l[izx][wid]->SetTextSize(0.02);
                }
                if (GeometryManager::Get().fChamber->wire_ch[izx+1][wid]!=-1){
                    text_cr_2r[izx][wid] = new TLatex(GeometryManager::Get().fChamber->wire_z[izx+1][wid][1]+10,GeometryManager::Get().fChamber->wire_x[izx+1][wid][1],Form("%d",wid));
                    text_cr_2r[izx][wid]->SetTextColor(color[wid]);
                    text_cr_2r[izx][wid]->SetTextSize(0.02);
                }
#ifdef PRINT_CROSSPOINTS
                for (int wjd = 0; wjd < NCEL; wjd++){
                    if (fabs(GeometryManager::Get().fChamber->wirecross_z[izx][wid][wjd])>300) continue;
                    point_cross_wire[izx][wid][wjd] = new TMarker(GeometryManager::Get().fChamber->wirecross_z[izx][wid][wjd],GeometryManager::Get().fChamber->wirecross_x[izx][wid][wjd],20);
                    point_cross_wire[izx][wid][wjd]->SetMarkerColor(color[wid]);
                    point_cross_wire[izx][wid][wjd]->SetMarkerSize(0.4);
                    text_cross_wire[izx][wid][wjd] = new TLatex(GeometryManager::Get().fChamber->wirecross_z[izx][wid][wjd],GeometryManager::Get().fChamber->wirecross_x[izx][wid][wjd]+5,Form("%d,%d",wid,wjd));
                    text_cross_wire[izx][wid][wjd]->SetTextColor(color[wid]);
                    text_cross_wire[izx][wid][wjd]->SetTextSize(0.02);
                }
#endif
            }
        }
    }

    //===================Prepare counters============================
    //  and to record the position of hits in each layer
    int    nHits_layer[NLAY];
    int    nHits_cell[NLAY][NCEL];
    int    iHit_cell[NLAY][NCEL][MULTI];
    std::map<int,bool> iHit_selected;
    // to record each dd & fd
    double dd_cell[NLAY][NCEL][MULTI];
    double fd_cell[NLAY][NCEL][MULTI];
    // to record track position on each layer
    double y_cell[NLAY][NCEL];
    for (int lid = 1; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            if (GeometryManager::Get().fChamber->wire_ch[lid][wid]!=-1){
                y_cell[lid][wid] = (GeometryManager::Get().fChamber->wire_y[lid][wid][1]+GeometryManager::Get().fChamber->wire_y[lid][wid][0])/2.;
            }
        }
    }

    //===================Loop in Events============================
    TString prefix = "";
    MyNamedLog("display","  Start loop");
    Long64_t N = InputOutputManager::Get().GetEntries();
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    if (m_iEntryStop<0||m_iEntryStart<0){m_iEntryStart = 0; m_iEntryStop=N-1;}
    for (Long64_t iEntry = m_iEntryStart; iEntry<=m_iEntryStop; iEntry++){
        MyNamedTrace("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
        if (iEntry%m_modulo == 0){
            MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
            std::cout<<iEntry<<std::endl;
        }
        InputOutputManager::Get().Reset();
        InputOutputManager::Get().GetEntry(iEntry);

        int    nHits = InputOutputManager::Get().nHits;
        int    nHitsG = InputOutputManager::Get().nHitsG;
        if (nHits<=0){
            MyNamedLog("display","No hits in event "<<iEntry<<", continue");
            continue;
        }

        // set prefix
        if (m_foundTrackFile){
            if (nHitsG<=6) prefix = "incom.";
            else if (nHitsG==7 ) prefix = "single.";
            else if (nHitsG==8 ) prefix = "single.";
            else if (nHitsG==9 ) prefix = "single.";
            else if (nHitsG==10) prefix = "n10.";
            else if (nHitsG==11) prefix = "n11.";
            else if (nHitsG==12) prefix = "n12.";
            else if (nHitsG>=13) prefix = "multi.";
        }

        // Reset counters
        for (int lid = 0; lid<NLAY; lid++){
            nHits_layer[lid] = 0;
            for (int wid = 0; wid<NCEL; wid++){
                nHits_cell[lid][wid] = 0;
            }
        }

        // count
        for (int ihit = 0; ihit<nHits; ihit++){
            int lid = InputOutputManager::Get().LayerID->at(ihit);
            int wid = InputOutputManager::Get().CellID->at(ihit);
            if (nHits_cell[lid][wid]<MULTI){ // don't count the rest
                iHit_cell[lid][wid][nHits_cell[lid][wid]] = ihit;
                nHits_cell[lid][wid]++;
            }
            nHits_layer[lid]++;
        }

        //===================Draw waveforms============================
        // clear objects from the previous event
        for (int bid = 0; bid<NBRD; bid++){
            for (int ch = 0; ch<NCHS; ch++){
                if (gr_waveForm[bid][ch]) delete gr_waveForm[bid][ch]; gr_waveForm[bid][ch] = 0;
            }
        }
        // Draw waveforms
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->cd();
            text_title->SetText(0.05,0.96,Form("Entry %lld, Trigger Number %d, %d TDCs, %d considered in tracking",iEntry,InputOutputManager::Get().triggerNumber,nHits,nHitsG));
            text_title->Draw();
            text_runsum->Draw();
            for (int ch = 0; ch<NCHS; ch++){
                // get waveform and draw
                int chg = bid*NCHS+ch;
                int lid = GeometryManager::Get().fChamber->wire_lid[bid][ch];
                int wid = GeometryManager::Get().fChamber->wire_wid[bid][ch];
                pad_WF[bid][ch]->cd();
                gr_waveForm[bid][ch] = new TGraph(NSAM,vSample,InputOutputManager::Get().adc[chg]);
                gr_waveForm[bid][ch]->SetTitle(Form("Layer %d Wire %d Board %d Channel %d ",lid,wid,bid,ch));
                gr_waveForm[bid][ch]->GetXaxis()->SetRangeUser(0,NSAM-1);
                gr_waveForm[bid][ch]->GetYaxis()->SetRangeUser(MIN_ADC,MAX_ADC);
                gr_waveForm[bid][ch]->GetXaxis()->SetTitle("Sample Index");
                gr_waveForm[bid][ch]->GetYaxis()->SetTitle("ADC");
                gr_waveForm[bid][ch]->SetMarkerStyle(20);
                gr_waveForm[bid][ch]->SetMarkerSize(0.3);
                gr_waveForm[bid][ch]->Draw("APL");
                // set title text
                //textWF[bid][ch]->SetText(1,MAX_ADC-50,Form("%d peaks, ADC sum %.0lf",InputOutputManager::Get().tdcNhit[chg],pk_aa[chg])); // TODO: how to get adc sum all for the given channel
                textWF[bid][ch]->SetText(1,MAX_ADC-50,Form("%d peaks",InputOutputManager::Get().tdcNhit[chg]));
                if (InputOutputManager::Get().tdcNhit[chg]>0)
                    textWF[bid][ch]->SetTextColor(kBlue); // has hit
                else
                    textWF[bid][ch]->SetTextColor(kGreen); // no hit
                textWF[bid][ch]->Draw();
            }
        }
        for (int ihit = 0; ihit<nHits; ihit++){ // update with hit information
            int lid = InputOutputManager::Get().LayerID->at(ihit);
            int wid = InputOutputManager::Get().CellID->at(ihit);
            int bid = GeometryManager::Get().fChamber->wire_bid[lid][wid];
            int ch = GeometryManager::Get().fChamber->wire_ch[lid][wid];
            int chg = bid*NCHS+ch;
            int ip = InputOutputManager::Get().iPeakInChannel->at(ihit);
            int clk = InputOutputManager::Get().TDCClock->at(ihit);
            int height = InputOutputManager::Get().ADCheight->at(ihit);
            ca_WF[bid]->cd();
            pad_WF[bid][ch]->cd();
            markerTDC[bid][ch][ip]->SetX(clk);
            markerTDC[bid][ch][ip]->SetY(height);
            textTDC[bid][ch][ip]->SetText(clk,height,Form("%d",(int)(InputOutputManager::Get().tdc[chg][ip])));
            // TODO: consider about this iPeakInChannel
            //int iPeak = InputOutputManager::Get().iPeakInChannel->at(iHit);
            // if (peakType!= TrackingPara::kAllPeaks&&iPeakOverCut!=0) good = false;
            if (isGood(ihit)){ // good hit
                markerTDC[bid][ch][ip]->SetMarkerColor(kRed);
                textTDC[bid][ch][ip]->SetTextColor(kRed);
                textWF[bid][ch]->SetTextColor(kRed);
                textWF[bid][ch]->Draw();
            }
            else{ // bad hit
                markerTDC[bid][ch][ip]->SetMarkerColor(14);
                textTDC[bid][ch][ip]->SetTextColor(14);
            }
            markerTDC[bid][ch][ip]->Draw();
            textTDC[bid][ch][ip]->Draw();
        }
        for (int bid = 0; bid<NBRD; bid++){
            ca_WF[bid]->SaveAs(prefix+Form("wf.%lld.b%d.pdf",iEntry,bid));
            ca_WF[bid]->SaveAs(prefix+Form("wf.%lld.b%d.png",iEntry,bid));
        }

        //===================Draw waveforms by layer============================
        ca_WFL->cd();
        text_title->Draw();
        text_runsum->Draw();
        for (int cell = 0; cell<NCELL4DR; cell++){
            // get waveform and draw
            int lid = 0;
            int wid = 0;
            cid4dr2lidwid(cell,lid,wid);
            int ch = GeometryManager::Get().fChamber->wire_ch[lid][wid];
            int bid = GeometryManager::Get().fChamber->wire_bid[lid][wid];
            pad_WFL[cell]->cd();
            gr_waveForm[bid][ch]->GetYaxis()->UnZoom();
            gr_waveForm[bid][ch]->Draw("APL");
            double maxadc  = gr_waveForm[bid][ch]->GetYaxis()->GetXmax();
            double minadc  = gr_waveForm[bid][ch]->GetYaxis()->GetXmin();
            double rangeadc = maxadc-minadc;
            textWF[bid][ch]->SetY(minadc+rangeadc/20);
            textWF[bid][ch]->Draw();
        }
        for (int ihit = 0; ihit<nHits; ihit++){ // update with hit information
            int lid = InputOutputManager::Get().LayerID->at(ihit);
            int wid = InputOutputManager::Get().CellID->at(ihit);
            int cell = 0;
            bool notincluded = lidwid2cid4dr(lid,wid,cell);
            if (notincluded) continue;
            if (cell<0||cell>=NCELL4DR) continue;
            int bid = GeometryManager::Get().fChamber->wire_bid[lid][wid];
            int ch = GeometryManager::Get().fChamber->wire_ch[lid][wid];
            int chg = bid*NCHS+ch;
            int ip = InputOutputManager::Get().iPeakInChannel->at(ihit);
            pad_WFL[cell]->cd();
            if (isGood(ihit)){
                markerTDC[bid][ch][ip]->SetMarkerColor(kRed);
                textTDC[bid][ch][ip]->SetTextColor(kRed);
            }
            else{
                markerTDC[bid][ch][ip]->SetMarkerColor(kBlue);
                textTDC[bid][ch][ip]->SetTextColor(kBlue);
            }
            // for test 180404
            //markerTDC[bid][ch][ip]->SetMarkerSize(1);
            markerTDC[bid][ch][ip]->Draw();
            textTDC[bid][ch][ip]->Draw();
        }
        ca_WFL->SaveAs(prefix+Form("wfl.%lld.pdf",iEntry));
        ca_WFL->SaveAs(prefix+Form("wfl.%lld.png",iEntry));

        //===================Draw xyADC plot============================
        ca_xyADC->cd();
        text_runsum->Draw();
        for (int iCand = 0; iCand<(m_foundTrackFile?InputOutputManager::Get().nCandidatesFound:0); iCand++){
            iHit_selected.clear();
            for (int lid = 0; lid<NLAY; lid++){
                int ihit = InputOutputManager::Get().hitIndexSelected[lid][iCand];
                if (ihit>=0){
                    iHit_selected[ihit] = true;
                }
            }
            //===================Draw XY in xyADC plot============================
            // update wire position
            int count_nwires = 0;
            for (int lid = 0; lid<NLAY; lid++){
                for (int wid = 0; wid<NCEL; wid++){
                    if (GeometryManager::Get().fChamber->wire_ch[lid][wid]!=-1){
                        double y = (GeometryManager::Get().fChamber->wire_y[lid][wid][0]+GeometryManager::Get().fChamber->wire_y[lid][wid][1])/2;
                        double x;
                        if (m_foundTrackFile){ // according to z-y relation from tracking
                            double z = InputOutputManager::Get().slopeZ[iCand]*(y-GeometryManager::Get().ReferenceY)+InputOutputManager::Get().interceptZ[iCand];
                            y = GeometryManager::Get().fChamber->GetY(lid,wid,z);
                            z = InputOutputManager::Get().slopeZ[iCand]*(y-GeometryManager::Get().ReferenceY)+InputOutputManager::Get().interceptZ[iCand];
                            x = GeometryManager::Get().fChamber->GetX(lid,wid,z);
                            y = GeometryManager::Get().fChamber->GetY(lid,wid,z);
                        }
                        else{
                            x = GeometryManager::Get().fChamber->GetX(lid,wid,0);
                            y = GeometryManager::Get().fChamber->GetY(lid,wid,0);
                        }
                        gr_wireCenter->SetPoint(count_nwires,x,y);
                        count_nwires++;
                    }
                }
            }
            // Draw the background graph for x-y plane
            pad_xyADC->cd();
            if (m_foundTrackFile)
                gr_wireCenter->SetTitle(Form("Ent %lld, nHitsG (%d)%d(%d), icom %d, isel %d, sl_{z}: %.2e->%.2e, #chi^{2}: %.2e->%.2e",iEntry,nHits,nHitsG,InputOutputManager::Get().nHitsS[iCand],InputOutputManager::Get().iCombination[iCand],InputOutputManager::Get().iSelection[iCand],InputOutputManager::Get().slopeZInput[iCand],InputOutputManager::Get().slopeZ[iCand],InputOutputManager::Get().chi2Input[iCand],InputOutputManager::Get().chi2[iCand]));
            else
                gr_wireCenter->SetTitle(Form("Entry %lld nHits = %d",iEntry,nHits));
            gr_wireCenter->Draw("AP");
            // Draw the hit circles
            for (int lid = 0; lid<NLAY; lid++){
                for (int wid = 0; wid<NCEL; wid++){
                    // get wire position
                    double wxro = GeometryManager::Get().fChamber->wire_x[lid][wid][1];
                    double wyro = GeometryManager::Get().fChamber->wire_y[lid][wid][1];
                    double wzro = GeometryManager::Get().fChamber->chamberLength/2;
                    double wxhv = GeometryManager::Get().fChamber->wire_x[lid][wid][0];
                    double wyhv = GeometryManager::Get().fChamber->wire_y[lid][wid][0];
                    double wzhv = -GeometryManager::Get().fChamber->chamberLength/2;
                    double wy = (wyro+wyhv)/2.;
                    double wx = (wxro+wxhv)/2.;
                    if (m_foundTrackFile){
                        // correct wx wy wz according to the track position
                        double wz = InputOutputManager::Get().slopeZ[iCand]*(wy-GeometryManager::Get().ReferenceY)+InputOutputManager::Get().interceptZ[iCand];
                        wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
                        wz = InputOutputManager::Get().slopeZ[iCand]*(wy-GeometryManager::Get().ReferenceY)+InputOutputManager::Get().interceptZ[iCand];
                        wx = ((wzro-wz)*wxhv+(wz-wzhv)*wxro)/(wzro-wzhv);
                        wy = ((wzro-wz)*wyhv+(wz-wzhv)*wyro)/(wzro-wzhv);
                    }
                    double resmin = 1e9;
                    double thefitd = 0;
                    for (int ip = 0; ip<nHits_cell[lid][wid]; ip++){
                        int ihit = iHit_cell[lid][wid][ip];
                        // Get hit information
                        double fitd = 0;
                        if (m_foundTrackFile){
                            fitd = GeometryManager::Get().GetDOCA(InputOutputManager::Get().LayerID->at(ihit),InputOutputManager::Get().CellID->at(ihit),InputOutputManager::Get().slopeX[iCand],InputOutputManager::Get().interceptX[iCand],InputOutputManager::Get().slopeZ[iCand],InputOutputManager::Get().interceptZ[iCand]);
                        }
                        int status;
                        // TODO: add t0offset
                        double dd = XTManager::Get().t2x(InputOutputManager::Get().DriftT->at(ihit),InputOutputManager::Get().LayerID->at(ihit),InputOutputManager::Get().CellID->at(ihit),fitd,status);
                        // set the hit circles
                        MyNamedInfo("display","Hit "<<ihit<<" ["<<lid<<"]["<<wid<<"]: dd = "<<dd<<" @ "<<wx<<" "<<wy);
                        circle_driftD[lid][wid][ip]->SetX1(wx);
                        circle_driftD[lid][wid][ip]->SetY1(wy);
                        circle_driftD[lid][wid][ip]->SetR1(fabs(dd));
                        circle_driftD[lid][wid][ip]->SetR2(fabs(dd));
                        if (!isGood(ihit)){ // bad hit
                            circle_driftD[lid][wid][ip]->SetLineStyle(2);
                            circle_driftD[lid][wid][ip]->SetLineColor(14);
                        }
                        else{
                            circle_driftD[lid][wid][ip]->SetLineStyle(1);
                            if (m_foundTrackFile&&iHit_selected[ihit]) // used for fitting
                                circle_driftD[lid][wid][ip]->SetLineColor(kRed);
                            else
                                circle_driftD[lid][wid][ip]->SetLineColor(kOrange);
                        }
                        // set the text of each hit
                        circle_driftD[lid][wid][ip]->Draw(); // from track fitting/finding

                        // get the min residual
                        if (fabs(fitd-dd)<fabs(resmin)&&(m_foundTrackFile&&isGood(ihit))){ // only print when t_XXX and good hit
                            resmin = fitd-dd;
                            thefitd = fitd;
                        }

                        // Get information for cross points
                        dd_cell[lid][wid][ip] = dd;
                        fd_cell[lid][wid][ip] = fitd;
                        double delta = dd*GeometryManager::Get().fChamber->chamberLength/2*2/sqrt(GeometryManager::Get().fChamber->chamberLength/2*GeometryManager::Get().fChamber->chamberLength/2*4+(GeometryManager::Get().fChamber->wire_x[lid][wid][0]-GeometryManager::Get().fChamber->wire_x[lid][wid][1])*(GeometryManager::Get().fChamber->wire_x[lid][wid][0]-GeometryManager::Get().fChamber->wire_x[lid][wid][1])); // correction for x position
                        if (m_drawZX&&m_foundTrackFile){
                            for (int ilr = 0; ilr<2; ilr++){
                                l_zx[lid][wid][ip][ilr]->SetY1(wxhv+(ilr?delta:-delta));
                                l_zx[lid][wid][ip][ilr]->SetY2(wxro+(ilr?delta:-delta));
                            }
                        }
                    }
                    if (resmin<1e9){ // found a hit in this wire
                        text_xyhit[lid][wid]->SetText(wx,wy,Form("%.3f",resmin));// draw the one with the smallest res
                        text_xyhit[lid][wid]->Draw();
                    }
                    y_cell[lid][wid] = wy;
                }
            }
            // draw the tracks on the x-y plane
            double xdown  = InputOutputManager::Get().interceptX[iCand] + (ydown-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeX[iCand];
            double xup    = InputOutputManager::Get().interceptX[iCand] + (yup-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeX[iCand];
            double xdowni = InputOutputManager::Get().interceptXInput[iCand] + (ydown-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeXInput[iCand];
            double xupi   = InputOutputManager::Get().interceptXInput[iCand] + (yup-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeXInput[iCand];
            MyNamedInfo("display","Track slx "<<InputOutputManager::Get().slopeX[iCand]<<" inx "<<InputOutputManager::Get().interceptX[iCand]<<" slz "<<InputOutputManager::Get().slopeZ[iCand]<<" inz "<<InputOutputManager::Get().interceptZ[iCand]);
            MyNamedInfo("display","Input slx "<<InputOutputManager::Get().slopeXInput[iCand]<<" inx "<<InputOutputManager::Get().interceptXInput[iCand]<<" slz "<<InputOutputManager::Get().slopeZInput[iCand]<<" inz "<<InputOutputManager::Get().interceptZInput[iCand]);
            if (m_foundTrackFile){
                l_itrack->SetX1(xupi);
                l_itrack->SetX2(xdowni);
                l_itrack->Draw();
                l_track->SetX1(xup);
                l_track->SetX2(xdown);
                l_track->Draw();
            }
            if (m_foundTrackFile){
                ca_xyADC->SaveAs(prefix+Form("xyADC.%lld.i%d.pdf",iEntry,iCand));
                ca_xyADC->SaveAs(prefix+Form("xyADC.%lld.i%d.png",iEntry,iCand));
            }
            else{
                ca_xyADC->SaveAs(prefix+Form("xyADC.%lld.pdf",iEntry));
                ca_xyADC->SaveAs(prefix+Form("xyADC.%lld.png",iEntry));
            }

            //===================Draw ZX plots============================
            // draw the z-x planes
            if (m_drawZX&&m_foundTrackFile){
                for (int izx = 1; izx<NZXP; izx++){ // z-x planes corresponding to the layerID of the lower layer counting from 1
                    ca_zx[izx]->cd();
                    // draw Background graph for z-x planes
                    gr_all[izx]->Draw("AP");
                    text_cr_1[izx]->Draw();
                    text_cr_2[izx]->Draw();
                    for (int wid = 0; wid<NCEL; wid++){
                        if (gr_wire[izx][wid]){
                            gr_wire[izx][wid]->Draw("PLSAME");
                            text_cr_1l[izx][wid]->Draw();
                        }
                        if (gr_wire[izx+1][wid]){
                            gr_wire[izx+1][wid]->Draw("PLSAME");
                            text_cr_2r[izx][wid]->Draw();
                        }
#ifdef PRINT_CROSSPOINTS
                        for (int wjd = 0; wjd < NCEL; wjd++){
                            if (!point_cross_wire[izx][wid][wjd]
                                    ||!text_cross_wire[izx][wid][wjd]) continue;
                            point_cross_wire[izx][wid][wjd]->Draw();
                            text_cross_wire[izx][wid][wjd]->Draw();
                        }
#endif
                    }
                    // draw the track and driftT lines on the z-x planes
                    for (int ilr = 0; ilr<2; ilr++){
                        for (int wid = 0; wid<NCEL; wid++){
                            for (int ip = 0; ip<nHits_cell[izx][wid]; ip++){
                                int ihit = iHit_cell[izx][wid][ip];
                                if (isGood(ihit))
                                    l_zx[izx][wid][ip][ilr]->SetLineColor(color[wid]);
                                else 
                                    l_zx[izx][wid][ip][ilr]->SetLineColor(14);
                                l_zx[izx][wid][ip][ilr]->Draw();
                            }
                            for (int ip = 0; ip<nHits_cell[izx+1][wid]; ip++){
                                int ihit = iHit_cell[izx+1][wid][ip];
                                if (isGood(ihit))
                                    l_zx[izx+1][wid][ip][ilr]->SetLineColor(color[wid]);
                                else 
                                    l_zx[izx+1][wid][ip][ilr]->SetLineColor(14);
                                l_zx[izx+1][wid][ip][ilr]->Draw();
                            }
                        }
                    }
                    // position of the track point
                    double z_track = 0;
                    double x_track = 0;
                    double z_itrack = 0;
                    double x_itrack = 0;
                    // is there a cross point?
                    if (nHits_layer[izx]>0&&nHits_layer[izx+1]>0){
                        for (int wid = 0; wid<NCEL; wid++){
                            for (int ip = 0; ip<nHits_cell[izx][wid]; ip++){
                                int ihit = iHit_cell[izx][wid][ip];
                                for (int wjd = 0; wjd<NCEL; wjd++){
                                    for (int jp = 0; jp<nHits_cell[izx+1][wjd]; jp++){
                                        int jhit = iHit_cell[izx+1][wjd][jp];
                                        int isel = 0;
                                        int jsel = 0;
                                        if (m_foundTrackFile){
                                            isel = iHit_selected[ihit];
                                            jsel = iHit_selected[jhit];
                                        }
                                        // position of the cross point
                                        for (int icombi = 0; icombi<4; icombi++){
                                            if (!(isGood(ihit))||(!isGood(jhit))){ // bad cross
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(14);
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerSize(0.35);
                                            }
                                            else if (isel&&jsel){ // selected cross
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(kBlack);
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerSize(0.55);
                                            }
                                            else{ // good cross
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerColor(kRed);
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetMarkerSize(0.55);
                                            }
                                            double dd1 = dd_cell[izx][wid][ip];
                                            double dd2 = dd_cell[izx+1][wjd][jp];
                                            if (icombi<2) dd1 = -dd1; // reverse izx when icombi is 0 or 1
                                            if (icombi%2==0) dd2 = -dd2; // reverse izx+1 when icombi is 0 or 2
                                            double theta1 = GeometryManager::Get().fChamber->wire_theta[izx][wid];
                                            double theta2 = GeometryManager::Get().fChamber->wire_theta[izx+1][wjd];
                                            double sintheta12 = sin(theta1-theta2);
                                            double zc_fix_slx = 0;
                                            double deltaY = y_cell[izx+1][wjd]-y_cell[izx][wid];
                                            if (m_foundTrackFile){
                                                zc_fix_slx = deltaY*InputOutputManager::Get().slopeX[iCand]/(tan(theta2)-tan(theta1));
                                            }
                                            else {
                                                zc_fix_slx = deltaY*m_slxini/(tan(theta2)-tan(theta1));
                                            }
                                            double xc = GeometryManager::Get().fChamber->wirecross_x[izx][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
                                            double zc = GeometryManager::Get().fChamber->wirecross_z[izx][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
                                            if (zc>-GeometryManager::Get().fChamber->chamberLength/2&&zc<GeometryManager::Get().fChamber->chamberLength/2){
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetX(zc);
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetY(xc);
                                                point_cross_zx[izx][wid][wjd][ip][jp][icombi]->Draw();
                                            }
                                            if (m_foundTrackFile){
                                                if (icombi==3){
                                                    double y_track = (y_cell[izx][wid]+y_cell[izx+1][wjd])/2.;
                                                    z_track = InputOutputManager::Get().interceptZ[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeZ[iCand];
                                                    x_track = InputOutputManager::Get().interceptX[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeX[iCand];
                                                    z_itrack = InputOutputManager::Get().interceptZInput[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeZInput[iCand];
                                                    x_itrack = InputOutputManager::Get().interceptXInput[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeXInput[iCand];
                                                    double fd1 = fd_cell[izx][wid][ip];
                                                    double fd2 = fd_cell[izx+1][wjd][jp];
                                                    double xcf = GeometryManager::Get().fChamber->wirecross_x[izx][wid][wjd]+fd1*sin(theta2)/(-sintheta12)+fd2*sin(theta1)/sintheta12;
                                                    double zcf = GeometryManager::Get().fChamber->wirecross_z[izx][wid][wjd]+fd1*cos(theta2)/(-sintheta12)+fd2*cos(theta1)/sintheta12+zc_fix_slx;
                                                    gr_all[izx]->SetTitle(Form("DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) DD_{[%d,%d]}: %.2lf mm (%.2lf #mum) #Delta_{x}: %.0lf #mum #Delta_{z}: %.0lf #mum",izx,wid,dd1,(fd1-dd1)*1000,izx+1,wjd,dd2,(fd2-dd2)*1000,(xc-x_track)*1000,(zc-z_track)*1000));
                                                }
                                            }
                                            else{
                                                gr_all[izx]->SetTitle(Form("Layer %d and Layer %d",izx,izx+1));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else{
                        int wid = 0;
                        for (; wid<NCEL; wid++){
                            if (nHits_cell[izx][wid]) break;
                        }
                        if (wid==NCEL) wid = NCEL/2;
                        int wjd = 0;
                        for (; wjd<NCEL; wjd++){
                            if (nHits_cell[izx+1][wjd]) break;
                        }
                        if (wjd==NCEL) wjd = NCEL/2;
                        double y_track = (y_cell[izx][wid]+y_cell[izx+1][wjd])/2.; // take the y value from a previous event
                        if (m_foundTrackFile){
                            z_track = InputOutputManager::Get().interceptZ[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeZ[iCand];
                            x_track = InputOutputManager::Get().interceptX[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeX[iCand];
                            z_itrack = InputOutputManager::Get().interceptZInput[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeZInput[iCand];
                            x_itrack = InputOutputManager::Get().interceptXInput[iCand]+(y_track-GeometryManager::Get().ReferenceY)*InputOutputManager::Get().slopeXInput[iCand];
                        }
                        gr_all[izx]->SetTitle(Form("Layer %d and Layer %d",izx,izx+1));
                    }
                    // draw the track point
                    if (m_foundTrackFile){
                        point_itrack_zx[izx]->SetX(z_itrack);
                        point_itrack_zx[izx]->SetY(x_itrack);
                        point_itrack_zx[izx]->Draw();
                        point_track_zx[izx]->SetX(z_track);
                        point_track_zx[izx]->SetY(x_track);
                        point_track_zx[izx]->Draw();
                    }
                    if (m_foundTrackFile){
                        ca_zx[izx]->SaveAs(prefix+Form("zx.%lld.i%d.l%d.pdf",iEntry,iCand,izx));
                        ca_zx[izx]->SaveAs(prefix+Form("zx.%lld.i%d.l%d.png",iEntry,iCand,izx));
                    }
                    else{
                        ca_zx[izx]->SaveAs(prefix+Form("zx.%lld.l%d.pdf",iEntry,izx));
                        ca_zx[izx]->SaveAs(prefix+Form("zx.%lld.l%d.png",iEntry,izx));
                    }
                }

                // create one more z-x plane with all points on it
                ca_zx_all->cd();
                gr_all[0]->Draw("AP");
                for (int izx = 1; izx<NZXP; izx++){
                    if (m_foundTrackFile)
                        point_track_zx[izx]->Draw();
                    if (nHits_layer[izx]>0&&nHits_layer[izx+1]>0){
                        for (int wid = 0; wid<NCEL; wid++){
                            for (int ip = 0; ip<nHits_cell[izx][wid]; ip++){
                                int ihit = iHit_cell[izx][wid][ip];
                                for (int wjd = 0; wjd<NCEL; wjd++){
                                    for (int jp = 0; jp<nHits_cell[izx+1][wjd]; jp++){
                                        int jhit = iHit_cell[izx+1][wjd][jp];
                                        // position of the cross point
                                        for (int icombi = 0; icombi<4; icombi++){
                                            if ((!isGood(ihit))||(!isGood(jhit))){ // bad cross
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextColor(14);
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextSize(0.008);
                                            }
                                            else{
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextColor(color[izx]);
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetTextSize(0.01);
                                            }
                                            double dd1 = dd_cell[izx][wid][ip];
                                            double dd2 = dd_cell[izx+1][wjd][jp];
                                            if (icombi<2) dd1 = -dd1; // reverse izx when icombi is 0 or 1
                                            if (icombi%2==0) dd2 = -dd2; // reverse izx+1 when icombi is 0 or 2
                                            double theta1 = GeometryManager::Get().fChamber->wire_theta[izx][wid];
                                            double theta2 = GeometryManager::Get().fChamber->wire_theta[izx+1][wjd];
                                            double sintheta12 = sin(theta1-theta2);
                                            double zc_fix_slx = 0;
                                            double deltaY = y_cell[izx+1][wjd]-y_cell[izx][wid];
                                            if (m_foundTrackFile){
                                                zc_fix_slx = deltaY*InputOutputManager::Get().slopeX[iCand]/(tan(theta2)-tan(theta1));
                                            }
                                            else{
                                                zc_fix_slx = deltaY*m_slxini/(tan(theta2)-tan(theta1));
                                            }
                                            double xc = GeometryManager::Get().fChamber->wirecross_x[izx][wid][wjd]+dd1*sin(theta2)/(-sintheta12)+dd2*sin(theta1)/sintheta12;
                                            double zc = GeometryManager::Get().fChamber->wirecross_z[izx][wid][wjd]+dd1*cos(theta2)/(-sintheta12)+dd2*cos(theta1)/sintheta12+zc_fix_slx;
                                            if (zc>-GeometryManager::Get().fChamber->chamberLength/2&&zc<GeometryManager::Get().fChamber->chamberLength/2){
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetX(zc);
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->SetY(xc);
                                                text_cross_zx[izx][wid][wjd][ip][jp][icombi]->Draw();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (m_foundTrackFile){
                    ca_zx_all->SaveAs(prefix+Form("zx.%lld.i%d.all.png",iEntry,iCand));
                    ca_zx_all->SaveAs(prefix+Form("zx.%lld.i%d.all.pdf",iEntry,iCand));
                }
                else{
                    ca_zx_all->SaveAs(prefix+Form("zx.%lld.all.png",iEntry));
                    ca_zx_all->SaveAs(prefix+Form("zx.%lld.all.pdf",iEntry));
                }
            }
        }

        MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    }
    MyNamedDebug("Memory","Memory size: @"<<__LINE__<<": "<<pMyProcessManager->GetMemorySize());
    return 0;
}

double tdc2t(int deltaTDC){
    return (deltaTDC)/0.96;
}

void cid4dr2lidwid(int cid, int & lid, int & wid){
    lid = cid/NWIRE4DR+IFLAY4DR; // skip the first layer (dummy) by counting from 1.
    wid = cid%NWIRE4DR+IFWIRE4DR; // skip first wire (on boundary)  by counting from 1.

    // for test 180404
    //	lid = 5;
    //	if (cid==0) wid=6;
    //	else if (cid==1) wid=7;
    //	else if (cid==2) wid=1;
    //	else if (cid==3) wid=2;
}

bool lidwid2cid4dr(int lid, int wid, int & cid){
    cid = (lid-IFLAY4DR)*NWIRE4DR+wid-IFWIRE4DR;
    if (lid<IFLAY4DR||lid-IFLAY4DR>=NLAY4DR||wid<IFWIRE4DR||wid-IFWIRE4DR>=NWIRE4DR) return true;
    else return false;

    // for test 180404
    //	if (lid!=5) return true;
    //	if (wid==6) cid=0;
    //	else if (wid==7) cid=1;
    //	else if (wid==1) cid=2;
    //	else if (wid==2) cid=3;
    //	else return true;
    //	return false;
}

bool isGood(int iHit){
    double sumCut = ParameterManager::Get().TrackingParameters.sumCut;
    double aaCut = ParameterManager::Get().TrackingParameters.aaCut;
    double tmin = ParameterManager::Get().TrackingParameters.tmin;
    double tmax = ParameterManager::Get().TrackingParameters.tmax;
    double aa = InputOutputManager::Get().ADCsumAll->at(iHit);
    double sum = InputOutputManager::Get().ADCsumPacket->at(iHit);
    double driftT = InputOutputManager::Get().DriftT->at(iHit);
    int ip = InputOutputManager::Get().iPeakInChannel->at(iHit);
    bool good = true;
    if (ip>0) good = false;
    if (aa<aaCut) good = false;
    if (sum<sumCut) good = false;
    if (driftT<tmin||driftT>tmax) good = false;
    return good;
}

void print_usage(char* prog_name){
    fprintf(stderr,"Usage %s [options] runName\n",prog_name);
    fprintf(stderr,"[options]\n");
    fprintf(stderr,"\t -D <name>=[error,severe,warn,debug,trace]\n");
    fprintf(stderr,"\t\t Change the named debug level\n");
    fprintf(stderr,"\t\t if the given name is \"general\" then set the default debug level\n");
    fprintf(stderr,"\t -V <name>=[quiet,log,info,verbose]\n");
    fprintf(stderr,"\t\t Change the named log level\n");
    fprintf(stderr,"\t\t if the given name is \"general\" then set the default log level\n");
    fprintf(stderr,"\t -M\n");
    fprintf(stderr,"\t\t Turning on memory debug mode\n");
    fprintf(stderr,"\t -C <file>\n");
    fprintf(stderr,"\t\t Set the configure file\n");
    fprintf(stderr,"\t -P <n>\n");
    fprintf(stderr,"\t\t Printing modulo set to n\n");
    fprintf(stderr,"\t -R <run>\n");
    fprintf(stderr,"\t\t Run number set to run\n");
    fprintf(stderr,"\t -L <lid>\n");
    fprintf(stderr,"\t\t Test layer set to lid\n");
    fprintf(stderr,"\t -B <n>\n");
    fprintf(stderr,"\t\t Starting entry index set to n\n");
    fprintf(stderr,"\t -E <n>\n");
    fprintf(stderr,"\t\t Stopping entry index set to n\n");
    fprintf(stderr,"\t -X <xtfile>\n");
    fprintf(stderr,"\t\t Instead of searching for the xt file from this run, use the provided one\n");
    fprintf(stderr,"\t -H <suf>\n");
    fprintf(stderr,"\t\t Add suffix to the hit file, like h_100SUFFIX.root\n");
}

void getRunTimeParameters(TString configureFile){
    if (configureFile!=""){
        ParameterManager::Get().ReadInputFile(configureFile,"",false,false);
        ParameterManager::Get().LoadParameters(ParameterManager::kTracking);
        ParameterManager::Get().LoadParameters(ParameterManager::kXTManager);
    }
}
