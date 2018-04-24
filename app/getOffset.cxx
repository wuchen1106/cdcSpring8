#include <vector>
#include <stdlib.h>
#include <math.h>

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

#include "XTAnalyzer.h"

#include "header.h"

void printUsage(char * name);
int getHitType(int type,bool isRight);

int main(int argc, char** argv){

	if (argc<6){
	    printUsage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
    TString prerunname = argv[2];
    TString runname = argv[3];
	int geoSetup = (int)strtol(argv[4],NULL,10); // 0: normal scintillator; 1: finger scintillator
	int wptype = (int)strtol(argv[5],NULL,10); // 0: don't generate new wiremap; 1: generate new one with offset
    double scale = 1; // normalize the step size with this scale when updating wiremap
    if (argc>6)
        scale = (double)strtod(argv[6],NULL);
    int debugLevel = 0;
    double stepSize = 0;
    if (argc>7)
        stepSize = (double)strtod(argv[7],NULL);
    double minslz = 0;
    double maxslz = 0;
    if (argc>9){
        minslz = (double)strtod(argv[8],NULL);
        maxslz = (double)strtod(argv[9],NULL);
    }
    double mininx = 0;
    double maxinx = 0;
    if (argc>11){
        mininx = (double)strtod(argv[10],NULL);
        maxinx = (double)strtod(argv[11],NULL);
    }
    double maxchi2 = 1;
    if (argc>12)
        maxchi2= (double)strtod(argv[12],NULL);
    if (argc>13)
        debugLevel = (int)strtol(argv[13],NULL,10);
    std::vector<int> wireIDs;
    if (argc>14){
        for (int i = 14; i<argc; i++){
            wireIDs.push_back((int)strtol(argv[i],NULL,10));
        }
    }
    else{
        for (int i = 0; i<NCEL; i++){
            wireIDs.push_back(i);
        }
    }
    printf("##############Input %d Parameters##################\n",argc);
    printf("runNo       = %d\n",runNo);
    printf("prerunname  = \"%s\"\n",prerunname.Data());
    printf("runname     = \"%s\"\n",runname.Data());
    printf("geoSetup:     %s\n",geoSetup==0?"normal scintillator":"finger scintillator");
    printf("wptype:       %s\n",wptype==0?"no new wiremap":"new wiremap");
    printf("scale       = %.3e\n",scale);
    printf("stepSize    = %.3e\n",stepSize);
    printf("minslz      = %.3e\n",minslz);
    printf("maxslz      = %.3e\n",maxslz);
    printf("mininx      = %.3e\n",mininx);
    printf("maxinx      = %.3e\n",maxinx);
    printf("maxchi2     = %.3e\n",maxchi2);
    printf("debug       = %d\n",debugLevel);
    printf("wires       : ");
    for (int i = 0; i<wireIDs.size(); i++){
        printf("%d ",wireIDs[i]);
    }
    printf("\n");
    fflush(stdout);

    TString HOME=getenv("CDCS8WORKING_DIR");

	// prepare offset
	TH1D * h_off[NLAY][NCEL];
	TH1D * h_slz[NLAY][NCEL];
	TH1D * h_inx[NLAY][NCEL];
	double v_off[NLAY][NCEL];
	double v_slz[NLAY][NCEL];
	double v_inx[NLAY][NCEL];
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			v_off[lid][wid] = 0;
			v_slz[lid][wid] = 0;
			v_inx[lid][wid] = 0;
			h_off[lid][wid] = new TH1D(Form("h_off_%d_%d",lid,wid),Form("Offset of wire [%d,%d]",lid,wid),128,-1,1);
			h_slz[lid][wid] = new TH1D(Form("h_slz_%d_%d",lid,wid),Form("slope Z of hits on wire [%d,%d]",lid,wid),128,-0.16,0.16);
			h_inx[lid][wid] = new TH1D(Form("h_inx_%d_%d",lid,wid),Form("intercetp X of hits on wire [%d,%d]",lid,wid),512,-50,50);
		}
	}

	// get slzmc inxmc of each wire
    double v_inxmc[NLAY][NCEL];
    double v_slzmc[NLAY][NCEL];
    for (int lid = 0; lid<NLAY; lid++){
        for (int wid = 0; wid<NCEL; wid++){
            v_inxmc[lid][wid] = 0;
            v_slzmc[lid][wid] = 0;
        }
    }
	TChain * ichain_beam = new TChain("t","t");
	ichain_beam->Add(Form("%s/Input/beammap.%d.root",HOME.Data(),runNo));
	if (ichain_beam->GetEntries()){
        double beam_inxmc;
        double beam_slzmc;
        int    beam_lid;
        int    beam_wid;
        ichain_beam->SetBranchAddress("inxmc",&beam_inxmc);
        ichain_beam->SetBranchAddress("slzmc",&beam_slzmc);
        ichain_beam->SetBranchAddress("lid",&beam_lid);
        ichain_beam->SetBranchAddress("wid",&beam_wid);
        for (int i = 0; i<ichain_beam->GetEntries(); i++){
            ichain_beam->GetEntry(i);
            if (beam_lid<0||beam_lid>=NLAY||beam_wid<0||beam_wid>=NCEL) continue;
            v_inxmc[beam_lid][beam_wid] = beam_inxmc;
            v_slzmc[beam_lid][beam_wid] = beam_slzmc;
        }
    }
    else{
        printf("Cannot find%s/Input/beammap.%d.root, will assume slzmc/inxmc center at 0\n",HOME.Data(),runNo);
    }

    // input file
    int triggerNumber;
    int nHits;
    int nHitsG;
    std::vector<int> *    i_layerID = 0;
    std::vector<int> *    i_wireID = 0;
    std::vector<double> * i_driftT = 0;
    std::vector<double> * i_driftDmc = 0;
    std::vector<int> *    i_type = 0;
    std::vector<int> *    i_np = 0;
    std::vector<int> *    i_ip = 0;
    std::vector<int> *    i_clk = 0;
    std::vector<int> *    i_width = 0;
    std::vector<int> *    i_peak = 0;
    std::vector<int> *    i_height = 0;
    std::vector<int> *    i_mpn = 0;
    std::vector<int> *    i_mpi = 0;
    std::vector<int> *    i_rank = 0;
    bool has_rank = false;
    std::vector<double> * i_aa = 0;
    std::vector<double> * i_ped = 0;
    bool has_ped = false;
    std::vector<double> * i_sum = 0;
    std::vector<double> * i_driftD[NCAND] = {0};
    std::vector<double> * i_calD[NCAND] = {0};
    int npairs[NCAND];
    int isel[NCAND];
    int icom[NCAND];
    double iinx[NCAND];
    double islx[NCAND];
    double iinz[NCAND];
    double islz[NCAND];
    double chi2x[NCAND];
    double chi2z[NCAND];
    double chi2i[NCAND];
    int nHitsS[NCAND];
    double inx[NCAND];
    double slx[NCAND];
    double inz[NCAND];
    double slz[NCAND];
    double chi2[NCAND];
    double chi2p[NCAND];
    double chi2a[NCAND];
    double chi2mc[NCAND];
    double chi2pmc[NCAND];
    double chi2amc[NCAND];
    double inxmc;
    double inzmc;
    double slxmc;
    double slzmc;
    std::vector<double> * i_fitD[NCAND] = {0};
    std::vector<int> * i_sel[NCAND] = {0};

    // Loop in layers
	for (int lid = 1; lid<NLAY; lid++){
        if (debugLevel>0) {printf("In Layer %d: preparing input TChain\n",lid);fflush(stdout);}
		TChain * ichain = new TChain("t","t");
		ichain->Add(Form("%s/root/t_%d.%s.layer%d.root",HOME.Data(),runNo,runname.Data(),lid));
        ichain->SetBranchAddress("triggerNumber",&triggerNumber);
        ichain->SetBranchAddress("nHits",&nHits);
        ichain->SetBranchAddress("nHitsG",&nHitsG);
        ichain->SetBranchAddress("layerID",&i_layerID);
        ichain->SetBranchAddress("wireID",&i_wireID);
        ichain->SetBranchAddress("driftT",&i_driftT);
        ichain->SetBranchAddress("type",&i_type);
        ichain->SetBranchAddress("np",&i_np);
        ichain->SetBranchAddress("ip",&i_ip);
        ichain->SetBranchAddress("clk",&i_clk);
        ichain->SetBranchAddress("width",&i_width);
        ichain->SetBranchAddress("peak",&i_peak);
        ichain->SetBranchAddress("height",&i_height);
        ichain->SetBranchAddress("mpn",&i_mpn);
        ichain->SetBranchAddress("mpi",&i_mpi);
        has_rank = (ichain->SetBranchAddress("rank",&i_rank)==0);
        ichain->SetBranchAddress("aa",&i_aa);
        has_ped = (ichain->SetBranchAddress("ped",&i_ped)==0);
		ichain->SetBranchAddress("sum",&i_sum);
		for (int iCand = 0; iCand<NCAND; iCand++){
			ichain->SetBranchAddress(Form("driftD%d",iCand),&(i_driftD[iCand]));
			ichain->SetBranchAddress(Form("npairs%d",iCand),&(npairs[iCand]));
			ichain->SetBranchAddress(Form("isel%d",iCand),&(isel[iCand]));
			ichain->SetBranchAddress(Form("icom%d",iCand),&(icom[iCand]));
			ichain->SetBranchAddress(Form("islx%d",iCand),&(islx[iCand]));
			ichain->SetBranchAddress(Form("islz%d",iCand),&(islz[iCand]));
			ichain->SetBranchAddress(Form("iinx%d",iCand),&(iinx[iCand]));
			ichain->SetBranchAddress(Form("iinz%d",iCand),&(iinz[iCand]));
			ichain->SetBranchAddress(Form("chi2x%d",iCand),&(chi2x[iCand]));
			ichain->SetBranchAddress(Form("chi2z%d",iCand),&(chi2z[iCand]));
			ichain->SetBranchAddress(Form("chi2i%d",iCand),&(chi2i[iCand]));
			ichain->SetBranchAddress(Form("calD%d",iCand),&(i_calD[iCand]));
			ichain->SetBranchAddress(Form("nHitsS%d",iCand),&(nHitsS[iCand]));
			ichain->SetBranchAddress(Form("slx%d",iCand),&(slx[iCand]));
			ichain->SetBranchAddress(Form("slz%d",iCand),&(slz[iCand]));
			ichain->SetBranchAddress(Form("inx%d",iCand),&(inx[iCand]));
			ichain->SetBranchAddress(Form("inz%d",iCand),&(inz[iCand]));
			ichain->SetBranchAddress(Form("chi2%d",iCand),&(chi2[iCand]));
			ichain->SetBranchAddress(Form("chi2p%d",iCand),&(chi2p[iCand]));
			ichain->SetBranchAddress(Form("chi2a%d",iCand),&(chi2a[iCand]));
			ichain->SetBranchAddress(Form("fitD%d",iCand),&(i_fitD[iCand]));
			ichain->SetBranchAddress(Form("sel%d",iCand),&(i_sel[iCand]));
		}

        // Loop in events
        ichain->GetEntries();
        Long64_t N = ichain->GetEntries();
        if (N==0){
            fprintf(stderr,"WARNING: \"%s/root/t_%d.%s.layer%d.root\" is empty! Will ignore this layer.\n",HOME.Data(),runNo,runname.Data(),lid);
            continue;
        }
        if (debugLevel>0) {printf("Processing %d events\n",N);fflush(stdout);}
        for ( int iEntry = 0 ; iEntry<N; iEntry++){
            if (iEntry%10000==0) printf("%d\n",iEntry);
            if (debugLevel>=20) printf("Entry%d: \n",iEntry);
            ichain->GetEntry(iEntry);

			// decide which candidate to use
			int theCand = 0;

            // ignore events with bad fitting
            if (nHitsS[theCand]<6) continue;
			if (chi2[theCand]>maxchi2) continue;
            //if (nHitsG>nHitsS[theCand]) continue;
            if (geoSetup==1){
                if (fabs(inz[theCand])>24) continue;
            }
            else{
                if (fabs(slz[theCand])>0.15) continue;
            }

            if (debugLevel>=20) printf("  Good Event! Looping in %d hits\n",nHits);
            // find the closest hit in the test layer
            double minres = 1e9;
            bool has = false;
            int wireID;
            double driftD, driftT, fitD;
            // FIXME: test more cut
            bool hasBadHit = false;
            for (int ihit = 0; ihit<nHits; ihit++){
                int tlayerID = (*i_layerID)[ihit];
                int twireID = (*i_wireID)[ihit];
                double tfitD = (*i_fitD[theCand])[ihit];
                double tdriftD = (*i_driftD[theCand])[ihit];
            	if ((*i_sel[theCand])[ihit]==1&&(fabs(tdriftD)<0.5||fabs(tdriftD)>7.5)) hasBadHit = true;
                if (tlayerID!=lid) continue;
				int ttype = getHitType((*i_type)[ihit],tfitD>=0);
                if (fabs(tfitD-tdriftD)<fabs(minres)){ // no cut for test layer!
                    minres = tfitD-tdriftD;
                    wireID = (*i_wireID)[ihit];
                    fitD = tfitD;
                    driftD = tdriftD;
                    driftT = (*i_driftT)[ihit];
                    has = true;
                }
            }
            if (!has) continue; // no hits found in test layer
            //if (hasBadHit) continue;

            if (debugLevel>=20) printf("  Found hit! pushing to XTAnalyzer\n");
			// tell analyzer a new data point
			if (fabs(driftD)>2&&fabs(driftD)<6&&abs(fitD-driftD)<1){ // only trust the body part
                h_off[lid][wireID]->Fill(fitD-driftD);
                h_slz[lid][wireID]->Fill(slz[theCand]);
                h_inx[lid][wireID]->Fill(inx[theCand]);
            }
        }
	}

	// output
	TFile * ofile = new TFile(Form("%s/info/offset.%d.%s.root",HOME.Data(),runNo,runname.Data()),"RECREATE");
	TTree * otree = new TTree("t","t");
	double o_off_delta;
	double o_off_slz;
	double o_off_inx;
	int o_off_lid;
	int o_off_wid;
	otree->Branch("d",&o_off_delta);
	otree->Branch("slz",&o_off_slz);
	otree->Branch("inx",&o_off_inx);
	otree->Branch("wid",&o_off_wid);
	otree->Branch("lid",&o_off_lid);
	for (int lid = 0; lid<NLAY; lid++){
		for (int wid = 0; wid<NCEL; wid++){
			h_off[lid][wid]->Write();
			h_slz[lid][wid]->Write();
			h_inx[lid][wid]->Write();
			if (h_off[lid][wid]->GetEntries()<100) continue;
			o_off_wid = wid;
			o_off_lid = lid;
			o_off_delta = h_off[lid][wid]->GetMean();
			v_off[lid][wid] = o_off_delta;
			o_off_slz = h_slz[lid][wid]->GetMean();
			v_slz[lid][wid] = o_off_slz;
			o_off_inx = h_inx[lid][wid]->GetMean();
			v_inx[lid][wid] = o_off_inx;
			otree->Fill();
		}
	}
	otree->Write();
	ofile->Close();

    if (wptype){
        //===================Get Wire Position============================
        TFile * TFile_wirepos = new TFile(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),runNo,prerunname.Data()));
        TTree * TTree_wirepos = (TTree*) TFile_wirepos->Get("t");
        int     wp_bid;
        int     wp_ch;
        int     wp_wid;
        int     wp_lid;
        double  wp_xro;
        double  wp_yro;
        double  wp_xc;
        double  wp_yc;
        double  wp_xhv;
        double  wp_yhv;
        TTree_wirepos->SetBranchAddress("b",&wp_bid);
        TTree_wirepos->SetBranchAddress("ch",&wp_ch);
        TTree_wirepos->SetBranchAddress("l",&wp_lid);
        TTree_wirepos->SetBranchAddress("w",&wp_wid);
        TTree_wirepos->SetBranchAddress("xhv",&wp_xhv);
        TTree_wirepos->SetBranchAddress("yhv",&wp_yhv);
        TTree_wirepos->SetBranchAddress("xro",&wp_xro);
        TTree_wirepos->SetBranchAddress("yro",&wp_yro);
        std::vector<int>     vwp_bid;
        std::vector<int>     vwp_ch;
        std::vector<int>     vwp_wid;
        std::vector<int>     vwp_lid;
        std::vector<double>  vwp_xro;
        std::vector<double>  vwp_yro;
        std::vector<double>  vwp_xc;
        std::vector<double>  vwp_yc;
        std::vector<double>  vwp_xhv;
        std::vector<double>  vwp_yhv;
        for (int i = 0; i<TTree_wirepos->GetEntries(); i++){
            TTree_wirepos->GetEntry(i);
            vwp_bid.push_back(wp_bid);
            vwp_ch.push_back(wp_ch);
            vwp_wid.push_back(wp_wid);
            vwp_lid.push_back(wp_lid);
            vwp_xro.push_back(wp_xro);
            vwp_yro.push_back(wp_yro);
            vwp_xc.push_back(wp_xc);
            vwp_yc.push_back(wp_yc);
            vwp_xhv.push_back(wp_xhv);
            vwp_yhv.push_back(wp_yhv);
        }
        TFile_wirepos->Close();

        //===================Set Wire Position============================
        TFile_wirepos = new TFile(Form("%s/info/wire-position.%d.%s.root",HOME.Data(),runNo,runname.Data()),"RECREATE");
        TTree_wirepos = new TTree("t","t");
        TTree_wirepos->Branch("b",&wp_bid);
        TTree_wirepos->Branch("ch",&wp_ch);
        TTree_wirepos->Branch("l",&wp_lid);
        TTree_wirepos->Branch("w",&wp_wid);
        TTree_wirepos->Branch("xhv",&wp_xhv);
        TTree_wirepos->Branch("yhv",&wp_yhv);
        TTree_wirepos->Branch("xc",&wp_xc);
        TTree_wirepos->Branch("yc",&wp_yc);
        TTree_wirepos->Branch("xro",&wp_xro);
        TTree_wirepos->Branch("yro",&wp_yro);
        for (int i = 0; i<vwp_bid.size(); i++){
            wp_bid = vwp_bid[i];
            wp_ch = vwp_ch[i];
            wp_wid = vwp_wid[i];
            wp_lid = vwp_lid[i];
            wp_xro = vwp_xro[i];
            wp_yro = vwp_yro[i];
            wp_xc = vwp_xc[i];
            wp_yc = vwp_yc[i];
            wp_xhv = vwp_xhv[i];
            wp_yhv = vwp_yhv[i];
            bool doOffset = false;
            for (int index = 0; index < wireIDs.size(); index++){
                if (wp_wid==wireIDs[index]){
                    doOffset = true;
                    break;
                }
            }
            if (doOffset&&wp_lid>0&&wp_lid<NLAY&&wp_wid>0&&wp_wid<NCEL){
                double theOff = v_off[wp_lid][wp_wid]*scale;
                double deltaSlz = v_slz[wp_lid][wp_wid]-v_slzmc[wp_lid][wp_wid];
                double deltaInx = v_inx[wp_lid][wp_wid]-v_inxmc[wp_lid][wp_wid];
                if (stepSize){
                    if(fabs(theOff)>stepSize) theOff = theOff>0?stepSize:-stepSize;
                }
                if ((minslz||maxslz)&&(deltaSlz>maxslz||deltaSlz<minslz)){
                    theOff = 0;
                }
                if ((mininx||maxinx)&&(deltaInx>maxinx||deltaInx<mininx)){
                    theOff = 0;
                }
                wp_xro += theOff;
                wp_xc += theOff;
                wp_xhv += theOff;
            }
            TTree_wirepos->Fill();
        }
        TTree_wirepos->Write();
        TFile_wirepos->Close();
    }

    return 0;
}

int getHitType(int type,bool isRight){
	int ttype = (type/10)%10;
	if (isRight){
		if (ttype==1||ttype==4) type-=ttype*10; // l- or l+
	}
	else{
		if (ttype==2||ttype==5) type-=ttype*10; // r- or r+
	}
	return type;
}

void printUsage(char * name){
    fprintf(stderr,"%s [runNo] [prerunname] [runname] [geoSetup: 0, normal;1, finger] [wptype: 0, no update; 1, update] <[scale (1)] [stepSize (0)] [minslz (0)] [maxslz (0)] [mininx (0)] [maxinx (0)] [maxchi2 (1)] [debug: (0)] [wireIDs (all)]>\n",name);
}
