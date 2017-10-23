#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"

#define NSAM 32
#define NCHS 48
#define NCHT 96
#define NBD 2
#define MIN_ADC 50
#define MAX_ADC 750
#define NBINS  256
#define MAX_WIDTH 10 // about merging. Depending on drift velocity: HV and Gas

int power2_15 = pow(2,15);

void print_usage(char* prog_name);
int main(int argc, char** argv){
	if (argc<2){
		print_usage(argv[0]);
		return 1;
	}
	int runNo = (int)strtol(argv[1],NULL,10);
	int nEventMax = 0;
	if (argc>=3) nEventMax = (int)strtol(argv[2],NULL,10);

    TString HOME=getenv("CDCS8WORKING_DIR");

	//===================Get input ROOT file============================
	TChain * c = new TChain("tree","tree");
	c->Add(HOME+Form("/root/run_%0.6d_built.root",runNo));
	int tdcNhit[NCHT];
	int clockNumberDriftTime[NCHT][NSAM];
	int adc[NCHT][NSAM];
	int driftTime[NCHT][NSAM];
	int triggerNumber;
	int triggerNumberMax = 0;
	c->SetBranchAddress("triggerNumber",&triggerNumber);
	c->SetBranchAddress("tdcNhit",tdcNhit);
	c->SetBranchAddress("clockNumberDriftTime",clockNumberDriftTime);
	c->SetBranchAddress("driftTime",driftTime);
	c->SetBranchAddress("adc",adc);

	//===================Prepare output ROOT file============================
	gROOT->ProcessLine(".L loader.C+");
	int o_nHits;
	double o_pedestal[NCHT];
	int o_nPeaks[NCHT];
	double o_areaall[NCHT];
	std::vector<std::vector<int> > * o_clk = 0;
	std::vector<std::vector<int> > * o_tdc = 0;
	std::vector<std::vector<int> > * o_peak = 0;
	std::vector<std::vector<double> > * o_sum = 0;
	std::vector<std::vector<int> > * o_width = 0;
	std::vector<std::vector<int> > * o_height = 0;
	std::vector<std::vector<int> > * o_mpn = 0;
	std::vector<std::vector<int> > * o_mpi = 0;
	TFile * f = new TFile(HOME+Form("/root/p_%d.root",runNo),"RECREATE");
	TTree * t = new TTree("t","t");
//	t->Branch("adc",adc,Form("adc[%d][%d]/I",NCHT,NSAM));
	t->Branch("triggerNumber",&triggerNumber);
	t->Branch("nh",&o_nHits,"nh/I");
	t->Branch("ped",o_pedestal,Form("ped[%d]/D",NCHT));
	t->Branch("np",o_nPeaks,Form("np[%d]/I",NCHT));
	t->Branch("aa",o_areaall,Form("aa[%d]/D",NCHT));
	t->Branch("clk",&o_clk);
	t->Branch("tdc",&o_tdc);
	t->Branch("peak",&o_peak);
	t->Branch("sum",&o_sum);
	t->Branch("width",&o_width);
	t->Branch("height",&o_height);
	t->Branch("mpn",&o_mpn);
	t->Branch("mpi",&o_mpi);
	double prepedestal[NCHT];
	for (int ch = 0; ch<NCHT; ch++){
		prepedestal[ch] = 210;
	}
    o_clk = new std::vector<std::vector<int> >;
    o_tdc = new std::vector<std::vector<int> >;
    o_peak = new std::vector<std::vector<int> >;
    o_sum = new std::vector<std::vector<double> >;
    o_width = new std::vector<std::vector<int> >;
    o_height = new std::vector<std::vector<int> >;
    o_mpn = new std::vector<std::vector<int> >;
    o_mpi = new std::vector<std::vector<int> >;
    std::vector<int> temp_clk;
    std::vector<int> temp_tdc;
    std::vector<double> temp_sum;
    std::vector<int> temp_peak;
    std::vector<int> temp_width;
    std::vector<int> temp_height;
    std::vector<int> temp_mpn;
    std::vector<int> temp_mpi;

	// Loop in events
	Long64_t N = c->GetEntries();
	if (nEventMax&&nEventMax<N) N = nEventMax;
	printf("Processing %d entries...\n",N);
	for (Long64_t i = 0;i<N; i++){
		//FIXME
		if (i%1000==0) printf("Entry %d, %.2e%%\n",i,(double)i/N*100);
		c->GetEntry(i);

        // reset
        o_nHits = 0;
        o_clk->clear();
        o_tdc->clear();
        o_peak->clear();
        o_sum->clear();
        o_width->clear();
        o_height->clear();
        o_mpn->clear();
        o_mpi->clear();

		// scan channels
		for(int ch = 0; ch<NCHT; ch++){
			int nPeaks = tdcNhit[ch];
			// reset
			o_pedestal[ch]=0;
			o_areaall[ch]=0;
            temp_clk.clear();temp_clk.resize(nPeaks);
            temp_tdc.clear();temp_tdc.resize(nPeaks);
            temp_peak.clear();temp_peak.resize(nPeaks);
            temp_sum.clear();temp_sum.resize(nPeaks);
            temp_width.clear();temp_width.resize(nPeaks);
            temp_height.clear();temp_height.resize(nPeaks);
            temp_mpn.clear();temp_mpn.resize(nPeaks);
            temp_mpi.clear();temp_mpi.resize(nPeaks);

			// get pedestal
			int clk;
			for(clk = 0; clk<(nPeaks<=0?NSAM:clockNumberDriftTime[ch][0]-1); clk++){
				o_pedestal[ch]+=adc[ch][clk];
			}
			if (clk==0){
				o_pedestal[ch] = prepedestal[ch];
			}
			else {
				o_pedestal[ch]/=clk;
				prepedestal[ch] = o_pedestal[ch];
			}
			// get peaks
			int pre_end_clk = 0;
			int the_mpi = 0;
			for ( int ip = 0; ip<nPeaks; ip++){
				if (driftTime[ch][ip]>0) driftTime[ch][ip]-=power2_15;
				temp_tdc[ip] = driftTime[ch][ip];
				temp_clk[ip] = clockNumberDriftTime[ch][ip];
				temp_height[ip] = adc[ch][temp_clk[ip]];
				if (ip>0&&pre_end_clk>temp_clk[ip]){ // should merge with previous peak in the same packet
				    the_mpi++;
				    temp_peak[ip] = temp_peak[ip-1];
				    temp_sum[ip] = temp_sum[ip-1];
				    temp_width[ip] = temp_width[ip-1];
				    temp_mpn[ip] = temp_mpn[ip];
				    temp_mpi[ip] = the_mpi;
				}
                else{ // new packet
                    the_mpi=0;
                    temp_peak[ip] = 0;;
                    temp_sum[ip] = 0;
                    int clk;
                    for(clk = clockNumberDriftTime[ch][ip]; clk<(clockNumberDriftTime[ch][ip]+MAX_WIDTH>=NSAM?NSAM:clockNumberDriftTime[ch][ip]+MAX_WIDTH); clk++){
                        if (adc[ch][clk]>temp_peak[ip]){
                            temp_peak[ip]=adc[ch][clk];
                        }
                        if (clk!=clockNumberDriftTime[ch][ip]&&adc[ch][clk]<o_pedestal[ch]+0.5) break; // FIXME: should think about the threshold, sigma of threshould is about 0.1 ADC
                        temp_sum[ip] += adc[ch][clk]-o_pedestal[ch];
                    }
                    int mpn = 1;
                    for (int jp = ip+1; jp<nPeaks; jp++){
                        if (clockNumberDriftTime[ch][jp]<clk) mpn++;
                    }
                    temp_width[ip] = clk-clockNumberDriftTime[ch][ip];
                    temp_mpn[ip] = mpn;
                    temp_mpi[ip] = the_mpi;
                    o_areaall[ch] += temp_sum[ip];
                    pre_end_clk = clk;
                }
			}
            o_nPeaks[ch] = nPeaks;
			o_clk->push_back(temp_clk);
			o_tdc->push_back(temp_tdc);
			o_peak->push_back(temp_peak);
			o_sum->push_back(temp_sum);
			o_width->push_back(temp_width);
			o_height->push_back(temp_height);
			o_mpn->push_back(temp_mpn);
			o_mpi->push_back(temp_mpi);
            if (o_nPeaks[ch]) o_nHits++;
		}
		t->Fill();
	}
	t->Write();

	return 0;
}

void print_usage(char* prog_name)
{
	fprintf(stderr,"\t%s [runNo] <[nEventMax]>\n",prog_name);
}
