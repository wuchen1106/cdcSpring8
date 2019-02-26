# cdcSpring8

### Introduction

### Work flow

#### Prepare raw data in ROOT format

Should follow the format of `run_XXXXXX.root`

#### Perform wave analysis

Usage:

```
getP runNo
```

Input file:

`root/raw/run_XXXXXX_built.root`

Output file:

`root/peaks/p_XXXX.root`

Tree structure of output file:
```
int triggerNumber;
int nh;
int np[NCHT];
double ped[NCHT];
double aa[NCHT];
std::vector<std::vector<int> > * clk; // the first level of vector has the same size as NCHT (total number of channels). The second level of vector has the size of number of TDCs in this channel in this event.
std::vector<std::vector<int> > * tdc;
std::vector<std::vector<int> > * peak;
std::vector<std::vector<double> > * sum;
std::vector<std::vector<int> > * width;
std::vector<std::vector<int> > * height;
std::vector<std::vector<int> > * mpn;
std::vector<std::vector<int> > * mpi;
```

#### Perform peak selection

Usage:

```
getH runNo
```

Input file:

`Input/run-info.root`  
`Input/wire-position.root`  
`root/raw/run_XXXXXX_built.root`  
`root/peaks/p_XXXX.root`  

Output file:

`root/hits/h_XXXX.root`  
and histograms

Tree structure of output file:  
```
int triggerNumber;
int nHits; // number of hits. Size is the sum of TDC numbers in all channels
int nLayers; // number of layers with at least one hit
std::vector<int> * layerID; // vector of hits. Sequence decided by DAQ.
std::vector<int> * wireID;
std::vector<int> * type; // denotes the region (R) in the chamber of the cell: 1: 4 dummy; 3 guard; 2 right boundary; 1 left boundary; 0 center part
std::vector<int> * np;
std::vector<int> * ip;
std::vector<int> * clk;
std::vector<int> * width;
std::vector<int> * peak;
std::vector<int> * height;
std::vector<int> * mpn;
std::vector<int> * mpi;
std::vector<int> * rank;
std::vector<double> * ped;
std::vector<double> * sum;
std::vector<double> * aa;
std::vector<double> * driftT;
```

#### Do the tracking

Usage:

```
Usage tracking [options] prerunname runname
[options]
	 -D <name>=[error,severe,warn,debug,trace]
		 Change the named debug level
	 -V <name>=[quiet,log,info,verbose]
		 Change the named log level
	 -C <file>
		 Set the configure file
	 -M <n>
		 Printing modulo set to n
	 -R <run>
		 Run number set to run
	 -B <n>
		 Starting entry index set to n
	 -E <n>
		 Stopping entry index set to n
	 -L <l>
		 Test layer set to l
	 -n <n>
		 Maximum number of hits cut set to n
	 -x <x>
		 T0 shift on board 0 set to x
	 -y <y>
		 T0 shift on board 1 set to y
	 -l <l>
		 Minimum time cut set to l
	 -u <u>
		 Maximum time cut set to u
	 -s <s>
		 ADC sum over peak cut set to s
	 -a <a>
		 ADC sum over all cut set to a
	 -g <g>
		 Geometry setup set to g
		 (0): normal; 1: finger
	 -i <i>
		 Input type set to i
		 (0) for data; 1 for MC
	 -p <p>
		 Peak type set to p
		 (0) only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
	 -w <w>
		 Work type set to w
		 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
```

Input file:

`Input/run-info.root`  
`Input/crosspoint.root`  
`info/wire-position.XXXX.PRERUNNAME.root` or `Input/wire-position.root`  
`info/xt.XXXX.PRERUNNAME.root` or `info/xt.GAS.HV.root`
`info/reso.XXXX.PRERUNNAME.root` or `info/reso.GAS.HV.root`
`info/offset.XXXX.RUNNAME.root` if it exists  
* inputType=0: `root/h_XXXX.root` for data
* inputType=2: `root/h_XXXX.MC.root` for MC without layer specified
* inputType=3: `root/h_XXXX.layerX.MC.root` for MC with layer specified

Output file:

`root/tracks/t_XXXX.RUNNAME.layerX.root`  

Tree structure of output file:
```
    int                   triggerNumber;
    // inheriting from raw hits
    int                   nHits;
    std::vector<int> *    layerID;
    std::vector<int> *    wireID;
    std::vector<double> * driftT; 
    std::vector<double> * driftDmc; // if MC sample
    std::vector<int> *    type; // in dec, [IMASTR].
    // I: peak index (only counting peaks over m_sumCut);
    // M: peak index in a packet;
    // A: smaller than aa cut?
    // S: smaller than sum cut?
    // T: 0 good; 1 l-r0; 2 l0r-; 3 l-r-; 4 l+r0; 5 l0r+; 6 l+r+; 7 else; l/r means left/right XT function. "-/+" means drift time is lower/upper to the limited region. "0" means in the region.
    // R: 0 center, 1 left, 2 right, 3 guard, 4 dummy
    std::vector<int> *    np;
    std::vector<int> *    ip;
    std::vector<int> *    clk;
    std::vector<int> *    width;
    std::vector<int> *    peak;
    std::vector<int> *    height;
    std::vector<int> *    mpn;
    std::vector<int> *    mpi;
    std::vector<int> *    rank;
    std::vector<double> * ped;
    std::vector<double> * sum;
    std::vector<double> * aa;
    // basic
    int                   nHitsG; // number of good hits in layers other than the test one: in t region and with good peak quality
    std::vector<double> * dxl;
    std::vector<double> * dxr;
    // track finding/fitting with different candidates;
    int                   nFit;
    int                   nFind;
    std::vector<double> * driftD[NCAND];
    std::vector<double> * calD[NCAND];
    std::vector<double> * fitD[NCAND];
    std::vector<int>    * sel[NCAND];
    int                   icom[NCAND];
    int                   isel[NCAND];
    int                   npairs[NCAND];
    double                iinx[NCAND];
    double                iinz[NCAND];
    double                islx[NCAND];
    double                islz[NCAND];
    double                chi2x[NCAND];
    double                chi2z[NCAND];
    double                chi2i[NCAND];
    double                chi2pi[NCAND];
    double                chi2ai[NCAND];
    int                   nHitsS[NCAND]; // number of hits selected from finding and fed to fitting
    double                inx[NCAND];
    double                inz[NCAND];
    double                slx[NCAND];
    double                slz[NCAND];
    double                chi2[NCAND];
    double                chi2p[NCAND];
    double                chi2a[NCAND];
    // if this is MC sample
    double                chi2mc[NCAND];
    double                chi2pmc[NCAND];
    double                chi2amc[NCAND];
    double                inxmc;
    double                inzmc;
    double                slxmc;
    double                slzmc;
```

#### analyze the reconstructed tracks and get XT relations

Usage:

```
Usage ana [options] prerunname runname [wires to be calibrated (all)]
[options]
	 -D <name>=[error,severe,warn,debug,trace]
		 Change the named debug level
	 -V <name>=[quiet,log,info,verbose]
		 Change the named log level
		 If equal sign is not found, set verbose level to the given value
	 -C <file>
		 Set the configure file
	 -M <n>
		 Printing modulo set to n
	 -R <run>
		 Run number set to run
	 -B <n>
		 Starting entry index set to n
	 -E <n>
		 Stopping entry index set to n
	 -L <l>
		 Default layer set to l
	 -H <h>
		 Histogram saving level set to h
	 -T
		 Disable in-trigger cut
	 -G
		 In getting offset and new XT, require all hits used in tracking are golden hits
	 -F <suffix>
		 Adding suffix to output files
	 -W
		 (false) Update wire map
	 -X <xt file>
		 Use given xt file instead of generating a new one
	 -n <n>
		 Maximum number of hits cut set to n
	 -f <f>
		 Minimum number of selected hits cut set to f
	 -c <c>
		 Maximum chi2 cut set to c
	 -v <v>
		 Minimum p-value cut set to v
	 -r <r>
		 Maximum resolution cut set to r
	 -z <z>
		 Maximum y-z slope cut set to z
	 -d <d>
		 Maximum fitD cut set to d
	 -o <o>
		 Maximum time range set to o
	 -s <s>
		 ADC sum over peak cut set to s
	 -a <a>
		 ADC sum over all cut set to a
	 -l <l>
		 Minimum time on axis set to l
	 -u <u>
		 Maximum time on axis set to u
	 -t <t>
		 Number of bins on time axis set to t
	 -m <m>
		 Number of bins on space axis set to m
	 -y <y>
		 Number of bins on resolution axis set to y
	 -g <g>
		 Geometry setup set to g
		 (0): normal; 1: finger
	 -i <i>
		 Input type set to i
		 (0) for data; 1 for MC
	 -x <x>
		 xt type set to x
		 XYZ (055) means polX for center, polY for middle and polZ for tail. If X is 0 then let middle function fit the center region.
	 -A
		 Use asymmetric XT
	 -S
		 Select candidate by:
		 ((O)riginal): the first one given by tracking (global chi2)
		 (G)lobalChi2: find the candidate with the smallest global chi2 including the test layer
		 (F)ittingChi2: find the candidate with the smallest fitting chi2
		 (L)eastLatePeak: find the candidate with the least late arrival peaks
	 -P
		 Use the peak with smallest residual to get driftT. By default just use the first one
```

Input file:

`info/wire-position.XXXX.PRERUNNAME.root` or `Input/wire-position.root`
`Input/crosspoint.root`
`info/xt.XXXX.PRERUNNAME.root`  
`root/tracks/t_XXXX.RUNNAME.layerX.root`  

Output file:

`info/resi.XXXX.RUNNAME.layerX.root`  
`info/offset.XXXX.RUNNAME.layerX.root`  
`info/wire-position.XXXX.RUNNAME.layerX.root` if the `m_UpdateWireMap` option is set  
`info/xt.XXXX.RUNNAME.root` if the `m_ExternalXT` option is not set  
`root/ana/ana_XXXX.layerX.root`  
 
 with branch structure
 ```
    int                   triggerNumber;
    // test layer relavant info
    double                res;
    bool                  isGood; // if this event meets all cuts
    int                   nHitsT[2]; // record 2 closest cells
    bool                  hashit[2];
    double                theDD;
    double                theDT;
    double                theFD[2];
    int                   theST[2];
    double                theAA[2];
    double                theCharge; // charge in the test layer
    int                   theWid;
    double                theSum;
    double                sum1st;
    double                dt1st;
    double                thePeak;
    double                theHeight;
    int                   theIp;
    int                   theMpi;
    int                   theCand;
    int                   highBid;
    int                   highCh;
    int                   highLid;
    int                   highWid;
    int                   highIp;
    double                highSum;
    double                highAA;
    double                highDT;
    // dE/dX related
    int                   nLayers; // number of layers used for chage on track
    double                chargeOnTrack[NLAY]; // charge along the track
    double                adcsumOnTrack[NLAY]; // ADC sum along the track
    double                chargeOnTrackIndex[NLAY]; // index of the corresponding hit along the track
    double                theGG;
    double                trackGG;
    // special hits count
    int                   nHitsSmallAll;
    int                   nHitsSmallSASD;
    int                   nSHits;
    int                   nLHits;
    int                   nSSHits;
    int                   nBHits;
    int                   nSBHits;
    // inheriting from raw hits
    int                   nHits;
    std::vector<int> *    layerID;
    std::vector<int> *    wireID;
    std::vector<int> *    channelID; // newly added info
    std::vector<int> *    boardID; // newly added info
    std::vector<double> * driftT; 
    std::vector<double> * driftDmc; // if MC sample
    std::vector<int> *    type; // in dec, [IMASTR]. I: peak index (only counting peaks over m_sumCut); M: peak index in a packet; A: smaller than aa cut? S: smaller than sum cut? T: -1 <m_tmin, 0 good, 1 >m_tmax; R: 0 center, 1 left, 2 right, 3 guard, 4 dummy
    std::vector<int> *    np;
    std::vector<int> *    ip;
    std::vector<int> *    clk;
    std::vector<int> *    width;
    std::vector<int> *    peak;
    std::vector<int> *    height;
    std::vector<int> *    mpn;
    std::vector<int> *    mpi;
    std::vector<int> *    rank; // if input has it
    std::vector<double> * ped; // if input has it
    std::vector<double> * sum;
    std::vector<double> * aa;
    // inheriting from tracking
    int                   nHitsG; // number of good hits in layers other than the test one: in t region and with good peak quality
    std::vector<int> *    driftDs;
    //    Currently we only choose one candidate to analyze, so the output is simpler
    std::vector<double> * driftD;
    std::vector<double> * driftD0; // old drift D from tracking
    std::vector<double> * calD;
    std::vector<double> * fitD;
    std::vector<int>    * sel;
    int                   icom;
    int                   isel;
    int                   npairs;
    double                iinx;
    double                iinz;
    double                islx;
    double                islz;
    double                chi2x;
    double                chi2z;
    double                chi2i;
    int                   nHitsS; // number of hits selected from finding and fed to fitting
    double                inx;
    double                inz;
    double                slx;
    double                slz;
    double                chi2;
    double                chi2p;
    double                chi2a;
    // if MC sample
    double                chi2mc;
    double                chi2pmc;
    double                chi2amc;
    double                inxmc;
    double                inzmc;
    double                slxmc;
    double                slzmc;
 ```
 
#### Test the residual distribution with MC simulation to get tracking error and spatial resolution

```
Usage updateRes [options] prerunname runname
[options]
	 -D <name>=[error,severe,warn,debug,trace]
		 Change the named debug level
	 -V <name>=[quiet,log,info,verbose]
		 Change the named log level
		 If equal sign is not found, set verbose level to the given value
	 -C <file>
		 Set the configure file
	 -M <n>
		 Printing modulo set to n
	 -R <run>
		 Run number set to run
	 -B <n>
		 Starting entry index set to n
	 -E <n>
		 Stopping entry index set to n
	 -L <l>
		 Default layer set to l
	 -A
		 Use averaged tracking error
	 -I
		 This is the initial iteration step
	 -n <n>
		 Maximum number of hits cut set to n
	 -c <c>
		 Maximum chi2 cut set to c
	 -d <d>
		 Maximum fitD cut set to d
	 -c <c>
		 Maximum chi2 cut set to c
	 -g <g>
		 Geometry setup set to g
		 (0): normal; 1: finger
	 -i <i>
		 Input type set to i
		 (0) for data; 1 for MC; 2 for MC with real driftD information
	 -x <x>
		 xt type set to x
		 (2) sym; 1 sym+offset; 0 no; 6 sym+offset+first OT peak; 7 sym+offset+first OT peak+2segments

```

Input file:

`info/resi.XXXX.RUNNAME.layerX.root`  
`info/xt.XXXX.ORIRUNNAME.root`  
`root/tracks/t_XXXX.RUNNAME.layerX.root`  

Output file:

`info/reso.XXXX.RUNNAME.layerX.root`  
`root/anamc/anamc_XXXX.layerX.MC.root`  

#### To do the iteration on real data

```
doIter.sh:  do the iteration!
Syntax:
    doIter.sh [options] runname
    [options]
    -h     display this help and exit
    -R [R] the run number
    -T [T] the staring thread index
    -N [N] the number of threads to ask for
    -I [I] start from iteration I
    -J [J] stop at iteration J
    -H [h] save histograms at this level (0)
    -L     (false) Take the last iteration as the final step and do the default complete checkings
    -W     (false) update the wire position map
    -U     (false) keep the xt curves unchanged.
    -S [S] set the start name (Garfield)
    -D [D] set this layer (4) as the default layer to save in XT file as fr/l_0
    -l [l1 (l2 ...)] Do the training of the given layers (4)
    -w [w1 (w2 ...)] Do the calibration of the given wires () in the ubove given layers (4)
    -c [c] Add configure file as C (default one "Para/default.dat" is always loaded first)
    -g [g] geometry setup (0). 0 ordinary scintillator; 1 finger scintillator
    -t [t] work type (0) for tracking. 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
    -a [a] aa cut (0)
    -s [s] sum cut (-10)
    -p [p] set peak type (0) for tracking. 0, only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
    -x [XYZ] xt type. XYZ (055) means polX for center, polY for middle and polZ for tail. If X is 0 then let middle function fit the center region.
    -n [n] maximum number (30) of hits to be used in ana
    -m [m] maximum number (15) of good hits to be used in tracking
```

#### To do the iteration to get spatial resolution with MC

```
doIter.res.sh runNo runTag iThreadStart nThreads iIterStart iIterStop runName [averageEtrack(1) maxChi2(2) minNhitsS(7) maxNhitsG(0) maxSlz(0.1)]
```

which generate an initial resolution file (`info/res.XXX`) first:

```
updateRes -R runNo -A -I averageEtrack StartName preRunName
```

and then loop from `iIterStart` to `iIterStop` to do:

```
mc2fitinput_Chen root/ana/ana_XXX info/res.XXX root/hits/h_XXX.MC.root 0 0 0 nLayers 1 1 iStart iStop maxChi2 minNhitsS maxNhitsG maxSlz
hadd root/hits/h_runNo.MC.root root/hits/h_runNo.iStart-iStop.MC.root ...
tracking
combine runNo currunname nEvtPerRun h_runNo.MC.root
updateRes -R runNo -A StartName preRunName
mv root/hits/h_runNo.MC.root root/hits/h_runNo.prerunname.MC.root
```
#### To loop in a list of runs

```
loopBestRuns.sh list runTag [averageEtrack(1) maxChi2(2) maxSlz(0.1)]
```
which will loop in runNo and runName in list and do the following (`maxNhitsG` is taken from runName "XXX.sXaXnX.iX", while `minNhitsS` is by default 7 for good runs and 5 fow low energy runs)

```
doIter.res.sh runNo runTag iThreadStart nThreads iIterStart iIterStop runName [averageEtrack(1) maxChi2(2) minNhitsS(7) maxNhitsG(0) maxSlz(0.1)]
```
