# cdcSpring8

### Introduction

### Work flow

#### Prepare raw data in ROOT format

Should follow the format of `run_XXXXXX.root`

#### Perform wave analysis

Input file:

`root/run_XXXXXX_built.root`

Output file:

`root/p_XXXX.root`

Usage:

```
getP runNo
```

#### Perform peak selection

Input file:

`Input/run-info.root`  
`Input/wire-position.root`  
`root/run_XXXXXX_built.root`  
`root/p_XXXX.root`  

Output file:

`root/h_XXXX.root`  
and histograms

Usage:

```
getH runNo
```

#### Do the tracking

Input file:

`Input/run-info.root`  
`Input/crosspoint.root`  
`info/wire-position.XXXX.PRERUNNAME.root` or `Input/wire-position.root`  
`info/xt.XXXX.PRERUNNAME.root`  
`info/offset.XXXX.RUNNAME.root` if it exists  
* inputType=1: `root/h_XXXX.root`
* inputType=2: `root/h_XXXX.MC.root`
* inputType=3: `root/h_XXXX.layerX.MC.root`

Output file:

`root/t_XXXX.RUNNAME.layerX.root`  
`info/xt.XXXX.RUNNAME.root`  
`root/ana_XXXX.layerX.MC.root`  

Usage:

```
Usage ./BinaryFiles/bin/tracking [options] prerunname runname
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

#### analyze the reconstructed tracks and get XT relations

Input file:

`info/wire-position.XXXX.PRERUNNAME.root` or `Input/wire-position.root`
`Input/crosspoint.root`
`info/xt.XXXX.PRERUNNAME.root`
* inputType=1: `root/h_XXXX.root`
* inputType=2: `root/h_XXXX.MC.root`
* inputType=3: `root/h_XXXX.layerX.MC.root`

Output file:

``

Usage:

```
Usage ./BinaryFiles/bin/getXT [options] prerunname runname
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
		 Default layer set to l
	 -n <n>
		 Maximum number of hits cut set to n
	 -c <c>
		 Maximum chi2 cut set to c
	 -g <g>
		 Geometry setup set to g
		 (0): normal; 1: finger
	 -i <i>
		 Input type set to i
		 (0) for data; 1 for MC
	 -p <p>
		 Peak type set to p
		 (0) only the first peak over threshold; 1, all peaks over threshold; 2, even including shaddowed peaks
	 -x <x>
		 xt type set to x
		 (2) sym; 1 sym+offset; 0 no; 6 sym+offset+first OT peak; 7 sym+offset+first OT peak+2segments
```

#### get residual distribution from the analyzed file

```
Usage ./BinaryFiles/bin/ana [options] runname
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
	 -f <f>
		 Minimum number of selected hits cut set to f
	 -c <c>
		 Maximum chi2 cut set to c
	 -p <p>
		 Minimum p-value cut set to p
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
	 -x <x>
		 Number of bins on space axis set to x
	 -y <y>
		 Number of bins on resolution axis set to y
	 -g <g>
		 Geometry setup set to g
		 (0): normal; 1: finger
```

#### Test the residual distribution with MC simulation to get tracking error and spatial resolution

```
loopBestRuns.sh list runTag [averageEtrack(1) maxChi2(2) maxSlz(0.1)]
```
which will loop in runNo and runName in list and do the following (`maxNhitsG` is taken from runName "XXX.sXaXnX.iX", while `minNhitsS` is by default 7 for good runs and 5 fow low energy runs)

```
doIter.res.sh runNo runTag iThreadStart nThreads iIterStart iIterStop runName [averageEtrack(1) maxChi2(2) minNhitsS(7) maxNhitsG(0) maxSlz(0.1)]
```

which generate an initial resolution file (`info/res.XXX`) first:

```
updateRes runNo StartName preRunName 1 averageEtrack
```

and then loop from `iIterStart` to `iIterStop` to do:

```
mc2fitinput_Chen root/ana_XXX info/res.XXX root/h_XXX.MC.root 0 0 0 nLayers 1 1 iStart iStop maxChi2 minNhitsS maxNhitsG maxSlz
hadd root/h_runNo.MC.root root/h_runNo.iStart-iStop.MC.root ...
tracking
combine runNo currunname nEvtPerRun h_runNo.MC.root
updateRes runNo StartName preRunName 0 averageEtrack
mv root/h_runNo.MC.root root/h_runNo.prerunname.MC.root
```
