# cdcSpring8

### Introduction

### Work flow

#### Prepare raw data in ROOT format

Should follow the format of `run_XXXXXX.root`

#### Perform wave analysis

```
getP runNo
```

#### Perform peak selection

```
getH runNo
```

#### Do the tracking

```
tracking runNo
```

#### analyze the reconstructed tracks and get XT relations

```
getXT runNo
```

#### get residual distribution from the analyzed file

```
ana runNo
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
