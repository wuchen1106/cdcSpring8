#!/bin/bash

QUEUE="h"
MSGDIR=/dev/null
#MSGDIR=$PWD/root
layers=""
geoSetup=0 # 1: finger; 0: normal
inputType=0 # 1: MC; 0: Real Data
workType=0 # 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
nHitsMax=13
t0shift=0
tmin=-20
tmax=800
sumCut=-20
aaCut=20
debug=-1

THISCMD=$(basename $0)
usages() {
cat << EOF
$THISCMD:  Generate simulation/reconstruction jobs and submit to pbs
Syntax:
	$THISCMD [options] runNo N prerunname runname [layers]
Options:
	-h  display this help and exit
	-g  geoSetup [$geoSetup] 1: finger; 0: normal
	-i  inputType [$inputType] 1: MC; 0: Real Data
	-w  workType [$workType] 0, fr/l_0; 1, even/odd; -1, even/odd reversed; others, all layers
	-n  nHitsMax [$nHitsMax]
	-t  t0shift [$t0shift]
	-a  tmin [$tmin]
	-z  tmax [$tmax]
	-s  sumCut [$sumCut]
	-q  aaCut [$aaCut]
	-d  debug [$debug]
EOF
}

while getopts ':hg:i:w:n:t:a:z:s:q:d:' optname
do
	case "$optname" in
	'h')
		usages
		exit 0
		;;
	'g')
		geoSetup="$OPTARG"
		;;
	'i')
		inputType="$OPTARG"
		;;
	'w')
		workType="$OPTARG"
		;;
	'n')
		nHitsMax="$OPTARG"
		;;
	't')
		t0shift="$OPTARG"
		;;
	'a')
		tmin="$OPTARG"
		;;
	'z')
		tmax="$OPTARG"
		;;
	's')
		sumCut="$OPTARG"
		;;
	'q')
		aaCut="$OPTARG"
		;;
	'd')
		debug="$OPTARG"
		;;
	'?')
		echo "Unknown option $OPTARG"
		usages
		exit -1
		;;
	':')
		echo "Need argument value for option $OPTARG"
		usages
		exit -1
		;;
	*)
		echo 'Unknown error while processing options'
		usages
		exit -1
		;;
	esac
done

if [ ! $(($#+1-$OPTIND)) -lt 4 ]
then
    runNo=${@:$OPTIND:1}
    N=${@:$OPTIND+1:1}
    prerunname=${@:$OPTIND+2:1}
    runname=${@:$OPTIND+3:1}
else
    usages
	exit -1
fi

for (( i=$(($OPTIND+4)); i<=$#; i++ ))
do
    layers=$layers" "${!i}
done

for i in $layers;
do
    testlayer=$i
    for (( j=0; j<N; j+=5000 ))
    do
        iEntryStart=$j
        iEntryStop=$((j+4999))
        if (( iEntryStop>=N ))
        then
            iEntryStop=$((N-1))
        fi
        trunname="${runname}.$iEntryStart-$iEntryStop"
        log="$PWD/root/t_${runNo}.${trunname}.layer${testlayer}.log"
        err="$PWD/root/t_${runNo}.${trunname}.layer${testlayer}.err"
        bsub -o $MSGDIR -e $MSGDIR -q $QUEUE "tracking $runNo $testlayer $prerunname $trunname $nHitsMax $t0shift $tmin $tmax $geoSetup $sumCut $aaCut $iEntryStart $iEntryStop $workType $inputType $debug > $log 2> $err"
    done
done
