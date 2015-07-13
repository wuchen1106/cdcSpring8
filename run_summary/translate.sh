#!/bin/bash
log="run.eff.ori"
sed -i "s/:/ /g" $log
sed -i "s/->/ /g" $log
sed -i "s/,/ /g" $log
sed -i "s/\(\<\w*\>\) *= *\(\<\w*\>\) *= *\(\<\w*\>\)/\1 \3 \2 \3/g" $log
sed -i "s/=/ /g" $log
sed -i "s/(/ /g" $log
sed -i "s/)/ /g" $log
sed -i "s/test run/test_run/g" $log
sed -i "s/beam on/beam_on/g" $log
sed -i "s/noise check/noise_check/g" $log
sed -i "s/check efficiency/check_efficiency/g" $log
sed -i "s/check HV/check_HV/g" $log
sed -i "s/\<\([0-9][0-9]*\)mV\>/\1/g" $log
sed -i "s/\<\([0-9][0-9]*\)V\>/\1/g" $log

sed -i "s/@TXY\>/@TXY1\&2/g" $log
sed -i "s/@TXY12\>/@TXY1\&2/g" $log
sed -i "s/@TX12\>/@TXY1\&2/g" $log
sed -i "s/@TY12\>/@TXY1\&2/g" $log
sed -i "s/@TX1\>/@TXY1/g" $log
sed -i "s/@TY1\>/@TXY1/g" $log
sed -i "s/@TX2\>/@TXY2/g" $log
sed -i "s/@TY2\>/@TXY2/g" $log

sed -i "s/@PT23\>/@PT2\&3/g" $log
sed -i "s/@PT\>/@PT2\&3/g" $log
sed -i "s/@prototype-II\>/@PT2/g" $log
sed -i "s/@prototype-III\>/@PT3/g" $log
sed -i "s/Vth for prototype-2&3 is/Vth@PT2\&3/g" $log
sed -i "s/Vth for all is/Vth@PT2\&3/g" $log
sed -i "s/HV of PT2 and PT3 *HV *\([0-9][0-9]*\)/HV HV@PT2\&3 \1/g" $log
sed -i "s/HV of PT2 and PT3 *\([0-9][0-9]*\)/HV HV@PT2\&3 \1/g" $log

sed -i "s/HV of PT2 and PT3 (HV\>/HV@PT2\&3/g" $log
sed -i "s/_PT2/@PT2/g" $log
sed -i "s/_PT3/@PT3/g" $log
sed -i "s/Vth Scan *\([0-9][0-9]*\)/Vth@PT2\&3 \1/g" $log
sed -i "s/\w* Scan$//g" $log

gawk 'BEGIN{tht1=3750;tht2=3750;thp3=3750;thp2=3750;hvt1=1800;hvt2=1800;hvp3=1700;hvp2=1700;}{ \
       for(i=0; i<NF; i++){ \
          if ($i=="HV@TXY1&2"){hvt1=$(i+1);hvt2=$(i+1);} \
          if ($i=="HV@TXY1"){hvt1=$(i+1)} \
          if ($i=="HV@TXY2"){hvt2=$(i+1)} \
          if($i=="HV@PT2&3"){hvp2=$(i+1);hvp3=$(i+1);}\
          if($i=="HV@PT2"){hvp2=$(i+1);}\
          if($i=="HV@PT3"){hvp3=$(i+1);}\
          if($i=="Vth@TXY1&2"){tht1=$(i+1);tht2=$(i+1);}\
          if($i=="Vth@TXY1"){tht1=$(i+1);}\
          if($i=="Vth@TXY2"){tht2=$(i+1);}\
          if($i=="Vth@PT2&3"){thp3=$(i+1);thp2=$(i+1);}\
          if($i=="Vth@PT3"){thp3=$(i+1);}\
          if($i=="Vth@PT2"){thp2=$(i+1);}\
       }\
       print $1" "$2" "$12" "hvt1" "hvt2" "hvp3" "hvp2" "tht1" "tht2" "thp3" "thp2;}'\
       $log
