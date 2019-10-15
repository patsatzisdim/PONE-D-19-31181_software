#!/bin/bash
#    comments
rm -rf ATrend
mkdir ATrend
#
DIR1="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
input1=$DIR1/paramet.i
input2=$DIR1/SUB_main.f
YYin=$DIR1/BLM_ICs.f 


#for (( i = 0; i < 20; i++ )); do
var_LacInit1=$(grep "c2=" $input1)
#B=1.0
#D=$(echo "($i+2) * $B" | bc -l)
#C=$(echo "($i+1) * $B" | bc -l)
#var_LacCha="${var_LacInit/$C/$D}"
#
#sed -i -e "s/$var_LacInit/$var_LacCha/" $input1
#done
#exit 1

#sed -i -e "s/$var_LacCha/$var_LacInit/" $input1





#
# loop over the change of initial conditions
for (( i = 0; i < 21; i++ )); do
var_LacInit=$(grep "c2=" $input1)
B=1.0
D=$(echo "($i+2) * $B" | bc -l)
C=$(echo "($i+1) * $B" | bc -l)
var_LacCha="${var_LacInit/$C/$D}"

# make dir for diagns
rm -rf ADiag
mkdir ADiag
#
#                                  First, the system attains it's SS for 11020 s.
rm *.o; rm *.dat; make; ./main
#
#                                  Non Steady-State (Do not cancel j8 from now on)
var_SS=$(grep "fac=" $input1)
var_NSS="${var_SS/1.0d0/1.0d0}" 
sed -i -e "s/$var_SS/$var_NSS/" $input1
#
#                                  Changes on Rates to sjow the validity of diagnostics
#var_RC=$(grep "FM=" $input1)
#var_RC1="${var_RC/1.0d0/1.2d0}" 
#sed -i -e "s/$var_RC/$var_RC1/" $input1
#                                  Oxidative (OCF) and glycolytic (GCF) pulse till 11040 s. 
var_tIn=$(grep "tend=" $input2)
var_tOut="${var_tIn/11.02d+3/2.0d+1}" 
varIn=$(grep "GCF=" $input1)
varOut="${varIn/+0.0d0/+0.3d-1}"
var1In=$(grep "OCF=" $input1)
var1Out="${var1In/+0.0d0/-1.58d-2}"
#
#var_NoEMin=$(grep "NoEM=9 " $input2)
#var_NoEMout="${var_NoEMin/=9 /=8 }" 
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$varIn/$varOut/" $input1
sed -i -e "s/$var1In/$var1Out/" $input1
sed -i -e "s/$var_tIn/$var_tOut/" $input2 
#sed -i -e "s/$var_NoEMin/$var_NoEMout/" $input2 
#
# Change the BLM_ICs.f to take new ICs
ICflagIn=$(grep "iflag=1" $YYin)
ICflagOut="${ICflagIn/1/2}"
sed -i -e "s/$ICflagIn/$ICflagOut/" $YYin
#
rm *.o; make; ./main
#
#                                  Change on glucolytic (GCF) pulse till 11060 s.
varOut1="${varOut/+0.3d-1/+1.5d-2}"
#
# Change paramet.i 
sed -i -e "s/$varOut/$varOut1/" $input1
#
rm *.o; make; ./main
#
#                                  Take out one part of oxidative (OCF) pulse till 11140 s.
var_tOut1="${var_tOut/2.0d+1/8.0d+1}" 
var1Out1="${var1Out/-1.58d-2/-2.8d-3}"
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$var1Out/$var1Out1/" $input1
sed -i -e "s/$var_tOut/$var_tOut1/" $input2 
#
rm *.o; make; ./main
#
#                                  Add a part on oxidative (OCF) pulse till 11240 s.
var_tOut2="${var_tOut1/8.0d+1/1.0d+2}" 
var1Out2="${var1Out1/-2.8d-3/+1.2d-3}"
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$var1Out1/$var1Out2/" $input1
sed -i -e "s/$var_tOut1/$var_tOut2/" $input2 
#
rm *.o; make; ./main
#
#                                  Take out a part of oxidative (OCF) pulse till 11320 s.
# var_tOut1 is the variable I want for time
var1Out3="${var1Out2/+1.2d-3/+4.0d-3}"
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$var1Out2/$var1Out3/" $input1
sed -i -e "s/$var_tOut2/$var_tOut1/" $input2 
#sed -i -e "s/$var_NoEMout/$var_NoEMin/" $input2 
#
rm *.o; make; ./main
#
#                                  Return to the initial values till 11340 s.
sed -i -e "s/$var_tOut1/$var_tOut/" $input2 
sed -i -e "s/$varOut1/$varIn/" $input1
sed -i -e "s/$var1Out3/$var1In/" $input1
#sed -i -e "s/$var_NoEMin/$var_NoEMout/" $input2 
#
rm *.o; make; ./main
#
#                                  Oxidative (OCF) and Glycolytic (GCF) pulse till 11460 s
var_tOut2="${var_tOut/2.0d+1/1.2d+2}" 
var1Out4="${var1Out3/+4.0d-3/-4.15d-2}"
varOut2="${varOut1/+1.5d-2/-2.9d-2}"
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$var1In/$var1Out4/" $input1
sed -i -e "s/$varIn/$varOut2/" $input1
sed -i -e "s/$var_tOut/$var_tOut2/" $input2 
#
rm *.o; make; ./main
#
#                                  Take out the Glycolytic (GCF) pulse and change Oxidative (OCF) till 11540 s
#var_tOut1 is the time I want
var1Out5="${var1Out4/-4.15d-2/+8.0d-3}"
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$var1Out4/$var1Out5/" $input1
sed -i -e "s/$varOut2/$varIn/" $input1
sed -i -e "s/$var_tOut2/$var_tOut1/" $input2 
#
rm *.o; make; ./main
#
#                                   Return to the initial values till 1600 s (end of domain)-I will leave it to 3000 s.
var_tOut3="${var_tOut2/1.2d+2/1.96d+3}" 
#
# Change paramet.i and SUB_main.f
sed -i -e "s/$var1Out5/$var1In/" $input1
sed -i -e "s/$var_tOut1/$var_tOut3/" $input2 
#
#                                   Steady-State (Cancel j8 from now on)
sed -i -e "s/$var_NSS/$var_SS/" $input1
#
#                                  Return to initial values of Rates
#sed -i -e "s/$var_RC1/$var_RC/" $input1
#
rm *.o; make; ./main
#
# Before close change to initial values for the next run
sed -i -e "s/$ICflagOut/$ICflagIn/" $YYin
sed -i -e "s/$var_tOut3/$var_tIn/" $input2
sed -i -e "s/$varOut2/$varIn/" $input1
sed -i -e "s/$var1Out5/$var1In/" $input1
#
#new condition
sed -i -e "s/$var_LacInit/$var_LacCha/" $input1

done
sed -i -e "s/$var_LacCha/$var_LacInit1/" $input1

#echo $input