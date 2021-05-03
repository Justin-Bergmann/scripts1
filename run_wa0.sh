#!/bin/bash
#script to setup wa=0 caclulation in subfolder wa0
rm -rf wa0
echo "wa=0 caluslation is tstartet"
wa0="0.0000;"
mkdir wa0
cp -p * wa0
cd wa0
wa=` grep -m 1 "wa=" mmrun2 `
wa1=${wa##*=}
#echo $wa1
echo $wa 
sed -i "s/$wa1/$wa0/g" mmrun2
wa=` grep -m 1 "wa=" mmrun2 `
echo $wa
#q=q_*

q=`ls * | grep "q_" ` 

#echo $q
q1=${q#*_} 
#echo $q1
q2="q_0$q1"
#echo $q2

mv $q $q2

sbatch $q2
echo "wa=0 calc startet"




