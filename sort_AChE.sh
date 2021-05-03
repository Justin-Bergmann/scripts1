#!/bin/bash


if [ -z "$1" ]
then
      pdb=sortet_files
else
      pdb=$1
fi

name=${pdb}_sort
echo $name
rm -r $name
rm -r ${name}.zip
mkdir $name
cd $name

m1=1
mmax=5 # loops through method m
method[1]='P00'
method[2]='P03'
method[3]='P05'
method[4]='P12'
method[5]='P13'
echo 0
pwd
for (( m=${m1}; m<=${mmax}; m++ ))
do
   mkdir -p ${method[$m]}
   echo 1
   pwd
   cd ..
   echo 2
   pwd
   cd ${method[$m]}
   mimic > pdb 
   pdb_name=${method[$m]}.pdb
   pdb_name_qm=${method[$m]}_qm.pdb
   mtz_name=${method[$m]}.mtz
   FO_name=${method[$m]}_2mFo-DFc_map.ccp4
   FC_name=${method[$m]}_mFo-DFc_map.ccp4
   cp bindividual1.pdb ../$name/${method[$m]}/$pdb_name
   cd Gz
   gzip -d qm.txt.gz comqum.pdb.gz
   gen_syst1.py qm.txt comqum.pdb
   cp label_vis.sh ..
   cd ..
   chmod +x label_vis.sh
   ./label_vis.sh
   cp pdb ../$name/${method[$m]}/$pdb_name_qm
   cp edstats/pdb_map_coeffs.mtz ../$name/${method[$m]}/$mtz_name
   cp edstats/pdb_2mFo-DFc_map.ccp4 ../$name/${method[$m]}/$FO_name
   cp edstats/pdb_mFo-DFc_map.ccp4 ../$name/${method[$m]}/$FC_name
   echo 3
   pwd
   cd ../$name
done
cd ..
zip_name=${name}.zip 
zip -r $zip_name $name
