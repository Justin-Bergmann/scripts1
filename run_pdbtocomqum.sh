#!/bin/bash


if [ -z "$1" ]
then
      proto=N
else
      proto=$1
fi


module load GCCcore/8.3.0
module load Python/3.7.4

changepdb << EOF
comqum.pdb
ren

w

q
EOF



/home/justin/snic2019-35-66/Justin/programs/scripts/gen_syst1.py qm.txt comqum.pdb $proto
pdbtocomqum << EOF
logfil
comqum.pdb
n

syst1
0
1000
N
N
0
c
12
new
EOF

comqumtoturbo

mimic > pdb




















