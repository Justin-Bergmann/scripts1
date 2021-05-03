#!/bin/bash

#!/bin/bash


if [ -z "$1" ]
then
      pdb=pdb
else
      pdb=$1
fi



ncns <generate.inp >generate.out
grep ERR generate.out
grep WRN generate.out
grep INFO generate.out
grep error generate.out
grep warning generate.out 


cns <alternate.inp >alternate.out


changepdb << EOF
alternate.pdb
al
$pdb
w
mm3.pdb
q
EOF
\mv alternate.mtf mm3.mtf

cp mm3.pdb comqum.pdb


