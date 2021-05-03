#!/bin/bash
comqumtocns <<EOF
x
mm1.par
all.par



EOF

changepdb <<EOF
mm3.pdb
cqf
comqum.dat

w
fix.pdb
q
EOF
  
