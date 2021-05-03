#!/bin/bash

echo "run makebs < makebs.in > makebs.out"
makebs < makebs.in > makebs.out

echo " "
echo "from control file:"
grep -A 1 '$alpha shells' control
grep -A 1 '$beta shells' control

echo " "
echo "from makebs:"
grep '$alpha shells  a' makebs.out
grep '$beta shells   a' makebs.out











