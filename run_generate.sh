#!/bin/bash
cns  <generate.inp >generate.out
grep ERR     generate.out
grep WRN     generate.out
grep INFO    generate.out
grep error   generate.out
grep warning generate.out
