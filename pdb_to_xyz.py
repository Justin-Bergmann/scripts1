#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### run :
# python3 d_hbond.py file.pdb



import time
## set timer
time_start = time.time()
import sys
##load mandatory modules
import numpy as np
import copy
from collections import OrderedDict
## define vdw parameter in AA from Bondi, A.V., 1964. van der Waals volumes and radii. The Journal of physical chemistry, 68(3), pp.441-451.

def read_pdb(input_pdb):
## reads in pdb coordinates hardcodedt 
##https://zhanglab.ccmb.med.umich.edu/BindProfX/pdb_atom_format.html
##COLUMNS        DATA TYPE       CONTENTS                            
##--------------------------------------------------------------------------------
 ##1 -  6        Record name     "ATOM  "                                            
## 7 - 11        Integer         Atom serial number.                   
##13 - 16        Atom            Atom name.                            
##17             Character       Alternate location indicator.         
##18 - 20        Residue name    Residue name.                         
##22             Character       Chain identifier.                     
##23 - 26        Integer         Residue sequence number.              
##27             AChar           Code for insertion of residues.       
##31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
##39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
##47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
##55 - 60        Real(6.2)       Occupancy.                            
##61 - 66        Real(6.2)       Temperature factor (Default = 0.0).                   
##73 - 76        LString(4)      Segment identifier, left-justified.   
##77 - 78        LString(2)      Element symbol, right-justified.      
##79 - 80        LString(2)      Charge on the atom.
        print(" read pdb file")
        read_pdb = open(input_pdb, "r")
        org_pdb =[]
        work_pdb =[]
        for line in read_pdb:
            idef = line[0:6].strip()
            if idef=="ATOM" or idef=="HETATM":

                l_pdb=[]
                a       = line[0:6].strip()
                a_num   = int(line[6:11].strip())
                a_name  = line[12:16].strip()
                x       = float(line[30:38].strip())
                y       = float(line[38:46].strip())
                z       = float(line[46:54].strip())

                l_pdb.append(a)
                l_pdb.append(a_num)
                l_pdb.append(a_name)
                l_pdb.append(x)
                l_pdb.append(y)
                l_pdb.append(z)

                org_pdb.append(l_pdb)
                work_pdb.append(l_pdb)
#               print(l_pdb)

        print("close file:", input_pdb)
        return org_pdb, work_pdb




def write_xyz(output, pdb):
        out=  open(output, "w")
        string= str(len(pdb))+" \n"+" \n"
        out.write(string)

        for i in range(len(pdb)):
           line=str(pdb[i][2])+" "+str(pdb[i][3])+" "+str(pdb[i][4])+" "+str(pdb[i][5])+" \n"
           out.write(line)

        print("xyz file is written")

        out.close()



 
################################################
##stript starts


if len(sys.argv)>=2:
        input_pdb=sys.argv[1]
        print("read from ", input_pdb)
else:
        input_pdb="pdb"
        print("read from defalt ", input_pdb)
output="pdb.xyz"
## [[0,   1,      2,    3,      4          5       6
## [[a, a_num, a_name, alter, res_name,  chain, res_num,
##   7    8  9 10  11   12      13     14    15
## inser, x, y, z, occ, b_fac, seg_i, ele, charge ]]        
org_pdb, pdb = read_pdb(input_pdb)
print(len(pdb))
write_xyz(output, pdb)




time_ende = time.time()
print("program ends normally after "'{:5.3f}s'.format(time_ende-time_start),
" or ", '{:5.2f}min'.format((time_ende-time_start)/60))

