#!/usr/bin/env python3


import time
## set timer
time_start = time.time()
import argparse, sys
##load mandatory modules
import numpy as np
#import scipy as sc

parser=argparse.ArgumentParser()
parser.add_argument('--pdb', help='input pdb file')
parser.add_argument('--list', help='input list')
parser.add_argument('--out', help='output file')
args=parser.parse_args()

if args.pdb != None:
   input_pdb=args.pdb
   print("input pdb file is:", input_pdb)
else:
         input_pdb="pdb"
         print("input pdb file is:", input_pdb)


if args.list != None:
           input_list=args.list
           print("read from ", input_list)
else:
           input_list="restart_anam"
           print("read from defalt ", input_list)                   

if args.out != None:
           out_file=args.out
           print("output_file ",out_file)
else:
           out_file="results"
           print("default out_file ", out_file)


#define rotines
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
        f_head=False
        head=[]
        for line in read_pdb:
            idef = line[0:6].strip()
            if idef=="ATOM" or idef=="HETATM":
                f_head=True
                l_pdb=[]
                a       = line[0:6].strip()
                a_num   = int(line[6:11].strip())
                a_name  = line[12:16].strip()
                alter   = line[16].strip()
                res_name= line[17:20].strip()
                chain   = line[21].strip()
                res_num = int(line[22:26].strip())
                inser   = line[26].strip()
                x       = float(line[30:38].strip())
                y       = float(line[38:46].strip())
                z       = float(line[46:54].strip())
                occ     = float(line[54:60].strip())
                b_fac   = float(line[60:66].strip())
                seg_i   = line[72:76].strip()
                ele     = line[76:78].strip()
                charge  = line[78:80].strip()

                l_pdb.append(a)
                l_pdb.append(a_num)
                l_pdb.append(a_name)
                l_pdb.append(alter)
                l_pdb.append(res_name)
                l_pdb.append(chain)
                l_pdb.append(res_num)
                l_pdb.append(inser)
                l_pdb.append(x)
                l_pdb.append(y)
                l_pdb.append(z)
                l_pdb.append(occ)
                l_pdb.append(b_fac)
                l_pdb.append(seg_i)
                l_pdb.append(ele)
                l_pdb.append(charge)
#[0,   1  ,   2   ,  3   ,    4    ,   5  ,    6,  ,  7   , 8, 9,10
#[a, a_num, a_name, alter, res_name, chain, res_num, inser, x, y, z
#11 ,  12  ,  13  , 14 ,  15
#occ, b_fac, seg_i, ele, charge]
                org_pdb.append(l_pdb)
                work_pdb.append(l_pdb)
#               print(l_pdb)
            if f_head==False:
                head.append(line)
        print("close file:", input_pdb)
        return org_pdb, work_pdb, head

def read_qm(qm):
## read in atoms or resiodues for the qm system
## inly actueal atoms trancuation ateoms are definem on it own
## ATOM in PDB formart
## RESI only one number
## ANAM atom name with chain and residue
## ANUM atom number



        print(" read qm file")
        read_qm_f = open(qm, "r")
        qm_atom = []
        qm_resi = []
        qm_anam = []
        qm_anum = []
        qm_rnum = []
        for line in read_qm_f:
            idef = line[0:6].strip()
            idef2= line.split()



            if idef=="ANAM" or idef2[0]=="ANAM":
                s_line=line.split()
                l_resi=[]
                l_resi.append(s_line[1])
                l_resi.append(int(s_line[2]))
                if len(s_line)<4:
                        s_line.append("  ")
                l_resi.append((s_line[3]))
                qm_anam.append(l_resi)


        read_qm_f.close()
        print("close file:", qm)
        return  qm_anam

def extract_qm(pdb,qm):
# extract heavy atoms of syst1
        syst1=[]
        for i in range(len(pdb)):
          for j in range(len(qm)):
              if qm[j][1]==pdb[i][6]:
                 if qm[j][2]==pdb[i][5]:
                   if qm[j][0]==pdb[i][2]:
                      syst1.append(pdb[i])

        return syst1



def write_pdb(pdb,head, output):
## wite out PDB file
        res = open(output, "w")
#       res= open("comqum.pdb","w")
        for i in range(len(head)):
            line=head[i].split()
            if line[0]== "CRYST1":
               res.write(head[i])
            if line[0]== "SCALE1":
               res.write(head[i])
            if line[0]== "SCALE2":
               res.write(head[i])
            if line[0]== "SCALE2":
               res.write(head[i])

#           res.write("\n")
        for i in range(len(pdb)):
               string = str('{:6}'.format(pdb[i][0]))
               string = string + str('{:5.0f}'.format(pdb[i][1]))
               string = string + "  "
               string = string + str('{:3s}'.format(str(pdb[i][2])))
               string = string + str('{:1}'.format(pdb[i][3]))
               string = string + str('{:3}'.format(pdb[i][4]))
               string = string + str('{:>2}'.format(pdb[i][5]))
               string = string + str('{:4}'.format(pdb[i][6]))
               string = string + str('{:1}'.format(pdb[i][7]))
               string = string + "   "
               string = string + str('{:8.3f}'.format(pdb[i][8]))
               string = string + str('{:8.3f}'.format(pdb[i][9]))
               string = string + str('{:8.3f}'.format(pdb[i][10]))
               string = string + str('{:6.2f}'.format(pdb[i][11]))
               string = string + str('{:6.2f}'.format(pdb[i][12]))
               string = string + str('{:>7}'.format(pdb[i][13]))
               string = string + str('{:>5}'.format(pdb[i][14]))


               string = string + "\n"
               res.write(string)
        res.write("END")





qm_anam = read_qm(input_list)
org_pdb, work_pdb, head = read_pdb(input_pdb)
syst1 = extract_qm(work_pdb,qm_anam)
print(qm_anam)
for i in range(len(syst1)):
    print(syst1[i])
write_pdb(syst1,head, out_file)
