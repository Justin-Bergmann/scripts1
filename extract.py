#!/usr/bin/env python3
#####!/software/sse/easybuild/prefix/software/Anaconda3/5.3.0-extras-nsc1/bin-wrapped/python 
# -*- coding: utf-8 -*-


## program to automatise teh generation 
#of the syst1 file for quantum refinment


import time
## set timer
time_start = time.time()
import sys
##load mandatory modules
import numpy as np
#import numpy as np
print(sys.argv)
if len(sys.argv)>=2:
        qm=sys.argv[1]
        print("read from line 1")
else:
        qm="restart_anam"


if len(sys.argv)>=3:
        input_pdb=sys.argv[2]
        print("read from line 2")
else:
        input_eds="edstats.out"

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



def read_eds(input_eds):
## reads oxigen positions for water from a PDB file
#       print(" read D positions from pdb file")
        f_eds = open(input_eds, "r")
        l_eds=[]
        for line in f_eds:
            l_split = line.split()
            l_eds.append(l_split)
        f_eds.close()
        print("close file:", input_eds)
        return l_eds


def find_match(l_eds, qm):
## fine residues which where inlcudet in syst1
        match = []
        for i in range(len(l_eds)):
            for j in range(len(qm)):
#               print(len(l_eds[i]), len(qm[j]))
                if len(l_eds[i])>=2 and len(qm[j])==3 :
#                  print("TEST")
#                  print(l_eds[i][1],qm[j][2],l_eds[i][2],qm[j][1])
                   if l_eds[i][40].strip()== qm[j][2].strip() and  int(l_eds[i][2])== int(qm[j][1]):
                      match.append(l_eds[i])
#                     print(l_eds[i])
        
        for k in range(len(match)-2,-1,-1):
              if match[k][0]==match[k+1][0] and  match[k][1]==match[k+1][1] and match[k][2]==match[k+1][2]:
                del match[k+1]
        return(match)



def write_table(match):
## write results in table for latex formart
        tab=open("tabel_edstats","w")
        tab.write("protonation_state % & % RES % & %  N_grid_points % & %  RSR %  & %  RSCC % & %  RSZD % \\\\   \n")
        for i in range(len(match)):
            string=str(match[i][2])+"-"+ str(match[i][0])+"-"+str(match[i][1])
            string= string+"%&%"+match[i][4] +"%&%"+match[i][29] +"%&%"+ match[i][34] +"%&%"+ match[i][36]+ "% \\\\  \n" 
            tab.write(string)

        tab.write("\n\n\n ")
        tab.write("%&% res_num  ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][40]))
        tab.write("%\\\\ \n")
        tab.write("%&%res_nam  ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][0]))
        tab.write("%\\\\ \n")


        tab.write("%&%chain  ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][40]))
        tab.write("%\\\\ \n")

        tab.write("%&%N_points  ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][4]))
        tab.write("%\\\\ \n")
        tab.write("%&%RSR  ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][29]))
        tab.write("%\\\\ \n")
        tab.write("%&%RSCC  ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][34]))
        tab.write("%\\\\ \n")
        tab.write("%&%RSZD   ")
        for i in range(len(match)):
            tab.write("% & % ")
            tab.write(str(match[i][36]))
        tab.write("%\\\\\\hline \n")
        


def write_table2(match):
## write results in table for latex formart
        tab=open("tabel_edstats2","w")
        tab.write("protonation_state   &   RES   &    N_grid_points   &    RSR    &    RSCC   &    RSZD    \n")
        for i in range(len(match)):
            string=str(match[i][2])+"-"+  str(match[i][0])+"-"+str(match[i][40])
            string= string+" & "+match[i][4] +" & "+match[i][29] +" & "+ match[i][34] +" & "+ match[i][36]+ "   \n" 
            tab.write(string)

        tab.write("\n\n\n ")
        tab.write(" &  res_num  ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][2]))
        tab.write("  \n")
        tab.write(" & res_nam  ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][0]))
        tab.write("\n")


        tab.write(" & chain  ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][40]))
        tab.write(" \n")

        tab.write(" & N_points  ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][4]))
        tab.write("  \n")
        tab.write(" & RSR  ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][29]))
        tab.write("  \n")
        tab.write(" & RSCC  ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][34]))
        tab.write("  \n")
        tab.write(" & RSZD   ")
        for i in range(len(match)):
            tab.write("  &   ")
            tab.write(str(match[i][36]))
        tab.write("  \n")
        






#########################################################################################################
##processing is starting first read in all infomration from files


qm_anam = read_qm(qm)

l_eds = read_eds(input_eds)

#for i in range(len(qm_anam)):
#        print(qm_anam[i])
#        print(l_eds[i])
match = find_match(l_eds, qm_anam)
#print(qm_anam)
for i in range(len(match)):
        print(match[i])
print(len(match))

write_table(match)

write_table2(match)
time_ende = time.time()
print("program ends normally after "'{:5.3f}s'.format(time_ende-time_start)," or ", '{:5.2f}min'.format((time_ende-time_start)/60))














