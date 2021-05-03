#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: justin
"""

import time
## set timer
time_start = time.time()
##load mandatory modules
import numpy as np
import argparse, sys
#import scipy as sc
import gemmi

mesh = 0.10
label_H="D"
##input variables:
    
parser=argparse.ArgumentParser()
parser.add_argument('--pdb', help='input pdb file')
parser.add_argument('--map', help='input ccp4 map')
parser.add_argument('--out', help='output file')
parser.add_argument('--wat', help='protonate water, yes or no')
args=parser.parse_args()



## process input
if args.pdb != None:
   input_pdb=args.pdb
   print("input pdb file is:", input_pdb)
else:
         input_pdb="pdb"
         print("input pdb file is:", input_pdb)

if args.map != None:
           input_grid=args.map
           map_ex=True
           print("read from ", input_grid)
else:
           input_grid="None"
           map_ex= False
           print("read from defalt ", input_grid)

if args.out != None:
           out_file=args.out
           print("output_file ",out_file)
else:
           out_file="results"
           print("default out_file ", out_file)

if args.wat != None:
           wat=args.wat
           if wat=="yes":
                wat_ex=True
                print("water is getting protonatet")
           else:
              wat_ex=False
              print("water is not getting protonatet")

else:
          wat_ex=False
          print("water is not getting protonatet")

## reading and writeing rotines
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

        print("close file:", input_pdb)
        return org_pdb, work_pdb

def atm_to_res(pdb_read):
## change pdb from list of atoms to list of res with atoms
        pdb=[]
        pdb.append([pdb_read[0]])
        for i in range(1,len(pdb_read)):
           if pdb_read[i][6] != pdb_read[i-1][6]:
              pdb.append([pdb_read[i]])
           if pdb_read[i][6] == pdb_read[i-1][6]:
              pdb[len(pdb)-1].append(pdb_read[i])

        return pdb

def read_grid(input_grid):
#read ccp4 map with gemmi
        print(input_grid)

        map= gemmi.read_ccp4_map(input_grid)
        print(gemmi.InfoMap())
#       print("head",map.header_str)
        print(map.__repr__())
#       print(help(map))
        print(" map is readet")
        return map
    
    
    
def write_pdb(pdb, output):
#writes out final pdb file with D atoms
        res = open(output, "w")
        atm_nr=0
        for i in range(len(pdb)):
            for k in range(len(pdb[i])):
               atm_nr = atm_nr + 1
               string = str('{:6}'.format(pdb[i][k][0]))
               string = string + str('{:5.0f}'.format(atm_nr))
               string = string + " "
               string = string + str('{:4s}'.format(str(pdb[i][k][2])))
               string = string + str('{:1}'.format(pdb[i][k][3]))
               string = string + str('{:3}'.format(pdb[i][k][4]))
               string = string + str('{:>2}'.format(pdb[i][k][5]))
               string = string + str('{:4}'.format(pdb[i][k][6]))
               string = string + str('{:1}'.format(pdb[i][k][7]))
               string = string + "   "
               string = string + str('{:8.4f}'.format(pdb[i][k][8]))
               string = string + str('{:8.4f}'.format(pdb[i][k][9]))
               string = string + str('{:8.4f}'.format(pdb[i][k][10]))
               string = string + str('{:6.2f}'.format(pdb[i][k][11]))
               string = string + str('{:6.2f}'.format(pdb[i][k][12]))
               string = string + str('{:>10}'.format(pdb[i][k][13]))               
               string = string + str('{:>2}'.format(pdb[i][k][14]))


               string = string + "\n"
               res.write(string)
    
    
    

## work on pdb file

def treat_res(pdb,map):
#work on each residio seperatly 
        alerts=[["CH"],["NH"],["OH"],["SH"],["flip"],["wat"],["disordert"],[]]
#       alerts =[CH,NH,OH,SH,flip]     
#       pdb, alerts = treat_n_nerm(pdb,map,alerts)
        for res_num in range(len(pdb)):
#         print("res",pdb[res_num][0][6])
          if pdb[res_num][0][4]=="GLY":
                   pdb, alerts = treat_GLY( pdb, res_num,alerts )
          if pdb[res_num][0][4]=="PRO":
                   pdb, alerts = treat_PRO( pdb, res_num , alerts)
          if pdb[res_num][0][4]=="PHE":
                   pdb, alerts = treat_PHE( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="ALA":
                   pdb, alerts = treat_ALA( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="VAL":
                   pdb, alerts = treat_VAL( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="LEU":
                   pdb, alerts = treat_LEU( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="ILE":
                   pdb, alerts = treat_ILE( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="MET":
                   pdb, alerts = treat_MET( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="ARG":
                   pdb, alerts = treat_ARG( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="TRP":
                   pdb, alerts = treat_TRP( pdb, res_num, alerts )
          if pdb[res_num][0][4]=="ASN":
                   pdb, alerts = treat_ASN( pdb, res_num,map, alerts )
          if pdb[res_num][0][4]=="GLN":
                   pdb, alerts = treat_GLN(pdb,res_num,  map, alerts)
          if pdb[res_num][0][4]=="ASP":
                   pdb, alerts = treat_ASP( pdb, res_num,map, alerts )
          if pdb[res_num][0][4]=="GLU":
                   pdb, alerts = treat_GLU( pdb, res_num,map, alerts )
          if pdb[res_num][0][4]=="HIS" or pdb[res_num][0][4]=="HID" or   \
          pdb[res_num][0][4]=="HIE" or pdb[res_num][0][4]=="HIP" :
                   pdb, alerts = treat_HIS(pdb,res_num,  map, alerts)
          if pdb[res_num][0][4]=="LYS":
                   pdb, alerts = treat_LYS( pdb, res_num, map, alerts )
          if pdb[res_num][0][4]=="SER":
                   pdb, alerts = treat_SER( pdb, res_num ,map, alerts)
          if pdb[res_num][0][4]=="THR":
                   pdb, alerts = treat_THR( pdb, res_num ,map, alerts)
          if pdb[res_num][0][4]=="TYR":
                   pdb, alerts = treat_TYR( pdb, res_num ,map, alerts)
          if pdb[res_num][0][4]=="CYS":
                   pdb, alerts = treat_CYS( pdb, res_num, map, alerts )
          if pdb[res_num][0][4]=="WAT" or pdb[res_num][0][4]=="HOH":
                   pdb, alerts = treat_WAT( pdb, res_num,map, alerts  )
                   
        return pdb, alerts
    
######################################################
##All H well-defined
######################################################
def treat_GLY(pdb, res_num,alerts):
## add deteria to GLY
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
#       print(pdb[res_num])
        pdb = ad_CH2(pdb, res_num,alter,"N","C","CA")
        
        return pdb, alerts
        
        
def treat_PRO(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at PRO
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CD","CB","CG")
        pdb = ad_CH2(pdb, res_num,alter,"N","CG","CD")
        return pdb, alerts
        
        
def treat_PHE(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at PHE
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH_ar(pdb, res_num,alter,"CG","CE2","CD2")
        pdb = ad_CH_ar(pdb, res_num,alter,"CD2","CZ","CE2")
        pdb = ad_CH_ar(pdb, res_num,alter,"CG","CE1","CD1")
        pdb = ad_CH_ar(pdb, res_num,alter,"CD1","CZ","CE1")
        pdb = ad_CH_ar(pdb, res_num,alter,"CE2","CE1","CZ")

        return pdb, alerts

##############################################################
##CH3 in principle free, but uninteresting
#############################################################


def treat_ALA(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at ALA
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH3(pdb, res_num,alter,"C","CA","CB")
        return pdb, alerts


def treat_VAL(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at VAL
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH_R3(pdb, res_num,alter,"CG2","CG1","CB","CA")
        pdb = ad_CH3(pdb, res_num,alter,"CA","CB","CG1")
        pdb = ad_CH3(pdb, res_num,alter,"CA","CB","CG2")

        return pdb, alerts


def treat_LEU(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at LEU
        alter, alerts = find_alter(pdb,res_num,alerts)
#       print("ALTER",alter)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH_R3(pdb, res_num,alter,"CD2","CD1","CG","CB")
        pdb = ad_CH3(pdb, res_num,alter,"CB","CG","CD1")
        pdb = ad_CH3(pdb, res_num,alter,"CB","CG","CD2")

        return pdb, alerts

def treat_ILE(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at ILE
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH_R3(pdb, res_num,alter,"CG2","CG1","CB","CA")
        pdb = ad_CH3(pdb, res_num,alter,"CA","CB","CG2")
        pdb = ad_CH2(pdb, res_num,alter,"CD1","CB","CG1")
        pdb = ad_CH3(pdb, res_num,alter,"CB","CG1","CD1")
        return pdb, alerts

def treat_MET(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at MET
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH2(pdb, res_num,alter,"SD","CB","CG")
        pdb = ad_CH3(pdb, res_num,alter,"CG","SD","CE")

        return pdb, alerts

##############################################################
##All H well-defined
##############################################################

def treat_ARG(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at ARG
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CD","CB","CG")
        pdb = ad_CH2(pdb, res_num,alter,"CB","NE","CD")
        pdb = ad_NH_s(pdb, res_num,alter,"CZ","CD","NE")
        pdb = ad_NH2(pdb, res_num,alter,"NE","CZ","NH1")
        pdb = ad_NH2(pdb, res_num,alter,"NE","CZ","NH2")

        return pdb, alerts


##############################################################
## 1 well-defined H
##############################################################
def treat_TRP(pdb,res_num, alerts):
# ad hydrogen to standard psoitions at TRP
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH_ar(pdb, res_num,alter,"CG","NE1","CD1")
        pdb = ad_NH_s(pdb, res_num,alter,"CE2","CD1","NE1")
        pdb = ad_CH_ar(pdb, res_num,alter,"CZ3","CD2","CE3")
        pdb = ad_CH_ar(pdb, res_num,alter,"CE3","CH2","CZ3")
        pdb = ad_CH_ar(pdb, res_num,alter,"CZ3","CZ2","CH2")
        pdb = ad_CH_ar(pdb, res_num,alter,"CE2","CH2","CZ2")
        return pdb, alerts

##############################################################
## Flip possible
##############################################################


def treat_ASN(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at ASN and checks if O and N are
# fliped
        alter, alerts = find_alter(pdb,res_num,alerts)
#       print("Alter",alter)
        if map_ex==True:
           pdb,alerts = check_flip(pdb,res_num,alter,map,"ND2","OD1",alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num, alter,"CA","CG","CB")
        pdb = ad_NH2(pdb, res_num,alter,"CB","CG","ND2")
        return pdb, alerts


def treat_GLN(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at GLN
#        print("GLN!",pdb[res_num][0])
        alter, alerts = find_alter(pdb,res_num,alerts)

        if map_ex==True:
           pdb,alerts = check_flip(pdb,res_num,alter,map,"NE2","OE1",alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CD","CB","CG")
        pdb = ad_NH2(pdb, res_num,alter,"CG","CD","NE2")

        return pdb, alerts
#rotines for flip in ANS and GLN
def check_flip(pdb,res_num,alter,map,N,O,alerts):
#checks if HIS is fliped based on the map
#for this the map is integratet around CD2, CE1 and ND1, NE2. If the sum of the integrated values for carbon are smaller than for nitrogen the HIS is getting flipped.
        for alt in range(len(alter)):
         com=ceck_bond_atoms_2(pdb,res_num,alter[alt],N,O)
         if com==False:
            continue
         N_coord=anam_to_coord(pdb,res_num,alter[alt],N)
         O_coord=anam_to_coord(pdb,res_num,alter[alt],O)
         int_N=int_around_point(N_coord,map,def_r_int("N"),mesh)
         int_O=int_around_point(O_coord,map,def_r_int("N"),mesh)
#        print("sum_C", sum_C)        
#        print("sum_N", sum_N)
         if int_O > int_N:
#            print(alter[alt])
#            print(pdb[res_num][0])
#            print(pdb[res_num][1])
            print("N and C are  getting fliped for" ,pdb[res_num][0][5],pdb[res_num][0][4],pdb[res_num][0][6])
            pdb = del_NH2(pdb,res_num, alter[alt],N_coord)
            string = str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])+" got flipped"
            alerts[4].append(string)
            for i in range(len(pdb[res_num])):
              if pdb[res_num][i][3]==alter[alt]:
                if pdb[res_num][i][2]==N:
#                  print("¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤")
#                  print(pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10])
                   pdb[res_num][i][8]=O_coord[0]
                   pdb[res_num][i][9]=O_coord[1]
                   pdb[res_num][i][10]=O_coord[2]
#                  print(pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10])
#                  print("##################")
                if pdb[res_num][i][2]==O:
                   pdb[res_num][i][8]=N_coord[0]
                   pdb[res_num][i][9]=N_coord[1]
                   pdb[res_num][i][10]=N_coord[2]
#               print(pdb[res_num][i])



        return pdb, alerts
def ceck_bond_atoms_2(pdb,res_num,alt, A1, A2):
## cecks if all needet atoms to ad a D are existing in residue
        com = True
        res = []
        for i in range(len(pdb[res_num])):
#           print(len(list(pdb[res_num][i][3])))
            if pdb[res_num][i][3]==alt or\
            len(list(pdb[res_num][i][3]))==0:
#              print(pdb[res_num][i][3],len(list(pdb[res_num][i][3])))
               res.append(pdb[res_num][i][2])
            elif len(list(alt))==0:
               if pdb[res_num][i][3]=="A":
                  res.append(pdb[res_num][i][2])
#       print("A1, A2, A3",A1, A2, A3)
#       print("res",res)
        if A1 not in res or A2 not in res:
                com=False

        return com


def del_NH2(pdb,res_num,alt,A_coord):
##find CH bond du A3
#        H_coord=[]
        ele1="N"
        atom1=A_coord
        cov1=ele_to_cov_rad(ele1)
        for i in range(len(pdb[res_num])):
          if  pdb[res_num][i][3]==alt or len(list(pdb[res_num][i][3]))==0:
            atom2 =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
            ele2 = pdb[res_num][i][14]
            cov2=ele_to_cov_rad(ele2)
            d = CDist2(atom1,atom2)
            dmax= cov1+cov2+0.4
            if 0.5 <= d <= dmax:
               if pdb[res_num][i][14] == "H" or pdb[res_num][i][14]=="D":
                  del pdb[res_num][i]
        return pdb

##############################################################
## No on charged, 4 posible position on neutral
##############################################################

def treat_ASP(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at ASP and chek if it is
# protonatet or unprotonatet based on the mp
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        if map_ex==True:
           av_CH, std_CH= int_CH(pdb,res_num,alter,map,mesh)
           pdb,alerts = ad_CO2H_map(pdb, res_num,alter,map,"OD1", "CG","OD2",av_CH, std_CH, alerts)
        return pdb, alerts





def treat_GLU(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at GLU
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CD","CB","CG")
        if map_ex==True:
           av_CH, std_CH= int_CH(pdb,res_num,alter,map,mesh)
           pdb,alerts=ad_CO2H_map(pdb,res_num,alter,map,"OE1","CD","OE2",av_CH, std_CH,alerts)

        return pdb, alerts

def ad_CO2H_map(pdb, res_num,alter,map, A1, A2, A3,av_CH,std_CH,alerts):
# chekc if CO2- group is protonatet
     for alt in range(len(alter)):
        com =   check(pdb,res_num, alter[alt],A1, A2, A3, 2)
        if com == False:
                return pdb
        r  = 0.997
        Ad_l=[106.2,114.9]
        Dd_l=[0,180]
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        H_coords=[]
        for k in range(len(Dd_l)):
            H_coords.append([ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad_l[k],Dd_l[k]),A3])
            H_coords.append([ZtoXYZ(A3xyz,A2xyz,A1xyz,r,Ad_l[k],Dd_l[k]),A1])
        r_int=def_r_int("H")
        int_H_l=[]
        for k in range(len(H_coords)):
            int_H_l.append([H_coords[k][0],int_around_point(H_coords[k][0],map,r_int,mesh),H_coords[k][1]])
        test_H=[0,0,0]
#       print("int_H_l",int_H_l)

        for i in range(len(int_H_l)):
            if int_H_l[i][1]> test_H[1]:
                test_H[0]=int_H_l[i][0]
                test_H[1]=int_H_l[i][1]
                test_H[2]=int_H_l[i][2]
#       print("av_CH std_CH test_H[1]",av_CH, std_CH, test_H[1])
        if (av_CH -1*std_CH) <= test_H[1] and test_H[1]> 0:
#          print("test_H[2]",test_H[1],av_CH,std_CH)
           nam=list(test_H[2])
           if len(nam)<=2:
              name=nam[1]
           else:
             name = nam[1]#+nam[2]

           name1 = label_H+name+"2"

#          name1 ="H"+test_H[2]
#           print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!",nam[2],"######################################")
           pdb = add_to_pdb(pdb,res_num,alter[alt], test_H[2],   test_H[0], name1,test_H[1])
#          
           if nam[2]=="1":
              for i in range(len(pdb[res_num])):
                  if pdb[res_num][i][2]==A1 and pdb[res_num][i][3]==alter[alt]:
                     O1=i
                  if pdb[res_num][i][2]==A3 and pdb[res_num][i][3]==alter[alt]:
                     O2=i
              pdb[res_num][O1][2]=A3
              pdb[res_num][O2][2]=A1
           print(name1 ,"is addet to ",test_H[2],res_num)
           string=str(name1)+" is addet to "+str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])
           alerts[2].append(string)
#       else:
#          string=str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])+" is charged"
#          alerts[2].append(string)
     return pdb, alerts

##############################################################
## Hid, Hie, Hip and flip
##############################################################
def treat_HIS(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at HIS
        alter, alerts = find_alter(pdb,res_num,alerts)
        if map_ex==True:
           pdb,alerts =check_flip_HIS(pdb,res_num,alter,map,alerts)
        proto = pdb[res_num][0][4]
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH_ar(pdb, res_num,alter,"ND1","NE2","CE1")
        pdb = ad_CH_ar(pdb, res_num,alter,"CG","NE2","CD2")
        if map_ex==True:
           av_CH,sig_CH= int_CH(pdb,res_num,alter,map,mesh)
           pdb,alerts = ad_NH_s_map(pdb, res_num,alter,map,"CE1","CD2","NE2",av_CH,sig_CH,alerts)
           pdb,alerts = ad_NH_s_map(pdb,           res_num,alter,map,"CE1","CG","ND1",av_CH,sig_CH,alerts)
        else:
           if proto == "HIS" or\
           proto == "HIE" or proto == "HIP":
                   pdb = ad_NH_s(pdb, res_num,alter,"CE1","CD2","NE2")
           if proto == "HID" or  proto == "HIP":
                   pdb = ad_NH_s(pdb, res_num,alter,"CE1","CG","ND1")

        return pdb, alerts

##############################################################
## 3 H on a ring; may be deprotonated
##############################################################
def treat_LYS(pdb,res_num, map, alerts):
# ad hydrogen to standard psoitions at LYS
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","CG","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CD","CB","CG")
        pdb = ad_CH2(pdb, res_num,alter,"CE","CG","CD")
        pdb = ad_CH2(pdb, res_num,alter,"NZ","CD","CE")
        if map_ex==True:
           av_CH,sig_CH= int_CH(pdb,res_num,alter,map,mesh)
           pdb = ad_NH3_map(pdb,res_num,alter,"CD","CE","NZ",map,av_CH,sig_CH)
        else:
           pdb = ad_NH3(pdb, res_num,alter,"CD","CE","NZ")
        return pdb, alerts

def treat_n_nerm(pdb,map,alerts):
## protonates NH3 as Nterminal
        res_num=0
        alter, alerts = find_alter(pdb,res_num,alerts)
        if map_ex==True:
           pdb = ad_NH3_map_ter(pdb,res_num,alter,"C","CA","N",map)
        else:
           pdb = ad_NH3(pdb, res_num,alter,"C","CA","N")



        return pdb, alerts
##############################################################
## 1 H on a ring
##############################################################
def treat_SER(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at SER
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","OG","CB")
        if map_ex==True:
           pdb, alerts = ad_OH_map(pdb, res_num,alter,"CA","CB","OG",map, alerts)
        else:
           pdb = ad_OH(pdb, res_num,alter,"CA","CB","OG")
        return pdb, alerts

def treat_THR(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at THR
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH_R3(pdb, res_num,alter,"CG2","OG1","CB","CA")
        pdb = ad_CH3(pdb, res_num,alter,"CA","CB","CG2")
        if map_ex==True:
           pdb, alerts = ad_OH_map(pdb, res_num,alter,"CA","CB","OG1",map, alerts)
        else:
           pdb = ad_OH(pdb, res_num,alter,"CA","CB","OG1")

        return pdb, alerts

def treat_TYR(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at TYR
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num, alter,"CA","CG","CB")
        pdb = ad_CH_ar(pdb, res_num,alter,"CE2","CG","CD2")
        pdb = ad_CH_ar(pdb, res_num,alter,"CZ","CD2","CE2")
        pdb = ad_CH_ar(pdb, res_num,alter,"CD1","CZ","CE1")
        pdb = ad_CH_ar(pdb, res_num,alter,"CG","CE1","CD1")
        if map_ex==True:
           pdb, alerts = ad_OH_map(pdb, res_num,alter,"CE1","CZ","OH",map, alerts)
        else:
           pdb = ad_OH(pdb, res_num,alter,"CE1","CZ","OH")

        return pdb, alerts
##############################################################
## 1 H on a ring; may be deprotonated
##############################################################
def treat_CYS(pdb,res_num, map, alerts):
# ad hydrogen to standard psoitions at CYS
        alter, alerts = find_alter(pdb,res_num,alerts)
        pdb = ad_NH(pdb,res_num,alter)
        pdb = ad_CH_R3(pdb, res_num,alter,"C","N","CA","CB")
        pdb = ad_CH2(pdb, res_num,alter,"CA","SG","CB")

        if map_ex==True:
           pdb, alerts = ad_SH_map(pdb, res_num,alter,"CA","CB","SG",map,alerts)
        else:
           pdb = ad_SH(pdb, res_num,alter,"CA","CB","SG")
        return pdb, alerts

##############################################################
## 2 H on a sphere
##############################################################
def treat_WAT(pdb,res_num,map, alerts):
# ad hydrogen to standard psoitions at WAT
        if map_ex==True and wat_ex==True:
           alter, alerts = find_alter(pdb,res_num,alerts)
           print("alter",alter)
           sphere=def_sphere(0.98, 20)
           r_int=def_r_int("H")
           coord_O=[pdb[res_num][0][8],pdb[res_num][0][9],pdb[res_num][0][10]]
#          print(pdb[res_num])
#          print(coord_O)
           int_sphere_l=int_sphere(sphere,map,coord_O,r_int)
           best_fit=find_best_fit(int_sphere_l,coord_O)
           if best_fit[2]<= best_fit[0][2]:
              print("integratet value is to small hydrogens are not aded")
              return pdb , alerts
           print("best_fit","sum",best_fit[2],"H1",best_fit[0][1],"H2",best_fit[1][1],"O",best_fit[0][2])
           pdb = add_to_pdb(pdb,res_num,alter[0], "O", best_fit[0][0], "H2",best_fit[1][1])
           pdb = add_to_pdb(pdb,res_num,alter[0], "O", best_fit[1][0], "H1",best_fit[0][1])


        return pdb, alerts
def find_best_fit(l_int, coord_O):
# finds atom pair with best fit in 104,5 +- 5 degree
        pairs=[]
#        bes_fit=[]
        g_opt=104.5
        g_sig=10
        int_sum=-1000
        for i in range(len(l_int)):
          v1=twoP_to_vec(coord_O,l_int[i][0])
          for k in range(i,len(l_int)):
            v2=twoP_to_vec(coord_O,l_int[k][0])
            ang=C2Angle(v1,v2)
            if g_opt-g_sig <= ang <= g_opt+g_sig:
               pairs.append([l_int[i],l_int[k],l_int[i][1]+l_int[k][1]])
        print("len(pairs)",len(pairs))
        for i in range(len(pairs)):
            if int_sum < pairs[i][2]:
               int_sum=pairs[i][2]
               best_fit=pairs[i]
        return best_fit


def int_sphere(sphere,map,coord_O,r_int):
#integrates for all points in list sphere on map 
        int_sphere_l=[]
#       r_int_O=ele_to_cov_rad("O")
        int_O=int_around_point(coord_O,map,r_int, mesh)

        for i in range(len(sphere)):
            H=[sphere[i][0]+coord_O[0],sphere[i][1]+coord_O[1],sphere[i][2]+coord_O[2]]
            int_H=int_around_point(H,map,r_int, mesh)
            int_sphere_l.append([H,int_H,int_O])
        return int_sphere_l


def def_sphere(r,deg):
#defines sphere wit radius r around 0,0,0
        sphere=[]
        for teth in range(0,180,deg):
            for phi in range(0,360,deg):
                x=r*np.sin(np.radians(teth))*np.cos(np.radians(phi))
                y=r*np.sin(np.radians(teth))*np.sin(np.radians(phi))
                z=r*np.cos(np.radians(teth))
                sphere.append([x,y,z])
        print("len(sphere)",len(sphere))


        return sphere

## rotiens to work on residues
def find_alter(pdb,res_num,alerts):
## finds alternative conformations in residue
        alter=[]
        for i in range(len(pdb[res_num])):
            alter.append(pdb[res_num][i][3])
        alter=list(dict.fromkeys(alter))
        if len(alter)>=2:
           string=str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])+" is disordert"
           alerts[6].append(string)
#       print("alter",alter)

        return alter, alerts

def ad_NH(pdb,res_num,alter):
## ad HN if a naiboring resiue is existing
#     print("ad_NH_print alter",alter,pdb[res_num][0][6])
     for alt in range(len(alter)):
#          res=pdb[res_num][0][6]
##      print("RES",res)
#        for i in range(len(pdb)):
#          for k in range(len(pdb[res_num])):
#            if pdb[i][k][6]== res-1 and pdb[i][k][2].strip()=="C" and  pdb[res_num][0][5]==pdb[i][k][5]:
                com = ceck_bond_atoms_N(pdb,res_num,res_num,alter[alt],"N","CA", "C")
#               if len(alter)>1:
#                  print("alter[alt],com",alter[alt],com)
                if com == False:
                    print("missing atoms to locate NH")
                    return pdb
#               print(pdb[res_num])
#               print("alter,res:num", alter[alt],res_num)
#                nr_b=nr_bonds(pdb,res_num,alter[alt],"N")
#                ext_H=False
#                print("alter",alter)
                ext_H=check_extH_N(pdb,res_num,"N",alter[alt])
                if ext_H == True:
                    print("NH not, H is already existing!!!!!!!!!!!!!!!!")
#                    print("at residue",pdb[res_num][0][6])
#                    if alter[alt]=="B":
#                        print("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
#                        print("Res_num",pdb[i][k][6])
                    
                    return pdb

                int_H=0
                Nxyz  = anam_to_coord(pdb,res_num,alter[alt],"N")
#               print(pdb[res_num])
                CAxyz  = anam_to_coord(pdb,res_num,alter[alt],"CA")
                Cxyz  = anam_to_coord(pdb,res_num-1,alter[alt],"C")
                r=1.016
                Ad=117.7
                Dd=180
#                if alter[alt]=="B":
 #                       print("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
 #                       print("Res_num",pdb[i][k][6])
#               if alter[alt]=="A":
#                       print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
#                       print("res:num",pdb[i][k][6])
                H = ZtoXYZ(Cxyz,CAxyz,Nxyz,r,Ad,Dd)
                if map_ex==True:
                   int_H=int_around_point(H,map,def_r_int("H"), mesh)
                pdb =  add_to_pdb(pdb,res_num,alter[alt],"N",H,label_H,int_H)
     return pdb
 
def check_extH_N(pdb,res_num,atom,alt):
# cheks if H of NH is existing
    ext=False
    atom1, ele1 =anam_to_coord_ele(pdb,res_num,alt,atom)
    cov1=ele_to_cov_rad(ele1)
#    print("atom1, ele1",atom1, ele1,alt)
#    print("find NH ex", pdb[res_num][0])
    for i in range(len(pdb[res_num])):
          if  pdb[res_num][i][3]==alt:# or len(list(pdb[res_num][i][3]))==0:
#           print(alt)
#           print("test alt",alt)
            if pdb[res_num][i][14]=="H" or pdb[res_num][i][14]=="D":
                atom2 =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
                ele2 = pdb[res_num][i][14]
                cov2=ele_to_cov_rad(ele2)
                d = CDist2(atom1,atom2)
                dmax= cov1+cov2+0.4
                if 0.5 <= d <= dmax:
                    ext=True
#                    print("ATOM2", pdb[res_num][i])

    return ext





def ad_CA(pdb,res_num,alter):
#add HCA at standard positons
     for alt in range(len(alter)):
        com= check(pdb,res_num,alter[alt],"C","N","CA", 3)
        print("ALTER",alter)
        print("len(alter",len(alter))
        if len(alter)>1:
           print("alter,com,QQQQQQQQQQQQ", alter,com)

        if com == False:
                return pdb
        if len(alter)>1:
         print("QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ")
        CAxyz  = anam_to_coord(pdb,res_num,alter[alt],"CA")
        Nxyz  = anam_to_coord(pdb,res_num,alter[alt],"N")
        Cxyz  = anam_to_coord(pdb,res_num,alter[alt],"C")

        r = 1.095
        Ad = 109
        Dd= -121.5
        HCA = ZtoXYZ(Nxyz,Cxyz,CAxyz,r,Ad,Dd)
        int_H=0
        if map_ex==True:
           int_H=int_around_point(HCA,map,def_r_int("H"), mesh)
        pdb = add_to_pdb(pdb,res_num,alter[alt],"CA",HCA,"HA",int_H)
     return pdb

def ad_CH3(pdb, res_num,alter, A1, A2, A3):
# ad 3H as CH3 group to A3
     for alt in range(len(alter)):
#       print("alter[alt]",alter[alt])
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
#       print(com)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r   = 1.095
        Ad  = 109
        Dd  = 180
        Dd1 = 60
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd1)
        adH3 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,-Dd1)
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name+"1"
        name2 = label_H+name+"2"
        name3 = label_H+name+"3"
        int_H1=int_H2=int_H3=0
        if map_ex==True:
           int_H1=int_around_point(adH1,map,def_r_int("H"), mesh)
           int_H2=int_around_point(adH2,map,def_r_int("H"), mesh)
           int_H3=int_around_point(adH3,map,def_r_int("H"), mesh)
#       print("adH1",adH1)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H1)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH2, name2,int_H2)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH3, name3,int_H3)
     return pdb

def ad_NH3(pdb, res_num,alter, A1, A2, A3):
# ad 3H as NH3+ group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb

#        print( pdb[res_num])
#        print("A2", A2)
        r   = 1.019
        Ad  = 109
        Dd  = 180
        Dd1 = 60
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd1)
        adH3 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,-Dd1)
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name+"1"
        name2 = label_H+name+"2"
        name3 = label_H+name+"3"
        int_H1=int_H2=int_H3=0
        
        if map_ex==True:
           int_H1=int_around_point(adH1,map,def_r_int("H"), mesh)
           int_H2=int_around_point(adH2,map,def_r_int("H"), mesh)
           int_H3=int_around_point(adH3,map,def_r_int("H"), mesh)

        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H1)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH2, name2,int_H2)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH3, name3,int_H3)
     return pdb



def ad_NH3_map(pdb, res_num,alter, A1, A2, A3,map,av_CH,sig_CH):
# ad 3H as NH3+ group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb

#        print( pdb[res_num])
#        print("A2", A2)
        r   = 1.019
        Ad  = 109
        Dd  = 180
        r_int=def_r_int("H")
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        int_l=[]
        int_H1=int_H2=int_H3=0
        for i in range(0,120,5):
           adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd+i)
           adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd-120+i)
           adH3 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd-(2*120)+i)
           int_H1=int_around_point(adH1,map,r_int, mesh)
           int_H2=int_around_point(adH2,map,r_int, mesh)
           int_H3=int_around_point(adH3,map,r_int, mesh)
           int_H12=int_H1+int_H2
           int_H13=int_H1+int_H3
           int_H23=int_H2+int_H3
           int_H123=int_H1+int_H2+int_H3           
           int_l.append([[adH1,int_H1],[adH2,int_H2],[adH3,int_H3],int_H123,int_H12,int_H13,int_H23,int_H23])
        test_int=-1000
#        print("av_CH,sig_CH",av_CH,sig_CH)
        for i in range(len(int_l)):
#           if (av_CH*3 -6*sig_CH) <= int_l[i][3] :
               if test_int < int_l[i][3]:
                  test_int=int_l[i][3]
                  adH1=int_l[i][0][0]
                  adH2=int_l[i][1][0]
                  adH3=int_l[i][2][0]
                  int_H1=int_l[i][0][1]
                  int_H2=int_l[i][1][1]
                  int_H3=int_l[i][2][1]
#        print("av_CH*3-6*sig_CH", av_CH*3-6*sig_CH)         
#        print("test_int",test_int/3)
#        print("######################")
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name+"1"
        name2 = label_H+name+"2"
        name3 = label_H+name+"3"
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H1)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH2, name2,int_H2)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH3, name3,int_H3)
     return pdb


def ad_NH3_map_ter(pdb, res_num,alter, A1, A2, A3,map):
# ad 3H as NH3+ group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb

#        print( pdb[res_num])
#        print("A2", A2)
        r   = 1.019
        Ad  = 109
        Dd  = 180
        r_int=def_r_int("H")
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        int_l=[]
        int_H1=int_H2=int_H3=0
        for i in range(0,120,5):
           adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd+i)
           adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd-120+i)
           adH3 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd-(2*120)+i)
#          print("ATOMS",i,adH1,adH2,adH3)

           int_H1=int_around_point(adH1,map,r_int, mesh)
           int_H2=int_around_point(adH2,map,r_int, mesh)
           int_H3=int_around_point(adH3,map,r_int, mesh)
           int_H12=int_H1+int_H2
           int_H13=int_H1+int_H3
           int_H23=int_H2+int_H3
           int_H123=int_H1+int_H2+int_H3
           int_l.append([[adH1,int_H1],[adH2,int_H2],[adH3,int_H3],int_H123,int_H12,int_H13,int_H23,int_H23])
        test_int=-1000
        for i in range(len(int_l)):
#           if (av_CH*3 -6*sig_CH) <= int_l[i][3] :
               if test_int < int_l[i][3]:
                  test_int=int_l[i][3]
                  adH1=int_l[i][0][0]
                  adH2=int_l[i][1][0]
                  adH3=int_l[i][2][0]
                  int_H1=int_l[i][0][1]
                  int_H2=int_l[i][1][1]
                  int_H3=int_l[i][2][1]
        nam=list(A3)
        if len(nam)<=2:
           name=nam[0]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name+"1"
        name2 = label_H+name+"2"
        name3 = label_H+name+"3"
#       print("TER",name1,adH1)
#       print("TER",name2,adH2)
#       print("TER",name3,adH3)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH3, name3,int_H3)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH2, name2,int_H2)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H1)
     return pdb

def ad_CH2(pdb, res_num, alter, A1, A2, A3):
# ad 2H as CH2 group to A3
        for alt in range(len(alter)):
           
           com =   check(pdb,res_num, alter[alt],A1, A2, A3, 2)
           if com == False:
                   return pdb
          
#          print( pdb[res_num])
#          print("A2", A2)
           r = 1.095
           Ad = 109
           Dd= -121.5
           A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
           A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
           A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
           adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
           adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,-Dd)
           nam=list(A3)
           if len(nam)<=2:
              name=nam[1]
           else:
              name = nam[1]+nam[2]
           int_H1=int_H2=0 
           if map_ex==True:
              int_H1=int_around_point(adH1,map,def_r_int("H"), mesh)
              int_H2=int_around_point(adH2,map,def_r_int("H"), mesh)           

           name1 = label_H+name+"1"
           name2 = label_H+name+"2"
           pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H1)
           pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH2, name2,int_H2)

        return pdb

def ad_NH2(pdb, res_num,alter, A1, A2, A3):
# ad 2H as NH2 group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb

#       print( pdb[res_num])
#       print("A2", A2)
        r   = 1.012
        Ad  = 120
        Dd  = 180
        Dd1 = 0
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd1)
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]
        int_H1=int_H2=0
        if map_ex==True:
           int_H1=int_around_point(adH1,map,def_r_int("H"), mesh)
           int_H2=int_around_point(adH2,map,def_r_int("H"), mesh)
        name1 = label_H+name+"1"
        name2 = label_H+name+"2"
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H1)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH2, name2,int_H2)

     return pdb



def ad_CH_ar(pdb, res_num,alter, A1, A2, A3):
# ad H as CH (SP2) group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num, alter[alt],A1, A2, A3, 2)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 1.014
        Ad = 120
        Dd = 180
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
#       name1 =label_H+A3
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name
        int_H=0
        if map_ex==True:
           int_H=int_around_point(adH1,map,def_r_int("H"), mesh)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H)

     return pdb

def ad_NH_s(pdb, res_num,alter, A1, A2, A3):
# ad H as NH (SP2) group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num, alter[alt],A1, A2, A3, 2)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 1.032
        Ad = 120
        Dd = 180
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
#       name1 ="H"+A3
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name
        int_H=0
        if map_ex==True:
           int_H=int_around_point(adH1,map,def_r_int("H"), mesh)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H)

     return pdb



def ad_NH_s_map(pdb, res_num,alter,map, A1, A2, A3,av_CH,std_CH,alerts):
# ad H as NH (SP2) group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num, alter[alt],A1, A2, A3, 2)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 1.032
        Ad = 120
        Dd = 180
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
#       print("a1,a2,a3", A1,A2,A3)
#       print(A1xyz,A2xyz,A3xyz)

        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        r_int=def_r_int("H")
#       print("adH1", adH1, r_int)
        int_H=int_around_point(adH1,map,r_int, mesh)        
#       print("int_NH",int_H,av_CH, std_CH)
        if (av_CH -2*std_CH) <= int_H :
           nam=list(A3)
           if len(nam)<=2:
              name=nam[1]
           else:
              name = nam[1]+nam[2]
           
           name1 = label_H+name

           pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H)
           string= A3+" in HIS "+str(pdb[res_num][0][6])+" is protonatet"
           alerts[1].append(string)
        else:
           string= A3+" in HIS "+str(res_num)+" is deprotonatet"
           alerts[1].append(string)

     return pdb,alerts


def check_flip_HIS(pdb,res_num,alter,map,alerts):
#checks if HIS is fliped based on the map
#for this the map is integratet around CD2, CE1 and ND1, NE2. If the sum of the integrated values for carbon are smaller than for nitrogen the HIS is getting flipped.
        for alt in range(len(alter)):
         
         CD2=anam_to_coord(pdb,res_num,alter[alt],"CD2")
         CE1=anam_to_coord(pdb,res_num,alter[alt],"CE1")
         ND1=anam_to_coord(pdb,res_num,alter[alt],"ND1")
         NE2=anam_to_coord(pdb,res_num,alter[alt],"NE2")
         int_CD2=int_around_point(CD2,map,def_r_int("N"),mesh)
         int_CE1=int_around_point(CE1,map,def_r_int("N"),mesh)
         int_ND1=int_around_point(ND1,map,def_r_int("N"),mesh)
         int_NE2=int_around_point(NE2,map,def_r_int("N"),mesh)
         sum_C=int_CD2 + int_CE1
         sum_N=int_ND1 + int_NE2
#        print("sum_C", sum_C)        
#        print("sum_N", sum_N)
         if sum_C > sum_N:
            string = str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])+" got flipped"
            alerts[4].append(string)
            print("HIS is getting fliped",pdb[res_num][0][5],pdb[res_num][0][4],pdb[res_num][0][6])
            for i in range(len(pdb[res_num])):
                if pdb[res_num][i][2]=="CD2":
                   print("¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤")
                   print(pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10])
                   pdb[res_num][i][8]=ND1[0]
                   pdb[res_num][i][9]=ND1[1]
                   pdb[res_num][i][10]=ND1[2]
                   print(pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10])
                   print("##################")
                if pdb[res_num][i][2]=="CE1":
                   pdb[res_num][i][8]=NE2[0]
                   pdb[res_num][i][9]=NE2[1]
                   pdb[res_num][i][10]=NE2[2]
                if pdb[res_num][i][2]=="ND1":
                   pdb[res_num][i][8]=CD2[0]
                   pdb[res_num][i][9]=CD2[1]
                   pdb[res_num][i][10]=CD2[2]
                if pdb[res_num][i][2]=="NE2":
                   pdb[res_num][i][8]=CE1[0]
                   pdb[res_num][i][9]=CE1[1]
                   pdb[res_num][i][10]=CE1[2]
                print(pdb[res_num][i])



        return pdb, alerts


def ad_OH(pdb, res_num,alter, A1, A2, A3):
# ad H as OH group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 0.978
        Ad = 106
        Dd = 180
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name

        int_H=0
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H)

     return pdb



def ad_OH_map(pdb, res_num,alter, A1, A2, A3,map, alerts):
# ad H as OH group to A3
#     print("res name and number",pdb[res_num][0][4],pdb[res_num][0][6])
     av_CH,sig_CH= int_CH(pdb,res_num,alter,map,mesh)
     r_int=def_r_int("H")
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 0.978
        Ad = 106
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        int_l=[]
        for Dd in range(0,360,5):
            adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
            int_H1=int_around_point(adH1,map,r_int, mesh)
            int_l.append([adH1,int_H1])

        test_int=-1000
        for i in range(len(int_l)):
#              print("test_int, int_l",test_int,int_l[i][1])
               if test_int < int_l[i][1]:
                  test_int=int_l[i][1]
                  adH1=int_l[i][0]
#        print("av_CH-sig_CH", av_CH,sig_CH)
#        print("test_int",test_int)
#        print("#########OH###########")
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]
        name1 = label_H+name
        if av_CH-2*sig_CH> test_int:
           string=str(name1)+" of "+str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])+" might be deprotonated"
           alerts[2].append(string)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1,  name1,test_int)

     return pdb, alerts




def ad_SH(pdb, res_num, alter, A1, A2, A3):
# ad H as SH group to A3
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 1.343
        Ad = 96
        Dd = 180
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name

        int_H=0
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H)

     return pdb

def ad_SH_map(pdb, res_num, alter, A1, A2, A3,map,alerts):
# ad H as SH group to A3 based on the map
#     print("res name and number",pdb[res_num][0][4],pdb[res_num][0][6])
     av_CH,sig_CH= int_CH(pdb,res_num,alter,map,mesh)
     r_int=def_r_int("H")
     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
        if com == False:
                return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        r  = 1.343
        Ad = 96
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        int_l=[]
        for Dd in range(0,360,5):
            adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
            int_H1=int_around_point(adH1,map,r_int, mesh)
            int_l.append([adH1,int_H1])


        test_int=-1000
        for i in range(len(int_l)):
#              print("test_int, int_l",test_int,int_l[i][1])
               if test_int < int_l[i][1]:
                  test_int=int_l[i][1]
                  adH1=int_l[i][0]
#        print("av_CH-sig_CH", av_CH,sig_CH)
#        print("test_int",test_int)
#        print("¤¤¤¤¤¤¤¤¤¤¤¤SH¤¤¤¤¤¤¤¤¤¤¤¤¤¤")
        
        if test_int <  av_CH-2*sig_CH or test_int<=0:
           string=str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])+" is deprotonated"
           alerts[3].append(string)
           print("retrunrn unprotonatet S")
           return pdb,alerts
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 =label_H+name
#        print("test_int",test_int)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1,    name1,test_int)
        string=str(name1)+" is addet to "+str(pdb[res_num][0][4])+" "+str(pdb[res_num][0][6])+" "+str(pdb[res_num][0][5])
        alerts[3].append(string)

     return pdb,alerts


def ad_CH_R3(pdb, res_num,alter, A1, A2, A3, A4):
# ad H as R3CH (SP3) group to A3

     for alt in range(len(alter)):
        com =   check(pdb,res_num,alter[alt],A1, A2, A3, 3)
        if com == False:
                return pdb
        com =   check(pdb,res_num,alter[alt],A4, A2, A3, 3)
        if com == False:
                return pdb
        r  = 1.032
#        Ad = 109
#        Dd = 116
        A1xyz  = anam_to_coord(pdb,res_num,alter[alt],A1)
        A2xyz  = anam_to_coord(pdb,res_num,alter[alt],A2)
        A3xyz  = anam_to_coord(pdb,res_num,alter[alt],A3)
        A4xyz  = anam_to_coord(pdb,res_num,alter[alt],A4)
        adH1 = H_tert(A1xyz, A2xyz, A3xyz, A4xyz,r)
        nam=list(A3)
        if len(nam)<=2:
           name=nam[1]
        else:
           name = nam[1]+nam[2]

        name1 = label_H+name

        int_H=0
        if map_ex==True:
           int_H=int_around_point(adH1,map,def_r_int("H"), mesh)
        pdb = add_to_pdb(pdb,res_num,alter[alt], A3, adH1, name1,int_H)

     return pdb
 
def H_tert(A1, A2, A3, A4,r):
## set terzier H to A3 in dinstace r
        v1=twoP_to_vec(A3,A1)
        v2=twoP_to_vec(A3,A2)
        v3=twoP_to_vec(A3,A4)
        v_sum=[0,0,0]
        for i in range(len(v_sum)):
            v_sum[i]=-(v1[i]+v2[i]+v3[i])

        len_v_sum=len3dvec(v_sum)
        H_vec=[0,0,0]
        for i in range(len(v_sum)):
            H_vec[i]=(r*v_sum[i])/len_v_sum
        H_pos=[0,0,0]
        for i in range(len(H_vec)):
                  H_pos[i]=A3[i]+H_vec[i]
        return H_pos
    
## H functions
def cov_rad():
#list of Covalent radii in  from  Cambridge Structural Database
        cov_rad=[["H",0.31],["D",0.31],["C",0.76],["N",0.71],["O",0.66],["S",1.05],["Fe",1.52],["dum",0.0]]
        return cov_rad
    
def anam_to_coord(pdb,res_num,alt,anam):
#strip coordinates from atom name and res number
        for i in range(len(pdb[res_num])):
           if pdb[res_num][i][2].strip()==anam:
            if pdb[res_num][i][3].strip()==alt \
            or len(list(pdb[res_num][i][3]))==0:
              xyz =[0,0,0]
              xyz =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
            elif len(list(alt))==0:
               if pdb[res_num][i][3]=="A":
                  xyz =[0,0,0]
                  xyz =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
        return xyz

def add_to_pdb(pdb,res_num,alt,old_atm,new_koord,new_nam,int_H):
## ad a new deuterium atom to behind an existing atom in the pdb list
        atom=0
        for i in range(len(pdb[res_num])):
#          print("ALT",alt)
#          atom=0
#          print(pdb[res_num][i][3])
#          print(pdb[res_num][i][2].strip())
           if pdb[res_num][i][2].strip()==old_atm and\
           pdb[res_num][i][3]==alt:
              atom=pdb[res_num][i]
              indice = i+1
           
        if atom == 0:
           return pdb
        l_pdb=[]       
        l_pdb.append(atom[0])
        l_pdb.append(atom[1])
        l_pdb.append(new_nam)
        l_pdb.append(atom[3])
        l_pdb.append(atom[4])
        l_pdb.append(atom[5])
        l_pdb.append(atom[6])
        l_pdb.append(atom[7])
        l_pdb.append(new_koord[0])
        l_pdb.append(new_koord[1])
        l_pdb.append(new_koord[2])
        l_pdb.append(atom[11])
        l_pdb.append(atom[12])
        l_pdb.append(atom[13])
        l_pdb.append("D")
        l_pdb.append(" ")
        if map_ex==True:
           l_pdb.append(int_H)
        pdb[res_num].insert(indice, l_pdb)        
#       pdb[res_num].append(l_pdb)
        return pdb
    

def CDist2(A,B):
#calculate distance betweenn two points
        dist = len3dvec(twoP_to_vec(A, B))
        return dist



def len3dvec(vec):
## calculates lengh of a 3D vecor
## input as list
        a = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        return a

def twoP_to_vec(A,B):
#creates vector between two points
        vec = np.array([B[0]-A[0], B[1]-A[1], B[2]-A[2]])

        return vec

def int_around_point(point, map, r, mesh):
#integrates in certan radius around point
        int_val=0
        l_val=[]
        vol=0
        ran_int=range_to_int(point,r,mesh)
        for i in range(ran_int[0],ran_int[1]):
          for j in range(ran_int[2],ran_int[3]):
            for k in range(ran_int[4],ran_int[5]):
                p=[i*mesh, j*mesh,k*mesh]
                dist=CDist2(point,p)
#               print(dist)
#               print(point)
#               print(p) 
                if dist <= r:
                  val=map.grid.interpolate_value(gemmi.Position(p[0],p[1],p[2]))
                  val=val*mesh**3
                  vol=vol+mesh**3
                  l_val.append(val)
        int_val=sum(l_val)#/vol
        return int_val

def range_to_int(point,r, mesh):
#set up range to integrate [x1,x2,y1,y2,z1,z2]
        range_int=[0,0,0,0,0,0]
        range_int[0]=int((point[0]-r)/mesh-1)
        range_int[1]=int((point[0]+r)/mesh+1)
        range_int[2]=int((point[1]-r)/mesh-1)
        range_int[3]=int((point[1]+r)/mesh+1)
        range_int[4]=int((point[2]-r)/mesh-1)
        range_int[5]=int((point[2]+r)/mesh+1)

        return range_int


def def_r_int(ele):
## creates integration radios out of the covalent radius
    cov=cov_rad()
    for i in range(len(cov)):
        if cov[i][0]==ele:
           r_int=cov[i][1]
    return r_int

def check(pdb,res_num, alt, A1, A2, A3, nr_b):
#performs all checks bevore continue
        test_cret = True
        atoms = ceck_bond_atoms(pdb,res_num,alt ,A1, A2, A3)
        if atoms == False:
           return atoms
        b_number=nr_bonds(pdb,res_num,alt ,A3)
        if b_number > nr_b:
           bond = False
           return bond
        return test_cret
    
    
def ceck_bond_atoms(pdb,res_num,alt, A1, A2, A3):
## cecks if all needet atoms to ad a D are existing in residue
        com = True
        res = []
        for i in range(len(pdb[res_num])):
#           print(len(list(pdb[res_num][i][3])))
            if pdb[res_num][i][3]==alt or\
            len(list(pdb[res_num][i][3]))==0:
#              print(pdb[res_num][i][3],len(list(pdb[res_num][i][3])))
               res.append(pdb[res_num][i][2])
            elif len(list(alt))==0:
               if pdb[res_num][i][3]=="A":
                  res.append(pdb[res_num][i][2])
#       print("A1, A2, A3",A1, A2, A3)
#       print("res",res)
        if A1 not in res or A2 not in res or A3 not in res:
                com=False
        
        return com    
    
def ceck_bond_atoms_N(pdb,res_num,res_num2,alt,A1, A2, A3):
## cecks if all needet atoms to ad a D are existing in residue for N
        com = True
        res = []
        res2=[]
#       print("alt!!",alt)
        for i in range(len(pdb[res_num])):
            if pdb[res_num][i][3]==alt or\
            len(list(pdb[res_num][i][3]))==0:
 #             print(pdb[res_num][i][3],len(list(pdb[res_num][i][3])))
               res.append(pdb[res_num][i][2])
            elif len(list(alt))==0:
               if pdb[res_num][i][3]=="A":
                  res.append(pdb[res_num][i][2])

        if res_num2 < 0:
           com = False 
           return com
        for i in range(len(pdb[res_num2])):
            if pdb[res_num2][i][3]==alt or\
            len(list(pdb[res_num2][i][3]))==0:
 #             print(pdb[res_num][i][3],len(list(pdb[res_num][i][3])))
               res2.append(pdb[res_num2][i][2])
            elif len(list(alt))==0:
               if pdb[res_num2][i][3]=="A":
                  res2.append(pdb[res_num2][i][2])

        if A1 not in res or A2 not in res:
                com=False
        if A3 not in res2:
                com=False
        return com
    
    
def nr_bonds(pdb,res_num,alt,A3):
##counts number of bonds for atom
#       cov = cov_rad()
        nr=0
        atom1, ele1 =anam_to_coord_ele(pdb,res_num,alt,A3)
        cov1=ele_to_cov_rad(ele1)
        for i in range(len(pdb[res_num])):
#         print(alt)
#         print(pdb[res_num][i][3])
          if  pdb[res_num][i][3]==alt or len(list(pdb[res_num][i][3]))==0:
#           print(alt)
#           print("test alt",alt)
            atom2 =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
            ele2 = pdb[res_num][i][14]
            cov2=ele_to_cov_rad(ele2)
            d = CDist2(atom1,atom2)
            dmax= cov1+cov2+0.4
            if 0.5 <= d <= dmax:
               nr= nr + 1
#              print(nr)
          elif len(list(alt))==0:
            if pdb[res_num][i][3]=="A":
               atom2 =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
               ele2 = pdb[res_num][i][14]
               cov2=ele_to_cov_rad(ele2)
               d = CDist2(atom1,atom2)
               dmax= cov1+cov2+0.4
               if 0.5 <= d <= dmax:
                  nr= nr + 1
#                 
        return nr    
    
    
def anam_to_coord_ele(pdb,res_num,alt,anam):
#strip coordinates from atom name and res number
#       print("res_alt",alt)
#       print("serch anam",anam)
#       print("pdb(res_num]",pdb[res_num])        
#       print("len(list(alt))",len(list(alt)))
        xyz=[0,0,0]
        ele="dum"
        for i in range(len(pdb[res_num])):
#          print(pdb[res_num][i])
           if pdb[res_num][i][2].strip()==anam and\
           pdb[res_num][i][3]==alt:
              xyz =[0,0,0]
              xyz =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
              ele = pdb[res_num][i][14]
           elif len(list(alt))==0:
               if pdb[res_num][i][3]=="A":
                  xyz =[0,0,0]
                  xyz =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
                  ele = pdb[res_num][i][14]

        return xyz, ele  
    
def ele_to_cov_rad(ele):
## givs covalend radius for elment
        cov = cov_rad()
        r=-100
        for i in range(len(cov)):
            if ele == cov[i][0]:
               r = cov[i][1]
        return r   
    
def ZtoXYZ(axyz,bxyz,cxyz,R,Ad,Dd):
## calculates coordinates for atom in distance r to cxyz with angle ad
## to bxyz with dihedral Dd to axyz 
        dxyz = [0,0,0]
        rtodeg = 57.2957795

        ## first check if atoms are linear (of yes something is wrong)
        tangle=CAngle(axyz,bxyz,cxyz)
        if abs(tangle)<= 0.1 or 179.9<= abs(tangle) <= 180.0:
                print("The atoms are collinear", tangle)
                sys.exit()
        
        # transforme vrom degree to rad
        A = Ad/rtodeg
        D = Dd/rtodeg
        
        ## Calculate the coordinates in a simple coordinate system
        dxyz[0] = (-R)* np.sin(A)*np.sin(D)
        dxyz[1] =   R * np.cos(A)
        dxyz[2] =   R * np.sin(A)*np.cos(D)
        b   = np.sqrt(CDist2(bxyz,cxyz))
        ab1 = np.sqrt(CDist2(axyz,bxyz))
        ang = CAngle(axyz,bxyz,cxyz)/rtodeg
        a2  = b - np.cos(ang)*ab1
        a3  = np.sin(ang)*ab1
        a1  = 0.0000

        ## Now, atom D is transformed into the original coordinate system
        ## 1st rotation
        tv = np.array([0.00, b, 0.00])
        bcv = np.array([bxyz[0]-cxyz[0], bxyz[1]-cxyz[1],bxyz[2]-cxyz[2]])
        rv = Cross(tv, bcv)
        rv = Normlz(rv)
        rv = np.array([rv[0],rv[1],rv[2]])
        zerov = ZeroVector()
        phi = CAngle(tv,zerov,bcv)/rtodeg
        an = [0,0,0]
        an[0]=(rv[0]*rv[0]+(1-rv[0])*(1+rv[0])*np.cos(phi))*a1+\
            (rv[0]*rv[1]*(1-np.cos(phi))-rv[2] *np.sin(phi))*a2+\
            (rv[0]*rv[2]*(1-np.cos(phi))+rv[1] *np.sin(phi))*a3

        an[1]=(rv[0]*rv[1]*(1-np.cos(phi))+rv[2] *np.sin(phi))*a1+\
            (rv[1]*rv[1]+(1-rv[1])*(1+rv[1])*np.cos(phi))*a2+\
            (rv[1]*rv[2]*(1-np.cos(phi))-rv[0] *np.sin(phi))*a3

        an[2]=(rv[0]*rv[2]*(1-np.cos(phi))-rv[1] *np.sin(phi))*a1+\
            (rv[1]*rv[2]*(1-np.cos(phi))+rv[0] *np.sin(phi))*a2+\
            (rv[2]*rv[2]+(1-rv[2])*(1+rv[2])*np.cos(phi))*a3

        dn = [0,0,0]
        dn[0]=(rv[0]*rv[0]+(1-rv[0])*(1+rv[0])*np.cos(phi))*dxyz[0]+   \
              (rv[0]*rv[1]*(1-np.cos(phi))-rv[2] *np.sin(phi))*dxyz[1]+\
              (rv[0]*rv[2]*(1-np.cos(phi))+rv[1] *np.sin(phi))*dxyz[2]

        dn[1]=(rv[0]*rv[1]*(1-np.cos(phi))+rv[2] *np.sin(phi))*dxyz[0]+\
              (rv[1]*rv[1]+(1-rv[1])*(1+rv[1])*np.cos(phi))*dxyz[1]+   \
              (rv[1]*rv[2]*(1-np.cos(phi))-rv[0] *np.sin(phi))*dxyz[2]

        dn[2]=(rv[0]*rv[2]*(1-np.cos(phi))-rv[1] *np.sin(phi))*dxyz[0]+\
              (rv[1]*rv[2]*(1-np.cos(phi))+rv[0] *np.sin(phi))*dxyz[1]+\
              (rv[2]*rv[2]+(1-rv[2])*(1+rv[2])*np.cos(phi))*dxyz[2]

        dxyz = [dn[0], dn[1],dn[2]]
        
 
        # 2nd rotation
        tv[0]=axyz[0]-cxyz[0]
        tv[1]=axyz[1]-cxyz[1]
        tv[2]=axyz[2]-cxyz[2]
        phi=CDihed(tv,bcv,zerov,an)/rtodeg
        bcv = Normlz(bcv)

        dn[0]=(bcv[0]*bcv[0]+(1-bcv[0])*(1+bcv[0])*np.cos(phi))*dxyz[0]+\
           (bcv[0]*bcv[1]*(1-np.cos(phi))-bcv[2]*np.sin(phi)) *dxyz[1]+\
           (bcv[0]*bcv[2]*(1-np.cos(phi))+bcv[1]*np.sin(phi)) *dxyz[2]
        
        dn[1]=(bcv[0]*bcv[1]*(1-np.cos(phi))+bcv[2]*np.sin(phi)) *dxyz[0]+\
           (bcv[1]*bcv[1]+(1-bcv[1])*(1+bcv[1])*np.cos(phi))*dxyz[1]+\
           (bcv[1]*bcv[2]*(1-np.cos(phi))-bcv[0]*np.sin(phi)) *dxyz[2]
        
        dn[2]=(bcv[0]*bcv[2]*(1-np.cos(phi))-bcv[1]*np.sin(phi)) *dxyz[0]+\
           (bcv[1]*bcv[2]*(1-np.cos(phi))+bcv[0]*np.sin(phi)) *dxyz[1]+\
           (bcv[2]*bcv[2]+(1-bcv[2])*(1+bcv[2])*np.cos(phi))*dxyz[2]

        dxyz = [dn[0], dn[1],dn[2]]

        #Final translation
        dxyz[0]=dxyz[0]+cxyz[0]
        dxyz[1]=dxyz[1]+cxyz[1]
        dxyz[2]=dxyz[2]+cxyz[2]


        return dxyz    



def C2Angle(x,y):
#Calculates the angle between x and y
#Answer in degrees
        #Calculate the angle between x and y
        rtodeg = 57.2957795
        C2angle=(ScalPr(x,y)/(len3dvec(x)*len3dvec(y)))
        if C2angle == 1 :
           C2angle=0
        elif C2angle== -1:
                C2angle=180
        else:
                C2angle=np.arccos(C2angle)*rtodeg
        return C2angle

def ScalPr(x,y):
#calculate the scalar product
        pro= x[0]*y[0]+x[1]*y[1]+x[2]*y[2]
        return pro

def CDihed(x,y,z,w):
#Calculate the dihedral angle x-y-z-w
#Answer in degrees between -180 and +180
        #Set v1=y-x, v2=z-y, v3=w-z
        v1=[0,0,0]
        v2=[0,0,0]
        v3=[0,0,0]
        for i in range(3):
           v1[i]=y[i]-x[i]
           v2[i]=z[i]-y[i]
           v3[i]=w[i]-z[i]

        #Calculate the normal vectors n1 and n2
        n1=Cross(v1,v2)
        n2=Cross(v2,v3)
        
        #Calculate the torsion angle;
        #The sign is determined by the sign of v1.n2
        CDihed=C2Angle(n1,n2)
        if ScalPr(v1,n2) < 0:
                CDihed=-CDihed
        return  CDihed


def Cross(x1,x2):
#Calculates the cross product x3 = x1 x x2
        x3 = [0,0,0] 
        x3[0]=x1[1]*x2[2]-x2[1]*x1[2]
        x3[1]=x1[2]*x2[0]-x2[2]*x1[0]
        x3[2]=x1[0]*x2[1]-x2[0]*x1[1]
        return x3

def ZeroVector():
## creates ZeroVector
        zerov = np.array([0,0,0])
        return zerov

def Normlz(xyz):
# Normalise xyz 
        temp = 1/len3dvec(xyz)
        for i in range(len(xyz)):
            xyz[i]=xyz[i]*temp

        return xyz



def CAngle(x,y,z):
# calculate angle between a,b,c in degree
    x1=[0,0,0]
    x2=[0,0,0]
    for i in range(len(x)):
        x1[i]=y[i]-x[i]
        x2[i]=y[i]-z[i]
    Cangle=C2Angle(x1,x2)
    return Cangle


def int_CH(pdb,res_num,alter,map,mesh):
#integrates for all CH atoms in list
        v_CH=[]
        l_H=[]
        av_CH=0
        a_list=[]
        r_int=def_r_int("H")
        for j in range(len(pdb[res_num])):
           if pdb[res_num][j][14]=="C":
             a_list.append(pdb[res_num][j][2])
        a_list=list(set(a_list))
#       print("a_list",a_list)
        for a in range(len(alter)):
          for i in range(len(a_list)):
             l_H = l_H + find_CH(pdb,res_num,alter[a],a_list[i]) 
        for i in range(len(l_H)):
         v_CH.append(int_around_point(l_H[i],map,r_int, mesh))

        av_CH=average(v_CH)
        sig_CH=std(v_CH, av_CH)
#       print("int_CH",av_CH,"sig_CH",sig_CH)
        return av_CH, sig_CH

def find_CH(pdb,res_num,alt,A3):
##find CH bond du A3
        H_coord=[]
        atom1, ele1 =anam_to_coord_ele(pdb,res_num,alt,A3)
        cov1=ele_to_cov_rad(ele1)
        for i in range(len(pdb[res_num])):
          if  pdb[res_num][i][3]==alt or len(list(pdb[res_num][i][3]))==0:
            atom2 =np.array([pdb[res_num][i][8],pdb[res_num][i][9],pdb[res_num][i][10]])
            ele2 = pdb[res_num][i][14]
            cov2=ele_to_cov_rad(ele2)
            d = CDist2(atom1,atom2)
            dmax= cov1+cov2+0.4
            if 0.5 <= d <= dmax:
               if pdb[res_num][i][14] == "H" or pdb[res_num][i][14]=="D":
                  H_coord.append(atom2)
        return H_coord


def std(i_list,mid):
## calculate standard deviation
        l_sqr=[]
        for i in range(len(i_list)):
            num=(i_list[i]-mid)**2
            l_sqr.append(num)
        s2=sum(l_sqr)/(len(l_sqr))
        s=np.sqrt(s2)
        return s


def average(i_list):
## calculate average of list
        ave=sum(i_list) / len(i_list)
        return ave




### statistics
def stat_H(pdb,alerts):
## do statistics on all addet hydrogen if existing
    if map_ex==False:
       return alerts

    for i in range(len(pdb)):
        int_H_l=[]
#       if i >0:
#          if pdb[i-1][0][6] == pdb[i][0][6]-1:
#            for k in range(len(pdb[i-1])):
#              if pdb[i-1][k][14]=="H" or pdb[i-1][k][14]=="D":
#                 if len(pdb[i-1][k])==17:
#                    print(pdb[i-1][k])
#                    int_H_l.append(pdb[i-1][k][16])
        for k in range(len(pdb[i])):
               if pdb[i][k][14]=="H" or pdb[i][k][14]=="D":
                  if len(pdb[i][k])==17:
                     int_H_l.append(pdb[i][k][16])
#       if i < len(pdb)-1:
#        if pdb[i+1][0][6] == pdb[i][0][6]+1:
#          for k in range(len(pdb[i+1])):
#              if pdb[i+1][k][14]=="H" or pdb[i+1][k][14]=="D":
#                 if len(pdb[i+1][k])==17:
#                   int_H_l.append(pdb[i+1][k][16])

        if len(int_H_l)>=2:
           av=average(int_H_l)
           sig=std(int_H_l,av)
#          print("AV",av,"sig",sig)
        
         
           for k in range(len(pdb[i])):
               if pdb[i][k][14]=="H" or pdb[i][k][14]=="D":
                  if len(pdb[i][k])==17:
                     if pdb[i][k][16] <= av-2.0*sig:
                       string=str(pdb[i][k][2])+" in res "+str(pdb[i][k][4])+"-"+str(pdb[i][k][6])+"-"+str(pdb[i][k][5])+" might not exisitng"
                       alerts[7].append(string)
                        



    return alerts
def vdw_rad():
    #vdw radie of elments

        l_vdw=[["H",1.20],["D",1.20],["C",1.70],["O",1.52],["N",1.55],["CL",1.75],["F",1.47],["BR",1.85],["dum",0.0]]
#       h_vdw=1.20
#       c_vdw=1.70
#       o_vdw=1.52
#       n_vdw=1.55
#       s_vdw=1.80
#       cl_vdw= 1.75
#       f_vdw=1.47
#       br_vdw= 1.85
        return l_vdw

def ele_to_vdw_rad(ele):
## givs covalend radius for elment
        vdw = vdw_rad()
        r=-100
        for i in range(len(vdw)):
            if ele == vdw[i][0]:
               r = vdw[i][1]
        return r


def find_clash(pdb):
#search for vdw and cov clashes beetween hydrogen atoms
        clash_cov=[]
        clash_vdw=[]
        H="H"
        cov_H=ele_to_cov_rad(H)
        vdw_H=ele_to_vdw_rad(H)
#        print("H",cov_H,vdw_H)
        for i in range(len(pdb)):
           for k in range(len(pdb[i])):
              if pdb[i][k][14]=="H" or pdb[i][k][14]=="D":
#                print(pdb[i][k])
                 for j in range(len(pdb)):
                   if j != i:
                     for l in range(len(pdb[j])):
                        if pdb[j][l][14]=="H" or pdb[j][l][14]=="D":
                           if pdb[j][l][3]==pdb[i][k][3] or  pdb[j][l][3]=="" or pdb[i][k][3] =="":
                               H1=[pdb[i][k][8],pdb[i][k][9],pdb[i][k][10]]
                               H2=[pdb[j][l][8],pdb[j][l][9],pdb[j][l][10]]
                               dist=CDist2(H1,H2)
                               if dist <= cov_H*2:
                                  clash_cov.append([pdb[i][k],pdb[j][l]])
                               if dist <= vdw_H*2:
                                  clash_vdw.append([pdb[i][k],pdb[j][l]])
        print("number of cov clashes ", int(len(clash_cov)/2))
        print("number of vdw clashes ", int(len(clash_vdw)/2))
        clash_cov_sort=[clash_cov[0]]
        
        for i in range(len(clash_cov)):
            ext=False
            for k in range(len(clash_cov_sort)):
                if clash_cov[i][0][1]==clash_cov_sort[k][1][1] and  clash_cov[i][1][1]==clash_cov_sort[k][0][1]:
                   ext=True
                if clash_cov[i][0][1]==clash_cov_sort[k][0][1] and  clash_cov[i][1][1]==clash_cov_sort[k][1][1]:
                   ext=True
            if ext==False:
                clash_cov_sort.append(clash_cov[i])
        clash_cov=clash_cov_sort
        for i in range(len(clash_cov)):
            print("atom1 ",clash_cov[i][0][2],clash_cov[i][0][3],clash_cov[i][0][4],clash_cov[i][0][6],clash_cov[i][0][5],clash_cov[i][0][1])
            print("atom2  ",clash_cov[i][1][2],clash_cov[i][1][3],clash_cov[i][1][4],clash_cov[i][1][6],clash_cov[i][1][5],clash_cov[i][1][1])

        return clash_cov,clash_vdw






    
    
## actual runn of script
pdb_read, work_pdb = read_pdb(input_pdb)
pdb = atm_to_res(pdb_read)
if map_ex==True:
        map= read_grid(input_grid)
        map.setup()

pdb, alerts = treat_res(pdb, map)

alerts=stat_H(pdb,alerts)

for i in range(len(alerts)):
   print(" ")
   for k in range(len(alerts[i])):
      print("alert", alerts[i][k])
   print(len(alerts[i]))
clash_cov,clash_vdw=find_clash(pdb)

##for i in range(len(pdb)):
##    for k in range(len(pdb[i])):
##        print(pdb[i][k])

write_pdb(pdb, out_file)
time_ende = time.time()
print("program ends normally after "'{:5.3f}s'.format(time_ende-time_start),
" or ", '{:5.2f}min'.format((time_ende-time_start)/60))
