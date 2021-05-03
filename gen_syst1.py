#!/usr/bin/env python3
#######!/software/sse/easybuild/prefix/software/Anaconda3/5.3.0-extras-nsc1/bin-wrapped/python 
# -*- coding: utf-8 -*-

# run with protonation: gen_syst1.py qm.tex pdb P
# run without protonation: gen_syst1.py qm.tex pdb 

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
        qm="qm.txt"


if len(sys.argv)>=3:
        input_pdb=sys.argv[2]
        print("read from line 2")
else:
        input_pdb="comqum.pdb"

f_p=False
if len(sys.argv)>=4:
        pq=sys.argv[3]
        print(pq)
        if pq == "P" or pq == "p":
                f_p=True
                print("protonation is on")




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
            
            if idef=="ATOM" or idef=="HETATM":

                l_pdb=[]
                a       = line[0:6]
                a_num   = int(line[6:11].strip())
                a_name  = line[12:16]
                alter   = line[16]
                res_name= line[17:20]
                chain   = line[21]
                res_num = int(line[22:26].strip())

                l_pdb.append(a)
                l_pdb.append(a_num)
                l_pdb.append(a_name)
                l_pdb.append(alter)
                l_pdb.append(res_name)
                l_pdb.append(chain)
                l_pdb.append(res_num)
                qm_atom.append(l_pdb)

            if len(idef2)>=1:
              if idef=="RESI" or idef2[0]=="RESI":
                s_line=line.split()
                l_resi=[]
                l_resi.append(s_line[1])
                l_resi.append(int(s_line[2]))
                qm_resi.append(l_resi)
                

              if idef=="ANAM" or idef2[0]=="ANAM":
                s_line=line.split()
                l_resi=[]
                l_resi.append(s_line[1])
                l_resi.append(int(s_line[2]))
                if len(s_line)<4:
                        s_line.append("  ")
                l_resi.append((s_line[3]))
                qm_anam.append(l_resi)

              if idef=="ANUM" or idef2[0]=="ANUM":
                s_line=line.split()
                l_resi=[]
                l_resi.append(int(s_line[1]))
                qm_anum.append(l_resi)

              if idef=="RNUM" or idef2[0]=="RNUM":
                s_line=line.split()
                l_resi=[]
                l_resi.append(int(s_line[1]))
                l_resi.append(int(s_line[2]))
                qm_rnum.append(l_resi)


        print("close file:", qm)
        return qm_atom, qm_resi, qm_anam, qm_anum, qm_rnum


def assign_atm_nr(qm_atom, qm_resi, qm_anum, qm_anam, qm_rnum ,pdb):
        atm_nr=[]
        atm_qm=[]
        if len(qm_atom)>0:
           for i in range(len(qm_atom)):
               res_num= qm_atom[i][6]
               chain  = qm_atom[i][5].strip()
               a_name = qm_atom[i][2].strip()
               for j in range(len(pdb)):
                 for k in range(len(pdb[j])):
                   res_num_pdb = pdb[j][k][6]
                   chain_pdb   = pdb[j][k][5].strip()
                   a_name_pdb  = pdb[j][k][2].strip()
                   if res_num == res_num_pdb and chain == chain_pdb and a_name == a_name_pdb:
                        atm_nr.append(pdb[j][k][1])
                        atm_qm.append(pdb[j][k])


        if len(qm_resi)>0:
            for i in range(len(qm_resi)):
               res_num= qm_resi[i][1]
               chain  = qm_resi[i][0].strip()
               for j in range(len(pdb)):
                 for k in range(len(pdb[j])):
                   res_num_pdb = pdb[j][k][6]
                   chain_pdb   = pdb[j][k][5].strip()
#                  print(res_num, res_num_pdb)
                   if res_num == res_num_pdb and chain == chain_pdb:
#                       print(work_pdb[j][1])
                        atm_nr.append(pdb[j][k][1])
                        atm_qm.append(pdb[j][k])
        
        if len(qm_anum)>0:
           for i in range(len(qm_anum)):
               a_num=qm_anum[i][0]
               for j in range(len(pdb)):
                 for k in range(len(pdb[j])):
                   a_num_pdb=pdb[j][k][1]
                   if a_num == a_num_pdb:
                        atm_nr.append(pdb[j][k][1])
                        atm_qm.append(pdb[j][k])


        if len(qm_anam)>0:
           for i in range(len(qm_anam)):
               res_num= qm_anam[i][1]
               chain  = qm_anam[i][2].strip()
               a_name = qm_anam[i][0].strip()
               for j in range(len(pdb)):
                 for k in range(len(pdb[j])):
                   res_num_pdb = pdb[j][k][6]
                   chain_pdb   = pdb[j][k][5].strip()
                   a_name_pdb  = pdb[j][k][2].strip()
                   if res_num == res_num_pdb and chain == chain_pdb and a_name == a_name_pdb:
                        atm_nr.append(pdb[j][k][1])
                        atm_qm.append(pdb[j][k])
        
        if len(qm_rnum)>0:
           for i in range(len(qm_rnum)):
               low = qm_rnum[i][0]
               up  = qm_rnum[i][1]
#              print(low,up)
               for j in range(len(pdb)):
                   a_num_pdb=pdb[j][k][1]
                 # print(a_num_pdb, low,up)
                   if low <= a_num_pdb <= up:
                        atm_nr.append(pdb[j][k][1])
                        atm_qm.append(pdb[j][k])



        atm_nr = list(set(atm_nr))
        atm_nr= sorted(atm_nr)
        atom_qm=[]
#       print(atm_nr)
#       print(atm_qm)
        for i in range(len(atm_nr)):
            for j in range(len(atm_qm)):
#               print(atm_qm[j])
                if atm_nr[i]== atm_qm[j][1]:
                   atom_qm.append(atm_qm[j])
                   break
        
        
        return atm_nr, atom_qm
def len3dvec(vec):
## calculates lengh of a 3D vecor
## input as list
        a = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        return a

def find_truncation(atm_nr, atom_qm, pdb):
## serche for trancuation atoms in radius aof 1.8 A around atoms
        trunc_nr = []
        trunc_atm= []
        print(atom_qm)
        for i in range(len(atom_qm)):
            vec1=[0,0,0]
            vec1[0] = atom_qm[i][8]
            vec1[1] = atom_qm[i][9]
            vec1[2] = atom_qm[i][10]
            for j in range(len(pdb)):
             for k in range(len(pdb[j])):
               vec2=[0,0,0]
               vec2[0] = pdb[j][k][8]
               vec2[1] = pdb[j][k][9]
               vec2[2] = pdb[j][k][10]
               vec=[0,0,0]
               for h in range(len(vec)):
                 vec[h]=vec1[h]-vec2[h]
               lenvec= len3dvec(vec)
               if 0.5 < lenvec <= 1.8:
                  trunc_nr.append(pdb[j][k][1])
                  trunc_atm.append(pdb[j][k])
#       print("trunc_nr",trunc_nr)

        trunc_nr=sorted(list(set(trunc_nr)))
#       print(len(trunc_nr),"trunc_nr2",trunc_nr)
#       print("atm_nr",atm_nr)
        trunc_nr = list(set(trunc_nr)-set(atm_nr))
        trunc_atom=[]
#       print("trunc_atm",trunc_atm)
#       print("trunc_nr",trunc_nr)
        for i in range(len(trunc_nr)):
            for j in range(len(trunc_atm)):
                if trunc_nr[i]==trunc_atm[j][1]:
                   trunc_atom.append(trunc_atm[j])
                   break
        return trunc_nr, trunc_atom



def sum_list(atm_nr, atom_qm, trunc_nr, trunc_atm):
## define list of QM atom and truncation atoms
        f_trunc=[]
        f_all= atm_nr + trunc_nr
        for i in range(len(trunc_atm)):
           if trunc_atm[i][14].strip() != "H":
                f_trunc.append(trunc_atm[i][1])
#          if trunc_atm[i][2].strip() != "H":
#               f_trunc.append(trunc_atm[i][1])
        
        f_all = sorted(list(set(f_all)))
        f_trunc=sorted(list(set(f_trunc)))
        f_qm=list(set(f_all)-set(f_trunc))
        f_qm = sorted(list(set(f_qm)))

        return f_all, f_qm, f_trunc


def write_syst1(f_all,f_trunc):

        syst1=open("syst1","w")
        syst1.write("coment  \n")
        f_comp=compress_list(f_all)
        for i in range(len(f_comp)):
            if len(f_comp[i])==1:
                string=str(f_comp[i][0])+" \n"
            if len(f_comp[i])==2:
                string = str(f_comp[i][0])+"-"+str(f_comp[i][1])+" \n"
            syst1.write(string)
        syst1.write(" \n")
        for i in range(len(f_trunc)):
            string=str(f_trunc[i])+" \n"
            syst1.write(string)
        syst1.write(" \n \n \n")
        syst1.close()
        print("syst1 file is ready")

def compress_list(list_in):
        list_out=[]
#       print("len(list_in)",len(list_in))
        second = False
        for i in range(len(list_in)):
#           print(second, i, )

            if i == 0:
               pm = -1
               pp = list_in[i+1]-1
               p = list_in[i]
            if 0 < i < len(list_in)-1:
               pm = list_in[i-1]+1
               pp = list_in[i+1]-1
               p = list_in[i]
            if  i == len(list_in)-1:
               pm = list_in[i-1]+1
               pp = -1
               p = list_in[i]


#           print(i, list_in[i],p,pm,pp, second)
            if second == False:
               if p != pm and p == pp:
                  list_out.append([list_in[i]])
                  second = True
               if p != pm and p != pp:
                  list_out.append([list_in[i]])
                  second = False
               continue
            if second == True :
               if p == pm and p == pp :
#                 print("###############")
                  continue
               if p == pm and p != pp :
#                 print("!!!!!!!!!!!")
                  list_out[len(list_out)-1].append(list_in[i])
                  second = False
               if p != pm and p != pp :
#                 print("¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤")
                  list_out.append([list_in[i]])
                  second = False

#       print("len(list_out)",len(list_out))
        return list_out



def write_restart_anam(f_qm,work_pdb):
## write input to restart this programm
        r_anam=open("restart_anam","w")
        w_pdb=open("pdb_syst1","w")
#       w_label=open("label_vis.sh","w")
        atm_nr=0
        for i in range(len(f_qm)):
            for j in range(len(work_pdb)):
                if f_qm[i]==work_pdb[j][1]:
#                  print(work_pdb[j])
                   string="ANAM "+ work_pdb[j][2]+" "+str(work_pdb[j][6])+" " + work_pdb[j][5] + "\n"
                   s_pdb, atm_nr= write_pdbi_line(work_pdb[j], atm_nr)
#                  s_label=prep_vis(work_pdb[j])
                   r_anam.write(string)
                   w_pdb.write(s_pdb)
#                  w_label.write(s_label)

def write_rename(f_all,f_trunc,work_pdb):
    w_label=open("label_vis.sh","w")
    for i in range(len(f_all)):
            for j in range(len(work_pdb)):
                if f_all[i]==work_pdb[j][1]:
                   ele=work_pdb[j][14]
                   if f_all[i] in f_trunc:
                      ele="H"

                   s_label=prep_vis(work_pdb[j],ele)
                   w_label.write(s_label)

def prep_vis(pdb,ele):
#preper bash file for the corect labeling of the pdb file from mimic
#[0,   1  ,   2   ,  3   ,    4    ,   5  ,    6,  ,  7   , 8, 9,10
#[a, a_num, a_name, alter, res_name, chain, res_num, inser, x, y, z
#11 ,  12  ,  13  , 14 ,  15
#occ, b_fac, seg_i, ele, charge]

        res_nam= pdb[4]
        res_num= pdb[6]
        chain=   pdb[5]
        string="sed -i \"0,/"
        string = string + str('{:2}'.format(ele))
#       string = string + "  ???     1/s//"
        string = string + "  ???     1/s/"
        string = string + str('{:2}'.format(ele))
        string = string + "  ???     1/"


        string = string + str('{:2}'.format(ele))
        string = string + "  "
        string = string + str(res_nam) + " "
        string = string + str(chain)#+ " "
        string = string + str('{:4}'.format(res_num))
        string = string + "/g\" pdb"
        string = string + "\n"
        return string

def write_pdbi_line(pdb, atm_nr):
## wite out PDB file
        atm_nr = atm_nr + 1
        string = str('{:6}'.format(pdb[0]))
        string = string + str('{:5.0f}'.format(atm_nr))
        string = string + "  "
        string = string + str('{:3s}'.format(str(pdb[2])))
        string = string + str('{:1}'.format(pdb[3]))
        string = string + str('{:3}'.format(pdb[4]))
        string = string + str('{:>2}'.format(pdb[5]))
        string = string + str('{:4}'.format(pdb[6]))
        string = string + str('{:1}'.format(pdb[7]))
        string = string + "   "
        string = string + str('{:8.3f}'.format(pdb[8]))
        string = string + str('{:8.3f}'.format(pdb[9]))
        string = string + str('{:8.3f}'.format(pdb[10]))
        string = string + str('{:6.2f}'.format(pdb[11]))
        string = string + str('{:6.2f}'.format(pdb[12]))
        string = string + str('{:>7}'.format(pdb[13]))
        string = string + str('{:>5}'.format(pdb[14]))


        string = string + "\n"
        return string, atm_nr


def atm_to_res(pdb_read):
# sort pdb in residues
        pdb=[]
        pdb.append([pdb_read[0]])
        for i in range(1,len(pdb_read)):
           if pdb_read[i][6] != pdb_read[i-1][6]:
#             print(len(pdb))
              pdb.append([pdb_read[i]])
           if pdb_read[i][6] == pdb_read[i-1][6]:
              pdb[len(pdb)-1].append(pdb_read[i])

        return pdb

def proto_qm(pdb,atom_qm,f_p):
# protonates qm system if protons are not existing
#       print("protonate syst1")
        if f_p== False:
           return pdb
        print("protonate syst1")
        for i in range(len(atom_qm)):
#          print("protonate syst1")

#          print(atom_qm[i])
           xyz=[atom_qm[i][8],atom_qm[i][9],atom_qm[i][10]]
#          print(xyz)
           ele=atom_qm[i][14]
#          print(ele)
#          print(atom_qm[i])
           if ele=="C":
              pdb=prot_c(pdb,xyz,ele,atom_qm[i])
           if ele=="N":
              pdb=prot_n(pdb,xyz,ele,atom_qm[i])
           if ele=="O":
              pdb=prot_o(pdb,xyz,ele,atom_qm[i])
        return pdb

def prot_n(pdb,xyz,ele,atom_qm):
## identyfy atom pype and protanes it
        bond_atom=f_bond_atoms(pdb,atom_qm)
        n_bonds=def_n_bonds(xyz,ele,bond_atom)
        if 1.5< n_bonds <= 2.5 and len(bond_atom)==2:
           pdb= ad_NH_s(pdb,bond_atom[0],bond_atom[1],atom_qm)
           return pdb
        if n_bonds==1 and len(bond_atom)==1:
           pdb=ad_NH2(pdb,bond_atom[0],atom_qm)
           return pdb
        return pdb


def ad_NH2(pdb,  A2, A3):
# ad 2H as NH2 group to A3
#    print("¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤")
     r   = 1.012
     Ad  = 120
     Dd  = 180
     Dd1 = 0
     A1=find_A1(pdb,A3,A2) 
     A3xyz  = [A3[8],A3[9],A3[10]]
     A2xyz  = [A2[8],A2[9],A2[10]]
     A1xyz  = [A1[8],A1[9],A1[10]]

     adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
     adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd1)
     nam=list(A3[2])
     name = str(nam[0])

     name1 = "H1"+name
     name2 = "H2"+name
     pdb = add_to_pdb(pdb, A3, adH1, name1)
     pdb = add_to_pdb(pdb, A3, adH2, name2)

     return pdb







def ad_NH_s(pdb, A1, A2, A3):
# ad H as NH (SP2) group to A3
#    print("atoms treatet")
#    print("A3",A3)
#    print("A2",A2)
#    print("A1",A1)
     r  = 1.032
     Ad = 120
     Dd = 180
     A3xyz  = [A3[8],A3[9],A3[10]]
     A2xyz  = [A2[8],A2[9],A2[10]]
     A1xyz  = [A1[8],A1[9],A1[10]]
     adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
     nam=list(A3[2])
     if len(nam)<=2:
        name=nam[0]
     else:
        name = nam[1]+nam[2]

     name1 = "H"+name

#    print("name1",name1)
#    print(adH1)
#    print(A3)
#    print("next step ad to pdb")
     pdb = add_to_pdb(pdb, A3, adH1, name1)
#    print("addet to pdb")
#    print("&&&&&&&&&&&&&&&TEST&&&&&&&&&&&&")
     return pdb


def prot_o(pdb,xyz,ele,atom_qm):
## identyfy atom pype and protanes it
        bond_atom=f_bond_atoms(pdb,atom_qm)
        n_bonds=def_n_bonds(xyz,ele,bond_atom)
#       print("atom_qm",n_bonds,len(bond_atom), atom_qm)
#       if atom_qm[4]=="HOH":
#          print(bond_atom)
        if n_bonds >= 2:
           return pdb
        if n_bonds==1 and len(bond_atom)==1:
           pdb=ad_OH(pdb,bond_atom[0],atom_qm)
        if n_bonds==0 and atom_qm[4]=="HOH":
           pdb = ad_OH_water(pdb,atom_qm)
        return pdb


def ad_OH_water(pdb,atom_qm):
## protonates water in direction of posible hydrogen bonds
        print(atom_qm)
        l_vdw=vdw()
        l_h_bonds=find_h_bond(atom_qm,pdb,l_vdw)
#       for i in range(len(l_h_bonds)):
#          print("l_h_bonds",l_h_bonds[i])
        if len(l_h_bonds)>=1:
#               print("######################adH1#######################")
                adH1=find_short_H(atom_qm,l_h_bonds)
                if adH1[0]==0:
                      return pdb
                name1 ="H1"+atom_qm[2]
                pdb = add_to_pdb(pdb, atom_qm, adH1, name1)
                if len(l_h_bonds)>=2: 
#                  print("______________________________adH2_______________________")
                   adH2=find_short_H(atom_qm,l_h_bonds)
                   if adH2[0]==0:
                      return pdb
                   name2 ="H2"+atom_qm[2]
                   pdb = add_to_pdb(pdb, atom_qm, adH2, name2)

        return pdb


def find_short_H(atom_qm,l_h_bonds):
## finds shortest istace in list of coordinates
        P1=[atom_qm[8],atom_qm[9],atom_qm[10]]
        d_min=1000000
        for i in range(len(l_h_bonds)):
            P2=[l_h_bonds[i][8],l_h_bonds[i][9],l_h_bonds[i][10]]
            dist=CDist2(P1,P2)
            if dist< d_min:
               exist=f_ex_h(pdb, P1, P2)
#              print("exist",exist)
#              if exist == True:
#                 print("l_h_bonds[i]",l_h_bonds[i])
               if exist == False:
                  d_min=dist
                  atom_min=l_h_bonds[i]

        P_s=[atom_min[8],atom_min[9],atom_min[10]]
        vec=twoP_to_vec(P1,P_s)
#       print("len3Dvec(vec)",len3dvec(vec))
        h_vec=resize_vec(vec,0.98)
#       print("len3Dvec(h_vec)",len3dvec(h_vec))
        
        coord=[P1[0]+h_vec[0],P1[1]+h_vec[1],P1[2]+h_vec[2]]
        if coord == P1:
           coord= [0]
        return coord

def f_ex_h(pdb, A1, A2):
#check if hydrogen is positiond in 30 deree to OHO distance
## false if proton not existing true if existing
     flag=False
     dh=1.5
     max_dist=0.6
     for j in range(len(pdb)):
      for k in range(len(pdb[j])):  
        if pdb[j][k][14]=="H" or pdb[j][k][14]=="D" :
           p0=[pdb[j][k][8],pdb[j][k][9],pdb[j][k][10]]
           d1= CDist2(A1,p0)
           d2=CDist2(A2,p0)
           ang=CAngle(A1,p0,A2)
           if d1 <= dh or d2 <=dh:
              if 120 <= ang <= 220:
#                print("ATOM H",pdb[j][k])
#                print("ang",ang)
                 flag=True


     return flag


def d_p_to_line(p0,p1,p2):
## calculate distance between a point0 and a line between p1 and p2
        if (p1[0] == p2[0] and p1[1] == p2[1] and p1[2] == p2[2]):
                d=0
        else:
                if (p2[0]-p1[0] != 0):
                        t=-((p1[0]-p0[0])*(p2[0]-p1[0]))/((abs(p2[0]-p1[0]))**2)
                elif (p2[1]-p1[2] != 0 ):
                       t=-((p1[1]-p0[1])*(p2[1]-p1[1]))/((abs(p2[1]-p1[1]))**2)
                elif ( p2[2]-p1[2] != 0):
                        t=-((p1[2]-p0[2])*(p2[2]-p1[2]))/((abs(p2[2]-p1[2]))**2)

                d2=((p1[0]-p0[0])+(p2[0]-p1[0])*t)**2+((p1[1]-p0[1])+(p2[1]-p1[1])*t)**2+((p1[2]-p0[2])+(p2[2]-p1[2])*t)**2
                d=d2**(0.5)
        return d




def resize_vec(vec,r):
        vec=Normlz(vec)
#       print("len3Dvec(vec)",len3dvec(vec))

        new_vec=[vec[0]*r,vec[1]*r,vec[2]*r]

        return new_vec

def find_h_bond(atom_qm,pdb,l_vdw):
        OO=l_vdw[3][1]+l_vdw[3][1]
        ON=l_vdw[3][1]+l_vdw[4][1]

        h_bond=[]
        O=[atom_qm[8],atom_qm[9],atom_qm[10]]
        for j in range(len(pdb)):
         for k in range(len(pdb[j])):
          if pdb[j][k][14]=="O" or pdb[j][k][14]=="N" :
#           print(pdb[k][6],l_water[i][j][6])
            if pdb[j][k][6]!=atom_qm[6]:
             acc=[pdb[j][k][8],pdb[j][k][9],pdb[j][k][10]]
             dist=CDist2(O,acc)
             if pdb[j][k][14]=="O":
                if dist <= OO:
                   h_bond.append(pdb[j][k])
             if pdb[j][k][14]=="N":
                if dist <= ON:
                   h_bond.append(pdb[j][k])

        return h_bond


def vdw():
        l_vdw=[["H",1.20],["D",1.20],["C",1.70],["O",1.52],["N",1.55],["CL",1.75],["F",1.47],["BR",1.85]]
#       h_vdw=1.20
#       c_vdw=1.70
#       o_vdw=1.52
#       n_vdw=1.55
#       s_vdw=1.80
#       cl_vdw= 1.75
#       f_vdw=1.47
#       br_vdw= 1.85
        return l_vdw


def ad_OH(pdb, A2, A3):
# ad H as OH group to A3
     r  = 0.978
     Ad = 106
     Dd = 180
     A1=find_A1(pdb,A3,A2)
     A3xyz  = [A3[8],A3[9],A3[10]]
     A2xyz  = [A2[8],A2[9],A2[10]]
     A1xyz  = [A1[8],A1[9],A1[10]]

     adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
     name1 ="H"+A3[2]
     pdb = add_to_pdb(pdb, A3, adH1, name1)

     return pdb



def prot_c(pdb,xyz,ele,atom_qm):
## identyfy atom pype and protanes it
        bond_atom=f_bond_atoms(pdb,atom_qm)

#       print(bond_atom)
#       print("xyz###",xyz)
        n_bonds=def_n_bonds(xyz,ele,bond_atom)
        res_num=atom_qm[6]
        alter=atom_qm[3]
        if n_bonds >= 4:
           return pdb
        if n_bonds==1:
           pdb = ad_CH3(pdb,res_num,alter,bond_atom[0],atom_qm)
        
        if n_bonds<=2.5 and len(bond_atom)==2:
           pdb = ad_CH2(pdb, res_num, alter, bond_atom[0], bond_atom[1], atom_qm) 
        if 2.5< n_bonds<= 3.5 and len(bond_atom)==2:
           pdb = ad_CH_ar(pdb,  bond_atom[0], bond_atom[1], atom_qm)
        if 2.5< n_bonds<= 3.5 and len(bond_atom)==3:
           pdb = ad_CH_R3(pdb,  bond_atom[0], bond_atom[1], atom_qm, bond_atom[2])

#       print("n_bonds",ele,n_bonds)
         
        return pdb

def ad_CH_R3(pdb, A1, A2, A3, A4):
        r  = 1.032
        Ad = 109
        Dd = 116
        A3xyz  = [A3[8],A3[9],A3[10]]
        A2xyz  = [A2[8],A2[9],A2[10]]
        A1xyz  = [A1[8],A1[9],A1[10]]
        A4xyz  = [A4[8],A4[9],A4[10]]

        adH1 = H_tert(A1xyz, A2xyz, A3xyz, A4xyz,r)
        name1 ="H"+A3[2]
        pdb = add_to_pdb(pdb, A3, adH1, name1)

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




def ad_CH_ar(pdb,  A1, A2, A3):
# ad H as CH (SP2) group to A3
     r  = 1.014
     Ad = 120
     Dd = 180
     A3xyz  = [A3[8],A3[9],A3[10]]
     A2xyz  = [A2[8],A2[9],A2[10]]
     A1xyz  = [A1[8],A1[9],A1[10]]
     adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
     name1 ="H"+A3[2]
     nam=list(A3[2])
     if len(nam)<=2:
        name=nam[1]
     else:
        name = nam[1]+nam[2]

     name1 = "H"+name
     pdb = add_to_pdb(pdb, A3, adH1, name1)

     return pdb









def ad_CH2(pdb, res_num, alter, A1, A2, A3):
# ad 2H as CH2 group to A3

           r = 1.095
           Ad = 109
           Dd= -121.5
           A3xyz  = [A3[8],A3[9],A3[10]]
           A2xyz  = [A2[8],A2[9],A2[10]]
           A1xyz  = [A1[8],A1[9],A1[10]]

           adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
           adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,-Dd)
           nam=list(A3[2])
           if len(nam)<=2:
              name=nam[1]
           else:
              name = nam[1]+nam[2]

           name1 = "1H"+name
           name2 = "2H"+name
           pdb = add_to_pdb(pdb, A3, adH1, name1)

           pdb = add_to_pdb(pdb, A3, adH2, name2)
           return pdb

def ad_CH3(pdb, res_num,alter, A2, A3):
# ad 3H as CH3 group to A3
#    for alt in range(len(alter)):
#       print("alter[alt]",alter[alt])
#       com =   check(pdb,res_num,alter[alt],A1, A2, A3, 1)
#       print(com)
#       if com == False:
#               return pdb
#       print( pdb[res_num])
#       print("A2", A2)
        A1=find_A1(pdb,A3,A2)
#       print("A1",A1)
        r   = 1.095
        Ad  = 109
        Dd  = 180
        Dd1 = 60
        A3xyz  = [A3[8],A3[9],A3[10]]
        A2xyz  = [A2[8],A2[9],A2[10]]
        A1xyz  = [A1[8],A1[9],A1[10]]
        adH1 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd)
        adH2 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,Dd1)
        adH3 = ZtoXYZ(A1xyz,A2xyz,A3xyz,r,Ad,-Dd1)
        nam=list(A3[2])
        if len(nam)<=2:
           name=nam[0]
        else:
           name = nam[0]+nam[1]

        name1 = "H1"+name
        name2 = "H2"+name
        name3 = "H3"+name
        pdb = add_to_pdb(pdb, A3, adH1, name1)
        pdb = add_to_pdb(pdb, A3, adH2, name2)
        pdb = add_to_pdb(pdb, A3, adH3, name3)
        return pdb

def add_to_pdb(pdb,old_atm,new_koord,new_nam):
## ad a new deuterium atom to behind an existing atom in the pdb list
        atom=0
        new_nam="H"
        for i in range(len(pdb)):
          for k in range(len(pdb[i])): 
#          print("ALT",alt)
#          atom=0
#          print(pdb[res_num][i][3])
#          print(pdb[res_num][i][2].strip())
           if pdb[i][k][2].strip()==old_atm[2].strip() and\
           pdb[i][k][3].strip()==old_atm[3].strip() and\
           pdb[i][k][5].strip()==old_atm[5].strip() and\
           pdb[i][k][6]==old_atm[6]:
              atom=pdb[i][k]
              indice = k+1
              res_num=i

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
        l_pdb.append("H")
        l_pdb.append(" ")
#       print("old_atom",old_atm)
#       print("atom",atom)
#       print("l_pdb",l_pdb)
#       print("???????????????????????????????????????????")
        pdb[res_num].insert(indice, l_pdb)
#       pdb[res_num].append(l_pdb)
        return pdb



def find_A1(pdb,A3,A2):
# find atom bond to A2 to define dihedral to A3
#       print("A3",A3)
#       print("A2",A2)
        bond_atoms=f_bond_atoms(pdb,A2)
#       print(bond_atoms)
        for i in range(len(bond_atoms)):
            if bond_atoms[i][2]!= A3[2]:
               A1=bond_atoms[i]


        return A1

def CAngle(x,y,z):
# calculate angle between a,b,c in degree
    x1=[0,0,0]
    x2=[0,0,0]
    for i in range(len(x)):
        x1[i]=y[i]-x[i]
        x2[i]=y[i]-z[i]
    Cangle=C2Angle(x1,x2)
    return Cangle
def C2Angle(x,y):
#Calculates the angle between x and y
#Answer in degrees
        #Calculate the angle between x and y
        rtodeg = 57.2957795
        C2angle=(ScalPr(x,y)/(len3dvec(x)*len3dvec(y)))
#       print("C2angle=",C2angle)
        if 0.999999999 <= C2angle <= 1.000000000001 :
           C2angle=0
        elif -0.999999999 >= C2angle >= -1.000000000001:
                C2angle=180
        else:
                C2angle=np.arccos(C2angle)*rtodeg
        return C2angle
def ScalPr(x,y):
#calculate the scalar product
        pro= x[0]*y[0]+x[1]*y[1]+x[2]*y[2]
        return pro


def Cross(x1,x2):
#Calculates the cross product x3 = x1 x x2
        x3 = [0,0,0]
        x3[0]=x1[1]*x2[2]-x2[1]*x1[2]
        x3[1]=x1[2]*x2[0]-x2[2]*x1[0]
        x3[2]=x1[0]*x2[1]-x2[0]*x1[1]
        return x3

def Normlz(xyz):
# Normalise xyz 
        temp = 1/len3dvec(xyz)
        for i in range(len(xyz)):
            xyz[i]=xyz[i]*temp

        return xyz
def ZeroVector():
## creates ZeroVector
        zerov = np.array([0,0,0])
        return zerov

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


def ZtoXYZ(axyz,bxyz,cxyz,R,Ad,Dd):
#defines coordinates of 4th atom fro coordinates of 3 atoms distance
#angel and dihedral 
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






def f_bond_atoms(pdb,atom_qm):
## find bondet atoms
       bond_atom=[]
       xyz=[atom_qm[8],atom_qm[9],atom_qm[10]]
       bond_atom=[]
#      print("xyz&&&",xyz)
       for i in range(len(pdb)):
          for k in range(len(pdb[i])):
              p2=[pdb[i][k][8],pdb[i][k][9],pdb[i][k][10]]
              dist=CDist2(xyz,p2)
       #      print(dist)
              if 0.5< dist< 2:
#                print(dist)
                 bond_atom.append(pdb[i][k])
       return bond_atom



def def_n_bonds(xyz,ele,bond_atoms):
# defines how many atoms are bond to the central atom with bond order
        nr_bonds=0
#       print("xyz!!!",xyz)
        for i in range(len(bond_atoms)):
            p2=[bond_atoms[i][8],bond_atoms[i][9],bond_atoms[i][10]]
            dist=CDist2(xyz,p2)
            l_ele=[ele,bond_atoms[i][14]]
#           print("dist",dist)
            if l_ele[0]=="C":
               if l_ele[1]=="C":
                  if 0.5 < dist < 1.25:
                     nr_bonds=nr_bonds+3
                  if 1.25<= dist <=1.37:
                     nr_bonds=nr_bonds+2
                  if 1.38<= dist<= 1.47:
                     nr_bonds=nr_bonds+1.5
                  if 1.48<= dist:
                     nr_bonds=nr_bonds+1
               elif l_ele[1]=="O":
                  if 0.5<= dist <=1.32:
                     nr_bonds=nr_bonds+2
                  if 1.33<= dist:
                     nr_bonds=nr_bonds+1
               elif l_ele[1]=="N":
                  if 0.5<= dist <=1.39:
                     nr_bonds=nr_bonds+1.5
                  if 1.4<= dist:
                     nr_bonds=nr_bonds+1
               elif l_ele[1]=="H":
                  if 0.5<= dist <=1.39:
                     nr_bonds=nr_bonds+1
               else:
                  if 0.5<= dist <=2.5:
                     nr_bonds=nr_bonds+1

            elif l_ele[0]=="N":
               if l_ele[1]=="C":
                  if 0.5<= dist <=1.39:
                     nr_bonds=nr_bonds+1
                  elif 1.4<= dist:
                     nr_bonds=nr_bonds+1
               elif l_ele[1]=="H":
                  if 0.5<= dist <=1.39:
                     nr_bonds=nr_bonds+1
               else:
                  if 0.5<= dist <=2.5:
                     nr_bonds=nr_bonds+1

            elif l_ele[0]=="O":
               if l_ele[1]=="C":
                  if 0.5<= dist <=1.32:
                     nr_bonds=nr_bonds+2
                  elif 1.32 < dist:
                     nr_bonds=nr_bonds+1
               elif l_ele[1]=="H":
                  if 0.5<= dist <=1.29:
                     nr_bonds=nr_bonds+1
               else:
                  if 1.55< dist <=2.5:
                     nr_bonds=nr_bonds+1
                  if 0.5<= dist <= 1.55:
                     nr_bonds=nr_bonds+2

        return nr_bonds



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


def renumber_pdb(pdb):
#renumbers pdb
        num=0
        for i in range(len(pdb)):
            for j in range(len(pdb[i])):
                num=num+1
                pdb[i][j][1]=num

        return pdb



def write_pdb(pdb,head, output):
## wite out PDB file
#       res = open(output, "w")
        res= open("comqum.pdb","w")
        atm_nr=0
        for i in range(len(head)):
            res.write(head[i])
#           res.write("\n")
        for i in range(len(pdb)):
            for k in range(len(pdb[i])):
               atm_nr = atm_nr + 1
               string = str('{:6}'.format(pdb[i][k][0]))
               string = string + str('{:5.0f}'.format(atm_nr))
               string = string + "  "
               string = string + str('{:3s}'.format(str(pdb[i][k][2])))
               string = string + str('{:1}'.format(pdb[i][k][3]))
               string = string + str('{:3}'.format(pdb[i][k][4]))
               string = string + str('{:>2}'.format(pdb[i][k][5]))
               string = string + str('{:4}'.format(pdb[i][k][6]))
               string = string + str('{:1}'.format(pdb[i][k][7]))
               string = string + "   "
               string = string + str('{:8.3f}'.format(pdb[i][k][8]))
               string = string + str('{:8.3f}'.format(pdb[i][k][9]))
               string = string + str('{:8.3f}'.format(pdb[i][k][10]))
               string = string + str('{:6.2f}'.format(pdb[i][k][11]))
               string = string + str('{:6.2f}'.format(pdb[i][k][12]))
               string = string + str('{:>7}'.format(pdb[i][k][13]))
               string = string + str('{:>5}'.format(pdb[i][k][14]))


               string = string + "\n"
               res.write(string)
        res.write("END")
#######################################################################################
#processing data
print(input_pdb)
print(qm)
qm_atom, qm_resi, qm_anam, qm_anum, qm_rnum= read_qm(qm)
print("Input")
print("ATOM",qm_atom)
print("RESI",qm_resi)
print("ANUM",qm_anum)
print("ANAM",qm_anam)
print("RNUM", qm_rnum)

org_pdb, work_pdb, head = read_pdb(input_pdb)


pdb=atm_to_res(work_pdb)

atm_nr, atom_qm = assign_atm_nr(qm_atom, qm_resi, qm_anum, qm_anam, qm_rnum ,pdb) 


pdb=proto_qm(pdb, atom_qm, f_p)
pdb =renumber_pdb(pdb)
atm_nr, atom_qm = assign_atm_nr(qm_atom, qm_resi, qm_anum, qm_anam, qm_rnum ,pdb) 
#print("atm_nr",atm_nr)
#for i in range(len(atom_qm)):
#    print(atom_qm[i])

trunc_nr, trunc_atm = find_truncation(atm_nr, atom_qm, pdb)
f_all, f_qm, f_trunc = sum_list(atm_nr, atom_qm, trunc_nr, trunc_atm)
print(" ")
print("f_all",f_all)
#print("f_qm",f_qm)
print("f_trunc",f_trunc)
write_syst1(f_all,f_trunc)
write_restart_anam(f_qm,work_pdb)
write_rename(f_all,f_trunc,work_pdb)
if f_p== True:
        write_pdb(pdb,head,input_pdb)


time_ende = time.time()
print("program ends normally after "'{:5.3f}s'.format(time_ende-time_start)," or ", '{:5.2f}min'.format((time_ende-time_start)/60))

