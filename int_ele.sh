#!/usr/bin/env python3
#programm to integrate electron densety around sulfulre atoms
# run ./int_s.sh --pdb=pdb --map=ccp4 --ele=ele --out=outfile
# --fout=formattet output

## new way of readin in grid files 
##https://gemmi.readthedocs.io/en/latest/grid.html


import time
## set timer
time_start = time.time()
import argparse, sys
##load mandatory modules
import numpy as np
import scipy as sc 
import gemmi

box_length = 5.2 #box oround atom
fin_res = 0.051 # minumum resolution the grit gets interpolate to
mesh = 0.05


parser=argparse.ArgumentParser()
parser.add_argument('--pdb', help='input pdb file')
parser.add_argument('--map', help='input ccp4 map')
parser.add_argument('--ele', help='element to integrat')
parser.add_argument('--out', help='output file')
parser.add_argument('--fout', help='list for relabeling ')






args=parser.parse_args()


if args.pdb != None:
   input_pdb=args.pdb
   print("input pdb file is:", input_pdb)
else:
         input_pdb="pdb"
         print("input pdb file is:", input_pdb)

if args.map != None:
           input_grid=args.map
           print("read from ", input_grid)
else:
           input_grid="map.dx"
           print("read from defalt ", input_grid)

if args.ele != None:
           ele=args.ele
           print("element ",ele)
else:
           ele="S"
           print("default ele ", ele)


if args.out != None:
           out_file=args.out
           print("output_file ",out_file)
else:
           out_file="results"
           print("default out_file ", out_file)

if args.fout != None:
           f_out=args.fout
           f_out_f=True
           print("read for formatted output ",f_out)
else:
           f_out="S"
           f_out_f=False
           print("no fromatet output ")














#########################################################################################
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








def str_to_float(list_in):
## converts list of str to list of float
        for i in range(len(list_in)):
                list_in[i]=float(list_in[i])
        return(list_in)



def def_box(box_length, delta):
#define box with dimensions of box lenght
        box = []
        for i in range(len(delta)):
            box.append(box_length * (1/delta[i]))
        ## create diameter in int
        for i in range(len(box)):
            box[i]= int(box[i])+1
        return box


def cov_rad():
#list of Covalent radii in  from  Cambridge Structural Database
        cov_rad=[["H",0.31],["D",0.31],["C",0.76],["N",0.71],["O",0.66],["S",1.05],["Fe",1.52],["Mo",1.54]]
        cov_rad.append(["FE",1.52])
        cov_rad.append(["MO",1.54])

        return cov_rad

def def_r_int(ele):
## creates integration radios out of the covalent radius
    cov=cov_rad()
    for i in range(len(cov)):
        if cov[i][0]==ele:
           r_int=cov[i][1]
    return r_int

def create_box(atom, box, origen, delta, coord):
## creates box orund atom
#       print("define box around oxygen atoms")

        pos=[0,0,0]
        for i in range(3):
            pos[i]=(coord[i]-origen[i])/delta[i]

        box_e = []
        offset_e_box=[]
        for i in range(len(box)):
            offset_e_box.append(((pos[i]-int(pos[i]-box[i]))-int(pos[i]-int(pos[i]-box[i])))*delta[i])
            box_e.append(int(pos[i]-box[i]))
            box_e.append(int(pos[i]+box[i]))
        return box_e, pos, offset_e_box


def  fill_box(box_e, data_mat, grid):
        ## (x,y,z)matrix arount atom at center
        ## x slow, y med, z fast
        grid_box = []
#       grid_box.append([box_o[0]])
        zeros=np.zeros((box_e[1]-box_e[0],box_e[3]-box_e[2],box_e[5]-box_e[4]))
#       grid_box.append(zeros)
        grid_box=zeros        
        cx = -1
        print("box_e",box_e)
        for x in range(box_e[0],box_e[1]):
                cx = cx + 1
                cy = -1
                for y in range(box_e[2],box_e[3]):
                        cy = cy +1
                        cz = -1
                        for z in range(box_e[4],box_e[5]):
                                cz=cz+1
                                if x < 0:
                                        x = grid[0]+x
                                if y < 0:
                                        y = grid[1]+y
                                if z < 0:
                                        z = grid[2]+z
                                grid_box[cx][cy][cz]=data_mat[x][y][z]
#       print("grid_box",grid_box)
#       print("box has ben filled with actual data")
        return grid_box



def interpol_to(grid, box_length, final_rel):
## interpolates to a final_rel resoluteion 

        res = dim_box(box_length,grid)
        while res[0] >= final_rel:
             grid = InterPol3D(grid)
             res = dim_box(box_length,grid)

        return grid, res[0]

def InterPol3D(matrix):
#interpolateds 3d grit
        R_mat = np.zeros((2*len(matrix)-1,2*len(matrix[0])-1,2*len(matrix[0][0])-1))

        for a in range(0,len(R_mat),2):
          for b in range(0,len(R_mat[0]),2):
            for c in range(0,len(R_mat[0][0]),2):
              R_mat[a][b][c] = matrix[int(a/2)][int(b/2)][int(c/2)]

        for a in range(0,len(R_mat),2):
          for b in range(0,len(R_mat[0]),2):
            for c in range(1,len(R_mat[0][0]),2):
              R_mat[a][b][c] = (R_mat[a][b][c-1]+R_mat[a][b][c+1])/2

        for a in range(0,len(R_mat),2):
          for b in range(1,len(R_mat[0]),2):
            for c in range(0,len(R_mat[0][0]),2):
                     R_mat[a][b][c] = (R_mat[a][b-1][c]+R_mat[a][b+1][c])/2

        for a in range(0,len(R_mat),2):
          for b in range(1,len(R_mat[0]),2):
            for c in range(1,len(R_mat[0][0]),2):
                     R_mat[a][b][c] = (R_mat[a][b-1][c-1]+R_mat[a][b+1][c-1]+R_mat[a][b-1][c+1]+R_mat[a][b+1][c+1])/4

        for a in range(1,len(R_mat),2):
          for b in range(0,len(R_mat[0]),2):
            for c in range(0,len(R_mat[0][0]),2):
                R_mat[a][b][c] = (R_mat[a-1][b][c]+R_mat[a+1][b][c])/2

        for a in range(1,len(R_mat),2):
          for b in range(0,len(R_mat[0]),2):
            for c in range(1,len(R_mat[0][0]),2):
                   R_mat[a][b][c] = (R_mat[a-1][b][c-1]+R_mat[a-1][b][c-1]+R_mat[a+1][b][c+1]+R_mat[a+1][b][c+1])/4

        for a in range(1,len(R_mat),2):
          for b in range(1,len(R_mat[0]),2):
            for c in range(0,len(R_mat[0][0]),2):
                     R_mat[a][b][c] = (R_mat[a-1][b-1][c]+R_mat[a-1][b+1][c]+R_mat[a+1][b-1][c]+R_mat[a+1][b+1][c])/4

        for a in range(1,len(R_mat),2):
          for b in range(1,len(R_mat[0]),2):
            for c in range(1,len(R_mat[0][0]),2):
                     R_mat[a][b][c] = (R_mat[a-1][b-1][c-1]+R_mat[a-1][b-1][c+1]+R_mat[a-1][b+1][c-1]+R_mat[a-1][b+1][c+1]+R_mat[a+1][b-1][c-1]+R_mat[a+1][b-1][c+1]+R_mat[a+1][b+1][c-1]+R_mat[a+1][b+1][c+1])/8

        return R_mat

def dim_box(box_length,grid):
#calculates the voxel lenght af a grid
        dim= [box_length/len(grid),box_length/len(grid[0]),box_length/len(grid[0][0])]
        return dim

def get_pos(grid_box, offset_e_box,box_length):
## gets position of element in box
    box=[len(grid_box),len(grid_box[0]),len(grid_box[0][0])]

    dim = dim_box(box_length, grid_box)
    pos_e_in_box=[0,0,0]
    for i in range(len(pos_e_in_box)):
                pos_e_in_box[i]=box[i]/2+(offset_e_box[i]/dim[i])

    return pos_e_in_box

def len3dvec(vec):
## calculates lengh of a 3D vecor
## input as list
        a = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        return a


def int_ball(grid_box,box_length,r_s,point,dim):
## integrades in radius r_s around point
        dim=dim_box(box_length, grid_box)
        print("box_length",box_length)
        print("dim",dim)
        print("r_s",r_s)
        v_len= vec_dim(r_s, dim)
        print("r_s",r_s)
        print("v_len",v_len)
        print("dim",dim)
        r_rad=[int(r_s/dim[0]) ,int(r_s/dim[1]), int(r_s/dim[2])]
        print(point,r_rad)
        p_int_low=[point[0]-r_rad[0],point[1]-r_rad[1],point[2]-r_rad[2]]
        p_int_high=[point[0]+r_rad[0],point[1]+r_rad[1],point[2]+r_rad[2]]

        int_high=[0,0,0]
        int_low=[0,0,0]
        for r in range(len(p_int_low)):
                if (p_int_low[r]< p_int_high[r]):
                        int_high[r] = p_int_high[r]
                        int_low[r]  = p_int_low[r]
                else:
                        int_high[r] = p_int_low[r]
                        int_low[r]  = p_int_high[r]

        for i in range(len(int_high)):
                int_low[i]=int(int_low[i])-2
                int_high[i]=int(int_high[i])+2




        p_to_integrate=[]
        p1=[0,0,0]
        for x in range(int_low[0], int_high[0]):
             for y in range(int_low[1], int_high[1]):
                 for z in range(int_low[2], int_high[2]):
                        p0=[x*dim[0],y*dim[1],z*dim[2]]
                        p1=[point[0]*dim[0],point[1]*dim[1],point[2]*dim[2]]
                        p0p1=[p0[0]-p1[0],p0[1]-p1[1],p0[2]-p1[2]]
#                       print(p0p1)
                        dist=len3dvec(p0p1)
#                       print(dist)
                        if dist <= r_s:
                                p_to_integrate.append([x,y,z])
#       print(len(p_to_integrate))
        counter=0
        for i in range(len(p_to_integrate)):
                counter=counter+grid_box[p_to_integrate[i][0]][p_to_integrate[i][1]][p_to_integrate[i][2]]
        inte_val=counter/len(p_to_integrate)
#       print("integration",inte_val)
        return inte_val




def vec_dim(length,dim):
## maximum of points fo a vector
        v_len = [length/dim[0], length/dim[1], length/dim[2]]
        return v_len





def fill_list(e_list,atom,int_ele):
## fill list with  atom number, atom name,chain, res_num, res, integration

        element=[]
        element.extend([atom[1],atom[2], atom[5], atom[6],atom[4]])
        element.extend([int_ele])
        e_list.append(element)

        return e_list


def average(i_list):
## calculate average of list
#       print(sum(i_list))
#       print(len(i_list))
        ave=sum(i_list) / len(i_list)
        return ave

def std(i_list,mid):
## calculate standard deviation
        l_sqr=[]
        for i in range(len(i_list)):
            num=(i_list[i]-mid)**2
            l_sqr.append(num)
        s2=sum(l_sqr)/(len(l_sqr)-1)
        s=np.sqrt(s2)
        return s



def ad_std(e_list, v_mid, v_std):
## add how manny stads it deviates from the average
        for i in range(len(e_list)):
            dif=(e_list[i][6]-v_mid)/v_std
            e_list[i].append(dif)

        return e_list


def stat(list):
## does some statistic om integration
        l_int=[]
#       print(list)
        for i in range(len(list)):
            l_int.append(list[i][6])
#       print(l_int)
        v_mid=average(l_int)
        v_std =  std(l_int,v_mid)
#       print("v_std",v_std) 
        list=ad_std(list, v_mid, v_std)

        return list , v_mid, v_std



def write_results(e_list,out_file,ave,std):
## writes results in formatet output file
        res = open(out_file, "w")
        head=" atom number| atom name | chain| res number and name| integradted value | difernence to average in sigma|"
        head = head +   str('{:3.3f}'.format(ave))+"  "+str('{:3.3f}'.format(std))+"\n"
        res.write(head)
        for i in range(len(e_list)):
           string = str('{:8}'.format(e_list[i][0]))+"  "
           string = string + str('{:5}'.format(e_list[i][1]))
           string = string + str('{:5}'.format(e_list[i][4]))
           string = string + str('{:3}'.format(e_list[i][3]))+"-"
           string = string + str('{:4}'.format(e_list[i][5]))
           string = string + str('{:9.4f}'.format(e_list[i][6]))
           string = string + str('{:9.4f}'.format(e_list[i][7]))
           string = string + "\n"
           res.write(string)

def write_f_results(e_list,out_file,ave,std,f_out):
#writes formatet output 
        f_o = open(f_out, "r")
        l_f_out=[]
        out_file= out_file + "_F"
        res = open(out_file, "w")
        head="atom_name   chain  integradted_value   difernence_to_average_in_sigma "
        head = head + " ave "+  str('{:3.3f}'.format(ave))+" sig  "+str('{:3.3f}'.format(std))+"\n"
        res.write(head)
        print("!!!!!!!!!!!!!",f_out)
        for line in f_o:
            l_f_out.append(line.split())
        for i in range(len(e_list)):
            for k in range(len(l_f_out)):
              if len(l_f_out[k])>=3:
                if int(e_list[i][5]) == int(l_f_out[k][2]):
                 if e_list[i][4] == l_f_out[k][1]:
                   string= str(l_f_out[k][0])+" "
                   string= string+ str(l_f_out[k][1])
                   string = string + str('{:9.4f}'.format(e_list[i][6]))
                   string = string + str('{:9.4f}'.format(e_list[i][7]))
                   string = string + "\n"
                   res.write(string)




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


def read_grid(input_grid):
#read ccp4 map
        print(input_grid)
        
        map= gemmi.read_ccp4_map(input_grid)
        print(gemmi.InfoMap())
#       print("head",map.header_str)
        print(map.__repr__())
#       print(help(map))
        print(" map is readet")
        return map


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
        print("#####################################")
        print("vol",vol)
        print("int_val",int_val)
        print("¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤")
        return int_val




def fill_l_int(l_int,atom,int_val):
#fill list with all important information 
#[atm_nr, atm_nam, ele, chain, res,integratet value]
        entry=[]
        entry.append(atom[1])
        entry.append(atom[2])
        entry.append(atom[14])
        entry.append(atom[4])
        entry.append(atom[5])
        entry.append(atom[6])
        entry.append(int_val)
        l_int.append(entry)
        return l_int


##################################################
#input_pdb, input_grid, out_file , ele, f_out, f_out_f = read_arg()

pdb_read, work_pdb = read_pdb(input_pdb)

map= read_grid(input_grid)
map.setup()
print(map.grid)
print("map",map)

arr = np.array(map.grid, copy=False)
#print(arr)

print("!!!!!!!!!", map.grid.get_value(1,1,1).conjugate())
print("################",map.grid.interpolate_value(gemmi.Position(1.5,1.5,1.5)))


r_int=def_r_int(ele) 
r_int=r_int/1
#r_int=0.5
l_int=[]
for atom in range(len(work_pdb)):
#   print(ele)
#   print(work_pdb[atom])
    if work_pdb[atom][14]== ele :
       print(work_pdb[atom])
       coord=[work_pdb[atom][8],work_pdb[atom][9],work_pdb[atom][10]]
       print("value at ",ele," position",map.grid.interpolate_value(gemmi.Position(coord[0],coord[1],coord[2])))
       int_val = int_around_point(coord,map,r_int, mesh)
       print("int_val",int_val)
       l_int=fill_l_int(l_int,work_pdb[atom],int_val)
e_list,ave,std=stat(l_int)

print("ave",ave,"std",std)
for i in range(len(e_list)):
       print(e_list[i])
write_results(e_list,out_file,ave,std)
if f_out_f==True:
   write_f_results(e_list,out_file,ave,std,f_out)
time_ende = time.time()
print("program ends normally after "'{:5.3f}s'.format(time_ende-time_start),
" or ", '{:5.2f}min'.format((time_ende-time_start)/60))





