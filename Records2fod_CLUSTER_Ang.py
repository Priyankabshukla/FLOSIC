### Records2fod_movie_CLUSTER.xyz (Angstrom Units)

import numpy as np
import os
import sys
root=os.getcwd()
path=root+'/'+'records'
print(path)
outpath=root+'/'+'fod_movie.xyz'
mol_path=root+'/'+'XMOL.xyz'

infile = open(path,'r')
outfile=open(outpath,'w+')
molfile=open(mol_path,'r')
found_ene=0
found_spins=0
atm=[]
frame=[]
count=0
natoms=0 #number of fods in records file
ene=[]
spin_updown=[]
mol_xyz=[] #pbs13
nmols=0 #number of atoms in XMOL.xyz file #pbs13
natomsplusnmols=0 #pbs13

######### pbs13 #####
for i, line in enumerate(molfile):
    if len(line.split())==1:
        nmols=line.split()[0]
        
    if len(line.split())==4:
        mol_xyz.append([str(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
#print(mol_xyz)
######### pbs13 #####

for i,line in enumerate(infile):
    if atm and count==natoms:
#        print(atm)
        
        frame.append(atm)
        
        count=0
        found_ene=0
        found_spins=0
        atm=[]
    
    
    if found_ene==1 and found_spins==1 and len(line.split())==5:
        
        count+=1
        atm.append([float(line.split()[0])*0.529177,float(line.split()[1])*0.529177,float(line.split()[2])*0.529177])
        
    if len(line.split())==1:
        found_ene=1
        found_ene_i=i
        temp_ene=line.split()
        
    if len(line.split())==4:
        spin_updown.append(line.split())

        
        if i-found_ene_i==1:
            found_spins=1
            ene.append(temp_ene)

        natoms=int(line.split()[0])+int(line.split()[1])
        natomsplusnmols=natoms+int(nmols) #pbs13
        

for i in range(len(frame)):
    fr_i=frame[i]
#     print(fr_i)

    outfile.write(str(natomsplusnmols)+'\n') #pbs13
    outfile.write(str(spin_updown[i][0])+' '+str(spin_updown[i][1])  +'\n')

    for j in range(len(fr_i)):
        if j<int(spin_updown[i][0]):
            outfile.write('He'+' '+str(fr_i[j][0])+' '+str(fr_i[j][1])+' '+ str(fr_i[j][2])+'\n')  
        else:
            outfile.write('X'+' '+str(fr_i[j][0])+' '+str(fr_i[j][1])+' '+ str(fr_i[j][2])+'\n')
    for k in range(len(mol_xyz)):
#         print(mol_xyz[k][0])
#         print(i)
        outfile.write(mol_xyz[k][0]+' '+ str(mol_xyz[k][1])+' '+str(mol_xyz[k][2])+' '+str(mol_xyz[k][3])+'\n')

    
outfile.close()

