### Records2fod_movie.xyz (Angstrom Units)

import numpy as np
import os
import sys
root=os.getcwd()
path=root+'/'+'records'
print(path)
outpath=root+'/'+'fod_movie.xyz'

infile = open(path,'r')
outfile=open(outpath,'w+')
found_ene=0
found_spins=0
atm=[]
frame=[]
count=0
natoms=0
ene=[]
spin_updown=[]
for i,line in enumerate(infile):
    if atm and count==natoms:
        
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
print(spin_updown[0][0])
for i in range(len(frame)):
    fr_i=frame[i]

    outfile.write(str(natoms)+'\n')
    outfile.write(str(spin_updown[i][0])+' '+str(spin_updown[i][1])  +'\n')

    for j in range(len(fr_i)):
        if j<int(spin_updown[i][0]):
            outfile.write('H'+' '+str(fr_i[j][0])+' '+str(fr_i[j][1])+' '+ str(fr_i[j][2])+'\n')  
        else:
            outfile.write('X'+' '+str(fr_i[j][0])+' '+str(fr_i[j][1])+' '+ str(fr_i[j][2])+'\n')

    
outfile.close()

