import numpy as np
import os
import sys

ene_input = sys.argv[1] 

print(ene_input)
root=os.getcwd()
path=root+'/'+'records'
print(path)
outpath=root+'/'+'movie_records'

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
        atm.append([line.split()[0],line.split()[1],line.split()[2]])
        
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
 
    
lsic_path=root+'/'+'FRMIDT'  
lsic_frmidt=open(f'{lsic_path}','w+')


for i in range(len(frame)):
    fr_i=frame[i]

    outfile.write(str(ene[i])+'\n')
    outfile.write(str(spin_updown[i][0])+' '+str(spin_updown[i][1])  +'\n')

    for j in range(len(fr_i)):
        outfile.write(str(fr_i[j][0])+' '+str(fr_i[j][1])+' '+ str(fr_i[j][2])+'\n')
    
    if round(float(ene[i][0]),12)== float(ene_input):  #input pzsic energy from fande.out to create FRMORB file
        print('yes')

        fr_i=frame[i]
        lsic_frmidt.write(str(spin_updown[i][0])+' '+str(spin_updown[i][1])  +'\n')
        for j in range(len(fr_i)):
            lsic_frmidt.write(str(fr_i[j][0])+' '+str(fr_i[j][1])+' '+ str(fr_i[j][2])+'\n')
        
        
        
        
    
    
outfile.close()
lsic_frmidt.close()

        
