import numpy as np
import matplotlib.pyplot as plt
import os
import sys

print_input = sys.argv[1] # give name of the print file to know it's evalues

root=os.getcwd()
path=root+'/'+f'{print_input}'

file=open(f'{path}','r')
lines=file.readlines()

found_lowden=1
found_lowsic=1
found_evaloc=1
evalues=[]
for line in lines:
    if len(line.split())==1 and line.split()[0]=='EVALOCC:' :
            found_evaloc=0
    if len(line.split())==3 and line.split()[2]=='LOWSIC' :
            found_lowsic=0
    if found_lowden==0 and found_lowsic==1 and found_evaloc==1:
        
    #    print(line.split())
        evalues.append(line.split())

   
    if len(line.split())==3 and line.split()[0]=='LOWDEN' and line.split()[1]=='OVERLAP':
        found_lowden=0
        found_lowsic=1
        found_evaloc=1
        
        
        
        
       
flat_list = []
 
for sublist in evalues:
    for item in sublist:
        flat_list.append(float(item))
        
if min(flat_list)<0.1:
    print('Bad eigenvalues: need another guess')
    print(min(flat_list))
else:
   print('Good')
