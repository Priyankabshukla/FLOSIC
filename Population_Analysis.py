import numpy as np
import os
import matplotlib.pyplot as plt
import numpy as np
import ase
import matplotlib.pyplot as plt
from ase.io import read, write
import math
from scipy.spatial.distance import pdist, squareform
from ase.neighborlist import NeighborList
import os
import pandas as pd

# Lowdin Charge Analysis from the FLOSIC output file
def Lowdin_charge(direc):
    f=open(f'{direc}/print.1')
    found=0
    stop=0
    count=0
    charge=[]
    for lines in f:
        if stop==1:
            break
        if found ==1:
            count+=1
            if len(lines.split())==4 or len(lines.split())==5:
                if count>=6:
                    charge.append(float(lines.split()[3]))
            if len(lines.split())==1 and count>6:
                stop=1

        if len(lines.split())==3:
            if lines.split()[0]=='LOWDIN':
                found=1
    return charge

# Lowdin Spin analysis and spin expectation value
def Lowdin_spin(direc):
    f=open(f'{direc}/print.1')
    found=0
    stop=0
    count=0
    charge=[]
    exp=[]
    for lines in f:
        if stop==1:
            break
        if found ==1:
            count+=1
            if len(lines.split()) ==3 and lines.split()[0]=='<S**2>':
                exp.append(float(lines.split()[2]))
                print(lines.split()[2])
            if len(lines.split())==4 or len(lines.split())==5:
                if count>=6:
                    charge.append(float(lines.split()[4]))
            if len(lines.split())==1 and count>6:
                stop=1

        if len(lines.split())==3:
            if lines.split()[0]=='LOWDIN':
                found=1
    return charge, exp


root = '/bgfs/kjohnson/pbs13/FLOSIC/Projects/Zeolites_SIC/Cu-SSZ-13/Final/fod_opt/all_result'

print("root",root)
for entry in os.listdir(root):
    os.chdir(root)
    if not entry.endswith('.py') and entry!='another_guess' and entry!='NRLMOL_INPUT.DAT' and entry!='nrlmolDft.opa.batch' and entry!='RUNS' and entry!='O2':
        path = os.path.join(root, entry)
        print(path)
        char_lda,exp_lda=Lowdin_spin(f'{path}/lda_pop/')
        atom_label=np.loadtxt(f'{path}/lda_pop/XMOL.xyz',skiprows=2, usecols=(0),dtype='str')
        charge_pzsic,exp_pzsic=Lowdin_spin(f'{path}/pop/')
        print('LDA <S**2>: ', round(float(exp_lda[0]),3))
        print('PZSIC <S**2>: ', round(float(exp_pzsic[0]),3))

        x = np.arange(len(char_lda))
        my_xticks = list(atom_label)
        plt.figure(figsize=(8,6))

        plt.plot(x,char_lda,'ro-',label=f'{entry} LDA',linewidth='3',markersize='10')
        plt.plot(x,charge_pzsic,'bo-',label=f'{entry} PZSIC',linewidth='3',markersize='10')
        plt.ylim(-0.4,0.4)
        plt.xticks(x, my_xticks,size=20)
        plt.yticks(size=20)
        plt.legend(frameon = False,prop={'size': 20})
        plt.xlabel('Atom Label',size=20)
        plt.ylabel('Spin',size=20)
        plt.show()




                

