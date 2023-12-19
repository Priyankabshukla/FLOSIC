import numpy as np
import ase
import matplotlib.pyplot as plt
from ase.io import read, write
import math
from scipy.spatial.distance import pdist, squareform
from ase.neighborlist import NeighborList
import os
import pandas as pd
import itertools


def e_sic(path):
    FLO,Charge,Asymp,e_coul,e_sic,exlocal,eclocal,exnon,ecnon = np.loadtxt(path,skiprows=1,unpack=True)
    plt.plot(e_sic,'ro')
    plt.show()

    plt.plot(e_sic[e_sic>-0.02],'o')
    plt.show()
    print(np.where(e_sic>-0.02))
    print(e_sic[e_sic>-0.02])

    print(np.where(e_sic>-0.01))
    print(e_sic[e_sic>-0.01])
    return np.where(e_sic>-0.01),len(e_sic)


def Print2Smallest(a,n,ind_list,fod_atom):
    VALUE = max(a) + 1
    idx_min0 = a.index(min(a))
    min0 = a[idx_min0]
    a[idx_min0] = VALUE
    idx_min1 = a.index(min(a))
    min1 = a[idx_min1]
    a[idx_min1] = VALUE
    a[idx_min0] = min0
    a[idx_min1] = min1
    
    
    return ([min0,min1],[idx_min0,idx_min1],[ind_list[idx_min0][0],ind_list[idx_min1][0]],
         [fod_atom.get_chemical_symbols()[ind_list[idx_min0][0]],
         fod_atom.get_chemical_symbols()[ind_list[idx_min1][0]]]
         )
    

def find_core(direc):
    Cu_list=[]
    fod_atom=read(f'{direc}')
    Distance=fod_atom.get_all_distances()
    # Find the indices of all atoms that match the desired element (carbon)
    H_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'H']
    O_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'O']
    Cu_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'Cu']
    N_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'N']
    Al_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'Al']
    He_indices=[i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'He']
    X_indices=[i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'X']
    cutoff = 1.2
    nl = NeighborList([cutoff/2]*len(fod_atom), self_interaction=False, bothways=True)
    nl.update(fod_atom)
    
    print(f'NOOH : {np.count_nonzero(Distance[Cu_indices]<0.6)-1}') # This counts the number of fods inside Cu
    Cu_list.append(np.where(Distance[Cu_indices[0]]<0.6))
    print(list(Cu_list[0][0]))
    print("Cu_inner",np.where(Distance[Cu_indices[0]]<0.6))
    print("Al_inner:",np.where(Distance[Al_indices[0]]<0.4))
    Al_list=[]
    Al_list.append(np.where(Distance[Al_indices[0]]<0.4))
    
    print(O_indices)
    O_list=[]
    O_list.append(np.where(Distance[O_indices]<0.4))
    print("*****",list(O_list[0][1]))
    print("H_indices",H_indices)
    H_list=[]
    H_list.append(np.where(Distance[H_indices]<0.4))
    print("*****",list(H_list[0][1]))
    N_list=[]
    N_list.append(np.where(Distance[N_indices]<0.4))
    print("*****",list(N_list[0][1]))
    
    for k in range(len(O_indices)):
        print(np.where(Distance[O_indices[k]]<0.4))
        
    for m in range(len(H_indices)):
        print(np.where(Distance[H_indices[m]]<0.4))
    for n in range(len(N_indices)):
        print(np.where(Distance[N_indices[n]]<0.4) )
        
        
    big_list=[]
    big_list.append([list(Cu_list[0][0]),list(Al_list[0][0]),list(O_list[0][1]),list(H_list[0][1]),list(N_list[0][1])])
    
    print("Final:",list(big_list))
        
        
    return big_list
       

    
fod_path=f'{path}'+'fod_xmol.xyz'
core_fod=find_core(fod_path)

a=np.array(core_fod,dtype=object)
b=list(itertools.chain(*core_fod))
c=list(itertools.chain(*b))
FLO,Charge,Asymp,e_coul,e_sic,exlocal,eclocal,exnon,ecnon = np.loadtxt(path+"SIC_ENERGY",skiprows=1,unpack=True)
e_sic
for i in j2:
    print(e_sic[i])  #SIC_ENERGY of all core shell fods
    
def fod_nearest_atom(path,orb_index,len_sic):
    Cu_list=[]
    Al_list=[]
    fod_atom=read(f'{path}')
    Distance=fod_atom.get_all_distances()
    # Find the indices of all atoms that match the desired element (carbon)
    H_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'H']
    O_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'O']
    Cu_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'Cu']
    N_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'N']
    Al_indices = [i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'Al']
    He_indices=[i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'He']
    X_indices=[i for i, symbol in enumerate(fod_atom.get_chemical_symbols()) if symbol == 'X']
    cutoff = 1.2
    nl = NeighborList([cutoff/2]*len(fod_atom), self_interaction=False, bothways=True)
    nl.update(fod_atom)

    # Create the connectivity matrix
    connectivity = np.zeros((len(fod_atom), len(fod_atom)), dtype=int)
    for i in range(len(fod_atom)):
        indices, offsets = nl.get_neighbors(i)
        for j in indices:
            if i != j:
                connectivity[i,j] = 1
    
    indices=np.argwhere((connectivity[orb_index]==1) & (np.arange(len(fod_atom)) > len_sic))
    
    ind_list=list(indices)
    ind_list
    ind_dis=[]
    for i in ind_list:
        ind_dis.append(Distance[orb_index,i[0]])
     
    a,b,atom_ind,atom_sym=Print2Smallest(ind_dis,len(ind_dis),ind_list,fod_atom)
    atom_ind
    j_m= Distance[atom_ind[0],atom_ind[1]]
    k_j=Distance[orb_index,atom_ind[0]]
    k_m=Distance[orb_index,atom_ind[1]]
    O_list=[]
    H_list=[]
    N_list=[]

    if j_m > k_j and j_m>k_m:
        if k_j < k_m:
#             print(f'{orb_index} belongs to {atom_sym[0]}')
            print(f'{orb_index} belongs to {atom_ind[0]}')
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[0])
#                 print(Cu_list)
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[0])
#                 print(Al_list)
            if atom_sym[0]=='O':
                O_list.append(atom_ind[0])
#                 print(O_list)
            if atom_sym[0]=='H':
                H_list.append(atom_ind[0])
#                 print(H_list)
                
            if atom_sym[0]=='N':
                N_list.append(atom_ind[0])
#                 print(N_list)
            
        else:
#             print(f'{orb_index} belongs to {atom_sym[1]}')
            print(f'{orb_index} belongs to {atom_ind[1]}')
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[1])
#                 print(Cu_list)
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[1])
#                 print(Al_list)
            if atom_sym[0]=='O':
                O_list.append(atom_ind[1])
#                 print(O_list)
            if atom_sym[0]=='H':
                H_list.append(atom_ind[1])
#                 print(H_list)
            if atom_sym[0]=='N':
                N_list.append(atom_ind[1])
#                 print(N_list)
                
    elif j_m<k_j or j_m<k_m:
        if k_j < k_m:
#             print(f'{orb_index} belongs to {atom_sym[0]}')
            print(f'{orb_index} belongs to {atom_ind[0]}')
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[0])
#                 print(Cu_list)
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[0])
#                 print(Al_list)
                
            if atom_sym[0]=='O':
                O_list.append(atom_ind[0])
#                 print(O_list)
            if atom_sym[0]=='H':
                H_list.append(atom_ind[0])
#                 print(H_list)
            if atom_sym[0]=='N':
                N_list.append(atom_ind[0])
#                 print(N_list)
        else:
#             print(f' {orb_index} belongs to {atom_sym[1]}')
            print(f'{orb_index} belongs to {atom_ind[1]}')
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[1])
#                 print(Cu_list)
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[1])
#                 print(Al_list)
            if atom_sym[0]=='O':
                O_list.append(atom_ind[1])
#                 print(O_list)
            if atom_sym[0]=='H':
                H_list.append(atom_ind[1])
#                 print(H_list)
            if atom_sym[0]=='N':
                N_list.append(atom_ind[1])
#                 print(N_list)
                
    return Cu_list, Al_list, O_list, H_list, N_list
            
    
    
    
