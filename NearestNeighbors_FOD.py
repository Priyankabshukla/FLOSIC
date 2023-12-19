import numpy as np
import ase
import matplotlib.pyplot as plt
from ase.io import read, write
import math
from scipy.spatial.distance import pdist, squareform
from ase.neighborlist import NeighborList
import os
import pandas as pd

def Print2Smallest(a,n,ind_list,fod_atom):
    VALUE = max(a) + 1
    idx_min0 = a.index(min(a))
    print("idx_min0: ",idx_min0)
    min0 = a[idx_min0]
    print("min0: ", min0)
    a[idx_min0] = VALUE
    print(a[idx_min0])

    idx_min1 = a.index(min(a))
    print("idx_min1:",idx_min1)
    min1 = a[idx_min1]
    a[idx_min1] = VALUE
    a[idx_min0] = min0
    a[idx_min1] = min1
    
    
    return ([min0,min1],[idx_min0,idx_min1],[ind_list[idx_min0][0],ind_list[idx_min1][0]],
         [fod_atom.get_chemical_symbols()[ind_list[idx_min0][0]],
         fod_atom.get_chemical_symbols()[ind_list[idx_min1][0]]]
         )
    

def fod_nearest_atom(path,orb_index,len_sic):
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

    # Print the connectivity matrix
    print("Connectivity matrix:")

    print(connectivity[orb_index])
    
    indices=np.argwhere((connectivity[orb_index]==1) & (np.arange(len(fod_atom)) > len_sic))
    print(indices)
    
    ind_list=list(indices)
    ind_list
    ind_dis=[]
    for i in ind_list:
        ind_dis.append(Distance[96,i[0]])
    #     print(Distance[96,i[0]])
     
    a,b,atom_ind,atom_sym=Print2Smallest(ind_dis,len(ind_dis),ind_list,fod_atom)
    atom_ind
    j_m= Distance[atom_ind[0],atom_ind[1]] #H-N distance
    print(j_m)
    k_j=Distance[orb_index,atom_ind[0]] #47-H distance
    print(k_j)
    k_m=Distance[orb_index,atom_ind[1]] #47-N distance
    print(k_m)

    if j_m > k_j and j_m>k_m:
        print(f'{atom_sym[0]}-{atom_sym[1]} is the hypoteneuse, then {orb_index} belongs to a bond between {atom_sym[0]},{atom_sym[1]}')
        if k_j < k_m:
            print(f'{orb_index} belongs to {atom_sym[0]}')
        else:
            print(f'{orb_index} belongs to {atom_sym[1]}')
    elif j_m<k_j or j_m<k_m:
        print(f'{atom_sym[0]}-{atom_sym[1]} is not the hypoteneuse, then {orb_index} is a lone pair to either {atom_sym[0]} or {atom_sym[1]}')
        if k_j < k_m:
            print(f'{orb_index} belongs to {atom_sym[0]}')
        else:
            print(f' {orb_index} belongs to {atom_sym[1]}')

    
    

    
    
    
    
