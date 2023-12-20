#Import required modules
import numpy as np
import ase
from ase.io import read, write
import matplotlib.pyplot as plt
import math
from scipy.spatial.distance import pdist, squareform
from ase.neighborlist import NeighborList
import os
import pandas as pd
import itertools

## Read SIC energies of each orbitals in a molecular and return SIC energies lower than 0.01 Ha. These orbitals will represent the delocalized orbitals in a system.
def e_sic(path):
    FLO,Charge,Asymp,e_coul,e_sic,exlocal,eclocal,exnon,ecnon = np.loadtxt(path,skiprows=1,unpack=True)
    plt.plot(e_sic,'ro')
    plt.show()
    print("Indices where SIC energies is smaller than 0.01 Ha: ",np.where(e_sic>-0.01))
    print("SIC energies smaller than 0.01 Ha:",e_sic[e_sic>-0.01])
    plt.plot(e_sic[e_sic>-0.01],'o')
    plt.show()
    return np.where(e_sic>-0.01),len(e_sic)
    
    
#Return the smallest two distances, indices and their chemical symbols. This will used to find the nereast two atoms to an electron-like FODs (Fermi-Lowdin Orbitals).
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
    
#Determine the core orbitals of H, O, Cu, and, Al atoms in a zeolite cluster. The cluster contains a total of 78 electrons. The other adsorbates such as HONO, NH3 and others can bind on the zeolite cluster to give more than 100 electrons in the system.
def find_core(direc):
    fod_atom=read(f'{direc}')
    Distance=fod_atom.get_all_distances() #Find distances of ase atom types. This gives a matrix of (N,N) where N represents the number of electrons and atoms in the system
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
    
    Cu_list=[]
    print(f'Core Cu FODs : {np.count_nonzero(Distance[Cu_indices]<0.6)-1}') # This counts the number of fods inside Cu excluding itself
    Cu_list.append(np.where(Distance[Cu_indices[0]]<0.6))
    print("Cu_inner FODs",np.where(Distance[Cu_indices[0]]<0.6))
    
    Al_list=[]
    Al_list.append(np.where(Distance[Al_indices[0]]<0.4))
    print("Al_inner FODs:",np.where(Distance[Al_indices[0]]<0.4))
    
    O_list=[]
    O_list.append(np.where(Distance[O_indices]<0.4))
    print("O_inner FODs:",list(O_list[0][1]))
    
    H_list=[]
    H_list.append(np.where(Distance[H_indices]<0.4))
    print("H_inner FODs:",list(H_list[0][1]))
    
    N_list=[]
    N_list.append(np.where(Distance[N_indices]<0.4))
    print("N_inner FODs:",list(N_list[0][1]))
    
    for k in range(len(O_indices)): # There are multiple O atoms in the zeolite cluster including the adsorbates
        print(np.where(Distance[O_indices[k]]<0.4))
    for m in range(len(H_indices)): # There are multiple H atoms in the zeolite cluster including the adsorbates
        print(np.where(Distance[H_indices[m]]<0.4))
    for n in range(len(N_indices)): # There are multiple N atoms in the zeolite cluster including the adsorbates
        print(np.where(Distance[N_indices[n]]<0.4) )
     
    #Returns a big list of the core electrons of atoms of the zeolite cluster
    big_list=[]
    big_list.append([list(Cu_list[0][0]),list(Al_list[0][0]),list(O_list[0][1]),list(H_list[0][1]),list(N_list[0][1])])
    return big_list
       
#The following function identifies the atom to which electron-like FODs (Fermi-Lowdin Orbitals) belong. This is done by creating a connectivity matrix and finding the indices where the value is 1 in the connectivity matrix. I then determine the two closest atoms (j,m) to a FOD and find if it's a bonding orbital or a lone pair orbital. The algorithm is if j_m (Distance between atom j and atom m) is less than k_j (Distance between FOD k and atom j) and j_m is less than k_m (Distance between FOD k and atom m), j_m is not the hypoteneuse of the triangle and FOD k is a lone pair of either atom j or m. If k_j is less than k_m, FOD k belongs to atom j. If j_m is greater than k_j and j_m is greater than k_m, then j_m is the hypoteneuse and FOD k belongs to a bonding orbital of either atom j or m. If k_j is less than k_m, FOD k belongs to atom j.
def fod_nearest_atom(path,orb_index,len_sic):
    atom_num=[]
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
    ind_dis=[]
    for i in ind_list:
        ind_dis.append(Distance[orb_index,i[0]])
     
    a,b,atom_ind,atom_sym=Print2Smallest(ind_dis,len(ind_dis),ind_list,fod_atom)
    j_m= Distance[atom_ind[0],atom_ind[1]]
    k_j=Distance[orb_index,atom_ind[0]]
    k_m=Distance[orb_index,atom_ind[1]]
    Cu_list=[]
    Al_list=[]
    O_list=[]
    H_list=[]
    N_list=[]

    if j_m > k_j and j_m>k_m:
        print(f'{atom_sym[0]}-{atom_sym[1]} is the hypoteneuse, then {orb_index} belongs to a bond between {atom_sym[0]},{atom_sym[1]}')
        print(f'{atom_ind[0]}-{atom_ind[1]} is the hypoteneuse, then {orb_index} belongs to a bond between {atom_ind[0]},{atom_ind[1]}')
        if k_j < k_m:
            print(f'{orb_index} belongs to {atom_ind[0]}')
            atom_num.append(atom_ind[0])
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[0])
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[0])
            if atom_sym[0]=='O':
                O_list.append(atom_ind[0])
            if atom_sym[0]=='H':
                H_list.append(atom_ind[0])
            if atom_sym[0]=='N':
                N_list.append(atom_ind[0])
            
        else:
            print(f'{orb_index} belongs to {atom_ind[1]}')
            atom_num.append(atom_ind[1])
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[1])
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[1])
            if atom_sym[0]=='O':
                O_list.append(atom_ind[1])
            if atom_sym[0]=='H':
                H_list.append(atom_ind[1])
            if atom_sym[0]=='N':
                N_list.append(atom_ind[1])
                
    elif j_m<k_j or j_m<k_m:
        print(f'{atom_sym[0]}-{atom_sym[1]} is not the hypoteneuse, then {orb_index} is a lone pair to either {atom_sym[0]} or {atom_sym[1]}')
        print(f'{atom_ind[0]}-{atom_ind[1]} is not the hypoteneuse, then {orb_index} is a lone pair to either {atom_ind[0]} or {atom_ind[1]}')
    
        if k_j < k_m:
            print(f'{orb_index} belongs to {atom_ind[0]}')
            atom_num.append(atom_ind[0])
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[0])
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[0])
            if atom_sym[0]=='O':
                O_list.append(atom_ind[0])
            if atom_sym[0]=='H':
                H_list.append(atom_ind[0])
            if atom_sym[0]=='N':
                N_list.append(atom_ind[0])
        else:
            print(f'{orb_index} belongs to {atom_ind[1]}')
            atom_num.append(atom_ind[1])
            if atom_sym[0]=='Cu':
                Cu_list.append(atom_ind[1])
            if atom_sym[0]=='Al':
                Al_list.append(atom_ind[1])
            if atom_sym[0]=='O':
                O_list.append(atom_ind[1])
            if atom_sym[0]=='H':
                H_list.append(atom_ind[1])
            if atom_sym[0]=='N':
                N_list.append(atom_ind[1])
                
    return atom_num
            
    

path=os.getcwd()
fod_path=f'{path}'+'fod_xmol.xyz'
ind_interest,len_sic=e_sic(path+'SIC_ENERGY')
core_fod=find_core(fod_path) # indices of core electrons of each atom type including the indices of each atom type

b=list(itertools.chain(*core_fod))
## List of fod numbers which are in the inner shell of Cu, Al, H, O and N. Atomic positions are also there
c=list(itertools.chain(*b))
# final list after removing atomic positions
j2 = [i for i in c if i <len_sic]

core_sic=[]
for i in j2:
    core_sic.append(e_sic[i])  #SIC_ENERGY of all core shell fods
    
np.savetxt('Core.txt',np.c_[j2,core_sic])

core_list=np.sort(np.array(j2))
full_list=list(np.arange(0,len_sic))
s = set(core_list)
valence_list = [x for x in full_list if x not in s]

valence_sic=[]
atom_belong=[]
for i in valence_list:
    atom_NN=fod_nearest_atom(fod_path,i,len_sic)
    valence_sic.append(e_sic[i])
    atom_belong.append(atom_NN[0])
    print(i,e_sic[i],atom_NN)
    
np.savetxt('Valence.txt',np.c_[atom_belong,valence_list,valence_sic])

            
    

