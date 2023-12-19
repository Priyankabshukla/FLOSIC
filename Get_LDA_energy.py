import numpy as np
import os

directory_to_check='/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result_final_cleandata/HTBH38/'
directory_NHTBH38='/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result_final_cleandata/NHTBH38'

energy=[]

T=[]
T_sub=[]
def myfile(path):
    global energy
    with open('SUMMARY','r') as f:
        lines=f.readlines()
        lda_energy=lines[-1].split()[2]
        energy.append(lda_energy)
        print(energy)
            
def myfunction(directory):
    global T,T_sub
    for filename in os.listdir(directory):
        if "SUMMARY" in filename:
            X=directory.split('/')
            T.append(X[10])
            T_sub.append(X[11])
            print(directory+"/"+filename)
            myfile(directory+"/"+filename)
            print("T_sub", T_sub)

directories = [os.path.abspath(x[0]) for x in os.walk(directory_NHTBH38)]


for i in directories:
    if "lda_pzsic_non_SCF" in i and "opt_FOD" not in i:
        os.chdir(i)
        myfunction(i)
