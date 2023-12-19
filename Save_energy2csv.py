import numpy as np
import os

import subprocess
directory_to_check = '/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result_final_cleandata/HTBH38'
T=[]
T_sub=[]

import subprocess
directory_to_check = '/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result_final_cleandata/HTBH38'
T=[]
T_sub=[]

def save_csv(global_name,mol_name,path):
    source='/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/LSIC_Refined_data'
    destination='/Users/priyankashukla/Desktop/FLOSIC/LSIC-BH76'
    global_reaction_path=root+'/'+global_name
    

################## Save LSIC sic energies from SIC-FLOSIC.txt to LSIC_Refined_data folder ########################################
def save_SIC(global_name,mol_name,path):
    root='/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/LSIC_Refined_data'
    os.chdir(root)
    global_reaction_path=root+'/'+global_name
    isexist_folder = os.path.exists(global_reaction_path)
    if not isexist_folder:
        os.mkdir(global_reaction_path)
    os.chdir(global_name)
    os.mkdir(mol_name)
    os.chdir(mol_name)
    os.system(f'cp {path}/SIC-FLOSIC.txt .')
    os.chdir('..')
    os.chdir('..')

###############Save SIC_ENERGY column 11 to SIC-FLOSIC.txt ##############################################################

def myfile(path):
            with open('SIC_ENERGY','r') as f:
                with open('SIC-FLOSIC.txt','wb') as f:
                    subprocess.call(["awk", '{print $5}', "SIC_ENERGY"], stdout=f)



############### Find SIC_ENERGY file path from lsic_final folder ##########################################
def myfunction(directory):
    global main_dir_reac, sub_dir_reac
    for filename in os.listdir(directory):
        if "SIC_ENERGY" in filename:
            X=directory.split('/')
            T.append(X[10])
            T_sub.append(X[11])
            print(directory+"/"+filename)
            myfile(directory+"/"+filename)
            
    

##################### Find path directory of lsic_final and call myfunction function #############################

directories = [os.path.abspath(x[0]) for x in os.walk(directory_to_check)]


for i in directories:
    if "lda_pzsic_non_SCF" in i:
        os.chdir(i)
        myfunction(i)
        
