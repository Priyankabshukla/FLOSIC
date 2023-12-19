import numpy as np
import os
import time

# Submit LDA DFT calculation
main_dir_path = '/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/non_hydrogen_transfer/TN9'

mol_name=['f-','ch3f','fch3fts']

for i in mol_name:
    full_path=main_dir_path+'/'+i
    os.chdir(main_dir_path)
    os.chdir(full_path)
    os.system('mkdir lda_pzsic_non_SCF')
    os.chdir('lda_pzsic_non_SCF')
    os.system('cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/T2/oh/lda_pzsic_non_SCF/GRPMAT .')
    os.system('cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/T2/oh/lda_pzsic_non_SCF/REPMAT .')
    os.system('cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/T2/oh/lda_pzsic_non_SCF/job.slurm .')
    os.system(f'cp {full_path}/ISYMGEN .')
    os.system(f'cp {full_path}/SYMBOL .')
    os.system(f'cp {full_path}/TMPTRE .')
    os.system(f'cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test/BH76-everything/all_result/non_hydrogen_transfer/RUNS .')
    os.system('sbatch job.slurm')
    os.system('cd ..')
    os.system('cd ..')
    os.system('cd ..')
    
    
#Submit LDA_SIC calculations
main_dir_path = '/bgfs/kjohnson/pbs13/FLOSIC/Projects/BH76_hypothesis_DFT_test'
mol_name=['c2h4','c2h5','c2h5ts','C2H6','c3h7','c3h7ts']
def func(mol_name):
    flag=0
    isExist = os.path.exists(main_dir_path) ##check if root directory exists
    full_path = path+'/'+mol_name
    print(full_path)
    isexist_folder = os.path.exists(full_path)
    if not isexist_folder:
        os.mkdir(full_path)
        
    os.chdir(full_path)
    os.system(f'cp {main_dir_path}/job.slurm .')
    run_ini_check_file = full_path+'/'+'SPNORB'  ### start nrlmol code first step
    isExists_SPNORB = os.path.exists(run_ini_check_file)
    if not isExists_SPNORB:
        os.system('sbatch job.slurm')
    
    time.sleep(180)
    source_dir = f'/bgfs/kjohnson/pbs13/PhD/FLOSIC/BH76/{mol_name}'

   
    DFT_run_file = full_path+'/'+'SUMMARY'
    isExists_SUMMARY = os.path.exists(DFT_run_file)
    if not isExists_SUMMARY:  ### Run LDA density with flosic
        flag=1
        os.system(f'cp {source_dir}/CLUSTER .')
        os.system('sbatch job.slurm')
    flag=0
    os.system(f'cp {main_dir_path}/NRLMOL_INPUT.DAT .')
    os.system(f'cp {source_dir}/FRMORB .')
    os.system('sbatch job.slurm')
    
        
    

    
    

    


