#!/usr/bin/env python3
import numpy as np
import os


os.system(f'cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/Zeolites_SIC/Cu-SSZ-13/Final/fod_opt/nrlmolDft.opa.batch .')
os.system(f'cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/Zeolites_SIC/Cu-SSZ-13/Final/fod_opt/RUNS .')
os.system(f'cp ../SYMBOL .')
os.system('cp ../ISYMGEN .')
os.system('cp ../GRPMAT .')
os.system('cp ../REPMAT .')
os.system('cp /bgfs/kjohnson/pbs13/FLOSIC/Projects/Zeolites_SIC/Cu-SSZ-13/Final/fod_opt/IE/Ge+/r2scan/TMPTRE .')


with open('SYMBOL', 'r') as file :
          filedata = file.read()

# Replace the target string
filedata = filedata.replace('LDA-PW91*LDA-PW91', 'MGGA-R2SCAN*MGGA-R2SCAN')

# Write the file out again
with open('SYMBOL', 'w') as file:
      file.write(filedata)

with open('nrlmolDft.opa.batch', 'r') as file :
          filedata = file.read()


   
