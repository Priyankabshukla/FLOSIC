import random
import numpy as np
import os

def perturb_xyz_file(input_file, output_file, perturbation_range):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read the number of atoms from the first line
        num_atoms = int(infile.readline())

        # Copy the comment line to the output file
        comment = infile.readline()
        outfile.write(f"{num_atoms}\n{comment}")

        # Loop through the lines with atomic coordinates
        for line in infile:
            atom_data = line.split()
            if len(atom_data) == 4:
                atom, x, y, z = atom_data
                x, y, z = float(x), float(y), float(z)

                # Perturb the coordinates
                dx = random.uniform(-pert_range, pert_range)
                print(dx,line)
                dy = random.uniform(-pert_range, pert_range)
                print(dy)
                dz = random.uniform(-pert_range, pert_range)
                print(dz)

                # Update the coordinates
                x += dx
                y += dy
                z += dz

                # Write the perturbed coordinates to the output file
                outfile.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")
            else:
                outfile.write(line)

# Input and output file names
input_file = "inputfod.xyz"  # Replace with your input XYZ file

# Range for perturbation (0 to 1 Angstrom)
pert_range = 0.5

# Perform the perturbation
for i in range(0,10):
    output_file = f"outputfod_{i}.xyz"  # Replace with your desired output XYZ file
    perturb_xyz_file(input_file, output_file, perturbation_range)

root=os.getcwd()
root

## Add atomic positions along with fod positions

for m in range(0,10):
    up_path=root+'/'+f'outputfod_{m}.xyz'
    XMOL=root+'/'+'XMOL.xyz'
    # n_fod=0
    up_infile=open(up_path,'r')
    xmol_infile=open(XMOL,'r')
    up_atm=[]
    mol_xyz=[]

    for i, line in enumerate(up_infile):
        if len(line.split())==1:
            n_fod=line.split()[0]

        if len(line.split())==4:
            up_atm.append([str(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])


    for i, line in enumerate(xmol_infile):
        if len(line.split())==1:
            nmols=line.split()[0]

        if len(line.split())==4:
            mol_xyz.append([str(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    tot_atm=int(nmols)+int(n_fod)
    #print(tot_atm)


    #         print(n_fod)

    outpath=root+'/'+f'output_fodatom_{m}.xyz'
    outfile=open(outpath,'w+')

    # outfile.write()
    outfile.write(str(tot_atm)+'\n') #pbs13
    outfile.write(  '\n')

    for j in range(len(up_atm)):
        outfile.write(up_atm[j][0]+' '+ str(up_atm[j][1])+' '+str(up_atm[j][2])+' '+str(up_atm[j][3])+'\n')


    for k in range(len(mol_xyz)):
        outfile.write(mol_xyz[k][0]+' '+ str(mol_xyz[k][1])+' '+str(mol_xyz[k][2])+' '+str(mol_xyz[k][3])+'\n')
    outfile.close()


