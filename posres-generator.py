# Generates position restraints from an existing coordinate file
# Parameters for the restraints are user-defined

# Dependencies
import tkinter as tk
from tkinter import filedialog
import os
from os import path

# User inputs for constants
restraint_force = input('Enter restraint force in kJ/mol nm^2:\n')
minmax = input('Enter upper and lower bounds in nm, separated by a comma:\n')

bounds_index = minmax.split(',')

# Determines bounds
if len(bounds_index) != 2 or bounds_index[0] == '' or float(bounds_index[0]) == float(bounds_index[1]):
    print('Error in bounds input!\nTerminating.')
    input('Press enter to exit.')
    exit()
elif float(bounds_index[0]) < float(bounds_index[1]):
    minZ = float(bounds_index[0])
    maxZ = float(bounds_index[1])
else:
    minZ = float(bounds_index[1])
    maxZ = float(bounds_index[0])

root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)

# Dialogue for selecting file
file_path = filedialog.askopenfilename(initialdir=os.getcwd(), title="Select Gromacs Coordinate File",
                                       filetypes=[("GRO Coordinate File", "*.gro")])

outfile_path = filedialog.asksaveasfilename(initialdir=os.getcwd(), title="Position Restraints File",
                                       filetypes=[("GRO Include Topology", "*.itp")])

if file_path == '':
    print('No GRO file selected!\nTerminating.')
    input('Press enter to exit.')
    exit()
else:
    print(os.path.basename(file_path), 'loaded!')

with open(file_path) as f:
    lines = [line.rstrip() for line in f]

del lines[0]
del lines[0]
del lines[-1]

atom_counter = 1

posres_index = []

reference_res_ID = lines[0][5:8].strip()

# Writes atom number
for i in lines:
    current_res_ID = i[5:8].strip()
    z_coord = (i[36:44].strip())
    if current_res_ID == reference_res_ID and (minZ < float(z_coord) < maxZ):
        posres_index.append(str(atom_counter))
        atom_counter += 1
    elif current_res_ID != reference_res_ID:
        reference_res_ID = current_res_ID
        atom_counter += 1
    else:
        atom_counter += 1

out_force = '1' + ' ' + restraint_force + ' ' + restraint_force + ' ' + restraint_force

for x in range(len(posres_index)):
    posres_index[x] = posres_index[x].rjust(5) + ' ' + out_force

posres_index.insert(0, '; atom type fx fy fz')
posres_index.insert(0, '[ position_restraints ]')

if path.exists('posres.itp'):
    os.remove('posres.itp')

with open(outfile_path, 'w') as outfile:
    for item in posres_index:
        outfile.write(("%s\n" % item))

print('Position restraint file written!')

input('Press enter to exit.')
