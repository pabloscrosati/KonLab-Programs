# PDB-ChainID_v1
# Written by Pablo M. Scrosati

# This program converts GRO files into PDB format
# If the GRO file contains multiple identical chains, the program will detect them
# and add chain IDs appropriately

# Dependencies
import os

# Default variables
res_number = []
res_name = []
atom_name = []
atom_number = []
x_coord = []
y_coord = []
z_coord = []
counter = 0
pdb_line = []
control_var = 0
iter_size = 0
char = 65
newFileName = ''
flag = 0
r = 0
file_name = ''
final_pdb = []
crys_x, crys_y, crys_z = 0, 0, 0

gro_list = []


# Finds all files with requested extension
def find_gro():
    for file in os.listdir('.'):
        if file.endswith('.gro'):
            gro_list.append(file)


# Clear functions for iterations
def clear_var():
    res_number.clear()
    res_name.clear()
    atom_name.clear()
    atom_number.clear()
    x_coord.clear()
    y_coord.clear()
    z_coord.clear()
    final_pdb.clear()
    pdb_line.clear()
    global counter, r, iter_size, control_var
    counter, r, iter_size, control_var = 0, 0, 0, 0


# Call to find GRO files
# Not really needed, but used in case additional functionality is needed later on
find_gro()

# Create output directory
if not os.path.exists('Output_PDB'):
    os.makedirs('Output_PDB')

# Main function, there is probably a better way to do this, but am leaving it as is for the current version
# Reads all GRO files and converts to PDB
# Finds identical chains
for line in gro_list:
    filename = line.split('.')
    # Retain filename
    newFileName = filename[0] + '.pdb'
    clear_var()
    with open(line) as f:
        lines = [line.rstrip() for line in f]

        # Store title and box coordinates, remove what is not needed
        first_line = lines[0]
        last_line = lines[-1]
        del lines[0]
        del lines[0]
        del lines[-1]

        # Parser for GRO file
        for i in lines:
            res_number.append(i[0:5].strip())
            res_name.append(i[5:10].strip())
            atom_name.append(i[10:15].strip())
            atom_number.append(i[15:20].strip())
            x_coord.append(i[20:28].strip())
            y_coord.append(i[28:36].strip())
            z_coord.append(i[36:44].strip())

        # Count only protein residues, may need an update or a reference list
        for i in res_name:
            if i == 'NA' or i == 'CL' or i == 'SOL':
                continue
            counter += 1

        # Detection of multiple chains, will report if only single chain detected
        for n in range(2, 101):
            if res_name[0:int(1 * counter / n)] == res_name[int(1 * counter / n):int(2 * counter / n)]:
                print(line, 'contains', n, 'protein chains.')
                break
            elif n == 100:
                print(line, 'does not appear to contain multiple chains. Will proceed assuming single chain.')
                n = 1
                break

        # PDB writer including chain ID
        iter_size = int(counter / n)
        print(iter_size)
        for i in range(len(res_name)):
            if control_var < iter_size and r < n:
                if flag == 1:
                    pdb_line.append('ATOM  ' + atom_number[i - 1].rjust(5) + ' ' + atom_name[i - 1].ljust(4)
                                    + ' ' + res_name[i - 1].rjust(3) + ' ' + chr((char + r)) + res_number[i - 1].rjust(
                        4)
                                    + '    ' + str("%.2f" % (float(x_coord[i - 1]) * 10)).rjust(8)
                                    + str("%.2f" % (float(y_coord[i - 1]) * 10)).rjust(8)
                                    + str("%.2f" % (float(z_coord[i - 1]) * 10)).rjust(8))
                    pdb_line.append('ATOM  ' + atom_number[i].rjust(5) + ' ' + atom_name[i].ljust(4)
                                    + ' ' + res_name[i].rjust(3) + ' ' + chr((char + r)) + res_number[i].rjust(4)
                                    + '    ' + str("%.2f" % (float(x_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(y_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(z_coord[i]) * 10)).rjust(8))
                    control_var += 1
                    flag = 0
                else:
                    pdb_line.append('ATOM  ' + atom_number[i].rjust(5) + ' ' + atom_name[i].ljust(4)
                                    + ' ' + res_name[i].rjust(3) + ' ' + chr((char + r)) + res_number[i].rjust(4)
                                    + '    ' + str("%.2f" % (float(x_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(y_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(z_coord[i]) * 10)).rjust(8))
                    control_var += 1
            elif r < n:
                control_var = 1
                r += 1
                flag = 1
            else:
                if flag == 1:
                    pdb_line.append('ATOM  ' + atom_number[i - 1].rjust(5) + ' ' + atom_name[i - 1].ljust(4)
                                    + ' ' + res_name[i - 1].rjust(3) + ' ' + ' ' + res_number[i - 1].rjust(4)
                                    + '    ' + str("%.2f" % (float(x_coord[i - 1]) * 10)).rjust(8)
                                    + str("%.2f" % (float(y_coord[i - 1]) * 10)).rjust(8)
                                    + str("%.2f" % (float(z_coord[i - 1]) * 10)).rjust(8))
                    pdb_line.append('ATOM  ' + atom_number[i].rjust(5) + ' ' + atom_name[i].ljust(4)
                                    + ' ' + res_name[i].rjust(3) + ' ' + ' ' + res_number[i].rjust(4)
                                    + '    ' + str("%.2f" % (float(x_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(y_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(z_coord[i]) * 10)).rjust(8))
                    control_var += 1
                    flag = 0
                else:
                    pdb_line.append('ATOM  ' + atom_number[i].rjust(5) + ' ' + atom_name[i].ljust(4)
                                    + ' ' + res_name[i].rjust(3) + ' ' + ' ' + res_number[i].rjust(4)
                                    + '    ' + str("%.2f" % (float(x_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(y_coord[i]) * 10)).rjust(8)
                                    + str("%.2f" % (float(z_coord[i]) * 10)).rjust(8))
                    control_var += 1
        pdb_line.append('END')

        # Formatting dependencies
        final_pdb.append('TITLE     ' + first_line)
        crys_x, crys_y, crys_z = ("%.2f" % (float((last_line.split())[0]) * 10)), (
                "%.2f" % (float((last_line.split())[1]) * 10)), ("%.2f" % (float((last_line.split())[2]) * 10))

        final_pdb.append('REMARK   1 PDB FILE WRITTEN WITH PDB-ChainID_v1')
        final_pdb.append('REMARK   2 PDB-ChainID_v1 WRITTEN BY PABLO M. SCROSATI')
        final_pdb.append('CRYST1' + str(crys_x.rjust(9)) + str(crys_y.rjust(9)) + str(
            crys_z.rjust(9)) + '  90.00  90.00  90.00 P 1           1')
        final_pdb.extend(pdb_line)

        # Output file portion
        with open(os.path.join('Output_PDB', newFileName), 'w') as outfile:
            for item in final_pdb:
                outfile.write("%s\n" % item)
