# This program removes solvent molecules from the inside of an SiO2 slab
# Additionally, this program will also update a topology file, if present

import os.path
from os import path

# find newest .gro file
files = [x for x in os.listdir('./') if x.endswith('.gro')]
newest = max(files, key=os.path.getctime)
print('Using', newest, 'as configuration file.')

boundaryAtomType = 'SIC3'

print('Processing coordinate file...')
with open(newest) as f:
    # read content of .gro file, strip newline character
    lines = [line.rstrip() for line in f]

    # define headers and footer to store
    firstLine = lines[0]
    secondLine = lines[1]
    lastLine = lines[-1]

    # other variables needed
    flag = 0
    SOLcount = 0
    ACNcount = 0

    # delete headers and footer to not clash with math operations
    del lines[0]
    del lines[0]
    del lines[-1]

    # other variables
    outlist = []
    lowValue = 99.0
    highValue = 0.0
    newList = lines
    listLength = range(len(newList))
    newTop = []
    # obtain z coords for min and max boundaries from atomtype
    for i in lines:
        ResName = i[10:15].strip()
        zCoord = i[36:44].strip()
        if ResName == boundaryAtomType:
            if float(zCoord) < float(lowValue):
                lowValue = float(zCoord)
            if float(zCoord) > float(highValue):
                highValue = float(zCoord)

    # Find solvent molecules that penetrate into the boundaries and remove
    for i in listLength:
        lineIndex = lines[i]
        extender = [lines[i]]
        ResID = lineIndex[5:8].strip()
        ResName = lineIndex[10:15].strip()
        zCoord = lineIndex[36:44].strip()
        if ResName == 'OW':
            if float(zCoord) > highValue or float(zCoord) < lowValue:
                flag = 0
            else:
                flag = 1
        elif ResName == 'HW1':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
            else:
                flag = 1
        elif ResName == 'HW2':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
                OW = [lines[(i - 2)]]
                HW1 = [lines[(i - 1)]]
                outlist.extend(OW)
                outlist.extend(HW1)
                outlist.extend(extender)
            else:
                SOLcount += 1

        elif ResName == 'C1' and ResID == 'ACN':
            if float(zCoord) > highValue or float(zCoord) < lowValue:
                flag = 0
            else:
                flag = 1
        elif ResName == 'H11' and ResID == 'ACN':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
            else:
                flag = 1
        elif ResName == 'H12' and ResID == 'ACN':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
            else:
                flag = 1
        elif ResName == 'H13' and ResID == 'ACN':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
            else:
                flag = 1
        elif ResName == 'C2' and ResID == 'ACN':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
            else:
                flag = 1
        elif ResName == 'N3' and ResID == 'ACN':
            if flag == 0 and (float(zCoord) > highValue or float(zCoord) < lowValue):
                flag = 0
                C1 = [lines[(i - 5)]]
                H11 = [lines[(i - 4)]]
                H12 = [lines[(i - 3)]]
                H13 = [lines[(i - 2)]]
                C2 = [lines[(i - 1)]]
                outlist.extend(C1)
                outlist.extend(H11)
                outlist.extend(H12)
                outlist.extend(H13)
                outlist.extend(C2)
                outlist.extend(extender)
            else:
                ACNcount += 1
        else:
            outlist.extend(extender)

# update coordinate file atom count
SOLatoms = SOLcount * 3
ACNatoms = ACNcount * 6
startAtoms = int(secondLine.strip())
newAtoms = startAtoms - SOLatoms - ACNatoms
secondLine = str(newAtoms).rjust(6)

finalList = [firstLine, secondLine]
finalList.extend(outlist)
finalList.append(lastLine)
fileName = newest.split('.')
newFileName = fileName[0] + '-desolv' + '.gro'

# write new coordinate file
with open(newFileName, 'w') as outfile:
    for item in finalList:
        outfile.write("%s\n" % item)
print('Coordinate file updated.')

print('Processing topology file...')
if path.exists('topol.top'):
    with open('topol.top') as topology:
        lines = [line.rstrip() for line in topology]
        for i in lines:
            ResName = i[0:5].strip()
            nMolLength = len(i[5:])
            nMolecules = i[5:].strip()
            if ResName == 'SOL':
                nMolecules = int(nMolecules) - SOLcount
                newline = ResName.ljust(5) + str(nMolecules).rjust(15)
                newTop.append(newline)
            elif ResName == 'ACN':
                nMolecules = int(nMolecules) - ACNcount
                newline = ResName.ljust(5) + str(nMolecules).rjust(15)
                newTop.append(newline)
            else:
                newTop.append(i)
    if path.exists('topol.top.backup'):
        os.remove('./topol.top.backup')
        os.rename(r'./topol.top', r'./topol.top.backup')
    else:
        os.rename(r'./topol.top', r'./topol.top.backup')
    with open('topol.top', 'w') as outfile:
        for item in newTop:
            outfile.write("%s\n" % item)
    print('Topology file updated.')
else:
    print('No topology file found, update topology manually.')
    print('Water (SOL) molecules removed:', SOLcount)
    print('Acetonitrile (ACN) molecules removed:', ACNcount)

print('Desolvator process complete!')
