"""
 locpot.py
 Average LOCPOT in one direction
 Latest update: 9th Feb. 2020
"""

from sys import argv
import numpy as np
import math

script, file, direction = argv

f = open(file, 'r')

f.readline()                                   # line 1: comment
line = f.readline()                            # line 2: unit length
unit = float(line)

vec = np.zeros((3, 3))
for i in range(3):
    line = f.readline(); tmp = line.split()        # line 3-5: cell vectors
    vec[i][0] = unit*float(tmp[0])
    vec[i][1] = unit*float(tmp[1])
    vec[i][2] = unit*float(tmp[2])

f.readline()                        # line 6: atomic species

line = f.readline(); tmp = line.split() # find the number of atoms

nion = 0
for i in range(len(tmp)):
    nion = nion + int(tmp[i])

for i in range(nion + 2) : f.readline()

line = f.readline(); tmp = line.split() # find the grid

grid = [int(tmp[0]),int(tmp[1]), int(tmp[2])]
tot_grid = int(grid[0]*grid[1]*grid[2])

raw_data = np.zeros(tot_grid)

for i in range(math.ceil(tot_grid/5)):
    line = f.readline(); tmp = line.split()
    for j in range(len(tmp)):
        n = 5*i + j
        raw_data[n] = float(tmp[j])

f.close()

data = np.reshape(raw_data, (grid[2], grid[1], grid[0]))

if direction.lower() == "x":
    idir = 0
    a = 1
    b = 2
elif direction.lower() == "y":
    a = 0
    idir = 1
    b = 2
else:
    a = 0
    b = 1
    idir = 2

average = np.zeros(grid[idir])
for ipt in range(grid[idir]):
    if direction.lower() == "x":
        average[ipt] = data[:,:,ipt].sum()
    elif direction.lower() == "y":
        average[ipt] = data[:,ipt,:].sum()
    else:
        average[ipt] = data[ipt,:,:].sum()

average = average / grid[a]*grid[b]

######## write output

out = file + '_AVE'
o = open(out, 'w')
#output.write('vacuum level : %f\n' % max(ave_xy))
#output.write('#position | averaged pot | e-field\n')
o.write('#position | averaged pot\n')

for i in range(grid[idir]):
    position = i * np.linalg.norm(vec[idir]) / grid[idir]
    o.write('%11f ' % position)
    o.write('%11f\n' % average[i])

o.close()
