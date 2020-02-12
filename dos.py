from math import *
import os

def read_poscar():
    # cell
    # [[cell[0][1], cell[0][1], cell[0][2]],
    #  [cell[1][0], cell[1][1], cell[1][2]],
    #  [cell[2][0], cell[2][1], cell[2][2]]]
    # atomic_coor
    # [[[1, 'atom1', x_coor, y_coor, z_coor], ..., [n, 'atom1', x_coor, y_coor, z_coor]]
    #  [[n+1, 'atom2', x_coor, y_coor, z_coor], ...]
    #  [...]]

    f = open('POSCAR', 'r')
    f.seek(0)

    f.readline()    # line 1 : Comment

    line = f.readline(); unit_len = float(line)   # line 2 : Unit Length
    cell = []
    for i in range(3):  # line 3-5 | Unit cell vectors
        line = f.readline(); tmp = line.split()
        cell.append([unit_len*float(tmp[0]), unit_len*float(tmp[1]), unit_len*float(tmp[2])])

    line = f.readline(); atoms = line.split()   # line 6 : Species
    line = f.readline(); tmp = line.split() # line 7 : Number of atoms
    n_atoms = []
    for i in range(len(tmp)):
        n_atoms.append(int(tmp[i]))
    #print(atoms, n_atoms)

    f.readline()    # line 8 : Direct|Cartesian

    atomic_coor = []
    n = 1
    for i, species in enumerate(atoms):
        atomic_coor_i = []
        for j in range(n_atoms[i]):
            line = f.readline(); tmp = line.split() # line 9- : Atomic positions
            atomic_coor_i.append([n, species, float(tmp[0]), float(tmp[1]), float(tmp[2])])
            n = n + 1
        atomic_coor.append(atomic_coor_i)

    f.close()

    return cell, atomic_coor

def write_tot_dos():

    f = open('DOSCAR', 'r')
    f.seek(0)

    for i in range(5):                   # line 1-5
        f.readline()

    line = f.readline(); tmp = line.split() # line 6: E(max), E(min), NEDOS,  E(fermi), 1.0000
    n_edos = int(tmp[2])
    e_fermi = float(tmp[3])

    tot_dos = []
    for i in range(n_edos):
        line = f.readline(); tmp = line.split() # line 7-: energy     DOS     integrated DOS
        tot_dos.append([float(tmp[0]), float(tmp[1]), float(tmp[2])])

    f.close()

    o = open('total.dos', 'w')

    o.write('#   Ef = %f\n' % e_fermi)
    o.write('#   E   DOS   integrated DOS\n')
    for i in range(n_edos):
        o.write('%12.3f%12.4E%12.4E\n' % (tot_dos[i][0], tot_dos[i][1], tot_dos[i][2]))

    o.close()

    return n_edos, e_fermi

def write_part_dos(atomic_coor):
    f = open('DOSCAR', 'r')
    f.seek(0)

    line = f.readline(); tmp = line.split() # line 1:
    n_atoms = int(tmp[1])                   # # of Ions (including empty spheres), # of Ions, 0 (no partial DOS) or 1 (incl. partial DOS), NCDIJ (currently not used)

    for i in range(2, 6):                   # line 2-5
        f.readline()
    
    line = f.readline(); tmp = line.split() # line 6: E(max), E(min), NEDOS,  E(fermi), 1.0000
    n_edos = int(tmp[2])
    e_fermi = float(tmp[3])

    for i in range(n_edos):
        f.readline() # lines for total dos

    for i in range(len(atomic_coor)):
        dos = []
        f. readline()

        out1 = str(atomic_coor[i][0][0]) + '_' + str(atomic_coor[i][0][1]) + '.dos'
        o1 = open(out1, 'w')

        for k in range(n_edos):
            line = f.readline(); tmp = line.split()
            dos.append([float(tmp[n]) for n in range(len(tmp))])
            o1.write('%12.3f' % float(tmp[0]))
            for j in range(len(tmp) - 1):
                o1.write('%12.4E' % float(tmp[j+1]))
            o1.write('\n')
        o1.close()

        for j in range(len(atomic_coor[i]) - 1):
            out2 = str(atomic_coor[i][j+1][0]) + '_' + str(atomic_coor[i][j+1][1]) + '.dos'
            o2 = open(out2, 'w')

            f. readline()
            for k in range(n_edos):
                line = f.readline(); tmp = line.split()
                for l in range(1, len(tmp)):
                    dos[k][l] = dos[k][l] + float(tmp[l])

                o2.write('%12.3f' % float(tmp[0]))
                for m in range(len(tmp) - 1):
                    o2.write('%12.4E' % float(tmp[m+1]))
                o2.write('\n')                  
            o2.close()
  
        out = atomic_coor[i][0][1] + '.sum.dos'
        o = open(out, 'w')

        for i in range(n_edos):
            o.write('%12.3f' % dos[i][0])
            for j in range(len(dos[0]) - 1):
                o.write('%12.4E' % dos[i][j+1])
            o.write('\n')
        
    return

def write_gnu(e_fermi):

    o = open('dos.gnu', 'w')

    o.write('set size ratio 0.5\n')
    o.write('#set key right center # font ",20" # {left|right|center} {top|bottom|center}\n')
    o.write('set zeroaxis\n')
    o.write('set tics out\n')
    o.write('set xtics 1\n')
    o.write('set mxtics 5\n')
    o.write('unset ytics\n')
    o.write('#set xra [-3:5]\n')
    o.write('set xlabel "Energy(eV)"\n')
    o.write('set ylabel "DOS(arb.)"\n')
    o.write('fermi = %f\n' % e_fermi)
    o.write('scale = %d\n' % 100)
    o.write('plot "total.dos" u ($1-fermi):($2/scale) w filledcurves lc rgb "gray" title "tot"\n')
    o.write('pause -1\n\n')

    o.write('## ISPIN = 1; LORBIT = 10; LNONCOLLINEAR = F\n')
    o.write('# energy s-DOS p-DOS d-DOS\n\n')
    
    o.write('## ISPIN = 1; LORBIT = 11; LNONCOLLINEAR = F\n')
    o.write('# energy  s  p_y p_z p_x d_{xy} d_{yz} d_{z2-r2} d_{xz} d_{x2-y2},...\n\n')
    
    o.write('## ISPIN = 2; LORBIT = 10; LNONCOLLINEAR = F\n')
    o.write('# energy s-DOS(up) s-DOS(down) p-DOS(up) p-DOS(dwn) d-DOS(up) d-DOS(dwn)\n\n')
    
    o.write('## ISPIN = 1; LORBIT = 10; LNONCOLLINEAR = T\n')
    o.write('# energy s-DOS(total) s-DOS(mx) s-DOS(my) s-DOS(mz) p-DOS(total) p-DOS(mx),...\n\n')
    
    o.close()

    return

#
# -- Main Body
#

cell, atomic_coor = read_poscar()
n_edos, e_fermi = write_tot_dos()
write_part_dos(atomic_coor)
write_gnu(e_fermi)
