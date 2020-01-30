#############################################################
#
# vband.py ver.2.1.1
# Code for band & fatband plot from VASP output files
# written by Jisook
# 
# 11th Mar, 2019
#
# JisookHong@lbl.gov
#
# input files: OUTCAR, KPOINTS, EIGENVAL, PROCAR
##############################################################

from numpy import *

def read_input():
    IN = open('vband.in', 'r')
    IN.seek(0)

    line = IN.readline()
    index = line.split()

    line = IN.readline()
    tmp = line.split()
    num_orb = int(tmp[0])


    jatom = [] ; jtype = [] ; jtype_index = []

    for i in range(num_orb):        
        line = IN.readline()
        tmp = line.split()
        jatom.append(int(tmp[0]))
        jtype.append(int(tmp[1]))
        jtype_index.append(tmp[3])

    IN.close()
    return index, num_orb, jatom, jtype, jtype_index

def read_info():
    OUTCAR = open('OUTCAR', 'r')
    KPOINTS = open('KPOINTS', 'r')

    OUTCAR.seek(0)
    KPOINTS.seek(0)
 
    KPOINTS.readline()
    line = KPOINTS.readline()
    tmp = line.split()
    kgrid = int(tmp[0])

    while True:
        line = OUTCAR.readline()
        if 'NKPTS' in line:
            tmp = line.split()
            nkpts = int(tmp[3])
            nbands = int(tmp[14])
            break
    while True:
        line = OUTCAR.readline()
        if 'LSORBIT' in line:
            tmp = line.split()
            soc = tmp[2]
            break
    while True:
        line = OUTCAR.readline()
        if 'E-fermi :' in line:
            tmp = line.split()
            fermi = float(tmp[2])
            break

    OUTCAR.close()
    KPOINTS.close()
    return kgrid, nkpts, nbands, soc, fermi

def k_list(nkpts):
    OUTCAR = open('OUTCAR', 'r')

    OUTCAR.seek(0)

    klist = zeros((nkpts, 5)) 
    k_vec = zeros((3,3))  # reciprocal lattice vectors in 1/Angstrom
    # k_vec[0] = k1, k_vec[1] = k2, k_vec[2] = k3

    while True:
        line = OUTCAR.readline()
        if 'k-points in units of 2pi/SCALE and weight' in line:
            for i in range(0, nkpts):
                line = OUTCAR.readline()
                tmp = line.split()
                klist[i][0] = float(tmp[0])
                klist[i][1] = float(tmp[1])
                klist[i][2] = float(tmp[2])
            break
            


#   klist[3][0] = 0.0
    for i in range(1, nkpts):
        klist[i][3] = klist[i-1][3] + sqrt((klist[i][0] - klist[i-1][0])**2 + (klist[i][1] - klist[i-1][1])**2 + (klist[i][2] - klist[i-1][2])**2)

    for i in range(0, nkpts):
        klist[i][4] = sqrt(klist[i][0]**2 + klist[i][1]**2 + klist[i][2]**2)

    OUTCAR.close()
    return klist

def energy_list(nkpts, nbands):
    EIGENVAL = open('EIGENVAL', 'r')
    EIGENVAL.seek(0)

    elist = zeros((nbands, nkpts))

    for i in range(0, 6):
        EIGENVAL.readline()

    for i in range(0, nkpts):
        EIGENVAL.readline()
        EIGENVAL.readline()
        for j in range(0, nbands):
            line = EIGENVAL.readline()
            tmp = line.split()
            elist[j][i] = float(tmp[1])
    EIGENVAL.close()
    return elist

def weight_list(jatom, jtype):
    PROCAR = open('PROCAR', 'r')
    
    # figure out how many lines in btw 'band n ...' and 'band n+1 ...'
    PROCAR.seek(0)
    for i in range(0, 3): PROCAR.readline()
    i = 0
    while True:
        line = PROCAR.readline()
        if 'band' in line:
            while True:
                line = PROCAR.readline()
                i = i + 1
                if 'band' in line:
                    break
            break
    b2b = i
    
    # start to record band character
    PROCAR.seek(0) #; nl = 0 ###

    PROCAR.readline() #; nl = nl + 1 ###
    line = PROCAR.readline() #; nl = nl + 1 ###
    tmp = line.split()
    nkpts = int(tmp[3])
    nbands = int(tmp[7])
    nions = int(tmp[11])

    wlist = zeros((nbands, nkpts))

    for i in range(0, nkpts):
        for n in range(0, 3): PROCAR.readline() #; nl = nl + 1 ###
        for j in range(0, nbands):
            for k in range(0, 2+jatom): PROCAR.readline() #; nl = nl + 1 ###
            line = PROCAR.readline() #; nl = nl + 1 ; print(nl) ### 
            tmp = line.split()
            wlist[j][i] = float(tmp[jtype])
            for k in range(0, b2b-3-jatom) : PROCAR.readline() #; nl = nl + 1 ###
        
    PROCAR.close()
    return wlist

def write_banddat(kgrid, nkpts, nbands, klist, elist):
    DAT = open('band.banddat', 'w')
    AXIS = open('band.axis', 'w')

    nkpath = int(nkpts / kgrid)

    e1 = min(elist[0])
    e2 = max(elist[-1])

    for i in range(0, nkpath - 1): 
        for j in range(0, 30):
            e = e1 - 10. + (e2 - e1 + 20.) / 29. * j 
            AXIS.write('%-11f' % klist[kgrid * (i + 1)][3])
            AXIS.write('%-11f\n' % e)
        AXIS.write('\n')
        AXIS.write('\n')

    DAT.write('#k-path | |k| | Energy\n')
    DAT.write('#1/Angstrom | eV\n')
    for i in range(0, nbands):
        DAT.write('# band %d\n' % (i + 1))
        for j in range(0, nkpts):
            DAT.write('%-13f' % klist[j][3])
            DAT.write('%-13f' % klist[j][4])
            DAT.write('%-13f\n' % elist[i][j])
        DAT.write('\n')
        DAT.write('\n')

    AXIS.close()
    DAT.close()


def write_fatbanddat(kgrid, nkpts, nbands, klist, elist, wlist):
    DAT = open('band.banddat', 'w')
    AXIS = open('band.axis', 'w')

    nkpath = nkpts / kgrid

    e1 = min(elist[0])
    e2 = max(elist[-1])

    for i in range(0, nkpath - 1):
        for j in range(0, 30):
            e = e1 - 10. + (e2 - e1 + 20.) / 29. * j
            AXIS.write('%-11f' % klist[kgrid * (i + 1)][3])
            AXIS.write('%-11f\n' % e)
        AXIS.write('\n')
        AXIS.write('\n')

    DAT.write('#k-path   |k|   Energy(eV)   weight\n')
    for i in range(0, nbands):
        DAT.write('# band %d\n' % (i + 1))
        for j in range(0, nkpts):
            DAT.write('%-13f' % klist[j][3])
            DAT.write('%-13f' % klist[j][4])
            DAT.write('%-13f' % elist[i][j])
            for k in range(len(wlist)):
                DAT.write('%-13f' % wlist[k][i][j])
            DAT.write('\n')
        DAT.write('\n')
        DAT.write('\n')

    AXIS.close()
    DAT.close()

def write_gnuplot(index, klist, fermi, jtype_index):       
    GNU = open('band.gnuplot', 'w')

    GNU.write('#set title "Band Plot"\n')
    GNU.write('#set terminal postscript eps enhanced color\n')
    GNU.write('#set output "band.eps"\n')
    GNU.write('set size ratio 1.5\n')
    GNU.write('set key left # font ",20" # {left|right|center} {top|bottom|center}\n')
    GNU.write('set zeroaxis\n')
    GNU.write('set ytics 1\n')
    GNU.write('set mytics 5\n')
    GNU.write('unset xtics\n')
    GNU.write('set xra [%f:%f]\n' % (klist[0][3], klist[-1][3]))
    GNU.write('set yra [-5:5]\n')
    GNU.write('set ylabel "Energy(eV)"\n')
    GNU.write('set xtics (')
    GNU.write('"%s" ' % index[0]) ; GNU.write('%f,' % klist[0][3])
    for i in range(0, len(index)-2):
        GNU.write('"%s" '% index[i+1]); GNU.write('%f,' % klist[kgrid*(i+1)-1][3])
    GNU.write('"%s" ' % index[-1]) ; GNU.write('%f' % klist[-1][3]) ; GNU.write(')\n')
    GNU.write('#set xtics 0.1\n')
    GNU.write('set mxtics 5\n')
    GNU.write('fermi = %f\n' % fermi)
    GNU.write('plot "band.banddat" u ($1):($3-fermi) w line lt rgb "red" notitle,\\\n')
    for i in range(len(jtype_index)):
        GNU.write('     "band.banddat" u ($1):($3-fermi):($%d*6) w points ps variable pt 7 title "%s",\\\n' % (i+4, jtype_index[i]))
    GNU.write('     "band.axis" u ($1):($2-fermi) w line lt rgb "black" notitle\n')
    GNU.write('pause -1')

    GNU.close()

#########################################################################
#                              Main                                     #
#########################################################################

kgrid, nkpts, nbands, soc, fermi = read_info()
nkpath = nkpts / kgrid
print('\n')
print('***The number of k-points in a k-path : %d' % kgrid)
print('***Total k-paths : %d' % nkpath)
print('***Total number of bands : %d' % nbands)
print('***Consider spin-orbit coupling : %s' % soc)
print('***Fermi level : %f eV' % fermi)
print('\n')

index, num_orb, jatom, jtype, jtype_index = read_input()

klist = k_list(nkpts)
elist = energy_list(nkpts, nbands)

if num_orb == 0:
    print('Plotting bands\n')
    write_banddat(kgrid, nkpts, nbands, klist, elist)
else:
    print('Plotting fatbands\n')
    wlist = []
    for i in range(num_orb):
        wlist_i = weight_list(jatom[i], jtype[i])
        wlist.append(wlist_i)

    write_fatbanddat(kgrid, nkpts, nbands, klist, elist, wlist)

write_gnuplot(index, klist, fermi, jtype_index)
