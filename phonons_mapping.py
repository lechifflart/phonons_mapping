#!/bin/python3
#
# Mapping of the phonons in the first Brillouin Zone to phonons at Gamma of supercells
#
#
# We take frequency and eigenvector obtained at some special q point for a 1D 2-atom unit cell
#     freq (    1) =      10.266743 [THz] =     342.461674 [cm-1]
# (  0.000000   0.000000    -0.000000   0.000000    -0.312398   0.000000   )
# (  0.000000   0.000000     0.000000   0.000000     0.949951   0.000000   )
#
#
import numpy as np
import cmath

# obtained from Matdyn()
#
Natoms = 2
Nfreq = 6
alat = 1.0

with open("dynmat_nosym4.eig",'r') as dynmat_eig:
    lines = dynmat_eig.readlines()

freq_pattern = "freq "
q_pattern = "q = "

eigenvectors=np.zeros((Nfreq,Natoms),dtype=list)


for i,line in enumerate(lines):
    if q_pattern in line:
        q_point = [2*np.pi/alat*float(el) for el in line.split()[-3:]]
    if freq_pattern in line:
        freqindex = int(line.split()[2].rstrip(")"))
        for j in range(Natoms):
            eivarray = [float(el) for el in  np.asarray(lines[i+j+1].split()[1:-1])]
            eigenvectors[freqindex-1][j] = np.asarray([complex(eivarray[i],eivarray[i+3]) for i in range(3)])
        


# We want to map to a supercell N=Nx*Ny*Nz
#
Nx , Ny, Nz = 2, 2, 1
N = Nx * Ny * Nz
# first, repeat N times the atoms in the unit cell
#
new_eiv = np.zeros((Nfreq,Natoms*N),dtype=list)


# define the special q point (head of the file here) and the translation vector 
#
# example of atom coordinates
# later obtained from Pwscf()
init_atoms = np.array([[0.6666,0.3333,0.0],[-0.6666,-0.3333,0.0]])
# in crystal unit

translation_vectors = np.zeros((Natoms*N,3), dtype = float)
big_box = [Nx,Ny,Nz]
count = 0

for z in range(Nz):
    for y in range(Ny):
        for x in range(Nx):
            cell = [x,y,z]
            for j in range(Natoms):
                for coord in range(3):
                    _ = init_atoms[j%Natoms][coord]
                    if _ < 0 :
                        _ +=1
                    if _ > 1 :
                        _ -=1
                    trans_vec = (_ + cell[coord])/float(big_box[coord])
                    translation_vectors[count][coord] = trans_vec
                    
                count += 1
#print(translation_vectors)


for w in range(Nfreq):
    for j in range(Natoms*N):
        # positions in crystal units
        # how to order the additional cells ?
        #  
        new_eiv[w][j] = eigenvectors[w][j%Natoms]*np.exp(
            complex(0,1)*np.dot(q_point,translation_vectors[j]))

print(new_eiv)

