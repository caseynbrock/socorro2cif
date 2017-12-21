#!/bin/env python
import numpy as np


def main():
    # read crystal file
    with open('crystal') as fin:
        crystal = fin.readlines()
    a,b,c,alpha,beta,gamma = get_unit_cell(crystal)
    atoms = get_atom_positions(crystal)

    #write cif file
    with open('newfile', 'w') as fout:
        write_cif_header(fout, a, b, c, alpha, beta, gamma)
        write_atom_positions(fout, atoms)


def get_unit_cell(crystal_text):
    ''' a, b, and c will need to be derived from scale and primitve vectors'''
    scale, avec, bvec, cvec = get_scale_and_lattice_vectors(crystal_text)
    a,b,c, = calc_abc(scale, avec, bvec, cvec)
    alpha, beta, gamma = calc_angles(avec, bvec, cvec)
    return a,b,c,alpha,beta,gamma

        
def get_scale_and_lattice_vectors(crystal_text):
    ''' scale is on 2nd line, vectors are on lines 3,4,5 '''
    scale = float(crystal_text[1])
    avec = [float(x) for x in crystal_text[2].split()]
    bvec = [float(x) for x in crystal_text[3].split()]
    cvec = [float(x) for x in crystal_text[4].split()]
    print scale, avec, bvec, cvec
    return scale, avec, bvec, cvec


def calc_abc(scale, avec, bvec, cvec):
    ''' calculates the magnitude of each lattice vector '''
    a = scale * np.sqrt(avec[0]**2 + avec[1]**2 + avec[2]**2) 
    b = scale * np.sqrt(bvec[0]**2 + bvec[1]**2 + bvec[2]**2) 
    c = scale * np.sqrt(cvec[0]**2 + cvec[1]**2 + cvec[2]**2) 
    print a,b,c
    return a, b, c


def calc_angles(avec, bvec, cvec):
    a_unscaled = np.sqrt(avec[0]**2 + avec[1]**2 + avec[2]**2)
    b_unscaled = np.sqrt(bvec[0]**2 + bvec[1]**2 + bvec[2]**2)
    c_unscaled = np.sqrt(cvec[0]**2 + cvec[1]**2 + cvec[2]**2)
    alpha = np.arccos(np.dot(bvec, cvec)/b_unscaled/c_unscaled) * 180./np.pi
    beta = np.arccos(np.dot(cvec, avec)/c_unscaled/a_unscaled) * 180./np.pi
    gamma = np.arccos(np.dot(avec, bvec)/a_unscaled/b_unscaled) * 180/np.pi
    print alpha, beta, gamma
    return alpha, beta, gamma


def get_atom_positions(crystal_text):
    ''' number of atoms on 7th line. Remaining lines are atoms '''
    N = int(crystal_text[6])
    print N
    atoms = []
    for i in range(N):
        atom_label = [crystal_text[i+7].split()[0]] # same as atom type for now
        atoms.append('  '.join(atom_label + crystal_text[i+7].split()))
    return atoms


def write_cif_header(file_handle, a, b, c, alpha, beta, gamma):
    header = ('data_global \n'
              '_audit_creation_method \'socorro2cif\' \n'
              '_cell_length_a ' + str(a) + '\n'
              '_cell_length_b ' + str(b) + '\n'
              '_cell_length_c ' + str(c) + '\n'
              '_cell_angle_alpha ' + str(alpha) + '\n'
              '_cell_angle_beta ' + str(beta) + '\n'
              '_cell_angle_gamma ' + str(gamma) + '\n'
              '_symmetry_space_group_name_H-M \'P 1\' \n'
              'loop_ \n'
              '_atom_site_label \n'
              '_atom_site_type_symbol \n'
              '_atom_site_fract_x \n'
              '_atom_site_fract_y \n'
              '_atom_site_fract_z \n')
    file_handle.write(header)


def write_atom_positions(file_handle, atoms):
    for atom in atoms:
        file_handle.write(atom+'\n')


if __name__=='__main__':
    main()
