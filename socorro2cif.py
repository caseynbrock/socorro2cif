#!/bin/env python
import numpy as np


def main():
    # read crystal file
    with open('crystal') as fin:
        crystal_text_0 = fin.readlines()
    crystal_0 = Crystal(crystal_text_0)

    ##write cif file
    #with open('newfile', 'w') as fout:
    #    write_cif_header(fout, a, b, c, alpha, beta, gamma)
    #    write_atom_positions(fout, atoms)


class Crystal(object):
    """
    INPUTS:
    crystal_text: text of a crystal in socorro format, as read by readlines()

    ATTRIBUTES:
    scale: lattice vector scale
    avec, bvec, cvec: lattice vectors specified in socorro crystal file (unscaled)
    a, b, c: magnitude of lattice vectors (after scaling)
    a;pha, beta, gamma: angles between lattice vectors

    METHODS:
    write_cif: writes crystal in cif format
    """
    def __init__(self, crystal_text):
        self.crystal_text = crystal_text
        self.scale, self.avec, self.bvec, self.cvec = self.get_scale_and_lattice_vectors()
        self.a, self.b, self.c = self.calc_abc()
        self.alpha, self.beta, self.gamma = self.calc_angles()
        self.atoms = self.get_atom_positions()
            
    def get_scale_and_lattice_vectors(self):
        """ scale is on 2nd line, vectors are on lines 3,4,5 """
        scale = float(self.crystal_text[1])
        avec = [float(x) for x in self.crystal_text[2].split()]
        bvec = [float(x) for x in self.crystal_text[3].split()]
        cvec = [float(x) for x in self.crystal_text[4].split()]
        print scale, avec, bvec, cvec
        return scale, avec, bvec, cvec
    
    def calc_abc(self):
        """ calculates the magnitude of each lattice vector """
        av = self.avec
        bv = self.bvec
        cv = self.cvec
        a = self.scale * np.sqrt(av[0]**2 + av[1]**2 + av[2]**2) 
        b = self.scale * np.sqrt(bv[0]**2 + bv[1]**2 + bv[2]**2) 
        c = self.scale * np.sqrt(cv[0]**2 + cv[1]**2 + cv[2]**2) 
        print a,b,c
        return a, b, c
    
    def calc_angles(self):
        """ calculates angles between lattice vectors """
        alpha = np.arccos(self.scale**2 * np.dot(self.bvec, self.cvec)/self.b/self.c) * 180./np.pi
        beta = np.arccos(self.scale**2 * np.dot(self.cvec, self.avec)/self.c/self.a) * 180./np.pi
        gamma = np.arccos(self.scale**2 * np.dot(self.avec, self.bvec)/self.a/self.b) * 180./np.pi
        print alpha, beta, gamma
        return alpha, beta, gamma
    
    def get_atom_positions(self):
        ''' number of atoms on 7th line. Remaining lines are atoms '''
        N = int(self.crystal_text[6])
        print N
        atoms = []
        for i in range(N):
            atom_label = [self.crystal_text[i+7].split()[0]] # same as atom type for now
            atoms.append('  '.join(atom_label + self.crystal_text[i+7].split()))
        return atoms
    
    #def write_cif_header(file_handle, a, b, c, alpha, beta, gamma):
    #    header = ('data_global \n'
    #              '_audit_creation_method \'socorro2cif\' \n'
    #              '_cell_length_a ' + str(a) + '\n'
    #              '_cell_length_b ' + str(b) + '\n'
    #              '_cell_length_c ' + str(c) + '\n'
    #              '_cell_angle_alpha ' + str(alpha) + '\n'
    #              '_cell_angle_beta ' + str(beta) + '\n'
    #              '_cell_angle_gamma ' + str(gamma) + '\n'
    #              '_symmetry_space_group_name_H-M \'P 1\' \n'
    #              'loop_ \n'
    #              '_atom_site_label \n'
    #              '_atom_site_type_symbol \n'
    #              '_atom_site_fract_x \n'
    #              '_atom_site_fract_y \n'
    #              '_atom_site_fract_z \n')
    #    file_handle.write(header)
    #
    #
    #def write_atom_positions(file_handle, atoms):
    #    for atom in atoms:
    #        file_handle.write(atom+'\n')


if __name__=='__main__':
    main()
