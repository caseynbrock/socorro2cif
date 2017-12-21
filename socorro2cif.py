#!/bin/env python
from __future__ import print_function
import numpy as np
import os


def main():
    """
    Converts crystal file to crystal_0.cif.
    Then, converts single new_crystal file to multiple cifs, 
    if new_crystal exists
    """
    # convert crystal file to cif
    with open('crystal') as fin:
        crystal_text_0 = fin.readlines()
    crystal_0 = Crystal(crystal_text_0)
    print('Lattice vector scale: ', crystal_0.scale)
    print('Lattice vectors (unscaled): ', crystal_0.avec, crystal_0.bvec, crystal_0.cvec)
    print('Lattice vector lengths (scaled): ', crystal_0.a, crystal_0.b, crystal_0.c)
    print('Lattice vector angles: ', crystal_0.alpha, crystal_0.beta, crystal_0.gamma)
    print('Number of atoms: ', crystal_0.num_atoms)
    crystal_0.write_cif('crystal_0.cif')
        
    # convert single new_crystal file to multiple cifs
    if os.path.isfile('new_crystal'):
        new_crystal_split = split_new_crystal()
        for i, crystal_text_i in enumerate(new_crystal_split):
            file_name = 'crystal_' + str(i+1) + '.cif'
            crystal_i = Crystal(crystal_text_i)
            crystal_i.write_cif(file_name)
    
def split_new_crystal():
    """ 
    Determines number of lines per crystal block by utilizing 
    repeated header text, the splits the new_crystal text into
    a list of crystal blocks.

    Each element of new_crystal_list is a full crystal.
    """
    with open('new_crystal') as fin:
        new_crystal_text = fin.readlines()
    num_lines = get_num_lines(new_crystal_text)
    num_blocks = len(new_crystal_text)/num_lines
    print('Number of crystal blocks in new_crystal: ', num_blocks)
    new_crystal_list = []
    for i in range(num_blocks):
        start_line = i*num_lines 
        stop_line = i*num_lines + num_lines
        new_crystal_list.append(new_crystal_text[start_line:stop_line])
    return new_crystal_list
    
def get_num_lines(new_crystal_text):
    """ determine number of lines *per block* in new_crystal file """
    header_text = new_crystal_text[0]
    for i, line in enumerate(new_crystal_text[1:]):
        if line==header_text:
            num_lines = i+1
            return num_lines

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
        self.num_atoms, self.atoms = self.get_atom_positions()
            
    def get_scale_and_lattice_vectors(self):
        """ scale is on 2nd line, vectors are on lines 3,4,5 """
        scale = float(self.crystal_text[1])
        avec = [float(x) for x in self.crystal_text[2].split()]
        bvec = [float(x) for x in self.crystal_text[3].split()]
        cvec = [float(x) for x in self.crystal_text[4].split()]
        return scale, avec, bvec, cvec
    
    def calc_abc(self):
        """ calculates the magnitude of each lattice vector """
        av = self.avec
        bv = self.bvec
        cv = self.cvec
        a = self.scale * np.sqrt(av[0]**2 + av[1]**2 + av[2]**2) 
        b = self.scale * np.sqrt(bv[0]**2 + bv[1]**2 + bv[2]**2) 
        c = self.scale * np.sqrt(cv[0]**2 + cv[1]**2 + cv[2]**2) 
        return a, b, c
    
    def calc_angles(self):
        """ calculates angles between lattice vectors """
        alpha = np.arccos(self.scale**2 * np.dot(self.bvec, self.cvec)/self.b/self.c) * 180./np.pi
        beta = np.arccos(self.scale**2 * np.dot(self.cvec, self.avec)/self.c/self.a) * 180./np.pi
        gamma = np.arccos(self.scale**2 * np.dot(self.avec, self.bvec)/self.a/self.b) * 180./np.pi
        return alpha, beta, gamma
    
    def get_atom_positions(self):
        ''' number of atoms on 7th line. Remaining lines are atoms '''
        N = int(self.crystal_text[6])
        atoms = []
        for i in range(N):
            atom_label = [self.crystal_text[i+7].split()[0]] # same as atom type for now
            atoms.append('  '.join(atom_label + self.crystal_text[i+7].split()))
        return N, atoms
    
    def write_cif(self, file_name):
        """ 
        fout is open file handle to write cif file 
        pass in a value for num to append _{num} to file name
        """
        with open(file_name, 'w') as fout:
            self.write_cif_header(fout)
            self.write_atom_positions(fout)

    def write_cif_header(self, file_handle):
        header = ('data_global \n'
                  '_audit_creation_method \'socorro2cif\' \n'
                  '_cell_length_a ' + str(self.a) + '\n'
                  '_cell_length_b ' + str(self.b) + '\n'
                  '_cell_length_c ' + str(self.c) + '\n'
                  '_cell_angle_alpha ' + str(self.alpha) + '\n'
                  '_cell_angle_beta ' + str(self.beta) + '\n'
                  '_cell_angle_gamma ' + str(self.gamma) + '\n'
                  '_symmetry_space_group_name_H-M \'P 1\' \n'
                  'loop_ \n'
                  '_atom_site_label \n'
                  '_atom_site_type_symbol \n'
                  '_atom_site_fract_x \n'
                  '_atom_site_fract_y \n'
                  '_atom_site_fract_z \n')
        file_handle.write(header)
    
    def write_atom_positions(self, file_handle):
        for atom in self.atoms:
            file_handle.write(atom+'\n')


if __name__=='__main__':
    main()
