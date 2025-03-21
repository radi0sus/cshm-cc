#!/usr/bin/env python3
# 
import requests                                       # download COD.cif
import tempfile                                       # temp dir for COD.cif download 
import argparse                                       # argument parser
import os                                             # filename handling
import sys                                            # exit
import re                                             # regex fo removing (x) from x.xx(x)
import gemmi                                          # cif handling
import numpy as np                                    # all sort of math
from itertools import permutations                    # permutations of vertices for CShM
from scipy.linalg import svd                          # SVD for CShM
from scipy.spatial import ConvexHull                  # polyhedral volume
from scipy.optimize import linear_sum_assignment      # Hungarian algorithm
from scipy.stats import special_ortho_group           # more evenly distributed 
                                                      # uniform random rotation matrix
from tabulate import tabulate                         # nice tables

# Definitions for several Shapes START ##################
# for CShM (Continuous Shape Measures)
# OC-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal Octahedron
# same as AB6
IDEAL_OC6 = np.array([
    [ 0.0           ,  0.0           , -1.080123449735],
    [ 1.080123449735,  0.0           ,  0.0           ],
    [ 0.0           ,  1.080123449735,  0.0           ],
    [-1.080123449735,  0.0,             0.0           ],
    [ 0.0           , -1.080123449735,  0.0           ],
    [ 0.0           ,  0.0           ,  1.080123449735],
    [ 0.0           ,  0.0           ,  0.0           ]
])


# for CShM (Continuous Shape Measures)
# TPR-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal trigonal prism
# same as APR_EQ
IDEAL_TPR6 = np.array([
    [ 0.816496580928,  0.0           , -0.707106781187],
    [-0.408248290464,  0.707106781187, -0.707106781187],
    [-0.408248290464, -0.707106781187, -0.707106781187],
    [ 0.816496580928,  0.0           ,  0.707106781187],
    [-0.408248290464,  0.707106781187,  0.707106781187],
    [-0.408248290464, -0.707106781187,  0.707106781187],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# JJPY-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the Johnson pentagonal pyramid (J2)
IDEAL_JPPY6 = np.array([
    [ 1.146281780821,  0.0           ,  0.101205871605],
    [ 0.354220550616,  1.090178757161,  0.101205871605],
    [-0.927361441027,  0.673767525738,  0.101205871605],
    [-0.927361441027, -0.673767525738,  0.101205871605],
    [ 0.354220550616, -1.090178757161,  0.101205871605],
    [ 0.0           ,  0.0           , -0.607235229628],
    [ 0.0           ,  0.0           ,  0.101205871605]
])

# for CShM (Continuous Shape Measures)
# HP-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the Hexagon
IDEAL_HP6 = np.array([
    [ 1.080123449735,  0.0           , 0.0],
    [ 0.540061724867,  0.935414346693, 0.0],
    [-0.540061724867,  0.935414346693, 0.0],
    [-1.080123449735,  0.0           , 0.0],
    [-0.540061724867, -0.935414346693, 0.0],
    [ 0.540061724867, -0.935414346693, 0.0],
    [ 0.0           ,  0.0           , 0.0],
])

# for CShM (Continuous Shape Measures)
# PPY-6 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the Pentagonal pyramid
IDEAL_PPY6 =  np.array([
    [ 0.0           ,  0.0           , -0.937042571332],
    [ 1.093216333220,  0.0           ,  0.156173761889],
    [ 0.337822425493,  1.039710517429,  0.156173761889],
    [-0.884430592103,  0.642576438232,  0.156173761889],
    [-0.884430592103, -0.642576438232,  0.156173761889],
    [ 0.337822425493, -1.039710517429,  0.156173761889],
    [ 0.0           ,  0.0           ,  0.156173761889]
])

# for CShM (Continuous Shape Measures)
# SPY-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal square pyramid
IDEAL_SPY5 = np.array([
    [ 0.0           ,  0.0           ,  1.095445115010],
    [ 1.060660171780,  0.0           , -0.273861278753],
    [ 0.0           ,  1.060660171780, -0.273861278753],
    [-1.060660171780,  0.0           , -0.273861278753],
    [ 0.0           , -1.060660171780, -0.273861278753],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# TBPY-5 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal trigonal bipyramid
# same as AB5_
IDEAL_TBPY5 = np.array([
    [ 0.0,             0.0,            -1.095445115010],
    [ 1.095445115010,  0.0,             0.0           ],
    [-0.547722557505,  0.948683298051,  0.0           ],
    [-0.547722557505, -0.948683298051,  0.0           ],
    [ 0.0,             0.0,             1.095445115010],
    [ 0.0,             0.0,             0.0           ],
])

# for CShM (Continuous Shape Measures)
# vOC-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal vacant octahedron (Johnson square pyramid, J1)
IDEAL_vOC5 = np.array([
    [ 0.0           ,  0.0           , -0.928476690885],
    [ 1.114172029062,  0.0           ,  0.185695338177],
    [ 0.0           ,  1.114172029062,  0.185695338177],
    [-1.114172029062,  0.0           ,  0.185695338177],
    [ 0.0           , -1.114172029062,  0.185695338177],
    [ 0.0           ,  0.0           ,  0.185695338177],
])

# for CShM (Continuous Shape Measures)
# TP-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal Pentagon
IDEAL_PP5 = np.array([
    [ 1.095445115010,  0.0           , 0.0],
    [ 0.338511156943,  1.041830214874, 0.0],
    [-0.886233714448,  0.643886483299, 0.0],
    [-0.886233714448, -0.643886483299, 0.0],
    [ 0.338511156943, -1.041830214874, 0.0],
    [ 0.0           ,  0.0           , 0.0]
])

# for CShM (Continuous Shape Measures)
# JTBPY-5 from 
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal Johnson trigonal bipyramid (J12)
IDEAL_JTBPY5 = np.array([
    [ 0.925820099773,  0.0           ,  0.0           ],
    [-0.462910049886,  0.801783725737,  0.0           ],
    [-0.462910049886, -0.801783725737,  0.0           ],
    [ 0.0           ,  0.0           ,  1.309307341416],
    [ 0.0           ,  0.0           , -1.309307341416],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# T-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal tetrahedron 
# same as AB4
IDEAL_T4 = np.array([
    [ 0.0           ,  0.912870929175, -0.645497224368],
    [ 0.0           , -0.912870929175, -0.645497224368],
    [ 0.912870929175,  0.0           ,  0.645497224368],
    [-0.912870929175,  0.0           ,  0.645497224368],
    [ 0.0           ,  0.0           ,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# SP-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal square 
# same as SQ5
IDEAL_SP4 = np.array([
    [ 1.118033988750,  0.0           , 0.0],
    [ 0.0           ,  1.118033988750, 0.0],
    [-1.118033988750,  0.0           , 0.0],
    [ 0.0           , -1.118033988750, 0.0],
    [ 0.0           ,  0.0           , 0.0],
])

# for CShM (Continuous Shape Measures)
# SS-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the ideal seesaw with center 
IDEAL_SS4 = np.array([
    [-0.235702260396, -0.235702260396, -1.178511301978],
    [ 0.942809041582, -0.235702260396,  0.0           ],
    [-0.235702260396,  0.942809041582,  0.0           ],
    [-0.235702260396, -0.235702260396,  1.178511301978],
    [-0.235702260396, -0.235702260396,  0.0           ]
])

# for CShM (Continuous Shape Measures)
# vTBPY-4 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# define the axially vacant trigonal bipyramid with center (trigonal pyramidal)
IDEAL_vTBPY4 = np.array([
    [ 0.0,             0.0,            -0.917662935482],
    [ 1.147078669353,  0.0,             0.229415733871],
    [-0.573539334676,  0.993399267799,  0.229415733871],
    [-0.573539334676, -0.993399267799,  0.229415733871],
    [ 0.0,             0.0,             0.229415733871]
])

# for CShM (Continuous Shape Measures)
# TP-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for trigonal planar
IDEAL_TP3 = np.array([
    [ 1.154700538379,  0.0, 0.0],
    [-0.577350269190,  1.0, 0.0],
    [-0.577350269190, -1.0, 0.0],
    [ 0.0           ,  0.0, 0.0]
])

# for CShM (Continuous Shape Measures)
# vT-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for Pyramid (vacant tetrahedron)
IDEAL_vT3 = np.array([
    [ 1.137070487230,  0.0,             0.100503781526],
    [-0.568535243615,  0.984731927835,  0.100503781526],
    [-0.568535243615, -0.984731927835,  0.100503781526],
    [ 0.0,             0.0,            -0.301511344578]
])

# for CShM (Continuous Shape Measures)
# fac-vOC-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for fac-Trivacant octahedron
IDEAL_facvOC3 = np.array([
    [ 1.0,            -0.333333333333, -0.333333333333],
    [-0.333333333333,  1.0,            -0.333333333333],
    [-0.333333333333, -0.333333333333,  1.0           ],
    [-0.333333333333, -0.333333333333, -0.333333333333]
])

# for CShM (Continuous Shape Measures)
# mer-vOC-3 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for mer-Trivacant octahedron (T-shape)
IDEAL_mvOC3 = np.array([
    [ 1.206045378311, -0.301511344578, 0.0],
    [ 0.0,             0.904534033733, 0.0],
    [-1.206045378311, -0.301511344578, 0.0],
    [ 0.0,            -0.301511344578, 0.0]
])
# Definitions for several Shapes END ##################

def download_cod_cif(cod_id):
    # download CIF from COD (Crystallography Open Database) 
    # https://www.crystallography.net/cod/
    # COD url
    url = 'https://www.crystallography.net/cod/'
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, f'{cod_id}.cif')
        # download CIF
            response = requests.get(f'{url}{cod_id}.cif')
            # raise error in case of problems
            response.raise_for_status()  
        # save CIF in temp dir
            with open(file_path, 'w', encoding = 'utf-8') as file:
                file.write(response.text)
            # to check for tmp folder
            #print(f'File saved in: {file_path}')
            try:
                doc = gemmi.cif.read_file(file_path)
            #file not found
            except IOError:
                print(f'{file_path} not found')
                sys.exit(1)
            #not a valid CIF
            except ValueError:
                print(f'{file_path} is not a valid CIF. Exit.')
                sys.exit(1)            
            return doc
    # CIF not found in COD
    except requests.exceptions.HTTPError:
        print(f'Warning! {cod_id}.cif not found in COD. Exit.')
        sys.exit(1)

def read_cif(cif):
    # read CIF
    try:
        doc = gemmi.cif.read_file(args.filename)
        return doc
    #file not found
    except IOError:
        print(f'{args.filename} not found')
        sys.exit(1)
    #not a valid cif
    except ValueError:
        print(f'{args.filename} is not a valid CIF. Exit.')
        sys.exit(1)

def get_metals(doc):
    # for automatic recognition of relevant central atoms (ca)
    # add or delete here
    tm = [
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
    ]
    
    # check if central atom is in CIF
    block = doc
    atoms_table = block.find(['_atom_site_label',
                              '_atom_site_type_symbol'])
                              
    atoms_table = [row['_atom_site_label'] 
                   for row in atoms_table 
                   if any(el in row['_atom_site_type_symbol'] 
                   for el in tm)]
    if atoms_table:
        return atoms_table
    else:
        # in a single combined CIF there could be a 'global' section
        # with general information
        # otherwise warn
        if doc.name != 'global':
            print(f"Warning! {doc.name} does not contain metal atoms or " 
                   "'_atom_site_type_symbol' is missing.\n")
        return None
        #sys.exit(1)

def get_bond_tables(doc, atom, exdist = None):
    # extract bonds from CIF
    block = doc
    
    # exit if no bond info is in CIF
    if not block.find(['_geom_bond_atom_site_label_1']):
        print('Warning! CIF does not contain distances. Exit.')
        sys.exit(1)
        
    # get all atom labels
    atoms_table = block.find(['_atom_site_label'])
    
    # check if a metal is in the atoms_table which contains all atom labels
    if atom not in [row[0] for row in atoms_table]:
        print(f'{atom} not found')
    
    # extract bonding information
    # there are some CIFs (from COD) where is no _geom_bond_site_symmetry_2
    # if it is missing, add it via gemmi loop 
    # the added value is always '.' (no symmetry relation)
    if block.find(['_geom_bond_site_symmetry_2']):
        bond_table = block.find(['_geom_bond_atom_site_label_1', 
                                 '_geom_bond_atom_site_label_2',
                                 '_geom_bond_distance',
                                 '_geom_bond_site_symmetry_2'])
    else:
        loop = block.find_loop('_geom_bond_distance').get_loop()
        loop.add_columns(['_geom_bond_site_symmetry_2'], value='.')
    
        bond_table = block.find(['_geom_bond_atom_site_label_1', 
                                 '_geom_bond_atom_site_label_2',
                                 '_geom_bond_distance',
                                 '_geom_bond_site_symmetry_2'])
    # check for the metal atom in the bonding information
    # and prevent symmetry related metal atoms in the bond_table
    if bond_table:
        bond_table = [row for row in bond_table if atom in row]
        bond_table = [row for row in bond_table 
                          if not (atom in row['_geom_bond_atom_site_label_2'] 
                          and '.' not in row['_geom_bond_site_symmetry_2'])]
        if bond_table:
            if exdist:
                # exclude atoms with a distance larger than exdist from bond_table
                bond_table = [row for row in bond_table 
                                  if abs(exdist) > float(re.match(r'([\d.]+)', 
                                  row['_geom_bond_distance']).group(1))]
            if bond_table:
                return bond_table
            else:
                # exit in case of too many exluded atoms
                print(f'Too many excluded atoms for {exdist} Å. Exit')
                sys.exit(1)
        else:
            # double check
            print(f'{atom} not found. Exit')
            #sys.exit(1)
    else:
        # no bonding information could be found at all
        print('CIF does not contain bonding informations. Exit')
        sys.exit(1)

def get_angle_tables(doc, atom, bond_table = None):
    block = doc
    
     # exit if no angle info is in CIF
    if not block.find(['_geom_angle_atom_site_label_1']):
        print('Warning!  CIF does not contain angles. Exit.')
        sys.exit(1)
        
    # a long way, to find a good comprehension
    # in case something is wrong, uncomment one of these
    #atms_from_btable = [row[0] and (row[1] + row[3]) for row in bond_table]
    #atms_from_btable = [row['_geom_bond_atom_site_label_2'] for row in bond_table]
    #atms_from_btable = [item for row in bond_table 
    #                         for item in (row[0] + row[3], row[1] + row[3])]
    # only atoms from bond_table should be considered for angles
    # also include symmetry info from row[3]
    # strip the '.' in case of no symmetry relation
    atms_from_btable = [item for row in bond_table 
                             for item in (row[0], row[1] + row[3].rstrip('.'))]
    
    # extract angle information
    # there are some CIFs (from COD) where is no _geom_angle_site_symmetry_1 & 3
    # if it is missing, add it via gemmi loop 
    # the added value is always '.' (no symmetry relation)    
    if block.find(['_geom_angle_site_symmetry_1']):
        angle_table = block.find(['_geom_angle_atom_site_label_1',
                                  '_geom_angle_atom_site_label_2',
                                  '_geom_angle_atom_site_label_3',
                                  '_geom_angle',
                                  '_geom_angle_site_symmetry_1',
                                  '_geom_angle_site_symmetry_3'])
    else:
        loop = block.find_loop('_geom_angle').get_loop()
        loop.add_columns(['_geom_angle_site_symmetry_1', 
                          '_geom_angle_site_symmetry_3'], value='.')
    
        angle_table = block.find(['_geom_angle_atom_site_label_1',
                                  '_geom_angle_atom_site_label_2',
                                  '_geom_angle_atom_site_label_3',
                                  '_geom_angle',
                                  '_geom_angle_site_symmetry_1',
                                  '_geom_angle_site_symmetry_3'])
    #for row in angle_table:
    #    print(row[0]+row[4]) 
    #    print(row[2]+row[5])
    #
    # the central atom (ca) is always in _geom_angle_atom_site_label_2; pos 2 in 1-2-3
    # further check for symmetry and not symmetry related atoms in pos 1 & pos 3 in 1-2-3
    # strip the '.' in case of no symmetry relation
    if angle_table:
        angle_table = [row for row in angle_table 
                           if atom in row['_geom_angle_atom_site_label_2']]
        angle_table = [row for row in angle_table 
                           if ((row[0]+row[4].rstrip('.')) and (row[2]+row[5].rstrip('.'))) 
                           in atms_from_btable] 
        if angle_table:
                return angle_table
        # exit in case of too many exluded atoms
        else:
            print(f'Too many excluded atoms. Exit')
            sys.exit(1)
    else:
        # no bonding information could be found at all
        print('CIF does not contain bonding informations. Exit')
        sys.exit(1)

def get_coordinates(doc, atom, bond_table, cn):
    # get the XYZ coordinates of the atoms for CShM calculation and the XYZ file
    # extract bond lengths and sort them to get the maximum and minmum distance
    bonds = [float(re.match(r'([\d.]+)', 
             row['_geom_bond_distance']).group(1)) 
             for row in bond_table]
    bonds.sort()
    # again list comprehension issues
    # atom labels that should be considered for coordinates and the XYZ file
    #atms_from_btable = [row[0] and row[1] for row in bond_table]
    atms_from_btable = set([item for row in bond_table for item in (row[0], row[1])])
    # prepare for gemmi neighbor search
    st = gemmi.make_small_structure_from_block(doc)
    ns = gemmi.NeighborSearch(st, 4).populate(include_h = False)
    # iterate over all sites (atoms) until the central atom is found
    for site in st.sites:
        if site.label == atom:
            # find neighbors of the central atom (ca) within min and max distance
            marks = ns.find_site_neighbors(site, min_dist = bonds[0]  - 0.05, 
                                                 max_dist = bonds[-1] + 0.05)
            # orthogonalize the fractional coordinates of the central atom (ca)
            cart_coord_ca = st.cell.orthogonalize(site.fract)
            break   
    
    # only coordination numers (cn) 3, 4, 5, and 6 are allowed
    if cn == 3 or cn == 4 or cn == 5 or cn == 6:
        # coordinates of the central atom for CShM; ca is at 0,0,0
        coordinates = np.array([0.0, 0.0, 0.0]) # the central atom is at 0, 0, 0
        # coordinates of the central atom for the XYZ file; ca is at 0,0,0
        xyz_list = [f'{site.element.name:<2} {0:>11.8f} {0:>11.8f} {0:>11.8f}']
        # list with central atom and ligand atoms
        ca_and_ligands = [atom]
        for mark in marks:
        # important: mark.pos gives position in unit cell, not outside
        # to_site and fract is useless in case of symmetry equivalents
        # pbc_position is the way to go 
        # check if the atom was excluded
            # element
            el_label = mark.to_site(st).element
            # atom label
            label = mark.to_site(st).label
            # if the atom label is in the list of considered atom labels
            if label in atms_from_btable:
                real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
                # coordinates for CShM
                neighbor_coordinate = np.array([[real_pos.x - cart_coord_ca.x, 
                                                real_pos.y - cart_coord_ca.y,
                                                real_pos.z - cart_coord_ca.z]])
                # coordinates for the XYZ file, which is: element x y z
                neighbor_xyz = (f'{el_label.name:<2}' 
                                f'{(real_pos.x - cart_coord_ca.x):>11.8f} ' 
                                f'{(real_pos.y - cart_coord_ca.y):>11.8f} ' 
                                f'{(real_pos.z - cart_coord_ca.z):>11.8f}')
                # add coordinates of neighbors
                xyz_list.append(neighbor_xyz)
                # list with central atom and coordinating ligand atoms
                ca_and_ligands.append(label)
                coordinates = np.vstack([coordinates, neighbor_coordinate])
    else:
        # if cn is different from 3, 4, 5, and 6 return None
        coordinates = np.array(None)
        xyz_list = None
        ca_and_ligands = None
    # return coordinates for CShM and the XYZ file
    return coordinates, xyz_list, ca_and_ligands
    
def calc_cn_number(bonds, angles, metal_atom, cif_name):
    # calculate the coordination number (cn) from number of bonds & angles
    if bonds == 3 and angles == 3:
        return 3
    elif bonds == 4 and angles == 6:
        return 4
    elif bonds == 5 and angles == 10:
        return 5
    elif bonds == 6 and angles == 15:
        return 6
    else:
        print(f'Warning! {metal_atom} ({cif_name}): CN not supported\n')
        return None

def calc_octahedricity(measured_angles):
    # O (rms angle deviation) = sqrt ((1/15)*sum(omega_measured - omega_ideal)^2)
    # ideal angle is 90 deg if angle is < 135 deg and 180 deg if angle is > 135 deg
    ideal_angles = [90 if angle <= 135 else 180 for angle in measured_angles]
    # deviation of the measured angles from ideal angles: delta omega
    # deviations = [(measured - ideal) for measured, ideal in zip(measured_angles, ideal_angles)]
    deviations = [(ideal - measured) for ideal, measured in zip(ideal_angles, measured_angles)]
    # squared deviations: delta omega^2
    squared_deviations = [dev**2 for dev in deviations]
    # octahedricity
    octahedricity = np.sqrt(sum(squared_deviations) / len(squared_deviations))
    return octahedricity

def calc_geom(angle_table, cn):
    # calculate  τ4, τ₄', τ5, O (octahedricity) from angles
    angles = [float(re.match(r'([\d.]+)', row['_geom_angle']).group(1)) 
              for row in angle_table]
              
    angles.sort(reverse = True)
    beta = angles[0]
    alpha = angles[1]
    
    if cn == 3:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = None
    elif cn == 4:
        tau4 = (f'{(360.0 - (alpha + beta )) / (360.0 - 2*109.5):.4f}')
        tau4impr = (f'{(beta - alpha) / (360.0 - 109.5) + (180.0 - beta) / (180.0 - 109.5):.4f}')
        tau5 = None
        O = None
    elif cn == 5:
        tau4 = None
        tau4impr = None
        tau5 = (f'{(beta - alpha) / 60.0:.4f}')
        O = None
    elif cn == 6:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = (f'{calc_octahedricity(angles):.4f}')
    else:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = None
    return [tau4,  tau4impr, tau5, O]

def normalize_structure(coordinates):
    # center and normalize the structure for CShM calculations
    centered_coords = coordinates - np.mean(coordinates, axis=0)
    norm = np.sqrt(np.mean(np.sum(centered_coords**2, axis=1)))
    return centered_coords / norm

def calc_cshm_fast(coordinates, ideal_shape, num_trials = 124):
    # faster Hungarian algorithm optimization
    # check number of trials, if it is to low, it calculates the
    # local and not the global minimum
    input_structure = normalize_structure(coordinates)
    ideal_sq_norms = np.sum(ideal_shape**2)
    # try different rotations first, then optimize assignment
    min_cshm = float('inf')
    
    # generate some initial rotations to avoid local minima
    for trial in range(num_trials):
        if trial == 0:
            # first trial with identity rotation
            R_init = np.eye(3)
        else:
            # random rotation matrix for subsequent trials
            # generate a random rotation matrix using QR decomposition
            
            #A = np.random.randn(3, 3)
            #Q, _ = np.linalg.qr(A)
            #R_init = Q
            # or 
            # more evenly distributed uniform random rotation matrix
            R_init = special_ortho_group.rvs(3)  
            
            # Ensure it's a proper rotation (det=1)
            if np.linalg.det(R_init) < 0:
                R_init[:, 0] *= -1
                
        # apply initial rotation to ideal shape
        rotated_ideal_init = np.dot(ideal_shape, R_init)
        
        # compute cost matrix based on squared Euclidean distances
        cost_matrix = np.linalg.norm(input_structure[:, None, :] - rotated_ideal_init[None, :, :], axis=2)
        
        # solve assignment problem (Hungarian algorithm)
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        # rearrange ideal_shape based on optimal assignment
        permuted_ideal = ideal_shape[col_ind]

        # compute optimal rotation using SVD
        H = np.dot(input_structure.T, permuted_ideal)
        U, _, Vt = svd(H)
        R = np.dot(Vt.T, U.T)

        rotated_ideal = np.dot(permuted_ideal, R)
        scale = np.sum(input_structure * rotated_ideal) / ideal_sq_norms
        cshm = np.mean(np.sum((input_structure - scale * rotated_ideal) ** 2, axis=1))
        
        min_cshm = min(min_cshm, cshm)

    return min_cshm * 100

def calc_cshm_exact(coordinates, ideal_shape, dummy):
    # calculation the continuous shape measures (CShM) parameter for a given structure
    # from the c++ code with some help of AI
    # https://github.com/continuous-symmetry-measure/shape
    # all permutations considered
    permut_list = list(permutations(range(len(coordinates))))
    ideal_sq_norms = np.sum(ideal_shape**2)
    input_structure = normalize_structure(coordinates)
    min_cshm = float('inf')
    for permuted_ideal in map(lambda p: ideal_shape[list(p)], permut_list):
        H = np.dot(input_structure.T, permuted_ideal)
        U, _, Vt = svd(H)
        R = np.dot(Vt.T, U.T)

        rotated_ideal = np.dot(permuted_ideal, R)
        scale = np.sum(input_structure * rotated_ideal) / ideal_sq_norms
        cshm = np.mean(np.sum((input_structure - scale * rotated_ideal)**2, axis=1))
        
        min_cshm = min(min_cshm, cshm)
        
    return min_cshm * 100

def highlight_min_cshm(values):
    # highlight lowest CShM with * * (italics in md)
    min_value = min(values)
    return [f"*{v:.4f}*" if v == min_value else f"{v:.4f}" for v in values]

def get_cshm(calc_cshm, coordinates, cn, num_trials):
    # CShM calculation with IDEAL shapes (from top) 
    # calc_cshm selects fast or slow calculation method
    # number of trials for fast method
    shape_cn3 = [None] * 4
    shape_cn4 = [None] * 4
    shape_cn5 = [None] * 5
    shape_cn6 = [None] * 6
    
    if cn == 3:
        TP3    = calc_cshm(coordinates, IDEAL_TP3, num_trials)
        vT3    = calc_cshm(coordinates, IDEAL_vT3, num_trials)
        fvOC3  = calc_cshm(coordinates, IDEAL_facvOC3, num_trials)
        mvOC3  = calc_cshm(coordinates, IDEAL_mvOC3, num_trials)
        
        shape_cn3_values = [TP3, vT3, fvOC3, mvOC3]
        shape_cn3 = highlight_min_cshm(shape_cn3_values)
        
    elif cn == 4:
        SP4    = calc_cshm(coordinates, IDEAL_SP4, num_trials)
        T4     = calc_cshm(coordinates, IDEAL_T4, num_trials)
        SS4    = calc_cshm(coordinates, IDEAL_SS4, num_trials)
        vTBPY4 = calc_cshm(coordinates, IDEAL_vTBPY4, num_trials)
        
        shape_cn4_values = [SP4, T4, SS4, vTBPY4]
        shape_cn4 = highlight_min_cshm(shape_cn4_values)
        
    elif cn == 5:
        PP5    = calc_cshm(coordinates, IDEAL_PP5, num_trials)
        vOC5   = calc_cshm(coordinates, IDEAL_vOC5, num_trials)
        TBPY5  = calc_cshm(coordinates, IDEAL_TBPY5, num_trials)
        SPY5   = calc_cshm(coordinates, IDEAL_SPY5, num_trials)
        JTBPY5 = calc_cshm(coordinates, IDEAL_JTBPY5, num_trials)
        
        shape_cn5_values = [PP5, vOC5, TBPY5, SPY5, JTBPY5]
        shape_cn5 = highlight_min_cshm(shape_cn5_values)
        
    elif cn == 6:
        HP6    = calc_cshm(coordinates, IDEAL_HP6, num_trials)
        PPY6   = calc_cshm(coordinates, IDEAL_PPY6, num_trials)
        OC6    = calc_cshm(coordinates, IDEAL_OC6, num_trials)
        TPR6   = calc_cshm(coordinates, IDEAL_TPR6, num_trials)
        JPPY6  = calc_cshm(coordinates, IDEAL_JPPY6, num_trials)
        
        shape_cn6_values = [HP6, PPY6, OC6, TPR6, JPPY6]
        shape_cn6 = highlight_min_cshm(shape_cn6_values)
        
    else:
        print('CN not supported')
        #sys.exit(1)
    return shape_cn3 + shape_cn4 + shape_cn5 + shape_cn6

def clean_table(table, prop):
    # clean the table, delete all empty rows
    if table:
        num_columns = len(next(iter(table.values())))
    
        # check rows for 'None'
        all_none_indices = [i for i in range(num_columns) if all(entry[i] is None for entry in table.values())]
    
        # remove rows with 'None' 
        cleaned_table = {
            f'**{key}**': [value[i] for i in range(num_columns) if i not in all_none_indices]
            for key, value in table.items()
        }
            
        cleaned_prop = [prop[i] for i in range(len(prop)) if i not in all_none_indices]
        cleaned_table = {'**compound**': cleaned_prop, **cleaned_table}
      
        return cleaned_table
    else:
        print('Warning! No values could be calculated. Exit.')
        sys.exit(1)
    
def print_bonds(bond_table):
    # print table with bonds
    # check if site symmetry is only '.'
    hide_3rd_column = all(row[3] == '.' for row in bond_table)
    
    # do not print 3rd column in case site symmetry is only '.'
    if hide_3rd_column:
        bond_lengths_table = [[f'{row[0]}-{row[1]}', row[2]] for row in bond_table]
        headers = ['**Atoms**', '**Bond length /Å**']
        col_align = 'center', 'right'
    else:
        bond_lengths_table = [[f'{row[0]}-{row[1]}', row[2], row[3]] for row in bond_table]
        headers = ['**Atoms**', '**Bond length /Å**', '**Site_Sym_2**']
        col_align = 'center', 'right', 'right'
    
    print(tabulate(bond_lengths_table, 
                  headers = headers, 
                  tablefmt='github', 
                  colalign = col_align),
                  '\n')

def print_angles(angle_table):
    # print table with angles
    # check if site symmetry is only '.'
    hide_4th_column = all(row[4] == '.' for row in angle_table)
    # check if site symmetry is only '.'
    hide_5th_column = all(row[5] == '.' for row in angle_table)
    
    # or do not print 4th or 5th or 4th and 5th column in case site symmetry is only '.'
    if hide_4th_column:
        angles_table = [[f'{row[0]}-{row[1]}-{row[2]}', row[3], row[5]] for row in angle_table]
        headers = ['**Atoms**', '**Angle /°**', '**Site_Sym_3**']
        col_align = 'center', 'right', 'right'
    if hide_5th_column:
        angles_table = [[f'{row[0]}-{row[1]}-{row[2]}', row[3], row[4]] for row in angle_table]
        headers = ['**Atoms**', '**Angle /°**', '**Site_Sym_1**']
        col_align = 'center', 'right', 'right'
    if hide_4th_column and hide_5th_column:
        angles_table = [[f'{row[0]}-{row[1]}-{row[2]}', row[3]] for row in angle_table]
        headers = ['**Atoms**', '**Angle /°**']
        col_align = 'center', 'right'
    else:
        angles_table = [[f'{row[0]}-{row[1]}-{row[2]}', row[3], row[4], row[5]] for row in angle_table]
        headers = ['**Atoms**', '**Angle /°**', '**Site_Sym_1**', '**Site_Sym_3**']
        col_align = 'center', 'right', 'right', 'right'
        
    print(tabulate(angles_table, 
                  headers = headers, 
                  tablefmt='github', 
                  colalign = col_align),
                  '\n')

def print_prop_table(prop_dict):
    # print the table with properties
    prop = ['CN', ' τ₄', "τ₄'", 'τ₅ ', 'O', 'V /Å³', ' ',
            'TP-3', 'vT-3', 'fvOC-3', 'mvOC-3', 
            'SP-4', 'T-4', 'SS-4', 'vTBPY-4', 
            'PP-5', 'vOC-5', 'TBPY-5', 'SPY-5', 'JTBPY-5', 
            'HP-6', 'PPY-6', 'OC-6', 'TPR-6', 'JPPY-6']
    # remove rows containing only 'None' 
    cleaned_table = clean_table(prop_dict, prop)
    
    print(tabulate(cleaned_table, 
                   headers = cleaned_table.keys(), 
                   tablefmt='github', 
                   stralign="right"))

def save_xyz(filename, xyz_dict):
    # save XYZ coordinates in a single file
    xyz_filename = f'{filename}.xyz'
    try:
        with open(xyz_filename, 'w', encoding='utf-8') as xyz_file:
            for key, values in xyz_dict.items():
                xyz_file.write(f'{len(values)}\n')
                xyz_file.write(f'{key}\n')
                for value in values:
                    xyz_file.write(f'{value}\n')
        print(f'\nXYZ file saved to {xyz_filename}')
    # file not found -> exit here
    except IOError:
        print("Write error. Exit.")
        sys.exit(1)    
  
# argument parser START ##################################################################
parser = argparse.ArgumentParser(prog='cshm-calc', 
        description = "Calculation of CShM and τ₄, τ₄', τ₅, and O geometry indices from "
                      "single and combined CIFs, or CIFs from the COD.")

#filename is required
parser.add_argument('filename',
    type = str,
    help = 'filename, CIF or COD number; e.g. mystructure.cif or 12345678')

#exclude atoms by distance
parser.add_argument('-n','--numtrials',
    type = int,
    default = 124,
    help = 'number of trials > 0 for fast calculation of CShM (default is 124 trials)')
    
parser.add_argument('-ex','--exact',
    action = 'store_true',
    help = 'slower CShM calculation that always finds the global minmum')

#exclude atoms by distance
parser.add_argument('-d','--dist',
    type = float,
    help = 'exclude atoms with distances larger than d in Å; e.g. -d 2.2')

parser.add_argument('-sxyz','--savexyz',
   action = 'store_true',
     help = 'save the XYZ coordinates of the central atom and its neighboring atoms, '
            'multiple entries are combined')
            
parser.add_argument('-v','--verbose',
   action = 'store_true',
     help = 'verbose output including distances and angles')

#parse arguments
args = parser.parse_args()

# argument parser END ####################################################################

# open file
# if the file has the extension CIF, read the file from hard disk
# if the file has no extension, assume it is a COD number
filename, file_extension = os.path.splitext(args.filename)
if file_extension == ('.cif' or '.CIF' or '.Cif'):
    cifs = read_cif(args.filename)
elif file_extension == '':
    cifs = download_cod_cif(args.filename)
else:
    print('Warning! File type not supported. Exit.')
    sys.exit(1)
    
# generate several dicts for tables
sum_table = {}
xyz_dict = {}
ca_and_l_dict = {}
bond_dict = {}
angle_dict = {}

# read a single CIF or a CIF with multiple entries or the CIF from the COD
for cif in cifs:
    # identify the central atom, the metal atom
    metal_atoms = get_metals(cif)
    # if there are metal atoms, continue and iterate over each metal atom
    if metal_atoms:
        for metal_atom in metal_atoms:
            # generate table with bond lengths
            bond_table = get_bond_tables(cif, metal_atom, exdist = args.dist)
            # generate table with angles
            angle_table = get_angle_tables(cif, metal_atom, bond_table)
            # calculate the coordination number (cn)
            cn = calc_cn_number(len(bond_table), len(angle_table), metal_atom, cif.name)
            # it the coordination number is 3, 4, 5, or 6, continue
            if cn:
                # calculate  τ4, τ₄', τ5, O (octahedricity) from angles
                geom = calc_geom(angle_table, cn)
                # generate coordinates, create content for the XYZ file,  
                # and prepare a list with the central atom along with ligand atoms  
                coordinates, xyz_list, ca_and_ligands = get_coordinates(cif, 
                                                                        metal_atom, 
                                                                        bond_table, cn)
                
                # exact calculation of CShM; slow, but always finds the global minmum
                if args.exact:
                    cshm = get_cshm(calc_cshm_exact, coordinates, cn, abs(args.numtrials))
                # fast but may miss the global minimum if the number of trials is insufficient  
                # around 100 trials should be sufficient (for cn 6)  
                else:
                    cshm = get_cshm(calc_cshm_fast, coordinates, cn, abs(args.numtrials))
                # table with CShM and geometry values, also calculate the polyhedral volume (ConvexHull) 
                sum_table[f'{metal_atom} ({cif.name})'] = [cn] + \
                                                          geom + \
                                                          [f'{ConvexHull(coordinates).volume:.4f}'] + \
                                                          [' '] + cshm
                # XYZ file
                xyz_dict[f'{cif.name}-{metal_atom}'] = xyz_list
                # table with bond lengths
                bond_dict[f'{cif.name}-{metal_atom}'] = bond_table 
                # table with angles
                angle_dict[f'{cif.name}-{metal_atom}'] = angle_table
                
                # in verbose mode, print bond lengths, angles, and other information  
                if args.verbose:
                    if file_extension:
                        # small info header 
                        # data from CIF
                        print(f'\nCIF: {cif.name}, '
                              f'M: {metal_atom}, '
                              f'CN: {cn}, Ligand atoms:', 
                              *ca_and_ligands[1:],'\n')
                    else:
                        # data from COD
                        print(f'\nCOD: {cif.name}, '
                              f'M: {metal_atom}, CN: {cn}, Ligand atoms:', 
                              *ca_and_ligands[1:],'\n')
                              
                    print_bonds(bond_table)
                    print_angles(angle_table)

# print the table with properties
print_prop_table(sum_table)

# save the XYZ file
if args.savexyz:
    save_xyz(filename, xyz_dict)
        
