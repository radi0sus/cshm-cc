#!/usr/bin/env python3
# 
import requests                                       # download COD.cif
import tempfile                                       # temp dir for COD.cif download 
import argparse                                       # argument parser
import os                                             # filename handling
import sys                                            # exit
import re                                             # regex 
import gemmi                                          # cif handling, coordinates
import numpy as np                                    # all sort of math
from itertools import permutations                    # permutations of vertices for CShM
from itertools import combinations                    # angle calculation
from itertools import cycle                           # ansi color cycle for bar graph
from scipy.linalg import svd                          # SVD for CShM
from scipy.spatial import ConvexHull,  QhullError     # polyhedral volume
from scipy.optimize import linear_sum_assignment      # Hungarian algorithm
from scipy.stats import special_ortho_group           # more evenly distributed 
                                                      # uniform random rotation matrix
from collections import Counter                       # finding duplicates
from tabulate import tabulate                         # nice tables

os.system("")                                         # for cmd.exe ANSI colors

# this will be considered as central atom
tm = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"
]

# radii collection START #################################################################    
# covalent radii from Alvarez (2008)
# DOI: 10.1039/b801115j
covalent_radii_alv = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 
    'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 
    'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 
    'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16, 
    'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 
    'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 
    'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 
    'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 
    'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 
    'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69
}

# covalent radii from Shelx
covalent_radii_shx = {
    'H': 0.32, 'He': 1.50, 'Li': 1.52, 'Be': 1.11, 'B': 0.82, 'C': 0.77, 
    'N': 0.70, 'O': 0.66, 'F': 0.64, 'Ne': 1.50, 'Na': 1.86, 'Mg': 1.60, 
    'Al': 1.25, 'Si': 1.17, 'P': 1.10, 'S': 1.03, 'Cl': 0.99, 'Ar': 1.50, 
    'K': 2.27, 'Ca': 1.97, 'Sc': 1.61, 'Ti': 1.45, 'V': 1.31, 'Cr': 1.24, 
    'Mn': 1.37, 'Fe': 1.24, 'Co': 1.25, 'Ni': 1.25, 'Cu': 1.28, 'Zn': 1.33, 
    'Ga': 1.26, 'Ge': 1.22, 'As': 1.21, 'Se': 1.17, 'Br': 1.14, 'Kr': 1.50, 
    'Rb': 2.48, 'Sr': 2.15, 'Y': 1.78, 'Zr': 1.59, 'Nb': 1.43, 'Mo': 1.36, 
    'Tc': 1.35, 'Ru': 1.33, 'Rh': 1.35, 'Pd': 1.38, 'Ag': 1.44, 'Cd': 1.49, 
    'In': 1.44, 'Sn': 1.40, 'Sb': 1.41, 'Te': 1.37, 'I': 1.33, 'Xe': 1.50, 
    'Cs': 2.65, 'Ba': 2.17, 'HF': 1.56, 'Ta': 1.43, 'W': 1.37, 'Re': 1.37, 
    'Os': 1.34, 'Ir': 1.36, 'Pt': 1.37, 'Au': 1.44, 'Hg': 1.50, 'Tl': 1.64, 
    'Pb': 1.60, 'Bi': 1.60, 'Po': 1.60, 'At': 1.60, 'Rn': 1.80, 'Fr': 2.80, 
    'Ra': 2.20, 'La': 1.87, 'Ce': 1.83, 'Pr': 1.82, 'Nd': 1.81, 'Pm': 1.81, 
    'Sm': 1.80, 'Eu': 2.00, 'Gd': 1.79, 'Tb': 1.76, 'Dy': 1.75, 'Ho': 1.74, 
    'Er': 1.73, 'Tm': 1.72, 'Yb': 1.94, 'Lu': 1.72, 'Ac': 1.90, 'Th': 1.85, 
    'Pa': 1.80, 'U': 1.80, 'Np': 1.80, 'Pu': 1.80, 'Am': 1.80, 'Cm': 1.80, 
    'Bk': 1.80, 'Cf': 1.80
}

# covalent radii from Jmol
covalent_radii_jmol = {
    'H': 0.23, 'He': 0.93, 'Li': 0.68, 'Be': 0.35, 'B': 0.83, 'C': 0.68, 
    'N': 0.68, 'O': 0.68, 'F': 0.64, 'Ne': 1.12, 'Na': 0.97, 'Mg': 1.1, 
    'Al': 1.35, 'Si': 1.2, 'P': 0.75, 'S': 1.02, 'Cl': 0.99, 'Ar': 1.57, 
    'K': 1.33, 'Ca': 0.99, 'Sc': 1.44, 'Ti': 1.47, 'V': 1.33, 'Cr': 1.35, 
    'Mn': 1.35, 'Fe': 1.34, 'Co': 1.33, 'Ni': 1.5, 'Cu': 1.52, 'Zn': 1.45, 
    'Ga': 1.22, 'Ge': 1.17, 'As': 1.21, 'Se': 1.22, 'Br': 1.21, 'Kr': 1.91, 
    'Rb': 1.47, 'Sr': 1.12, 'Y': 1.78, 'Zr': 1.56, 'Nb': 1.48, 'Mo': 1.47, 
    'Tc': 1.35, 'Ru': 1.4, 'Rh': 1.45, 'Pd': 1.5, 'Ag': 1.59, 'Cd': 1.69, 
    'In': 1.63, 'Sn': 1.46, 'Sb': 1.46, 'Te': 1.47, 'I': 1.4, 'Xe': 1.98, 
    'Cs': 1.67, 'Ba': 1.34, 'La': 1.87, 'Ce': 1.83, 'Pr': 1.82, 'Nd': 1.81, 
    'Pm': 1.8, 'Sm': 1.8, 'Eu': 1.99, 'Gd': 1.79, 'Tb': 1.76, 'Dy': 1.75, 
    'Ho': 1.74, 'Er': 1.73, 'Tm': 1.72, 'Yb': 1.94, 'Lu': 1.72, 'Hf': 1.57, 
    'Ta': 1.43, 'W': 1.37, 'Re': 1.35, 'Os': 1.37, 'Ir': 1.32, 'Pt': 1.5, 
    'Au': 1.5, 'Hg': 1.7, 'Tl': 1.55, 'Pb': 1.54, 'Bi': 1.54, 'Po': 1.68, 
    'At': 1.7, 'Rn': 2.4, 'Fr': 2, 'Ra': 1.9, 'Ac': 1.88, 'Th': 1.79, 
    'Pa': 1.61, 'U': 1.58, 'Np': 1.55, 'Pu': 1.53, 'Am': 1.51, 'Cm': 1.5, 
    'Bk': 1.5, 'Cf': 1.5, 'Es': 1.5, 'Fm': 1.5, 'Md': 1.5, 'No': 1.5, 
    'Lr': 1.5, 'Rf': 1.6, 'Db': 1.6, 'Sg': 1.6, 'Bh': 1.6, 'Hs': 1.6, 
    'Mt': 1.6
}

# max radii for each element from Alvarez, Shelx and Jmol 
# used a default
covalent_radii_max = {
    'H': 0.32, 'D': 0.32, 'He': 1.50, 'Li': 1.52, 'Be': 1.11, 'B': 0.84, 'C': 0.77, 
    'N': 0.71, 'O': 0.68, 'F': 0.64, 'Ne': 1.50, 'Na': 1.86, 'Mg': 1.60, 
    'Al': 1.35, 'Si': 1.20, 'P': 1.10, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.57, 
    'K': 2.27, 'Ca': 1.97, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 
    'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.50, 'Cu': 1.52, 'Zn': 1.45, 
    'Ga': 1.26, 'Ge': 1.22, 'As': 1.21, 'Se': 1.22, 'Br': 1.21, 'Kr': 1.91, 
    'Rb': 2.48, 'Sr': 2.15, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.45, 'Pd': 1.50, 'Ag': 1.59, 'Cd': 1.69, 
    'In': 1.63, 'Sn': 1.46, 'Sb': 1.46, 'Te': 1.47, 'I': 1.40, 'Xe': 1.98, 
    'Cs': 2.65, 'Ba': 2.17, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 
    'Pm': 1.99, 'Sm': 1.98, 'Eu': 2.00, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.94, 'Lu': 1.87, 'Hf': 1.75, 
    'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.50, 
    'Au': 1.50, 'Hg': 1.70, 'Tl': 1.64, 'Pb': 1.60, 'Bi': 1.60, 'Po': 1.68, 
    'At': 1.70, 'Rn': 2.40, 'Fr': 2.80, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.80, 
    'Bk': 1.80, 'Cf': 1.80
}

# radii collection END#T ################################################################# 

# definitions for several shapes START ###################################################

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

# for CShM (Continuous Shape Measures)
# L-2 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for Linear
IDEAL_L2 = np.array([
    [ 1.224744871392, 0.0, 0.0],
    [-1.224744871392, 0.0, 0.0],
    [ 0.0,            0.0, 0.0]
])

# for CShM (Continuous Shape Measures)
# vT-2 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for Divacant tetrahedron (V-shape, 109.47°)
IDEAL_vT2 = np.array([
    [ 0.801783725737,  0.801783725737,  0.267261241912],
    [-0.801783725737, -0.801783725737,  0.267261241912],
    [ 0.0,             0.0,            -0.534522483825]
])

# for CShM (Continuous Shape Measures)
# vOC-2 from
# https://github.com/GrupEstructuraElectronicaSimetria/
#         cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml
# definition for Tetravacant octahedron (L-shape, 90.00°)
IDEAL_vOC2 = np.array([
    [  1.0, -0.5, 0.0],
    [ -0.5,  1.0, 0.0],
    [ -0.5, -0.5, 0.0]
])
# definitions for several hhapes END #####################################################

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
            # print(f'File saved in: {file_path}')
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

def get_metals(doc, tm):
    # for automatic recognition of relevant central atoms (ca)
    # check if central atom is in CIF
    block = doc
    # check for cif items
    atoms_table = block.find(['_atom_site_label',
                              '_atom_site_type_symbol'])
    # compare with 'tm' list                         
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
            print(f"{doc.name}: Warning! No metal atoms or " 
                   "'_atom_site_type_symbol' is missing.  ")
        return None

def get_atm_prop(doc, atom, arg_dist = 0.0, inc_hydro = True, enlarge_bond = 10.0):
    # get ligand atoms that coordinate the central atom
    # doc           = cif
    # atom          = central atom
    # arg_dist      = in case the d parameter (args.distance) has been invoked
    # inc_hydro     = include hydrogen atoms (as possible ligand atoms)
    # enlarge_bonds = ad +x % to the bond lengths
    st = gemmi.make_small_structure_from_block(doc)
    ns = gemmi.NeighborSearch(st, 4).populate(include_h = inc_hydro)
    max_dist = 0.0
    min_dist = 100
    # sympops not used yet
    symops = st.symops

    for site in st.sites:
        if site.label == atom:
            # find neighbors of the central atom (ca) within min and max distance
            marks = ns.find_site_neighbors(site, max_dist = 4)
            # orthogonalize the fractional coordinates of the central atom (ca)
            cart_coord_ca = st.cell.orthogonalize(site.fract)
            break  
    # get the covalent radius of the central atom from dict of radii
    ca_radius = covalent_radii_max[site.element.name]
    
    # store several properties in a dict; here the central atom (ca)
    atoms_prop = [{'label': site.label,                         # label   Fe1
                  'element': site.element,                      # element Fe
                  'real_pos': cart_coord_ca,                    # cart. coordinates
                  # rel. position for the central atom is 0,0,0 
                  'relative_pos': np.array([0.0, 0.0, 0.0]),    
                  # distance to the central atom (is zero)
                  'distance': 0.0,                              
                  'mark_idx': 0,                         # gemmi sym op number
                  'sym_code': '.',                       # general sym code (1_555)
                  'sym_op': symops[0],                   # sym op x,y,z
                  'pbc_shift': (0, 0, 0)}]               # shift 1, 1, 1 for 1+x, 1+y, 1+z
    # now the ligand atoms
    for mark in marks:
        el_label = mark.to_site(st).element
        label = mark.to_site(st).label
        idx = mark.image_idx
        real_pos = st.cell.find_nearest_pbc_position(cart_coord_ca, mark.pos, 0)
        # convert to NumPy array
        real_pos_array = np.array([real_pos.x, real_pos.y, real_pos.z])  
        # position relative to the central atom (ca)
        relative_pos = np.array([real_pos.x - cart_coord_ca.x, 
                                 real_pos.y - cart_coord_ca.y,
                                 real_pos.z - cart_coord_ca.z])
        # get the covalent radius of the ligand atom from dict of radii
        # exit if element is not in dict
        try:
            l_radius = covalent_radii_max[el_label.name]
        except KeyError:
            print(f'{site.label} ({doc.name}): Warning! {el_label.name} is not in the '
                   'list of elements. Exit.  ') 
            sys.exit(1)
        # bond length: radius(ligand atom) + radius(central atom) + x% (10% default)
        bond = (l_radius + ca_radius) + (enlarge_bond/100.0)*(l_radius + ca_radius) 
        # distance of the central atom to the ligand atom from gemmi 
        distance = real_pos.dist(cart_coord_ca)
        
        if arg_dist:
            # if args.dist, the d parameter is invoked, use this as bond length
            bond = abs(arg_dist)
        # all distances smaller than the bond length (which ist the sum of the radii)  
        # are considered     
        if distance < bond:
            image = st.cell.find_nearest_pbc_image(cart_coord_ca, mark.pos, mark.image_idx)
            # distance > 0 because distance from ca to ca is zero
            if distance > 0.0:
                # bonds to symmetry equivalent and not symmetry equivalent central 
                # atoms are not considered, the central atom or other central atoms_prop
                # are no ligand atoms
                if label != site.label and site.element != mark.to_site(st).element:
                    # more or less same as above for central atom
                    atoms_prop.append ({'label':mark.to_site(st).label, 
                    'element': mark.to_site(st).element, 
                    'real_pos': real_pos,
                    'relative_pos': relative_pos,
                    'distance': distance,
                    'mark_idx': mark.image_idx,
                    'sym_code': image.symmetry_code(True),
                    'sym_op': symops[int(image.symmetry_code(True).split('_')[0])-1],
                    'pbc_shift': image.pbc_shift})
                    # max and min distance will be printed later
                    max_dist = max(distance, max_dist)
                    min_dist = min(distance, min_dist)
                else:
                    # usually for symmetry equivalent and not symmetry equivalent central 
                    # atoms
                    print(f'{site.label} ({doc.name}): Warning! {label} '
                          f'({image.symmetry_code(True)}) has been excluded '
                           'from coordinating atoms.  ')
    # print max and min distances and CN, CIF name and central atom name
    if max_dist and min_dist:
        print(f'{site.label} ({doc.name}): CN = {len(atoms_prop)-1}, '
              f'min dist. = {min_dist:.4f} Å, max dist. = {max_dist:.4f} Å  ')
    else:
        print(f'{site.label} ({doc.name}): Warning! No atomic distances calculated.  ') 
    # check duplicated atoms; gemmi is not taking all sym_ops or translations into account
    # maybe there will be a solution later
    # now, no idea how to circumvent
    # best is to exit
    duplicates = [atom['real_pos'].x + atom['real_pos'].y + atom['real_pos'].z for atom in atoms_prop]
    
    if any(count > 1 for count in Counter(duplicates).values()):
        print(f'{label} ({doc.name}): Warning! '
               'Symmetry related atoms on same positions. Exit.  ')
        sys.exit(1)
    
    return atoms_prop 
    
def normalize_structure(coordinates):
    # center and normalize the structure for CShM calculations
    centered_coords = coordinates - np.mean(coordinates, axis=0)
    norm = np.sqrt(np.mean(np.sum(centered_coords**2, axis=1)))
    return centered_coords / norm
    
def calc_cshm_fast(coordinates, ideal_shape, num_trials):
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
            # generate a random rotation matrix 
            
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

def get_cshm(calc_cshm, coordinates, cn, cif, metal, num_trials):
    # CShM calculation with IDEAL shapes (from top) 
    # calc_cshm selects fast or slow calculation method
    # number of trials for fast method
    # returns None for not calculated CShM for table generation
    shape_cn2 = [None] * 3
    shape_cn3 = [None] * 4
    shape_cn4 = [None] * 4
    shape_cn5 = [None] * 5
    shape_cn6 = [None] * 5
    
    if cn == 2:
        L2   = calc_cshm(coordinates, IDEAL_L2, num_trials)
        vT2  = calc_cshm(coordinates, IDEAL_vT2, num_trials)
        vOC2 = calc_cshm(coordinates, IDEAL_vOC2, num_trials)
        
        shape_cn2_values = [L2, vT2, vOC2]
        shape_cn2 = shape_cn2_values
        shape_cn2 = highlight_min_cshm(shape_cn2_values)

    elif cn == 3:
        TP3    = calc_cshm(coordinates, IDEAL_TP3, num_trials)
        vT3    = calc_cshm(coordinates, IDEAL_vT3, num_trials)
        fvOC3  = calc_cshm(coordinates, IDEAL_facvOC3, num_trials)
        mvOC3  = calc_cshm(coordinates, IDEAL_mvOC3, num_trials)
        
        shape_cn3_values = [TP3, vT3, fvOC3, mvOC3]
        shape_cn3 = shape_cn3_values
        shape_cn3 = highlight_min_cshm(shape_cn3_values)
        
    elif cn == 4:
        SP4    = calc_cshm(coordinates, IDEAL_SP4, num_trials)
        T4     = calc_cshm(coordinates, IDEAL_T4, num_trials)
        SS4    = calc_cshm(coordinates, IDEAL_SS4, num_trials)
        vTBPY4 = calc_cshm(coordinates, IDEAL_vTBPY4, num_trials)
        
        shape_cn4_values = [SP4, T4, SS4, vTBPY4]
        shape_cn4 = shape_cn4_values
        shape_cn4 = highlight_min_cshm(shape_cn4_values)
        
    elif cn == 5:
        PP5    = calc_cshm(coordinates, IDEAL_PP5, num_trials)
        vOC5   = calc_cshm(coordinates, IDEAL_vOC5, num_trials)
        TBPY5  = calc_cshm(coordinates, IDEAL_TBPY5, num_trials)
        SPY5   = calc_cshm(coordinates, IDEAL_SPY5, num_trials)
        JTBPY5 = calc_cshm(coordinates, IDEAL_JTBPY5, num_trials)
        
        shape_cn5_values = [PP5, vOC5, TBPY5, SPY5, JTBPY5]
        shape_cn5 = shape_cn5_values
        shape_cn5 = highlight_min_cshm(shape_cn5_values)
        
    elif cn == 6:
        HP6    = calc_cshm(coordinates, IDEAL_HP6, num_trials)
        PPY6   = calc_cshm(coordinates, IDEAL_PPY6, num_trials)
        OC6    = calc_cshm(coordinates, IDEAL_OC6, num_trials)
        TPR6   = calc_cshm(coordinates, IDEAL_TPR6, num_trials)
        JPPY6  = calc_cshm(coordinates, IDEAL_JPPY6, num_trials)
        
        shape_cn6_values = [HP6, PPY6, OC6, TPR6, JPPY6]
        shape_cn6 = shape_cn6_values
        shape_cn6 = highlight_min_cshm(shape_cn6_values)
        
    else:
        print(f'{metal} ({cif}): Warning! CN = {cn} is not supported for CShM.  ')
        
    return shape_cn2 + shape_cn3 + shape_cn4 + shape_cn5 + shape_cn6

def get_angles(atoms):
    # calculate angles 
    angles = []
    atom2_lbl = atoms[0]['label']
    atom2_pos = atoms[0]['real_pos']
    for atom1, atom3 in combinations(atoms[1:], 2):
        atom1_pos = atom1['real_pos']
        atom3_pos = atom3['real_pos']
        atom1_lbl = atom1['label']
        atom3_lbl = atom3['label']
        atom1_sym = atom1['sym_code']
        atom3_sym = atom3['sym_code']
        
        angle = np.degrees(gemmi.calculate_angle(atom1_pos, atom2_pos, atom3_pos))
        
        angles.append ({'atom1_label' : atom1_lbl, # atom label 1 for angle 1-2-3
                        'atom2_label' : atom2_lbl, # atom label 2 for angle 1-2-3
                        'atom3_label' : atom3_lbl, # atom label 3 for angle 1-2-3
                        # replace sym op 1_5555 with '.''
                        'atom1_sym': re.sub(r'\b1_555\b', '.', atom1_sym),
                        'atom3_sym': re.sub(r'\b1_555\b', '.', atom3_sym),
                        'angle_float' : angle})    # value of the angle
    if angles:
        return angles
    else:
        return None

def get_bonds(atoms):
    # summarize already calculated distances / bond lengths
    bonds = []
    for atom in atoms[1:]:
        bonds.append({'atom1_label' : atoms[0]['label'], 
                      'atom2_label' : atom['label'], 
                      'atom2_sym': re.sub(r'\b1_555\b', '.', atom['sym_code']),
                      'distance'    : atom['distance']})
                             
    if bonds:
        return bonds
    else:
        return None

def get_xyz(atoms):
    # summarize already calculated relative positions aka x, y, z coordinates
    xyz = []
    for atom in atoms:
        x, y, z = atom['relative_pos']
        xyz.append(f"{atom['element'].name:<2} {x:>11.8f} {y:>11.8f} {z:>11.8f}")
    if xyz:
        return xyz
    else:
        return None                         
        
def calc_octahedricity(measured_angles, cif, metal):
    # O (rms angle deviation) = sqrt ((1/15)*sum(omega_measured - omega_ideal)^2)
    # ideal angle is 90 deg if angle is < 135 deg and 180 deg if angle is > 135 deg
    cis_ang = sum(1 for angle in measured_angles if angle < 135)
    trans_ang = sum(1 for angle in measured_angles if angle > 135)
    # check if it is some sort of octahedra at all
    if cis_ang != 12 or trans_ang !=3:
        print(f'{metal} ({cif}): Warning! The number of cis or trans angles is '
               'inconsistent with an octahedron.')
        return None
     
    ideal_angles = [90 if angle <= 135 else 180 for angle in measured_angles]
    # deviation of the measured angles from ideal angles: delta omega
    deviations = [(ideal - measured) for ideal, measured in zip(ideal_angles, measured_angles)]
    # squared deviations: delta omega^2
    squared_deviations = [dev**2 for dev in deviations]
    # octahedricity
    octahedricity = np.sqrt(sum(squared_deviations) / len(squared_deviations))
    
    if octahedricity:
        return f'{octahedricity:.4f}'
    else:
        return None

def calc_geom(angle_table, cn, cif, metal):
    # calculate  τ₄, τ₄', τ₅, O (octahedricity) from angles
    # not calculated values will become None for table generation
    # get the coordination number also from angles
    angles = [angle['angle_float'] for angle in angle_table]
    cn_from_angles = len(angles)
    if len(angles) < 2:
        print(f'{metal} ({cif}): Warning! Number of angles is less than 2.  ')
        tau4 = None
        tau4impr = None
        tau5 = None
        O = None
        return [tau4,  tau4impr, tau5, O]
        
    angles.sort(reverse = True)
    beta = angles[0]
    alpha = angles[1]
    
    if cn == 3 and cn_from_angles == 3:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = None
    elif cn == 4 and cn_from_angles == 6:
        tau4 = (f'{(360.0 - (alpha + beta )) / (360.0 - 2*109.5):.4f}')
        tau4impr = (
        f'{(beta - alpha) / (360.0 - 109.5) + (180.0 - beta) / (180.0 - 109.5):.4f}'
        )
        tau5 = None
        O = None
    elif cn == 5 and cn_from_angles == 10:
        tau4 = None
        tau4impr = None
        tau5 = (f'{(beta - alpha) / 60.0:.4f}')
        O = None
    elif cn == 6 and cn_from_angles == 15:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = calc_octahedricity(angles, cif, metal)
    elif cn > 6:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = None
        # if CN > 6
        print(f'{metal} ({cif}): Warning! CN = {cn} is not supported for geometry indices.')
    else:
        tau4 = None
        tau4impr = None
        tau5 = None
        O = None
        # if coordination number from the number of ligand atoms
        # does not fit to the coordination number calc. from angles
        print(f'{metal} ({cif}): Warning! There is a mismatch of the number of '
               'coordinating atoms and the coordination number determined from '
               'the number of angles!  ')
    return [tau4,  tau4impr, tau5, O]

def calc_polyhedral_volume(coordinates, cif, metal):
    # calc polyhedral volume with ConvexHull
    # needs more than 3 atoms
    if len(coordinates) > 3:
        try:
            return f'{ConvexHull(coordinates).volume:.4f}'
        except QhullError:
            # sometimes it fails; catch the error
            print(f'{metal} ({cif}): Warning! Convex hull calculation failed.  ')
            return None  
    else:
        print(f'{metal} ({cif}): Warning! Insufficient number of atoms '
               'to calculate the polyhedral volume.  ')
        return None

def clean_table(table, prop):
    # clean the table, delete all empty rows
    if table:
        num_columns = len(next(iter(table.values())))
        all_none_indices = [i for i in range(num_columns) 
                              if all(entry[i] is None 
                              for entry in table.values())]
        # remove rows with 'None' 
        cleaned_table = {f'**{key}**': [value[i] for i in range(num_columns) 
                                    if i not in all_none_indices]
                                    for key, value in table.items()}
            
        cleaned_prop = [prop[i] for i in range(len(prop)) if i not in all_none_indices]
        cleaned_table = {'**compound**': cleaned_prop, **cleaned_table}
      
        return cleaned_table
    else:
        print('Warning! No values could be calculated. Exit.')
        sys.exit(1)

def print_prop_table(prop_dict):
    # print the table with properties
    prop = ['CN', ' τ₄', "τ₄'", 'τ₅ ', 'O', 'V /Å³', ' ',
            'L-2', 'vT-2', 'vOC-2', 
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
                     
def print_bonds(bonds):
    # print table with bonds
    # check if site symmetry is only '.'
    hide_3rd_column = all(bond['atom2_sym'] == '.' for bond in bonds)
    
    # do not print 3rd column in case site symmetry is only '.' (1_555)
    if hide_3rd_column:
        bond_lengths_table = [[f"{bond['atom1_label']}-{bond['atom2_label']}", 
                               f"{bond['distance']:.4f}"] 
                               for bond in bonds]
                               
        headers = ['**Atoms**', '**Bond length /Å**']
        col_align = 'center', 'right'
        
    else:
        bond_lengths_table = [[f"{bond['atom1_label']}-{bond['atom2_label']}", 
                               f"{bond['distance']:.4f}", 
                               bond['atom2_sym']] 
                               for bond in bonds]
                               
        headers = ['**Atoms**', '**Bond length /Å**', '**Site_Sym_2**']
        col_align = 'center', 'right', 'right'
    
    print(tabulate(bond_lengths_table, 
                  headers = headers, 
                  tablefmt='github', 
                  disable_numparse = True,
                  colalign = col_align),
                  '\n')
                    
def print_angles(angles):
    # print table with angles
    # check if site symmetry is only '.' (1_555)
    hide_3rd_column = all(angle['atom1_sym'] == '.' for angle in angles)
    # check if site symmetry is only '.'
    hide_4th_column = all(angle['atom3_sym'] == '.' for angle in angles)
    # or do not print 4th or 5th or 4th and 5th column in case site symmetry is only '.'
    if hide_3rd_column:
        angles_table = [[f"{angle['atom1_label']}-{angle['atom2_label']}-{angle['atom3_label']}",
                         f"{angle['angle_float']:.2f}", 
                         angle['atom3_sym']] 
                         for angle in angles]
                         
        headers = ['**Atoms**', '**Angle /°**', '**Site_Sym_3**']
        col_align = 'center', 'right', 'right'
        
    if hide_4th_column:
        angles_table = [[f"{angle['atom1_label']}-{angle['atom2_label']}-{angle['atom3_label']}", 
                         f"{angle['angle_float']:.2f}", 
                         angle['atom1_sym']] 
                         for angle in angles]
                         
        headers = ['**Atoms**', '**Angle /°**', '**Site_Sym_1**']
        col_align = 'center', 'right', 'right'
        
    if hide_3rd_column and hide_4th_column:
        angles_table = [[f"{angle['atom1_label']}-{angle['atom2_label']}-{angle['atom3_label']}", 
                         f"{angle['angle_float']:.2f}"] 
                         for angle in angles]
                         
        headers = ['**Atoms**', '**Angle /°**']
        col_align = 'center', 'right'
        
    if not hide_3rd_column and not hide_4th_column:
        angles_table = [[f"{angle['atom1_label']}-{angle['atom2_label']}-{angle['atom3_label']}", 
                         f"{angle['angle_float']:.2f}", 
                         angle['atom1_sym'], 
                         angle['atom3_sym']] 
                         for angle in angles]
                         
        headers = ['**Atoms**', '**Angle /°**', '**Site_Sym_1**', '**Site_Sym_3**']
        col_align = 'center', 'right', 'right', 'right'
        
    print(tabulate(angles_table, 
                  headers = headers, 
                  tablefmt='github', 
                  disable_numparse = True,
                  colalign = col_align),
                  '\n')
                  
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

def positive_int(numtrials):
    # custom function to ensure the argument is a positive integer
    # for the number of trials
    try:
        ntrials = int(numtrials)
    except ValueError:
        raise argparse.ArgumentTypeError('number of trials must be a positive integer > 0')
    if ntrials <= 0:
        raise argparse.ArgumentTypeError('number of trials must be a positive integer > 0')
    return ntrials

def plot_bars(val_dict, bar_labels, max_bar_length = 40):
    # plot bar graph for CShM in the terminal
    if not any(val is not None for sublist in val_dict.values() for val in sublist):
        print('\nWarning! No data to plot.  ')
        return 
        
    # scale for terminal representation
    max_value = max(float(val.strip('*')) 
                    for val in val_dict.values() 
                    for val in val if val is not None)
    for metal_atom, val_dict_values in val_dict.items():
        if any(val_dict_values):
            print(f'\n{metal_atom}:  ')

        # pair keys with their values and filter out None values
        paired_data = [(key, float(val.strip('*'))) 
                        for key, val in zip(bar_labels, val_dict_values) 
                        if val is not None]

        # sort the paired data by numerical value
        sorted_data = sorted(paired_data, key=lambda x: x[1])

        for key, numeric_value in sorted_data:
            bar_length = max(1, int((numeric_value / max_value) * max_bar_length))
            bar = '░' * bar_length
            print(f'{key:7}: {bar} {numeric_value:.2f}  ')

def plot_barscol(val_dict, bar_labels, max_bar_length = 40):
    # same as above but in color
    ANSI_COLORS = [
        '\033[91m',  # Red
        '\033[92m',  # Green
        '\033[93m',  # Yellow
        '\033[94m',  # Blue
        '\033[95m',  # Magenta
        '\033[96m',  # Cyan
    ]
    RESET = '\033[0m'
    
    if not any(val is not None for sublist in val_dict.values() for val in sublist):
        print('\nWarning! No data to plot.  ')
        return 
        
    # scale for terminal representation
    max_value = max(float(val.strip('*')) 
                    for val in val_dict.values() 
                    for val in val if val is not None)
    for metal_atom, val_dict_values in val_dict.items():
        if any(val_dict_values):
            print(f'\n{metal_atom}:  ')

        # pair keys with their values and filter out None values
        paired_data = [(key, float(val.strip('*'))) 
                        for key, val in zip(bar_labels, val_dict_values) 
                        if val is not None]

        # sort the paired data by numerical value
        sorted_data = sorted(paired_data, key=lambda x: x[1])
        ANSI_COLOR_CYCLE = cycle(ANSI_COLORS)
        
        for key, numeric_value in sorted_data:
            bar_length = max(1, int((numeric_value / max_value) * max_bar_length))
            color = next(ANSI_COLOR_CYCLE)
            bar = color + '█' * bar_length + RESET
            print(f'{key:7}: {bar} {numeric_value:.2f}  ')
            
# argument parser START ##################################################################

parser = argparse.ArgumentParser(prog='cshm-cc', 
        description = "Calculation of CShM and τ₄, τ₄', τ₅, and O geometry indices from "
                      "single and combined CIFs, or CIFs from the COD.")

#filename is required
parser.add_argument('filename',
    type = str,
    help = 'filename (CIF) or COD ID; e.g. mystructure.cif or 12345678')

#exclude atoms by distance
parser.add_argument('-n','--numtrials',
    type = positive_int,
    default = 124,
    help = 'number of trials > 0 for fast calculation of CShM (default is 124 trials)')
    
parser.add_argument('-ex','--exact',
    action = 'store_true',
    help = 'slower CShM calculation that always finds the global minimum')

#exclude atoms by distance
parser.add_argument('-d','--dist',
    type = float,
    help = 'exclude atoms with distances larger than d in Å; e.g. -d 2.2')

parser.add_argument('-eh','--exhydro',
     action = 'store_false',
       help = 'exclude all hydrogen atoms')

parser.add_argument('-sxyz','--savexyz',
   action = 'store_true',
     help = 'save the XYZ coordinates of the central atom and its neighboring atoms, '
            'multiple entries are combined')
            
parser.add_argument('-v','--verbose',
   action = 'store_true',
     help = 'verbose output including distances and angles')
     
parser.add_argument('-p','--plot',
   action = 'store_true',
     help = 'plot bar graphs of CShM in terminal')
     
parser.add_argument('-pc','--plotcolor',
   action = 'store_true',
     help = 'plot colored bar graphs of CShM in terminal')

#parse arguments
args = parser.parse_args()

# argument parser END ####################################################################

# main ###################################################################################

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

# generate some empty dicts
sum_table  = {}     # the main table
bond_dict  = {}     # table with bond lengths
angle_dict = {}     # table witha angles
xyz_dict   = {}     # for the XYZ file
cshm_dict  = {}     # CShM for bar graph plot

for cif in cifs:
    # iterate over cif entries and central atoms (ca or metal atoms)
    metal_atoms = get_metals(cif, tm)
    if metal_atoms:
        for metal_atom in metal_atoms:
            # get properties of atoms; coordinates, etc.
            atoms_prop = get_atm_prop(cif, metal_atom, args.dist, args.exhydro)
            # coordination number from the number of coordinating ligand atoms
            cn = len(atoms_prop)-1
            # coordinates for CShM calculation
            coordinates = np.vstack([atom['relative_pos'] for atom in atoms_prop])
            # exact calculation of CShM; slow, but always finds the global minmum
            if args.exact:
                cshm = get_cshm(calc_cshm_exact, 
                                coordinates, cn, 
                                cif.name, metal_atom,
                                args.numtrials)
             # fast but may miss the global minimum if the number of trials is insufficient  
             # around 100 trials should be sufficient (for cn 6)  
            else:
                cshm = get_cshm(calc_cshm_fast, 
                                coordinates, cn, 
                                cif.name, metal_atom, 
                                args.numtrials)
            # calculate the polyhedral volume
            ph_vol = calc_polyhedral_volume(coordinates, cif.name, metal_atom)
            # calculate angles and transfer the results
            angle_table = get_angles(atoms_prop)
            # transfer bond lengths / distances
            bond_table  = get_bonds(atoms_prop)
            # transfer XYZ coordinates
            xyz_list    = get_xyz(atoms_prop)
            
            if angle_table:
                # calculate  τ₄, τ₄', τ₅, O (octahedricity) from angles
                geom = calc_geom(angle_table, cn, cif.name, metal_atom)
                # collect results in a table
                # ca name, cif name, coordination number, geometry indices, 
                # polhedral volume, CShM 
                sum_table[f'{metal_atom} ({cif.name})'] = [cn] + \
                                                          geom + \
                                                          [ph_vol] + \
                                                          [' '] + cshm
                # table with bonds                                                
                bond_dict[f'{cif.name}-{metal_atom}'] = bond_table 
                # table with angles
                angle_dict[f'{cif.name}-{metal_atom}'] = angle_table
                # xyz coordinates
                xyz_dict[f'{cif.name}-{metal_atom}'] = xyz_list
                # CShM values for bar graphs
                cshm_dict[f'{metal_atom} ({cif.name}) CN = {cn}'] = cshm
                
                if args.verbose:
                    if file_extension:
                       # small info header 
                       # data from CIF
                        print(f'\nCIF: {cif.name}, '
                              f'M: {metal_atom}, '
                              f'CN: {cn}, Ligand atoms:', 
                            *[f"{atom['label']} ({atom['sym_code']})".replace('(1_555)','') 
                                for atom in atoms_prop[1:]],'\n')
                    else:
                       # data from COD
                       # small info header
                        print(f'\nCOD: {cif.name}, '
                              f'M: {metal_atom}, '
                              f'CN: {cn}, Ligand atoms:', 
                            *[f"{atom['label']} ({atom['sym_code']})".replace('(1_555)','') 
                                for atom in atoms_prop[1:]],'\n')
                    # print table with bond lengths
                    print_bonds(bond_table)
                    # print table with angles
                    print_angles(angle_table)

print('')
# print the table with all important values
print_prop_table(sum_table)

if args.savexyz:
# save the XYZ file
    save_xyz(filename, xyz_dict)

if args.plot:
# plot CShM as bars in bw
    plot_bars(cshm_dict,
            ['L-2', 'vT-2', 'vOC-2', 
            'TP-3', 'vT-3', 'fvOC-3', 'mvOC-3', 
            'SP-4', 'T-4', 'SS-4', 'vTBPY-4', 
            'PP-5', 'vOC-5', 'TBPY-5', 'SPY-5', 'JTBPY-5', 
            'HP-6', 'PPY-6', 'OC-6', 'TPR-6', 'JPPY-6'])
            
if args.plotcolor:
# plot CShM as bars in color
    plot_barscol(cshm_dict,
            ['L-2', 'vT-2', 'vOC-2', 
            'TP-3', 'vT-3', 'fvOC-3', 'mvOC-3', 
            'SP-4', 'T-4', 'SS-4', 'vTBPY-4', 
            'PP-5', 'vOC-5', 'TBPY-5', 'SPY-5', 'JTBPY-5', 
            'HP-6', 'PPY-6', 'OC-6', 'TPR-6', 'JPPY-6'])
