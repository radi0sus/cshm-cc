#!/usr/bin/env python3
# 
import numpy as np                                    # all sort of math
from itertools import permutations                    # permutations of vertices for CShM
from scipy.linalg import svd                          # SVD for CShM
from scipy.optimize import linear_sum_assignment      # Hungarian algorithm
from scipy.stats import special_ortho_group           # more evenly distributed 
                                                      # uniform random rotation matrix
import matplotlib.pyplot as plt                       # plot
                                                      
IDEAL_OC6 = np.array([
    [ 0.0           ,  0.0           , -1.080123449735],
    [ 1.080123449735,  0.0           ,  0.0           ],
    [ 0.0           ,  1.080123449735,  0.0           ],
    [-1.080123449735,  0.0,             0.0           ],
    [ 0.0           , -1.080123449735,  0.0           ],
    [ 0.0           ,  0.0           ,  1.080123449735],
    [ 0.0           ,  0.0           ,  0.0           ]
])

IDEAL_TPR6 = np.array([
    [ 0.816496580928,  0.0           , -0.707106781187],
    [-0.408248290464,  0.707106781187, -0.707106781187],
    [-0.408248290464, -0.707106781187, -0.707106781187],
    [ 0.816496580928,  0.0           ,  0.707106781187],
    [-0.408248290464,  0.707106781187,  0.707106781187],
    [-0.408248290464, -0.707106781187,  0.707106781187],
    [ 0.0           ,  0.0           ,  0.0           ]
])

IDEAL_JPPY6 = np.array([
    [ 1.146281780821,  0.0           ,  0.101205871605],
    [ 0.354220550616,  1.090178757161,  0.101205871605],
    [-0.927361441027,  0.673767525738,  0.101205871605],
    [-0.927361441027, -0.673767525738,  0.101205871605],
    [ 0.354220550616, -1.090178757161,  0.101205871605],
    [ 0.0           ,  0.0           , -0.607235229628],
    [ 0.0           ,  0.0           ,  0.101205871605]
])

IDEAL_HP6 = np.array([
    [ 1.080123449735,  0.0           , 0.0],
    [ 0.540061724867,  0.935414346693, 0.0],
    [-0.540061724867,  0.935414346693, 0.0],
    [-1.080123449735,  0.0           , 0.0],
    [-0.540061724867, -0.935414346693, 0.0],
    [ 0.540061724867, -0.935414346693, 0.0],
    [ 0.0           ,  0.0           , 0.0],
])

IDEAL_PPY6 =  np.array([
    [ 0.0           ,  0.0           , -0.937042571332],
    [ 1.093216333220,  0.0           ,  0.156173761889],
    [ 0.337822425493,  1.039710517429,  0.156173761889],
    [-0.884430592103,  0.642576438232,  0.156173761889],
    [-0.884430592103, -0.642576438232,  0.156173761889],
    [ 0.337822425493, -1.039710517429,  0.156173761889],
    [ 0.0           ,  0.0           ,  0.156173761889]
])

def normalize_structure(coordinates):
    centered_coords = coordinates - np.mean(coordinates, axis=0)
    norm = np.sqrt(np.mean(np.sum(centered_coords**2, axis=1)))
    return centered_coords / norm

def calc_cshm_exact(coordinates, ideal_shape):
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

def calc_cshm_fast(coordinates, ideal_shape, num_trials):
    input_structure = normalize_structure(coordinates)
    ideal_sq_norms = np.sum(ideal_shape**2)
    
    # try different rotations first, then optimize assignment
    min_cshm = float('inf')
    
    # generate some initial rotations to avoid local minima
    # this simulates checking multiple permutations like in the exhaustive approach
    for trial in range(num_trials):
        if trial == 0:
            # first trial with identity rotation
            R_init = np.eye(3)
        else:
            # generate a more evenly distributed uniform random rotation matrix
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

# test coordinates
coordinates = np.array([
    [ 0.00000000,  0.00000000,  0.00000000], # Fe
    [-2.04223974,  0.42694767, -0.37811642], # N
    [-0.02645702,  1.68133365, -1.47795507], # N
    [ 0.02645702, -1.68133365, -1.47795507], # N
    [ 2.04223974, -0.42694767, -0.37811642], # N
    [-1.08796929, -1.40363590,  1.24288269], # N
    [ 1.08796929,  1.40363590,  1.24288269]  # N
])

# first, calculate the values with the exact approach
OC6_ex = calc_cshm_exact(coordinates, IDEAL_OC6)
TPR6_ex = calc_cshm_exact(coordinates, IDEAL_TPR6)
JPP6_ex = calc_cshm_exact(coordinates, IDEAL_JPPY6)
PPY6_ex = calc_cshm_exact(coordinates, IDEAL_PPY6)
HP6_ex = calc_cshm_exact(coordinates, IDEAL_HP6)

# Store values for all runs and trials
oc6_values = []
tpr6_values = []
jpp6_values = []
ppy6_values = []
hp6_values = []
trial_numbers = []

for run in range(1, 101):
    for trial in range(1, 101):
        OC6 = calc_cshm_fast(coordinates, IDEAL_OC6, trial)
        TPR6 = calc_cshm_fast(coordinates, IDEAL_TPR6, trial)
        JPP6 = calc_cshm_fast(coordinates, IDEAL_JPPY6, trial)
        PPY6 = calc_cshm_fast(coordinates, IDEAL_PPY6, trial)
        HP6 = calc_cshm_fast(coordinates, IDEAL_HP6, trial)
        
        oc6_values.append(OC6)
        tpr6_values.append(TPR6)
        jpp6_values.append(JPP6)
        ppy6_values.append(PPY6)
        hp6_values.append(HP6)
        trial_numbers.append(trial)

# scatter plot with frequency-based marker sizes
plt.figure(figsize=(10, 5))
colors = ['blue', 'green', 'red', 'purple', 'orange']
labels = ['OC-6', 'TPR-6', 'JPP-6', 'PPY-6', 'HP-6']
data_sets = [oc6_values, tpr6_values, jpp6_values, ppy6_values, hp6_values]

for data, color, label in zip(data_sets, colors, labels):
    plt.scatter(trial_numbers, data, color=color, alpha=0.1, label=label)

# draw constant value lines with matching colors
# constant_values = [5.9068, 7.7281, 23.6662, 19.5521, 29.0969]
constant_values = [OC6_ex, TPR6_ex, JPP6_ex, PPY6_ex, HP6_ex]
for value, color in zip(constant_values, colors):
    plt.axhline(y=value, color=color, linestyle='--', linewidth=1)

# set the right y-axis ticks at the levels of the constant lines
plt.yticks(constant_values)

plt.xlim(0, 101)
plt.xticks(np.arange(0, 101, 10))  # set x-ticks every 10
plt.xlabel('Number of Trials')
plt.grid(True, which='both', axis='x', linestyle='dashed')
plt.ylabel('CShM ')
plt.title('CShM vs. Number of Trials (intensity indicates frequency) ')
leg = plt.legend()
for lh in leg.legend_handles: 
    lh.set_alpha(1)
plt.show()

