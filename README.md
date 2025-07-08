# cshm-cc

A Python 3 script that calculates CShM (Continuous Shape Measures) values for various shapes, geometry indices, including τ₄, τ₄', τ₅, *O* (octahedricity) of automatically selected atoms (transition metal atoms by default) from a crystallographic information file (CIF). The CIF may contain one or more entries. It can also calculate these values from COD (Crystallography Open Database) entries by simply entering the COD ID. It also saves the XYZ coordinates of the central atom and its neighboring atoms, including those from CIFs with multiple entries or COD entries and calculates the polyhedral volume. Optionally, tables in Markdown format containing bond lengths and angles can be generated.

## Introduction
This script calculates the Continuous Shape Measures (CShM) and geometry indices (τ₄, τ₄', τ₅, and *O* (octahedricity)) to assign coordination geometries to two-, three-, four-, five and six-coordinate transition metal atoms. The indices τ₄ and τ₅ are used to determine whether a compound adopts tetrahedral, trigonal pyramidal, square planar, or seesaw geometry (for four-coordinated compounds), and square pyramidal or trigonal bipyramidal geometry (for five-coordinated compounds). These assignments rely on the two largest angles enclosing the central atom. The octahedricity index, *O*, is calculated based on experimental X-M-X angles and assesses how close the geometry is to an ideal octahedron. The CShM value approaches zero when the shape closely matches the ideal geometry. Additionally, the polyhedral volume is calculated using the Convex Hull Algorithm (via SciPy and the Qhull Library). For more information, refer to the associated papers and references.

#### Equations and values for the geometry indices τ₄, τ₄', τ₅, and *O*:

$\boldsymbol{\tau_4} = \frac{360°-(\alpha+\beta)}{141°}$ 

Square planar geometry: $\tau_4 = 0$; Tetrahedral geometry:  $\tau_4 = 1$; Seesaw geometry: $\tau_4 \approx 0.43$

$\boldsymbol{\tau'_4} = \frac{\beta-\alpha}{250.5°}+\frac{180-\beta}{70.5°}$

Square planar geometry: $\tau'_4 = 0$; Tetrahedral geometry:  $\tau'_4 = 1$; Seesaw geometry: $\tau'_4 \approx 0.24$

$\boldsymbol{\tau_5} = \frac{\beta-\alpha}{60°}$

Square pyramidal geometry: $\tau_5 = 0$; Trigonal bipyramidal geometry: $\tau_5 = 1$

$\alpha$ and $\beta$ are the two largest angles, with $\beta > \alpha$.  

$\boldsymbol{O} = \sqrt{\frac{1}{15}\sum_{i=1}^{15}(\hat{\theta_i} - \theta)^2}$

$\hat{\theta_i}$ = 180° for *trans* X-M-X angles and 90° for *cis* X-M-X angles\
$\theta$ = experimental X-M-X angles

$O$ is close to zero for an almost ideal octahedron

## External modules
 `gemmi >= 0.7.3`, `numpy`, `scipy`, `tabulate`, `requests` 
 
## Quick start
For local CIFs with one or more entries start the script with:
```console
python3 cshm-cc.py example.cif
```
The following output will be printed: 

    𝗰𝗲𝗻𝘁𝗿𝗮𝗹 𝗮𝘁𝗼𝗺   𝗰𝗼𝗼𝗿𝗱𝗶𝗻𝗮𝘁𝗶𝗼𝗻 𝗻𝘂𝗺𝗯𝗲𝗿 
     ↓                ↓
    Fe1 (example): CN = 5, min dist. = 1.9799 Å, max dist. = 2.3664 Å 
           ↑                  ↑                     ↑
    𝗖𝗜𝗙 𝗲𝗻𝘁𝗿𝘆 𝗼𝗿 𝗖𝗢𝗗 𝗜𝗗       𝗺𝗶𝗻𝗶𝗺𝘂𝗺 𝗱𝗶𝘀𝘁𝗮𝗻𝗰𝗲        𝗺𝗮𝘅𝗶𝗺𝘂𝗺 𝗱𝗶𝘀𝘁𝗮𝗻𝗰𝗲 (𝘁𝗼 𝗰𝗲𝗻𝘁𝗿𝗮𝗹 𝗮𝘁𝗼𝗺)       
    
    |   **compound** |   **Fe1 (example)**     |  ← 𝗖𝗲𝗻𝘁𝗿𝗮𝗹 𝗮𝘁𝗼𝗺 (𝗖𝗜𝗙 𝗲𝗻𝘁𝗿𝘆) 
    |----------------|-------------------------|
    |             CN |                       5 |  ← 𝗖𝗼𝗼𝗿𝗱𝗶𝗻𝗮𝘁𝗶𝗼𝗻 𝗻𝘂𝗺𝗯𝗲𝗿 (𝗖𝗡)
    |             τ₅ |                  0.3053 |  ← 𝗚𝗲𝗼𝗺𝗲𝘁𝗿𝘆 𝗶𝗻𝗱𝗲𝘅 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟱)
    |          V /Å³ |                  7.2032 |  ← 𝗣𝗼𝗹𝘆𝗵𝗲𝗱𝗿𝗮𝗹 𝘃𝗼𝗹𝘂𝗺𝗲
    |                |                         | 
    |           PP-5 |                 31.2864 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟱)
    |          vOC-5 |                  2.1704 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟱)
    |         TBPY-5 |                  3.1110 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟱)
    |          SPY-5 |                *0.9603* |  ← 𝗖𝗦𝗵𝗠 (𝗹𝗼𝘄𝗲𝘀𝘁 𝘃𝗮𝗹𝘂𝗲 𝗵𝗶𝗴𝗵𝗹𝗶𝗴𝗵𝘁𝗲𝗱)
    |        JTBPY-5 |                  6.3056 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟱)

Or to retrieve structural data from the COD (Crystallography Open Database), for example:
```console
python3 cshm-cc.py 4110517
```
The following output will be printed: 

    Fe1 (4110517): CN = 6, min dist. = 1.8997 Å, max dist. = 1.9875 Å  
    Co1 (4110517): Warning! Co1 (1_655) has been excluded from coordinating atoms.  
    Co1 (4110517): CN = 6, min dist. = 2.0824 Å, max dist. = 2.1582 Å 

    |   **compound** |   **Fe1 (4110517)** |   **Co1 (4110517)** |  ← 𝗖𝗲𝗻𝘁𝗿𝗮𝗹 𝗮𝘁𝗼𝗺𝘀 (𝗖𝗢𝗗 𝗜𝗗) 
    |----------------|---------------------|---------------------|
    |             CN |                   6 |                   6 |  ← 𝗖𝗼𝗼𝗿𝗱𝗶𝗻𝗮𝘁𝗶𝗼𝗻 𝗻𝘂𝗺𝗯𝗲𝗿𝘀 (𝗖𝗡)
    |              O |              5.1303 |              0.1262 |  ← 𝗚𝗲𝗼𝗺𝗲𝘁𝗿𝘆 𝗶𝗻𝗱𝗶𝗰𝗲𝘀 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟲)
    |          V /Å³ |              9.7103 |             12.5312 |  ← 𝗣𝗼𝗹𝘆𝗵𝗲𝗱𝗿𝗮𝗹 𝘃𝗼𝗹𝘂𝗺𝗲
    |                |                     |                     |
    |           HP-6 |             30.9116 |             33.2437 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟲)
    |          PPY-6 |             26.5838 |             30.2095 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟲)
    |           OC-6 |            *0.4782* |            *0.0261* |  ← 𝗖𝗦𝗵𝗠 (𝗹𝗼𝘄𝗲𝘀𝘁 𝘃𝗮𝗹𝘂𝗲 𝗵𝗶𝗴𝗵𝗹𝗶𝗴𝗵𝘁𝗲𝗱)
    |          TPR-6 |             14.1660 |             16.7105 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟲)
    |         JPPY-6 |             30.0476 |             33.6764 |  ← 𝗖𝗦𝗵𝗠 (𝘀𝘂𝗶𝘁𝗮𝗯𝗹𝗲 𝗳𝗼𝗿 𝗖𝗡 = 𝟲)

The central atom(s) or transition metal atom(s) and their coordination number(s) will be automatically determined, 
and all related values will be calculated. There is no manual choice of atom(s); only a maximum bonding distance can be set (`-d` option).

Start the script with:
```console
python3 cshm-cc.py example.cif > example.md
```
will save the output in markdown format.

Convert markdown to docx (install [PANDOC](https://pandoc.org) first):
```console
pandoc example.md -o example.docx
```
This will convert the markdown file to a docx file. Open it with your favorite
word processor. Convert the file to even more formats such as HTML, PDF or TeX with PANDOC.

## Command-line options
- `filename or COD ID` (required): A CIF (Crystallographic Information File) with the .cif extension or a COD ID (an integer without an extension), e.g., `example.cif` or `12345678`.
- `-d` `N` (optional): Excludes atoms with bond lengths larger than `N` Å from the central atom in the calculation (e.g., `-d 2.1`).
- `-eh` (optional): Exclude all hydrogen atoms.
- `-n` `N` (optional): Specifies the number of trials (`N` > 0) for the fast CShM calculation (default: 124).
- `-ex` (optional): Uses a slower CShM calculation that guarantees finding the global minimum.
- `-sxyz` (optional): Saves the XYZ coordinates of the central atom and its neighboring atoms. If multiple entries are provided, they are combined (filename: `cif_name.xyz` or `cod_id.xyz`).
- `-v` (optional): Enables verbose output, printing all bond lengths and angles of the central atom with its neighboring atoms.
- `-p` (optional): Plot bar graphs of CShM metrics for easier visual comparison.
- `-pc` (optional): Plot color bar graphs of CShM metrics for easier visual comparison.

## Remarks
- All parameter calculations are based on the estimated coordination number(s).
- Parameters like angles and distances are calculated from the atomic coordinates using the `gemmi` module. Therefore, symmetry operations for bonded atoms may differ from those in the CIF.
  Estimated standard deviations (e.s.d.) on bond lengths and angles are not calculated.
- The advantage of this method, compared to the one in the classic version, is that it now recognizes a much larger number of CIFs, especially those from the COD.
  A bonding section in the CIF file is no longer required.
- Possible bonds are now determined based on covalent radii. Several reported radii have been tested with varying success.
  In the end, the maximum radius of each element is used, with 10% added to the sum of the radii of the central and ligand atoms.
  This approach successfully identifies nearly all bonds to ligand atoms but also detects numerous metal-metal bonds,
  including those between the central atom and itself or other metals of the same element.
  These metal-metal bonds are excluded from further calculations.
- The central atom(s) or transition metal atom(s) and their coordination number(s) are determined automatically. There is no manual selection of atoms possible.
- Manual selection of atoms is possible with the script [tau-calc](https://github.com/radi0sus/tau-calc).
- The XYZ coordinates of neighboring atoms are provided relative to the central atom, which is positioned at [0, 0, 0].
- The XYZ file (option: `-sxyz`) can be used for further analysis of coordination geometry.
- The polyhedral volume should match the value calculated by [Olex2](https://www.olexsys.org/olex2/).
- Reference shapes for CShM are from
[here (cosymlib)](https://github.com/GrupEstructuraElectronicaSimetria/cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml).
- The script's results should match those of the [online calculator](https://csm.ouproj.org.il/molecule) or the [Shape program](https://www.ee.ub.edu/downloads/) when using default options.
- The CShM calculation employs a fast optimization process using the Hungarian algorithm. When only a small number of random trials is performed, the result may converge to a local minimum.
  For recommended number of trials, [see below](#regarding-the-number-of-trials).
- The `-ex` option enables a slower CShM calculation, which ensures finding the global minimum.
- The script downloads the CIF from the COD into a temporary folder. The CIF is deleted after the script finishes.

## Known Issues
- In rare cases, the `gemmi` module does not return the correct or expected coordinates for symmetry-equivalent positions. When this happens, the script exits. This should be fixed with gemmi >= 0.7.3.
- If problems arise, the only available option is to reduce the number of bonds considered using the `-d` option.
- The CShM method is rewritten from the [C++ code](https://github.com/continuous-symmetry-measure/shape) and may still contain errors.

## Regarding the number of trials
The CShM calculation employs a fast optimization process using the Hungarian algorithm. When only a small number of random trials are performed, the result may converge to a local minimum.
To determine the optimal number of trials, 100 runs were conducted, each with trial counts ranging from 1 to 100 (a total of 100 × 100 calculations, see [statistics.py](https://github.com/radi0sus/cshm-cc/blob/main/statistics.py)). The results are shown in the figure below.

<img src='examples\Figure_1.png' width='1024' alt='CShM vs. Number of Trials' align='center'>

The values for the global minimum solutions are displayed on the y-axis. Except for $\textcolor{blue}{\textrm{OC-6}}$ (the lowest CShM value) and $\textcolor{orange}{\textrm{HP-6}}$, the algorithm occasionally converges to local minima, which are higher than the global minimum. As seen from the color intensity of the points, the tendency to optimize local minima decreases significantly after approximately 20 trials and disappears completely after 40 trials. Thus, a "number of trials" of around 100 should be sufficient.

## Polyhedra and Shape Reference
<a href='https://creativecommons.org/licenses/by-nc-sa/4.0/'><img src='examples\all_polyhedra-6.png' alt='Polyhedra and Shape Reference' width=800 align='center'></a>
  
## Examples

### Example 1

```console
python3 cshm-cc.py combined.cif
```
    Ru1 (I): CN = 6, min dist. = 2.0645 Å, max dist. = 2.0656 Å  
    Cu1 (anko4): CN = 3, min dist. = 1.8636 Å, max dist. = 1.9991 Å  
    7718332: Warning! No metal atoms or '_atom_site_type_symbol' is missing.  
    Fe1 (av12_25): Warning! Fe1 (2_666) has been excluded from coordinating atoms.  
    Fe1 (av12_25): CN = 4, min dist. = 1.9890 Å, max dist. = 2.2038 Å  
    
    |   **compound** |   **Ru1 (I)** |   **Cu1 (anko4)** |   **Fe1 (av12_25)** |
    |----------------|---------------|-------------------|---------------------|
    |             CN |             6 |                 3 |                   4 |
    |             τ₄ |               |                   |              0.9035 |
    |            τ₄' |               |                   |              0.8976 |
    |             τ₅ |               |                   |                     |
    |              O |        7.1226 |                   |                     |
    |          V /Å³ |       11.5171 |            0.0415 |              4.5706 |
    |                |               |                   |                     |
    |           TP-3 |               |          *2.1792* |                     |
    |           vT-3 |               |            4.8766 |                     |
    |         fvOC-3 |               |           12.5001 |                     |
    |         mvOC-3 |               |            6.0484 |                     |
    |           SP-4 |               |                   |             29.3605 |
    |            T-4 |               |                   |            *1.2134* |
    |           SS-4 |               |                   |              8.3243 |
    |        vTBPY-4 |               |                   |              3.0150 |
    |           PP-5 |               |                   |                     |
    |          vOC-5 |               |                   |                     |
    |         TBPY-5 |               |                   |                     |
    |          SPY-5 |               |                   |                     |
    |        JTBPY-5 |               |                   |                     |
    |           HP-6 |       28.5605 |                   |                     |
    |          PPY-6 |       27.0546 |                   |                     |
    |           OC-6 |      *0.9516* |                   |                     |
    |          TPR-6 |       13.6312 |                   |                     |
    |         JPPY-6 |       30.8553 |                   |                     |

### Example 2

```console
python3 cshm-cc.py 1100756 -v -sxyz
```
    Fe1 (1100756): CN = 6, min dist. = 1.9321 Å, max dist. = 2.1414 Å  
    
    COD: 1100756, M: Fe1, CN: 6, Ligand atoms: O1  O2  O5  O7  N1  N2  
    
    |  **Atoms**  |   **Bond length /Å** |
    |-------------|----------------------|
    |   Fe1-O1    |               1.9394 |
    |   Fe1-O2    |               1.9321 |
    |   Fe1-O5    |               2.0936 |
    |   Fe1-O7    |               1.9975 |
    |   Fe1-N1    |               2.1147 |
    |   Fe1-N2    |               2.1414 | 
    
    |  **Atoms**  |   **Angle /°** |
    |-------------|----------------|
    |  O1-Fe1-O2  |          92.13 |
    |  O1-Fe1-O5  |         174.32 |
    |  O1-Fe1-O7  |          92.81 |
    |  O1-Fe1-N1  |          85.28 |
    |  O1-Fe1-N2  |          96.52 |
    |  O2-Fe1-O5  |          93.54 |
    |  O2-Fe1-O7  |          98.78 |
    |  O2-Fe1-N1  |         167.91 |
    |  O2-Fe1-N2  |          84.63 |
    |  O5-Fe1-O7  |          86.65 |
    |  O5-Fe1-N1  |          89.10 |
    |  O5-Fe1-N2  |          83.70 |
    |  O7-Fe1-N1  |          93.15 |
    |  O7-Fe1-N2  |         169.95 |
    |  N1-Fe1-N2  |          83.93 | 
    
    Fe2 (1100756): CN = 6, min dist. = 1.9019 Å, max dist. = 2.1397 Å  
    
    COD: 1100756, M: Fe2, CN: 6, Ligand atoms: O6  O3  N6  O4  N3  N4  
    
    |  **Atoms**  |   **Bond length /Å** |
    |-------------|----------------------|
    |   Fe2-O6    |               2.1298 |
    |   Fe2-O3    |               1.9019 |
    |   Fe2-N6    |               2.1397 |
    |   Fe2-O4    |               1.9100 |
    |   Fe2-N3    |               2.1042 |
    |   Fe2-N4    |               2.1254 | 
    
    |  **Atoms**  |   **Angle /°** |
    |-------------|----------------|
    |  O6-Fe2-O3  |          93.74 |
    |  O6-Fe2-N6  |          75.28 |
    |  O6-Fe2-O4  |         164.47 |
    |  O6-Fe2-N3  |          86.67 |
    |  O6-Fe2-N4  |          86.82 |
    |  O3-Fe2-N6  |          90.19 |
    |  O3-Fe2-O4  |          95.39 |
    |  O3-Fe2-N3  |          87.07 |
    |  O3-Fe2-N4  |         172.35 |
    |  N6-Fe2-O4  |          92.13 |
    |  N6-Fe2-N3  |         161.53 |
    |  N6-Fe2-N4  |          97.32 |
    |  O4-Fe2-N3  |         106.30 |
    |  O4-Fe2-N4  |          85.80 |
    |  N3-Fe2-N4  |          85.35 | 
    
    
    |   **compound** |   **Fe1 (1100756)** |   **Fe2 (1100756)** |
    |----------------|---------------------|---------------------|
    |             CN |                   6 |                   6 |
    |              O |              6.1844 |              9.2859 |
    |          V /Å³ |             11.1317 |             11.2190 |
    |                |                     |                     |
    |           HP-6 |             32.8092 |             32.9834 |
    |          PPY-6 |             25.3548 |             23.8228 |
    |           OC-6 |            *0.4960* |            *1.3337* |
    |          TPR-6 |             13.2723 |             10.3584 |
    |         JPPY-6 |             29.1733 |             27.1195 |
        
    XYZ file saved to 1100756.xyz

  XYZ file content:
  
    7
    1100756-Fe1
    Fe  0.00000000  0.00000000  0.00000000
    O  -0.61659181 -1.06174100 -1.50121523
    O   0.95997426  1.21200823 -1.15862156
    O   0.57158964  1.00798188  1.74371741
    O  -1.72418244  1.00058396  0.12684006
    N  -0.70442517 -1.59688382  1.19403794
    N   1.89597425 -0.97050214  0.22175439
    7
    1100756-Fe2
    Fe  0.00000000  0.00000000  0.00000000
    O  -0.48516594 -2.07125725  0.10297420
    O  -1.20716320  0.47924102  1.38941634
    N  -1.64079225 -0.24189918 -1.35192126
    O   0.31465477  1.78925625 -0.58946897
    N   1.41570551 -0.38239341  1.50903978
    N   1.54651075 -0.55071842 -1.34996022

### Example 3

```console
python3 cshm-cc.py 1100600 -p
```
    Cd1 (1100600): Warning! Cd2 (1_555) has been excluded from coordinating atoms.  
    Cd1 (1100600): Warning! Cd4 (1_555) has been excluded from coordinating atoms.  
    Cd1 (1100600): Warning! Cd3 (1_555) has been excluded from coordinating atoms.  
    Cd1 (1100600): CN = 4, min dist. = 2.2531 Å, max dist. = 2.3998 Å  
    Cd2 (1100600): Warning! Cd4 (1_555) has been excluded from coordinating atoms.  
    Cd2 (1100600): Warning! Cd1 (1_555) has been excluded from coordinating atoms.  
    Cd2 (1100600): Warning! Cd3 (1_555) has been excluded from coordinating atoms.  
    Cd2 (1100600): CN = 4, min dist. = 2.2682 Å, max dist. = 2.3968 Å  
    Cd3 (1100600): Warning! Cd2 (1_555) has been excluded from coordinating atoms.  
    Cd3 (1100600): Warning! Cd4 (1_555) has been excluded from coordinating atoms.  
    Cd3 (1100600): Warning! Cd1 (1_555) has been excluded from coordinating atoms.  
    Cd3 (1100600): CN = 4, min dist. = 2.2747 Å, max dist. = 2.3846 Å  
    Cd4 (1100600): Warning! Cd2 (1_555) has been excluded from coordinating atoms.  
    Cd4 (1100600): Warning! Cd1 (1_555) has been excluded from coordinating atoms.  
    Cd4 (1100600): Warning! Cd3 (1_555) has been excluded from coordinating atoms.  
    Cd4 (1100600): CN = 4, min dist. = 2.2334 Å, max dist. = 2.4005 Å    
    
    |   **compound** |   **Cd1 (1100600)** |   **Cd2 (1100600)** |   **Cd3 (1100600)** |   **Cd4 (1100600)** |
    |----------------|---------------------|---------------------|---------------------|---------------------|
    |             CN |                   4 |                   4 |                   4 |                   4 |
    |             τ₄ |              0.7378 |              0.7348 |              0.7131 |              0.7341 |
    |            τ₄' |              0.7273 |              0.7064 |              0.6972 |              0.7251 |
    |          V /Å³ |              5.5171 |              5.5015 |              5.5687 |              5.4022 |
    |                |                     |                     |                     |                     |
    |           SP-4 |             33.2281 |             31.6257 |             33.2891 |             33.0562 |
    |            T-4 |            *3.6158* |            *3.7698* |            *3.8506* |            *3.5041* |
    |           SS-4 |              7.5277 |              6.5626 |              7.2498 |              7.3999 |
    |        vTBPY-4 |              4.4485 |              4.5457 |              4.3667 |              4.1798 |
    
    Cd1 (1100600) CN = 4:  
    T-4    : ░░░░ 3.62  
    vTBPY-4: ░░░░░ 4.45  
    SS-4   : ░░░░░░░░░ 7.53  
    SP-4   : ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 33.23  
    
    Cd2 (1100600) CN = 4:  
    T-4    : ░░░░ 3.77  
    vTBPY-4: ░░░░░ 4.55  
    SS-4   : ░░░░░░░ 6.56  
    SP-4   : ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 31.63  
    
    Cd3 (1100600) CN = 4:  
    T-4    : ░░░░ 3.85  
    vTBPY-4: ░░░░░ 4.37  
    SS-4   : ░░░░░░░░ 7.25  
    SP-4   : ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 33.29  
    
    Cd4 (1100600) CN = 4:  
    T-4    : ░░░░ 3.50  
    vTBPY-4: ░░░░░ 4.18  
    SS-4   : ░░░░░░░░ 7.40  
    SP-4   : ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 33.06 

Or color bar graphs using:

```console
python3 cshm-cc.py 1100600 -pc
```
    
<img src='examples\color_bars_2.png' width='460' alt='Color bar graphs' align='center'>

### Example 4

<img src='examples\show-use.gif' alt='Example 4' width=1024 align='center'>

## References
If you use τ<sub>4</sub>, τ<sub>5</sub>, the *O* index or CShM to describe the coordination geometry of your compounds, please cite one or more of the following articles:

**τ<sub>4</sub>**:
> "Structural variation in copper(i) complexes with pyridylmethylamide ligands: 
>  structural analysis with a new four-coordinate geometry index, τ<sub>4</sub>"
>  
> Lei Yang, Douglas R. Powell, Robert P. Houser,
> *Dalton Trans.* **2007**, 955-964.
> 
> DOI: https://doi.org/10.1039/B617136B

**τ<sub>4</sub>' (τ<sub>4</sub> improved)**:
> "Coordination polymers and molecular structures among complexes of 
>  mercury(II) halides with selected 1-benzoylthioureas"
> 
> Andrzej Okuniewski, Damian Rosiak, Jarosław Chojnacki, Barbara Becker,
> *Polyhedron* **2015**, *90*, 47–57.
> 
> DOI: https://doi.org/10.1016/j.poly.2015.01.035

**τ<sub>5</sub>**:
> "Synthesis, structure, and spectroscopic properties of copper(II) compounds containing 
>  nitrogen–sulphur donor ligands; the crystal and molecular structure of 
>  aqua[1,7-bis(N-methylbenzimidazol-2′-yl)-2,6-dithiaheptane]copper(II) perchlorate"
>  
> Anthony W. Addison, T. Nageswara Rao, Jan Reedijk, Jacobus van Rijn, Gerrit C. Verschoor, 
> *J. Chem. Soc., Dalton Trans.* **1984**, 1349-1356.
> 
> DOI: https://doi.org/10.1039/DT9840001349

***O***:
> "Structural, electrochemical and photophysical behavior of Ru(II) complexes with 
>  large bite angle sulfur-bridged terpyridyl ligands"
>  
> Christopher M. Brown, Nicole E. Arsenault, Trevor N. K. Cross, Duane Hean, Zhen Xu, 
> Michael O. Wolf, 
> *Inorg. Chem. Front.* **2020**, *7*, 117-127.
> 
> DOI: https://doi.org/10.1039/C9QI01009B

**CShM**:
> "Continuous Symmetry Measures. 5. The Classical Polyhedra"
>  
> Mark Pinsky, David Avnir, 
> *Inorg. Chem.* **1998**, *37*, 5575–5582.
> 
> DOI: https://doi.org/10.1021/ic9804925
> 
> "Shape maps and polyhedral interconversion paths in transition metal chemistry"
>  
> Santiago Alvarez, Pere Alemany, David Casanova, Jordi Cirera, Miquel Llunell, David Avnir,
> *Coord. Chem. Rev.*, **2005**, *249*, 1693–1708.
> 
> DOI: https://doi.org/10.1016/j.ccr.2005.03.031

The script uses the **Gemmi** library for CIF processing:
> "GEMMI: A library for structural biology"
> 
> Marcin Wojdyr,
> *Journal of Open Source Software* **2022**, *7 (73)*, 4200.
>
> DOI: https://doi.org/10.1021/ic9804925

https://gemmi.readthedocs.io/en/latest/

https://github.com/project-gemmi/gemmi

The **Crystallography Open Database**:

https://www.crystallography.net/cod/
