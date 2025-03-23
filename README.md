# cshm-cc

A Python 3 script that calculates CShM (Continuous Shape Measures) values for various shapes, geometry indices, including τ₄, τ₄', τ₅, *O* (octahedricity) of automatically selected atoms (transition metal atoms by default) from a crystallographic information file (CIF). The CIF may contain one or more entries. It can also calculate these values from COD (Crystallography Open Database) entries by simply entering the COD ID. It also saves the XYZ coordinates of the central atom and its neighboring atoms, including those from CIFs with multiple entries or COD entries and calculates the polyhedral volume. Optionally, tables in Markdown format containing bond lengths and angles can be generated.

## Introduction
This script calculates the Continuous Shape Measures (CShM) and geometry indices (τ₄, τ₄', τ₅, and *O* (octahedricity)) to assign coordination geometries to three-, four-, five and six-coordinate transition metal atoms. The indices τ₄ and τ₅ are used to determine whether a compound adopts tetrahedral, trigonal pyramidal, square planar, or seesaw geometry (for four-coordinated compounds), and square pyramidal or trigonal bipyramidal geometry (for five-coordinated compounds). These assignments rely on the two largest angles enclosing the central atom. The octahedricity index, *O*, is calculated based on experimental X-M-X angles and assesses how close the geometry is to an ideal octahedron. The CShM value approaches zero when the shape closely matches the ideal geometry. Additionally, the polyhedral volume is calculated using the Convex Hull Algorithm (via SciPy and the Qhull Library). For more information, refer to the associated papers and references.

#### Equations and values for the geometry indices τ₄, τ₄', τ₅, and *O*:

$\tau_4 = \frac{360°-(\alpha+\beta)}{141°}$ 

Square planar geometry: $\tau_4 = 0$; Tetrahedral geometry:  $\tau_4 = 1$; Seesaw geometry: $\tau_4 \approx 0.43$

$\tau'_4 = \frac{\beta-\alpha}{250.5°}+\frac{180-\beta}{70.5°}$

Square planar geometry: $\tau'_4 = 0$; Tetrahedral geometry:  $\tau'_4 = 1$; Seesaw geometry: $\tau_4 \approx 0.24$

$\tau_5 = \frac{\beta-\alpha}{60°}$

Square pyramidal geometry: $\tau_5 = 0$; Trigonal bipyramidal geometry: $\tau_5 = 1$

$O = \sqrt{\frac{1}{15}\sum_{i=1}^{15}(\hat{\theta_i} - \theta)^2}$

$\hat{\theta_i}$ = 180° for *trans* X-M-X angles and 90° for *cis* X-M-X angles\
$\theta$ = experimental X-M-X angles

$O$ is close to zero for an almost ideal octahedron

## External modules
 `gemmi`, `numpy`, `scipy`, `tabulate`, `requests` 
 
## Quick start
For local CIFs with one or more entries start the script with:
```console
python3 cshm-cc.py example.cif
```
The following output will be printed:

 
    |   **compound** |   **Fe1 (example)**     |   <- Central atom (CIF entry) 
    |----------------|-------------------------|
    |             CN |                       5 |   <- Coordination number (CN)
    |             τ₅ |                  0.3053 |   <- Geometry index (suitable for CN = 5)
    |          V /Å³ |                  7.2032 |   <- Polyhedral volume
    |                |                         | 
    |           PP-5 |                 31.2864 |   <- CShM (suitable for CN = 5)
    |          vOC-5 |                  2.1704 |   <- CShM (suitable for CN = 5)
    |         TBPY-5 |                  3.1110 |   <- CShM (suitable for CN = 5)
    |          SPY-5 |                *0.9603* |   <- CShM (lowest value highlighted)
    |        JTBPY-5 |                  6.3056 |   <- CShM (suitable for CN = 5)

Or to retrieve structural data from the COD (Crystallography Open Database), for example:
```console
python3 cshm-cc.py 4110517
```
The following output will be printed:

    |   **compound** |   **Fe1 (4110517)** |   **Co1 (4110517)** |   <- Central atoms (COD ID)  
    |----------------|---------------------|---------------------|
    |             CN |                   6 |                   6 |   <- Coordination numbers (CN)
    |              O |              5.1303 |              0.1262 |   <- Geometry indices (suitable for CN = 6)
    |          V /Å³ |              9.7103 |             12.5312 |   <- Polyhedral volumes
    |                |                     |                     |
    |           HP-6 |             30.9116 |             33.2437 |   <- CShM (suitable for CN = 6)
    |          PPY-6 |             26.5838 |             30.2095 |   <- CShM (suitable for CN = 6)
    |           OC-6 |            *0.4782* |            *0.0261* |   <- CShM (lowest values highlighted)
    |          TPR-6 |             14.1660 |             16.7105 |   <- CShM (suitable for CN = 6)
    |         JPPY-6 |             30.0476 |             33.6764 |   <- CShM (suitable for CN = 6)

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
- `-n` `N` (optional): Specifies the number of trials (`N` > 0) for the fast CShM calculation (default: 124).
- `-ex` (optional): Uses a slower CShM calculation that guarantees finding the global minimum.
- `-sxyz` (optional): Saves the XYZ coordinates of the central atom and its neighboring atoms. If multiple entries are provided, they are combined (filename: `cif_name.xyz` or `cod_id.xyz`).
- `-v` (optional): Enables verbose output, printing all bond lengths and angles of the central atom with its neighboring atoms.

## Remarks
- All parameter calculations are based on the estimated coordination number(s).
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
- The script determines ligand atoms based on the bonding section in the CIF file; therefore, this section must be present.
- Additionally, all bond angles involving the central atom must be included in the CIF file.
- The script may fail if bonding information is incomplete; for example, if symmetry-equivalent positions are missing or if manual entries have been made.
- Metal-metal bonding can also cause issues in some cases.
- If problems arise, the only available option is to reduce the number of bonds considered using the `-d` option.
- The script can only remove atoms from the coordination sphere, not add them. Ensure that the connectivity list is appropriate.
- Hydrogen atoms are generally ignored. As a result, hydrogen atoms in metal hydrides will not be considered.
- The CShM method is rewritten from the [C++ code](https://github.com/continuous-symmetry-measure/shape) and may still contain errors.

## Regarding the number of trials

The CShM calculation employs a fast optimization process using the Hungarian algorithm. When only a small number of random trials are performed, the result may converge to a local minimum.
To determine the optimal number of trials, 100 runs were conducted, each with trial counts ranging from 1 to 100 (a total of 100 × 100 calculations). The results are shown in the figure below.

<img src='examples\Figure_1.png' width='1024' alt='CShM vs. Number of Trials' align='center'>

The values for the global minimum solutions are displayed on the y-axis. Except for $\textcolor{blue}{\textrm{OC-6}}$ (the lowest CShM value) and $\textcolor{orange}{\textrm{HP-6}}$, the algorithm occasionally converges to local minima, which are higher than the global minimum. As seen from the color intensity of the points, the tendency to optimize local minima decreases significantly after approximately 20 trials and disappears completely after 40 trials. Thus, a "number of trials" of around 100 should be sufficient.

## Polyhedra and Shape Reference

<a href='https://creativecommons.org/licenses/by-nc-sa/4.0/'><img src='examples\all_polyhedra5.png' alt='Polyhedra and Shape Reference' width=800 align='center'></a>
  
## Examples

### Example 1

```console
python3 cshm-cc.py combined.cif
```
    |   **compound** |   **Ru1 (I)** |   **Cu1 (anko4)** |   **Fe1 (av12_25)** |
    |----------------|---------------|-------------------|---------------------|
    |             CN |             6 |                 3 |                   4 |
    |             τ₄ |               |                   |              0.9035 |
    |            τ₄' |               |                   |              0.8976 |
    |              O |        7.1201 |                   |                     |
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
    |           HP-6 |       28.5605 |                   |                     |
    |          PPY-6 |       27.0546 |                   |                     |
    |           OC-6 |      *0.9516* |                   |                     |
    |          TPR-6 |       13.6312 |                   |                     |
    |         JPPY-6 |       30.8553 |                   |                     |

### Example 2

```console
python3 cshm-cc.py 1100756 -v -sxyz
```
    COD: 1100756, M: Fe1, CN: 6, Ligand atoms: O1 O2 O5 O7 N1 N2 
    
    |  **Atoms**  |   **Bond length /Å** |
    |-------------|----------------------|
    |   Fe1-O2    |           1.9321(17) |
    |   Fe1-O1    |           1.9393(18) |
    |   Fe1-O7    |             1.998(2) |
    |   Fe1-O5    |           2.0937(18) |
    |   Fe1-N1    |             2.115(2) |
    |   Fe1-N2    |             2.141(2) | 
    
    |  **Atoms**  |   **Angle /°** |
    |-------------|----------------|
    |  O2-Fe1-O1  |       92.13(7) |
    |  O2-Fe1-O7  |       98.79(8) |
    |  O1-Fe1-O7  |       92.81(8) |
    |  O2-Fe1-O5  |       93.54(7) |
    |  O1-Fe1-O5  |      174.31(7) |
    |  O7-Fe1-O5  |       86.65(8) |
    |  O2-Fe1-N1  |      167.89(8) |
    |  O1-Fe1-N1  |       85.28(8) |
    |  O7-Fe1-N1  |       93.16(9) |
    |  O5-Fe1-N1  |       89.09(8) |
    |  O2-Fe1-N2  |       84.63(8) |
    |  O1-Fe1-N2  |       96.52(8) |
    |  O7-Fe1-N2  |      169.95(8) |
    |  O5-Fe1-N2  |       83.70(8) |
    |  N1-Fe1-N2  |       83.92(9) | 
    
    
    COD: 1100756, M: Fe2, CN: 6, Ligand atoms: O6 O3 N6 O4 N3 N4 
    
    |  **Atoms**  |   **Bond length /Å** |
    |-------------|----------------------|
    |   Fe2-O3    |           1.9020(17) |
    |   Fe2-O4    |           1.9099(18) |
    |   Fe2-N3    |             2.104(2) |
    |   Fe2-N4    |             2.126(2) |
    |   Fe2-O6    |           2.1298(18) |
    |   Fe2-N6    |             2.140(2) | 
    
    |  **Atoms**  |   **Angle /°** |
    |-------------|----------------|
    |  O3-Fe2-O4  |       95.39(7) |
    |  O3-Fe2-N3  |       87.08(8) |
    |  O4-Fe2-N3  |      106.30(8) |
    |  O3-Fe2-N4  |      172.34(8) |
    |  O4-Fe2-N4  |       85.80(8) |
    |  N3-Fe2-N4  |       85.33(8) |
    |  O3-Fe2-O6  |       93.73(7) |
    |  O4-Fe2-O6  |      164.47(7) |
    |  N3-Fe2-O6  |       86.67(7) |
    |  N4-Fe2-O6  |       86.82(7) |
    |  O3-Fe2-N6  |       90.19(8) |
    |  O4-Fe2-N6  |       92.14(8) |
    |  N3-Fe2-N6  |      161.53(8) |
    |  N4-Fe2-N6  |       97.33(8) |
    |  O6-Fe2-N6  |       75.28(7) | 
    
    |   **compound** |   **Fe1 (1100756)** |   **Fe2 (1100756)** |
    |----------------|---------------------|---------------------|
    |             CN |                   6 |                   6 |
    |              O |              6.1891 |              9.2875 |
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
