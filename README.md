# cshm-cc

A Python 3 script that calculates CShM (Continuous Shape Measures) values for various shapes, geometry indices, including τ₄, τ₄', τ₅, *O* (octahedricity) of automatically selected atoms (transition metal atoms by default) from a crystallographic information file (CIF). The CIF may contain one or more entries. It can also calculate these values from COD (Crystallography Open Database) entries by simply entering the COD ID. It also saves the XYZ coordinates of the central atom and its neighboring atoms, including those from CIFs with multiple entries or COD entries and calculates the polyhedral volume. Optionally, tables in Markdown format containing bond lengths and angles can be generated.

## Introduction
This script calculates the Continuous Shape Measures (CShM) and geometry indices (τ₄, τ₄', τ₅, and *O* (octahedricity)) to assign coordination geometries to three-, four-, five and six-coordinate transition metal atoms. The indices τ₄ and τ₅ are used to determine whether a compound adopts tetrahedral, trigonal pyramidal, square planar, or seesaw geometry (for four-coordinated compounds), and square pyramidal or trigonal bipyramidal geometry (for five-coordinated compounds). These assignments rely on the two largest angles enclosing the central atom. The octahedricity index, *O*, is calculated based on experimental X-M-X angles and assesses how close the geometry is to an ideal octahedron. The CShM value approaches zero when the shape closely matches the ideal geometry. Additionally, the polyhedral volume is calculated using the Convex Hull Algorithm (via SciPy and the Qhull Library). For more information, refer to the associated papers and references.

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
- The XYZ coordinates of neighboring atoms are provided relative to the central atom, which is positioned at [0, 0, 0].
- The XYZ file (option: `-sxyz`) can be used for further analysis of coordination geometry.
- The polyhedral volume should match the value calculated by [Olex2](https://www.olexsys.org/olex2/).
- Reference shapes for CShM are from
[here (cosymlib)](https://github.com/GrupEstructuraElectronicaSimetria/cosymlib/blob/master/cosymlib/shape/ideal_structures_center.yaml).
- The script's results should match those of the [online calculator](https://csm.ouproj.org.il/molecule) or the [Shape program](https://www.ee.ub.edu/downloads/) when using default options.
- The CShM calculation employs a fast optimization process using the Hungarian algorithm. When only a small number of random trials is performed, the result may converge to a local minimum.
- For recommended number of trials, see below.
- The `-ex` option enables a slower CShM calculation, which ensures finding the global minimum.

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
The CShM calculation employs a fast optimization process using the Hungarian algorithm. When only a small number of random trials is performed, the result may converge to a local minimum.
To determine the optimal number of trials, 100 runs were conducted with trial counts ranging from 1 to 100. The results are shown in the figure below:

<img src='examples\Figure_1.png' width='1024' alt='CShM vs. Number of Trials' align='center'>

The values for the global minimum solutions are displayed on the y-axis. Except for $\textcolor{blue}{\textrm{OC-6}}$ (the lowest CShM value) and $\textcolor{orange}{\textrm{HP-6}}$, the algorithm occasionally converges to local minima, which are higher than the global minimum. As seen from the color intensity of the points, the tendency to optimize local minima decreases significantly after approximately 20 trials and disappears completely after 40 trials. Thus, a "number of trials" of around 100 should be sufficient.

## Polyhedra and Shape Reference

<a href='https://creativecommons.org/licenses/by-nc-sa/4.0/'><img src='examples\all_polyhedra5.png' alt='Polyhedra and Shape Reference' width=800 align='center'></a>
  
## Examples

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
