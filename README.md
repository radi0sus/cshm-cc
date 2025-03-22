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
