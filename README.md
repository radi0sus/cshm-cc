# cshm-cc

A Python 3 script that calculates CShM (Continuous Shape Measures) values for various shapes, geometry indices, including τ₄, τ₄', τ₅, *O* (octahedricity) of selected atoms from a crystallographic information file (CIF). The CIF may contain one or more entries. It can easily calculate these values from COD (Crystallography Open Database) entries by simply entering the COD number. It also saves the XYZ coordinates of the central atom and its neighboring atoms, including those from CIFs with multiple entries or COD entries and calculates the polyhedral volume. Optionally, tables in Markdown format containing bond lengths and angles can be generated.

## Introduction

This script calculates the Continuous Shape Measures (CShM) and geometry indices (τ₄, τ₄', τ₅, and *O* (octahedricity)) to assign coordination geometries to three-, four-, five and six-coordinate compounds. The indices τ₄ and τ₅ are used to determine whether a compound adopts tetrahedral, trigonal pyramidal, square planar, or seesaw geometry (for four-coordinated compounds), and square pyramidal or trigonal bipyramidal geometry (for five-coordinated compounds). These assignments rely on the two largest angles enclosing the central atom. The octahedricity index, *O*, is calculated based on experimental X-M-X angles and assesses how close the geometry is to an ideal octahedron. The CShM value approaches zero when the shape closely matches the ideal geometry. Additionally, the polyhedral volume is calculated using the Convex Hull Algorithm (via SciPy and the Qhull Library). For more information, refer to the associated papers and references.

## External modules
 `gemmi`, `numpy`, `scipy`, `tabulate` 
 
## Quick start
For local CIFs with one or more entries start the script with:
```console
python3 cshm-cc.py filename.cif
```
or to retrieve structural data from the COD (Crystallography Open Database), for example:
```console
python3 cshm-cc.py 1546808
```
 
