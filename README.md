# symplectic-basis
A Symplectic Basis for Triangulated 3-Manifolds. 

https://arxiv.org/abs/2208.06969

## Setup

Install SnapPy and symplectic-basis

> python -m pip install snappy symplectic-basis
> 
> python -m snappy.app

In SnapPy, import symplectic-basis. 

> import symplectic_basis
> 
> M = Manifold('4_1')
> 
> symplectic-basis.symplectic_basis(M)

## Developement

``` symp_src ``` contains the source files for computing oscillating curves.
``` snappea_src ``` contains the full SnapPea kernel source files for running ``` symp_src ``` directly in C using the CMakeLists.txt.
``` snappy_src ``` contains a reduced SnapPea kernel which provides the minimal amount of source code to run ``` read_triangulation_from_string ``` which is used to get the C Triangulation struct from the SnapPy manifold python object.
It also contains ``` peripheral_curves.c ``` and ``` get_gluing_equations.c ``` which are used to construct the symplectic basis.

## Credits

```snappy_src``` and ```snappea_src``` contain parts of the SnapPy source code (M. Culler, N. M. Dunfield, M. Goerner, and J. R. Weeks, SnapPy, a computer program for studying the geometry and topology of 3-manifolds, http://snappy.computop.org)
