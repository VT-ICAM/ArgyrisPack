ArgyrisPack
===========
Overview
--------
ArgyrisPack is a library for doing things with the Argyris-5 finite element, a
C1 finite element with 21 degrees of freedom. The Argyris-5 element spans the
space of 5th degree polynomials on each element and has O(h^4) convergence in
the H2 norm. The main goals of ArgyrisPack are numerical correctness,
performance, and portability. ArgyrisPack is mostly written in C and can be called
from C, Matlab, Julia, or Python.

ArgyrisPack is available under the BSD-3 License. See the file `License.txt` for
more details.

Design & Usage Philosophy
-------------------------
ArgyrisPack assumes that it is faster to work with precomputed values of
reference basis functions. Therefore it is up to the user to set up appropriate
tables of reference data (function values, derivative values, presumably second
derivative values) and to use these when setting up the typical per-element
matricies.

Since performance and portability are both goals, the C code does not offer
much abstraction; one passes in pointers to arrays and the functions fill the
arrays appropriately. The non-C interfaces take care of wrapping output values
and checking bounds instead.

Compilation Details
-------------------
ArgyrisPack is written in C99 and may be compiled to use either row-major or
column-major storage.  Macros (contained in `order_logic.h`) ensure that the C
code does 'the right things' (correct BLAS calls, correct order traversal, etc)
in both cases. This is done by passing either `-DUSE_ROW_MAJOR` or
`-DUSE_COL_MAJOR` to the compiler. In particular, with the makefile, compile the
library with

    make STORAGE_ORDER=USE_ROW_MAJOR

to pass one symbol or another. To compile the mex-files, run

    make

at the MATLAB prompt. MATLAB (to the best of my knowledge) assumes MEX files on
*nix platforms will use the gcc, but you may need to tweak your local MEX
settings to compile the binaries correctly. Julia just needs `libargyris_pack.so`
in the current path to work. Python needs for the .so to be in the working path
as well.

ArgyrisPack currently targets GNU/Linux and OS X. GCC 4.2 and later should work;
please let us know if this is not the case. In principle, there is nothing
preventing a port to other operating systems, but all the provided makefiles
assume a *nix environment.

What can ArgyrisPack do?
------------------------
ArgyrisPack consists of four levels of C and Python functions/classes:
* Level 1: evaluation of the Argyris basis functions (and derivatives) at
  arbitrary points on the reference triangle. These are the functions beginning
  with `ref_`.
* Level 2: evaluation of the required coordinate transformations (Dominguez and
  Sayas' C, B, b, and Theta (Th) matricies) and evaluation of Argyris basis
  functions on arbitrary triangles. These are the functions beginning with
  `physical_`.
* Level 3: evaluation of a few common bilinear forms (the classic stiffness and
  mass matricies as well as the 'biharmonic' matrix resulting from
  discretization of the biharmonic operator). These functions begin with
  `matrix_`.
* Level 4: mesh generation. Argyris elements have 21 nodes (5 on each corner
  and one at the midpoint of each triangle edge). ArgyrisPack contains mesh
  parsing and creation classes for a variety of textual representations of
  meshes.

There are a few additional files; we wrote a 'multiply by a diagonal matrix'
routine, wrappers to make `dgemm` work with row or column order, as well as a
symbolic (`symbolic.py` and `symbolic.m`) version of the Argyris element to
verify the numerical computations.

Running the Tests
-----------------
The tests are currently implemented in Python. They require SAGE to be installed
and used as a library. This is because the tests work by interpolating the
Argyris polynomials on a physical element and then comparing the symbolic
results with the numerical results. To run the results from SAGE, go to the root
directory of ArgyrisPack and execute

    import ap.test; ap.test.run()

The meshing software has tests, but these are not linked to the numerical tests
at the moment.

Authors
-------
* Erich Foster, Virginia Tech
* Traian Iliescu, Virginia Tech
* David Wells, Virginia Tech

Relevant Papers
---------------
ArgyrisPack is a new implementation of older ideas. We are greatly indebted
to those who have done the necessary mathematical research for making this
library possible.

* Dominguez and Sayas: A simple Matlab implementation of the Argyris element,
  2006
* Argyris, Fried and Scharpf: The TUBA Family of Plate Elements for the Matrix
  Displacement Method, 1968.
