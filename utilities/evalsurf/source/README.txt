This program is a "stand alone" for generating energies, gradients, and
derivative couplings from surfgen.x input files:
 hd.data
 irrep.in
 coord.in
 surfgen.in

Please refer to surfgen.x documentation for information on these files.

Requires openmp libraries.

Build with:
make FC={Fortran compiler} BLAS_LIB={BLAS libraries + openmp libraries}

Run:
evalsurf.x [input geometry file]
