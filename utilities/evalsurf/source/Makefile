all:
	gfortran -c evalsurf.f90
	gfortran -o evalsurf.x evalsurf.o libsurfgen.a $(BLAS_LIB) -lgomp 
