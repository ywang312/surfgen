all: ioutil.o mabutil.o mabeval.o
	gfortran -g -o mabeval.x ioutil.o mabutil.o mabeval.o $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp
	cp mabeval.x ../bin/

ghz_position: ioutil.o ghz_position.o
	gfortran -g -o ghz_position.x ioutil.o ghz_position.o $(BLAS_LIB) -lgomp
	cp ghz_position.x ../bin/

compute_geomdist: ioutil.o compute_geomdist.o
	gfortran -g -o compute_geomdist.x ioutil.o compute_geomdist.o $(BLAS_LIB) -lgomp
	cp compute_geomdist.x ../bin/

clean:
	rm -rf ./*.o ./*.mod

%.o:%.f90
	gfortran -g -fbacktrace -c -o $@ $<
