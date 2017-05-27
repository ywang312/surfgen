
OBJECTS=evalsurf.o
LIBOBJS=hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o \
	io.o potlib.o getver.o

FC:=gfortran
FCOPT:=-g

all:	lib $(OBJECTS)
	$(FC) $(FCOPT) -c evalsurf.f90
	$(FC) $(FCOPT) -o evalsurf.x evalsurf.o libsurfgen.a $(BLAS_LIB) 

lib:	$(LIBOBJS)
	$(AR) -r -v libsurfgen.a $(LIBOBJS)
	@echo " Cleaning up..."
	rm -f $(LIBOBJS)
	rm -f ./*.mod

./%.o : ./lib/%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./lib/%.F90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

clean :
	rm -f ./*.mod ./*.o
