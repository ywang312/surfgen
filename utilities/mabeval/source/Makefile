MABEVALOBJS := ioutil.o mabutil.o mabeval.o
INTCDISPOBJS := ioutil.o cart2int_util.o mabutil.o intc_mabutil.o intc_displacements.o
INTCMABFOBJS := ioutil.o cart2int_util.o mabutil.o intc_mabutil.o intc_mab_f.o
INTCMABEVALOBJS := ioutil.o cart2int_util.o mabutil.o intc_mabutil.o intc_mabeval.o
GHZPOSITIONOBJS := ioutil.o ghz_position.o
GHZDISPLACEMENTOBJS := ioutil.o ghz_displacement.o
COMPGDISTOBJS := ioutil.o compute_geomdist.o
MKGHINTCOBJS := ioutil.o cart2int_util.o mabutil.o intc_mabutil.o mkgh_intc.o
COLIBX:=colib.a
UNIXDIR:=./UNIX
COLIBDIR:=./colib

# Objects for colib library 
COLIBO:=blaswrapper.o colib1.o colib2.o colib3.o colib4.o colib5.o colib6.o \
	colib7.o colib8.o colib9.o colib10.o
# Source files for colib lirary
COLIBF:=blaswrapper.f colib1.f colib2.f colib3.f colib4.f colib5.f colib6.f \
	colib7.f colib8.f colib9.f colib10.f 
UNIXO := fdate.o falloc.o fwtime.o hostnm.o flushstdout.o fsize.o
UNIXC := fdate.c falloc.c fwtime.c hostnm.c flushstdout.c fsize.c
COLIBOBJS:= $(addprefix $(COLIBDIR)/,$(COLIBO))
COLIBSRCF:= $(addprefix $(COLIBDIR)/,$(COLIBF))
UNIXOBJS := $(addprefix $(UNIXDIR)/,$(UNIXO))
UNIXSRCC := $(addprefix $(UNIXDIR)/,$(UNIXC))

#use predefined fortran compiler
ifdef FC 
	FC := $(FC)
else
	FC:=gfortran
	CC:=cc
endif

ifneq ($(findstring ifort, $(shell $(FC) --version)),)
	FFLAGS:=-g -i8 -O3 -auto -assume byterecl -traceback
	CC:=icc
else
	FFLAGS:=-g -fdefault-integer-8 -frecord-marker=4 -O3 -fbacktrace \
		-fno-align-commons
endif

# Cflags are the same
CFLAGS:=-Wall -std=c11 -march=native -funroll-loops -ffast-math \
	          -fomit-frame-pointer -fstrict-aliasing -O3
CCDEF  := -DFLUSH -DINT64 -DEXTNAME
AR:= ar rv
RANL:= ranlib
RM:=rm -rf

all: colib intc_mab_f ghz_position ghz_displacement compute_geomdist mabeval mkgh_intc

mabeval: $(MABEVALOBJS)
	$(FC) $(FFLAGS) -o mabeval.x $(MABEVALOBJS) $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp
	cp mabeval.x ../bin/

intc_displacements: colib $(INTCDISPOBJS)
	$(FC) $(FFLAGS) -o intc_displacements.x $(INTCDISPOBJS) $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp ./colib.a
	cp intc_displacements.x ../bin/

intc_mabeval: colib $(INTCMABEVALOBJS)
	$(FC) $(FFLAGS) -o intc_mabeval.x $(INTCMABEVALOBJS) $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp ./colib.a
	cp intc_mabeval.x ../bin/

intc_mab_f: colib $(INTCMABFOBJS)
	$(FC) $(FFLAGS) -o intc_mab_f.x $(INTCMABFOBJS) $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp ./colib.a
	cp intc_mab_f.x ../bin/

ghz_position: $(GHZPOSITIONOBJS)
	$(FC) $(FLAGS) -o ghz_position.x $(GHZPOSITIONOBJS) $(BLAS_LIB) -lgomp
	cp ghz_position.x ../bin/

ghz_displacement: $(GHZDISPLACEMENTOBJS)
	$(FC) $(FLAGS) -o ghz_displacement.x $(GHZDISPLACEMENTOBJS) $(BLAS_LIB) -lgomp
	cp ghz_displacement.x ../bin/

compute_geomdist: $(COMPGDISTOBJS)
	$(FC) $(FLAGS) -o compute_geomdist.x $(COMPGDISTOBJS) $(BLAS_LIB) -lgomp
	cp compute_geomdist.x ../bin/

mkgh_intc: colib $(MKGHINTCOBJS)
	$(FC) $(FLAGS) -o mkgh_intc.x $(MKGHINTCOBJS) $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp ./colib.a
	cp mkgh_intc.x ../bin/

colib: $(COLIBOBJS) $(UNIXOBJS)
	$(AR) $(COLIBX) $(COLIBOBJS) $(UNIXOBJS)
	$(RANL) $(COLIBX)
	chmod +x $(COLIBX)

clean:
	rm -rf ./colib/colib*.o ./*.o ./*.mod ./*.x ./*.a

$(UNIXDIR)/%.o: $(UNIXDIR)/%.c
	$(CC) -c -DFLUSH -DINT64 -DEXTNAME -o $@ $< $(CFLAGS)

$(COLIBDIR)/%.o: $(COLIBDIR)/%.f
	$(FC) -g -c -o $@ $< $(CPOPS) $(FFLAGS)


%.o:%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o:%.f
	$(FC) $(FFLAGS) -c -o $@ $<

