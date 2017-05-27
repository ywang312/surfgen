INTCMABFOBJS := ioutil.o cart2int_util.o intc_mabutil.o intc_mab_f.o
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

FC:=gfortran
FFLAGS:=-g -fdefault-integer-8 -frecord-marker=4 -O3
CC:=cc
CFLAGS:=-Wall -std=c11 -march=native -funroll-loops -ffast-math \
	          -fomit-frame-pointer -fstrict-aliasing -O3
CCDEF  := -DFLUSH -DINT64 -DEXTNAME
AR:= ar rv
RANL:= ranlib
RM:=rm -rf

all: colib intc_mab_f

intc_mab_f: colib $(INTCMABFOBJS)
	$(FC) $(FFLAGS) -fbacktrace -o intc_mab_f.x $(INTCMABFOBJS) $(SURFGEN)/lib/libsurfgen.a $(BLAS_LIB) -lgomp ./colib.a

colib: $(COLIBOBJS) $(UNIXOBJS)
	$(AR) $(COLIBX) $(COLIBOBJS) $(UNIXOBJS)
	$(RANL) $(COLIBX)
	chmod +x $(COLIBX)

clean:
	rm -rf ./colib/colib*.o ./*.o ./*.mod ./*.x ./*.a

$(UNIXDIR)/%.o: $(UNIXDIR)/%.c
	$(CC) -c -DFLUSH -DINT64 -DEXTNAME -o $@ $< $(CFLAGS) -Wno-implicit-function-declaration

$(COLIBDIR)/%.o: $(COLIBDIR)/%.f
	$(FC) -g -fbacktrace -c -o $@ $< $(CPOPS) $(FFLAGS)


%.o:%.f90
	$(FC) $(FFLAGS) -fbacktrace -c -o $@ $<

%.o:%.f
	$(FC) $(FFLAGS) -fbacktrace -fno-align-commons -c -o $@ $<

