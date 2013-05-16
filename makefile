# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define compiler and platform
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#fc      = ifort
fc      = gfortran
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Don't change below ~~~~~~~~~~~~~~~
ifeq ($(fc),gfortran)
  FCFLAGS   = -g -O3 -fdefault-real-8 -ffree-line-length-none
else
  FCFLAGS   = -g -O3 -r8
endif

OBJECTS     = \
   modchem.o \
   inputchem_simple.o \
   inputchem_mozart.o \
   bulk_chemistry.o  \
   iter_simple.o \
   iter_mozart.o
   

MXLCH_SOA : $(OBJECTS)
	$(fc) $(FCFLAGS) -o $@ $(OBJECTS)

modchem.o  : modchem.f90
	$(fc) $(FCFLAGS) -c  modchem.f90

bulk_chemistry.o   : bulk_chemistry.f90 
	$(fc) $(FCFLAGS) -c  bulk_chemistry.f90

iter_simple.o : iter_simple.f90 
	$(fc) $(FCFLAGS) -c  iter_simple.f90
	
inputchem_simple.o : inputchem_simple.f90
	$(fc) $(FCFLAGS) -c  inputchem_simple.f90

iter_mozart.o : iter_mozart.f90 
	$(fc) $(FCFLAGS) -c  iter_mozart.f90
	
inputchem_mozart.o : inputchem_mozart.f90
	$(fc) $(FCFLAGS) -c  inputchem_mozart.f90

clean:
	rm -f *.o *.mod MXLCH_SOA



