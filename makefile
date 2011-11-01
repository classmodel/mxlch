#
#   makefile to create postproc.x
#
#acc     = -fdefault-real-8
#deb     = -g
#nil     = -O2
vec     =

#FFLAGS  = -c $(deb) $(nil) $(vec) $(acc) -static -col120
#FFLAGS  = -c $(deb) $(nil) $(vec) $(acc) -static -extend_source
FFLAGS  = -c $(deb) $(nil) $(vec) $(acc) -ffree-line-length-none
OFLAGS  =    $(deb) $(nil) $(vec) $(acc) 


OBJECTS     = \
   modchem.o \
   inputchem_simple.o \
   inputchem_mozart.o \
   bulk_chemistry.o  \
   iter_simple.o \
   iter_mozart.o
   

expanded_MXL.exe : $(OBJECTS)
	gfortran $(OFLAGS) -o $@ $(OBJECTS)

modchem.o  : modchem.f90
	gfortran $(FFLAGS)  modchem.f90

bulk_chemistry.o   : bulk_chemistry.f90 
	gfortran $(FFLAGS)  bulk_chemistry.f90

iter_simple.o : iter_simple.f90 
	gfortran $(FFLAGS)  iter_simple.f90
	
inputchem_simple.o : inputchem_simple.f90
	gfortran $(FFLAGS)  inputchem_simple.f90

iter_mozart.o : iter_mozart.f90 
	gfortran $(FFLAGS)  iter_mozart.f90
	
inputchem_mozart.o : inputchem_mozart.f90
	gfortran $(FFLAGS)  inputchem_mozart.f90
