#
#   makefile to create postproc.x
#
acc     = -r8
#deb     = -g
#nil     = -O2
vec     =

#FFLAGS  = -c $(deb) $(nil) $(vec) $(acc) -static -col120
#FFLAGS  = -c $(deb) $(nil) $(vec) $(acc) -static -extend_source
FFLAGS  = -c $(deb) $(nil) $(vec) $(acc) -col132 
OFLAGS  =    $(deb) $(nil) $(vec) $(acc)


OBJECTS     = \
   modchem.o \
   inputchem_simple.o \
   inputchem_mozart.o \
   bulk_chemistry.o  \
   iter_simple.o \
   iter_mozart.o
   

expanded_MXL : $(OBJECTS)
	ifort $(OFLAGS) -o $@ $(OBJECTS)

modchem.o  : modchem.f90
	ifort -c  modchem.f90

bulk_chemistry.o   : bulk_chemistry.f90 
	ifort -c  bulk_chemistry.f90

iter_simple.o : iter_simple.f90 
	ifort -c  iter_simple.f90
	
inputchem_simple.o : inputchem_simple.f90
	ifort -c  inputchem_simple.f90

iter_mozart.o : iter_mozart.f90 
	ifort -c  iter_mozart.f90
	
inputchem_mozart.o : inputchem_mozart.f90
	ifort -c  inputchem_mozart.f90
