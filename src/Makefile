#*****************************************************
#  General Makefile
#
#*****************************************************

#**************************************************************************
# Macros defining the Fortran, C/C++ compiler and linker.

CC = mpif90
LINK = mpif90
FFLAGS = -O3 -fallow-argument-mismatch
#FFLAGS = -g -fallow-argument-mismatch

#**************************************************************************
# List of .o files that EXENAME depends on.  Edit as appropriate for MP.

OBJS = \
	src/DataStruct/NumConst.o src/DataStruct/PhysConst.o src/DataStruct/Pgrid.o \
	src/DataStruct/Data.o \
	src/Func/Timer.o src/Func/Transpose.o src/Func/Fldmger.o src/Func/Ptclmger.o  \
	src/Func/FFT.o  src/Func/math.o \
        src/Appl/BPM.o src/Appl/CCL.o src/Appl/CCDTL.o src/Appl/DTL.o src/Appl/SC.o src/Appl/Multipole.o \
	src/Appl/DriftTube.o src/Appl/Quadrupole.o src/Appl/ConstFoc.o src/Appl/SolRF.o \
	src/Appl/Sol.o src/Appl/Dipole.o src/Appl/EMfld.o src/Appl/EMfldCart.o src/Appl/EMfldCyl.o src/Appl/EMfldAna.o \
	src/Appl/BeamLineElem.o \
	src/Appl/CompDom.o src/Appl/BeamBunch.o src/Appl/Depositor.o src/Appl/Ranger.o src/Appl/Field.o src/Appl/Distribution.o \
	src/Contrl/Input.o src/Contrl/Output.o src/Contrl/AccSimulator.o src/Contrl/main.o

OBJS2 = \
	NumConst.o PhysConst.o Pgrid.o Data.o \
        Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o math.o \
	BPM.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o Multipole.o ConstFoc.o \
	SolRF.o Sol.o Dipole.o EMfld.o EMfldCart.o EMfldCyl.o EMfldAna.o \
	BeamLineElem.o CompDom.o BeamBunch.o Depositor.o Ranger.o Field.o Distribution.o \
	Input.o Output.o AccSimulator.o main.o	
#**************************************************************************
# Change this line if you don't like 'a.out'.

EXENAME = ImpactT.exe

#************************************************************************
# disable predefined suffixes and define your own set of allowed suffixes
 .SUFFIXES:
 .SUFFIXES: .o .f .F .c .f90 .F90

#*************************************************************************
# inference rules (how to compile object files that have no explicit rules)
#  $* = name part of target
#  $@ = full target name
#  $< = dependent name

.f90.o:
	$(CC) -c $(FFLAGS) $<

#**************************************************************************
# Rules for building EXENAME from OBJS and OBJS from your source.

$(EXENAME):  $(OBJS) 
	$(LINK) -o $(EXENAME) $(OBJS2) 

#************************************************************************
# if you wish to compile a certain object with different flags
# or in some special way, then specify the target & dependency explicitly
# the following line say Timer.o is depended on Timer.f90
#Timer.o: Timer.f90
#	$(CC) -c -O3 Timer.f90

	cp  AccSimulator.o main.o Input.o Output.o src/Contrl
	cp  BPM.o CCL.o CCDTL.o DTL.o SC.o DriftTube.o Quadrupole.o \
	    ConstFoc.o BeamLineElem.o BeamBunch.o Field.o CompDom.o \
	    Distribution.o SolRF.o Sol.o Dipole.o Multipole.o EMfld.o \
	    EMfldCart.o EMfldCyl.o EMfldAna.o Depositor.o Ranger.o src/Appl
	cp  Timer.o Transpose.o Fldmger.o Ptclmger.o FFT.o math.o src/Func	
	cp  NumConst.o PhysConst.o Data.o Pgrid.o src/DataStruct 
#***********************************************************************
clean:
	-rm -f *.o *.mod src/*/*.o 
