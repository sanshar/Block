#Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
#Copyright (c) 2012, Garnet K.-L. Chan                                                                                                                     
#This program is integrated in Molpro with the permission of 
#Sandeep Sharma and Garnet K.-L. Chan



#specify boost include file
BOOSTINCLUDE = /usr/include/boost

#specify boost and lapack-blas library locations
BOOSTLIB = -L/usr/lib/ -lboost_serialization -lboost_system -lboost_filesystem
LAPACKBLAS = -L/opt/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#use these variable to set if we will use mpi or not 
USE_MPI = yes

AR=ar
ARFLAGS=-qs
RANLIB=ranlib

# use this variable to set if we will use integer size of 8 or not.
# molpro compilation requires I8, since their integers are long
I8_OPT = yes
MOLPRO = no
OPENMP = no

ifeq ($(I8_OPT), yes)
	I8 = -DI8
endif

EXECUTABLE = block.spin_adapted

# change to icpc for Intel
CXX = g++
MPICXX = mpic++
HOME = .
NEWMATINCLUDE = $(HOME)/newmat10/
INCLUDE1 = $(HOME)/include/
INCLUDE2 = $(HOME)/
NEWMATLIB = $(HOME)/newmat10/
.SUFFIXES: .C .cpp

   

ifeq ($(MOLPRO), yes)
   MOLPROINCLUDE=-I $(HOME)/../
   MOLPRO_BLOCK= -DMOLPRO
endif
FLAGS =  -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOSTINCLUDE) -I$(HOME)/modules/twopdm/ -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm $(MOLPROINCLUDE)
LIBS =  -L$(NEWMATLIB) -lnewmat $(BOOSTLIB) $(LAPACKBLAS) -lgomp
MPI_OPT = -DSERIAL


# Intel compiler
ifeq ($(notdir $(firstword $(CXX))),icpc)
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -openmp -D_OPENMP 
   endif
	OPT = -DNDEBUG -O3 -funroll-loops $(OPENMP_FLAGS)
#	OPT = -g
#	CXX = icc
endif
# GNU compiler
ifeq ($(notdir $(firstword $(CXX))),g++)
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -fopenmp -D_OPENMP 
   endif
	OPT = -DNDEBUG -O3
#	OPT = -g
endif

OPT+=$(OPENMP_FLAGS) -DBLAS -DUSELAPACK $(MPI_OPT) $(I8) -DFAST_MTP $(MOLPRO_BLOCK)

ifeq ($(USE_MPI), yes)
	MPI_OPT = 
	MPI_LIB = -L$(BOOSTINCLUDE)/lib/ -lboost_mpi
   LIBS += $(MPI_LIB)
	CXX = $(MPICXX)
endif








SRC_genetic = genetic/CrossOver.C genetic/Evaluate.C genetic/GAInput.C genetic/GAOptimize.C genetic/Generation.C genetic/Mutation.C genetic/RandomGenerator.C genetic/ReadIntegral.C

SRC_spin_adapted =  dmrg.C fiedler.C least_squares.C set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C


SRC_spin_library =  fiedler.C IrrepSpace.C least_squares.C dmrg.C readinput.C save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C spinblock.C StateInfo.C set_spinblock_components.C op_components.C Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C linear.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C  modules/onepdm/sweep.C modules/onepdm/onepdm.C  modules/generate_blocks/sweep.C fci.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C


OBJ_spin_adapted=$(SRC_spin_adapted:.C=.o)
OBJ_spin_library=$(SRC_spin_library:.C=.o)

.C.o :
	$(CXX) -fPIC $(FLAGS) $(OPT) -c $< -o $@
.cpp.o :
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all	: $(EXECUTABLE) libqcdmrg.a

library : libqcdmrg.a $(NEWMATLIB)/libnewmat.a

libqcdmrg.a : $(OBJ_spin_library)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

$(EXECUTABLE) : $(OBJ_spin_adapted) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  $(EXECUTABLE) $(OBJ_spin_adapted) $(LIBS) -lnewmat

$(NEWMATLIB)/libnewmat.a : 
	cd $(NEWMATLIB) && $(MAKE) -f makefile libnewmat.a

clean:
	rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm/*.o $(NEWMATLIB)*.o libqcdmrg.so $(EXECUTABLE) $(NEWMATLIB)/libnewmat.a genetic/gaopt genetic/*.o

# DO NOT DELETE
