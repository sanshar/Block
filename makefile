#Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
#Copyright (c) 2012, Garnet K.-L. Chan                                                                                                                     
#This program is integrated in Molpro with the permission of 
#Sandeep Sharma and Garnet K.-L. Chan



#specify boost include file
#BOOSTINCLUDE = /usr/include/boost
BOOSTINCLUDE = /home/mark/libs/boost_1_50_0/

EIGENINCLUDE = /home/mark/libs/eigen/


#specify boost and lapack-blas library locations
#BOOSTLIB = -L/usr/lib/ -lboost_serialization -lboost_system -lboost_filesystem
#LAPACKBLAS = -L/opt/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
BOOSTLIB = -L/home/mark/libs/boost_1_50_0/stage/lib/ -lboost_serialization -lboost_system -lboost_filesystem  #maw
LAPACKBLAS = -L/usr/lib -lblas -llapack_atlas -llapack

#use these variable to set if we will use mpi or not 
USE_MPI = yes
#USE_MPI = no

AR=ar
ARFLAGS=-qs
RANLIB=ranlib

# use this variable to set if we will use integer size of 8 or not.
# molpro compilation requires I8, since their integers are long
I8_OPT = yes
MOLPRO = no
OPENMP = no

##MAW FIXME THIS BREAKS boost/spirit stuff!!!
#MAWifeq ($(I8_OPT), yes)
#MAW	I8 = -DI8
#MAWendif

EXECUTABLE = block.spin_adapted

# change to icpc for Intel
#CXX = g++
MPICXX = mpic++ -std=c++0x
CXX=g++ -std=c++0x
OMPI_CXX=g++ -std=c++0x

HOME = .
NEWMATINCLUDE = $(HOME)/newmat10/
INCLUDE1 = $(HOME)/include/
INCLUDE2 = $(HOME)/
NEWMATLIB = $(HOME)/newmat10/
.SUFFIXES: .C .cpp

   

ifeq ($(MOLPRO), yes)
   MOLPROINCLUDE=$(HOME)/../
   MOLPRO_BLOCK= -DMOLPRO
endif
FLAGS =  -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOSTINCLUDE) -I$(EIGENINCLUDE) -I$(HOME)/modules/npdm/ -I$(HOME)/modules/two_index_ops/ -I$(HOME)/modules/three_index_ops/ -I$(HOME)/modules/four_index_ops/ -I$(HOME)/modules/twopdm_old/ -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm -I$(MOLPROINCLUDE)
LIBS =  -L$(NEWMATLIB) -lnewmat $(BOOSTLIB) $(LAPACKBLAS) -lgomp
MPI_OPT = -DSERIAL


ifeq ($(notdir $(firstword $(CXX))),icpc)
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -openmp -D_OPENMP 
   endif
# Intel compiler
	OPT = -O3 -funroll-loops $(OPENMP_INC) -DBLAS -DUSELAPACK  $(MPI_OPT) $(I8) -DFAST_MTP $(MOLPRO_BLOCK)
#	OPT = -ggdb3 -DBLAS -DUSELAPACK  $(MPI_OPT) -DFAST_MTP 
	CXX = icc
endif

ifeq ($(notdir $(firstword $(CXX))),g++)
   ifeq ($(OPENMP), yes)
      OPENMP_INC= -fopenmp -D_OPENMP 
   endif
# GNU compiler
	OPT = -O3 $(OPENMP_FLAGS) -DBLAS -DFAST_MTP -DUSELAPACK $(MPI_OPT) $(I8) $(MOLPRO_BLOCK)
#	OPT = -ggdb3 -fopenmp   -DBLAS -DFAST_MTP -DUSELAPACK $(MPI_OPT)
endif

ifeq ($(USE_MPI), yes)
	MPI_OPT = -g
	MPI_LIB = -L$(BOOSTINCLUDE)/lib/ -lboost_mpi
   LIBS += $(MPI_LIB)
	CXX = $(MPICXX)
endif








SRC_genetic = genetic/CrossOver.C genetic/Evaluate.C genetic/GAInput.C genetic/GAOptimize.C genetic/Generation.C genetic/Mutation.C genetic/RandomGenerator.C genetic/ReadIntegral.C

SRC_npdm = modules/npdm/npdm.cpp modules/npdm/npdm_driver.cpp modules/npdm/npdm_patterns.cpp modules/npdm/npdm_expectations.cpp modules/npdm/npdm_wrapper_selector.cpp modules/npdm/npdm_spin_adaptation.cpp modules/npdm/npdm_expectations_engine.cpp modules/npdm/npdm_spin_ops.cpp  modules/two_index_ops/two_index_op_components.cpp modules/two_index_ops/two_index_ops.cpp modules/two_index_ops/two_index_wrappers.cpp  modules/three_index_ops/three_index_op_components.cpp modules/three_index_ops/three_index_compound_ops.cpp modules/three_index_ops/three_index_ops.cpp modules/three_index_ops/three_index_wrappers.cpp modules/three_index_ops/build_3index_ops.cpp modules/four_index_ops/four_index_op_components.cpp modules/four_index_ops/four_index_compound_ops.cpp modules/four_index_ops/four_index_ops.cpp modules/four_index_ops/four_index_wrappers.cpp modules/four_index_ops/build_4index_ops.cpp modules/npdm/onepdm_container.cpp modules/npdm/twopdm_container.cpp modules/npdm/threepdm_container.cpp modules/npdm/fourpdm_container.cpp modules/npdm/npdm_permutations.cpp modules/npdm/nevpt2_npdm_driver.cpp modules/npdm/nevpt2_A16_matrix.cpp modules/npdm/nevpt2_npdm.cpp

SRC_spin_adapted =  dmrg.C least_squares.C set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm_old/sweep.C modules/twopdm_old/twopdm_old.C modules/twopdm_old/twopdm_2.C $(SRC_genetic) $(SRC_npdm)

SRC_spin_library =  IrrepSpace.C least_squares.C dmrg.C readinput.C save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C spinblock.C StateInfo.C set_spinblock_components.C op_components.C Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C linear.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C anglib.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/generate_blocks/sweep.C fci.C $(SRC_genetic) $(SRC_npdm)


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
	rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm_old/*.o modules/two_index_ops/*.o modules/three_index_ops/*.o modules/four_index_ops/*.o modules/npdm/*.o $(NEWMATLIB)*.o libqcdmrg.so $(EXECUTABLE) $(NEWMATLIB)/libnewmat.a genetic/gaopt genetic/*.o

# DO NOT DELETE
