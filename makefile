#Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
#Copyright (c) 2012, Garnet K.-L. Chan                                                                                                                     
#This program is integrated in Molpro with the permission of 
#Sandeep Sharma and Garnet K.-L. Chan

##BOOSTINCLUDE = /home/sandeep/Work/Programs/boost_1_54_0/
#specify boost include file
BOOSTINCLUDE = /home/juny/boost_1_55_0/install/include

#specify boost and lapack-blas library locations
BOOSTLIB =  -L/home/juny/boost_1_55_0/install/lib -lboost_serialization -lboost_system -lboost_filesystem
LAPACKBLAS = -L/opt/intel/composer_xe_2013_sp1.0.080/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#use these variable to set if we will use mpi or not 
USE_MPI = yes
USE_MKL = no

ifeq ($(USE_MKL), yes)
MKLLIB = /opt/intel/composer_xe_2013_sp1.0.080/mkl/lib/intel64/
LAPACKBLAS = -L${MKLLIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
MKLFLAGS = /opt/intel/composer_xe_2013_sp1.0.080/mkl/include
MKLOPT = -D_HAS_INTEL_MKL
else
MKLFLAGS=.
endif

RUN_UNITTEST=no
ifeq ($(RUN_UNITTEST), yes)
UNITTEST = -DUNITTEST
endif

AR=ar
ARFLAGS=-qs
RANLIB=ranlib

# use this variable to set if we will use integer size of 8 or not.
# molpro compilation requires I8, since their integers are long
I8_OPT = no
MOLPRO = no
OPENMP = no

ifeq ($(I8_OPT), yes)
	I8 = -DMolpro_I8
endif

EXECUTABLE = block.spin_adapted

# change to icpc for Intel
CXX =  g++
MPICXX = mpicxx
BLOCKHOME = .
HOME = .
NEWMATINCLUDE = $(BLOCKHOME)/newmat10/
INCLUDE1 = $(BLOCKHOME)/include/
INCLUDE2 = $(BLOCKHOME)/
NEWMATLIB = $(BLOCKHOME)/newmat10/
BTAS = $(BLOCKHOME)/btas
.SUFFIXES: .C .cpp

   
MOLPROINCLUDE=.
ifeq ($(MOLPRO), yes)
   MOLPROINCLUDE=$(HOME)/../
   MOLPRO_BLOCK= -DMOLPRO
endif

FLAGS =  -I${MKLFLAGS} -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOSTINCLUDE) -I$(MOLPROINCLUDE) \
         -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm -I$(HOME)/modules/twopdm/ \
         -I$(HOME)/modules/npdm -I$(HOME)/modules/two_index_ops -I$(HOME)/modules/three_index_ops -I$(HOME)/modules/four_index_ops -std=c++0x \
	 -I$(HOME)/modules/ResponseTheory 

LIBS +=  -L$(NEWMATLIB) -lnewmat $(BOOSTLIB) $(LAPACKBLAS) -lgomp 
MPI_OPT = -DSERIAL

ifeq ($(notdir $(firstword $(CXX))),icpc)
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -openmp -D_OPENMP 
   endif
# Intel compiler
	OPT = -DNDEBUG -O3 -funroll-loops 
#	OPT = -g 
	CXX = icc
endif

ifeq ($(notdir $(firstword $(CXX))),g++)
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -fopenmp -D_OPENMP 
   endif
# GNU compiler
	OPT = -DNDEBUG -O3 -fPIC
#	OPT = -g
endif

ifeq ($(USE_MPI), yes)
	MPI_OPT = 
	MPI_LIB = -lboost_mpi
        LIBS += $(MPI_LIB)
	CXX = $(MPICXX)
endif


OPT	+= $(OPENMP_FLAGS) -DBLAS -DUSELAPACK $(MPI_OPT) $(I8) $(MOLPRO_BLOCK)  -DFAST_MTP -D_HAS_CBLAS -D_HAS_INTEL_MKL ${MKLOPT} ${UNITTEST}

SRC_genetic = genetic/CrossOver.C genetic/Evaluate.C genetic/GAInput.C genetic/GAOptimize.C genetic/Generation.C genetic/Mutation.C genetic/RandomGenerator.C genetic/ReadIntegral.C

SRC_npdm = modules/npdm/npdm.C modules/npdm/npdm_driver.C modules/npdm/npdm_patterns.C modules/npdm/npdm_expectations.C modules/npdm/npdm_expectations_engine.C  \
           modules/npdm/npdm_permutations.C modules/npdm/npdm_spin_adaptation.C modules/npdm/npdm_operator_selector.C modules/npdm/npdm_spin_ops.C \
           modules/npdm/npdm_array_buffer.C modules/npdm/onepdm_container.C modules/npdm/twopdm_container.C modules/npdm/threepdm_container.C modules/npdm/fourpdm_container.C  \
           modules/two_index_ops/two_index_wrappers.C modules/three_index_ops/three_index_wrappers.C modules/four_index_ops/four_index_wrappers.C  \
           modules/three_index_ops/three_index_compound_ops.C modules/four_index_ops/four_index_compound_ops.C  \
           modules/three_index_ops/three_index_op_components.C modules/four_index_ops/four_index_op_components.C  \
           modules/three_index_ops/three_index_ops.C modules/four_index_ops/four_index_ops.C  \
           modules/three_index_ops/build_3index_ops.C modules/four_index_ops/build_4index_ops.C \
           modules/npdm/nevpt2_npdm_driver.C modules/npdm/nevpt2_A16_container.C modules/npdm/nevpt2_npdm_matrices.C modules/npdm/externalsort.C

SRC_spin_adapted =  modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C Schedule.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) 

SRC_OH = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) 
OBJ_OH= OverlapHelement.o

SRC_spin_library =  modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C fciqmchelper.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C Schedule.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) 

OBJ_OH+=$(SRC_OH:.C=.o)
OBJ_spin_adapted=$(SRC_spin_adapted:.C=.o)
OBJ_spin_library=$(SRC_spin_library:.C=.o)

.C.o :
	$(CXX)  $(FLAGS) $(OPT) -c $< -o $@
.cpp.o :
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all	: $(EXECUTABLE) libqcdmrg.a OH

library : libqcdmrg.a $(NEWMATLIB)/libnewmat.a libqcdmrg.so

libqcdmrg.a : $(OBJ_spin_library)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

libqcdmrg.so : $(OBJ_spin_library)
	$(CXX) -shared -o $@ $^ $(LIBS)

$(EXECUTABLE) : $(OBJ_spin_adapted) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  $(EXECUTABLE) $(OBJ_spin_adapted) $(LIBS)

OH : $(OBJ_OH) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  OH $(OBJ_OH) $(LIBS)

$(NEWMATLIB)/libnewmat.a : 
	cd $(NEWMATLIB) && $(MAKE) -f makefile libnewmat.a

clean:
	rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm/*.o modules/npdm/*.o $(NEWMATLIB)*.o libqcdmrg.so $(EXECUTABLE) $(NEWMATLIB)/libnewmat.a genetic/gaopt genetic/*.o btas/lib/*.o btas/lib/libbtas.a modules/two_index_ops/*.o modules/three_index_ops/*.o modules/four_index_ops/*.o modules/ResponseTheory/*.o

# DO NOT DELETE


