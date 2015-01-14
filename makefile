#Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
#Copyright (c) 2012, Garnet K.-L. Chan                                                                                                                     
#This program is integrated in Molpro with the permission of 
#Sandeep Sharma and Garnet K.-L. Chan

##BOOSTINCLUDE = /home/sandeep/Work/Programs/boost_1_54_0/
#specify boost include file
BOOSTINCLUDE = /home/juny/libs/boost_1_55_0-gcc-4_8_2-install/include

#specify boost and lapack-blas library locations
BOOSTLIB = -L/home/juny/libs/boost_1_55_0-gcc-4_8_2-install/lib/ -lboost_serialization -lboost_system -lboost_filesystem
#BOOSTLIB = -lboost_serialization -lboost_system -lboost_filesystem
LAPACKBLAS = -lblas -llapack

USE_BOOST56 = no
ifeq ($(USE_BOOST56), yes)
	B56 = -DBOOST_1_56_0
endif

#use these variable to set if we will use mpi or not 
USE_MPI = yes
USE_MKL = yes

# use this variable to set if we will use integer size of 8 or not.
# molpro compilation requires I8, since their integers are long
I8_OPT = no
MOLPRO = no
OPENMP = no

# add Molcas interface to libqcdmrg.so
# molcas compilation w/ -64 option requires I8 as well
MOLCAS = no

ifeq ($(USE_MKL), yes)
MKLLIB = /opt/intel/composer_xe_2013_sp1.0.080/mkl/lib/intel64/
LAPACKBLAS = -L${MKLLIB} -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
MKLFLAGS = /opt/intel/composer_xe_2013_sp1.0.080/mkl/include
MKLOPT = -D_HAS_INTEL_MKL
else
MKLFLAGS = .
endif

RUN_UNITTEST=no
ifeq ($(RUN_UNITTEST), yes)
UNITTEST = -DUNITTEST
endif

AR=ar
ARFLAGS=-qs
RANLIB=ranlib

ifeq ($(I8_OPT), yes)
	I8 = -D_FORTINT_64
endif

EXECUTABLE = block.spin_adapted

# change to icpc for Intel
CXX =  g++
MPICXX = /usr/local/openmpi/1.6.5/devtoolset21/x86_64/bin/mpicxx
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
   MOLPROINCLUDE=$(BLOCKHOME)/../
   MOLPRO_BLOCK= -DMOLPRO
endif

FLAGS =  -I${MKLFLAGS} -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOSTINCLUDE) -I$(MOLPROINCLUDE) \
         -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm -I$(HOME)/modules/twopdm/ \
         -I$(HOME)/modules/npdm -I$(HOME)/modules/two_index_ops -I$(HOME)/modules/three_index_ops -I$(HOME)/modules/four_index_ops -std=c++0x \
	 -I$(HOME)/modules/ResponseTheory -I$(HOME)/modules/nevpt2 -I$(HOME)/molcas

LIBS +=  -L$(NEWMATLIB) -lnewmat $(BOOSTLIB) $(LAPACKBLAS) -lgomp 
MPI_OPT = -DSERIAL




ifeq (icpc, $(CXX))
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -openmp -D_OPENMP 
   endif
# Intel compiler
   OPT = -DNDEBUG -O3 -funroll-loops  
#  OPT = -g -fPIC
   ifeq ($(USE_MPI), no) 
      CXX = icc
   endif
endif

ifeq (g++, $(CXX))
   ifeq ($(OPENMP), yes)
      OPENMP_FLAGS= -fopenmp -D_OPENMP 
   endif
# GNU compiler
      OPT = -DNDEBUG -O3 
      GLIBS = -L/home/juny/libs/gcc-4.8.3-install/lib64 
      LIBS += $(GLIBS)
#   OPT = -g -fPIC
endif

ifeq ($(USE_MPI), yes)
     MPI_OPT = 
     MPI_LIB = -lboost_mpi -L/usr/local/openmpi/1.6.5/devtoolset21/x86_64/lib64 -lmpi -lmpi_cxx
     LIBS += $(MPI_LIB)
     CXX = $(MPICXX)
endif


OPT	+= $(OPENMP_FLAGS) -DBLAS -DUSELAPACK $(MPI_OPT) $(I8) $(B56) $(MOLPRO_BLOCK)  -DFAST_MTP -D_HAS_CBLAS -D_HAS_INTEL_MKL ${MKLOPT} ${UNITTEST} -fPIC

SRC_genetic = genetic/CrossOver.C genetic/Evaluate.C genetic/GAInput.C genetic/GAOptimize.C genetic/Generation.C genetic/Mutation.C genetic/RandomGenerator.C genetic/ReadIntegral.C

SRC_npdm = modules/npdm/npdm.C modules/npdm/npdm_driver.C modules/npdm/npdm_patterns.C modules/npdm/npdm_expectations.C modules/npdm/npdm_expectations_engine.C  \
           modules/npdm/npdm_permutations.C modules/npdm/npdm_spin_adaptation.C modules/npdm/npdm_operator_selector.C modules/npdm/npdm_spin_ops.C \
           modules/npdm/npdm_array_buffer.C modules/npdm/onepdm_container.C modules/npdm/twopdm_container.C modules/npdm/threepdm_container.C modules/npdm/fourpdm_container.C  \
           modules/two_index_ops/two_index_wrappers.C modules/three_index_ops/three_index_wrappers.C modules/four_index_ops/four_index_wrappers.C  \
           modules/three_index_ops/three_index_compound_ops.C modules/four_index_ops/four_index_compound_ops.C  \
           modules/three_index_ops/three_index_op_components.C modules/four_index_ops/four_index_op_components.C  \
           modules/three_index_ops/three_index_ops.C modules/four_index_ops/four_index_ops.C  \
           modules/three_index_ops/build_3index_ops.C modules/four_index_ops/build_4index_ops.C modules/npdm/pairpdm_container.C \
           modules/npdm/nevpt2_npdm_driver.C modules/npdm/nevpt2_A16_container.C modules/npdm/nevpt2_npdm_matrices.C modules/npdm/externalsort.C

SRC_spin_adapted =  modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C hfOccGenerator.C Schedule.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2)

SRC_OH = modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2)
OBJ_OH= OverlapHelement.o

SRC_CSFOH = CSFOverlapHelement.C modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2)


SRC_COEF= modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C wrapper.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C Schedule.C input.C hfOccGenerator.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2)
OBJ_COEF= Coefficients.o

SRC_spin_library =  modules/ResponseTheory/sweepResponse.C modules/ResponseTheory/sweepCompress.C fciqmchelper.C pario.C dmrg.C fiedler.C least_squares.C sweep_mps.C set_spinblock_components.C linear.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C hfOccGenerator.C Schedule.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C new_anglib.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C $(SRC_genetic) SpinSpace.C include/IntegralMatrix.C $(SRC_npdm) $(SRC_nevpt2) 

SRC_nevpt2 = modules/nevpt2/nevpt2.C modules/nevpt2/nevpt2_info.C modules/nevpt2/nevpt2_mpi.C \
             modules/nevpt2/nevpt2_opconstruct.C modules/nevpt2/nevpt2_operators.C \
             modules/nevpt2/nevpt2_pal.C modules/nevpt2/nevpt2_renormalize.C \
             modules/nevpt2/nevpt2_util.C modules/nevpt2/ripdm.C modules/nevpt2/sweep_gen_nevpt2.C \
             modules/nevpt2/sweep_nevpt2.C

SRC_molcas = molcas/block_calldmrg.C molcas/molpro_fcidump.C molcas/loadNpdm.C molcas/sortNpdm.C molcas/tranNpdm.C

OBJ_OH+=$(SRC_OH:.C=.o)
OBJ_CSFOH+=$(SRC_CSFOH:.C=.o)
OBJ_COEF+=$(SRC_COEF:.C=.o)
OBJ_spin_adapted=$(SRC_spin_adapted:.C=.o)
OBJ_spin_library=$(SRC_spin_library:.C=.o)
OBJ_nevpt2=$(SRC_nevpt2:.C=.o)
ifeq ($(MOLCAS), yes)
	OBJ_spin_library+=$(SRC_molcas:.C=.o)
endif

.C.o :
	$(CXX)  $(FLAGS) $(OPT) -c $< -o $@
.cpp.o :
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all	: library $(EXECUTABLE) OH COEF CSFOH

library : libqcdmrg.a $(NEWMATLIB)/libnewmat.a libqcdmrg.so

libqcdmrg.a : $(OBJ_spin_library) $(OBJ_molcas)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

libqcdmrg.so : $(OBJ_spin_library) $(OBJ_molcas)
	$(CXX) -shared -o $@ $^ $(LIBS)

$(EXECUTABLE) : $(OBJ_spin_adapted) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  $(EXECUTABLE) $(OBJ_spin_adapted) $(LIBS)

OH : $(OBJ_OH) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  OH $(OBJ_OH) $(LIBS)

CSFOH : $(OBJ_CSFOH) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  CSFOH $(OBJ_CSFOH) $(LIBS)

COEF : $(OBJ_COEF) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  COEF $(OBJ_COEF) $(LIBS)

$(NEWMATLIB)/libnewmat.a :
	cd $(NEWMATLIB) && $(MAKE) -f makefile libnewmat.a

clean:
	rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm/*.o modules/npdm/*.o $(NEWMATLIB)*.o libqcdmrg.a libqcdmrg.so $(EXECUTABLE) $(NEWMATLIB)/libnewmat.a genetic/gaopt genetic/*.o btas/lib/*.o btas/lib/libbtas.a modules/two_index_ops/*.o modules/three_index_ops/*.o modules/four_index_ops/*.o modules/ResponseTheory/*.o modules/nevpt2/*.o molcas/*.o

# DO NOT DELETE


