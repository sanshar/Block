#specify boost include file
BOOSTINCLUDE = /home/sharma/boost_1_47_0/

#specify boost and lapack-blas library locations
BOOSTLIB = -L/home/sharma/boost_1_47_0/lib/ -lboost_serialization -lboost_system -lboost_filesystem
LAPACKBLAS = -L/srv/usr/local/opt/intel/mkl/10.2.1.017/lib/em64t/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#use these variable to set if we will use mpi or not 
USE_MPI = yes

#use this variable to set if we will intel compiler or not
INTEL = no


CXX = g++
HOME = .
NEWMATINCLUDE = $(HOME)/newmat10/
INCLUDE1 = $(HOME)/include/
INCLUDE2 = $(HOME)/
NEWMATLIB = $(HOME)/newmat10/
.SUFFIXES: .C .cpp

FLAGS =  -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOSTINCLUDE) -I$(HOME)/modules/twopdm/ -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm
LIBS =  -L$(NEWMATLIB) -lnewmat $(BOOSTLIB) $(LAPACKBLAS)  


ifeq ($(INTEL), yes)
	OPT = -O3 -funroll-loops -openmp  -DBLAS -DUSELAPACK  $(MPI_OPT) -DFAST_MTP 
#	OPT = -g -openmp  -DBLAS -DUSELAPACK  $(MPI_OPT) #-DFAST_MTP 
	CXX = icc
else
	OPT = -O3 -fopenmp   -DBLAS -DFAST_MTP -DUSELAPACK $(MPI_OPT)
#	OPT = -g -fopenmp   -DBLAS -DFAST_MTP -DUSELAPACK $(MPI_OPT)
endif


ifeq ($(USE_MPI), yes)
	MPI_OPT = 
	MPI_LIB = -L$(BOOST)/lib/ -lboost_mpi
	CXX = mpic++
else
	MPI_FLAG =
	MPI_LIB =
	MPI_OPT = -DSERIAL
endif

LIBS += $(MPI_LIB)




SRC_spin_adapted =  dmrg.C set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C anglib.C diis.C diis_updateError.C diis_updateHamiltonian.C diis_transformBlock.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C modules/onepdm/sweep.C modules/onepdm/onepdm.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C


SRC_spin_library =  dmrg.C readinput.C save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C spinblock.C StateInfo.C set_spinblock_components.C op_components.C Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C linear.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C anglib.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C  modules/onepdm/sweep.C modules/onepdm/onepdm.C  modules/generate_blocks/sweep.C diis.C diis_updateError.C diis_updateHamiltonian.C diis_transformBlock.C fci.C


OBJ_spin_adapted=$(SRC_spin_adapted:.C=.o)
OBJ_spin_library=$(SRC_spin_library:.C=.o)

.C.o :
	$(CXX) -fPIC $(FLAGS) $(OPT) -c $< -o $@
.cpp.o :
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@

all	: block.spin_adapted libqcdmrg.so

libqcdmrg.so : $(OBJ_spin_library)
	$(CXX) $(FLAGS) $(OPT) -shared -o libqcdmrg.so  $(OBJ_spin_library)

block.spin_adapted : $(OBJ_spin_adapted) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  block.spin_adapted $(OBJ_spin_adapted) $(LIBS) -lnewmat

$(NEWMATLIB)/libnewmat.a : 
	cd $(NEWMATLIB) && $(MAKE) -f makefile libnewmat.a

clean:
	rm *.o include/*.o modules/generate_blocks/*.o modules/onepdm/*.o modules/twopdm/*.o $(NEWMATLIB)*.o libqcdmrg.so block.spin_adapted

