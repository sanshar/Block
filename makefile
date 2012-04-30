CXX = mpic++
F90 = gfortran
HOME = .
NEWMATINCLUDE = $(HOME)/newmat10/
INCLUDE1 = $(HOME)/include/
INCLUDE2 = $(HOME)/
BOOST = /home/sharma/boost_1_47_0/
NEWMATLIB = $(HOME)/newmat10/
.SUFFIXES: .C .cpp

FLAGS =  -I$(INCLUDE1) -I$(INCLUDE2) -I$(NEWMATINCLUDE) -I$(BOOST) -I$(HOME)/modules/twopdm/ -I$(HOME)/modules/generate_blocks/ -I$(HOME)/modules/onepdm -I$(GSL)/ 
LIBS =  -L$(NEWMATLIB) -lnewmat -L$(BOOST)/lib/ -lboost_serialization -L/srv/usr/local/opt/intel/mkl/10.2.1.017/lib/em64t/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core


#use this variable to set if we will use mpi or not
USE_MPI = yes
INTEL = no

ifeq ($(USE_MPI), yes)
	MPI_FLAG = $(BOOST)lib/libboost_mpi.a
	MPI_LIB = -L$(BOOST)/lib/ -lboost_mpi
	MPI_OPT = 
else
	MPI_FLAG =
	MPI_LIB =
	MPI_OPT = -DSERIAL
endif

ifeq ($(INTEL), yes)
	OPT = -O3 -funroll-loops -openmp  -DBLAS -DUSELAPACK  $(MPI_OPT) -DFAST_MTP 
#	OPT = -g -openmp  -DBLAS -DUSELAPACK  $(MPI_OPT) #-DFAST_MTP 
	CXX = mpic++
	LIBS += $(MPI_LIB)
else
	OPT = -O3 -fopenmp   -DBLAS -DFAST_MTP -DUSELAPACK $(MPI_OPT)
#	OPT = -g -fopenmp   -DBLAS -DFAST_MTP -DUSELAPACK $(MPI_OPT)
	LIBS += $(MPI_LIB)
endif



SRC_spin_adapted =  set_spinblock_components.C linear.C main.C readinput.C  save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C StateInfo.C  Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C anglib.C diis.C diis_updateError.C diis_updateHamiltonian.C diis_transformBlock.C fci.C spinblock.C op_components.C IrrepSpace.C modules/generate_blocks/sweep.C 
#SRC_spin_adapted =  main.C readinput.C save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C StateInfo.C set_spinblock_components.C op_components.C Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C linear.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C anglib.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C  modules/onepdm/sweep.C modules/onepdm/onepdm.C  modules/generate_blocks/sweep.C diis.C diis_updateError.C diis_updateHamiltonian.C diis_transformBlock.C fci.C spinblock.C

SRC_spin_library =  dmrg.C readinput.C save_load_block.C timer.C SpinQuantum.C Symmetry.C input.C orbstring.C slater.C csf.C spinblock.C StateInfo.C set_spinblock_components.C op_components.C Operators.C BaseOperator.C screen.C MatrixBLAS.C operatorfunctions.C opxop.C wavefunction.C solver.C linear.C davidson.C sweep_params.C sweep.C initblocks.C guess_wavefunction.C density.C rotationmat.C renormalise.C couplingCoeffs.C distribute.C anglib.C modules/twopdm/sweep.C modules/twopdm/twopdm.C modules/twopdm/twopdm_2.C  modules/onepdm/sweep.C modules/onepdm/onepdm.C  modules/generate_blocks/sweep.C diis.C diis_updateError.C diis_updateHamiltonian.C diis_transformBlock.C fci.C


OBJ_spin_adapted=$(SRC_spin_adapted:.C=.o)
OBJ_spin_library=$(SRC_spin_library:.C=.o)

.C.o :
	$(CXX) -fPIC $(FLAGS) $(OPT) -c $< -o $@
.cpp.o :
	$(CXX) $(FLAGS) $(OPT) -c $< -o $@
.f90.o :
	$(F90) -c $< -o $@ 

all	: block.spin_adapted libqcdmrg.so

libqcdmrg.so : $(OBJ_spin_library)
	$(CXX) $(FLAGS) $(OPT) -shared -o libqcdmrg.so  $(OBJ_spin_library)
#	mv libqcdmrg.so /home/sharma/.

block.spin_adapted : $(OBJ_spin_adapted) $(NEWMATLIB)/libnewmat.a
	$(CXX)   $(FLAGS) $(OPT) -o  block.spin_adapted $(OBJ_spin_adapted) $(LIBS) -lnewmat

$(NEWMATLIB)/libnewmat.a : 
	cd $(NEWMATLIB) && $(MAKE) -f makefile libnewmat.a

clean:
	rm *.o include/*.o modules/generate_blocks/*.o modules/twopdm/*.o modules/twopdm/*.o $(NEWMATLIB)*.o libqcdmrg.so block.spin_adapted

