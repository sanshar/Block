#CXX=/homec/naokin/gnu/gcc/4.8.2/bin/g++
CXX=mpic++
CXXFLAGS=-fPIC -O2

BOOSTROOT=/homec/naokin/boost/1.54.0
BLOCKROOT=/homec/naokin/devel/Block
INTELROOT=/home100/opt/intel

INCLUDEFLAGS=-I. -I$(BLOCKROOT) -I$(BLOCKROOT)/include -I$(BOOSTROOT)/include -I$(INTELROOT)/mkl/include

LIBS=-L$(BLOCKROOT) -lqcdmrg -L$(BLOCKROOT)/newmat10 -lnewmat -L$(BOOSTROOT)/lib -lboost_serialization -lboost_system -lboost_filesystem -lboost_mpi -L$(INTELROOT)/lib/intel64 -L$(INTELROOT)/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread

SRC_LIBDMRG = block_calldmrg.C molpro_fcidump.C loadNpdm.C sortNpdm.C tranNpdm.C
OBJ_LIBDMRG = $(SRC_LIBDMRG:.C=.o)

.C.o:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: libdmrg.so

libdmrg.so : $(OBJ_LIBDMRG)
	$(CXX) -shared -Wl,-soname=libqcdmrg.so.1 -o libqcdmrg.so $(OBJ_LIBDMRG) $(LIBS)
	mv libqcdmrg.so ..

clean:
	rm *.so; rm *.o; rm *.x; rm *.a;
