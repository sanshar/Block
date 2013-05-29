CXX=mpic++
CXXFLAGS=-g -O3 -fPIC

BOOSTDIR=/usr
BOOSTINCLUDE=$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib

NEWMATLIB=../newmat10

INCLUDES= -I. -I$(BOOSTINCLUDE) -I$(NEWMATLIB) -I.. -I../include
LIBS=$(BOOSTLIB) -lboost_serialization -lboost_mpi

FLAGS=$(CXXFLAGS)

OBJ_gaopt = main.o fiedler.o CrossOver.o Evaluate.o GAInput.o GAOptimize.o Generation.o Mutation.o RandomGenerator.o ReadIntegral.o

.C.o:
	$(CXX) $(FLAGS) $(INCLUDES) $(LIBS) -c $*.C

all     : gaopt

gaopt: $(OBJ_gaopt)
	$(CXX) $(FLAGS) -o $@ $(OBJ_gaopt) $(LIBS) -L$(NEWMATLIB) -lnewmat

clean:
	rm *.o; rm *.x; rm *.a
#
