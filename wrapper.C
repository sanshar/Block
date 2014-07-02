#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "global.h"
#include "fciqmchelper.h"
#include "input.h"
#include <pario.h>
#include "spinblock.h"
#include "wrapper.h"

void ReadInput(char* conf);
namespace SpinAdapted{
MPS globalMPS;
}
using namespace SpinAdapted;

void ReadInputFromC(char* conf, int outputlevel) {
  ReadInput(conf);
  dmrginp.setOutputlevel() = outputlevel;
  dmrginp.initCumulTimer();
}

void readMPSFromDiskAndInitializeStaticVariables(int mpsindex) {
  if(!dmrginp.spinAdapted())
    MPS::sweepIters = dmrginp.last_site()/2-2;
  else
    MPS::sweepIters = dmrginp.last_site()-2;
  MPS::spinAdapted = false;
  for (int i=0; i<MPS::sweepIters+2; i++)
    MPS::siteBlocks.push_back(SpinBlock(i, i, false)); //alway make transpose operators as well

#ifndef SERIAL
  if (mpigetrank() == 0)
#endif
    globalMPS = ::MPS(mpsindex);

#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, globalMPS, 0);
#endif

}

void test()
{
  MPS statea(0);
  double o, h;
  calcHamiltonianAndOverlap(statea, statea, h, o);
  cout << o<<"  "<<h<<endl;
  MPS stateb(1);
  calcHamiltonianAndOverlap(stateb, statea, h, o);
  cout << o<<"  "<<h<<endl;
  stateb.normalize();
  calcHamiltonianAndOverlap(stateb, stateb, h, o);
  cout << o<<"  "<<h<<endl;
  calcHamiltonianAndOverlap(stateb, statea, h, o);
  cout << o<<"  "<<h<<endl;
}

void evaluateOverlapAndHamiltonian(long *occ, int length, double* o, double* h) {
  MPS dmrgc(occ, length);
  calcHamiltonianAndOverlap(globalMPS, dmrgc, *h, *o);
}


