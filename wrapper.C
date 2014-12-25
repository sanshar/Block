#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "global.h"
#include "fciqmchelper.h"
#include "input.h"
#include "spinblock.h"
#include "wrapper.h"
#include "rotationmat.h"
#include <sstream>


void ReadInput(char* conf);
namespace SpinAdapted{
MPS globalMPS;
}
using namespace SpinAdapted;

void initBoostMPI(int argc, char* argv[]) {
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
#endif
}

void ReadInputFromC(char* conf, int outputlevel) {
  ReadInput(conf);
  dmrginp.setOutputlevel() = outputlevel;
  dmrginp.initCumulTimer();
}

void initializeGlobalMPS(int mpsindex) {
  SpinAdapted::globalMPS = MPS(mpsindex);
}

void readMPSFromDiskAndInitializeStaticVariables() {
  if(!dmrginp.spinAdapted())
    MPS::sweepIters = dmrginp.last_site()/2-2;
  else
    MPS::sweepIters = dmrginp.last_site()-2;
  MPS::spinAdapted = false;
  for (int i=0; i<MPS::sweepIters+2; i++) {
    if (i==0 && dmrginp.spinAdapted() && dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 ) {
      SpinBlock s(i, i, 0, false);
      SpinQuantum sq = dmrginp.molecule_quantum();
      sq = SpinQuantum(sq.get_s().getirrep(), sq.get_s(), IrrepSpace(0));
      int qs = 1, ns = 1;
      StateInfo addstate(ns, &sq, &qs); 
      SpinBlock dummyblock(addstate, 0);
      SpinBlock newstartingBlock;
      newstartingBlock.set_integralIndex() = 0;
      newstartingBlock.default_op_components(false, s, dummyblock, true, true, false);
      newstartingBlock.setstoragetype(LOCAL_STORAGE);
      newstartingBlock.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, s, dummyblock);
      MPS::siteBlocks.push_back(newstartingBlock); //alway make transpose operators as well
    }
    else
      MPS::siteBlocks.push_back(SpinBlock(i, i, 0, false)); //alway make transpose operators as well
  }
}

void writeFullMPS()
{
  int ncore = 2;
  int nactive = 8;
  int nvirt = 18;

  Matrix m(1,1); m(1,1)=1.0;
  std::vector<Matrix> corerotMat, virtcorMat, firstactMat;
  corerotMat.resize(3); virtcorMat.resize(3);firstactMat.resize(3);
  corerotMat[2] = m; virtcorMat[0] = m;
  firstactMat[0] = m; firstactMat[1] = m; firstactMat[2] = m;

  std::vector<Matrix> rotMat_init;rotMat_init.resize(6);
  rotMat_init[5] = m;

  std::vector<int> readSites(2,0); 
  std::vector< std::vector<Matrix> > activeRots; activeRots.resize(nactive);
  for (int i=ncore+1; i<ncore+nactive-1; i++) {
    readSites[1]++;
    LoadRotationMatrix(readSites, activeRots[i-ncore-1], 0); 
  }


  std::vector<int> sites(2,0); sites[1] = 1;
  SaveRotationMatrix(sites, rotMat_init, 0);

  //write ncore rotation matrices
  for (int i=2; i<ncore; i++) {
    sites[1] = i;
    SaveRotationMatrix(sites, rotMat_init, 0);
  }


  sites[1] = ncore;
  SaveRotationMatrix(sites, firstactMat, 0);


  //write ncore rotation matrices
  for (int i=ncore+1; i<ncore+nactive-1; i++) {
    sites[1] = i;
    SaveRotationMatrix(sites, activeRots[i-ncore-1], 0);
  }

  sites[1] = ncore+nactive-1;
  SaveRotationMatrix(sites, firstactMat, 0);

  //write ncore rotation matrices
  for (int i=ncore+nactive; i<ncore+nactive+nvirt; i++) {
    sites[1] = i;
    SaveRotationMatrix(sites, virtcorMat, 0);
  }


}

void test()
{
  setbuf(stdout, NULL);
  pout.precision(12);
  MPS statea(0);
  double o, h;
  calcHamiltonianAndOverlap(statea, statea, h, o);

  if (mpigetrank() == 0)
    printf("<0|0> = %18.10f   <0|H|0> = %18.10f\n", o, h);
  double e0 = h;

  for (int i=1; i<60; i++) {
    MPS statec(i);
    double h1,o1;
    calcHamiltonianAndOverlap(statea, statec, h, o);
    calcHamiltonianAndOverlap(statec, statec, h1, o1);
    if (mpigetrank() == 0)
      printf("<0|%i> = %18.10e   <0|H|%i> = %18.10e  <%i|%i> = %18.10e  <%i|H|%i> = %18.10e  perturb=%18.10e\n", i, o, i, h, i, i, o1, i, i, h1, h1-e0*o1);
  }
  /*
  MPS statec(2);
  calcHamiltonianAndOverlap(statea, statec, h, o);
  if (mpigetrank() == 0)
    printf("<0|1> = %18.10f   <0|H|1> = %18.10f\n", o, h);


  //MPS stateb(1);
  //calcHamiltonianAndOverlap(stateb, stateb, h, o);
  //pout << o<<"  "<<h<<endl;

  calcHamiltonianAndOverlap(statec, statec, h, o);
  if (mpigetrank() == 0)
    printf("<1|1> = %18.10f   <1|H|1> = %18.10f\n", o, h);
  */
}

void evaluateOverlapAndHamiltonian(unsigned long *occ, int length, double* o, double* h) {
  MPS dmrgc(occ, length);
  calcHamiltonianAndOverlap(SpinAdapted::globalMPS, dmrgc, *h, *o);
}


void intFromString(unsigned long &occ, const char* s) {
  occ = 0;
  long temp = 1;
  string ss(s);
  stringstream stream(ss);
  int n, i=0;
  while (stream >>n) {
    if (n==1)
      occ = occ | temp <<(63-i);
    i++;
  }
  return;
}
