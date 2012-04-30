#include "IntegralMatrix.h"
#include <fstream>
#include "input.h"
#include "pario.h"
#include "global.h"
#include "orbstring.h"
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <include/communicate.h>
#include "omp.h"
//the following can be removed later
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "spinblock.h"
#include "StateInfo.h"
#include "operatorfunctions.h"
#include "wavefunction.h"
#include "solver.h"
#include "davidson.h"
#include "guess_wavefunction.h"
#include "rotationmat.h"
#include "density.h"
#include "sweep.h"
#include "BaseOperator.h"

void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
//void test(double sweep_tol);


Timer globaltimer(false);
bool DEBUGWAIT = false;
bool DEBUG_MEMORY = false;
bool restartwarm = false;
OneElectronArray v_1;
TwoElectronArray v_2;
Input dmrginp;
int MAX_THRD = 1;


int main(int argc, char* argv [])
{
#ifndef SERIAL
  mpi::environment env(argc, argv);
  MPI_Comm mpiparent;
  MPI_Comm_get_parent(&mpiparent); 

  mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  pout << " communication topology " << endl;
  vector<int> sendlist;
  //makesendlist(sendlist);

  for (int i = 0; i < sendlist.size(); ++i)
    pout << "processor " << rank << " to " << sendlist[i] << endl;
  pout << "\t\t\t proc " << rank << " of " << size << endl;
#endif
  globaltimer.start();
  int randomseed = 243;
  bool RESTART = false;
  bool reset_iter = false;
  srand(randomseed);
  pout << "\t\t\t random seed " <<  randomseed << endl;  



  double sweep_tol = 1e-7;
  ifstream oneElectronIntegralFile;
  ifstream twoElectronIntegralFile;
  v_1.rhf= true; 
  v_2.rhf=true;
  v_2.permSymm = true;
  std::string configFile;

  int i = 1;
  configFile = argv [i];

  pout << "About to read Input File : "<< configFile<<endl;
  //read the config file
  dmrginp = Input(configFile);
  sweep_tol = dmrginp.get_sweep_tol();
  RESTART = dmrginp.get_restart();
  restartwarm = dmrginp.get_restart_warm();
  reset_iter = dmrginp.get_reset_iterations();
  oneElectronIntegralFile.open(dmrginp.get_oneintegral().c_str(), ios::in);
  twoElectronIntegralFile.open(dmrginp.get_twointegral().c_str(),ios::in);

  pout << "About to read integrals"<<endl;
  //read integrals
  if (mpigetrank() == 0)
  {
    cout << "v2bin "<<v_2.bin<<endl;
    v_1.ReadFromDumpFile(oneElectronIntegralFile);
    pout << "finished v1read" << endl;
    v_2.ReadFromDumpFile(twoElectronIntegralFile, v_1.NOrbs());
    pout << "finished v2read" << endl;
  }
  
#ifndef SERIAL
  mpi::broadcast(world,v_1,0);
  mpi::broadcast(world,v_2,0);
#endif

  //read the config file
  dmrginp = Input(configFile);

#ifndef SERIAL
  MAX_THRD=dmrginp.thrds_per_node()[rank];
#else
  MAX_THRD=dmrginp.thrds_per_node()[0];
#endif

  omp_set_num_threads(MAX_THRD);


  //initialise the size of all Slater determinants equal to orbsize
  Orbstring::init(dmrginp.last_site());

  bool direction;
  int restartsize;
  SweepParams sweepParams;
  sweepParams.restorestate(direction, restartsize);
  Sweep::do_one(sweepParams, false, !direction, false, 0);
  Sweep::do_one(sweepParams, false, direction, false, 0);

#ifndef SERIAL
  if(mpiparent!=MPI_COMM_NULL){
    int itmp;
    MPI_Comm_size(MPI_COMM_WORLD, &itmp);
    if(mpigetrank() == 0) MPI_Send(&itmp, 1, MPI_INT, 0, 0, mpiparent);
    MPI_Barrier(mpiparent);
    MPI_Comm_free(&mpiparent);
  }  
#endif

  return 0;

}


