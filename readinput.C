#include "IntegralMatrix.h"
#include <fstream>
#include "input.h"
#include "pario.h"
#include "global.h"
#include "orbstring.h"
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

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace SpinAdapted;
void ReadInput(char* conf)
{
#ifndef SERIAL
  mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  if (dmrginp.outputlevel() != 0) 
    pout << " communication topology " << endl;
  vector<int> sendlist;
  //makesendlist(sendlist);

  if (dmrginp.outputlevel() != 0) {
    for (int i = 0; i < sendlist.size(); ++i)
      pout << "processor " << rank << " to " << sendlist[i] << endl;
    pout << "\t\t\t proc " << rank << " of " << size << endl;
  }
#endif
  globaltimer.start();
  int randomseed = 243;
  srand(randomseed);
  //pout << "\t\t\t random seed " <<  randomseed << endl;  



  ifstream oneElectronIntegralFile;
  ifstream twoElectronIntegralFile;
  std::string configFile(conf);

  //pout << "About to read Input File : "<< configFile<<endl;
  //read the config file
  pout << "INPUT FILE"<<endl;
  dmrginp = Input(configFile);

  v_1.rhf= true; 
  v_2.rhf=true;
  if (sym != "dinfh")
    v_2.permSymm = true;

  RESTART = dmrginp.get_restart();
  FULLRESTART = dmrginp.get_fullrestart();
  restartwarm = dmrginp.get_restart_warm();
  reset_iter = dmrginp.get_reset_iterations();
  oneElectronIntegralFile.open(dmrginp.get_oneintegral().c_str(), ios::in);
  twoElectronIntegralFile.open(dmrginp.get_twointegral().c_str(),ios::in);

  //pout << "About to read integrals"<<endl;
  //read integrals
  if (mpigetrank() == 0)
  {
    //cout << "v2bin "<<v_2.bin<<endl;
    v_1.ReadFromDumpFile(oneElectronIntegralFile);
    //pout << "finished v1read" << endl;
    v_2.ReadFromDumpFile(twoElectronIntegralFile, v_1.NOrbs());
    //pout << "finished v2read" << endl;
  }

#ifndef SERIAL
  mpi::broadcast(world,v_1,0);
  mpi::broadcast(world,v_2,0);
#endif

#ifndef SERIAL
  MAX_THRD=dmrginp.thrds_per_node()[rank];
#else
  MAX_THRD=dmrginp.thrds_per_node()[0];
#endif
  //pout << "before omp "<<MAX_THRD<<endl;

  omp_set_num_threads(MAX_THRD);
  //pout << "after omp"<<endl;

  //initialise the size of all Slater determinants equal to orbsize
  Orbstring::init(dmrginp.slater_size());

}

