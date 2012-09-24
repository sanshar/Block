/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


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
#include <boost/filesystem.hpp>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace SpinAdapted;
void CheckFileExistance(string filename, string filetype)
{
  boost::filesystem::path p(filename);
  if (boost::filesystem::exists(p)) {
    if (!boost::filesystem::is_regular_file(p)) {
      pout << filetype<<" "<<filename<<" is not a regular file."<<endl;
      abort();
    }
  }
  else {
    pout << filetype<<" "<<filename<<" is not present."<<endl;
    abort();
  }
}

void ReadInput(char* conf)
{
#ifndef SERIAL
  mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  if (dmrginp.outputlevel() > 0) 
    pout << " communication topology " << endl;
  vector<int> sendlist;

  if (dmrginp.outputlevel() > 0) {
    for (int i = 0; i < sendlist.size(); ++i)
      pout << "processor " << rank << " to " << sendlist[i] << endl;
      pout << "\t\t\t proc " << rank << " of " << size << endl;
  }
#endif
  globaltimer.start();
  int randomseed = 243;
  srand(randomseed);


  std::string configFile(conf);

  CheckFileExistance(conf, "Input file ");
  //read the config file
  dmrginp = Input(configFile);

  RESTART = dmrginp.get_restart();
  FULLRESTART = dmrginp.get_fullrestart();
  restartwarm = dmrginp.get_restart_warm();
  reset_iter = dmrginp.get_reset_iterations();


#ifndef SERIAL
  MAX_THRD=dmrginp.thrds_per_node()[rank];
#else
  MAX_THRD=dmrginp.thrds_per_node()[0];
#endif

  omp_set_num_threads(MAX_THRD);

  //initialise the size of all Slater determinants equal to orbsize
  Orbstring::init(dmrginp.slater_size());

}

