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
#ifdef _OPENMP
#include "omp.h"
#endif
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
void CheckFileExistence(string filename, string filetype)
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
void CheckFileInexistence(string filename, string filetype){
   boost::filesystem::path p(filename);
   if (boost::filesystem::exists(p)) {
      p2out << filetype<<" "<<filename<<" is present."<<endl;
      if (filetype=="genetic algorithm reorder")
         pout << "You already ran the genetic algorithm reordering and the ordering is present in the location above." << endl;
      abort();
   }
}



void ReadInput(char* conf)
{
#ifndef SERIAL
  mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  p3out << " communication topology " << endl;
  vector<int> sendlist;

  for (int i = 0; i < sendlist.size(); ++i)
    p3out << "processor " << rank << " to " << sendlist[i] << endl;
    p3out << "\t\t\t proc " << rank << " of " << size << endl;
#endif
  globaltimer.start();
  int randomseed = 243;
  srand(randomseed);


  std::string configFile(conf);

  CheckFileExistence(conf, "Input file ");
  //read the config file
  dmrginp = Input(configFile);

  RESTART = dmrginp.get_restart();
  FULLRESTART = dmrginp.get_fullrestart();
  BACKWARD = dmrginp.get_backward();
  restartwarm = dmrginp.get_restart_warm();
  reset_iter = dmrginp.get_reset_iterations();


#ifndef SERIAL
  MAX_THRD=dmrginp.thrds_per_node()[rank];
#else
  MAX_THRD=dmrginp.thrds_per_node()[0];
#endif

#ifdef _OPENMP
  omp_set_num_threads(MAX_THRD);
#endif

  //initialise the size of all Slater determinants equal to orbsize
  Orbstring::init(dmrginp.slater_size());

}

