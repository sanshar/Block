/*                                                                           
Developed by Sandeep Sharma, Roberto Olivares-Amaya and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma, Garnet K.-L. Chan and Roberto Olivares-Amaya
*/


#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <pario.h>
#include <cstdlib>
#include "dmrg_wrapper.h"
//#include "input.h"
//#include "global.h"

using namespace std;
int calldmrg(char*, char*);

#ifndef UNITTEST
int main(int argc, char* argv []) 
{
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  if (world.rank() == 0) {
    pout << "Runing with " << world.size() << " processors" << endl;
  }
#endif

  //This needs to be added
  //SpinAdapted::dmrginp.initCumulTimer();

  if (argc == 1 && mpigetrank() == 0) {
    pout << "No command line argument found, expects an input file"<<endl;
    abort();
  } else if ((!strncmp(argv[1], "-v", 2)) ||
             (!strncmp(argv[1], "--version", 9))) {
    pout << "Block 1.1.1"<<endl;
    pout << "Copyright (C) 2012  Garnet K.-L. Chan"<<endl;
    pout << "This program comes with ABSOLUTELY NO WARRANTY; for details see license file."<<endl;
    pout << "This is free software, and you are welcome to redistribute it"<<endl;
    pout << "under certain conditions; see license file for details."<<endl;
    return 0;
  }
  char* output = 0;
  return calldmrg(argv[1], output);
}
#endif
