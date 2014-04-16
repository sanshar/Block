#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "global.h"
#include "fciqmchelper.h"
#include "input.h"
#include <pario.h>

void ReadInput(char* conf);
int SpinAdapted::MPS::sweepIters ;
bool SpinAdapted::MPS::spinAdapted ;

int main(int argc, char* argv []) {

#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
#endif

  if (argc == 1 && mpigetrank() == 0) {
    cout << "No command line argument found, expects an input file"<<endl;
    abort();
  }
  char* output = 0;
  ReadInput(argv[1]);
  
  int statea = 0, stateb = 0;


  //Initialize the stateinfo, overlap for site1+2 
  MPS::sweepIters = dmrginp.last_site()/2-2;
  MPS::spinAdapted = false;

  MPS dmrga(statea);
  MPS dmrgb(stateb);

  std::vector<bool> occ(dmrginp.last_site(),0);
 
  occ[0] = 1; occ[1] = 1;
  //occ[4] = 1; occ[9] = 1; ///
  occ[5] = 1; occ[8] = 1; ///
  occ[14] = 1; occ[15] = 1;

  //occ[4] = 1; occ[5] = 1;
  //occ[6] = 1; occ[7] = 1;
  //occ[8] = 1; occ[9] = 1;
  //occ[10] = 1; occ[11] = 1;
  //occ[12] = 1; occ[13] = 1;
  //occ[14] = 1; occ[15] = 1;
  //occ[16] = 1; occ[17] = 1;
  //occ[18] = 1; occ[19] = 1;
  //occ[20] = 1; occ[21] = 1;
  //occ[22] = 1; occ[23] = 1;

  for (int i=0; i<dmrginp.last_site(); i++)
    cout << occ[i]<<"  ";
  cout << endl;

  MPS dmrgc(occ);

  cout <<"overlap = "<< calculateOverlap(dmrga, dmrgc)<<endl;
  //cout <<"overlap = "<< calculateOverlap(dmrga, dmrgb)<<endl;
  return 0;
}
