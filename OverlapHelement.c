#include "wrapper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
#include "stdio.h"

int main(int argc, char* argv []) {

#ifndef SERIAL
  MPI_Init(&argc, &argv);
#endif

  ReadInputFromC(argv[1], -1);
  int mpsstate=0;
  readMPSFromDiskAndInitializeStaticVariables(mpsstate);

  double overlap, hvalue;

  long temp=1, occ;
  
  occ = temp<<63 | temp<<62 | temp<<61 | temp<<60 | temp<<59 | temp<<58 ;

  evaluateOverlapAndHamiltonian(&occ, 1, &overlap, &hvalue);
 
  printf("overlap = %10.5e\n", overlap);
  printf("helement = %10.5e\n", hvalue);
  return 0;
}
