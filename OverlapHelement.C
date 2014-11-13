#include "wrapper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
#include "stdio.h"
#include "stdlib.h"
#include <sstream>

int main(int argc, char* argv []) {

  int rank=0, size=1;
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  initBoostMPI(argc, argv);
  ReadInputFromC(argv[1], -1);

  int mpsstate=0;
  
  readMPSFromDiskAndInitializeStaticVariables(mpsstate);


  test();
  //exit(0);
  //
  //
  /*
  double overlap, hvalue;

  long temp=1, occ=0;

  string occstring = "11 00 11 00 11 00 00 11";
  stringstream stream(occstring);
  int n, i=0;
  while (strem >>n) {
    if (n==1)
      occ = occ | temp <<(63-i);
    i++;
  }

  evaluateOverlapAndHamiltonian(&occ, 1, &overlap, &hvalue);
 
  printf("overlap = %15.8e  %i\n", overlap, rank);
  printf("helement = %15.8e  %i\n", hvalue, rank);
  */

#ifndef SERIAL
  MPI_Finalize();
#endif
  return 0;
}
