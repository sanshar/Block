#include "wrapper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
#include "stdio.h"

int main(int argc, char* argv []) {

  int rank=0, size=1;
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  ReadInputFromC(argv[1], -1);
  int mpsstate=0;
  
  readMPSFromDiskAndInitializeStaticVariables(mpsstate);
  

  double overlap, hvalue;

  long temp=1, occ;
  
  //occ = temp<<63 | temp<<62 | temp<<61 | temp<<60 | temp<<59 | temp<<58 ;
  occ = temp<<63 | temp<<62 | temp<<61 | temp<<60 | temp<<59 | temp<<58 | temp<<57 | temp<<56 | temp<<55 | temp<<54 ;

  evaluateOverlapAndHamiltonian(&occ, 1, &overlap, &hvalue);
 
  printf("overlap = %10.5e  %i\n", overlap, rank);
  printf("helement = %10.5e  %i\n", hvalue, rank);

#ifndef SERIAL
  MPI_Finalize();
#endif
  return 0;
}
