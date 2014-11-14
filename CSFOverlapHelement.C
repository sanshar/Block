#include "wrapper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "fciqmchelper.h"
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "global.h"

using namespace std;

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
  
  dmrginp.add_noninteracting_orbs() = false;
  readMPSFromDiskAndInitializeStaticVariables();
  //initializeGlobalMPS(mpsstate);


  double overlap, hvalue;

  unsigned long temp=1, occ=0;
  int msgsize=1000;
  char msgctr[msgsize];

  ifstream file("determinants");
  file.getline(msgctr, msgsize);
  string s( msgctr);
  vector<string> tok;
  boost::split(tok, s, is_any_of(", \t"), token_compress_on);
  int ndets = atoi(tok[0].c_str());


  std::vector<MPS> dets;
  for (int i=0; i<ndets; i++) {
    file.getline(msgctr, msgsize);
    //string occstring = "1 1 0 0 1 1 0 0 1 1 0 0 0 0 1 1";
    intFromString(occ, msgctr);
    dets.push_back(MPS(&occ, 1));
  } 

  /*
  for (int i=0; i<ndets; i++)
    for (int j=i; j<ndets; j++) {
      calcHamiltonianAndOverlap(dets[i], dets[j], hvalue, overlap);
      printf("<D%i|H|D%i> = %18.10f \n", i, j, hvalue);
    }
  */
  printf ("Interaction Matrix, Row(i,j) = <Di|H|Dj>\n");
  for (int i=0; i<ndets; i++) {
    for (int j=0; j<ndets; j++) {
      calcHamiltonianAndOverlap(dets[i], dets[j], hvalue, overlap);
      printf(" %12.3e  ", hvalue);
    }
    printf("\n");
  }


#ifndef SERIAL
  MPI_Finalize();
#endif
  return 0;
}
