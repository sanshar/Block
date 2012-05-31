#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <pario.h>
#include <cstdlib>

using namespace std;
int callDmrg(char*, char*);

int main(int argc, char* argv []) 
{
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
#endif
  if (argc == 1) {
    pout << "No command line argument found, expects an input file"<<endl;
    abort();
  }
  char* output = 0;
  return callDmrg(argv[1], output);
}
