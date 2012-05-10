#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

int callDmrg(char*, char*);

int main(int argc, char* argv []) 
{
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
#endif
  char* output = 0;
  return callDmrg(argv[1], output);
}
