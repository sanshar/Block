
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <stdio.h>
#include "nevpt2_mpi.h"

int mpi_world_size(){
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
  return size;
#else
  return 1;
#endif  
}

int mpi_rank(){
#ifndef SERIAL
  boost::mpi::communicator world;
  int rank = world.rank();
  return rank;
#else
  return 0;
#endif  
}

void mpi_barrier(){
#ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
#else
  
#endif  
}

void mpi_divide_loop(int &start_loc, int &stop_loc, int start_global, int stop_global){
#ifndef SERIAL
    boost::mpi::communicator world;
    if (world.size()==1){
      start_loc = start_global;
      stop_loc   = stop_global;
    }
    else{
      int loopsize = stop_global - start_global + 1;
      int numproc = world.size();
      int quot = loopsize / numproc;
      int iproc = world.rank();
      start_loc = iproc * quot;
      stop_loc = (iproc+1) * quot;
      if (world.rank()==numproc-1){
        stop_loc = stop_global;
      }
    }
#else
    start_loc = start_global;
    stop_loc   = stop_global;
#endif
  }

void mpi_debug_break(){
#ifndef SERIAL
  int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  char msg[512];
  sprintf(msg,"\nPID %d on %s ready to be attached to\n", getpid(), hostname);
  printf("%s",msg);
  fflush(stdout);
  while (0 == i) sleep(5);
  mpi_barrier();
#endif
}

void mpi_message(char *msg){
  int rank = mpi_rank();
  char MSG[512];
  sprintf(MSG,"%s rank=%i",msg,rank);
  printf("%s",MSG);
  fflush(stdout);
  }

  
void PALDivideLoop(int &start_loc, int &stop_loc, int start_global, int stop_global){
#ifndef SERIAL
  boost::mpi::communicator world;
  if (world.size()==1){
    start_loc = start_global;
    stop_loc   = stop_global;
  }
  else{
    int loopsize = stop_global - start_global + 1;
    int numproc = world.size();
    int quot = loopsize / numproc;
    int iproc = world.rank();
    start_loc = iproc * quot;
    stop_loc = (iproc+1) * quot;
    if (numproc>=loopsize){
      if (iproc>=loopsize){
        start_loc = -1;
        stop_loc = -1;
      }
      else{
        start_loc = iproc;
        stop_loc = iproc + 1;
      }
    }
    if (world.rank()==numproc-1){
      stop_loc = stop_global;
    }
      
  }
#else
  start_loc = start_global;
  stop_loc   = stop_global;
#endif
}


//determine the (maximum) global dimension
void LoopSynchronizer::DetermineGlobalDim(){
  //if this is a serial calculation, set the local dimension as global maximum
  if (mpi_world_size()==1){
    MaxDimGlobal = LocalDim;
    HaveGlobalDim = true;
    return;
  }
  //if not, determine the maximum global dimension
  int tmp_recv;
#ifndef SERIAL
  boost::mpi::communicator world;
  int world_size = world.size();
  int rank = world.rank();
  if (rank==0){
    //set the local dimension as the maximum value
    int MaxDim = LocalDim;
    //gather all local dimensions and determine their maximum value
    for (int iproc=1;iproc<world_size;iproc++){
      world.recv(iproc, iproc, tmp_recv);
      if (tmp_recv>MaxDim) MaxDim=tmp_recv;
    }//iproc
    MaxDimGlobal = MaxDim;
  }//master
  else{
    world.send(0, rank, LocalDim);
  }//slave
  //broadcast the result
  boost::mpi::broadcast(world,MaxDimGlobal,0);
  HaveGlobalDim = true;
#endif
}

