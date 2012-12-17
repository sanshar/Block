/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_DISTRIBUTE_HEADER_H
#define SPIN_DISTRIBUTE_HEADER_H
#include <iostream>
#include <communicate.h>
#include "timer.h"
#include "pario.h"
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <vector>

namespace SpinAdapted
{
  class SparseMatrix;


  // tree mpi::broadcast algorithm ...
  void makesendlist(std::vector<int>& tolist, int offsetproc = 0);

  int receivefrom(int offsetproc = 0);

  namespace Broadcastsettings
  {
    const int npyramid = 2;
  };
#ifdef SERIAL
  template<class T> void distributedaccumulate(T& component){;}
#else
  template<class T> void plus_equals(T& lhs, const T& rhs)
    {
      lhs += rhs;
    }
  template<class T> void distributedaccumulate(T& component)
    {
      Timer distributetimer;
      boost::mpi::communicator world;
      int size = world.size();
      int rank = world.rank();
      if (size > 1)
	{
	  // i.  look at our sendlist
	  // ii. receive stuff from all the nodes in our sendlist
	  // iii. accumulate
	  // iv. send to receivefromnode
	  // v. mpi::broadcast
      
	  std::vector<int> sendlist;
	  makesendlist(sendlist);
      
	  // receive and accumulate
	  int listsize = sendlist.size();
	  for (int i = 0; i < listsize; ++i)
	    {
	      T part;
	      receiveobject(part, sendlist[i]);
	      plus_equals(component, part);
	      //component = component + part;	  
	    }
      
	  // send to root (if not already root)
	  if (rank != 0)
	    sendobject(component, receivefrom());

	  boost::mpi::broadcast(world,component,0);
	}
      //     pout << "distribution time " << distributetimer.elapsedwalltime() << " " << distributetimer.elapsedcputime() << endl;
    }
#endif
#ifdef SERIAL
  template<class T> void accumulateto(T& component, int proc){;}
#else
  template<class T> void accumulateto(T& component, int proc)
    {
      Timer distributetimer;
      boost::mpi::communicator world;
      int size = world.size();
      int rank = world.rank();
      if (size > 1)
	{
	  // i.  look at our sendlist
	  // ii. receive stuff from all the nodes in our sendlist
	  // iii. accumulate
	  // iv. send to receivefromnode
	  // v. mpi::broadcast
      
	  std::vector<int> sendlist;
	  makesendlist(sendlist);
      
	  // receive and accumulate
	  int listsize = sendlist.size();
	  for (int i = 0; i < listsize; ++i)
	    {
	      T part;
	      receiveobject(part, sendlist[i]);
	      component += part;	  
	    }
      
	  // send to root (if not already root)
	  if (rank != 0)
	    sendobject(component, receivefrom());

	  // send back out to proc from root
	  if (rank == 0)
	    sendobject(component, proc);
	  else if (rank == proc)
	    receiveobject(component, 0);
	}
      //  pout << "distribution time " << distributetimer.elapsedwalltime() << " " << distributetimer.elapsedcputime() << endl;
    }



#endif


  //template<> void distributedaccumulate(SparseMatrix& component);
  //template<> void accumulateto(SparseMatrix& component, int proc);

  template<class T> void initiateMultiThread( T* op, T* &op_array, T* &op_distributed, int MAX_THRD)
    {
#ifndef SERIAL
      boost::mpi::communicator world;
      int size = world.size();

      if (size == 1 && MAX_THRD == 1) {
	op_array = op;
	op_distributed = op;
      }
      else if (size > 1 && MAX_THRD == 1) {
	op_array = op;
	op_distributed = new T;
	*op_distributed = *op;
      }
      else if (size == 1 && MAX_THRD > 1) {
	op_array = new T[MAX_THRD];
	for (int i=0; i<MAX_THRD; i++)
	  op_array[i] = *op;
	op_distributed = op_array;
      }
      else {
	//op_array = new T[MAX_THRD];
	//for (int i=0; i<MAX_THRD; i++)
	//op_array[i] = *op;
	op_distributed = new T[MAX_THRD];
	for (int i=0; i<MAX_THRD; i++)
	  op_distributed[i] = *op;
      }
#else
      if (MAX_THRD == 1) {
	op_array = op;
	op_distributed = op;
      }
      else
	{
	  op_array = new T[MAX_THRD];
	  for (int i=0; i<MAX_THRD; i++)
	    op_array[i] = *op;
	  op_distributed = op_array;
	}
#endif
    }

  template<class T> void accumulateMultiThread(T* op, T* &op_array, T* &op_distributed, int MAX_THRD, int toproc = 0)
    {
#ifndef SERIAL
      boost::mpi::communicator world;
      int size = world.size();
  
      if (size == 1 && MAX_THRD == 1)
	return;
      else if (size > 1 && MAX_THRD == 1) { //only mpi
	if (toproc == 0) {
	  distributedaccumulate(op_distributed[0]);
	  *op += op_distributed[0];
	}
	else {
	  accumulateto(op_distributed[0], toproc);
	  if(mpigetrank() == toproc) *op += op_distributed[0];
	}
	delete op_distributed;
      }
      else if (size == 1 && MAX_THRD > 1) {  //only multithreaded
	for (int i=0; i<MAX_THRD; i++)
	  *op += op_array[i];
	delete [] op_array;
      }
      else {  //multithreaded and mpi
	//for (int i=0; i<MAX_THRD; i++)
	//*op += op_array[i];
	//delete [] op_array;

	for (int i=1; i<MAX_THRD; i++)
	  op_distributed[0] += op_distributed[i];

	if (toproc == 0) {
	  distributedaccumulate(op_distributed[0]);
	  *op += op_distributed[0];
	}
	else {
	  accumulateto(op_distributed[0], toproc);
	  if (mpigetrank() == toproc) *op += op_distributed[0];
	}
	delete [] op_distributed;    
      }
#else
      if ( MAX_THRD == 1)
	return;
      else {  //only multithreaded
	for (int i=0; i<MAX_THRD; i++)
	  *op += op_array[i];
	delete [] op_array;
      }
#endif

    }
}
#endif
