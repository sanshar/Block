/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef COMMUNICATE_HEADER_H
#define COMMUNICATE_HEADER_H
#ifndef SERIAL
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>

inline int mpigetrank() { boost::mpi::communicator world;  return world.rank(); }
inline int mpigetsize() { boost::mpi::communicator world; return world.size(); }
#else
inline int mpigetrank() { return 0; }
inline int mpigetsize() { return 1; }
#endif

#ifdef SERIAL 
template<class T> void sendobject(const T& object, int to) {}
#else
template<class T> void sendobject(const T& object, int to)  
{
  boost::mpi::communicator world;
  int tag = 0;
  world.send(to,tag,object);
}
#endif
// default argument + template specialization bug workaround
//template<class T> void receiveobject(T& object, int from, int tag = 0)
#ifdef SERIAL
template<class T> void receiveobject(T& object, int from) {}
#else
template<class T> void receiveobject(T& object, int from)
{
  boost::mpi::communicator world;
  int tag = 0;
  world.recv(from,tag,object);
}
#endif

#endif
