/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_OPERATOR_LOOPS_HEADER_H
#define SPIN_OPERATOR_LOOPS_HEADER_H
#include <vector>
#include <iostream>
#include <communicate.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <boost/shared_ptr.hpp>
#include "pario.h"


/**
 * Distributed loops to be used functors on OperatorArrays
 * 
 */

/**
 * Loop over all local (i.e. stored on current processor) in array
 * 
 */


namespace SpinAdapted{
class SpinBlock;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//loop over all local operators and build them
// different version include multithread build, single thread build and single thread build from csf

template<class A> void singlethread_build(A& array, SpinBlock& b, std::vector< Csf >& s, vector< vector<Csf> >& ladders)
{
#ifdef _OPENMP
#pragma omp parallel default(shared)
#pragma omp for schedule(guided) nowait
#endif
  for (int i = 0; i < array.get_size(); ++i) {
    //typedef typename A::OpType Op;
    std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
    for (int j=0; j<vec.size(); j++) {
      vec[j]->buildUsingCsf(b, ladders, s);
      //pout << *vec[j]<<endl;
    }
  }
}

template<class A> void singlethread_build(A& array, SpinBlock& b)
{
#ifdef _OPENMP
#pragma omp parallel default(shared)
#pragma omp for schedule(guided) nowait
#endif
  for (int i = 0; i < array.get_size(); ++i) {
    //typedef typename A::OpType Op;
    //std::vector<boost::shared_ptr<Op> >& vec = array.get_local_element(i);
    std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
    for (int j=0; j<vec.size(); j++) {
      vec[j]->build(b);
      //pout << *vec[j]<<endl;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//execute a function on all elements of an array
//single thread and multithread versions of the code
template<typename T2, class A> void for_all_singlethread(A& array, const T2& func)
{
  int i;
  {
    for (i = 0; i < array.get_size(); ++i) {
      std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
      func(vec);
    }
  }
}

template<typename T2, class A> void for_all_multithread(A& array, const T2& func)
{
#ifdef _OPENMP
#pragma omp parallel default(shared)
#pragma omp for schedule(guided) nowait
#endif
    for (int i = 0; i < array.get_size(); ++i) {
      std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
      func(vec);
    }
}
 
template<typename T2, class A> void for_all_operators_multithread(A& array, const T2& func)
{
  int i;
#ifdef _OPENMP
  #pragma omp parallel default(shared) private(i)
#endif
  {
#ifdef _OPENMP
    #pragma omp for schedule(guided) nowait
#endif
    for (i = 0; i < array.get_size(); ++i) {
      std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
      for (int j=0; j<vec.size(); j++){
	func(*vec[j]);
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Used with functions designed to build operators in core, but here we actually write to disk instead
template<typename T2, class A> void for_all_operators_to_disk(A& array, const SpinBlock& b, std::ofstream& ofs, const T2& func)
{     
  for (int i = 0; i < array.get_size(); ++i) {
    std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
    for (int j=0; j<vec.size(); j++){
      // Don't build if already built!
      assert( ! vec.at(j)->get_built() );
      assert( ! vec.at(j)->get_built_on_disk() );
        
      // Apply function to operator
      func( *(vec.at(j)) );
      
      // Store on disk
      vec.at(j)->set_built_on_disk() = true;
      boost::archive::binary_oarchive save_op(ofs);
      save_op << *(vec.at(j));
     
      // Deallocate memory for operator representation
      vec.at(j)->set_built() = false;
//FIXME is Clear() sufficient to deallocate?
      vec.at(j)->Clear();
    }
  }
}   

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
#endif


