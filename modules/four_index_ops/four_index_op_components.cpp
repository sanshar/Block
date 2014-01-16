/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
//
//  This is an extension of op_components.C (i.e. used with op_components.h header file)
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

#include <boost/format.hpp>
#include <fstream>
#include <stdio.h>
#include "spinblock.h"
#include "op_components.h"
//#include "screen.h"

namespace SpinAdapted {
  
//===========================================================================================================================================================
// Choose 4-index tuples on this MPI process such that 2-index are available to build them
// FIXME IMPLEMENT FOR NON-PARALLEL VERSION FIRST
// There are several cases????

//FIXME screening?
std::map< std::tuple<int,int,int,int>, int > get_4index_tuples(SpinBlock& b)
{
  std::map< std::tuple<int,int,int,int>, int > tuples;

//FIXME
//  if ( b.get_leftBlock() != NULL ) {
//    // Generate only mpi local tuples for compound block, consistent with existing operators on sys and dot
//    tuples = get_local_3index_tuples(b);    
//  }

//FIXME DEBUG ONLY
//  // (l == k == j == i)
//  std::vector<int> sites = b.get_sites();
//  for (int i = 0; i < sites.size(); ++i)
//    for (int j = i; j <= i; ++j)
//      for (int k = i; k <= i; ++k)
//        for (int l = i; l <= i; ++l)
//          tuples[ std::make_tuple(sites[i], sites[j], sites[k], sites[l]) ] = -1;


  // Generate all tuples such that (l <= k <= j <= i) and let para_array assign them to local processes as necessary
  std::vector<int> sites = b.get_sites();
  for (int i = 0; i < sites.size(); ++i)
    for (int j = 0; j <= i; ++j)
      for (int k = 0; k <= j; ++k) {
        for (int l = 0; l <= k; ++l) {
//          if ( b.get_leftBlock() != NULL ) {
//            // The -2 here means that this should be assigned to global_indices only (i.e. shouldn't be on this MPI thread)
//            if ( tuples.find(std::make_tuple(sites[i], sites[j], sites[k], sites[l])) == tuples.end() )
//              tuples[ std::make_tuple(sites[i], sites[j], sites[k], sites[l]) ] = -2;
//          }
          // The -1 here means there's no constraint on which MPI process (i.e. let para_array choose if it should belong to this one)
          tuples[ std::make_tuple(sites[i], sites[j], sites[k], sites[l]) ] = -1;
        }
      }

  return tuples;
}

//===========================================================================================================================================================
// RI_4_INDEX skeleton class
//----------------------------

template<>
string Op_component<RI4index>::get_op_string() const {
  return "RI_4_INDEX";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<RI4index>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return;

  // Set up 4-index (i,j,k,l) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int,int>, int > tuples = get_4index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

//  // Allocate new set of operators for each set of spatial orbitals
//  std::vector<int> orbs(4);
//  for (int i = 0; i < m_op.local_nnz(); ++i) {
//    orbs = m_op.unmap_local_index(i);
//  }
}

//===========================================================================================================================================================
// (Cre,Cre,Des,Des)
//-------------------

template<> 
string Op_component<CreCreDesDes>::get_op_string() const {
  return "CreCreDesDes";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<CreCreDesDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 4-index (i,j,k,l) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int,int>, int > tuples = get_4index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(4);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
pout << "New set of CCDD operators:  " << i << std::endl;
pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << " " << orbs[3] << std::endl;
    std::vector<boost::shared_ptr<CreCreDesDes> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));
    SpinQuantum spin4 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[3]));

    // Create (CC)(DD) structure
    //----------------------------------------
    std::vector< std::vector<SpinQuantum> > cc_dd_quantum_ladder;
    // (DD)
    std::vector<SpinQuantum> spinvec34 = spin3 + spin4;
    for (int p=0; p < spinvec34.size(); p++) {
      assert( spinvec34[p].get_s() == 2*p );
      assert( spinvec34[p].particleNumber == 2 );
      // (CC)
      std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
      for (int q=0; q < spinvec12.size(); q++) {
        std::vector<SpinQuantum> spinvec1234 = spinvec12[q] - spinvec34[p];
        for (int t=0; t < spinvec1234.size(); t++) {
          // Store (CC) first, then (DD), then 4-index spin quantums
          std::vector<SpinQuantum> tmp = { spinvec12[q], spinvec34[p], spinvec1234[t] };
cout << "p, spin12[q]   = " << q << ",  " << 0.5*spinvec12[q].get_s() << endl;
cout << "q, spin34[p]   = " << p << ",  " << 0.5*spinvec34[p].get_s() << endl;
cout << "t, spin1234[t] = " << t << ",  " << 0.5*spinvec1234[t].get_s() << endl;
          cc_dd_quantum_ladder.push_back( tmp );
          assert( spinvec1234[t].particleNumber == 0 );
        }
      }
    }
    assert( cc_dd_quantum_ladder.size() == 6 );

    // Create ((CC)D)(D) structure
    //----------------------------------------
    std::vector< std::vector<SpinQuantum> > cc_d__d_quantum_ladder;
    // (CC)
    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (CC)D
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] - spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> spinvec1234 = spinvec123[q] - spin4;
        // Store (CC) first, then (CC)D, then 4-index spin quantums
        for (int t=0; t < spinvec1234.size(); t++) {
          std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q], spinvec1234[t] };
//cout << "p, spin12[p]   = " << p << ",  " << 0.5*spinvec12[p].get_s() << endl;
//cout << "q, spin123[p]  = " << q << ",  " << 0.5*spinvec123[q].get_s() << endl;
//cout << "t, spin1234[t] = " << t << ",  " << 0.5*spinvec1234[t].get_s() << endl;
          cc_d__d_quantum_ladder.push_back( tmp );
          assert( spinvec1234[t].particleNumber == 0 );
        }
      }
    }
    assert( cc_d__d_quantum_ladder.size() == 6 );


    // Create (C(CD))(D) structure
    //----------------------------------------
    std::vector< std::vector<SpinQuantum> > c_cd__d_quantum_ladder;
    // (CD)
    std::vector<SpinQuantum> spinvec23 = spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // C(CD)
      std::vector<SpinQuantum> spinvec123 = spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> spinvec1234 = spinvec123[q] - spin4;
        // Store (CD) first, then C(CD), then 4-index spin quantums
        for (int t=0; t < spinvec1234.size(); t++) {
          std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q], spinvec1234[t] };
          c_cd__d_quantum_ladder.push_back( tmp );
          assert( spinvec1234[t].particleNumber == 0 );
        }
      }
    }
    assert( c_cd__d_quantum_ladder.size() == 6 );

    // Create (C)(C(DD)) structure
    //----------------------------------------
    std::vector< std::vector<SpinQuantum> > c__c_dd_quantum_ladder;
    // (DD)
    for (int p=0; p < spinvec34.size(); p++) {
      // C(DD)
      std::vector<SpinQuantum> spinvec234 = spin2 - spinvec34[p];
      for (int q=0; q < spinvec234.size(); q++) {
        std::vector<SpinQuantum> spinvec1234 = spin1 + spinvec234[q];
        // Store (DD) first, then C(DD), then 4-index spin quantums
        for (int t=0; t < spinvec1234.size(); t++) {
          std::vector<SpinQuantum> tmp = { spinvec34[p], spinvec234[q], spinvec1234[t] };
          c__c_dd_quantum_ladder.push_back( tmp );
          assert( spinvec1234[t].particleNumber == 0 );
        }
      }
    }
    assert( c__c_dd_quantum_ladder.size() == 6 );

    // Create (C)((CD)D) structure
    //----------------------------------------
    std::vector< std::vector<SpinQuantum> > c__cd_d_quantum_ladder;
    // (CD)
    for (int p=0; p < spinvec23.size(); p++) {
      // (CD)D
      std::vector<SpinQuantum> spinvec234 = spinvec23[p] - spin4;
      for (int q=0; q < spinvec234.size(); q++) {
        std::vector<SpinQuantum> spinvec1234 = spin1 + spinvec234[q];
        // Store (CD) first, then (CD)D, then 4-index spin quantums
        for (int t=0; t < spinvec1234.size(); t++) {
          std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec234[q], spinvec1234[t] };
          c__cd_d_quantum_ladder.push_back( tmp );
          assert( spinvec1234[t].particleNumber == 0 );
        }
      }
    }
    assert( c__cd_d_quantum_ladder.size() == 6 );


    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < cc_dd_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<CreCreDesDes>(new CreCreDesDes) );
      boost::shared_ptr<CreCreDesDes> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((CC)(DD))"] = cc_dd_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C)((CD)D))"] = c__cd_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C)(C(DD)))"] = c__c_dd_quantum_ladder.at(q);
      op->set_quantum_ladder()["(((CC)D)(D))"] = cc_d__d_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C(CD))(D))"] = c_cd__d_quantum_ladder.at(q);
      // Set default value, which is changed according to current build_pattern // FIXME .at(2) is brittle here!  
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(2);
    }

    assert( m_op.get_local_element(i).size() == 6);
  }
}

//===========================================================================================================================================================

}
