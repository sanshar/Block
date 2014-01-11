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
assert(false);
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

    std::vector< std::vector<SpinQuantum> > cc__dd_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c__cd_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c__c_dd_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > cc_d__d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_cd__d_quantum_ladder;

//??CONTINUE HERE
//////    // Create (CC)(DD) structure
//////    //----------------------------------------
//////    // DD => CC => minus sign and transpose and commute (factor -1??)
//////    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
//////    for (int p=0; p < spinvec12.size(); p++) {
//////      assert( spinvec12[p].get_s() == 2*p );
//////      assert( spinvec12[p].particleNumber == 2 );
//////      // DD_Cre => minus, plus sign ?? FIXME
//////      std::vector<SpinQuantum> spinvec123 = -spinvec12[p] + spin3;
//////      for (int q=0; q < spinvec123.size(); q++) {
//////        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
//////        dd_c_quantum_ladder.push_back( tmp );
//////        assert( spinvec123[q].particleNumber == -1 );
//////      }
//////    }
//////    assert( dd_c_quantum_ladder.size() == 3 );
//////
//////    // Create D(DC) structure
//////    //----------------------------------------
//////    // DC => minus, plus sign
//////    std::vector<SpinQuantum> spinvec23 = -spin2 + spin3;
//////    for (int p=0; p < spinvec23.size(); p++) {
//////      assert( spinvec23[p].get_s() == 2*p );
//////      assert( spinvec23[p].particleNumber == 0 );
//////      // D_(DC) => minus, minus sign ?? //FIXME
//////      std::vector<SpinQuantum> spinvec123 = -spin1 - spinvec23[p];
//////      for (int q=0; q < spinvec123.size(); q++) {
//////        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
//////        d_dc_quantum_ladder.push_back( tmp );
//////        assert( spinvec123[q].particleNumber == -1 );
//////      }
//////    }
//////    assert( d_dc_quantum_ladder.size() == 3 );
//////
//////    // Allocate new operator for each spin component
//////    //------------------------------------------------
//////    spin_ops.clear();
//////    for (int q=0; q < dd_c_quantum_ladder.size(); q++) {
//////      spin_ops.push_back( boost::shared_ptr<CreCreDesDes>(new CreCreDesDes) );
//////      boost::shared_ptr<CreCreDesDes> op = spin_ops.back();
//////      op->set_orbs() = orbs;
//////      op->set_initialised() = true;
//////      // 3-index ops are fermionic 
//////      op->set_fermion() = true;
//////      op->set_quantum_ladder()["((DD)C)"] = dd_c_quantum_ladder.at(q);
//////      op->set_quantum_ladder()["(D(DC))"] = d_dc_quantum_ladder.at(q);
////////FIXME  This should be updated when the build_pattern changes
////////FIXME  Set default
//////      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
//////    }
//////
//////    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================

}
