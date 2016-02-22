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
#include "screen.h"
#include "spinblock.h"
#include "op_components.h"
#include "build_3index_ops.h"
#include "pario.h"

namespace SpinAdapted {
  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
//FIXME update SCREEN.C
//
//std::vector<std::tuple<int,int,int> > screened_ccd_indices(const vector<int, std::allocator<int> >& sites,
//                                                           const vector<int, std::allocator<int> >& interactingix,
//                                                           const TwoElectronArray& twoe, double thresh)
//{
//  std::vector<std::tuple<int,int,int> > screened_indices;
//
//  if ( true ) {
//    // Set up site indices such that (k <= j <= i)
//  //FIXME should this stride order match para_array??
//    for (int i = 0; i < sites.size(); ++i) {
//      for (int j = 0; j <= i; ++j) {
//        for (int k = 0; k <= j; ++k) {
//  //      if (dmrginp.use_partial_two_integrals()) {
//          screened_indices.push_back(std::make_tuple(sites[i], sites[j], sites[k]));
//  //      }
//  //FIXME
//  //      else {
//  //        if (screen_ccd_interaction(indices[i], indices[j], interactingix, twoe, thresh))
//  //          screened_indices.push_back(make_pair(indices[i], indices[j]));
//  //      }
//        }
//      }
//    }
//  }
//  else {
//    // i=k=k
//    for (int i = 0; i < sites.size(); ++i) 
//      screened_indices.push_back(std::make_tuple(sites[i], sites[i], sites[i]));
//  }
//
//  return screened_indices;
//}
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
// Choose 3-index tuples on this MPI process such that 2-index are available to build them

std::map< std::tuple<int,int,int>, int > get_local_3index_tuples(SpinBlock& b)
{
  std::map< std::tuple<int,int,int>, int > tuples;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
  
  assert( dotBlock != NULL );
  if(dmrginp.spinAdapted())
  assert( dotBlock->get_sites().size() == 1 );
  if(!dmrginp.spinAdapted())
  assert( dotBlock->get_sites().size() == 2 );

  bool forward = true;
  if ( sysBlock->get_sites()[0] > dotBlock->get_sites()[0] ) forward = false;
  int dot = dotBlock->get_sites()[0];

  // 3 on dot (the -1 means there's no constraint on which MPI process the tuple is assigned to)
  //----------
  tuples[ std::make_tuple(dot, dot, dot) ] = -1;

  // 2 on dot
  //----------
//FIXME assume CRE and DES same (THEY'RE NOT, CRE is duplicated sometimes in save_load_block.C so use DES)
  std::vector< std::vector<int> > op_array = sysBlock->get_op_array(DES).get_array();
  for (auto op = op_array.begin(); op != op_array.end(); ++op) {
    int i = (*op)[0];
    if ( forward ) {
      if ( sysBlock->get_op_array(DES).is_local() )
        // When 1-index is duplicated on all ranks we don't want 3-index being duplicated too
        tuples[ std::make_tuple( dot, dot, i) ] = -1; 
      else
        tuples[ std::make_tuple( dot, dot, i) ] = mpigetrank();
    } 
    else {
      if ( sysBlock->get_op_array(DES).is_local() )
        // When 1-index is duplicated on all ranks we don't want 3-index being duplicated too
        tuples[ std::make_tuple( i, dot, dot) ] = -1;
      else
        tuples[ std::make_tuple( i, dot, dot) ] = mpigetrank();
    }
  }

  // 1 on dot
  //----------
//FIXME we assume CRE_CRE is representative of all 2-index ops (pass opType in general???)
  std::vector< std::vector<int> > ij_array = sysBlock->get_op_array(CRE_CRE).get_array();
  for (auto ij = ij_array.begin(); ij != ij_array.end(); ++ij) {
    int i = (*ij)[0];
    int j = (*ij)[1];
    assert( i >= j );
    if ( forward ) {
      if ( sysBlock->get_op_array(CRE_CRE).is_local() )
//FIXME not convinced this reliably removes duplicate 3-index ops
        // When 2-index is duplicated on all ranks we don't want 3-index being duplicated too
        tuples[ std::make_tuple( dot, i, j) ] = -1; 
      else
        tuples[ std::make_tuple( dot, i, j) ] = mpigetrank();
    } 
    else {
      if ( sysBlock->get_op_array(CRE_CRE).is_local() )
        // When 2-index is duplicated on all ranks we don't want 3-index being duplicated too
        tuples[ std::make_tuple( i, j, dot) ] = -1;
      else
        tuples[ std::make_tuple( i, j, dot) ] = mpigetrank();
    }
  }

  // 0 on dot
  //----------
  std::vector< std::vector<int> > ijk_array = sysBlock->get_op_array(RI_3INDEX).get_array();
  for (auto ijk = ijk_array.begin(); ijk != ijk_array.end(); ++ijk) {
    int i = (*ijk)[0];
    int j = (*ijk)[1];
    int k = (*ijk)[2];
    assert( i >= j );
    assert( j >= k );
    if ( sysBlock->get_op_array(RI_3INDEX).is_local() )
      // When 1-site 3-index is duplicated on all ranks we don't want multi-site 3-index being duplicated too
      tuples[ std::make_tuple( i, j, k) ] = -1;
    else
      tuples[ std::make_tuple( i, j, k) ] = mpigetrank();
  }

  return tuples;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
// Choose 3-index tuples on this MPI process such that 2-index are available to build them
// There are three cases:
//   (1) 3-index must be done on this thread  (e.g. built from 2-index on this thread)
//   (2) 3-index could be done on this thread if load-balancing wants it 
//   (3) 3-index must not be done on this thread (it has to be done on another and we don't want duplicates)

//FIXME screening?
std::map< std::tuple<int,int,int>, int > get_3index_tuples(SpinBlock& b)
{
  std::map< std::tuple<int,int,int>, int > tuples;
  std::vector<int> sites = b.get_sites();

  //add a special case for when rightblock is a dummyblock
  if (b.get_rightBlock() != NULL) 
    if (b.get_rightBlock()->get_sites().size() == 0) {
      tuples[ std::make_tuple(sites[0], sites[0], sites[0]) ] = -1;
      return tuples;
    }

  if ( b.get_leftBlock() != NULL ) {
    // Generate only mpi local tuples for compound block, consistent with existing operators on sys and dot
    tuples = get_local_3index_tuples(b);    
  }

  // Generate all tuples such that (k <= j <= i) and let para_array assign them to local processes as necessary
  for (int i = 0; i < sites.size(); ++i)
    for (int j = 0; j <= i; ++j)
      for (int k = 0; k <= j; ++k) {
        if ( b.get_leftBlock() != NULL ) {
          // The -2 here means that this should be assigned to global_indices only (i.e. shouldn't be on this MPI thread)
          if ( tuples.find(std::make_tuple(sites[i], sites[j], sites[k])) == tuples.end() )
            tuples[ std::make_tuple(sites[i], sites[j], sites[k]) ] = -2;
        }
        else {
          // The -1 here means there's no constraint on which MPI process (i.e. let para_array choose if it should belong to this one)
          tuples[ std::make_tuple(sites[i], sites[j], sites[k]) ] = -1;
        } 
      }

  return tuples;
}

//===========================================================================================================================================================

//FIXME the 3-index routines below have a lot of common ground that could be re-implemented more clearly

//===========================================================================================================================================================
// RI_3_INDEX skeleton class
//----------------------------

template<> 
string Op_component<RI3index>::get_op_string() const {
  return "RI_3_INDEX";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<RI3index>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  // This is only a skeleton class, so actual operators should never be built
  return;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<RI3index>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                             const std::vector<Matrix>& rightMat, const StateInfo *ket)
{ 
  //FIXME
  //No need to renormalise RI3index, since it is just gotten by product of two index operator and one operattor on single site.
  //This is just for compatibility of transition pdm.
  //abort(); 
  return;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<RI3index>::renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo* s)
{
  // This is only a skeleton class, so actual operators should never be built
  return;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<RI3index>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

//  // Initialize empty vector
//  for (int i = 0; i < m_op.local_nnz(); ++i) {
//    std::vector<boost::shared_ptr<RI3index> >& spin_ops = m_op.get_local_element(i);
//    spin_ops.clear();
//  }
}

//===========================================================================================================================================================
// (Des,Des,Des) 
//-------------------
// We build this as an alternative to using Tr(CCC)

template<> 
string Op_component<DesDesDes>::get_op_string() const {
  return "DesDesDes";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesDesDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( DES_DES_DES, b, DES, DES_DES, DES, DES_DES, rotateMatrix, stateinfo );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesDesDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{ Three_index_ops::build_3index_ops( DES_DES_DES, b, DES, DES_DES, DES, DES_DES, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<DesDesDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of DDD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<DesDesDes> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > dd_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > d_dd_quantum_ladder;

    // Create (DD)D structure
    //----------------------------------------
    // (DD)
    std::vector<SpinQuantum> spinvec12 = -spin1 - spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (DD)D
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] - spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        dd_d_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -3 );
      }
    }
    assert( dd_d_quantum_ladder.size() == 3 );

    // Create D(DD) structure
    //----------------------------------------
    // (DD)
    std::vector<SpinQuantum> spinvec23 = -spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // D(DD)
      std::vector<SpinQuantum> spinvec123 = -spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        d_dd_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -3 );
      }
    }
    assert( d_dd_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < dd_d_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<DesDesDes>(new DesDesDes) );
      boost::shared_ptr<DesDesDes> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((DD)(D))"] = dd_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["((D)(DD))"] = d_dd_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Cre,Cre)
//-------------------

template<> 
string Op_component<CreCreCre>::get_op_string() const {
  return "CreCreCre";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreCreCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( CRE_CRE_CRE, b, CRE, CRE_CRE, CRE, CRE_CRE, rotateMatrix, stateinfo );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreCreCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( CRE_CRE_CRE, b, CRE, CRE_CRE, CRE, CRE_CRE, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<CreCreCre>::build_iterators(SpinBlock& b)
{
//FIXME don't build CCC operators where all indices are the same!! (May screw up load-balancing, but not much...)
//
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
//  const double screen_tol = dmrginp.screen_tol();
//  assert( dmrginp.screen_tol() == 0 ); //FIXME otherwise some 1 or 2-index ops will be absent when we try to build the 3-index guy
//  std::vector< std::tuple<int,int,int> > tuples = screened_ccd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
//pout << "CCC indices\n";
//for (auto it = tuples.begin(); it != tuples.end(); ++it) {
//pout << "p" << mpigetrank() << ": " << std::get<0>(it->first) << "," << std::get<1>(it->first) << "," << std::get<2>(it->first) << " ; mode = " << it->second << endl;
//}
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
//FIXME remove this
//pout << "New set of CCC operators: p" << mpigetrank() << "; local size = " << m_op.local_nnz() << "; global size = " << m_op.global_nnz() << "; is local " << m_op.is_local() << std::endl;
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "p" << mpigetrank() << "; Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreCreCre> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    // Quantum ladder components give spin of each sub-contraction
    // e.g. cc_c_quantum_ladder[0] is spin of 2-index op
    //      cc_c_quantum_ladder[1] is spin of 3-index op 
    std::vector< std::vector<SpinQuantum> > cc_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_cc_quantum_ladder;

    // Create (CC)C structure
    //----------------------------------------
    // (C1C2)
    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (C1C2)C3
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] + spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        cc_c_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 3 );
      }
    }
    assert( cc_c_quantum_ladder.size() == 3 );

    // Create C(CC) structure
    //----------------------------------------
    // (C2C3)
    std::vector<SpinQuantum> spinvec23 = spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // C1(C2C3)
      std::vector<SpinQuantum> spinvec123 = spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        c_cc_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 3 );
      }
    }
    assert( c_cc_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < cc_c_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<CreCreCre>(new CreCreCre) );
      boost::shared_ptr<CreCreCre> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((CC)(C))"] = cc_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C)(CC))"] = c_cc_quantum_ladder.at(q);
      // This is updated when the build_pattern changes
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Cre,Des)
//-------------------

template<> 
string Op_component<CreCreDes>::get_op_string() const {
  return "CreCreDes";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreCreDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( CRE_CRE_DES, b, CRE, CRE_CRE, DES, CRE_DES, rotateMatrix, stateinfo );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreCreDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( CRE_CRE_DES, b, CRE, CRE_CRE, DES, CRE_DES, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<CreCreDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
//pout << "New set of CCD operators: p" << mpigetrank() << "; local size = " << m_op.local_nnz() << "; global size = " << m_op.global_nnz() << "; is local " << m_op.is_local() << std::endl;
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CCD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreCreDes> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > cc_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_cd_quantum_ladder;

    // Create (CC)D structure
    //----------------------------------------
    // (CC)
    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (CC)D
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] - spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        cc_d_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 1 );
      }
    }
    assert( cc_d_quantum_ladder.size() == 3 );

    // Create C(CD) structure
    //----------------------------------------
    // (CD)
    std::vector<SpinQuantum> spinvec23 = spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // C(CD)
      std::vector<SpinQuantum> spinvec123 = spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        c_cd_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 1 );
      }
    }
    assert( c_cd_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < cc_d_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<CreCreDes>(new CreCreDes) );
      boost::shared_ptr<CreCreDes> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((CC)(D))"] = cc_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C)(CD))"] = c_cd_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Des,Des)
//-------------------

template<> 
string Op_component<CreDesDes>::get_op_string() const {
  return "CreDesDes";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreDesDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( CRE_DES_DES, b, CRE, CRE_DES, DES, DES_DES, rotateMatrix, stateinfo );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreDesDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( CRE_DES_DES, b, CRE, CRE_DES, DES, DES_DES, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<CreDesDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CDD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreDesDes> >& spin_ops = m_op.get_local_element(i);
//pout << "before size() = " << m_op.get_local_element(i).size() << std::endl;

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > cd_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_dd_quantum_ladder;

    // Create (CD)D structure
    //----------------------------------------
    // (CD)
    std::vector<SpinQuantum> spinvec12 = spin1 - spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (CD)D
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] - spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        cd_d_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( cd_d_quantum_ladder.size() == 3 );

    // Create C(DD) structure
    //----------------------------------------
    // (DD)
    std::vector<SpinQuantum> spinvec23 = -spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // C(DD)
      std::vector<SpinQuantum> spinvec123 = spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        c_dd_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( c_dd_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < cd_d_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<CreDesDes>(new CreDesDes) );
      boost::shared_ptr<CreDesDes> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((CD)(D))"] = cd_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C)(DD))"] = c_dd_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Des,Cre)
//-------------------

template<> 
string Op_component<CreDesCre>::get_op_string() const {
  return "CreDesCre";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreDesCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( CRE_DES_CRE, b, CRE, CRE_DES, CRE, DES_CRE, rotateMatrix, stateinfo );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<CreDesCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( CRE_DES_CRE, b, CRE, CRE_DES, CRE, DES_CRE, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<CreDesCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CDC operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreDesCre> >& spin_ops = m_op.get_local_element(i);
//pout << "before size() = " << m_op.get_local_element(i).size() << std::endl;

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > cd_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_dc_quantum_ladder;

    // Create (CD)C structure
    //----------------------------------------
    // (CD)
    std::vector<SpinQuantum> spinvec12 = spin1 - spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (CD)C
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] + spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        cd_c_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 1 );
      }
    }
    assert( cd_c_quantum_ladder.size() == 3 );

    // Create C(DC) structure
    //----------------------------------------
    // (DC)
    std::vector<SpinQuantum> spinvec23 = -spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // C(DC)
      std::vector<SpinQuantum> spinvec123 = spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        c_dc_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 1 );
      }
    }
    assert( c_dc_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < cd_c_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<CreDesCre>(new CreDesCre) );
      boost::shared_ptr<CreDesCre> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((CD)(C))"] = cd_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["((C)(DC))"] = c_dc_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// 4PDM operators
//===========================================================================================================================================================

//===========================================================================================================================================================
// (Des,Cre,Des)
//-------------------

template<> 
string Op_component<DesCreDes>::get_op_string() const {
  return "DesCreDes";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesCreDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( DES_CRE_DES, b, DES, DES_CRE, DES, CRE_DES, rotateMatrix, stateinfo );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesCreDes>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( DES_CRE_DES, b, DES, DES_CRE, DES, CRE_DES, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<DesCreDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of DCD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<DesCreDes> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > dc_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > d_cd_quantum_ladder;

    // Create (DC)D structure
    //----------------------------------------
    // (DC)
    std::vector<SpinQuantum> spinvec12 = -spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (DC)D
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] - spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        dc_d_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( dc_d_quantum_ladder.size() == 3 );

    // Create D(CD) structure
    //----------------------------------------
    // (CD)
    std::vector<SpinQuantum> spinvec23 = spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // D(CD)
      std::vector<SpinQuantum> spinvec123 = - spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        d_cd_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( d_cd_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < dc_d_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<DesCreDes>(new DesCreDes) );
      boost::shared_ptr<DesCreDes> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((DC)(D))"] = dc_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["((D)(CD))"] = d_cd_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Des,Des,Cre)
//-------------------

template<> 
string Op_component<DesDesCre>::get_op_string() const {
  return "DesDesCre";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesDesCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( DES_DES_CRE, b, DES, DES_DES, CRE, DES_CRE, rotateMatrix, stateinfo  );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesDesCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( DES_DES_CRE, b, DES, DES_DES, CRE, DES_CRE, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<DesDesCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of DDC operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<DesDesCre> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > dd_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > d_dc_quantum_ladder;

    // Create (DD)C structure
    //----------------------------------------
    // (DD)
    std::vector<SpinQuantum> spinvec12 = -spin1 - spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (DD)C
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] + spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        dd_c_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( dd_c_quantum_ladder.size() == 3 );

    // Create D(DC) structure
    //----------------------------------------
    // (DC)
    std::vector<SpinQuantum> spinvec23 = -spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // D(DC)
      std::vector<SpinQuantum> spinvec123 = -spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        d_dc_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( d_dc_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < dd_c_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<DesDesCre>(new DesDesCre) );
      boost::shared_ptr<DesDesCre> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((DD)(C))"] = dd_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["((D)(DC))"] = d_dc_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Des,Cre,Cre)
//-------------------

template<> 
string Op_component<DesCreCre>::get_op_string() const {
  return "DesCreCre";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesCreCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) 
{
  Three_index_ops::build_3index_ops( DES_CRE_CRE, b, DES, DES_CRE, CRE, CRE_CRE, rotateMatrix, stateinfo  );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<>
void Op_component<DesCreCre>::build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, 
                                                                                              const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  Three_index_ops::build_3index_ops( DES_CRE_CRE, b, DES, DES_CRE, CRE, CRE_CRE, leftMat, bra, rightMat, ket);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<DesCreCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  std::map< std::tuple<int,int,int>, int > tuples = get_3index_tuples(b);
  m_op.set_tuple_indices( tuples, dmrginp.last_site() );

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of DCC operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<DesCreCre> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
    SpinQuantum spin2 = getSpinQuantum(orbs[1]);
    SpinQuantum spin3 = getSpinQuantum(orbs[2]);

    std::vector< std::vector<SpinQuantum> > dc_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > d_cc_quantum_ladder;

    // Create (DC)C structure
    //----------------------------------------
    // (DC)
    std::vector<SpinQuantum> spinvec12 = -spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      // (DC)C
      std::vector<SpinQuantum> spinvec123 = spinvec12[p] + spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        dc_c_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 1 );
      }
    }
    assert( dc_c_quantum_ladder.size() == 3 );

    // Create D(CC) structure
    //----------------------------------------
    // (CC)
    std::vector<SpinQuantum> spinvec23 = spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      // D(CC)
      std::vector<SpinQuantum> spinvec123 = -spin1 + spinvec23[p];
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec23[p], spinvec123[q] };
        d_cc_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == 1 );
      }
    }
    assert( d_cc_quantum_ladder.size() == 3 );

    // Allocate new operator for each spin component
    //------------------------------------------------
    spin_ops.clear();
    for (int q=0; q < dc_c_quantum_ladder.size(); q++) {
      spin_ops.push_back( boost::shared_ptr<DesCreCre>(new DesCreCre) );
      boost::shared_ptr<DesCreCre> op = spin_ops.back();
      op->set_orbs() = orbs;
      op->set_initialised() = true;
      op->set_quantum_ladder()["((DC)(C))"] = dc_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["((D)(CC))"] = d_cc_quantum_ladder.at(q);
      op->set_deltaQuantum(1, op->get_quantum_ladder().at( op->get_build_pattern() ).at(1) );
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================

}
