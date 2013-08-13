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

std::vector<std::tuple<int,int,int> > get_local_3index_tuples(SpinBlock& b)
{
  std::vector<std::tuple<int,int,int> > tuples;

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();
  
  assert( dotBlock != NULL );
  assert( dotBlock->get_sites().size() == 1 );

  // Forward or backwards sweep 
  bool forward = true;
  if ( sysBlock->get_sites()[0] > dotBlock->get_sites()[0] ) forward = false;
  cout << "MAW CCC iterators forward? p" << mpigetrank() << " = " << forward << std::endl;
  int dot = dotBlock->get_sites()[0];
  pout << "dot = " << dot << std::endl;

  // 3 on dot
  tuples.push_back( std::make_tuple(dot, dot, dot) );
  // 2 on dot
  for (auto i = sysBlock->get_sites().begin(); i != sysBlock->get_sites().end(); ++i) {
     if ( forward )
       tuples.push_back( std::make_tuple(dot, dot, *i) );
     else
       tuples.push_back( std::make_tuple(*i, dot, dot) );
  }
  // 1 on dot
//FIXME check that these are local!
  std::vector< std::vector<int> > ij_array = sysBlock->get_op_array(CRE_CRE).get_array();
  for (auto ij = ij_array.begin(); ij != ij_array.end(); ++ij) {
     pout << "ij = " << (*ij)[0] << "," << (*ij)[1] << std::endl;
     assert( (*ij)[0] >= (*ij)[1] );
     if ( forward )
       tuples.push_back( std::make_tuple(dot, (*ij)[0], (*ij)[1]) );
     else
       tuples.push_back( std::make_tuple((*ij)[0], (*ij)[1], dot) );
  }
  std::vector< std::vector<int> > ijk_array = sysBlock->get_op_array(CRE_CRE_CRE).get_array();
  // 0 on dot
  for (auto ijk = ijk_array.begin(); ijk != ijk_array.end(); ++ijk) {
     assert( (*ijk)[0] >= (*ijk)[1] );
     assert( (*ijk)[1] >= (*ijk)[2] );
     tuples.push_back( std::make_tuple((*ijk)[0], (*ijk)[1], (*ijk)[2]) );
  }

  return tuples;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
// Choose 3-index tuples on this MPI process such that 2-index are available to build them

//FIXME screening?
std::vector<std::tuple<int,int,int> > get_3index_tuples(SpinBlock& b, bool& global_is_local)
{
  std::vector<std::tuple<int,int,int> > tuples;

  if ( b.get_leftBlock() != NULL ) {
    // Generate only mpi local tuples for compound block, consistent with existing operators on sys and dot
    global_is_local = true;
    tuples = get_local_3index_tuples(b);    
  }
  else {
    // Generate all tuples such that (k <= j <= i) and let para_array assign them to local processes as necessary
    global_is_local = false;
    std::vector<int> sites = b.get_sites();
    for (int i = 0; i < sites.size(); ++i)
      for (int j = 0; j <= i; ++j)
        for (int k = 0; k <= j; ++k)
          tuples.push_back(std::make_tuple(sites[i], sites[j], sites[k]));
  }

  return tuples;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
//
//bool is_forward(SpinBlock& b) 
//{
//  SpinBlock* sysBlock = b.get_leftBlock();
//  SpinBlock* dotBlock = b.get_rightBlock();
//  
//  assert( dotBlock != NULL );
//  assert( dotBlock->get_sites().size() == 1 );
//  bool forward = true;
//  if ( sysBlock->get_sites()[0] > dotBlock->get_sites()[0] ) forward = false;
//  cout << "MAW CCC iterators forward? p" << mpigetrank() << " = " << forward << std::endl;
//
//  return forward;
//}
//
//===========================================================================================================================================================

//FIXME the 3-index routines below have a lot of common ground that could be re-implemented more clearly

//===========================================================================================================================================================
// (Cre,Cre,Cre)
//-------------------

template<> 
string Op_component<CreCreCre>::get_op_string() const {
  return "CRECRECRE";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<CreCreCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
//  const double screen_tol = dmrginp.screen_tol();
//  std::vector< std::tuple<int,int,int> > tuples = screened_ccd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
  bool global_is_local = false;
  std::vector< std::tuple<int,int,int> > tuples = get_3index_tuples(b, global_is_local);
  m_op.set_tuple_indices(tuples, dmrginp.last_site(), global_is_local);      

  // Allocate new set of operators for each set of spatial orbitals
//FIXME check that we have load balancing
cout << "New set of CCC operators: p" << mpigetrank() << "; size = " << m_op.local_nnz() << std::endl;
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreCreCre> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    // Quantum ladder components give spin of each sub-contraction
    // e.g. cc_c_quantum_ladder[0] is spin of cc 2-index op
    //      cc_c_quantum_ladder[1] is spin of c 1-index op
    std::vector< std::vector<SpinQuantum> > cc_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_cc_quantum_ladder;

    // Create (CC)C structure
    //----------------------------------------
    // Cre_Cre => plus sign
    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      assert( spinvec12[p].particleNumber == 2 );
      // X_Cre => plus sign
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
    // Cre_Cre => plus sign
    std::vector<SpinQuantum> spinvec23 = spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      assert( spinvec23[p].get_s() == 2*p );
      assert( spinvec23[p].particleNumber == 2 );
      // X_(CC) => plus sign
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
      // 3-index ops are fermionic 
      op->set_fermion() = true;
      op->set_quantum_ladder()["((CC)C)"] = cc_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["(C(CC))"] = c_cc_quantum_ladder.at(q);
//FIXME  This should be updated when the build_pattern changes
//FIXME  Set default
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Cre,Des)
//-------------------

template<> 
string Op_component<CreCreDes>::get_op_string() const {
  return "CRECREDES";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<CreCreDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  bool global_is_local = false;
  std::vector< std::tuple<int,int,int> > tuples = get_3index_tuples(b, global_is_local);
  m_op.set_tuple_indices(tuples, dmrginp.last_site(), global_is_local);      

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CCD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreCreDes> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    std::vector< std::vector<SpinQuantum> > cc_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_cd_quantum_ladder;

    // Create (CC)D structure
    //----------------------------------------
    // Cre_Cre => plus sign
    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      assert( spinvec12[p].particleNumber == 2 );
      // X_Des => minus sign
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
    // Cre_Des => minus sign
    std::vector<SpinQuantum> spinvec23 = spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      assert( spinvec23[p].get_s() == 2*p );
      assert( spinvec23[p].particleNumber == 0 );
      //FIXME minus or plus sign???
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
      // 3-index ops are fermionic 
      op->set_fermion() = true;
      op->set_quantum_ladder()["((CC)D)"] = cc_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["(C(CD))"] = c_cd_quantum_ladder.at(q);
//FIXME  This should be updated when the build_pattern changes
//FIXME  Set default
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Des,Des)
//-------------------

template<> 
string Op_component<CreDesDes>::get_op_string() const {
  return "CREDESDES";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<CreDesDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  bool global_is_local = false;
  std::vector< std::tuple<int,int,int> > tuples = get_3index_tuples(b, global_is_local);
  m_op.set_tuple_indices(tuples, dmrginp.last_site(), global_is_local);      

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CDD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreDesDes> >& spin_ops = m_op.get_local_element(i);
//pout << "before size() = " << m_op.get_local_element(i).size() << std::endl;

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    std::vector< std::vector<SpinQuantum> > cd_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_dd_quantum_ladder;

    // Create (CD)D structure
    //----------------------------------------
    // Cre_Des => minus sign
    std::vector<SpinQuantum> spinvec12 = spin1 - spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      assert( spinvec12[p].particleNumber == 0 );
      // X_Des => minus sign
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
    // Des_Des => Cre_Cre => plus sign and transpose and commute (factor -1??)
    std::vector<SpinQuantum> spinvec23 = spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      assert( spinvec23[p].get_s() == 2*p );
//FIXME
      assert( spinvec23[p].particleNumber == 2 );
//FIXME minus or plus sign???  FACTOR -1 ??
      std::vector<SpinQuantum> spinvec123 = spin1 - spinvec23[p];
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
      // 3-index ops are fermionic 
      op->set_fermion() = true;
      op->set_quantum_ladder()["((CD)D)"] = cd_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["(C(DD))"] = c_dd_quantum_ladder.at(q);
//FIXME  This should be updated when the build_pattern changes
//FIXME  Set default
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================
// (Cre,Des,Cre)
//-------------------

template<> 
string Op_component<CreDesCre>::get_op_string() const {
  return "CREDESCRE";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

template<> 
void Op_component<CreDesCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  bool global_is_local = false;
  std::vector< std::tuple<int,int,int> > tuples = get_3index_tuples(b, global_is_local);
  m_op.set_tuple_indices(tuples, dmrginp.last_site(), global_is_local);      

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CDC operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreDesCre> >& spin_ops = m_op.get_local_element(i);
//pout << "before size() = " << m_op.get_local_element(i).size() << std::endl;

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    std::vector< std::vector<SpinQuantum> > cd_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_dc_quantum_ladder;

    // Create (CD)C structure
    //----------------------------------------
    // Cre_Des => minus sign
    std::vector<SpinQuantum> spinvec12 = spin1 - spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      assert( spinvec12[p].particleNumber == 0 );
      // X_Cre => plus sign
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
    // Des_Cre => Cre_Des => minus sign and transpose and commute (factor -1??)
//FIXME how to do DC??
    std::vector<SpinQuantum> spinvec23 = -spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      assert( spinvec23[p].get_s() == 2*p );
//FIXME
      assert( spinvec23[p].particleNumber == 0 );
//FIXME minus or plus sign???  FACTOR -1 ??
      std::vector<SpinQuantum> spinvec123 = spin1 - spinvec23[p];
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
      // 3-index ops are fermionic 
      op->set_fermion() = true;
      op->set_quantum_ladder()["((CD)C)"] = cd_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["(C(DC))"] = c_dc_quantum_ladder.at(q);
//FIXME  This should be updated when the build_pattern changes
//FIXME  Set default
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}
  
//===========================================================================================================================================================
// (Des,Cre,Des)
//-------------------

template<> 
string Op_component<DesCreDes>::get_op_string() const {
  return "DESCREDES";
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
template<> 
void Op_component<DesCreDes>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  bool global_is_local = false;
  std::vector< std::tuple<int,int,int> > tuples = get_3index_tuples(b, global_is_local);
  m_op.set_tuple_indices(tuples, dmrginp.last_site(), global_is_local);      

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of DCD operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<DesCreDes> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    std::vector< std::vector<SpinQuantum> > dc_d_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > d_cd_quantum_ladder;

    // Create (DC)D structure
    //----------------------------------------
    // Des_Cre => minus, plus
    std::vector<SpinQuantum> spinvec12 = -spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      assert( spinvec12[p].particleNumber == 0 );
      // X_Des => minus sign
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
    // Cre_Des => plus sign
    std::vector<SpinQuantum> spinvec23 = spin2 - spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      assert( spinvec23[p].get_s() == 2*p );
      assert( spinvec23[p].particleNumber == 0 );
      // D_(CD) => minus, plus sign //FIXME ???
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
      // 3-index ops are fermionic 
      op->set_fermion() = true;
      op->set_quantum_ladder()["((DC)D)"] = dc_d_quantum_ladder.at(q);
      op->set_quantum_ladder()["(D(CD))"] = d_cd_quantum_ladder.at(q);
//FIXME  This should be updated when the build_pattern changes
//FIXME  Set default
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
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
void Op_component<DesDesCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  // Set up 3-index (i,j,k) spatial operator indices for this SpinBlock
  bool global_is_local = false;
  std::vector< std::tuple<int,int,int> > tuples = get_3index_tuples(b, global_is_local);
  m_op.set_tuple_indices(tuples, dmrginp.last_site(), global_is_local);      

  // Allocate new set of operators for each set of spatial orbitals
  std::vector<int> orbs(3);
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of DDC operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<DesDesCre> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    std::vector< std::vector<SpinQuantum> > dd_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > d_dc_quantum_ladder;

    // Create (DD)C structure
    //----------------------------------------
    // DD => CC => minus sign and transpose and commute (factor -1??)
    std::vector<SpinQuantum> spinvec12 = spin1 + spin2;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      assert( spinvec12[p].particleNumber == 2 );
      // DD_Cre => minus, plus sign ?? FIXME
      std::vector<SpinQuantum> spinvec123 = -spinvec12[p] + spin3;
      for (int q=0; q < spinvec123.size(); q++) {
        std::vector<SpinQuantum> tmp = { spinvec12[p], spinvec123[q] };
        dd_c_quantum_ladder.push_back( tmp );
        assert( spinvec123[q].particleNumber == -1 );
      }
    }
    assert( dd_c_quantum_ladder.size() == 3 );

    // Create D(DC) structure
    //----------------------------------------
    // DC => minus, plus sign
    std::vector<SpinQuantum> spinvec23 = -spin2 + spin3;
    for (int p=0; p < spinvec23.size(); p++) {
      assert( spinvec23[p].get_s() == 2*p );
      assert( spinvec23[p].particleNumber == 0 );
      // D_(DC) => minus, minus sign ?? //FIXME
      std::vector<SpinQuantum> spinvec123 = -spin1 - spinvec23[p];
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
      // 3-index ops are fermionic 
      op->set_fermion() = true;
      op->set_quantum_ladder()["((DD)C)"] = dd_c_quantum_ladder.at(q);
      op->set_quantum_ladder()["(D(DC))"] = d_dc_quantum_ladder.at(q);
//FIXME  This should be updated when the build_pattern changes
//FIXME  Set default
      op->set_deltaQuantum() = op->get_quantum_ladder().at( op->get_build_pattern() ).at(1);
    }

    assert( m_op.get_local_element(i).size() == 3);
  }
}

//===========================================================================================================================================================

}
