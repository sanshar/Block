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

#include "spinblock.h"
#include "op_components.h"
//#include "screen.h"

namespace SpinAdapted {
  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
//FIXME update SCREEN.C

std::vector<std::tuple<int,int,int> > screened_ccd_indices(const vector<int, std::allocator<int> >& sites,
                                                           const vector<int, std::allocator<int> >& interactingix,
                                                           const TwoElectronArray& twoe, double thresh)
{
  // Set up site indices such that (k <= j <= i)
  std::vector<std::tuple<int,int,int> > screened_indices;
//FIXME  only i=j=k at the moment!!
//  for (int i = 0; i < sites.size(); ++i) {
//    int j =i ;
//    int k =i ;
//    {
//    {
//FIXME should this stride order match para_array??
  for (int i = 0; i < sites.size(); ++i) {
    for (int j = 0; j <= i; ++j) {
      for (int k = 0; k <= j; ++k) {
//      if (dmrginp.use_partial_two_integrals()) {
        screened_indices.push_back(std::make_tuple(sites[i], sites[j], sites[k]));
//pout << i << " " << j << " " << k << std::endl;

//      }
//FIXME
//      else {
//        if (screen_ccd_interaction(indices[i], indices[j], interactingix, twoe, thresh))
//          screened_indices.push_back(make_pair(indices[i], indices[j]));
//      }
      }
    }
  }
  return screened_indices;
}

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
  const double screen_tol = dmrginp.screen_tol();
//FIXME
  std::vector< std::tuple<int,int,int> > screened_ccc_ix = screened_ccd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
  m_op.set_tuple_indices(screened_ccc_ix, dmrginp.last_site());      
  std::vector<int> orbs(3);

  // Allocate new set of operators for each set of spatial orbitals (For 3-index this part cannot be stored all in core, even locally ??)
  for (int i = 0; i < m_op.local_nnz(); ++i) {
    orbs = m_op.unmap_local_index(i);
//pout << "New set of CCC operators:  " << i << std::endl;
//pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreCreCre> >& spin_ops = m_op.get_local_element(i);

    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));

    std::vector< std::vector<SpinQuantum> > cc_c_quantum_ladder;
    std::vector< std::vector<SpinQuantum> > c_cc_quantum_ladder;

    // Create (CC)C structure first
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

    // Create C(CC) structure next
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
  const double screen_tol = dmrginp.screen_tol();
//FIXME
  std::vector< std::tuple<int,int,int> > screened_ccd_ix = screened_ccd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
  m_op.set_tuple_indices(screened_ccd_ix, dmrginp.last_site());      
  std::vector<int> orbs(3);

  // Allocate new set of operators for each set of spatial orbitals (For 3-index this part cannot be stored all in core, even locally ??)
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

    // Create (CC)D structure first
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

    // Create C(CD) structure next
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
  const double screen_tol = dmrginp.screen_tol();
//FIXME
  std::vector< std::tuple<int,int,int> > screened_cdd_ix = screened_ccd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
  m_op.set_tuple_indices(screened_cdd_ix, dmrginp.last_site());      
  std::vector<int> orbs(3);

  // Allocate new set of operators for each set of spatial orbitals (For 3-index this part cannot be stored all in core, even locally ??)
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

    // Create (CD)D structure first
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

    // Create C(DD) structure next
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
  const double screen_tol = dmrginp.screen_tol();
//FIXME
  std::vector< std::tuple<int,int,int> > screened_cdc_ix = screened_ccd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
  m_op.set_tuple_indices(screened_cdc_ix, dmrginp.last_site());      
  std::vector<int> orbs(3);

  // Allocate new set of operators for each set of spatial orbitals (For 3-index this part cannot be stored all in core, even locally ??)
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

    // Create (CD)C structure first
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

    // Create C(DC) structure next
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

}

