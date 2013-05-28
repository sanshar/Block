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
  // Set up site indices such that (i <= j <= k)  (Note Block mostly sets up with j<=i, but that's not very convenient for us here).
  std::vector<std::tuple<int,int,int> > screened_indices;
//FIXME  only i=j=k at the moment!!
  for (int k = 0; k < sites.size(); ++k) {
    for (int j = k; j <= k; ++j) {
      for (int i = j; i <= j; ++i) {
//      if (dmrginp.use_partial_two_integrals()) {
        screened_indices.push_back(std::make_tuple(sites[i], sites[j], sites[k]));
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
pout << "New set of 3-index operators:  " << i << std::endl;
pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreCreDes> >& spin_ops = m_op.get_local_element(i);
pout << "before size() = " << m_op.get_local_element(i).size() << std::endl;
    spin_ops.clear();

    // Combine first 2 indices first
    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    // Note Cre_Cre => plus sign
    std::vector<SpinQuantum> spinvec12 = spin1+spin2;

    // Now combine with 3rd index
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));
pout << "spins = " << spin1.get_s()/2.0 << " " << spin2.get_s()/2.0 << " " << spin3.get_s()/2.0 << std::endl;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      // Note X_Des => minus sign
      std::vector<SpinQuantum> spinvec123 = spinvec12[p]-spin3;

      // Allocate new operator for each spin component
      for (int q=0; q < spinvec123.size(); q++) {
        spin_ops.push_back( boost::shared_ptr<CreCreDes>(new CreCreDes) );
        SparseMatrix& op = *spin_ops[spin_ops.size()-1];
        op.set_orbs() = orbs;
        op.set_initialised() = true;
//FIXME!!! why false or true??
        op.set_fermion() = true;
        op.set_deltaQuantum() = spinvec123[q];      
        op.set_quantum_ladder() = { spinvec12[p], spinvec123[q] };
        assert( spinvec12[p].particleNumber  == 2 );
        assert( spinvec123[p].particleNumber == 1 );
pout << "3-index operators spin composition:\n";
pout << spinvec12[p].get_s()/2.0 << "  " << spinvec123[q].get_s()/2.0 << std::endl;
      }
    }
pout << "after size() = " << m_op.get_local_element(i).size() << std::endl;
  }

}
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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
pout << "New set of 3-index operators:  " << i << std::endl;
pout << "Orbs = " << orbs[0] << " " << orbs[1] << " " << orbs[2] << std::endl;
    std::vector<boost::shared_ptr<CreDesDes> >& spin_ops = m_op.get_local_element(i);
pout << "before size() = " << m_op.get_local_element(i).size() << std::endl;
    spin_ops.clear();

    // Combine first 2 indices first
    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
    // Note Cre_Des => minus sign
    std::vector<SpinQuantum> spinvec12 = spin1-spin2;

    // Now combine with 3rd index
    SpinQuantum spin3 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[2]));
pout << "spins = " << spin1.get_s()/2.0 << " " << spin2.get_s()/2.0 << " " << spin3.get_s()/2.0 << std::endl;
    for (int p=0; p < spinvec12.size(); p++) {
      assert( spinvec12[p].get_s() == 2*p );
      // Note X_Des => minus sign
      std::vector<SpinQuantum> spinvec123 = spinvec12[p]-spin3;

      // Allocate new operator for each spin component
      for (int q=0; q < spinvec123.size(); q++) {
        spin_ops.push_back( boost::shared_ptr<CreDesDes>(new CreDesDes) );
        SparseMatrix& op = *spin_ops[spin_ops.size()-1];
        op.set_orbs() = orbs;
        op.set_initialised() = true;
//FIXME!!! why false or true??
        op.set_fermion() = true;
        op.set_deltaQuantum() = spinvec123[q];      
        op.set_quantum_ladder() = { spinvec12[p], spinvec123[q] };
        assert( spinvec12[p].particleNumber  == 2 );
        assert( spinvec123[p].particleNumber == 1 );
pout << "3-index operators spin composition:\n";
pout << spinvec12[p].get_s()/2.0 << "  " << spinvec123[q].get_s()/2.0 << std::endl;
      }
    }
pout << "after size() = " << m_op.get_local_element(i).size() << std::endl;
  }

}
  
//===========================================================================================================================================================

}

