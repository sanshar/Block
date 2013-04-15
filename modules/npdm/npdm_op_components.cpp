/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
//
//  This is an extension of op_components.C
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

#include "spinblock.h"
#include "op_components.h"
//#include "screen.h"

namespace SpinAdapted {
  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
//FIXME update SCREEN.C

std::vector<std::tuple<int,int,int> > screened_ccc_indices(const vector<int, std::allocator<int> >& indices,
                                                           const vector<int, std::allocator<int> >& interactingix,
                                                           const TwoElectronArray& twoe, double thresh)
{
  std::vector<std::tuple<int,int,int> > screened_indices;
//FIXME
//  for (int i = 0; i < indices.size(); ++i)
//    for (int j = 0; j <= i; ++j) {
//      if (dmrginp.use_partial_two_integrals()) {
//        screened_indices.push_back(make_pair(indices[i], indices[j]));
//      }
//      else {
//        if (screen_ccc_interaction(indices[i], indices[j], interactingix, twoe, thresh))
//          screened_indices.push_back(make_pair(indices[i], indices[j]));
//      }
//    }
  return screened_indices;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------  
// (Cre,Cre,Cre)
//-------------------

template<> 
string Op_component<CreCreCre>::get_op_string() const {
  return "CRECRECRE";
}

template<> 
void Op_component<CreCreCre>::build_iterators(SpinBlock& b)
{
  // Blank construction (used in unset_initialised() Block copy construction, for use with STL)
  if (b.get_sites().size () == 0) return; 

  const double screen_tol = dmrginp.screen_tol();
//FIXME
  std::vector< std::tuple<int,int,int> > screened_ccc_ix = screened_ccc_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
//FIXME
//  m_op.set_tuple_indices(screened_ccc_ix, dmrginp.last_site());      
  std::vector<int> orbs(3);

  for (int i = 0; i < m_op.local_nnz(); ++i) {
//FIXME    std::tuple<int,int,int> tuple = m_op.unmap_local_index(i);
//    orbs[0] = std::get<0>(tuple);
//    orbs[1] = std::get<1>(tuple);
//    orbs[2] = std::get<2>(tuple);
  
    std::vector<boost::shared_ptr<CreCreCre> >& vec = m_op.get_local_element(i);
    SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
    SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
  
    std::vector<SpinQuantum> spinvec = spin1-spin2;
  
    vec.resize(spinvec.size());
    for (int j=0; j<spinvec.size(); j++) {
      vec[j]=boost::shared_ptr<CreCreCre>(new CreCreCre);
      SparseMatrix& op = *vec[j];
      op.set_orbs() = orbs;
      op.set_initialised() = true;
      op.set_fermion() = false;
      op.set_deltaQuantum() = spinvec[j];      
    }
  
  }
}
  
//-------------------------------------------------------------------------------------------------------------------------------------------------------------  

}

