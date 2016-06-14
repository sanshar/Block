/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "spinblock.h"
#include "op_components.h"
#include "screen.h"

namespace SpinAdapted {
 
  // -------------------- C_S1 ---------------------------  
  template<> string Op_component<Cre>::get_op_string() const {
    return "CRE";
  }

  template<> void Op_component<Cre>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      std::vector<int> screened_ix;
      if(dmrginp.calc_type() == MPS_NEVPT)
      {
        if(b.nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2))
          screened_ix =  screened_cdd_c_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
        else
          screened_ix =  screened_ccd_c_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);

      }
      else{

        int integralIndex = b.get_integralIndex();
        screened_ix = (dmrginp.hamiltonian() == BCS) ? 
          screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) : 
          screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      }
      m_op.set_indices(screened_ix, dmrginp.last_site());  
      std::vector<int> orbs(1);
      
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<Cre>(new Cre);
	  SparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  op.set_deltaQuantum(1, getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
	  //op.set_deltaQuantum() = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOf(orbs[0]));      
     op.set_quantum_ladder()["(C)"] = { op.get_deltaQuantum(0) };
	}
      
    }
  
  
  
  template<> std::vector< std::vector<int> > Op_component<Cre>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }
  
  template<> void Op_component<Cre>::add_local_indices(int i, int j , int k) 
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<Cre> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<Cre>(new Cre);
    }

  //usually not needed, because it can be calculated as a transpose of C, but
  //when the bra and ket state in the block are different than transpose cannot be used
  // -------------------- D_S1 ---------------------------  
  template<> string Op_component<Des>::get_op_string() const {
    return "DES";
  }

  template<> void Op_component<Des>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      int integralIndex = b.get_integralIndex();
      const double screen_tol = dmrginp.oneindex_screen_tol();
      std::vector<int> screened_ix;
      if(dmrginp.calc_type() == MPS_NEVPT)
      {
        if(b.nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2))
          screened_ix =  screened_cdd_d_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
        else
          screened_ix =  screened_ccd_d_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);

      }
      else{

        screened_ix = (dmrginp.hamiltonian() == BCS) ? 
          screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) : 
          screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      }
      m_op.set_indices(screened_ix, dmrginp.last_site());  
      std::vector<int> orbs(1);
      
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<Des>(new Des);
	  SparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  op.set_deltaQuantum(1, -getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
     op.set_quantum_ladder()["(D)"] = { op.get_deltaQuantum(0) };
	}
      
    }
  
  
  
  template<> std::vector< std::vector<int> > Op_component<Des>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }
  

  template<> void Op_component<Des>::add_local_indices(int i, int j , int k) 
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<Des> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<Des>(new Des);
    }


  
  // -------------------- Cd_ ---------------------------  
  template<> string Op_component<CreDes>::get_op_string() const {
    return "CREDES";
  }

  template<> void Op_component<CreDes>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ? 
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<CreDes> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1-spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<CreDes>(new CreDes);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
       op.set_quantum_ladder()["(CD)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  template<> void Op_component<CreDes>::add_local_indices(int i, int j , int k) {};

  // -------------------- dC_ ---------------------------  
  template<> string Op_component<DesCre>::get_op_string() const {
    return "DESCRE";
  }

  template<> void Op_component<DesCre>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ? 
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<DesCre> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<DesCre>(new DesCre);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
       op.set_quantum_ladder()["(DC)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  template<> void Op_component<DesCre>::add_local_indices(int i, int j , int k) {};
  
  
  // -------------------- Cc_ ---------------------------  
  template<> string Op_component<CreCre>::get_op_string() const {
    return "CRECRE";
  }
  template<> void Op_component<CreCre>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :        
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<CreCre> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<CreCre>(new CreCre);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
       op.set_quantum_ladder()["(CC)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  template<> void Op_component<CreCre>::add_local_indices(int i, int j , int k) {};
  
  // -------------------- Dd_ ---------------------------  
  template<> string Op_component<DesDes>::get_op_string() const {
    return "DESDES";
  }
  template<> void Op_component<DesDes>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :        
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<DesDes> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = -getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = -getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<DesDes>(new DesDes);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
       op.set_quantum_ladder()["(DD)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  template<> void Op_component<DesDes>::add_local_indices(int i, int j , int k) {};
  
  // -------------------- Cdcomp_ ---------------------------  
  template<> string Op_component<CreDesComp>::get_op_string() const {
    return "CREDES_COMP";
  }
  template<> void Op_component<CreDesComp>::build_iterators(SpinBlock& b)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.twoindex_screen_tol();
    vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ?
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
    m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
    
    std::vector<int> orbs(2);
    for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<CreDesComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<CreDesComp>(new CreDesComp);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
        if (dmrginp.hamiltonian() == BCS) {
          op.resize_deltaQuantum(3);
          op.set_deltaQuantum(0) = spinvec[j];
          op.set_deltaQuantum(1) = SpinQuantum(2, spinvec[j].get_s(), spinvec[j].get_symm());
          op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
        } else {
          op.set_deltaQuantum(1, spinvec[j]);
        }
      }
	}
  }

  template<> void Op_component<CreDesComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<CreDesComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin2-spin1;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) 
	vec[j]=boost::shared_ptr<CreDesComp>(new CreDesComp);
    }
  
  // -------------------- dCcomp_ ---------------------------  
  template<> string Op_component<DesCreComp>::get_op_string() const {
    return "DESCRE_COMP";
  }
  template<> void Op_component<DesCreComp>::build_iterators(SpinBlock& b)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.twoindex_screen_tol();
    vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ?
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
    m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
    
    std::vector<int> orbs(2);
    for (int i = 0; i < m_op.local_nnz(); ++i)
    {
	   orbs = m_op.unmap_local_index(i);
      std::vector<boost::shared_ptr<DesCreComp> >& vec = m_op.get_local_element(i);
      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
      SpinQuantum spin2 = getSpinQuantum(orbs[1]);
      std::vector<SpinQuantum> spinvec = spin1-spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<DesCreComp>(new DesCreComp);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
        if (dmrginp.hamiltonian() == BCS) {
          op.resize_deltaQuantum(3);
          op.set_deltaQuantum(0) = spinvec[j];
          op.set_deltaQuantum(1) = SpinQuantum(2, spinvec[j].get_s(), spinvec[j].get_symm());
          op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
        } else {
          op.set_deltaQuantum(1, spinvec[j]);
        }
      }
    }
  }
  
  
  template<> void Op_component<DesCreComp>::add_local_indices(int i, int j , int k) 
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<DesCreComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1-spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) 
	vec[j]=boost::shared_ptr<DesCreComp>(new DesCreComp);
    }
  
  // -------------------- Ddcomp_ ---------------------------  
  template<> string Op_component<DesDesComp>::get_op_string() const {
    return "DESDES_COMP";
  }
  template<> void Op_component<DesDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	   orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<DesDesComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<DesDesComp>(new DesDesComp);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
        
        if (dmrginp.hamiltonian() == BCS) {
          op.resize_deltaQuantum(3);          
          op.set_deltaQuantum(0) = -spinvec[j];
          op.set_deltaQuantum(1) = -SpinQuantum(0, spinvec[j].get_s(), spinvec[j].get_symm());
          op.set_deltaQuantum(2) = -SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
        } else {
	      op.set_deltaQuantum(1, -spinvec[j]);
        }  
	  }
	}
      
    }
  
  template<> void Op_component<DesDesComp>::add_local_indices(int i, int j , int k) 
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<DesDesComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1+spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) 
	    vec[j]=boost::shared_ptr<DesDesComp>(new DesDesComp);
    }
  
  
  // -------------------- CCcomp_ ---------------------------  
  template<> string Op_component<CreCreComp>::get_op_string() const {
    return "CRECRE_COMP";
  }
  template<> void Op_component<CreCreComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<CreCreComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<CreCreComp>(new CreCreComp);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
        if (dmrginp.hamiltonian() == BCS) {
          op.resize_deltaQuantum(3);          
          op.set_deltaQuantum(0) = spinvec[j];
          op.set_deltaQuantum(1) = SpinQuantum(0, spinvec[j].get_s(), spinvec[j].get_symm());
          op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
        } else {
	      op.set_deltaQuantum(1, spinvec[j]);
        }
	  }
	}
      
    }
  
  template<> void Op_component<CreCreComp>::add_local_indices(int i, int j , int k) 
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<CreCreComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1+spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) 
	vec[j]=boost::shared_ptr<CreCreComp>(new CreCreComp);
    }
  
  
  // -------------------- Ccdcomp_ ---------------------------  
  template<> string Op_component<CreCreDesComp>::get_op_string() const {
    return "CRECREDES_COMP";
  }
  template<> void Op_component<CreCreDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      vector< int > screened_cdd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      m_op.set_indices(screened_cdd_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<CreCreDesComp>(new CreCreDesComp);
	  SparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
      if (dmrginp.hamiltonian() == BCS) {
        op.resize_deltaQuantum(4);
        SpinQuantum qorb = getSpinQuantum(orbs[0]);
        op.set_deltaQuantum(0) = qorb;
        op.set_deltaQuantum(1) = SpinQuantum(3, qorb.get_s(), qorb.get_symm());
        op.set_deltaQuantum(2) = SpinQuantum(-1, qorb.get_s(), qorb.get_symm());
        op.set_deltaQuantum(3) = SpinQuantum(-3, qorb.get_s(), qorb.get_symm());
      } else {
	    op.set_deltaQuantum(1, getSpinQuantum(orbs[0]));
      }
	}
    }
  
  
  template<> std::vector<std::vector<int> > Op_component<CreCreDesComp>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }

  template<> void Op_component<CreCreDesComp>::add_local_indices(int i, int j , int k) {};

  //usually not needed, because it can be calculated as a transpose of CCDcomp, but
  //when the bra and ket state in the block are different than transpose cannot be used
  // -------------------- Cddcomp_ ---------------------------  
  template<> string Op_component<CreDesDesComp>::get_op_string() const {
    return "CREDESDES_COMP";
  }
  template<> void Op_component<CreDesDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      vector< int > screened_cdd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc, v_cccc, v_cccd, screen_tol) :
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      m_op.set_indices(screened_cdd_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<CreDesDesComp>(new CreDesDesComp);
	  SparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
      if (dmrginp.hamiltonian() == BCS) {
        op.resize_deltaQuantum(4);
        SpinQuantum qorb = getSpinQuantum(orbs[0]);
        op.set_deltaQuantum(0) = -qorb;
        op.set_deltaQuantum(1) = -SpinQuantum(3, qorb.get_s(), qorb.get_symm());
        op.set_deltaQuantum(2) = -SpinQuantum(-1, qorb.get_s(), qorb.get_symm());
        op.set_deltaQuantum(3) = -SpinQuantum(-3, qorb.get_s(), qorb.get_symm());
      } else {
	    op.set_deltaQuantum(1, -getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
      }     
	}
    }
  
  
  template<> std::vector<std::vector<int> > Op_component<CreDesDesComp>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }

  template<> void Op_component<CreDesDesComp>::add_local_indices(int i, int j , int k) {};
  
  // -------------------- HAM ---------------------------  
  template<> string Op_component<Ham>::get_op_string() const {
    return "HAM";
  }
  template<> void Op_component<Ham>::build_iterators(SpinBlock& b)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<Ham>(new Ham);
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = false;
      if (dmrginp.hamiltonian() == BCS) {
        m_op(0)[0]->resize_deltaQuantum(5);
        for (int i = 0; i <5; ++i) {
          m_op(0)[0]->set_deltaQuantum(i) = SpinQuantum(2*(i-2), SpinSpace(0), IrrepSpace(0) );
        }    
      } else {
        m_op(0)[0]->set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
      }      
    }
  
  template<> std::vector<std::vector<int> > Op_component<Ham>::get_array() const 
    {
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      return ret_val;
    }

  template<> void Op_component<Ham>::add_local_indices(int i, int j , int k) {};
  
  // -------------------- Overlap ---------------------------  
  template<> string Op_component<Overlap>::get_op_string() const {
    return "OVERLAP";
  }
  template<> void Op_component<Overlap>::build_iterators(SpinBlock& b)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<Overlap>(new Overlap);
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = false;

      m_op(0)[0]->set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
            
    }
  
  template<> std::vector<std::vector<int> > Op_component<Overlap>::get_array() const 
    {
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      return ret_val;
    }
  
  template<> void Op_component<Overlap>::add_local_indices(int i, int j , int k) {};

// three-index opertors
  template<> void Op_component<CreCreDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreDesDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreDesCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreCreCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<DesCreDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<DesDesCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<DesCreCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<DesDesDes>::add_local_indices(int i, int j , int k) {};

// four-index opertors
  template<> void Op_component<CreCreDesDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreDesCreDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreDesDesCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreDesDesDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreCreCreDes>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreCreDesCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreDesCreCre>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<CreCreCreCre>::add_local_indices(int i, int j , int k) {};

// NPDM
  template<> void Op_component<RI3index>::add_local_indices(int i, int j , int k) {};
  template<> void Op_component<RI4index>::add_local_indices(int i, int j , int k) {};

// Added by NN, 2016-06-14
  template<> void Op_component<CDD_sum>::add_local_indices(int, int, int) { }
  template<> void Op_component<CDD_CreDesComp>::add_local_indices(int, int, int) { }
  template<> void Op_component<CDD_DesDesComp>::add_local_indices(int, int, int) { }
  template<> void Op_component<CCD_sum>::add_local_indices(int, int, int) { }
  template<> void Op_component<CCD_CreDesComp>::add_local_indices(int, int, int) { }
  template<> void Op_component<CCD_CreCreComp>::add_local_indices(int, int, int) { }

}
