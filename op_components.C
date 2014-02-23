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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
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
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  template<> void Op_component<Cre>::add_local_indices(int i, int j)
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<Cre> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<Cre>(new Cre);
    }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------- C_S1 ---------------------------  

  template<> string Op_component<Cre>::get_op_string() const {
    return "CRE";
  }

  template<> void Op_component<Cre>::build_iterators(SpinBlock& b)
  {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.screen_tol();
      std::vector<int> screened_c_ix = screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1, *b.get_twoInt(), screen_tol); 
      m_op.set_indices(screened_c_ix, dmrginp.last_site());  
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
        op.set_deltaQuantum() = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
        //op.set_deltaQuantum() = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOf(orbs[0]));      
        //MAW
        op.set_quantum_ladder()["(C)"] = { op.get_deltaQuantum() };
      }
      
  }
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------- C_S1 ---------------------------  

  template<> string Op_component<Des>::get_op_string() const {
    return "DES";
  }

  template<> void Op_component<Des>::build_iterators(SpinBlock& b)
  {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.screen_tol();
      std::vector<int> screened_c_ix = screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1, *b.get_twoInt(), screen_tol); 
      m_op.set_indices(screened_c_ix, dmrginp.last_site());  
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
        // NOte minus sign
        op.set_deltaQuantum() = - SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
        //op.set_deltaQuantum() = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOf(orbs[0]));      
        //MAW
        op.set_quantum_ladder()["(D)"] = { op.get_deltaQuantum() };
      }
      
  }
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------- Cd_ ---------------------------  

  template<> string Op_component<CreDes>::get_op_string() const {
    return "CREDES";
  }

  template<> void Op_component<CreDes>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.screen_tol();
      vector< pair<int, int> > screened_cd_ix = screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
//MAW	  pair<int, int> opair = m_op.unmap_local_index(i);
//MAW	  orbs[0] = opair.first; orbs[1] = opair.second;
	  orbs = m_op.unmap_local_index(i);
//MAW <<<<<<
	  std::vector<boost::shared_ptr<CreDes> >& vec = m_op.get_local_element(i);
assert( vec.size() == 0);
	  SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1-spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<CreDes>(new CreDes);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum() = spinvec[j];      
//MAW
       op.set_quantum_ladder()["(CD)"] = { op.get_deltaQuantum() };
       assert( op.get_deltaQuantum().particleNumber == 0 );
	  }
	}
    }
  
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------- Cc_ ---------------------------  

  template<> string Op_component<CreCre>::get_op_string() const {
    return "CRECRE";
  }

  template<> void Op_component<CreCre>::build_iterators(SpinBlock& b) {

    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.screen_tol();
    
    vector< pair<int, int> > screened_dd_ix = screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
    m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
    //cout << "New set of CC operators: p" << mpigetrank() << "; local size = " << m_op.local_nnz() << "; global size = " << m_op.global_nnz() << std::endl;
    std::vector<int> orbs(2);
    for (int i = 0; i < m_op.local_nnz(); ++i) {
      orbs = m_op.unmap_local_index(i);
      //cout << "p" << mpigetrank() << "; Orbs = " << orbs[0] << " " << orbs[1] << std::endl;
      std::vector<boost::shared_ptr<CreCre> >& vec = m_op.get_local_element(i);
      SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
      SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
      std::vector<SpinQuantum> spinvec = spin1+spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) {
        vec[j]=boost::shared_ptr<CreCre>(new CreCre);
        SparseMatrix& op = *vec[j];
        op.set_orbs() = orbs;
        op.set_initialised() = true;
        op.set_fermion() = false;
        op.set_deltaQuantum() = spinvec[j];      
        //MAW NPDM
        op.set_quantum_ladder()["(CC)"] = { op.get_deltaQuantum() };
        assert( op.get_deltaQuantum().particleNumber == 2 );
      }
    }
  }
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------- Cdcomp_ ---------------------------  

  template<> string Op_component<CreDesComp>::get_op_string() const {
    return "CREDES_COMP";
  }
  template<> void Op_component<CreDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.screen_tol();
      vector< pair<int, int> > screened_cd_ix = screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
//MAW	  pair<int, int> opair = m_op.unmap_local_index(i);
//MAW	  orbs[0] = opair.first; orbs[1] = opair.second;
	  orbs = m_op.unmap_local_index(i);
//MAW <<<<<<
	  std::vector<boost::shared_ptr<CreDesComp> >& vec = m_op.get_local_element(i);
	  //SpinQuantum spin1 = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOfSpatialOrb(orbs[0]));
	  //SpinQuantum spin2 = SpinQuantum(1, SpinOf(orbs[1]), SymmetryOfSpatialOrb(orbs[1]));
	  SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<CreDesComp>(new CreDesComp);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum() = spinvec[j];      
//MAW
       op.set_quantum_ladder()["CDcomp"] = { op.get_deltaQuantum() };
	  }
	}
    }
  
  template<> void Op_component<CreDesComp>::add_local_indices(int i, int j)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<CreDesComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin2-spin1;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) 
	vec[j]=boost::shared_ptr<CreDesComp>(new CreDesComp);
    }
  
  
  
  // -------------------- Ddcomp_ ---------------------------  
  template<> string Op_component<DesDesComp>::get_op_string() const {
    return "DESDES_COMP";
  }
  template<> void Op_component<DesDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.screen_tol();
      vector< pair<int, int> > screened_dd_ix = screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
//MAW	  pair<int, int> opair = m_op.unmap_local_index(i);
//MAW	  orbs[0] = opair.first; orbs[1] = opair.second;
	  orbs = m_op.unmap_local_index(i);
//MAW <<<<<<
	  std::vector<boost::shared_ptr<DesDesComp> >& vec = m_op.get_local_element(i);
	  //SpinQuantum spin1 = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOfSpatialOrb(orbs[0]));
	  //SpinQuantum spin2 = SpinQuantum(1, SpinOf(orbs[1]), SymmetryOfSpatialOrb(orbs[1]));
	  SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<DesDesComp>(new DesDesComp);
	    SparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum() = -spinvec[j];      
//MAW
       op.set_quantum_ladder()["DDcomp"] = { op.get_deltaQuantum() };
	  }
	}
      
    }
  
  template<> void Op_component<DesDesComp>::add_local_indices(int i, int j)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<DesDesComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1+spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) 
	vec[j]=boost::shared_ptr<DesDesComp>(new DesDesComp);
    }
  
  
  // -------------------- Ccdcomp_ ---------------------------  
  template<> string Op_component<CreCreDesComp>::get_op_string() const {
    return "CRECREDES_COMP";
  }
  template<> void Op_component<CreCreDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.screen_tol();
      vector< int > screened_cdd_ix = screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1, *b.get_twoInt(), screen_tol);
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
	  //op.set_deltaQuantum() = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOfSpatialOrb(orbs[0]) );      
	  op.set_deltaQuantum() = SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );      
//MAW
     op.set_quantum_ladder()["CCDcomp"] = { op.get_deltaQuantum() };
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
      m_op(0)[0]->set_deltaQuantum() = SpinQuantum(0, 0, IrrepSpace(0) );      
//MAW
      m_op(0)[0]->set_quantum_ladder()["Ham"] = { m_op(0)[0]->get_deltaQuantum() };
    }
  
  template<> std::vector<std::vector<int> > Op_component<Ham>::get_array() const 
    {
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      return ret_val;
    }
  


}
