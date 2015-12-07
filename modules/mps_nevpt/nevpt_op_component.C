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

  
  // -------------------- CDD_sum ---------------------------  
  template<> string Op_component<CDD_sum>::get_op_string() const {
    return "CDD_SUM";
  }
  template<> void Op_component<CDD_sum>::build_iterators(SpinBlock& b)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<CDD_sum>(new CDD_sum);
      //m_op(0)[0]->set_orbs() = std::vector<int>(1,b.nonactive_orb(0));
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = true;
      m_op(0)[0]->set_deltaQuantum(1, -getSpinQuantum(b.nonactive_orb(0)));
    }
  
  // -------------------- CDD_CreDescomp_ ---------------------------  
  template<> string Op_component<CDD_CreDesComp>::get_op_string() const {
    return "CDD_CREDES_COMP";
  }
  template<> void Op_component<CDD_CreDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      double screen_tol = dmrginp.oneindex_screen_tol();
    //screen_tol = 0.0;
      std::vector<int> screened_d_ix = screened_cdd_d_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
      m_op.set_indices(screened_d_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	    {
	      orbs[0] = m_op.get_local_indices()[i];
	      std::vector<boost::shared_ptr<CDD_CreDesComp> >& vec = m_op.get_local_element(i);
	      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	      SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
	      std::vector<SpinQuantum> spinvec = spin2-spin1;
        vec.resize(spinvec.size());
        for(int j=0; j < spinvec.size(); j++)
        {
	        vec[j]=boost::shared_ptr<CDD_CreDesComp>(new CDD_CreDesComp);
	        SparseMatrix& op = *m_op.get_local_element(i)[j];
	        op.set_orbs() = orbs;
	        op.set_initialised() = true;
	        op.set_fermion() = false;
	        op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
        }
	    }
    }
  
  // -------------------- CDD_DesDescomp_ ---------------------------  
  template<> string Op_component<CDD_DesDesComp>::get_op_string() const {
    return "CDD_DESDES_COMP";
  }
  template<> void Op_component<CDD_DesDesComp>::build_iterators(SpinBlock& b)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    double screen_tol = dmrginp.oneindex_screen_tol();
    //screen_tol = 0.0;
    std::vector<int> screened_c_ix = screened_cdd_c_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
    m_op.set_indices(screened_c_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    for (int i = 0; i < m_op.local_nnz(); ++i)
  	{
  	  orbs[0] = m_op.get_local_indices()[i];
  	  std::vector<boost::shared_ptr<CDD_DesDesComp> >& vec = m_op.get_local_element(i);
  	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
  	  SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
  	  std::vector<SpinQuantum> spinvec = spin2+spin1;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
  	    vec[j]=boost::shared_ptr<CDD_DesDesComp>(new CDD_DesDesComp);
  	    SparseMatrix& op = *m_op.get_local_element(i)[j];
  	    op.set_orbs() = orbs;
  	    op.set_initialised() = true;
  	    op.set_fermion() = false;
  	    op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
      }
  	}
  }
  
  
  // -------------------- CCD_sum ---------------------------  
  template<> string Op_component<CCD_sum>::get_op_string() const {
    return "CCD_SUM";
  }
  template<> void Op_component<CCD_sum>::build_iterators(SpinBlock& b)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<CCD_sum>(new CCD_sum);
      //m_op(0)[0]->set_orbs() = std::vector<int>(1,b.nonactive_orb(0));
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = true;
      m_op(0)[0]->set_deltaQuantum(1, getSpinQuantum(b.nonactive_orb(0)));
    }
  
  // -------------------- CCD_CreDescomp_ ---------------------------  
  template<> string Op_component<CCD_CreDesComp>::get_op_string() const {
    return "CCD_CREDES_COMP";
  }
  template<> void Op_component<CCD_CreDesComp>::build_iterators(SpinBlock& b)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      std::vector<int> screened_c_ix = screened_ccd_c_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);
      m_op.set_indices(screened_c_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	    {
	      orbs[0] = m_op.get_local_indices()[i];
	      std::vector<boost::shared_ptr<CCD_CreDesComp> >& vec = m_op.get_local_element(i);
	      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	      SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
	      std::vector<SpinQuantum> spinvec = spin1-spin2;
        vec.resize(spinvec.size());
        for(int j=0; j < spinvec.size(); j++)
        {
	        vec[j]=boost::shared_ptr<CCD_CreDesComp>(new CCD_CreDesComp);
	        SparseMatrix& op = *m_op.get_local_element(i)[j];
	        op.set_orbs() = orbs;
	        op.set_initialised() = true;
	        op.set_fermion() = false;
	        op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
        }
	    }
    }
  
  // -------------------- CCD_CreCrecomp_ ---------------------------  
  template<> string Op_component<CCD_CreCreComp>::get_op_string() const {
    return "CCD_DesDES_COMP";
  }
  template<> void Op_component<CCD_CreCreComp>::build_iterators(SpinBlock& b)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.oneindex_screen_tol();
    std::vector<int> screened_d_ix = screened_ccd_d_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);
    m_op.set_indices(screened_d_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    for (int i = 0; i < m_op.local_nnz(); ++i)
  	{
  	  orbs[0] = m_op.get_local_indices()[i];
  	  std::vector<boost::shared_ptr<CCD_CreCreComp> >& vec = m_op.get_local_element(i);
  	  SpinQuantum spin1 = -getSpinQuantum(orbs[0]);
  	  SpinQuantum spin2 = -getSpinQuantum(b.nonactive_orb(0));
  	  std::vector<SpinQuantum> spinvec = spin2+spin1;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
  	    vec[j]=boost::shared_ptr<CCD_CreCreComp>(new CCD_CreCreComp);
  	    SparseMatrix& op = *m_op.get_local_element(i)[j];
  	    op.set_orbs() = orbs;
  	    op.set_initialised() = true;
  	    op.set_fermion() = false;
  	    op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
      }
  	}
  }
  
}
