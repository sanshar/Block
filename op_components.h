/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_OP_COMPONENTS_H
#define SPIN_OP_COMPONENTS_H

#include <boost/function.hpp>
#include <boost/functional.hpp>
#include <boost/format.hpp>
#include <para_array.h>
#include <boost/shared_ptr.hpp>
#include <list>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include "Operators.h"
#include "operatorloops.h"
#include <string>
#include "three_index_ops.h"
#include "four_index_ops.h"
#include "para_array_3d.h"
#include "para_array_4d.h"

namespace SpinAdapted{
class SpinBlock;

//===========================================================================================================================================================
//choose the type of array for different types of Operators

template <class T> struct ChooseArray {
  typedef para_array_1d<std::vector<boost::shared_ptr<SparseMatrix> > > ArrayType;
};
template <> struct ChooseArray<Cre> {
  typedef para_array_1d<std::vector<boost::shared_ptr<Cre> > > ArrayType; // Cre, CreDes, etc. are sparse matrices: <a|a_i^\dagger|b>
};
template <> struct ChooseArray<Des> {
  typedef para_array_1d<std::vector<boost::shared_ptr<Des> > > ArrayType;
};
template <> struct ChooseArray<CreDes> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreDes> > > ArrayType;
};
template <> struct ChooseArray<DesCre> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<DesCre> > > ArrayType;
};
template <> struct ChooseArray<CreCre> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreCre> > > ArrayType;
};
template <> struct ChooseArray<DesDes> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<DesDes> > > ArrayType;
};

template <> struct ChooseArray<CreDesComp> {
    typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreDesComp> > > ArrayType;
};
template <> struct ChooseArray<DesCreComp> {
    typedef para_array_triang_2d<std::vector<boost::shared_ptr<DesCreComp> > > ArrayType;
};
template <> struct ChooseArray<DesDesComp> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<DesDesComp> > > ArrayType;
};
template <> struct ChooseArray<CreCreComp> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreCreComp> > > ArrayType;
};
template <> struct ChooseArray<CreCreDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CreCreDesComp> > > ArrayType;
};
template <> struct ChooseArray<CreDesDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CreDesDesComp> > > ArrayType;
};

template <> struct ChooseArray<Ham> {
  typedef para_array_0d<std::vector<boost::shared_ptr<Ham> > > ArrayType;
};
template <> struct ChooseArray<Overlap> {
  typedef para_array_0d<std::vector<boost::shared_ptr<Overlap> > > ArrayType;
};

template <> struct ChooseArray<RI3index> {
  typedef para_array_3d<std::vector<boost::shared_ptr<RI3index> > > ArrayType;
};  
template <> struct ChooseArray<RI4index> {
  typedef para_array_4d<std::vector<boost::shared_ptr<RI4index> > > ArrayType;
};
// 3PDM
template <> struct ChooseArray<CreCreDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<CreCreDes> > > ArrayType;
};
template <> struct ChooseArray<CreDesDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<CreDesDes> > > ArrayType;
};
template <> struct ChooseArray<CreDesCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<CreDesCre> > > ArrayType;
};
template <> struct ChooseArray<CreCreCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<CreCreCre> > > ArrayType;
};
// 4PDM
template <> struct ChooseArray<DesCreDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<DesCreDes> > > ArrayType;
};
template <> struct ChooseArray<DesDesCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<DesDesCre> > > ArrayType;
};
template <> struct ChooseArray<DesCreCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<DesCreCre> > > ArrayType;
};
template <> struct ChooseArray<DesDesDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<DesDesDes> > > ArrayType;
};
template <> struct ChooseArray<CreCreDesDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreCreDesDes> > > ArrayType;
};
template <> struct ChooseArray<CreDesCreDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreDesCreDes> > > ArrayType;
};
template <> struct ChooseArray<CreDesDesCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreDesDesCre> > > ArrayType;
};
template <> struct ChooseArray<CreDesDesDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreDesDesDes> > > ArrayType;
};
template <> struct ChooseArray<CreCreCreDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreCreCreDes> > > ArrayType;
};
template <> struct ChooseArray<CreCreDesCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreCreDesCre> > > ArrayType;
};
template <> struct ChooseArray<CreDesCreCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreDesCreCre> > > ArrayType;
};
template <> struct ChooseArray<CreCreCreCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<CreCreCreCre> > > ArrayType;
};
// mps_nevpt
template <> struct ChooseArray<CDD_sum> {
  typedef para_array_0d<std::vector<boost::shared_ptr<CDD_sum> > > ArrayType;
};
template <> struct ChooseArray<CDD_CreDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CDD_CreDesComp> > > ArrayType;
};
template <> struct ChooseArray<CDD_DesDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CDD_DesDesComp> > > ArrayType;
};
template <> struct ChooseArray<CCD_sum> {
  typedef para_array_0d<std::vector<boost::shared_ptr<CCD_sum> > > ArrayType;
};
template <> struct ChooseArray<CCD_CreDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CCD_CreDesComp> > > ArrayType;
};
template <> struct ChooseArray<CCD_CreCreComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CCD_CreCreComp> > > ArrayType;
};
//===========================================================================================================================================================

class Op_component_base
{

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & m_core & m_deriv;
  }

 protected:
  bool m_core;
  bool m_deriv;

 public:
  virtual void build_iterators(SpinBlock& b)=0;
  virtual void build_operators(SpinBlock& b)=0;
  virtual void build_csf_operators(std::vector< Csf >& dets, std::vector< std::vector<Csf> >& ladders, SpinBlock& b) = 0;
  virtual void renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo* s) =0;
  virtual void renormalise_transform(const std::vector<Matrix>& leftMat, const StateInfo* bra, const std::vector<Matrix>& rightMat, const StateInfo* ket) =0;
  virtual void build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) =0;
  virtual void build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket) =0;

  //virtual string type_name() = 0;
  virtual int get_size() const =0;
  virtual int size() const=0;
  virtual void clear() =0;
  virtual std::vector<boost::shared_ptr<SparseMatrix> > get_local_element(int i) =0;
  virtual std::vector<boost::shared_ptr<SparseMatrix> > get_global_element(int i)=0;
  const bool &is_core() const {return m_core;}
  const bool &is_deriv() const {return m_deriv;}
  void set_core(bool is_core) {m_core = is_core;}
  virtual void add_local_indices(int i, int j=-1, int k=-1) {};
  virtual bool is_local() const = 0;
  virtual bool& set_local() = 0; 
  virtual std::vector< std::vector<int> > get_array() const =0;
  virtual const std::vector<boost::shared_ptr<SparseMatrix> > get_element(int i, int j=-1, int k=-1, int l=-1) const = 0;
  virtual std::vector<boost::shared_ptr<SparseMatrix> > get_element(int i, int j=-1, int k=-1, int l=-1) = 0;
  virtual bool has(int i, int j=-1, int k=-1, int l=-1) const = 0;
  virtual bool has_local_index(int i, int j=-1, int k=-1, int l=-1) const = 0;
  virtual boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) = 0;
  virtual const boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) const = 0;
  virtual boost::shared_ptr<SparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1) = 0;
  virtual const boost::shared_ptr<SparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1) const = 0;
  virtual std::string get_op_string() const = 0;
  virtual std::string get_filename() const = 0;
  virtual ~Op_component_base() {}  

};

//===========================================================================================================================================================

template <class Op> class Op_component : public Op_component_base
{

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Op_component_base>(*this);
    ar.register_type(static_cast<Op *>(NULL));
    ar & m_op & uniqueID;
  }

 protected:
  typedef typename ChooseArray<Op>::ArrayType paraarray;
  typedef Op OpType; 
  paraarray m_op;
  int uniqueID;
  // Use for unique filename for NPDM disk-based operator storage 
  static int nIDgenerator;

 public:
  Op_component() { m_deriv=false; uniqueID = nIDgenerator++; }
  Op_component(bool core) { m_core=core; m_deriv=false; uniqueID = nIDgenerator++; }
  std::string get_op_string() const;
  bool& set_local() {return m_op.set_local();}
  bool is_local() const {return m_op.is_local();}
  int get_size() const {return m_op.local_nnz();}
  int size() const  {return m_op.global_nnz();}
  bool has(int i, int j=-1, int k=-1, int l=-1) const {return m_op.has(i, j, k, l);}
  bool has_local_index(int i, int j=-1, int k=-1, int l=-1) const {return m_op.has_local_index(i, j, k, l);}
  virtual void add_local_indices(int i, int j=-1, int k=-1);
  void clear(){m_op.clear();}

  void build_iterators(SpinBlock& b);
  void build_operators(SpinBlock& b) 
    { singlethread_build(*this, b); }

  // Note for NPDM higher-index operators, there are template specializations for these functions
  void build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo)
    {for_all_operators_multithread(*this, bind(&SparseMatrix::build_and_renormalise_transform, _1, &b, boost::ref(ot), boost::ref(rotateMatrix), stateinfo));}

  void build_and_renormalise_operators(SpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
    {for_all_operators_multithread(*this, bind(&SparseMatrix::build_and_renormalise_transform, _1, &b, boost::ref(ot), boost::ref(leftMat), bra, boost::ref(rightMat), ket));}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Use for unique filename for NPDM disk-based operator storage
  std::string get_filename() const
  {
    std::string file;
    file = str( boost::format("%s%s%s%s%d%s%d%s") % dmrginp.load_prefix() % "/" % get_op_string() % "_" % uniqueID % "_p" % mpigetrank() % ".tmp" );
    return file;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void build_csf_operators(std::vector< Csf >& c, vector< vector<Csf> >& ladders, SpinBlock& b) 
  {
    for_all_operators_multithread( *this, bind(&SparseMatrix::buildUsingCsf, _1, boost::ref(b), boost::ref(ladders), boost::ref(c)) );
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo* s)
  {
    if ( m_op.num_indices() > 2 ) {
      // renormalise_transform_on_disk( rotateMatrix, s );
      assert(false);
    }
    else {
      // For operators built in core
      for_all_operators_multithread( *this, bind(&SparseMatrix::renormalise_transform, _1, boost::ref(rotateMatrix), s) );
    }
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void renormalise_transform( const std::vector<Matrix>& leftMat, const StateInfo* bra, const std::vector<Matrix>& rightMat, const StateInfo* ket)
  {
    if ( m_op.num_indices() > 2 ) {
      // renormalise_transform_on_disk( ... )
      assert(false);
    }
    else {
      // For operators built in core
      for_all_operators_multithread( *this, bind(&SparseMatrix::renormalise_transform, _1, boost::ref(leftMat), bra, boost::ref(rightMat), ket) );
    }
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<boost::shared_ptr<SparseMatrix> > get_local_element(int i) 
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op.get_local_element(i).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op.get_local_element(i)[l]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<boost::shared_ptr<SparseMatrix> > get_global_element(int i)
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op.get_global_element(i).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op.get_global_element(i)[l]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector< std::vector<int> > get_array() const 
  {
    std::vector<int> orbs(2);
    std::vector< std::vector<int> > ret_val(m_op.local_nnz());
    for (int i=0; i<m_op.local_nnz(); i++)
      ret_val[i] = m_op.unmap_local_index(i);
    return ret_val;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  const std::vector<boost::shared_ptr<SparseMatrix> >  get_element(int i, int j=-1, int k=-1, int l=-1) const 
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op(i,j,k,l).size());
    for (int p=0; p<vec.size(); p++)
      vec[p] = m_op(i,j,k,l)[p]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<boost::shared_ptr<SparseMatrix> >  get_element(int i, int j=-1, int k=-1, int l=-1)
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op(i,j,k,l).size());
    for (int p=0; p<vec.size(); p++)
      vec[p] = m_op(i,j,k,l)[p]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1)
  {
    Op* o = 0;
    std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++) {
      if (s == vec[p]->get_deltaQuantum())
	    return m_op(i,j,k,l)[p];
    }
    return boost::shared_ptr<Op>(o);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  const boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) const
  {
    Op* o = 0;
    const std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++)
      if (s == vec[p]->get_deltaQuantum())
	    return m_op(i,j,k,l)[p];
    return boost::shared_ptr<Op>(o);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  boost::shared_ptr<SparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1)
  {
    Op* o = 0;
    std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++) {
      if (s == vec[p]->get_quantum_ladder())
	    return m_op(i,j,k,l)[p];
    }
    return boost::shared_ptr<Op>(o);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  const boost::shared_ptr<SparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1) const
  {
    Op* o = 0;
    const std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++)
      if (s == vec[p]->get_quantum_ladder())
	    return m_op(i,j,k,l)[p];
    return boost::shared_ptr<Op>(o);
  }

};

//===========================================================================================================================================================
 
template <class Op> int Op_component<Op>::nIDgenerator = 1;

}

#endif
