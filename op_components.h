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
#include <para_array.h>
#include <boost/shared_ptr.hpp>
#include <list>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include "Operators.h"
#include "operatorloops.h"
#include <string>

namespace SpinAdapted{
class SpinBlock;

//********************************
//choose the type of array for different types of Operators
template <class T> struct ChooseArray {
  typedef para_array_1d<std::vector<boost::shared_ptr<SparseMatrix> > > ArrayType;
};
template <> struct ChooseArray<Cre> {
  typedef para_array_1d<std::vector<boost::shared_ptr<Cre> > > ArrayType; // Cre, CreDes, etc. are sparse matrices: <a|a_i^\dagger|b>
};
template <> struct ChooseArray<CreDes> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreDes> > > ArrayType;
};
template <> struct ChooseArray<CreCre> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreCre> > > ArrayType;
};
template <> struct ChooseArray<CreDesComp> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<CreDesComp> > > ArrayType;
};
template <> struct ChooseArray<DesDesComp> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<DesDesComp> > > ArrayType;
};
template <> struct ChooseArray<CreCreDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<CreCreDesComp> > > ArrayType;
};
template <> struct ChooseArray<Ham> {
  typedef para_array_0d<std::vector<boost::shared_ptr<Ham> > > ArrayType;
};
//*************************************
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
  virtual void build_operators(SpinBlock& b)=0;
  virtual void build_csf_operators(std::vector< Csf >& dets, std::vector< std::vector<Csf> >& ladders, SpinBlock& b) = 0;
  virtual void build_iterators(SpinBlock& b)=0;
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
  //virtual std::vector<SparseMatrix*> get_element(int i, int j=-1, int k=-1) = 0;
  virtual const std::vector<boost::shared_ptr<SparseMatrix> > get_element(int i, int j=-1, int k=-1) const = 0;
  virtual std::vector<boost::shared_ptr<SparseMatrix> > get_element(int i, int j=-1, int k=-1) = 0;
  virtual bool has(int i, int j = -1, int k = -1) const = 0;
  virtual bool has_local_index(int i, int j=-1, int k=-1) const = 0;
  virtual std::vector< std::vector<int> > get_array() const =0;
  virtual boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1) = 0;
  virtual const boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1) const = 0;
  virtual std::string get_op_string() const = 0;
  virtual ~Op_component_base() {}  
};

template <class Op> class Op_component : public Op_component_base
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Op_component_base>(*this);
    ar.register_type(static_cast<Op *>(NULL));
    ar & m_op;
  }
 protected:
  typedef typename ChooseArray<Op>::ArrayType paraarray;
  typedef Op OpType; 
  paraarray m_op;
 public:
  Op_component() {m_deriv=false;}
  Op_component(bool core) {m_core=core;m_deriv=false;}
  bool& set_local() {return m_op.set_local();}
  bool is_local() const {return m_op.is_local();}
  int get_size() const {return m_op.local_nnz();}
  int size() const  {return m_op.global_nnz();}
  bool has(int i, int j=-1, int k=-1) const {return m_op.has(i, j, k);}
  bool has_local_index(int i, int j=-1, int k=-1) const {return m_op.has_local_index(i, j, k);}
  virtual void add_local_indices(int i, int j=-1, int k=-1){};
  std::vector<boost::shared_ptr<SparseMatrix> > get_local_element(int i) 
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op.get_local_element(i).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op.get_local_element(i)[l]; 
    return vec;
  }
  std::vector<boost::shared_ptr<SparseMatrix> > get_global_element(int i)
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op.get_global_element(i).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op.get_global_element(i)[l]; 
    return vec;
  }

  std::vector< std::vector<int> > get_array() const 
  {
    std::vector<int> orbs(2);
    std::vector< std::vector<int> > ret_val(m_op.local_nnz());
    for (int i=0; i<m_op.local_nnz(); i++)
      {
	pair<int, int> opair = m_op.unmap_local_index(i);
	orbs[0] = opair.first; orbs[1] = opair.second;
	ret_val[i] = orbs;
      }
    return ret_val;
  }

  void clear(){m_op.clear();}
  void build_iterators(SpinBlock& b);
  void build_operators(SpinBlock& b) {singlethread_build(*this, b);}
  void build_csf_operators(std::vector< Csf >& c, vector< vector<Csf> >& ladders, SpinBlock& b) {singlethread_build(*this, b, c, ladders);}
  const std::vector<boost::shared_ptr<SparseMatrix> >  get_element(int i, int j=-1, int k=-1) const 
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op(i,j,k).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op(i,j,k)[l]; 
    return vec;
  }
  std::vector<boost::shared_ptr<SparseMatrix> >  get_element(int i, int j=-1, int k=-1)
  {
    std::vector<boost::shared_ptr<SparseMatrix> > vec(m_op(i,j,k).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op(i,j,k)[l]; 
    return vec;
  }
  boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1)
  {
    Op* o = 0;
    std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k);
    for (int l=0; l<vec.size(); l++) {
      if (boost::iequals(s, vec[l]->get_deltaQuantum()))
	    return m_op(i,j,k)[l];
    }
    return boost::shared_ptr<Op>(o);
  }
  const boost::shared_ptr<SparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1) const
  {
    Op* o = 0;
    const std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k);
    for (int l=0; l<vec.size(); l++)
      if (boost::iequals(s, vec[l]->get_deltaQuantum()))
	    return m_op(i,j,k)[l];
    return boost::shared_ptr<Op>(o);
  }
  std::string get_op_string() const;

};


 
}


#endif
