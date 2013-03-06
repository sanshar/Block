#ifndef SPIN_OPERATORS_HEADER
#define SPIN_OPERATORS_HEADER
#include "BaseOperator.h"
#include <boost/function.hpp>
#include <boost/functional.hpp>


typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor;

typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&,
                              std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor2;

namespace SpinAdapted {

//
// FIXME:
// Since SparseMatrix::build(const SpinBlock&) and Op_component::build_operators(SpinBlock&) both depend on state index,
// Op_component also has state index, and follows must be the copy of it (but this design is not good for OOP).
//

class GenericOperator : public SpinAdapted::SparseMatrix
{
private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<SparseMatrix>(*this);
    ar & m_state_index;
  }
protected:
  opTypes m_state_index;
public:
  GenericOperator(opTypes index = 0) { orbs.resize(1); fermion = true; m_state_index = index & GENERIC_MASK; }
  virtual void build(const SpinBlock& b) = 0;
  virtual boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) = 0;
  virtual double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b) = 0;
  virtual boost::shared_ptr<SparseMatrix> deepCopy() const = 0;
};

class Cre: public SpinAdapted::GenericOperator
{
public:
  Cre(opTypes index = 0) : GenericOperator(index) { orbs.resize(1); fermion = true; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new Cre(*this)); }
};


class CreDes: public SpinAdapted::GenericOperator
{
public:
  CreDes(opTypes index = 0) : GenericOperator(index) { orbs.resize(2); fermion = false; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new CreDes(*this)); }
};

class CreCre: public SpinAdapted::GenericOperator
{
public:
  CreCre(opTypes index = 0) : GenericOperator(index) { orbs.resize(2); fermion = false; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new CreCre(*this)); }
};


class CreDesComp: public SpinAdapted::GenericOperator
{
public:
  CreDesComp(opTypes index = 0) : GenericOperator(index) { orbs.resize(2); fermion = false; }
  void build_fromcd(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new CreDesComp(*this)); }
};


class DesDesComp: public SpinAdapted::GenericOperator
{
public:
  DesDesComp(opTypes index = 0) : GenericOperator(index) { orbs.resize(2); fermion = false; }
  void build_fromcc(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new DesDesComp(*this)); }
};


class CreCreDesComp: public SpinAdapted::GenericOperator
{
public:
  CreCreDesComp(opTypes index = 0) : GenericOperator(index) { orbs.resize(1); fermion = true; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new CreCreDesComp(*this)); }
};

class Ham: public SpinAdapted::GenericOperator
{
public:
  Ham(opTypes index = 0) : GenericOperator(index) { fermion = false; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
  boost::shared_ptr<SparseMatrix> deepCopy() const { return boost::shared_ptr<SparseMatrix>(new Ham(*this)); }
};

}


#endif
