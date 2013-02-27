#ifndef SPIN_OPERATORS_HEADER
#define SPIN_OPERATORS_HEADER
#include "BaseOperator.h"
#include <boost/function.hpp>
#include <boost/functional.hpp>


typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor;

typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&,
                              std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor2;

namespace SpinAdapted{


class Cre: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  Cre(opTypes index = 0) { orbs.resize(1); fermion = true; m_state_index = index & GENERIC_MASK;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreDes: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  CreDes(opTypes index = 0) { orbs.resize(2); fermion = false; m_state_index = index & GENERIC_MASK;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CreCre: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  CreCre(opTypes index = 0) { orbs.resize(2); fermion = false; m_state_index = index & GENERIC_MASK;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreDesComp: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  CreDesComp(opTypes index = 0) { orbs.resize(2); fermion = false; m_state_index = index & GENERIC_MASK;}
  void build_fromcd(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class DesDesComp: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  DesDesComp(opTypes index = 0) { orbs.resize(2); fermion = false; m_state_index = index & GENERIC_MASK;}
  void build_fromcc(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreCreDesComp: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  CreCreDesComp(opTypes index = 0) { orbs.resize(1); fermion = true; m_state_index = index & GENERIC_MASK;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class Ham: public SpinAdapted::SparseMatrix
{
 private:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & m_state_index; }
 private:
  opTypes m_state_index;
 public:
  Ham(opTypes index = 0) { fermion = false; m_state_index = index & GENERIC_MASK;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

}


#endif
