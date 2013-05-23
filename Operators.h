#ifndef SPIN_OPERATORS_HEADER
#define SPIN_OPERATORS_HEADER
#include "BaseOperator.h"
#include <boost/function.hpp>
#include <boost/functional.hpp>


typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::SparseMatrix> >&)> Functor;

namespace SpinAdapted{


class Cre: public SpinAdapted::SparseMatrix
{
 public:
  Cre() { orbs.resize(1); fermion = true;}
//MAW
  void build_in_csf_space(const SpinBlock& b) {assert(false);}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreDes: public SpinAdapted::SparseMatrix
{
 public:
  CreDes() { orbs.resize(2); fermion = false;}
//MAW
  void build_in_csf_space(const SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CreCre: public SpinAdapted::SparseMatrix
{
 public:
  CreCre() { orbs.resize(2); fermion = false;}
//MAW
  void build_in_csf_space(const SpinBlock& b) ;
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CreDesComp() { orbs.resize(2); fermion = false;}
//MAW
  void build_in_csf_space(const SpinBlock& b) {assert(false);}
  void build_fromcd(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class DesDesComp: public SpinAdapted::SparseMatrix
{
 public:
  DesDesComp() { orbs.resize(2); fermion = false;}
//MAW
  void build_in_csf_space(const SpinBlock& b) {assert(false);}
  void build_fromcc(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreCreDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CreCreDesComp() { orbs.resize(1); fermion = true;}
//MAW
  void build_in_csf_space(const SpinBlock& b) {assert(false);}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class Ham: public SpinAdapted::SparseMatrix
{
 public:
  Ham() { fermion = false;}
//MAW
  void build_in_csf_space(const SpinBlock& b) {assert(false);}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

}


#endif
