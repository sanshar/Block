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
  Cre() { orbs.resize(1); fermion = true; build_pattern = "(C)"; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class Des: public SpinAdapted::SparseMatrix
{
 public:
  Des() { orbs.resize(1); fermion = true; build_pattern = "(D)"; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreDes: public SpinAdapted::SparseMatrix
{
 public:
  CreDes() { orbs.resize(2); fermion = false; build_pattern = "(CD)"; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class DesCre: public SpinAdapted::SparseMatrix
{
 public:
  DesCre() { orbs.resize(2); fermion = false; build_pattern = "(DC)"; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CreCre: public SpinAdapted::SparseMatrix
{
 public:
  CreCre() { orbs.resize(2); fermion = false; build_pattern = "(CC)"; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class DesDes: public SpinAdapted::SparseMatrix
{
 public:
  DesDes() { orbs.resize(2); fermion = false; build_pattern = "(DD)"; }
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CreDesComp() { orbs.resize(2); fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class DesCreComp: public SpinAdapted::SparseMatrix
{
 public:
  DesCreComp() { orbs.resize(2); fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class DesDesComp: public SpinAdapted::SparseMatrix
{
 public:
  DesDesComp() { orbs.resize(2); fermion = false;}
  void build_fromcc(SpinBlock& b);
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CreCreComp: public SpinAdapted::SparseMatrix
{
 public:
  CreCreComp() { orbs.resize(2); fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


class CreCreDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CreCreDesComp() { orbs.resize(1); fermion = true;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CreDesDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CreDesDesComp() { orbs.resize(1); fermion = true;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class Ham: public SpinAdapted::SparseMatrix
{
 public:
  Ham() { fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class Overlap: public SpinAdapted::SparseMatrix
{
 public:
  Overlap() { fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//CDD_sum is used in V_a subspace of nevpt2.
//It is the sum of perturber operator on one block.
class CDD_sum: public SpinAdapted::SparseMatrix
{
 public:
  CDD_sum() { fermion = true;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CDD_CreDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CDD_CreDesComp() { fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CDD_DesDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CDD_DesDesComp() { fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//CCD_sum is used in V_i subspace of nevpt2.
//It is the sum of perturber operator on one block.
class CCD_sum: public SpinAdapted::SparseMatrix
{
 public:
  CCD_sum() { fermion = true;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CCD_CreDesComp: public SpinAdapted::SparseMatrix
{
 public:
  CCD_CreDesComp() { fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

class CCD_CreCreComp: public SpinAdapted::SparseMatrix
{
 public:
  CCD_CreCreComp() { fermion = false;}
  void build(const SpinBlock& b) ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};


double vcccc_4idx_asymm(int i, int j, int k, int l);

}


#endif
