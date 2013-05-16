/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_BASEOPERATOR_HEADER
#define SPIN_BASEOPERATOR_HEADER
#include <cmath>
#include <boost/function.hpp>
#include <boost/functional.hpp>
#include <boost/bind.hpp>
#include <boostutils.h>
#include "SpinQuantum.h"
#include "ObjectMatrix.h"
#include "csf.h"
#include "StateInfo.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

class TensorOp;
namespace SpinAdapted{
class SpinBlock;

//enum opTypes{ HAM, CRE, CRE_CRE, DES_DESCOMP, CRE_DES, CRE_DESCOMP, CRE_CRE_DESCOMP};
//MAW
enum opTypes{ HAM, CRE, CRE_CRE, DES_DESCOMP, CRE_DES, CRE_DESCOMP, CRE_CRE_DESCOMP, CRE_CRE_DES };

enum CompType{CD, DD, CCD, C};


template<class T> class Baseoperator
{
 public:
  virtual bool get_fermion() const = 0;
  virtual int nrows() const = 0;
  virtual int ncols() const = 0;
  virtual bool get_initialised() const = 0;
  virtual const std::vector<int>& get_orbs() const = 0;
  virtual int get_orbs(int i) const = 0;
  virtual const char& allowed(int i, int j) const = 0;
  virtual char& allowed(int i, int j) = 0;
  virtual T& operator_element(int i, int j) = 0;
  virtual const T& operator_element(int i, int j) const = 0;
  virtual T& operator()(int i, int j) = 0;
  virtual const T& operator()(int i, int j) const = 0;
  virtual SpinQuantum get_deltaQuantum() const = 0;
  virtual char conjugacy() const = 0;
  virtual ~Baseoperator() {};
  virtual int get_spin() const = 0;
  virtual IrrepSpace get_symm() const = 0;
  virtual double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const = 0;
  Baseoperator() {};
};


class SparseMatrix : public Baseoperator<Matrix>
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & orbs \
	& deltaQuantum \
	& fermion \
	& initialised \
	& built \
	& allowedQuantaMatrix \
	& operatorMatrix & Sign;
    }

 protected:
  std::vector<int> orbs;
  bool fermion;
  ObjectMatrix<char> allowedQuantaMatrix;
  bool initialised;
  bool built;
  ObjectMatrix<Matrix> operatorMatrix;
  SpinQuantum deltaQuantum;
  int Sign;
 public:
 SparseMatrix() : orbs(2), initialised(false), built(false), Sign(1), fermion(false){};
  virtual ~SparseMatrix(){};
  int nrows() const { return allowedQuantaMatrix.nrows(); }
  int ncols() const { return allowedQuantaMatrix.ncols(); }
  int nrows(char conj) const { return allowedQuantaMatrix.Nrows(conj); }
  char conjugacy() const { return 'n'; }
  int get_sign() const {return Sign;}
  int ncols(char conj) const { return allowedQuantaMatrix.Ncols(conj); }
  bool get_initialised() const { return initialised; }
  bool &set_initialised() { return initialised; }
  bool get_fermion() const { return fermion; }
  bool &set_fermion() { return fermion; }
  const char& allowed(int i, int j) const { return allowedQuantaMatrix(i, j); }
  char& allowed(int i, int j) { return allowedQuantaMatrix(i, j); }
  SpinQuantum &set_deltaQuantum() { return deltaQuantum; }
  int get_spin()const  { return deltaQuantum.get_s();}
  IrrepSpace get_symm()const  { return deltaQuantum.get_symm();}
  SpinQuantum get_deltaQuantum() const { return deltaQuantum; }
  int get_orbs(int i) const 
  { 
    if(i >= orbs.size())
      return -1;
    else
      return orbs[i]; 
  }
  double memoryUsed(const SpinBlock& b);
  const std::vector<int>& get_orbs() const { return orbs; }
  std::vector<int>& set_orbs() { return orbs; }
  const bool& get_built() const { return built; }
  bool& set_built() { return built; }  
  double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const {return 1.0;}

  void resize(int n, int c) { operatorMatrix.ReSize(n, c); allowedQuantaMatrix.ReSize(n, c); }
  const Matrix& operator_element(int i, int j) const { return operatorMatrix(i, j); }
  const Matrix& operator()(int i, int j) const { return operatorMatrix(i, j); }
  Matrix& operator_element(int i, int j) { return operatorMatrix(i, j); }
  Matrix& operator()(int i, int j) { return operatorMatrix(i, j); }
  void allocate(const StateInfo& s);
  void allocate(const StateInfo& sr, const StateInfo& sc);
  void allocate(const SpinBlock& b);
  virtual boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) =0;
  virtual void build(const SpinBlock& b) =0;
  void buildUsingCsf(const SpinBlock& b, vector< vector<Csf> >& ladders, std::vector< Csf >& s) ;
  virtual double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b=0)=0;
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const TwoElectronArray& v_2);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2);
  bool nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb);
  std::vector<double> calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2);
  friend ostream& operator<<(ostream& os, const SparseMatrix& a);
  void Randomise();
  void Normalise(int* success);
  void Clear();
  void CleanUp();

  void Save(std::ofstream &ofs) const
  {
    boost::archive::binary_oarchive save_op(ofs);
    save_op << *this;
  }
  void Load(std::istream &ifs)
  {
    boost::archive::binary_iarchive load_op(ifs);
    load_op >> *this;
  }
  void OperatorMatrixReference(ObjectMatrix<Matrix*>& m, const std::vector<int>& oldToNewStateI, const std::vector<int>& oldToNewStateJ);

  void renormalise_transform(const std::vector<Matrix>& rotate_matrix, const StateInfo *stateinfo);
  void build_and_renormalise_transform(SpinBlock *big, const opTypes &ot, const std::vector<Matrix>& rotate_matrix, 
				       const StateInfo *newStateInfo);
  SparseMatrix& operator+=(const SparseMatrix& other);
};

class Transposeview : public SparseMatrix
{
private:
  boost::shared_ptr<SparseMatrix> opdata; 
public:
  Transposeview(const boost::shared_ptr<SparseMatrix>& opptr) : opdata(opptr) {}
  Transposeview(SparseMatrix& op) { opdata = boost::shared_ptr<SparseMatrix>(&op, boostutils::null_deleter());}
  SpinQuantum get_deltaQuantum() const { return -opdata->get_deltaQuantum(); }  
  bool get_fermion() const { return opdata->get_fermion(); }
  bool get_initialised() const { return opdata->get_initialised(); }
  int nrows() const { return opdata->ncols(); }
  int ncols() const { return opdata->nrows(); }
  const char &allowed(int i, int j) const { return opdata->allowed(j, i); }
  char &allowed(int i, int j) { return opdata->allowed(j, i); }
  const Matrix& operator_element(int i, int j) const { return opdata->operator_element(j, i); }
  Matrix& operator_element(int i, int j) { return opdata->operator_element(j, i); }
  int get_spin()const  { return opdata->get_deltaQuantum().get_s();}
  IrrepSpace get_symm()const  { return -opdata->get_deltaQuantum().get_symm();}
  int get_orbs(int i) const {return opdata->get_orbs(i);}
  const std::vector<int>& get_orbs() const { return opdata->get_orbs(); }
  const Matrix& operator()(int i, int j) const { return opdata->operator()(j, i); }
  Matrix& operator()(int i, int j) { return opdata->operator()(j, i); }
  char conjugacy() const { if (opdata->conjugacy() == 'n') return 't'; else return 'n';}
  double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const 
  {
    int rights = rightq.get_s(), lefts = leftq.get_s();
    int leftirrep = leftq.get_symm().getirrep(), rightirrep = rightq.get_symm().getirrep();
    int lS, rS, lirrep, rirrep;
    double scale = 1.0;
    if (conjugacy() == 'n') {return 1.0;}
    else {lS = rights; rS = lefts; lirrep = rightirrep; rirrep = leftirrep;}
    if (opdata->get_deltaQuantum().get_s() == 2)
    {
      if (lS == rS)
	scale = 1.0;
      else if (lS == rS+2)
	scale = -1.0*sqrt((rS+3.0)/(rS+1.0));
      else if (rS == lS+2)
	scale = -1.0*sqrt((lS+1.0)/(lS+3.0));
      else
	{cout<< "inappropriate left and right spin states"<<endl;exit(0);}
    }
    else if (opdata->get_deltaQuantum().get_s() == 0)
      scale = 1.0;
    else if (opdata->get_deltaQuantum().get_s() == 1)
    {
      if (lS+1 == rS)
	scale = -1.0*sqrt((lS+1.0)/(lS+2.0));
      else if (lS == rS+1)
	scale = sqrt((rS+2.0)/(rS+1.0));
      else
	{cout<< "inappropriate left and right spin states "<<endl;exit(0);}
    }
    if (Symmetry::sizeofIrrep(opdata->get_deltaQuantum().get_symm().getirrep()) != 1) {
      if (Symmetry::sizeofIrrep(rirrep) == 1)
	scale = scale* sqrt(2.0);
      else if (Symmetry::sizeofIrrep(lirrep) == 1)
	scale = scale/sqrt(2.0);
    }
    return scale;
  }
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) {return opdata;}
  void build(const SpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b){return 0.0;}
}; 

const Transposeview Transpose(SparseMatrix& op);

double getCommuteParity(SpinQuantum a, SpinQuantum b, SpinQuantum c);
void Normalise(SparseMatrix& a, int* success = 0);
void ScaleAdd(double d, const SparseMatrix& a, SparseMatrix& b);
double DotProduct(const SparseMatrix& lhs, const SparseMatrix& rhs);
void Scale(double d, SparseMatrix& a);
void assignloopblock(SpinBlock*& loopblock, SpinBlock*& otherblock, SpinBlock* leftSpinBlock, SpinBlock* rightSpinBlock);
void copy(const ObjectMatrix<Matrix>& a, ObjectMatrix<Matrix>& b);
void copy(const Matrix& a, Matrix& b);
}

#endif
