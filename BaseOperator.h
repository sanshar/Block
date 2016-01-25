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
//FIXME
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
//<<
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

class TensorOp;
namespace SpinAdapted{
class SpinBlock;

//when bra and ket states are different then DES cannot be made by just taking the transpose of
//CRE. In such a case we need to explicitly make the operators in the second row. 
enum opTypes{ HAM, CRE, CRE_CRE, DES_DESCOMP, CRE_DES, CRE_DESCOMP, CRE_CRE_DESCOMP, 
              DES, DES_DES, CRE_CRECOMP, DES_CRE, DES_CRECOMP, CRE_DES_DESCOMP, OVERLAP,
              // Extra (empty) operator classes for RI-approx NPDM
              RI_3INDEX, RI_4INDEX,
              // Extra 3PDM operators
              CRE_CRE_CRE, CRE_CRE_DES, CRE_DES_DES, CRE_DES_CRE,
              // Extra 4PDM operators
              DES_CRE_DES, DES_DES_CRE, DES_CRE_CRE, DES_DES_DES,
              CRE_CRE_DES_DES, CRE_DES_CRE_DES, CRE_DES_DES_CRE, CRE_DES_DES_DES,
              CRE_CRE_CRE_DES, CRE_CRE_DES_CRE, CRE_DES_CRE_CRE, CRE_CRE_CRE_CRE ,CDD_SUM, CDD_CRE_DESCOMP, CDD_DES_DESCOMP, CCD_SUM, CCD_CRE_DESCOMP, CCD_CRE_CRECOMP};

enum CompType{CD, DD, CCD, C, CDD, CDD_CD, CCD_CD};

template<class T> class Baseoperator  // The abstract class of an operator
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
  virtual int get_deltaQuantum_size() const = 0;
  virtual std::vector<SpinQuantum> get_deltaQuantum() const = 0;
  virtual SpinQuantum get_deltaQuantum(int i) const = 0;  
  virtual char conjugacy() const = 0;
  virtual ~Baseoperator() {};
  virtual SpinSpace get_spin(int i=0) const = 0;
  virtual IrrepSpace get_symm(int i=0) const = 0;
  virtual double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const = 0;
  Baseoperator() {};
};


class SparseMatrix : public Baseoperator<Matrix>  // the sparse matrix representation of the operator
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & orbs \
         & deltaQuantum \
         & quantum_ladder \
         & build_pattern \
         & fermion \
         & initialised \
         & built \
         & built_on_disk \
         & allowedQuantaMatrix \
         & quantum_ladder \
         & build_pattern \
         & operatorMatrix \
         & Sign;
    }

 protected:
  std::vector<int> orbs;
  bool fermion;
  ObjectMatrix<char> allowedQuantaMatrix;  // some whether it is allowed?
  bool initialised;
  bool built;
  bool built_on_disk;
  ObjectMatrix<Matrix> operatorMatrix;  // put the dense block in the place it should be
  std::vector<SpinQuantum> deltaQuantum;    // allowed quantum
  int Sign;
  // With N-index ops (N>2), there are several ways to build them...
  std::string build_pattern;
  // ...and for each way we record the spin ladder components
  std::map< std::string, std::vector<SpinQuantum> > quantum_ladder;

 public:
  SparseMatrix() : fermion(false), orbs(2), initialised(false), built(false), built_on_disk(false), Sign(1){};
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
  std::vector<SpinQuantum> &set_deltaQuantum() { return deltaQuantum; }
  SpinQuantum &set_deltaQuantum(int i) { return deltaQuantum[i]; }
  void resize_deltaQuantum(int i) { deltaQuantum.resize(i); }
  void set_deltaQuantum(int i, const SpinQuantum s) { deltaQuantum.assign(i, s); }
  int get_deltaQuantum_size() const { return deltaQuantum.size(); }
  SpinSpace get_spin(int i=0) const  { return deltaQuantum[i].get_s();}
  IrrepSpace get_symm(int i=0) const  { return deltaQuantum[i].get_symm();}
  std::vector<SpinQuantum> get_deltaQuantum() const { return deltaQuantum; }
  std::map< std::string, std::vector<SpinQuantum> >  get_quantum_ladder() const { return quantum_ladder; }
  std::map< std::string, std::vector<SpinQuantum> >& set_quantum_ladder() { return quantum_ladder; }
  std::string  get_build_pattern() const { return build_pattern; }
  std::string& set_build_pattern() { return build_pattern; }
  SpinQuantum get_deltaQuantum(int i) const { return deltaQuantum[i]; }
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
  const bool& get_built_on_disk() const { return built_on_disk; }
  bool& set_built_on_disk() { return built_on_disk; }
  double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const {return 1.0;}

  void resize(int n, int c) { operatorMatrix.ReSize(n, c); allowedQuantaMatrix.ReSize(n, c); }
  const Matrix& operator_element(int i, int j) const { return operatorMatrix(i, j); }
  const Matrix& operator()(int i, int j) const { return operatorMatrix(i, j); }
  Matrix& operator_element(int i, int j) { return operatorMatrix(i, j); }
  Matrix& operator()(int i, int j) { return operatorMatrix(i, j); }
  void allocate(const StateInfo& s);
  void allocate(const StateInfo& sr, const StateInfo& sc);
  void allocate(const SpinBlock& b);
  void makeIdentity(const StateInfo& s); 
  virtual boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) =0;
  virtual void build(const SpinBlock& b) =0;
  void buildUsingCsf(const SpinBlock& b, vector< vector<Csf> >& ladders, std::vector< Csf >& s) ;
  virtual double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b=0)=0;
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const TwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const PerturbTwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const OneElectronArray& v_1, int integralIndex);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const CCCCArray& vcccc);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const CCCDArray& vcccd);  
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const PerturbTwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const CCCCArray& vcccc);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const CCCDArray& vcccd);
  bool nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb);
  std::vector<double> calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2);
  friend ostream& operator<<(ostream& os, const SparseMatrix& a);
  void Randomise();
  void SymmetricRandomise();
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
  void renormalise_transform(const std::vector<Matrix>& leftrotate_matrix, const StateInfo* leftstateinfo, const std::vector<Matrix>& rightrotate_matrix, const StateInfo *rightstateinfo);
  void build_and_renormalise_transform(SpinBlock *big, const opTypes &ot, const std::vector<Matrix>& leftrotate_matrix, const StateInfo *leftstateinfo, 
				       const std::vector<Matrix>& rightrotate_matrix,  const StateInfo *newStateInfo);
  SparseMatrix& operator+=(const SparseMatrix& other);
};

class Transposeview : public SparseMatrix
{
private:
  boost::shared_ptr<SparseMatrix> opdata; 
public:
  Transposeview(const boost::shared_ptr<SparseMatrix>& opptr) : opdata(opptr) {}
  Transposeview(SparseMatrix& op) { opdata = boost::shared_ptr<SparseMatrix>(&op, boostutils::null_deleter());}
  int get_deltaQuantum_size() const { return opdata->get_deltaQuantum_size(); }  
  SpinQuantum get_deltaQuantum(int i) const {return -opdata->get_deltaQuantum(i);}
  std::vector<SpinQuantum> get_deltaQuantum() const {
    std::vector<SpinQuantum> q;
    for (int i = 0; i < opdata->get_deltaQuantum_size(); ++i) {
      q.push_back(-opdata->get_deltaQuantum(i));
    }
    return q;
  }
  bool get_fermion() const { return opdata->get_fermion(); }
  bool get_initialised() const { return opdata->get_initialised(); }
  int nrows() const { return opdata->ncols(); }
  int ncols() const { return opdata->nrows(); }
  const char &allowed(int i, int j) const { return opdata->allowed(j, i); }
  char &allowed(int i, int j) { return opdata->allowed(j, i); }
  const Matrix& operator_element(int i, int j) const { return opdata->operator_element(j, i); }
  Matrix& operator_element(int i, int j) { return opdata->operator_element(j, i); }
  SpinSpace get_spin(int i=0) const  { return opdata->get_deltaQuantum(i).get_s();}
  IrrepSpace get_symm(int i=0) const  { return -opdata->get_deltaQuantum(i).get_symm();}
  int get_orbs(int i) const {return opdata->get_orbs(i);}
  const std::vector<int>& get_orbs() const { return opdata->get_orbs(); }
  const Matrix& operator()(int i, int j) const { return opdata->operator()(j, i); }
  Matrix& operator()(int i, int j) { return opdata->operator()(j, i); }
  char conjugacy() const { if (opdata->conjugacy() == 'n') return 't'; else return 'n';}
  double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const ;
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) {return opdata;}
  void build(const SpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b){return 0.0;}
}; 

class SubSparseMatrix : public SparseMatrix
{
private:
  boost::shared_ptr<SparseMatrix> opdata;
  int section;
  ObjectMatrix<char> SuballowedQuantaMatrix;
public:
  SubSparseMatrix(SparseMatrix& op, int sec, const StateInfo& s);  
  SubSparseMatrix(const boost::shared_ptr<SparseMatrix>& opptr, int sec, const StateInfo& s);
  SubSparseMatrix(const boost::shared_ptr<SparseMatrix>& opptr, int sec, const StateInfo& sr, const StateInfo& sc);
  std::vector<SpinQuantum> get_deltaQuantum() const {
    std::vector<SpinQuantum> deltaQuantum(1, opdata->get_deltaQuantum(section));
    return deltaQuantum;
  }
  int get_deltaQuantum_size() const {return 1;}
  SpinQuantum get_deltaQuantum(int i) const {
    return get_deltaQuantum()[0]; // this i is dummy
  }
  bool get_fermion() const { return opdata->get_fermion(); }
  bool get_initialised() const { return opdata->get_initialised(); }
  int nrows() const { return opdata->nrows(); }
  int ncols() const { return opdata->ncols(); }
  const char &allowed(int i, int j) const { return SuballowedQuantaMatrix(i, j); }  
  char &allowed(int i, int j) { return SuballowedQuantaMatrix(i, j); }
  const Matrix& operator_element(int i, int j) const { return opdata->operator_element(i, j); }
  Matrix& operator_element(int i, int j) { return opdata->operator_element(i, j); }
  SpinSpace get_spin(int i=0) const  { return opdata->get_deltaQuantum(section).get_s();}
  IrrepSpace get_symm(int i=0) const  { return -opdata->get_deltaQuantum(section).get_symm();}
  int get_orbs(int i) const {return opdata->get_orbs(i);}
  const std::vector<int>& get_orbs() const { return opdata->get_orbs(); }
  const Matrix& operator()(int i, int j) const { return opdata->operator()(i, j); }
  Matrix& operator()(int i, int j) { return opdata->operator()(i, j); }
  char conjugacy() const { return opdata->conjugacy(); }
  double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const {  return 1.0; }
  boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) {return opdata;}
  void build(const SpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b){return 0.0;}
};

const Transposeview Transpose(SparseMatrix& op);

double getCommuteParity(SpinQuantum a, SpinQuantum b, SpinQuantum c);
void Normalise(SparseMatrix& a, int* success = 0);
void ScaleAdd(double d, const SparseMatrix& a, SparseMatrix& b);
double DotProduct(const SparseMatrix& lhs, const SparseMatrix& rhs);
double trace(const SparseMatrix& lhs);
void Scale(double d, SparseMatrix& a);
void assignloopblock(SpinBlock*& loopblock, SpinBlock*& otherblock, SpinBlock* leftSpinBlock, SpinBlock* rightSpinBlock);
void copy(const ObjectMatrix<Matrix>& a, ObjectMatrix<Matrix>& b);
void copy(const Matrix& a, Matrix& b);
} 


#endif
