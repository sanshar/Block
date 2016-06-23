#ifndef FCIQMC_HELPER_HEADER_H
#define FCIQMC_HELPER_HEADER_H

#include "global.h"
#include <vector>
#include "boost/shared_ptr.hpp"
#include "StateInfo.h"
#include "Operators.h"
#include "wavefunction.h"
#include "ObjectMatrix.h"

typedef unsigned long int ulong;

namespace SpinAdapted{


class QSTensor {
private:
  friend class boost::serialization::access;
  template<class Archive> 
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & allowedQuanta & leftQuanta & rightQuanta & data & nl & nr;
    }  
  ObjectMatrix<char> allowedQuanta;
  std::vector<SpinQuantum> leftQuanta, rightQuanta;
  ObjectMatrix<Matrix> data;
  int nl, nr;
public:
  QSTensor():nl(0), nr(0) {}
  QSTensor(int l, int r):nl(l), nr(r) {allowedQuanta.ReSize(nl, nr); data.ReSize(nl, nr);}
  QSTensor(const std::vector<SpinQuantum>& lq, const std::vector<SpinQuantum>& rq):leftQuanta(lq), rightQuanta(rq), nl(lq.size()), nr(rq.size()) {
    allowedQuanta.ReSize(nl, nr); data.ReSize(nl, nr);
  }
  void Init(int l, int r) {
    nl = l;
    nr = r;
    allowedQuanta.ReSize(nl, nr); data.ReSize(nl, nr);
  }
  void Init(const std::vector<SpinQuantum>& lq, const std::vector<SpinQuantum>& rq) {
    leftQuanta = lq;
    rightQuanta = rq;
    nl = lq.size();
    nr = rq.size();
    allowedQuanta.ReSize(nl, nr); data.ReSize(nl, nr);
  }
  bool allowed(int i, int j) const {return allowedQuanta(i,j);}
  ObjectMatrix<char>& allowedMatrix() {return allowedQuanta;}
  Matrix& operator() (int i, int j) {return data(i,j);}
  const Matrix& operator() (int i, int j) const {return data(i,j);}
  std::vector<SpinQuantum>& lQuanta() {return leftQuanta;}
  const std::vector<SpinQuantum>& lQuanta() const {return leftQuanta;}
  std::vector<SpinQuantum>& rQuanta() {return rightQuanta;}
  const std::vector<SpinQuantum>& rQuanta() const {return rightQuanta;}
  int lsize() const { return nl;}
  int rsize() const { return nr;}
  void CleanUp() {
    data.ReSize(0,0);
    allowedQuanta.ReSize(0,0);
    leftQuanta.clear();
    rightQuanta.clear();
    nl = 0;
    nr = 0;
  }
  void remove_empty();
  friend ostream& operator<< (ostream& os, const QSTensor& q);
};

QSTensor TensorProduct(const QSTensor& A, const QSTensor& B);

void SaveQSTensor(const int& site, const QSTensor& m, int state);
void LoadQSTensor(const int& site, QSTensor& m, int state);

//the MPS is stored in the left canonical form
//LLLLL..LC 

class MPS{

 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & SiteTensors \
	& w & MPSrep;
    }

  std::vector< std::vector<Matrix> > SiteTensors; //these are the L matrices
  std::vector<std::vector<QSTensor>> MPSrep;
  Wavefunction w; //the last wavefunction

  void Init(std::vector<bool>& occ);
 public:
  static int sweepIters;
  static bool spinAdapted;
  static std::vector<SpinBlock> siteBlocks;
  // It is the site block with implicit transpose.

  MPS() {};
  MPS(int stateindex); 
  MPS(std::vector<bool>& occ);
  MPS(ulong* occnum, int length);
  void buildMPSrep();
  std::vector<Matrix>& getSiteTensors(int i) {return SiteTensors[i];}
  const std::vector<Matrix>& getSiteTensors(int i) const {return SiteTensors[i];}
  const Wavefunction& getw() const {return w;}
  void scale(double r) {Scale(r, w);}
  void normalize() {int success; w.Normalise(&success);}
  double get_coefficient(const vector<bool>& occ_strings);
  void writeToDiskForDMRG(int state, bool writeStateAverage=false);
};


 //statea is multiplied with Operator O|Mpsa> and then we compress it to get stateb
 //void compressOperatorTimesMPS(const MPS& statea, MPS& stateb);

 //calculate overlap between a and b <Mpsa|Mpsb>
 double calculateOverlap (const MPS& a, const MPS& b);

 //calculate hamiltonian matrix between a and b <Mpsa|H|Mpsb>
 void calcHamiltonianAndOverlap(const MPS& statea, const MPS& stateb, double& h, double& o, bool sameStates=false) ;

}

#endif
