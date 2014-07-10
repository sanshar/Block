/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_CSF_HEADER
#define SPIN_CSF_HEADER
#include "orbstring.h"
#include <iostream>
#include "global.h"
#include <map>
#include "slater.h"
#include "Symmetry.h"
#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include "IrrepVector.h"
#include "tensor_operator.h"

namespace SpinAdapted{
struct Csf
{
  map<Slater, double> det_rep;
  
  int n;
  SpinSpace S;
  int Sz;
  IrrepVector irrep;

public:
  inline Csf () {}
  Csf( const map<Slater, double>& p_dets, const int p_n, const SpinSpace p_S, const int p_Sz, const IrrepVector pirrep);

  Csf (const Csf& s) : det_rep(s.det_rep), n(s.n), S(s.S), Sz(s.Sz), irrep(s.irrep){} 
  //void operator= (const Csf& s) { if (this != &s) {det_rep =s.det_rep; n = s.n; S = s.S; Sz = s.Sz; sym=s.sym;}}

  // accessors
  inline int size () const { return det_rep.size (); }
  inline int n_is () const { return n; }
  inline SpinSpace S_is () const { return S; }
  inline int Sz_is () const { return Sz; }
  inline int row() const {return irrep.getrow();}
  inline IrrepSpace sym_is() const {
    return IrrepSpace(irrep.getirrep());
  }
  
  void set_det_rep(map<Slater, double> p_det, SpinSpace pS, IrrepVector pirrep){  
    det_rep = p_det;
    S = pS;
    map<Slater, double>::iterator it = det_rep.begin();
    n = it->first.n_is();
    Sz = it->first.Sz_is();
    irrep = pirrep;
  }

  void set_n(int p_n){n = p_n;}
  void set_S(SpinSpace p_S){S = p_S;}
  void set_Sz(int p_Sz){Sz = p_Sz;}
  void set_irrep(IrrepVector p_irrep) {irrep = p_irrep;}

  inline map<Slater, double> c (int i)
  { 
    Slater s;
    double d;
    map<Slater, double> dets;
    for (map<Slater, double>::iterator it = det_rep.begin(); it!= det_rep.end(); it++) {
      s = (*it).first;
      d = (*it).second;
      s.c(i);
      dets[s] = d;
    }
    return dets;
  }
  
  inline map<Slater, double> d (int i)
  {
    Slater s;
    double d;
    map<Slater, double> dets;
    for (map<Slater, double>::iterator it = det_rep.begin(); it!= det_rep.end(); it++) {
      s = (*it).first;
      d = (*it).second;
      s.d(i);
      dets[s] = d;
    }
    return dets;
  }

  void applySplus(Csf& output);
  void applySminus(Csf& output);
  void applyRowminus(Csf& output, int orbL);
  vector<Csf> spinLadder(int k);
  void outerProd(const Csf& csf, double factor, map<Slater, double>& output) const;

  void normalize()
  {
    double norm=0;
    map<Slater, double>::iterator it = det_rep.begin();
    for (; it!= det_rep.end(); it++) 
      norm += pow((*it).second,2);
    for (it = det_rep.begin(); it!= det_rep.end(); it++) 
      (*it).second /= sqrt(norm);
  }

  double norm()
  {
    double norm=0;
    map<Slater, double >::iterator it = det_rep.begin();
    for (; it!= det_rep.end(); it++)
      norm += pow((*it).second,2);
    return norm;
  }

     
  bool isempty()
  {
    for (map<Slater, double>::iterator it = det_rep.begin(); it!= det_rep.end(); it++)
      if (!(*it).first.alpha.isempty())
	return false;
    return true;
  }

  bool operator== (const Csf& s) const
  {
    if (SpinQuantum(n, S, sym_is()) == SpinQuantum(s.n, s.S, s.sym_is()))
      return true;
    else
      return false;
  }



  bool operator< (const Csf& s) const;
  friend ostream& operator<< (ostream& os, const Csf& s)
  {
    os <<"n: "<<s.n<<" S: "<<s.S<<" Sz: "<<s.Sz<<"  Irrep: "<<s.irrep<<","<<s.row()<<endl;
    map<Slater, double>::const_iterator it = s.det_rep.begin();    
    for (; it!= s.det_rep.end(); it++)
      os<<(*it).second<<" "<<(*it).first;
    return os;
  }

  static std::vector< Csf > distribute (const int n, const int s, const IrrepVector &sym, const int left, const int right, const int edge, int integralIndex);
  static std::vector<Csf> distributeNonSpinAdapted (const int n, const int sp, const IrrepVector &sym, const int left, const int right, const int edge, int integralIndex);
    
};


namespace CSFUTIL {
  std::vector< Csf > spinfockstrings(const std::vector<int>& orbs, std::vector<std::vector<Csf> >& ladders);  
  std::vector< Csf > spinfockstrings(const std::vector<int>& orbs);  
  void TensorProduct(Csf& lhs, vector<Csf>& lhs_csfs, Csf& rhs, vector<Csf>& rhs_csfs, vector< Csf >& output, vector< vector<Csf> >& outputladder);
  void TensorProduct(Csf& lhs, Csf& rhs, vector< Csf >& output);
  Csf applyTensorOp(const TensorOp& newop, int spinL);
  vector< vector<int> > generate_partitions(int k);
}

struct Csfcompare
{
  Csfcompare(){};
  bool operator() (const boost::shared_ptr<Csf>& c1, const boost::shared_ptr<Csf>& c2) const {return (*c1 < *c2);}
};

 double csf_energy (const Csf& s, int integralIndex);
}
#endif
