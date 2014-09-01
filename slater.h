/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SLATER_HEADER
#define SPIN_SLATER_HEADER
#include "orbstring.h"
#include <iostream>
#include "global.h"
#include <map>
#include "Symmetry.h"
#include <boost/functional/hash.hpp>


namespace SpinAdapted{
  //this is the slater orbital index
inline int SymmetryOfSpinOrb (const int i)
{
  return dmrginp.spin_orbs_symmetry()[i];
}

  //this is the irrep orbital index
inline IrrepSpace SymmetryOfSpatialOrb(const int i)
{
  int symm = dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[i]];
  if (sym == "dinfh")
    symm = abs(symm);

  return IrrepSpace(symm);
}

inline IrrepSpace SymmetryOfOrb(const int i)
{
  if(dmrginp.spinAdapted()) return SymmetryOfSpatialOrb(i);

  int symm = dmrginp.spin_orbs_symmetry()[i];
  if (sym == "dinfh")
    symm = abs(symm);

  return IrrepSpace(symm);
}

inline int SpinOf (const int i)
{
  return i%2==0?1:-1;
}
inline int SzOf(const int i)
{
  return i%2==0? 1 : -1;
}

inline SpinQuantum getSpinQuantum(const int i)
{
  if(dmrginp.spinAdapted())
    return SpinQuantum(1, SpinSpace(1), SymmetryOfSpatialOrb(i));
  else
    return SpinQuantum(1, SpinSpace(SpinOf(i)), IrrepSpace(abs(dmrginp.spin_orbs_symmetry()[i])));
}

void ConvertList (std::vector<int>& a, const std::vector<int>& b);

struct Slater
{
  Orbstring alpha;

  int n;
  int Sz;

public:
  Slater (): n(0), Sz(0) {}
  Slater (const Orbstring& a);
  Slater (const vector<bool>& occ_rep, int sign);
  Slater (const Slater& s) : alpha (s.alpha), n (s.n), Sz (s.Sz) {}
  void operator= (const Slater& s) { if (this != &s) {alpha = s.alpha; n = s.n; Sz = s.Sz;} }

  // accessors
  inline int size () const { return alpha.size (); }
  inline int n_is () const { return n; }
  inline int Sz_is () const { return Sz; }
  inline Slater& c (int i)
  { 
    this->alpha.c (i);
    ++n;
    Sz += SzOf (i);
    return (*this);
  }
  
  boost::shared_ptr<Slater> getLeftSlater(int index);
  boost::shared_ptr<Slater> getRightSlater(int index);

  inline Slater& d (int i)
  {
    this->alpha.d (i);
    --n;
    Sz -= SzOf (i);
    return (*this);
  }

  bool isempty()
  {
    return alpha.isempty();
  }
  inline int parity(int i)
  {
    return alpha.parity(i);
  }

  inline void setSign(int i)
  {
    alpha.setSign(i);
  }

  inline int getSign()
  {
    return alpha.getSign();
  }

  inline int trace (const Slater& s)
  {
    return this->alpha.trace (s.alpha);
  }
  void connect (const Slater& s, std::vector<int>& cv, std::vector<int>& dv) const;

  bool operator< (const Slater& s) const;
  bool operator== (const Slater& s) const ;
  double hash_value() const;

  void outerProd(const Slater& s, Slater& output) const;

  friend ostream& operator<< (ostream& os, const Slater& s)
  {
    os << s.alpha << endl;
    return os;
  }
};


inline IrrepSpace AbelianSymmetryOf (const Slater& s)
{
  if(sym == "dinfh") {perr <<"dinfh is not abelian"<<endl;exit(0);}
  IrrepSpace symmetry = IrrepSpace(0), tempsym;
  for (int i=0; i<s.alpha.size(); ++i) {
    if (s.alpha[i]) {
      tempsym = (symmetry + IrrepSpace(SymmetryOfSpinOrb(i)))[0];
      symmetry = tempsym;
    }
  }
  return symmetry;
}

extern void Convert (Slater& s);


inline int flipindex(int i)
{
  return 0;
}
double det_energy (const Slater& s, int integralIndex);
}
#endif
