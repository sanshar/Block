/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "SpinQuantum.h"
#include "global.h"

namespace SpinAdapted{

SpinQuantum::SpinQuantum () : particleNumber (0), totalSpin (0), orbitalSymmetry (0) {}
SpinQuantum::SpinQuantum (const int p, const int s, const IrrepSpace orbS) : particleNumber (p), totalSpin (s), orbitalSymmetry(orbS) {}

int SpinQuantum::insertionNum(const SpinQuantum& ql, const SpinQuantum& qr) const
{
  int index = 0;
  for (int i=abs(ql.totalSpin- qr.totalSpin); i<= ql.totalSpin+qr.totalSpin; i+=2) {
    vector<IrrepSpace> vec = ql.orbitalSymmetry+qr.orbitalSymmetry;
    for (int j=0; j<vec.size(); j++) {
      if (*this == SpinQuantum(ql.particleNumber+qr.particleNumber, i, vec[j]))
	return index;
      else 
	index++;
    }
  }
  return -1;
  
}

vector<SpinQuantum> SpinQuantum::operator+ (const SpinQuantum q) const
{
  vector<SpinQuantum> quanta;
  for (int i=abs(q.totalSpin- totalSpin); i<= q.totalSpin+totalSpin; i+=2) {
    vector<IrrepSpace> vec = orbitalSymmetry+q.orbitalSymmetry;
    for (int j=0; j<vec.size(); j++) {
      quanta.push_back(SpinQuantum(particleNumber+q.particleNumber, i, vec[j]));
    }
  }
  return quanta;
}

vector<SpinQuantum> SpinQuantum::operator- (const SpinQuantum q) const
{
  vector<SpinQuantum> quanta;
  for (int i=abs(q.totalSpin- totalSpin); i<= q.totalSpin+totalSpin; i+=2) {
    IrrepSpace negativeQSym = -q.orbitalSymmetry;
    vector<IrrepSpace> vec = orbitalSymmetry+negativeQSym;
    for (int j=0; j<vec.size(); j++)
      quanta.push_back(SpinQuantum(particleNumber-q.particleNumber, i, vec[j]));
  }
  return quanta;
}

SpinQuantum SpinQuantum::operator-() const
{
  return SpinQuantum(-particleNumber, totalSpin, -orbitalSymmetry);
}

bool SpinQuantum::allow(const SpinQuantum s1, const SpinQuantum s2) const
{
  if (particleNumber != s1.get_n() + s2.get_n()) return false;
  if ( totalSpin < abs(s1.get_s() - s2.get_s()) || totalSpin > (s1.get_s() + s2.get_s())) return false;
  if ( (totalSpin - abs(s1.get_s() - s2.get_s()))%2 != 0) return false;

  IrrepSpace sym1 = s1.get_symm(), sym2 = s2.get_symm();
  std::vector<IrrepSpace> sym3 = sym1+sym2;

  for (int i=0 ;i<sym3.size(); i++)
    if (sym3[i] == get_symm())
      return true;

  return false;
}
  
bool SpinQuantum::operator== (const SpinQuantum q) const
{
  return ((particleNumber == q.particleNumber) && (totalSpin == q.totalSpin) && (orbitalSymmetry == q.orbitalSymmetry));
}

bool SpinQuantum::operator!= (const SpinQuantum q) const
{
  return ((particleNumber != q.particleNumber) || (totalSpin != q.totalSpin) || (orbitalSymmetry != q.orbitalSymmetry));
}

bool SpinQuantum::operator< (const SpinQuantum q) const
{
  return ((particleNumber < q.particleNumber) || ((particleNumber == q.particleNumber) && (totalSpin < q.totalSpin)) ||
	  ((particleNumber == q.particleNumber) && (totalSpin == q.totalSpin) && (orbitalSymmetry < q.orbitalSymmetry)));
}

bool IsFermion (const SpinQuantum q)
{
  return (q.particleNumber & 1);
}
ostream& operator<< (ostream& os, const SpinQuantum q)
{
  os << q.particleNumber << ":" << q.totalSpin << ":" << q.orbitalSymmetry;
  return os;
}

bool SpinQuantum::can_complement (SpinQuantum q)
{
  int q_try = dmrginp.total_particle_number() - q.get_n();
  int s_try = abs(dmrginp.total_spin_number() - q.get_s());
  return (abs (s_try) <= q_try);
}

void SpinQuantum::complementize ()
{
  particleNumber = dmrginp.total_particle_number() - particleNumber;
  totalSpin = abs(dmrginp.total_spin_number() - totalSpin);
  orbitalSymmetry = (dmrginp.total_symmetry_number() + orbitalSymmetry)[0];
}

vector<SpinQuantum> SpinQuantum::get_complement () const
{
  vector<SpinQuantum> quanta;
  int n = dmrginp.total_particle_number() - particleNumber;
  vector<IrrepSpace> vec = dmrginp.total_symmetry_number() + orbitalSymmetry;
  for (int i=abs(dmrginp.total_spin_number()- totalSpin); i<= dmrginp.total_spin_number()+totalSpin; i+=2) {
    if (i <= n) {
      for (int j=0; j<vec.size(); j++)
	quanta.push_back(SpinQuantum(n, i, vec[j]));
    }
  }
  return quanta;
}

}
