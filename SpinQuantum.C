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
SpinQuantum::SpinQuantum (const int p, const SpinSpace s, const IrrepSpace orbS) : particleNumber (p), totalSpin (s), orbitalSymmetry(orbS) {}

int SpinQuantum::insertionNum(const SpinQuantum& ql, const SpinQuantum& qr) const
{
  int index = 0;
  vector<SpinSpace> spins = ql.totalSpin+qr.totalSpin;
  vector<IrrepSpace> vec = ql.orbitalSymmetry+qr.orbitalSymmetry;
  for (int i=0; i< spins.size(); i++) {
    for (int j=0; j<vec.size(); j++) {
      if (*this == SpinQuantum(ql.particleNumber+qr.particleNumber, spins[i], vec[j]))
	return index;
      else 
	index++;
    }
  }
  return -1;
  
}

vector<SpinQuantum> SpinQuantum::spinToNonSpin() const
{
  vector<SpinQuantum> quanta;
  if(totalSpin.getirrep() < 0 ) {
    pout << "The quanta is already nonspinAdapted"<<endl;
    exit(0);
  }
  if (dmrginp.spinAdapted()) {
    quanta.push_back(*this);
    return quanta;
  }
  else {
    for (int i=-totalSpin.getirrep(); i<=totalSpin.getirrep(); i+=2)
      quanta.push_back(SpinQuantum(particleNumber, SpinSpace(i), orbitalSymmetry));
    return quanta;
  }
}

vector<SpinQuantum> SpinQuantum::operator+ (const SpinQuantum& q) const
{
  vector<SpinQuantum> quanta;
  vector<SpinSpace> spins = totalSpin+q.totalSpin;
  vector<IrrepSpace> vec = orbitalSymmetry+q.orbitalSymmetry;
  for (int i=0; i< spins.size(); i++) 
  for (int j=0; j<vec.size(); j++) 
    quanta.push_back(SpinQuantum(particleNumber+q.particleNumber, spins[i], vec[j]));
  
  
  return quanta;
}

vector<SpinQuantum> SpinQuantum::operator+ (const vector<SpinQuantum>& q) const
{
  vector<SpinQuantum> quanta;
  for (int i=0; i< q.size(); i++){
    vector<SpinQuantum> tmp=*this+q[i] ;
    quanta.reserve(quanta.size()+tmp.size());
    quanta.insert(quanta.end(),tmp.begin(),tmp.end());
  }
  return quanta;
}

vector<SpinQuantum> SpinQuantum::operator- (const SpinQuantum& q) const
{
  SpinQuantum negative = -q;
  return *this+negative;
}

vector<SpinQuantum> SpinQuantum::operator- (const vector<SpinQuantum>& q) const
{
  vector<SpinQuantum> quanta;
  for (int i=0; i< q.size(); i++){
    vector<SpinQuantum> tmp=*this-q[i] ;
    quanta.reserve(quanta.size()+tmp.size());
    quanta.insert(quanta.end(),tmp.begin(),tmp.end());
  }
  return quanta;
}

SpinQuantum SpinQuantum::operator-() const
{
  return SpinQuantum(-particleNumber, -totalSpin, -orbitalSymmetry);
}

bool SpinQuantum::allow(const SpinQuantum s1, const SpinQuantum s2) const
{
  std::vector<SpinQuantum> sumQ = s1+s2;
  for (int i=0; i<sumQ.size(); i++)
    if (*this == sumQ[i])
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
  int s_try = abs(dmrginp.total_spin_number().getirrep() - q.get_s().getirrep());
  return (abs(s_try) <= q_try);
}

vector<SpinQuantum> SpinQuantum::get_complement () const
{
  vector<SpinQuantum> quanta;
  vector<int> ns;
  if (dmrginp.hamiltonian() == BCS) {
    for (int n = 0; n <= dmrginp.total_particle_number(); n+=2) {
      if (n-particleNumber < 0) continue;
      ns.push_back(n-particleNumber);
    }
  } else {
    ns.push_back(dmrginp.total_particle_number() - particleNumber);
  }
  vector<SpinSpace> spins = dmrginp.total_spin_number() + (- totalSpin);
  vector<IrrepSpace> vec = dmrginp.total_symmetry_number() + (- orbitalSymmetry);
  for (int n_idx = 0; n_idx < ns.size(); ++n_idx) {
    int n = ns[n_idx];
    for (int i=0; i< spins.size(); i++) {
      if (abs(i) <= n) {
        for (int j=0; j<vec.size(); j++)
      quanta.push_back(SpinQuantum(n, spins[i], vec[j]));
      }
    }
  }
  return quanta;
}

}
