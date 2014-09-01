/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "slater.h"
#include "global.h"

bool SpinAdapted::Slater::operator< (const Slater& s) const
{
  if (n < s.n) return true;
  for (int i=alpha.size()-1; i>=0; i--)
  {
    if (alpha[i] == s.alpha[i])
      continue;
    else if(alpha[i] > s.alpha[i])
      return false;
    else if (alpha[i] < s.alpha[i])
      return true;
  }
  return false;
}

double SpinAdapted::Slater::hash_value() const
{
  double j = 0.0;
  boost::hash<double> hasher;
  if (alpha.isempty())
    return hasher(j);
  for (int i=0; i<alpha.size(); i++)
    j += alpha[i]*pow(2.0,i);
  return j;
}

bool SpinAdapted::Slater::operator==(const Slater& s2) const
{
  if (alpha.isempty() && s2.alpha.isempty())
    return true;
  if ( (!alpha.isempty() && s2.alpha.isempty()) || (alpha.isempty() && !s2.alpha.isempty()))
    return false;

  if ((alpha.get_occ_rep().size() != s2.alpha.get_occ_rep().size()) ) return false;
  for (int i=0; i<alpha.get_occ_rep().size(); i++) 
    if (alpha[i] != s2.alpha[i])
      return false;
  return true;
}

SpinAdapted::Slater::Slater (const vector<bool>& occrep, int sign):alpha(occrep, sign)
{
  
  n = 0;
#ifdef DEC  
  count (alpha.get_occ_rep().begin (), alpha.get_occ_rep().end (), 1, n);
#else
  n = count (alpha.get_occ_rep().begin (), alpha.get_occ_rep().end (), 1);
#endif
  Sz = 0;
  for (int i = 0; i < alpha.get_occ_rep().size (); ++i)
    if (alpha.get_occ_rep() [i]) Sz += SzOf(i);

}


SpinAdapted::Slater::Slater (const Orbstring& a)
{
  alpha = a;
  n = 0;
#ifdef DEC  
  count (a.get_occ_rep().begin (), a.get_occ_rep().end (), 1, n);
#else
  n = count (a.get_occ_rep().begin (), a.get_occ_rep().end (), 1);
#endif
  Sz = 0;
  for (int i = 0; i < a.get_occ_rep().size (); ++i)
    if (a.get_occ_rep() [i]) Sz += SzOf(i);

}

boost::shared_ptr<SpinAdapted::Slater> SpinAdapted::Slater::getLeftSlater (int index)
{
  std::vector<bool>::const_iterator it = alpha.get_occ_rep().begin();
  std::vector<bool> occ_rep_left;
  occ_rep_left.insert(occ_rep_left.begin(), it, it+index);
  Orbstring o(occ_rep_left);
  boost::shared_ptr<Slater> leftSlater(new Slater(o));
  return leftSlater;
}

boost::shared_ptr<SpinAdapted::Slater> SpinAdapted::Slater::getRightSlater (int index)
{
  std::vector<bool>::const_iterator it = alpha.get_occ_rep().begin();
  std::vector<bool> occ_rep_left;
  occ_rep_left.insert(occ_rep_left.begin(), it+index, alpha.get_occ_rep().end());
  Orbstring o(occ_rep_left);
  boost::shared_ptr<Slater> rightSlater(new Slater(o));
  return rightSlater;
}


void SpinAdapted::ConvertList (std::vector<int>& a, const std::vector<int>& b)
{
  assert (a.size () == b.size ());
  std::vector<int> tmp (a.size ());
  for (int i = 0; i < a.size (); ++i)
    if (a [i] != 0)
      tmp [b [i]] = 1;
  a = tmp;
}

double SpinAdapted::det_energy (const Slater& s, int integralIndex)
{
  // energy of a single slater determinant, used in truncation                                                                            
  double energy = 0;
  Slater bra;
  Slater ket;

  for (int i = 0; i < s.size (); ++i)
    {
      ket = s;
      bra = s;
      energy += ket.trace (bra.d(i).c(i)) *
        v_1[integralIndex] (i,i);
    }

  // diagonal                                                                                                                             
  for (int i = 0; i < s.size (); ++i)
    for (int j = 0; j < s.size (); ++j)
      {
	ket = s;
	bra = s;
	energy += .5 * ket.trace (bra.d(j).d(i).c(j).c(i))
	  * (v_2[integralIndex] (i,j,j,i) - v_2[integralIndex] (i,j,i,j));
      }

  return energy;
}


void SpinAdapted::Slater::connect (const Slater& s, std::vector<int>& cv, std::vector<int>& dv) const
{
  assert ((cv.size() == 0) && (dv.size() == 0));
  this->alpha.connect (s.alpha, cv, dv);
}



void SpinAdapted::Slater::outerProd(const Slater& s, Slater& output) const
{
  vector<bool> occ(Slater().size());

  for (int i=0; i<Slater().size(); i++)
  {
    if (s.alpha.get_occ_rep()[i] == 1 && alpha.get_occ_rep()[i]==1) {
      pout <<"cannot get outerprod of slater determinants\ndet1: "<<s<<"\ndet2: "<<*this<<endl;
      throw 20;
    }
    occ[i] = s.alpha.get_occ_rep()[i] + alpha.get_occ_rep()[i];
  }
  Orbstring o(occ, s.alpha.getSign()*alpha.getSign());
  Slater temp(o);
  output = temp;
}
