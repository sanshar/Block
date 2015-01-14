/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "couplingCoeffs.h"
#include "global.h"
#include "pario.h"

namespace SpinAdapted{

ninejCoeffs& ninejCoeffs::getinstance()
{
  static ninejCoeffs nj(2);
  return nj;
}

void ninejCoeffs::init(int maxj_)
{
  if (maxj != maxj_) {
    maxj = maxj_;
    buildArray();
  }
  else
    return;

}

void ninejCoeffs::buildArray()
{
  p2out << "Building Array with maxj: "<<maxj<<endl;
  coeffs.resize(10);
  initarray(0, 0, 0, 0);
  initarray(1, 1, 1, 0);
  initarray(2, 2, 2, 0);
  initarray(3, 1, 1, 2);
  initarray(4, 2, 0, 2);
  initarray(5, 0, 2, 2);
  initarray(6, 1, 2, 1);
  initarray(7, 2, 1, 1);
  initarray(8, 0, 1, 1);
  initarray(9, 1, 0, 1);
}

ninejCoeffs::ninejCoeffs(int maxj_): maxj(maxj_)
{
  buildArray();
}


double ninejCoeffs::operator() (int ja, int jb, int jc, int jd, int je, int jf, int jg, int jh, int ji) const
{
  if (!dmrginp.spinAdapted()) return 1.0;
  int a = jd, b=je, c=jf;
  int index = arrayindex(a, b, c);
  if (index == -1)
    return Ninej(ja, jb, jc, jd, je, jf, jg, jh, ji);
  if (ja >= maxj || jb >=maxj || jc >=maxj) {
    return Ninej(ja, jb, jc, jd, je, jf, jg, jh, ji);
  }
  int imax = a+1, jmax = b+1, kmax = c+1;
  int ii = (jg+a-ja)/2, jj = (jh+b-jb)/2, kk = (ji+c-jc)/2;

  return coeffs[index](ja*imax+ii, jb*jmax+jj, jc*kmax + kk);
}

void ninejCoeffs::initarray(int arrayindex, int a, int b, int c)
{
  int imax = a+1, jmax = b+1, kmax = c+1;
  coeffs[arrayindex].resize(imax*maxj, jmax*maxj, kmax*maxj);
  for (int i=0; i<maxj; i++)
    for (int j=0; j<maxj; j++)
      for (int k=0; k<maxj; k++)
	for (int ii = 0; ii<imax; ii++)
	  for (int jj = 0; jj<jmax; jj++)
	    for (int kk = 0; kk<kmax; kk++)
	      if (i+2*ii-a < 0 ||  j+2*jj-b<0 || k+2*kk-c <0)
		coeffs[arrayindex](i*imax + ii, j*jmax +jj, k*kmax +kk) = 0.0;
	      else
		coeffs[arrayindex](i*imax + ii, j*jmax +jj, k*kmax +kk) = Ninej(i, j, k,  a, b, c, i+2*ii-a, j+2*jj-b, k+2*kk-c);
		
}

int ninejCoeffs::arrayindex(int a, int b, int c) const
{
  if (a == 0 && b == 0 && c == 0)
    return 0;
  else if (a == 1 && b == 1 && c == 0)
    return 1;
  else if (a == 2 && b == 2 && c == 0)
    return 2;
  else if (a == 1 && b == 1 && c == 2)
    return 3;
  else if (a == 2 && b == 0 && c == 2)
    return 4;
  else if (a == 0 && b == 2 && c == 2)
    return 5;
  else if (a == 1 && b == 2 && c == 1)
    return 6;
  else if (a == 2 && b == 1 && c == 1)
    return 7;
  else if (a == 0 && b == 1 && c == 1)
    return 8;
  else if (a == 1 && b == 0 && c == 1)
    return 9;
  else
    {
      return -1;
    }
}


}
