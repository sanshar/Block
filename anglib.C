/*
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012
Copyright (c) 2012, Garnet K.-L. Chan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/*This has been converted from original anglib.f90 code to a c code. The 
the license agreement is reproduced below

!    anglib.f90: angular momentum coupling coefficients in Fortran 90
!    Copyright (C) 1998  Paul Stevenson
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
!    02110-1301,USA
*/
#include <stdio.h>
#include <stdlib.h>
#include "anglib.h"
#include <cmath>
#include <algorithm>

using namespace std;

double cleb_(int j1, int m1, int j2, int m2, int j, int m)
{
  double cleb,factor,sum;
  int par,z,zmin,zmax;


  if (2*(j1/2)-int(2*(j1/2.0)) != 2*(abs(m1)/2)-int(2*(abs(m1)/2.0)) ||
      2*(j2/2)-int(2*(j2/2.0)) != 2*(abs(m2)/2)-int(2*(abs(m2)/2.0)) || 
         2*(j/2)-int(2*(j/2.0)) != 2*(abs(m)/2)-int(2*(abs(m)/2.0)) ||
         j1<0 || j2<0 || j<0 || abs(m1)>j1 || abs(m2)>j2 ||
      abs(m)>j || j1+j2<j || abs(j1-j2)>j || m1+m2!=m)
    cleb= 0.0;
  else
  {  
    factor = 0.0;
    factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2);
    factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2);
    factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2);
    factor = sqrt(factor);
       
    zmin = max(max(0,j2+(j1-m1)/2-(j1+j2+j)/2),j1+(j2+m2)/2-(j1+j2+j)/2);
    zmax = min(min((j1+j2-j)/2,(j1-m1)/2),(j2+m2)/2);
       
    sum=0.0;
    for (z = zmin; z<=zmax; z++) {
      par=1;
      if(2*(z/2)-int(2*(z/2.0)) != 0) 
	par=-1;
      sum=sum+par*binom((j1+j2-j)/2,z)*binom((j1-j2+j)/2,(j1-m1)/2-z)*binom((-j1+j2+j)/2,(j2+m2)/2-z);
    }
    cleb = factor*sum;
  }
  return cleb;
}


double sixj_(int a, int b, int c, int d, int e, int f)
{
  double sixj;
  int nlo, nhi, n;
  double outfactors, sum, sumterm;

  sixj=0.0;
  if((a+b)%2 != c%2) return sixj;
  if((c+d)%2 != e%2) return sixj;
  if((a+e)%2 != f%2) return sixj;
  if((b+d)%2 != f%2) return sixj;
  if(abs(a-b)>c || a+b<c) return sixj;
  if(abs(c-d)>e || c+d<e) return sixj;
  if(abs(a-e)>f || a+e<f) return sixj;
  if(abs(b-d)>f || b+d<f) return sixj;

  outfactors = angdelta(a,e,f)/angdelta(a,b,c);
  outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e);

  nlo = max(max(max( (a+b+c)/2, (c+d+e)/2), (b+d+f)/2), (a+e+f)/2 );
  nhi = min(min( (a+b+d+e)/2, (b+c+e+f)/2), (a+c+d+f)/2);

  sum=0.0;
  for(n=nlo; n<=nhi; n++) {
    if (n%2 == 1)
      sumterm = -1.0;
    else
      sumterm = 1.0;

    sumterm = sumterm * binom(n+1,n-(a+b+c)/2);
    sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2);
    sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2);
    sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2);
    sum=sum+sumterm;
  }

  sixj = sum * outfactors;
  return sixj;
}

double angdelta(int a,int b,int c)
{
  double angdelta_, scr1;

  scr1= factorial((a+b-c)/2);
  scr1=scr1/factorial((a+b+c)/2+1);
  scr1=scr1*factorial((a-b+c)/2);
  scr1=scr1*factorial((-a+b+c)/2);
  angdelta_=sqrt(scr1);
  return angdelta_;
}


double ninej_(int a, int b, int c, int d, int e, int f, int g, int h, int i)
{
  double ninej, sum;
  int xlo, xhi, term;
  int x;

  ninej=0.0;

  if(abs(a-b)>c  || a+b<c) return ninej;
  if(abs(d-e)>f  || d+e<f) return ninej;
  if(abs(g-h)>i  || g+h<i) return ninej;
  if(abs(a-d)>g  || a+d<g) return ninej;
  if(abs(b-e)>h  || b+e<h) return ninej;
  if(abs(c-f)>i  || c+f<i) return ninej;
  
  xlo = max(max(abs(b-f),abs(a-i)),abs(h-d));
  xhi = min(min(b+f,a+i),h+d);
    
  sum=0.0;
  for (int x=xlo; x<=xhi; x=x+2)
  {
    if (x%2 == 1)
      term = -1;
    else
      term = 1;

    sum=sum+term*(x+1)*sixj_(a,b,c,f,i,x)*sixj_(d,e,f,b,x,h)*sixj_(g,h,i,x,a,d);
  }
  ninej=sum;
  return ninej;
}

 

double factorial(int n) {

  double res;

  res = 1.0;
  if (n==0 || n ==1)
    return res;
  for (int i=2; i<=n; i++)
    res *= i;
  return res;
}

double binom(int n, int r)
{
  double res;

  if(n==r || r==0) 
    res = 1.0;
  else if (r==1)
    res = n;
  else
    res = 1.0*n/(n-r)*binom(n-1,r);
  return res;
}

