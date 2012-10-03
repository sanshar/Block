
#include <stdio.h>
#include <stdlib.h>
#include "new_anglib.h"
#include "anglib.h"
#include <cmath>
#include <algorithm>

#include <iostream>
using namespace std;

/*
double six-j(int j1, int j2, int j3, int l1, int l2, int l3) {
	double sixj;
	
   sixj=delta(j1, j2, j3)*delta(j1, l2, l3) * delta(l1, j2, j3) *
         delta(l1,j2,j3) * w6j(j1, j2, j3, l1, l2, l3); 
}
//end six-j

double w6j(j1, j2, j3, l1, l2, l3) {
}
//end w6j
*/

double three_j(int j1, int j2, int j3, int m1, int m2, int m3) {
   double cleb =0.0;
	double threej = 0.0;
   double fj1, fj2, fj3, fm3;
   fj1 = j1/2.;
   fj2 = j2/2.;
   fj3 = j3/2.;
   fm3 = m3/2.;
   cleb = clebsch(j1, m1, j2, m2, j3, m3);
   threej = mone(fj1-fj2+fm3)*cleb/sqrt(2*fj3+1);
   return threej;
}

double clebsch(int nj1, int nm1, int nj2, int nm2, int nj3, int nm3) {
   double j1, j2, j3;
   double m1, m2, m3;

   //Converting to half its value
   j1=nj1/2.;
   j2=nj2/2.;
   j3=nj3/2.;
   m1=nm1/2.;
   m2=nm2/2.;
   m3=nm3/2.;

   double cleb=0.0;
   if ( j1 < 0 || j2 < 0 || j3 < 0 || abs(m1) > j1 || abs(m2) > j2 ||
      abs(m3) > j3 || j1 + j2 < j3 || abs(j1-j2) > j3 || m1 + m2 != m3) {

      cleb=0.0;
   }
   else
   {
      double factor = 0.0;
      double sum = 0.0;
      int t;

      double num1 = pow(2*j3+1,2);
      double num2 = fbinom(j1+j2+j3+1, j1+j2-j3);
      double num3 = fbinom(2*j3, j3+m3);
      double den1 = (2*j1+1);
      double den2 = (2*j2+1);
      double den3 = fbinom(j1+j2+j3+1, j1-j2+j3);
      double den4 = fbinom(j1+j2+j3+1, j2-j1+j3);
      double den5 = fbinom(2*j1, j1+m1);
      double den6 = fbinom(2*j2, j2+m2);

      double num = num1*num2*num3;
      double den = den1*den2*den3*den4*den5*den6;
      factor = sqrt(num/den);

      double mint = max(max(0., j1-m1-(j3-m3)), j2 + m2 - (j3 + m3)); 
      double maxt = min(min(j1-m1, j2+m2),j1+j2-j3);
      
      cout << "mint " << mint << endl;
      cout << "maxt " << maxt << endl;
      double bin1;
      double bin2;
      double bin3;
      for (t=mint; t<=maxt; t++) {
         bin1=fbinom(j1+j2-j3, t);
         bin2=fbinom(j3-m3,     j1-m1-t);
         bin3=fbinom(j3+m3,     j2+m2-t);
         sum = sum + mone(t)*bin1*bin2*bin3;

         cout << "t " << t << endl;
         cout << "sum " << sum << endl;
         cout << "bin1 " << bin1 << endl;
         cout << "bin2 " << bin2 << endl;
         cout << "bin3 " << bin3 << endl;
      }

      cleb = factor*sum;
      cout << "factor: " << factor << endl;
      cout << "sum: " << sum << endl;
      cout << "Clebsch: " << cleb << endl;
   }
      return cleb;
}

double facto(double n) {
	double fac;
   int nint;
   nint = get_cast(n);
	fac=1.0;
	int i;
	if (n==0 || n==1)
		return fac;
	for (i=2; i<=nint; i++) 
		fac *= i;
	return fac;
}

int mone (double n) {
	int value;
   int nint;
   nint = get_cast(n);
   //cout << "nint %2 " << nint %2 << endl;
	if (nint % 2 == 0)
		value = 1;
	else
		value = -1;
	return (value);
}	
		
int get_cast(double x) {
   int i;
   //i = (x / (int) x >= 1) ? (int) x : (int) x + 1 ;
   i = (int) x;
   //cout << "x " << x << " i " << i << endl;
   return i;
}

double fbinom(double dn, double dr)
{
  double res;
  int n = get_cast(dn);
  int r = get_cast(dr);


  if(n==r || r==0) 
  {
    res = 1.0;
  }
  else if (r==1)
    res = n;
  else
    res = 1.0*n/(n-r)*fbinom((double)n-1,(double)r);
//    cout << n << " " << r<< " -> " << res << endl;
  return res;
}

















