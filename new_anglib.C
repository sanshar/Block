
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



/*
double clebsch(int j1, int m1, int j2, int m2, int j3, int m3){
   // calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
   // arguments are integer and twice the true value. 
   double threej = 0.0;
   threej= three_j(j1, j2, j3, m1, m2, -m3);
   double cleb=0.0;
   double dj1, dj2, dj3;
   double dm3;

   dj1=j1/2.;
   dj2=j2/2.;
   dj3=j3/2.;
   dm3=m3/2.;
   //cout << "-j1 +j2 -m3 "  << -j1 +j2 -m3 << endl;
   //cout << "-dj1 +dj2 -dm3 "  << -dj1 +dj2 -dm3 << endl;
   //cout << "mone " << mone(-dj1 +dj2 -dm3) << endl;
   //cout << "dj3 " << dj3 << " j3 " << j3 << endl;
   cleb = mone(-dj1 +dj2 -dm3)*sqrt(2*dj3+1)*threej;
   cout << "threej " << threej << endl;;
   //cout << "Clebsch " << cleb << endl;;
   return cleb;
}
*/
/*
double three_j(int j1, int j2, int j3, int m1, int m2, int m3) {
	double threej;
	int expo;
   int prefac;
   double delta=0.0;
   double dobleu=0.0;
   double dj1, dj2, dj3;
   double dm1, dm2, dm3;

   //cout << "Entering three_j new" << endl;

	//if (j1<0 || j2 < 0 || j3 < 0 || abs(m1) < j1 || abs(m2) < j2 || 
	//	abs(m3) < j3 || j1 + j2 < j3 || abs(j1 -j2) > j3 || m1 + m2 != m3)

   if (j1 < 0 || j2 < 0 || j3 < 0 || abs(m1) > j1 || abs(m2) > j2 || abs(m3) > j3 ||
         j1 + j2 -j3 < 0 || j2 + j3 -j1 < 0 || j3 + j1 -j2 < 0 || m1 +m2 + m3 != 0) {

      //cout << j1 << endl;
      //cout << j2 << endl;
      //cout << j3 << endl;
      //cout << m1 << endl;
      //cout << m2 << endl;
      //cout << m3 << endl;
		threej = 0.0;
   }
	else {
      dj1=j1/2.;
      dj2=j2/2.;
      dj3=j3/2.;
      dm1=m1/2.;
      dm2=m2/2.;
      dm3=m3/2.;

		threej=0.0;
		expo= dj1-dj2-dm3;
		prefac=mone(expo);
      delta=cj_delta(dj1, dj2, dj3);
      dobleu = w3j(dj1, dj2, dj3, dm1, dm2, dm3);
		//threej=prefac * cj_delta(j1, j2, j3) * w3j(j1,j2,j3, m1, m2, m3);
		threej=prefac * delta * dobleu;
      //cout << "factor: " << delta << endl;
      //cout << "w3j: " << dobleu << endl;
      //cout << "threej: " << threej << endl;
   }
   return threej;
}
*/
//end three_j

double w3j(double j1, double j2, double j3, double m1, double m2, double m3) {

   int v;
   int vmin;
   int vmax;
   double prod=0.0;
   double facall=0.0;
   double fac1=0.0;
   double fac2=0.0;
   double fac3=0.0;
   double fac4=0.0;
   double fac5=0.0;
   double fac6=0.0;

   //vmin = max(max(0., j3-j2+2*m1), j3-j1-m2);
   //vmax = min(min(j1-m1, j2+m2), j1+j2-j3);

   vmin = min(j3-j1-m2, j3-j2+m1);
   //cout << "j3 j1 m2, j3-j1-m2 " << j3 << " " << j1 << " " << m2 << ", " << j3-j1-m2 << endl;
   //cout << "j3 j2 m1, j3-j2+m1 " << j3 << " " << j2 << " " << m1 << ", " << j3-j2+m1 << endl;
   vmin = -vmin;
   vmin = max(vmin, 0);

   vmax = min(min(j1-m1, j2+m2), j1+j2-j3);
   cout << "vmin: " << vmin << endl;
   cout << "vmax: " << vmax << endl;
   //cout << "j1 m1, j1-m1 " << j1 << " " << m1 << " " << j1-m1 << endl;
   //cout << "j2 m2, j2+m2 " << j2 << " " << m2 << " " << j2+m2 << endl;
   //cout << "j1 j2 j3, j1+j2-j3 " << j1 << " " << j2 << " " << j3 << ", " << j1+j2-j3 << endl;
   
   double prefac = 0.0;

   int j1m1_1 = facto(j1 + m1);
   int j1m1_2 = facto(j1 - m1);
   int j2m2_1 = facto(j2 + m2);
   int j2m2_2 = facto(j2 - m2);
   int j3m3_1 = facto(j3 + m3);
   int j3m3_2 = facto(j3 - m3);

   prefac = j1m1_1*j1m1_2*j2m2_1*j2m2_2*j3m3_1*j3m3_2;
   prefac = sqrt(prefac);

   for (v=vmin; v<=vmax; v++) { //  for (v=0; v<20; v++) {
  //    if ((j1+j2-j3-v < 0) && (j1-m1-v<0) &&  (j2+m2-v < 0 ) && (j3-j2+m1-v < 0) && (j3-j1-m2-v>0))
   //   {

      fac1=facto(j1+j2-j3-v);
      fac2=facto(j1-m1-v);
      fac3=facto(j2+m2-v);
      fac4=facto(j3-j2+m1+v);
      fac5=facto(j3-j1-m2+v);
      fac6 =facto(v);

      //cout << "1 " << j1+j2-j3-v << endl;
      //cout << "2 " << j1-m1-v << endl;
      //cout << "3 " << j2+m2-v << endl;
      //cout << "4 " << j3-j2+m1+v << endl;
      //cout << "5 " << j3-j1-m2+v << " fac 5 " << fac5 << endl;
      //cout << "6 " << v << endl;

      facall = fac1*fac2*fac3*fac4*fac5*fac6;
      prod = prod + mone(v)/facall;
      //cout << "v " << v << " prod: " << prod << endl;
   }
   //}
   cout << "prefac: " << prefac << endl;
   cout << "produc: " << prod << endl;

   prod = prefac*prod;
   cout << "prod: " << prod << endl;
   return prod;
}
//end w3j

double cj_delta(double a, double b, double c) {
	double prefac = 0.0;
	double num1 = facto(a+b-c);
   //cout << "a+b-c " << a+ b -c << endl;
	double num2 = facto(a-b+c);
   //cout << "a-b+c " << a - b + c << endl;
	double num3 = facto(-a+b+c);
   //cout << "-a+b+c " << - a + b + c << endl;
	double den1 = facto( a + b + c + 1);
   //cout << "a+b+c+1 " << a + b + c + 1 << endl;

	prefac=num1*num2*num3/den1;
	prefac=sqrt(prefac);
	return prefac;
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

















