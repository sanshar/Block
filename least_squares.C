#include <cmath>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <fstream>
//#include <vector>
#include <boost/serialization/vector.hpp>

#include "least_squares.h"
#include "pario.h"

using namespace std;

vector<double> dw;
vector<double> energy;


void least_squares(vector<double> x, vector<double> y){//, double* a){
   double avgx = 0.0;
   double avgy = 0.0;
   double avgxx = 0.0; 
   double avgyy = 0.0; 
   double avgxy = 0.0; 
   double sserr = 0.0; 
   double sstot = 0.0; 
   double beta;
   double alpha;
   double r2;
   int n; //size of n
   int i;

   double f;
   n = x.size();



   for (i=0; i<n; i++) {
      avgx+=x[i];
      avgy+=y[i];
      avgxx+=x[i]*x[i];
      avgxy+=x[i]*y[i];
      avgyy+=y[i]*y[i];
   }
   avgx/=n;
   avgy/=n;
   avgxx/=n;
   avgxy/=n;
   avgyy/=n;

   beta=(avgxy-avgx*avgy)/(avgxx-avgx*avgx);
   alpha = avgy-beta*avgx;

   //*a = alpha;

   for (i=0; i<n; i++) {
      f=alpha + beta*x[i];
      sserr += pow(y[i]-f, 2);
      sstot += pow(y[i]-avgy, 2);
   }

   r2 = 1. - sserr/sstot;
   //cout << "beta " << beta << " alpha " << alpha << " r2 " << r2 << endl;
   //pout << "Extrapolated energy: " << alpha << " a.u." << endl;
#ifndef MOLPRO
   printf("\n\t\t\tExtrapolated Energy = %20.10f a.u.\n",alpha);
   cout << "\n\t\t\tExtrapolated Energy = " << fixed << setprecision(10) << alpha << " a.u." << endl << endl;
#else 
   xout << "\n\t\t\tExtrapolated Energy = " << fixed << setprecision(10) << alpha << " a.u." << endl << endl;
#endif
}
