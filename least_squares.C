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


void least_squares(vector<double> x, vector<double> y){
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
   int max;
   int n; //size of n
   int n_low; //n- max
   int i;

   double f;
   max=3; //take the max number of values
   n = x.size();
   n_low = n - max;



   for (i=n_low; i<n; i++) {

      avgx+=x[i];
      avgy+=y[i];
      avgxx+=x[i]*x[i];
      avgxy+=x[i]*y[i];
      avgyy+=y[i]*y[i];
   }
   avgx/=max;
   avgy/=max;
   avgxx/=max;
   avgxy/=max;
   avgyy/=max;

   beta=(avgxy-avgx*avgy)/(avgxx-avgx*avgx);
   alpha = avgy-beta*avgx;

   //*a = alpha;

   for (i=n_low; i<n; i++) {
      f=alpha + beta*x[i];
      sserr += pow(y[i]-f, 2);
      sstot += pow(y[i]-avgy, 2);
   }

   r2 = 1. - sserr/sstot;
   //pout << "beta " << beta << " alpha " << alpha << " r2 " << r2 << endl;
   //pout << "Extrapolated energy: " << alpha << " a.u." << endl;
   pout << "\n\t\t\tExtrapolated Energy = " << fixed << setprecision(10) << alpha << " a.u." << endl << endl;
}
