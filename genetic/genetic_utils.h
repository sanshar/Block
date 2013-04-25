#ifndef GENETIC_UTILS_H
#define GENETIC_UTILS_H

#include <vector>
#include <newmat.h>

namespace genetic
{
  // random number generator
  inline int irand(const int& n) { return (int) (((double)(rand())) * n / RAND_MAX); }

  inline int irand(const int& nmin, const int& nmax)
  {
    int n = nmax - nmin;
    return nmin + irand(n);
  }

  inline double drand(const double& alpha) { return (alpha * (((double)(rand())) / RAND_MAX)); }

  // random bit
  inline int brand(const double& alpha)
  {
    if(drand(1.0) < alpha) return 1;
    else                   return 0;
  }

  vector<int> RandomBitSequence(const double& alpha, const int& nSize);
  vector<int> RandomSequence(const int& nSize);
  vector<int> ReadSequence(const int& nSize);
  inline vector<int> RandomBitSequence(const int& nSize) { return RandomBitSequence(0.5, nSize); }

  // cross over functions
  vector<int> PMXU(const vector<int>& a, const vector<int>& b);
  vector<int> CrossOver(const vector<int>& a, const vector<int>& b);

  // mutation functions
  vector<int> PointMutation(const vector<int>& a);
  vector<int> GlobalMutation(const vector<int>& a);

};

#endif
