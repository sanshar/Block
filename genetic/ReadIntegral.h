#ifndef READ_INTEGRAL_H
#define READ_INTEGRAL_H

#include <fstream>
#include <newmat.h>
using namespace std;

namespace genetic
{
  void ReadIntegral(ifstream& fdump, Matrix& K);
  void ReadIntegral_nevpt(ifstream& fdump, Matrix& K, int nact);
  void ReadIntegral_BCS(ifstream& fdump, Matrix& K);
};

#endif
