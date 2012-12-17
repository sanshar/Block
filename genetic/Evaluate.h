#ifndef EVALUATE_H
#define EVALUATE_H

#include <vector>
#include <newmat.h>
#include "Gene.h"

namespace genetic
{
  void Permute(const Matrix& K, const vector<int>& sequence, Matrix& Kp);
  double Evaluate(const double& scale, const double& np, const Gene& gene, const Matrix& K);
};

#endif
