#include <cmath>
#include "Evaluate.h"

void genetic::Permute(const Matrix& K, const vector<int>& sequence, Matrix& Kp)
{
  int nSize = sequence.size();
  Kp.ReSize(nSize, nSize); Kp = 0.0;
  for(int i = 0; i < nSize; ++i)
  {
    int ip = sequence[i];
    for(int j = i + 1; j < nSize; ++j)
    {
      int jp = sequence[j];
      Kp.element(i, j) = K.element(ip, jp);
    }
  }
}

double genetic::Evaluate(const double& scale, const double& np, const Gene& gene, const Matrix& K)
{
  Matrix Kp;
  Permute(K, gene.Sequence(), Kp);

  double weight = 0.0;
  for(int i = 0; i < Kp.Nrows(); ++i)
    for(int j = i + 1; j < Kp.Ncols(); ++j)
    {
      double d = j - i;
//    weight += Kp.element(i, j) * pow(d, np);
      weight += Kp.element(i, j) * d * d;
    }

  return scale * weight;
}
