#include <vector>
#include <algorithm>
#include "genetic_utils.h"
using namespace std;

vector<int> genetic::PointMutation(const vector<int>& original)
{
  int nSize = original.size();
  vector<int> mutated(original);
  int nTry = irand(3) + 1;
  for(int i = 0; i < nTry; ++i)
  {
    int j1 = irand(nSize); int j2 = irand(nSize);
    if(j1 != j2) swap(mutated[j1], mutated[j2]);
  }
  return mutated;
}

vector<int> genetic::GlobalMutation(const vector<int>& original)
{
  int nSize = original.size();
  vector<int> mutated(nSize, 0);
  vector<int> np(4, 0);
  for(int i = 0; i < 4; ++i) np[i] = irand(nSize);
  sort(np.begin(), np.end());
  int n = 0;
  for(int i = 0;     i < np[0]; ++i) mutated[n++] = original[i];
  for(int i = np[2]; i < np[3]; ++i) mutated[n++] = original[i];
  for(int i = np[1]; i < np[2]; ++i) mutated[n++] = original[i];
  for(int i = np[0]; i < np[1]; ++i) mutated[n++] = original[i];
  for(int i = np[3]; i < nSize; ++i) mutated[n++] = original[i];
  return mutated;
}
