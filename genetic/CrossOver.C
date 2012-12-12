#include <map>
#include "genetic_utils.h"
using namespace std;

// ordered crossover ( Partically Mapped Crossover-Uniform, PMXU )
vector<int> genetic::PMXU(const vector<int>& a, const vector<int>& b)
{
  int nSize = a.size();
  map<int, int> aMap;
  for(int i = 0; i < nSize; ++i) aMap.insert(make_pair(a[i], i));
  map<int, int> bMap;
  for(int i = 0; i < nSize; ++i) bMap.insert(make_pair(b[i], i));

  vector<int> mask(RandomBitSequence(nSize));
  vector<int> c(nSize, -1);

  vector<int> bNew(b);
  for(int i = 0; i < nSize; ++i)
  {
    if(mask[i] == 1)
    {
      c[i] = a[i];
      int k = aMap.find(b[i])->second;
      if(mask[k] == 0)
      {
        int j = bMap.find(a[i])->second;
        while(mask[j]) j = bMap.find(a[j])->second;
        swap(bNew[i], bNew[j]);
      }
    }
  }
  for(int i = 0; i < nSize; ++i) if(mask[i] == 0) c[i] = bNew[i];

  return c;
}

vector<int> genetic::CrossOver(const vector<int>& a, const vector<int>& b)
{
  vector<int> gene(PMXU(a, b));
  return gene;
}
