#include <iostream>
#include <map>
#include <cstdlib>
#include "genetic_utils.h"
using namespace std;

// Random Bit Sequence Generator
vector<int> genetic::RandomBitSequence(const double& alpha, const int& nSize)
{
  vector<int> sequence(nSize, 0);
  for(int i = 0; i < nSize; ++i) sequence[i] = brand(alpha);
  return sequence;
}

// Random Sequence Generator
vector<int> genetic::RandomSequence(const int& nSize)
{
  vector<int> sequence(nSize, 0);
  multimap<double, int> SeqMap;
  for(int i = 0; i < nSize; ++i) SeqMap.insert(make_pair(drand(1.0), i));
  multimap<double, int>::iterator it = SeqMap.begin();
  for(int i = 0; i < nSize; ++i) sequence[i] = (it++)->second;
  return sequence;
}

vector <int> genetic::ReadSequence(const int& nSize) {
  vector<int> sequence(nSize, 0);
  multimap<double, int> SeqMap;
  for(int i = 0; i < nSize; ++i) SeqMap.insert(make_pair(drand(1.0), i));
  multimap<double, int>::iterator it = SeqMap.begin();
  //ROA: need to continue
  return sequence;
}
