#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "ReadIntegral.h"
using namespace std;
using namespace boost::algorithm;

void genetic::ReadIntegral(ifstream& fdump, Matrix& K)
{
  // save file pointer
  ifstream::pos_type fp = fdump.tellg();
  // rewind
  fdump.seekg(0, ios::beg);

  char sbuf[256];
  fdump.getline(sbuf, 256);
  string entry(sbuf);

  vector<string> fields;
  split(fields, entry, is_any_of("=, \t"), token_compress_on);

  int nOrbs = 0;
  for(int i = 0; i < fields.size(); ++i)
  {
    if(fields[i] == "NORB")
    {
      nOrbs = atoi(fields[i+1].c_str());
      break;
    }
  }

  K.ReSize(nOrbs, nOrbs); K = 0.0;

  while(fdump >> entry) if((entry == "&END") || (entry == "/")) break;

  int i, j, k, l;
  double v;
  while(fdump >> v >> i >> j >> k >> l)
  {
    if(i == 0 && j == 0 && k == 0 && l == 0) continue;
    //* Read by Mulliken Notation
    if(i == k && j == l)
    {
      i--; j--;
      K.element(i, j) += fabs(v);
      K.element(j, i) = K.element(i, j);
    }

    if(k == 0 && l == 0)
    {
      i--; j--;
      K.element(i, j) += 1.0e-7 * fabs(v);
      K.element(j, i)  = K.element(i, j);
    }
  }
  // set file pointer to original position
  fdump.clear();
  fdump.seekg(fp);
}

void genetic::ReadIntegral_nevpt(ifstream& fdump, Matrix& K, int nact)
{
  // save file pointer
  ifstream::pos_type fp = fdump.tellg();
  // rewind
  fdump.seekg(0, ios::beg);

  char sbuf[256];
  fdump.getline(sbuf, 256);
  string entry(sbuf);

  vector<string> fields;
  split(fields, entry, is_any_of("=, \t"), token_compress_on);

  int nOrbs = 0;
  for(int i = 0; i < fields.size(); ++i)
  {
    if(fields[i] == "NORB")
    {
      nOrbs = atoi(fields[i+1].c_str());
      break;
    }
  }

  K.ReSize(nact, nact); K = 0.0;

  while(fdump >> entry) if((entry == "&END") || (entry == "/")) break;

  int i, j, k, l;
  double v;
  while(fdump >> v >> i >> j >> k >> l)
  {
    if(i == 0 && j == 0 && k == 0 && l == 0) break;
    if (i > nact) break;
    //* Read by Mulliken Notation
    if(i == k && j == l)
    {
      i--; j--;
      K.element(i, j) += fabs(v);
      K.element(j, i) = K.element(i, j);
    }

    if(k == 0 && l == 0)
    {
      i--; j--;
      K.element(i, j) += 1.0e-7 * fabs(v);
      K.element(j, i)  = K.element(i, j);
    }
  }
  // set file pointer to original position
  fdump.clear();
  fdump.seekg(fp);
}

void genetic::ReadIntegral_BCS(ifstream& fdump, Matrix& K)
{
  // save file pointer
  ifstream::pos_type fp = fdump.tellg();
  // rewind
  fdump.seekg(0, ios::beg);

  char sbuf[256];
  fdump.getline(sbuf, 256);
  string entry(sbuf);

  vector<string> fields;
  split(fields, entry, is_any_of("=, \t"), token_compress_on);

  int nOrbs = 0;
  for(int i = 0; i < fields.size(); ++i)
  {
    if(fields[i] == "NORB")
    {
      nOrbs = atoi(fields[i+1].c_str());
      break;
    }
  }

  K.ReSize(nOrbs, nOrbs); K = 0.0;  
  while(fdump >> entry) if((entry == "&END") || (entry == "/")) break;

  int i, j, k, l;
  double v;
  while(fdump >> v >> i >> j >> k >> l)
  {
    if(i == 0 && j == 0 && k == 0 && l == 0) continue;
    //* Read by Mulliken Notation
    if(i == k && j == l)
    {
      i--; j--;
      K.element(i, j) += fabs(v);
      K.element(j, i) = K.element(i, j);
    }
    if(i == l && j == k)
    {
      i--; j--;
      K.element(i, j) += fabs(v);
      K.element(j, i) = K.element(i, j);
    }
    if(k == 0 && l == 0)
    {
      i--; j--;
      K.element(i, j) += 1.0e-7 * fabs(v);
      K.element(j, i)  = K.element(i, j);
    }
  }
  // set file pointer to original position
  fdump.clear();
  fdump.seekg(fp);
}
