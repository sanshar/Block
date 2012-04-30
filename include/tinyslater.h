#ifndef TINY_SLATER_H
#define TINY_SLATER_H

#include <vector>
#include <iostream>
#include <utils.h>
#include <cstdio>
#include <fstream>
class tiny_slater
{
  int sign;
  bool empty;
  vector<int> occ_rep;

public:
  tiny_slater() : sign(0), empty(0) {}
  tiny_slater(const vector<int>& occ, const int sn = 1) : occ_rep (occ), sign (sn), empty(false)
  {
    for (int i=0; i<occ_rep.size(); ++i) occ_rep[i] = (occ_rep[i]) ? occ_rep[i] : 0;
  }
  
  tiny_slater(const int* occ, const int sz, const int sn = 1) : sign (sn), empty (false)
  {
    occ_rep.resize (sz);
    for (int i = 0; i < sz; ++i)
      {
	assert ((occ[i] == 1) || (occ[i] == 0));
	occ_rep[i] = occ [i];
      }
  }

  tiny_slater(const int sz, const int sn = 1) : sign (sn), empty (false)
  {
    occ_rep.resize(sz);
  }

  // accessors
  int& operator[](int i) { return occ_rep[i]; }
  int operator[](int i) const { return occ_rep[i]; }
  int size() const { return occ_rep.size(); }

  // c is create, d is destroy
  tiny_slater& c(int i)
  {
    // applies creation operator i to ket
    {
      if (empty)
	return *this;
      else if (occ_rep[i])
	{
	  empty = true;
	  return *this;
	}
      else 
	{
	  sign *= parity(i);
	  occ_rep[i] = 1;
	}
      return *this;
    }
  }
  tiny_slater& d(int i)
  {
    // applies destruction operator i to ket
    
    if (empty)
      return *this;
    else if (!occ_rep[i])
      {
	empty = true;
	return *this;
    }
    else
      {
	sign *= parity(i);
	occ_rep[i] = 0;
      }
    return *this;
  }

  int parity (int i)
  // determines number of interchanges by number of
  // occ states between end of ket and orb i
  {
    int p = 0;
    for (int j = occ_rep.size() - 1; j > i; --j)
      if (occ_rep[j]) ++p;
    //    return (p & 1) ? -1 : 1;
    return (p % 2) ? -1 : 1;
  }

  int trace(const tiny_slater& s)
  {
    if (empty || s.empty)
      return 0;
    else
      return sign * s.sign * (occ_rep == s.occ_rep);
  }

  void clear()
  {
    sign = 1;
    for (int i=0; i<occ_rep.size (); ++i) occ_rep [i] = 0;
  }

  // returns a list of creation cv and destruction operators dv
  // to produce state s from this (with arbitrary parity)
  void connect (const tiny_slater& s, vector<int>& cv, vector<int>& dv)
  {
    for (int i=0; i<occ_rep.size(); ++i)
      {
	if (occ_rep[i] && !s.occ_rep[i])
	  cv.push_back(i);
	else if (s.occ_rep[i] && !occ_rep[i])
	  dv.push_back(i);
      }
  }

  void savetext(ofstream& f)
  {
    f << empty << " " <<  sign << " " << occ_rep.size() << " " << occ_rep << endl;
  }
  void loadtext(ifstream& f)
  {
    f >> empty >> " " >> sign;
    int adhoc;
    f >> adhoc;
    occ_rep.resize(adhoc);
    for (int i = 0; i << adhoc; ++i)
      f >> occ_rep[i];
  }
  void save(FILE* f)
  {
    fwrite(&empty, sizeof(empty), 1, f);
    fwrite(&sign, sizeof(sign), 1, f);
    int adhoc = occ_rep.size();
    fwrite(&adhoc, sizeof(int), 1, f);
    fwrite(&occ_rep[0], sizeof(int), adhoc, f);
  }
  void load(FILE* f)
  {
    fread(&empty, sizeof(empty), 1, f);
    fread(&sign, sizeof(sign), 1, f);
    int adhoc;
    fread(&adhoc, sizeof(int), 1, f);
    occ_rep.resize(adhoc);
    fread(&occ_rep[0], sizeof(int), adhoc, f);
  } 
  friend ostream& operator<< (ostream& os, const tiny_slater& s)
  {
    if (s.empty) 
    os << "empty" << endl;
    else
      {
	if (s.sign == 1)
	  os << "+" << " ";
	else
	  os << "-" << " ";
	
	for (int i=0; i<s.occ_rep.size(); ++i)
	  {
	    if (s.occ_rep[i])
	      os << 1 << " ";
	    else
	      os << 0 << " ";
	  }
	os << endl;
      }
    return os;
  }
};  

#endif
