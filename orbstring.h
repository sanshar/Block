/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_ORBSTRING_HEADER
#define SPIN_ORBSTRING_HEADER
#include <iostream>
#include <vector>
#include <algorithm>
#include <assert.h>

using namespace std;
namespace SpinAdapted{


class Orbstring
{
private:
  std::vector<bool> occ_rep;
  int sign;
  bool empty;

public:
  static void init (const int o);
  Orbstring ();
  inline Orbstring (const std::vector<bool>& occ, const int sn = 1) 
    : occ_rep (occ), sign (sn), empty (false)
  {
    // a little fix for IBMs xlC, see notes in Orbstring::distribute 
    for (int i = 0; i < occ_rep.size (); ++i) occ_rep [i] = abs (occ_rep [i]);   
  }
  
  inline Orbstring (const bool* occ, const int sz, const int sn = 1) 
    : sign (sn), empty (false)
  {
    occ_rep.resize (sz);
    for (int i = 0; i < sz; ++i)
      {
	assert ((occ[i] == 1) || (occ[i] == 0));
	occ_rep [i] = occ [i];
      }
  }

  inline Orbstring (const int* occ, const int sz, const int sn = 1) 
    : sign (sn), empty (false)
  {
    occ_rep.resize (sz);
    for (int i = 0; i < sz; ++i)
      {
	assert ((occ[i] == 1) || (occ[i] == 0));
	occ_rep [i] = occ [i]? 0 : 1;
      }
  }

  bool isempty() const {return empty;}
  int getSign() const {return sign;}
  void setSign(int i) {sign = i;} 
  
  const std::vector<bool>& get_odd_rep()const {return occ_rep;}
  
  inline Orbstring (const int sz, const int sn = 1)
    : sign (sn), empty (false)
  {
    occ_rep.resize (sz);
  }

  const std::vector<bool> &get_occ_rep() const { return occ_rep; }

  Orbstring (const Orbstring& o) : sign (o.sign), empty (o.empty), occ_rep (o.occ_rep) {}
  void operator= (const Orbstring& o) { if (this != &o) {sign = o.sign; empty = o.empty; occ_rep = o.occ_rep;} }

  inline vector<bool>::reference operator[] (int i) { return occ_rep[i]; }
  inline bool operator[] (int i) const { return occ_rep[i]; }

  int size () const;

  // c is create, d is destroy
  Orbstring& c (int i);
  Orbstring& d (int i);  // applies destruction operator i to ket

  inline int parity (int i)
  {
    int p = 0;
    for (int j = 0; j < i; j++)
      if (occ_rep[j]) ++p;
    return (p%2) ? -1 : 1;
  }

  inline int trace (const Orbstring& s)
  {
    if (this->empty || s.empty)
      return 0;
    else
      return this->sign * s.sign * (occ_rep == s.occ_rep);
  }

  void clear ()
  {
    sign = 1;
    for (int i = 0; i < occ_rep.size (); ++i) occ_rep [i] = 0;
  }

  // returns a list of creation cv and destruction operators dv
  // to produce this state from s (with arbitrary parity)
  void connect (const Orbstring& s, std::vector<int>& cv, std::vector<int>& dv) const
  {
    for (int i = 0; i < occ_rep.size(); ++i)
      {
	if (occ_rep[i] && !s.occ_rep[i])
	  cv.push_back(i);
	else if (s.occ_rep[i] && !occ_rep[i])
	  dv.push_back(i);
      }
  }

  friend ostream& operator<< (ostream& os, const Orbstring& s){
    if (s.empty) 
      os << "empty" << endl;
    else
      {
	if (s.sign == 1)
	  os << "+" << " ";
	else
	  os << "-" << " ";
	
	for (int i = 0; i < s.occ_rep.size(); ++i)
	  {
	    if (s.occ_rep[i])
	      os << 1 << " ";
	    else
	      os << 0 << " ";
	  }
      }
    return os;
  }

};
}
  
#endif
