/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_MATRIX_DECLARATION_HEADER
#define SPIN_MATRIX_DECLARATION_HEADER
#include "pario.h"
#include <newmat.h>
#include <newmatap.h>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include "blas_calls.h"
#include <boost/pool/poolfwd.hpp>
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/serialization/serialization.hpp>

// basic templated matrix class 
// with simple deep copy semantics
// intended for use when access to matrices or copying is infrequent

namespace SpinAdapted{
template<class T> class ObjectMatrix
{
public:
  //  typedef Elementtype T;
//  vector<T, boost::pool_allocator<T> > rep;
  std::vector<T> rep;
  int nrs;
  int ncs;
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & nrs & ncs & rep;
  }

public:
  ObjectMatrix (void) : nrs (0), ncs (0), rep () {}
  ObjectMatrix (int s, int t) : nrs (s), ncs (t) { rep.resize (s * t); }
  ObjectMatrix (const ObjectMatrix<T>& m) : nrs (m.nrs), ncs (m.ncs) { rep = m.rep; }
  T& operator() (int s, int t) 
  { 
    assert ((s >= 0) && (t >= 0) && (s < nrs) && (t < ncs));
    return rep.at(s * ncs + t);
  }
  const T& operator() (int s, int t) const 
  { 
    assert ((s >= 0) && (t >= 0) && (s < nrs) && (t < ncs));
    return rep.at(s * ncs + t);
  }

  T& operator () (int s, int t, char conj)
  {
    if (conj == 'n')
      assert ((s >= 0) && (t >= 0) && (s < nrs) && (t < ncs));
    else
      assert ((s >= 0) && (t >= 0) && (t < nrs) && (s < ncs));
    return (conj == 'n') ? rep.at(s * ncs + t) : rep.at(t * ncs + s);
  }

  const T& operator () (int s, int t, char conj) const
  {
    if (conj == 'n')
      assert ((s >= 0) && (t >= 0) && (s < nrs) && (t < ncs));
    else
      assert ((s >= 0) && (t >= 0) && (t < nrs) && (s < ncs));
    return (conj == 'n') ? rep.at(s * ncs + t) : rep.at(t * ncs + s);
  }

  ObjectMatrix<T>& operator= (const ObjectMatrix<T>& m)
    {
      if (this != &m)
	{
	  rep = m.rep; nrs = m.nrs; ncs = m.ncs;
	}
      return *this;
    }

  int ncols () const { return ncs; }
  int nrows () const { return nrs; }
  int Ncols () const { return ncs; }
  int Nrows () const { return nrs; }
  int Ncols (char t) const { return (t == 'n') ? ncs : nrs; }
  int Nrows (char t) const { return (t == 'n') ? nrs : ncs; }
  int capacity () const { return rep.capacity (); }
 
  void resize (int s, int t) { rep.resize (s * t); nrs = s; ncs = t; } 
  void ReSize (int s, int t) { rep.resize (s * t); nrs = s; ncs = t; }
  friend ostream& operator<< (ostream& os, const ObjectMatrix<T>& m) 
  { 
    for (int i = 0; i < m.nrs; ++i)
      {
	for (int j = 0; j < m.ncs; ++j)
	  os << m (i,j) << "\t";
	os << endl;
      }
    return os;
  }


  void Clear () { rep.clear (); nrs = 0; ncs = 0; }

};



template<class T> class ObjectLowerTriangularMatrix
{
  std::vector< std::vector<T> > rep;
  int n;
public:
  ObjectLowerTriangularMatrix (void) : n (0) {}
  ObjectLowerTriangularMatrix (int s) : n (s) 
  { 
    rep.resize (s); 
    for (int i = 0; i < s; ++i) rep [i].resize (i+1);
  }
  ObjectLowerTriangularMatrix (const ObjectLowerTriangularMatrix<T>& m) : rep (m.rep), n (m.n) {}
  T& operator() (int s, int t) 
  { 
    assert (s >= 0 && t >= 0 && s < n && t < n &&  s >= t);
    return rep [s][t];
  }
  const T& operator() (int s, int t) const 
  { 
    assert (s >= 0 && t >= 0 && s < n && t < n && s >= t);
    return rep [s][t];
  }
  
  /* void operator= (const ObjectMatrix<T>& m) */
  /*   { */
  /*     if (this != &m) */
  /* 	{ */
  /* 	  rep = m.rep; n = m.n; */
  /* 	} */
  /*     return *this; */
  /*   } */
  ObjectLowerTriangularMatrix<T>& operator= (const ObjectMatrix<T>& m)
    {
      if (this != &m)
	{
	  rep = m.rep; n = m.n;
	}
      return *this;
    }

  int ncols () const { return n; }
  int nrows () const { return n; }
  // non-destructive resizing
  void Clear () { rep.clear (); n = 0; }
  void resize (int s) 
  { 
    rep.resize (s); 
    for (int i = 0; i < s; ++i) rep [i].resize (i+1); 
    n = s;
  }
};

template<class T> class ObjectMatrix3D
{
  std::vector< std::vector< std::vector<T> > > rep;
  int n0;
  int n1;
  int n2;

public:
  ObjectMatrix3D (void) : n0 (0), n1 (0), n2 (0) {}
  ObjectMatrix3D (int m0, int m1, int m2) : n0 (m0), n1 (m1), n2 (m2) { ReSize (n0, n1, n2); }
  
  void ReSize (int m0, int m1, int m2)
  {
    n0 = m0; n1 = m1; n2 = m2;
    rep.resize (n0);
    for (int i = 0; i < n0; ++i)
      {
	rep [i].resize (n1);
	for (int j = 0; j < n1; ++j)
	  rep [i][j].resize (n2);
      }
  }
  
  T& operator() (int i, int j, int k) { return rep [i][j][k]; }
  const T& operator() (int i, int j, int k) const { return rep [i][j][k]; }

  int NDim0 () const { return n0; }
  int NDim1 () const { return n1; }
  int NDim2 () const { return n2; }
  int Nrows0 () { return n0; }
  int Nrows1 () { return n1; }
  int Nrows2 () { return n2; }
};

template<class T> class ObjectMatrix4D
{
  std::vector< std::vector< std::vector<std::vector<T> > > > rep;
  int n0;
  int n1;
  int n2;
  int n3;

public:
  ObjectMatrix4D (void) : n0 (0), n1 (0), n2 (0), n3(0) {}
  ObjectMatrix4D (int m0, int m1, int m2, int m3) : n0 (m0), n1 (m1), n2 (m2), n3 (m3) { ReSize (n0, n1, n2, n3); }
  
  void ReSize (int m0, int m1, int m2, int m3)
  {
    n0 = m0; n1 = m1; n2 = m2; n3 = m3;
    rep.resize (n0);
    for (int i = 0; i < n0; ++i)
      {
	rep [i].resize (n1);
	for (int j = 0; j < n1; ++j)
        {
	  rep [i][j].resize (n2);
          for( int k =0; k< n2; ++k)
	    rep [i][j][k].resize (n3);
        }
      }
  }
  
  T& operator() (int i, int j, int k, int l) { return rep [i][j][k][l]; }
  const T& operator() (int i, int j, int k, int l) const { return rep [i][j][k][l]; }

  int NDim0 () const { return n0; }
  int NDim1 () const { return n1; }
  int NDim2 () const { return n2; }
  int NDim3 () const { return n3; }
  int Nrows0 () { return n0; }
  int Nrows1 () { return n1; }
  int Nrows2 () { return n2; }
  int Nrows3 () { return n3; }
};
}
#endif
