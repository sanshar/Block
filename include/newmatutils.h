/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef NEWMATUTILS_HEADER_H
#define NEWMATUTILS_HEADER_H
#include <newmat.h>
#include <fstream>
#include <cstdio>

using namespace std;

inline void savenewmat(char* filename, Matrix& m)
{
  ofstream f(filename);
  f.precision(20);
  f << m.Nrows() << " " << m.Ncols() << endl;
  for (int i = 0; i < m.Nrows(); ++i)
    for (int j = 0; j < m.Ncols(); ++j)
      f << i << " " << j << " " << m.element(i, j) << endl;
}
inline void loadnewmat(char* filename, Matrix& m)
{
  ifstream f(filename);
  int rows, cols;
  f >> rows  >> cols;
  m.ReSize(rows, cols);
  int dum;
  for (int i = 0; i < m.Nrows(); ++i)
    for (int j = 0; j < m.Ncols(); ++j)
      f >> dum >> dum >> m.element(i, j); 
}
inline void loadbinarynewmat(char* filename, Matrix& m)
{
  FILE* file = fopen(filename, "rb");
  int nrs, ncs;
  fread(&nrs, sizeof(int), 1, file);
  fread(&ncs, sizeof(int), 1, file);
  m.ReSize(nrs, ncs);
  fread(m.Store(), sizeof(double), m.Storage(), file);
}

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, SymmetricMatrix& a, const unsigned int version)
{
  int dim = a.Nrows();
  ar & dim;
  if(dim != a.Nrows())
    a.ReSize(dim);
  for(int i=0;i<a.Storage();++i)
    ar & a.Store()[i];
}

template<class Archive>
void serialize(Archive & ar, DiagonalMatrix& a, const unsigned int version)
{
  int dim = a.Nrows();
  ar & dim;
  if(dim != a.Nrows())
    a.ReSize(dim);
  for(int i=0;i<a.Storage();++i)
    ar & a.Store()[i];
}

template<class Archive>
void serialize(Archive & ar, Matrix& a, const unsigned int version)
{
  int Nrs = a.Nrows();
  int Ncs = a.Ncols();
  ar & Nrs;
  ar & Ncs;

  if(a.Nrows() != Nrs || a.Ncols() != Ncs)
    {
      a.ReSize(Nrs, Ncs);
    }
  for(int i=0;i<a.Storage();++i) {
    ar & a.Store()[i];
  }
}



} }

#endif
