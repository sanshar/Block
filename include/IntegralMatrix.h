/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef INTEGRAL_ARRAY_HEADER
#define INTEGRAL_ARRAY_HEADER
#define WANT_STREAM
#include <iostream>
#include <fstream>
#include <cmath>
//#include "parser.h"
#include <newmat.h>
#include <newmatio.h>
#include <boost/filesystem.hpp>
#include <multiarray.h>
#include <newmatutils.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "blas_calls.h"
#include <boost/algorithm/string.hpp>
#include <string>
#include <algorithm>
#include "input.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "pario.h"
#include "ObjectMatrix.h"

using namespace boost;

namespace SpinAdapted{
  
class OneElectronArray
  {
  private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & rep;
      ar & dim & rhf & dummyZero & bin & Occnum;
    }
#ifdef NO_SYM_2E
    Matrix rep;
#else
    SymmetricMatrix rep;
#endif
    int dim;
    double dummyZero; /* this variable holds the value 0, for use with rhf integrals */
  public:
    bool rhf;
    bool bin; //if the file is binary
    vector<double> Occnum;
  OneElectronArray() : dim(0), dummyZero(0.0), rhf(false), bin(false)
    {}
  OneElectronArray(int n, bool rhf_=false, bool bin_=false):dummyZero(0.0), rhf(rhf_), bin(bin_)
    {
      ReSize(n);
    }

    double& operator()(int i, int j);
    double operator()(int i, int j) const;
    void ReSize(int n);
    
    int NOrbs() const {
      return dim;
    }

#ifdef NO_SYM_2E
      Matrix& GetRepresentation()
        {
          return rep;
        }
#else
        SymmetricMatrix& GetRepresentation()
        {
          return rep;
        }
#endif

    void ReadOccNum(ifstream& inFile);
    
    void ReadFromDumpFile(ifstream& dumpFile, int norbs);

    void DumpToFile(ofstream& dumpFile) const;

    friend ostream& operator<<(ostream& os, const OneElectronArray& integral)
    {
      os << integral.rep;
      return os;
    }
  };

// It is an interface to unify TwoElectronArray in active space and TwoElectronArray between active space and nonactivespace in nevpt2.
class GeneralTwoElectronArray
{
  public:
    virtual double& operator()(int i, int j, int k, int l) =0;
    virtual double operator () (int i, int j, int k, int l) const =0;

};

//class VaTwoElectronArray : public GeneralTwoElectronArray // 2e array, store the two electron integral in Va subspace of nevpt2.
//{
//  private :
//    ObjectMatrix4D
//
//}

class PerturbTwoElectronArray : public GeneralTwoElectronArray // 2e array, store the two electron integral in Va subspace of nevpt2.
{
  std::vector< std::vector< std::vector<std::vector<double> > > > rep;
  int n0;
  int n1;
  int n2;
  int n3;
  double dummyZero; // dummy variable to hold zero for RHF integrals

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & rep;
    ar & n0 & n1 & n2 & n3;
  }

public:
  PerturbTwoElectronArray () : n0 (0), n1 (0), n2 (0), n3(0) {}
  PerturbTwoElectronArray (int m0, int m1, int m2, int m3) : n0 (m0), n1 (m1), n2 (m2), n3 (m3) { ReSize (n0, n1, n2, n3); }
  
  void ReSize (int m0, int m1, int m2, int m3)
  {
    assert(m0%2==0);
    assert(m1%2==0);
    assert(m2%2==0);
    assert(m3%2==0);
    n0 = m0; n1 = m1; n2 = m2; n3 = m3;
    rep.resize (n0/2);
    for (int i = 0; i < n0/2; ++i)
      {
	rep[i].resize (n1/2);
	for (int j = 0; j < n1/2; ++j)
        {
	  rep[i][j].resize (n2/2);
          for( int k =0; k< n2/2; ++k)
	    rep[i][j][k].resize (n3/2);
        }
      }
  }
  
  virtual double& operator() (int i, int j, int k, int l)
  {
    assert(i >= 0 && i < n0 && j >= 0 && j < n1 && k >= 0 && k < n2 && l >= 0 && l < n3);
    bool is_odd_i = (i & 1);
    bool is_odd_j = (j & 1);
    bool is_odd_k = (k & 1);
    bool is_odd_l = (l & 1);

    bool zero = !((is_odd_i == is_odd_k) && (is_odd_j == is_odd_l));
    if (zero) {
      dummyZero = 0.0;
      return dummyZero;
    }

    i=i/2;
    j=j/2;
    k=k/2;
    l=l/2;
    return rep[i][j][k][l];
  }

  virtual double operator() (int i, int j, int k, int l) const
  {
    assert(i >= 0 && i < n0 && j >= 0 && j < n1 && k >= 0 && k < n2 && l >= 0 && l < n3);
    bool is_odd_i = (i & 1);
    bool is_odd_j = (j & 1);
    bool is_odd_k = (k & 1);
    bool is_odd_l = (l & 1);

    bool zero = !((is_odd_i == is_odd_k) && (is_odd_j == is_odd_l));
    if (zero)
      return 0.0;

    i=i/2;
    j=j/2;
    k=k/2;
    l=l/2;
    return rep[i][j][k][l];
  }

  int NDim0 () const { return n0; }
  int NDim1 () const { return n1; }
  int NDim2 () const { return n2; }
  int NDim3 () const { return n3; }
  int Nrows0 () { return n0; }
  int Nrows1 () { return n1; }
  int Nrows2 () { return n2; }
  int Nrows3 () { return n3; }
};

//class PerturbTwoElectronArray : public GeneralTwoElectronArray // 2e array, store the two electron integral in Va subspace of nevpt2.
//{
//  private:
//    ObjectMatrix4D<double> rep;
//  public:
//    void ReSize(int i, int j, int k, int l){ rep.ReSize(i/2,j/2,k/2,l/2) };
//    T& operator() (int i, int j, int k, int l) {
//      assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
//      return rep [i][j][k][l]; }
//    const T& operator() (int i, int j, int k, int l) const { return rep [i][j][k][l]; }
//
//
//    
//}
//
// ***************************************************************************
// class TwoElectronArray
// ***************************************************************************
/*!
@brief Class for real Hermitian two-electron quantities such as two-e integrals and two-particle r.d.m.

- Uses (12|1'2') format.		

- Hermitian symmetry (ij|kl) = (kl|ij) is AUTOMATICALLY enforced in the underlying storage.

- Additional permutational symmetry (appropriate for 2-e integrals, is
used by default e.g. (i(1)j(2)|k(1)l(2)) = (k(1)j(2)|i(1)l(2)) but
can be turned off through the flag permSymm (default is true).

- The convention for indexing and spins is: odd-index = alpha, even-index = beta.
e.g. twoe(0, 2, 0, 2) would be a all-beta type quantity.

- For spin-restricted quantities, the underlying storage (in rep) can be reduced
by setting the rhf flag. However, to the external user the indexing still behaves
as unrestricted. Default is false.
*/
class TwoElectronArray : public GeneralTwoElectronArray // 2e integral, notation (12|12), symmetric matrix
  {
  private:

    SymmetricMatrix rep; //!< underlying storage enforcing real hermitian symmetry.
    int dim;  //!< #spinorbs
    array_2d<int> indexMap;
    double dummyZero; // dummy variable to hold zero for RHF integrals

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & rep;
      ar & dim;
      ar & indexMap;
      ar & rhf;
      ar & bin;
      ar & permSymm;
    }


  public:
    // default is unrestrictedPermSymm
    enum TwoEType { unrestrictedPermSymm, restrictedPermSymm,
                    unrestrictedNonPermSymm, restrictedNonPermSymm };
    bool rhf; //!< reduced storage for restricted quantities.
    bool permSymm; //!< enforce additional 2e- permutational symmetry.
    bool bin; //if binary file is to be read

  public:

    TwoElectronArray() : dim(0), rhf(false), permSymm(true), bin(false), dummyZero(0.0)
    {}

    TwoElectronArray(TwoEType twoetype);
    
    explicit TwoElectronArray(int n, TwoEType twoetype)
    {
      *this = TwoElectronArray(twoetype);
      ReSize(n);
    }

    int NOrbs() const
      {
        return dim;
      }
    array_2d<int>& GetMap()
    {
      return indexMap;
    }
    SymmetricMatrix& GetRepresentation()
    {
      return rep;
    }

    void ReSize(int n);

    int MapIndices(int n);

    virtual void Load(std::string prefix, int index) {}
    virtual void Save(std::string prefix, int index) {}
    virtual double& operator()(int i, int j, int k, int l);
    virtual double operator () (int i, int j, int k, int l) const;

    void ReadFromDumpFile(ifstream& dumpFile, int norbs);
    
    void DumpToFile(ofstream& dumpFile) const;

    friend ostream& operator<<(ostream& os, const TwoElectronArray& integral)
    {
      os << integral.rep;
      return os;
    }
  };



class PartialTwoElectronArray : public TwoElectronArray // 2e integral, notation (OrbIndex2|12), where OrbIndex is given
{
  private:

    std::vector<std::vector<double> > rep; //!< underlying storage 
    std::vector<int> OrbIndex; //this is the orb i
    int dim;  //!< #spinorbs
    double dummyZero;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & rep & dim & OrbIndex;
    }


  public:

    PartialTwoElectronArray(std::vector<int>& pOrbIndex) : OrbIndex(pOrbIndex)
    {}

    PartialTwoElectronArray(std::vector<int>& pOrbIndex, int pdim) : OrbIndex(pOrbIndex), dim(pdim)
    {}

    void populate(TwoElectronArray& v2);

    void setOrbIndex(std::vector<int>& pOrbIndex) {OrbIndex = pOrbIndex;}

    void Load(std::string prefix, int index=0);

    void Save(std::string prefix, int index=0);

    double operator () (int i, int j, int k, int l) const;
 
    double& operator () (int i, int j, int k, int l);

  };

// ****************************************************
// class PairArray
// ****************************************************
/*!
@brief Class for electron-pairing quantities such as pairing potential and pairing matrix in particle number nonconserving calculations. Only singlet pairing is included

- pairing part in Hamiltonian is given by \Delta_{ij}a_{ia}^\dagger a_{jb}^\dagger+c.c.

- pairing potential \Delta can be symmetric or nonsymmetric, depends on spin-restricted or unrestricted calculation

- default is unrestricted

*/


class PairArray {
  private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive &ar, const unsigned int version) {
      ar & rep;
      ar & srep;
      ar & dim;
      ar & rhf;
      ar & bin;
    }
    Matrix rep;
    SymmetricMatrix srep;
    int dim;
    double dummyZero;
  public:
    bool rhf;
    bool bin;

  public:
    PairArray(): dim(0), rhf(false), bin(false), dummyZero(0.0) {}
    PairArray(int n, bool rhf_=false, bool bin_=false): rhf(rhf_), bin(bin_), dummyZero(0.0) {
      ReSize(n);
    }

    double operator() (int i, int j) const;    
    double& operator() (int i, int j);
    
    void ReSize(int n);

    int NOrbs() const {
      return dim;
    }

    Matrix& GetRepresentation() {
      cerr << "PairArray::GetRepresentation not implemented yet!";
      abort();
    }

    void ReadFromDumpFile(ifstream& dumpFile, int norbs) {
      cerr << "PairArray::ReadFromDumpFile not implemented yet!";
      abort();
    }
    
    void DumpToFile(ofstream& dumpFile) const {
      cerr << "PairArray::DumpToFile not implemented yet!";
      abort();
    }

    friend ostream& operator<<(ostream& os, const PairArray& integral)
    {
      if (integral.rhf) {
        os << integral.srep;        
      } else {
        os << integral.rep;
      }
      return os;
    }
  };

// ****************************************************
// class CCCCArray
// ****************************************************
/*!
@brief Class for CCCC type quantities

- The Hamiltonian is given by
    0.25*w_{ijkl}C_{ia}C_{ja}C_{kb}C_{lb} + c.c.
  where C means creation operator

- Antisymmetry is enforced:
    w_{ijkl}=-w_{jikl}=-w_{ijlk}=w_{jilk}
  therefore, the underlying storage has i>j, l>k, composed indice ij=i*(i-1)/2+j, lk=l*(l-1)/2+k, and the storage is a matrix with (ij, lk) as indice

- Spin symmetry is controled by rhf option. if rhf = true, additionaly
    w_{ijkl}=w_{lkji}
  the underlying storage bocomes a symmetric matrix wrt ij and lk
*/


class CCCCArray {
  private:
    SymmetricMatrix srep;
    Matrix rep;
    int dim; // dim is number of spin orbitals, i.e. 2*norbs
    array_2d<int> indexMap;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, const unsigned int version)
    {
      ar & srep;
      ar & rep;
      ar & dim;
      ar & indexMap;
      ar & rhf;
      ar & bin;
    }

  public:
    bool rhf;
    bool bin;

    // constructors
    CCCCArray(): dim(0), rhf(false), bin(false) {}
    CCCCArray(bool _rhf): dim(0), rhf(_rhf), bin(false) {}
    explicit CCCCArray(int n, bool _rhf) {
      *this = CCCCArray(_rhf);
      ReSize(n);
    }

    int NOrbs() const { return dim;}
    array_2d<int>& GetMap() {  return indexMap;}
    Matrix& GetRepresentation() {
      cerr << "CCCCArray::GetRepresentation not implemented yet!" << endl;
      abort();
    }

    void ReSize(int n);

    int MapIndices(int n);
    
    virtual double operator() (int i, int j, int k, int l) const;
    virtual void set(int i, int j, int k, int l, double value);
    

    void ReadFromDumpFile(ifstream& dumpFile, int norbs) {
      cerr << "CCCCArray::ReadFromDumpFile not implemented yet!";
      abort();
    }
    
    void DumpToFile(ofstream& dumpFile) const {
      cerr << "CCCCArray::DumpToFile not implemented yet!";
      abort();
    }
    
    friend ostream& operator<<(ostream& os, const CCCCArray& integral)
        {
          if (integral.rhf) {
            os << integral.srep;
          } else {
            os << integral.rep;
          }
          return os;
        }
      };

// ****************************************************
// class CCCDArray
// ****************************************************
/*!
@brief Class for CCCD type quantities

- The Hamiltonian is given by
    1/6*w_{ijkl}C_iC_jC_kD_l
  where C means creation operator and D means destruction operator

- Antisymmetry is enforced:
    w_{ijkl}=-w_{ikjl}=-w_{jikl}=w_{jkil}=-w_{kjil}=w_{kijl}
  Considering spin, we only store w_{ijkl} where i>j,k,l and spin of i,j,l the same, while k the opposite. Those which cannot be deduced from this type is 0. composed indice ij=i*(i-1)/2+j, kl=k*n+l and the storage is two matrices with (ij, kl) as indice

- Spin symmetry is controled by rhf option. if rhf = true, additionaly
    w_{ia,ja,kb,la}=-w_{ib,jb,ka,lb}
  the underlying storage bocomes only one matrix
*/


class CCCDArray {
  private:
    Matrix repA, repB;
    int dim;
    array_2d<int> indexMap;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, const unsigned int version)
    {
      ar & repA;
      ar & repB;
      ar & dim;
      ar & indexMap;
      ar & rhf;
      ar & bin;
    }
  
  public:
    bool rhf;
    bool bin;

    // constructors
    CCCDArray(): dim(0), rhf(false), bin(false) {}
    CCCDArray(bool _rhf): dim(0), rhf(_rhf), bin(false) {}
    explicit CCCDArray(int n, bool _rhf) {
      *this = CCCDArray(_rhf);
      ReSize(n);
    }

    int NOrbs() const { return dim;}
    array_2d<int>& GetMap() {  return indexMap;}
    Matrix& GetRepresentation() {
      cerr << "CCCDArray::GetRepresentation not implemented yet!";
      abort();
    }

    void ReSize(int n);
    
    int MapIndices(int n);

    virtual double operator() (int i, int j, int k, int l) const;    
    virtual void set(int i, int j, int k, int l, double value);

    void ReadFromDumpFile(ifstream& dumpFile, int norbs) {
      cerr << "CCCDArray::ReadFromDumpFile not implemented yet!";
      abort();
    }
    
    void DumpToFile(ofstream& dumpFile) const {
      cerr << "CCCDArray::DumpToFile not implemented yet!";
      abort();
    }
    
    friend ostream& operator<<(ostream& os, const CCCDArray& integral)
    {
      os << integral.repA;
      if (!integral.rhf) {  
        os << integral.repB;
      }
      return os;
    }
  };

 
//  enum OnePerturbType{ Vi_1=0, Va_1, Vai_1, OnePertEnd};
//  enum TwoPerturbType{ Vij=0, Vi, Vab, Va, Vabij, Vai, Vabi, Vaij,TwoPertEnd};
} // namespace

#endif
