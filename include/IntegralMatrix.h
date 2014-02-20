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

using namespace boost;

namespace SpinAdapted{
  
class OneElectronArray
  {
  private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & rep;
      ar & dim;
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
  OneElectronArray() : dim(0), rhf(false), bin(false), dummyZero(0.0)
    {}
  OneElectronArray(int n, bool rhf_=false, bool bin_=false): rhf(rhf_), bin(bin_), dummyZero(0.0)
    {
      ReSize(n);
    }

    double& operator()(int i, int j)
    {
      assert(i >= 0 && i < dim && j >= 0 && j < dim);
      if (rhf)
        {
          bool is_odd_i = (i & 1);
          bool is_odd_j = (j & 1);

          bool zero = !(is_odd_i == is_odd_j);
          if (zero)
            return dummyZero;

          i=i/2;
          j=j/2;
        }
      return rep(i + 1, j + 1);
    }

    double operator()(int i, int j) const
      {
        assert(i >= 0 && i < dim && j >= 0 && j < dim);
        if (rhf)
          {
            bool is_odd_i = (i & 1);
            bool is_odd_j = (j & 1);

            bool zero = !(is_odd_i == is_odd_j);
            if (zero)
              return 0.0;
            i=i/2;
            j=j/2;
          }
        return rep(i + 1, j + 1);
        return rep(i + 1, j + 1);
      }

    void ReSize(int n)
    {
      dim = n;
      // dim is the dimension of the spinorbital-space
      if (rhf)
        {
#ifdef NO_SYM_2E
          rep.ReSize(n/2, n/2);
#else

          rep.ReSize(n/2);
#endif

        }
      else
        {
#ifdef NO_SYM_2E
          rep.ReSize(n, n);
#else

          rep.ReSize(n);
#endif

        }
      rep = 0.;
    }

    int NOrbs() const
      {
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

    void ReadOccNum(ifstream& inFile)
    {
      int n;
      inFile >> n;
      if (rhf)
	{
	  n=2*n;
	}
      Occnum.resize(n);
      double occ;
      int i=0;
      while (inFile >> occ )
	{
	  Occnum.at(2*i) = occ;
	  Occnum.at(2*i+1) = occ;
	  i++;
	}
    }
    
    void ReadFromDumpFile(ifstream& dumpFile, int norbs)
    {
      if (bin) {
        dumpFile.seekg (0, ios::end);
	//int size = dumpFile.tellg();
	double size = dumpFile.tellg();
        dumpFile.seekg (0, ios::beg);
        int nmo = rhf ? static_cast<int>(2*sqrt(size / (sizeof(double)))) : static_cast<int>(sqrt(size / (sizeof(double))));
	ReSize(nmo);
        if (rhf) nmo /= 2;
        char buffer[nmo*nmo*sizeof(double)] ;
        dumpFile.read(buffer, nmo*nmo*sizeof(double));
	Matrix Aoints, Moints;
	Aoints.ReSize(nmo,nmo);
	Moints.ReSize(nmo, nmo);
	Aoints = 0;
        for(int i=0;i<nmo;++i)
          for(int j=0; j<nmo; ++j)
	    {
	      int a=i,b=j;
	      Aoints(a+1,b+1) = ((double*)(&buffer[(i*nmo +j)*sizeof(double)]))[0];
	    }

	//the above are the ao integrals...need mo integrals
	//first read the mo coefficients
	ifstream moCoeff;
	moCoeff.open("42.0", ios::binary);
       
	if (rhf) {
	  Matrix CoeffMatrix;
	  CoeffMatrix.ReSize(nmo, nmo); 
	  char coeffchars[nmo*nmo*sizeof(double)];
	  moCoeff.read(coeffchars, nmo*nmo*sizeof(double));
	  double* coeffs = ((double*)(coeffchars));

	  for (int i=0; i<nmo; i++)
	    for (int j=0; j<nmo; j++) {
	      CoeffMatrix(i+1,j+1) = coeffs[i*nmo+j];
	    }

	  moCoeff.read(coeffchars, nmo*sizeof(double));
	  double* occnums = ((double*)(coeffchars));
	  Occnum.resize(2*nmo);
	  for (int i=0; i<nmo; i++) {
	    Occnum.at(2*i) = occnums[i];
	    Occnum.at(2*i+1) = occnums[i];
	  }

	  double scale=1.0, cfactor=0.0;
	  double* inter = new double[nmo*nmo];
	  char n='n', t='t';
	  DGEMM(n, n, nmo, nmo, nmo, scale, Aoints.Store(), nmo, CoeffMatrix.Store (), nmo, cfactor, inter, nmo);
	  DGEMM(t, n, nmo, nmo, nmo, scale, CoeffMatrix.Store (), nmo, inter, nmo, cfactor, Moints.Store (), nmo);
	  delete [] inter;


	  for(int i=0;i<nmo;++i)
	    for(int j=0; j<nmo; ++j)
	      {
		int a=i,b=j;
		if (rhf)
		  {
		    a*=2;
		    b*=2;
		  }
		(*this)(a,b) = Moints(a/2+1, b/2+1);
	      }

	}
      }
      else {
	int n = 0;
	string msg; int msgsize = 5000;
	Input::ReadMeaningfulLine(dumpFile, msg, msgsize);
	vector<string> tok;
	boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
	if (tok.size() != 1) {
	  cerr << "The first line of one electron integral file should be number of orbitals"<<endl;
	  cerr << "Error at line :"<<msg<<endl;
	  abort();
	}
	if (atoi(tok[0].c_str()) != norbs) {
	  cerr << "Number of orbitals in one electron integral file should be equal to one given in input file"<<endl;
	  cerr << "# orbs in input file : "<<norbs<<endl;
	  cerr << "# orbs in one electron integral file : "<<atoi(tok[0].c_str())/2<<endl;
	  abort();
	}
	n = norbs;
	
	if (rhf)
	  {
	    n=2*n;
	  }

	ReSize(n);
	int i, j;
	
	Input::ReadMeaningfulLine(dumpFile, msg, msgsize);
	while (msg.size() != 0)
	{
	  boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
	  if (tok.size() != 3) {
	    cerr<< "The format of one electron integral file incorrect"<<endl;
	    cerr <<"error at this line: "<<msg<<endl;
	    abort();
	  }
	  i = atoi(tok[0].c_str());
	  j = atoi(tok[1].c_str());
	  if (i >= n || j >= n) {
	    cerr << "index of orbitals in one electron integral file cannot be bigger than "<<n<<endl;
	    cerr<< "error at this line: "<<msg<<endl;
	    abort();
	  }
	  if (rhf)
	  {
	    i=2*i;
	    j=2*j;
	  }
	  (*this)(i, j) = atof(tok[2].c_str());

	  msg.resize(0);
	  Input::ReadMeaningfulLine(dumpFile, msg, msgsize);

	}
      }
    }

    void DumpToFile(ofstream& dumpFile) const
      {
        dumpFile.precision(20);
        dumpFile << dim << endl;
        for (int i=0; i<dim; ++i)
          for (int j=0; j<dim; ++j)
            {
              if (fabs((*this)(i, j)) > 1.e-20)
                dumpFile << i << " " << j << " " << (*this)(i, j) << endl;
            }
      }

    friend ostream& operator<<(ostream& os, const OneElectronArray& integral)
    {
      os << integral.rep;
      return os;
    }
  };


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
class TwoElectronArray // 2e integral, notation (12|12), symmetric matrix
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

    TwoElectronArray(TwoEType twoetype) : dummyZero(0.0)
    {
      switch (twoetype)
        {
        case unrestrictedPermSymm:
          rhf = false;
          permSymm = true;
          break;
        case restrictedPermSymm:
          rhf = true;
          permSymm = true;
          break;
        case unrestrictedNonPermSymm:
          rhf = false;
          permSymm = false;
          break;
        case restrictedNonPermSymm:
          rhf = true;
          permSymm = false;
          break;
        default:
          abort();
        }
    }

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

    void ReSize(int n)
    {
      //cout << "resizing 2e by n " << n << endl;
      dim = n;
      // dim is the dimension of the spinorbital-space

      int matDim;
      if (rhf)
        {
          indexMap.resize(n/2, n/2);
          matDim = MapIndices(n/2);
        }
      else
        {
          indexMap.resize(n, n);
          matDim = MapIndices(n);
        }
      //    rep.ReSize(matDim, matDim);
      rep.ReSize(matDim);
      rep = 0.;
    }

    int MapIndices(int n)
    {
      int count = 0;
      if (permSymm)
        {
          for (int i=0; i<n; ++i)
            for (int j=0; j<=i; ++j)
              {
                indexMap(i, j) = count;
                indexMap(j, i) = count;
                ++count;
              }
        }
      else
        {
          for (int i=0; i<n; ++i)
            for (int j=0; j<n; ++j)
              {
                indexMap(i, j) = count;
                ++count;
              }
        }

      return count;
    }

    virtual void Load(std::string prefix) {}
    virtual void Save(std::string prefix) {}
    virtual double& operator()(int i, int j, int k, int l)
    {
      assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
      if (rhf)
        {
          bool is_odd_i = (i & 1);
          bool is_odd_j = (j & 1);
          bool is_odd_k = (k & 1);
          bool is_odd_l = (l & 1);

          bool zero = !((is_odd_i == is_odd_k) && (is_odd_j == is_odd_l));
          if (zero)
            return dummyZero;

          i=i/2;
          j=j/2;
          k=k/2;
          l=l/2;
        }
      int n = indexMap(i, k);
      int m = indexMap(j, l);
      return rep(n + 1, m + 1);
    }

    virtual double operator () (int i, int j, int k, int l) const
      {
        assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
        if (rhf)
          {
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
          }
        int n = indexMap(i, k);
        int m = indexMap(j, l);

        return rep(n + 1, m + 1);
      }

    void ReadFromDumpFile(ifstream& dumpFile, int norbs)
    {
      int n = 0;
      string msg; int msgsize = 5000;
      Input::ReadMeaningfulLine(dumpFile, msg, msgsize);
      vector<string> tok;
      boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
	if (tok.size() != 1) {
	  cerr << "The first line of two electron integral file should be number of orbitals"<<endl;
	  cerr << "Error at line :"<<msg<<endl;
	  abort();
	}
	if (atoi(tok[0].c_str()) != norbs) {
	  cerr << "Number of orbitals in two electron integral file should be equal to one given in input file"<<endl;
	  cerr << "# orbs in input file : "<<norbs<<endl;
	  cerr << "# orbs in two electron integral file : "<<atoi(tok[0].c_str())<<endl;
	  abort();
	}

	n = norbs;

      if (rhf)
	{
	  n=2*n;
	}
      ReSize(n);  // this resizes the integral matrix with spin-convention
      int a, b, c, d;
      Input::ReadMeaningfulLine(dumpFile, msg, msgsize);
      while (msg.size() != 0)
      {
	boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
	if (tok.size() != 5) {
	  cerr<< "The format of two electron integral file incorrect"<<endl;
	  cerr <<"error at this line: "<<msg<<endl;
	  abort();
	}
	a = atoi(tok[0].c_str());
	b = atoi(tok[1].c_str());
	c = atoi(tok[2].c_str());
	d = atoi(tok[3].c_str());
	if (a >= n || b >= n || c >= n || d >= n) {
	  cerr << "index of orbitals in two electron integral file cannot be bigger than "<<n<<endl;
	  cerr<< "error at this line: "<<msg<<endl;
	  abort();
	}
	
	if (rhf)
	  {
	    a=2*a;
	    b=2*b;
	    c=2*c;
	    d=2*d;
	  }
	(*this)(a, b, c, d) = atof(tok[4].c_str());
	msg.resize(0);
	Input::ReadMeaningfulLine(dumpFile, msg, msgsize);

      }
      
    }
    void DumpToFile(ofstream& dumpFile) const
      {
        dumpFile.precision(20);
        dumpFile << dim << endl;
        for (int i=0; i<dim; ++i)
          for (int j=0; j<dim; ++j)
            for (int k=0; k<dim; ++k)
              for (int l=0; l<dim; ++l)
                if (fabs((*this)(i, j, k, l)) > 1.e-20)
                  dumpFile << i << " " << j << " " << k << " " << l << " " << (*this)(i, j, k, l) << endl;
      }

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

    void populate(TwoElectronArray& v2)
    {
      dim = v2.NOrbs()/2;
      rep.resize(OrbIndex.size());
      for (int orb = 0; orb<OrbIndex.size(); orb++) {
	rep[orb].resize(dim*dim*dim);
	for (int i=0; i<dim; i++)
	for (int j=0; j<dim; j++)
	for (int k=0; k<dim; k++)
	  rep.at(orb)[i*dim*dim+j*dim+k] = v2(2*OrbIndex[orb], 2*i, 2*j, 2*k);
      }
    }

    void setOrbIndex(std::vector<int>& pOrbIndex) {OrbIndex = pOrbIndex;}

    void Load(std::string prefix)
    {
      char file [5000];
      sprintf (file, "%s%s%d%s%d%s%d%s", prefix.c_str(), "/integral-", OrbIndex[0],"-",OrbIndex[OrbIndex.size()-1], ".", mpigetrank(), ".tmp");
      if(mpigetrank() == 0) {
	cout << "\t\t\t Reading Integral file "<<file <<endl;
	std::ifstream ifs(file, std::ios::binary);
	boost::archive::binary_iarchive load_integral(ifs);
	load_integral >> *this ;
      }
    }

    void Save(std::string prefix)
    {
      char file [5000];
      sprintf (file, "%s%s%d%s%d%s%d%s", prefix.c_str(), "/integral-", OrbIndex[0],"-",OrbIndex[OrbIndex.size()-1], ".", mpigetrank(), ".tmp");
      if(mpigetrank() == 0) {
	cout << "\t\t\t Saving Integral file "<<file <<endl;
	std::ofstream ofs(file, std::ios::binary);
	boost::archive::binary_oarchive save_integral(ofs);
	save_integral << *this ;
      }
    }



    double operator () (int i, int j, int k, int l) const
      {
        assert(i >= 0 && i < 2*dim && j >= 0 && j < 2*dim && k >= 0 && k < 2*dim && l >= 0 && l < 2*dim);
	
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

	for (int orb=0; orb<OrbIndex.size(); orb++) {
	  
	  if (i == OrbIndex[orb])
	    return rep[orb][j*dim*dim+k*dim+l];
	  else if (j == OrbIndex[orb])
	    return rep[orb][i*dim*dim+l*dim+k];
	  else if (k == OrbIndex[orb])
	    return rep[orb][l*dim*dim+i*dim+j];
	  else if (l == OrbIndex[orb])
	    return rep[orb][k*dim*dim+j*dim+i];
	}
	cout << "OrbIndex = "<< OrbIndex[0]<< " does not match any of the indices "<<i<<"  "<<j<<"  "<<k<<"  "<<l<<endl;
	abort();
      }

    double& operator () (int i, int j, int k, int l) 
      {
        assert(i >= 0 && i < 2*dim && j >= 0 && j < 2*dim && k >= 0 && k < 2*dim && l >= 0 && l < 2*dim);
	
	bool is_odd_i = (i & 1);
	bool is_odd_j = (j & 1);
	bool is_odd_k = (k & 1);
	bool is_odd_l = (l & 1);
	
	bool zero = !((is_odd_i == is_odd_k) && (is_odd_j == is_odd_l));
	if (zero)
	  return dummyZero;
	
	i=i/2;
	j=j/2;
	k=k/2;
	l=l/2;

	for (int orb=0; orb<OrbIndex.size(); orb++) {
	  
	  if (i == OrbIndex[orb])
	    return rep[orb][j*dim*dim+k*dim+l];
	  else if (j == OrbIndex[orb])
	    return rep[orb][i*dim*dim+l*dim+k];
	  else if (k == OrbIndex[orb])
	    return rep[orb][l*dim*dim+i*dim+j];
	  else if (l == OrbIndex[orb])
	    return rep[orb][k*dim*dim+j*dim+i];
	}
	cout << "OrbIndex = "<< OrbIndex[0]<< " does not match any of the indices "<<i<<"  "<<j<<"  "<<k<<"  "<<l<<endl;
	abort();
      }

  };


}

#endif
