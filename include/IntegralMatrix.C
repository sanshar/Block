/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <IntegralMatrix.h>
#include "pario.h"

using namespace boost;

double& SpinAdapted::OneElectronArray::operator()(int i, int j) {
  assert(i >= 0 && i < dim && j >= 0 && j < dim);
  if (rhf) {
      bool is_odd_i = (i & 1);
      bool is_odd_j = (j & 1);

      bool zero = !(is_odd_i == is_odd_j);
      if (zero) {
        dummyZero = 0.0;
        return dummyZero;
      }

      i=i/2;
      j=j/2;
  }
  return rep(i + 1, j + 1);
}

double SpinAdapted::OneElectronArray::operator()(int i, int j) const {
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
}

void SpinAdapted::OneElectronArray::ReSize(int n) {
  dim = n;
  // dim is the dimension of the spinorbital-space
  if (rhf) {
#ifdef NO_SYM_2E
    rep.ReSize(n/2, n/2);
#else
    rep.ReSize(n/2);
#endif
  } else {
#ifdef NO_SYM_2E
  rep.ReSize(n, n);
#else
  rep.ReSize(n);
#endif
  }
  rep = 0.;
}

void SpinAdapted::OneElectronArray::ReadOccNum(ifstream& inFile) {
  int n;
  inFile >> n;
  if (rhf) {
    n=2*n;
  }
  Occnum.resize(n);
  double occ;
  int i=0;
  while (inFile >> occ ) {
    Occnum.at(2*i) = occ;
    Occnum.at(2*i+1) = occ;
    i++;
  }
}

void SpinAdapted::OneElectronArray::ReadFromDumpFile(ifstream& dumpFile, int norbs) {
  if (bin) {
    dumpFile.seekg (0, ios::end);
	//int size = dumpFile.tellg();
	double size = dumpFile.tellg();
        dumpFile.seekg (0, ios::beg);
        FORTINT nmo = rhf ? static_cast<int>(2*sqrt(size / (sizeof(double)))) : static_cast<int>(sqrt(size / (sizeof(double))));
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
	  dgemm_ (&n, &n, &nmo, &nmo, &nmo, &scale, Aoints.Store(), &nmo, CoeffMatrix.Store (), &nmo, &cfactor, inter, &nmo);
	  dgemm_ (&t, &n, &nmo, &nmo, &nmo, &scale, CoeffMatrix.Store (), &nmo, inter, &nmo, &cfactor, Moints.Store (), &nmo);
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
	  perr << "The first line of one electron integral file should be number of orbitals"<<endl;
	  perr << "Error at line :"<<msg<<endl;
	  abort();
	}
	if (atoi(tok[0].c_str()) != norbs) {
	  perr << "Number of orbitals in one electron integral file should be equal to one given in input file"<<endl;
	  perr << "# orbs in input file : "<<norbs<<endl;
	  perr << "# orbs in one electron integral file : "<<atoi(tok[0].c_str())/2<<endl;
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
	    perr<< "The format of one electron integral file incorrect"<<endl;
	    perr <<"error at this line: "<<msg<<endl;
	    abort();
	  }
	  i = atoi(tok[0].c_str());
	  j = atoi(tok[1].c_str());
	  if (i >= n || j >= n) {
	    perr << "index of orbitals in one electron integral file cannot be bigger than "<<n<<endl;
	    perr<< "error at this line: "<<msg<<endl;
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

void SpinAdapted::OneElectronArray::DumpToFile(ofstream& dumpFile) const {
  dumpFile.precision(20);
  dumpFile << dim << endl;
  for (int i=0; i<dim; ++i)
    for (int j=0; j<dim; ++j) {
      if (fabs((*this)(i, j)) > 1.e-20)
        dumpFile << i << " " << j << " " << (*this)(i, j) << endl;
    }
}

SpinAdapted::TwoElectronArray::TwoElectronArray(TwoEType twoetype) : dummyZero(0.0) {
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

void SpinAdapted::TwoElectronArray::ReSize(int n) {
  //pout << "resizing 2e by n " << n << endl;
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

int SpinAdapted::TwoElectronArray::MapIndices(int n) {
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

double& SpinAdapted::TwoElectronArray::operator()(int i, int j, int k, int l) {
  assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
  if (rhf)
    {
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
    }
  int n = indexMap(i, k);
  int m = indexMap(j, l);
  return rep(n + 1, m + 1);
}

double SpinAdapted::TwoElectronArray::operator()(int i, int j, int k, int l) const {
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

void SpinAdapted::TwoElectronArray::ReadFromDumpFile(ifstream& dumpFile, int norbs) {
  int n = 0;
  string msg; int msgsize = 5000;
  Input::ReadMeaningfulLine(dumpFile, msg, msgsize);
  vector<string> tok;
  boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
  if (tok.size() != 1) {
    perr << "The first line of two electron integral file should be number of orbitals"<<endl;
    perr << "Error at line :"<<msg<<endl;
    abort();
  }
  if (atoi(tok[0].c_str()) != norbs) {
    perr << "Number of orbitals in two electron integral file should be equal to one given in input file"<<endl;
    perr << "# orbs in input file : "<<norbs<<endl;
    perr << "# orbs in two electron integral file : "<<atoi(tok[0].c_str())<<endl;
    abort();
  }

  n = norbs;

  if (rhf) {
    n=2*n;
  }
  ReSize(n);  // this resizes the integral matrix with spin-convention
  int a, b, c, d;
  Input::ReadMeaningfulLine(dumpFile, msg, msgsize);
  while (msg.size() != 0) {
      boost::split(tok, msg, is_any_of(" \t"), token_compress_on);
    if (tok.size() != 5) {
      perr<< "The format of two electron integral file incorrect"<<endl;
      perr <<"error at this line: "<<msg<<endl;
      abort();
    }
    a = atoi(tok[0].c_str());
    b = atoi(tok[1].c_str());
    c = atoi(tok[2].c_str());
    d = atoi(tok[3].c_str());
    if (a >= n || b >= n || c >= n || d >= n) {
      perr << "index of orbitals in two electron integral file cannot be bigger than "<<n<<endl;
      perr<< "error at this line: "<<msg<<endl;
      abort();
    }
    
    if (rhf) {
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
  
void SpinAdapted::TwoElectronArray::DumpToFile(ofstream& dumpFile) const {
    dumpFile.precision(20);
    dumpFile << dim << endl;
    for (int i=0; i<dim; ++i)
      for (int j=0; j<dim; ++j)
        for (int k=0; k<dim; ++k)
          for (int l=0; l<dim; ++l)
            if (fabs((*this)(i, j, k, l)) > 1.e-20)
              dumpFile << i << " " << j << " " << k << " " << l << " " << (*this)(i, j, k, l) << endl;
}

void SpinAdapted::PartialTwoElectronArray::populate(TwoElectronArray& v2) {
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

void SpinAdapted::PartialTwoElectronArray::Load(std::string prefix, int index) {
  char file [5000];
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", prefix.c_str(), "/integral-", OrbIndex[0],"-",OrbIndex[OrbIndex.size()-1], ".", index, ".", mpigetrank(), ".tmp");
  if(mpigetrank() == 0) {
    pout << "\t\t\t Reading Integral file "<<file <<endl;
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_integral(ifs);
    load_integral >> *this ;
  }
}

void SpinAdapted::PartialTwoElectronArray::Save(std::string prefix, int index) {
  char file [5000];
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", prefix.c_str(), "/integral-", OrbIndex[0],"-",OrbIndex[OrbIndex.size()-1], ".", index, ".", mpigetrank(), ".tmp");
  if(mpigetrank() == 0) {
    p1out << "\t\t\t Saving Integral file "<<file <<endl;
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save_integral(ofs);
    save_integral << *this ;
  }
}

double SpinAdapted::PartialTwoElectronArray::operator () (int i, int j, int k, int l) const {
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
	pout << "OrbIndex = "<< OrbIndex[0]<< " does not match any of the indices "<<i<<"  "<<j<<"  "<<k<<"  "<<l<<endl;
	abort();
}

double& SpinAdapted::PartialTwoElectronArray::operator () (int i, int j, int k, int l) {
  assert(i >= 0 && i < 2*dim && j >= 0 && j < 2*dim && k >= 0 && k < 2*dim && l >= 0 && l < 2*dim);
  
  bool is_odd_i = (i & 1);
  bool is_odd_j = (j & 1);
  bool is_odd_k = (k & 1);
  bool is_odd_l = (l & 1);
  
  bool zero = !((is_odd_i == is_odd_k) && (is_odd_j == is_odd_l));
  if (zero) {
    dummyZero = 0.;
    return dummyZero;
  }
  
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
  pout << "OrbIndex = "<< OrbIndex[0]<< " does not match any of the indices "<<i<<"  "<<j<<"  "<<k<<"  "<<l<<endl;
  abort();
}

double& SpinAdapted::PairArray::operator() (int i, int j) {
  assert(i >= 0 && i < dim && j >= 0 && j < dim);
  bool is_odd_i = (i & 1);
  bool is_odd_j = (j & 1);
  bool zero = is_odd_i || (!is_odd_j);
  if (zero) {
    dummyZero = 0.;
    return dummyZero;
  }
  i=i/2;
  j=j/2;
  if (rhf) {
    return srep(i+1, j+1);
  } else {
    return rep(i+1, j+1);
  }
}

double SpinAdapted::PairArray::operator() (int i, int j) const {
  assert(i >= 0 && i < dim && j >= 0 && j < dim);
  bool is_odd_i = (i & 1);
  bool is_odd_j = (j & 1);
  bool zero = is_odd_i || (!is_odd_j);
  if (zero)
    return 0.0;
  i=i/2;
  j=j/2;
  if (rhf) {
    return srep(i+1, j+1);
  } else {
    return rep(i+1, j+1);
  }
}

void SpinAdapted::PairArray::ReSize(int n) {
  dim = n;
  if (rhf) {
    srep.ReSize(n/2);
    srep = 0.;
  } else {
    rep.ReSize(n/2, n/2);
    rep = 0.;
  }
}

void SpinAdapted::CCCCArray::ReSize(int n) {
  dim = n;
  int matDim;
  indexMap.resize(n/2, n/2);
  matDim = MapIndices(n/2);
  if (rhf) {
    srep.ReSize(matDim);
    srep = 0.;
  } else {
    rep.ReSize(matDim, matDim);
    rep = 0.;
  }
}

int SpinAdapted::CCCCArray::MapIndices(int n) {
  int count = 0;
  for (int i=0; i<n; ++i) {
    indexMap(i, i) = 0;
    for (int j=0; j<i; ++j) {
      ++count;
      indexMap(i, j) = count;
      indexMap(j, i) = -count;
    }
  }
  return count;
}

void SpinAdapted::CCCCArray::set(int i, int j, int k, int l, double value) {
  assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
  bool is_odd_i = (i & 1);
  bool is_odd_j = (j & 1);
  bool is_odd_k = (k & 1);
  bool is_odd_l = (l & 1);
  
  bool zero = is_odd_i || is_odd_j || (!is_odd_k) || (!is_odd_l) || (i==j) || (k==l);
  if (zero) {
    return;
  }

  i=i/2;
  j=j/2;
  k=k/2;
  l=l/2;

  int n = indexMap(i, j);
  int m = indexMap(l, k);
  
  int sign = (n>0) == (m>0) ? 1:-1;

  if (rhf) {
    srep(abs(n), abs(m)) = sign * value;
  } else {
    rep(abs(n), abs(m)) = sign * value;
  }
}

double SpinAdapted::CCCCArray::operator()(int i, int j, int k, int l) const {
  assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);

  int spin_i = 1 - (i % 2)*2;
  int spin_j = 1 - (j % 2)*2;
  int spin_k = 1 - (k % 2)*2;
  int spin_l = 1 - (l % 2)*2;
  bool zero = (spin_i+spin_j+spin_k+spin_l != 0);
  if (zero) {
    return 0.;
  }

  i=i/2;
  j=j/2;
  k=k/2;
  l=l/2;

  int factor = 1.;
  int n = 0, m = 0;
  if (spin_i == spin_j) {
    n = indexMap(i, j);
    m = indexMap(l, k);
  } else if (spin_i == spin_k) {
    factor = -1;
    n = indexMap(i, k);
    m = indexMap(l, j);
  } else {
    n = indexMap(i, l);
    m = indexMap(k, j);
  }
  
  if (n == 0 || m == 0) {
    return 0.;
  }

  if (spin_i == -1) {
    int temp = n;
    n = m;
    m = temp;
  }

  factor *= (n>0) == (m>0) ? 1:-1;
  if (rhf) {
    return factor * srep(abs(n), abs(m));
  } else {
    return factor * rep(abs(n), abs(m));
  }
}

void SpinAdapted::CCCDArray::ReSize(int n) {
  dim = n;
  int matDim;
  indexMap.resize(n/2, n/2);
  matDim = MapIndices(n/2);
  repA.ReSize(matDim, dim*dim/4);
  repA = 0;
  if (!rhf) {
    repB.ReSize(matDim, dim*dim/4);
    repB = 0.;        
  }
}

int SpinAdapted::CCCDArray::MapIndices(int n) {
  int count = 0;
  for (int i=0; i<n; ++i) {
    indexMap(i, i) = 0;
    for (int j=0; j<i; ++j) {
      ++count;
      indexMap(i, j) = count;
      indexMap(j, i) = -count;
    }
  }
  return count;
}

void SpinAdapted::CCCDArray::set(int i, int j, int k, int l, double value) {
  assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
  bool is_odd_i = (i & 1);
  bool is_odd_j = (j & 1);
  bool is_odd_k = (k & 1);
  bool is_odd_l = (l & 1);

  bool zero = !((is_odd_i == is_odd_j) && (is_odd_i == is_odd_l) &&(is_odd_i != is_odd_k)) || (i==j);
  if (zero) {
    return;
  }

  i=i/2;
  j=j/2;
  k=k/2;
  l=l/2;

  int n = indexMap(i, j);
  int m = k*(dim/2)+l+1;

  if (rhf) {
    int sign = ((n>0) == is_odd_i)? -1:1;
    repA(abs(n), m) = sign * value;
  } else {
    int sign = (n>0) ? 1:-1;    
    if (is_odd_i) {
      repB(abs(n), m) = sign * value;
    } else {
      repA(abs(n), m) = sign * value;
    }
  }
}

double SpinAdapted::CCCDArray::operator()(int i, int j, int k, int l) const {
  assert(i >= 0 && i < dim && j >= 0 && j < dim && k >= 0 && k < dim && l >= 0 && l < dim);
  int spin_i = 1 - (i % 2)*2;
  int spin_j = 1 - (j % 2)*2;
  int spin_k = 1 - (k % 2)*2;
  int spin_l = 1 - (l % 2)*2;

  bool zero = (spin_i+spin_j+spin_k-spin_l != 0);
  if (zero) {
    return 0.;
  }
  i=i/2;
  j=j/2;
  k=k/2;
  l=l/2;

  int factor = 1.;
  int idx1, idx2, idx3;
  int spin;
  if (spin_i != spin_l) {
    idx1 = j;idx2 = k; idx3 = i; spin = spin_i;
  } else if (spin_j != spin_l) {
    factor = -1; idx1 = i;idx2 = k; idx3 = j; spin = spin_j;
  } else {
    idx1 = i;idx2 = j; idx3 = k; spin = spin_k;
  }

  int n = indexMap(idx1, idx2);
  if (n == 0) {
    return 0.;
  }
  int m = idx3*(dim/2)+l+1;
  if (rhf) {
    factor *= ((n>0) == ( spin == -1)) ? 1:-1;
    return factor * repA(abs(n), m);
  } else {
    factor *= (n>0) ? 1:-1;
    if (spin == -1) {
      return factor * repA(abs(n), m);
    } else {
      return factor * repB(abs(n), m);
    }
  }
}
