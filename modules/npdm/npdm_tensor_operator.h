/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_TENSOR_OPERATOR_H
#define NPDM_TENSOR_OPERATOR_H

#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include <stdio.h>
//#include "Symmetry.h"
#include "new_anglib.h"
#include <stdlib.h>
#include <cmath>
#include "pario.h"

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

using namespace std;
//using namespace SpinAdapted;

struct IrrepSpace {
  int irrep;
IrrepSpace() : irrep(0) {};
  IrrepSpace(int p) : irrep(p) {};
  int getirrep() {return 0;}
  friend IrrepSpace operator-(IrrepSpace lhs) {return lhs;} 
  friend std::vector<IrrepSpace> operator+(IrrepSpace lhs, IrrepSpace rhs)
  {
    return std::vector<IrrepSpace>(1);
  }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//FIXME
//namespace NpdmSymmetry {
//  int sizeofIrrep(int i) {return 1;}
//};


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class TensorOp {

 public:
  int Spin; //total spin
  int irrep;
  int rows;
  bool empty;
  std::vector< std::vector<double> > Szops; //there are 2*Spin+1 Szops, which is the outer vector
  //the inner vector is the coefficients for each product of indices eg, [c(ij), c(ij+1), c(i+1j), c(i+1j+1)]
  //where i represents a single creation or destruction op with sz=-1/2 and i+1 with sz=1/2  

  std::vector< std::vector<int> > opindices; //instead of chars it stores the index order of operators
  std::vector<int> optypes; //cre corresponds to 1 and des corresponds to -1
 public:
//FIXME
  int npdm_sizeofIrrep(int i) {return 1;}
  const std::vector< std::vector<double> >& get_Szops() const {return Szops;}
  TensorOp():Spin(0), irrep(0), empty(true) {}
  TensorOp(int k, int sign) :empty(false) 
  {
    int K = 2*k; //convert spatial id to spin id because slaters need that
    int Kirrep = 0;

    if (npdm_sizeofIrrep(Kirrep)>1 && sign <0) {
      int ind1[] = {K+3, K+2, K+1, K+0};
      *this = TensorOp(-1, Kirrep, vector<int>(ind1, ind1+4));
    }
    else if (npdm_sizeofIrrep(Kirrep)==1 && sign <0){
      int ind1[] = {K+1, K+0};
      *this = TensorOp(-1, Kirrep, vector<int>(ind1, ind1+2));
    }
    else if (npdm_sizeofIrrep(Kirrep)>1 && sign > 0) {
      int ind1[] = {K+0, K+1, K+2, K+3};
      *this = TensorOp(1, Kirrep, vector<int>(ind1, ind1+4));
    }
    else {
      int ind1[] = {K+0, K+1};
      *this = TensorOp(1, Kirrep, vector<int>(ind1, ind1+2));
    }

  }
 
 TensorOp(int sgn, int pirrep, std::vector<int> ind) :empty(false)
  {
    if (!(sgn == 1 || sgn == -1)) {
      perr<<"sign not correct";exit(0);}
    Spin = 1;
    if (sgn == -1)
      irrep = (-IrrepSpace(pirrep)).getirrep();
    else
      irrep = pirrep;
    optypes.push_back(sgn);
    if ( npdm_sizeofIrrep(irrep) == 1) {
      std::vector<double> vec1(2,0);
      std::vector<double> vec2(2,0);
      vec1[0] = sgn*1.0; vec2[1] = 1.0;
      Szops.push_back(vec1); Szops.push_back(vec2);
      opindices.resize(2);
      opindices[0].push_back(ind[0]); opindices[1].push_back(ind[1]);
      rows = 1;
    }
    else {
      std::vector<double> vec1(4,0), vec2(4,0), vec3(4,0), vec4(4,0);

      vec1[0] = sgn*1.0; vec2[1] = 1.0; vec3[2] = sgn*1.0; vec4[3] = 1.0;
      Szops.push_back(vec1); Szops.push_back(vec2);
      Szops.push_back(vec3); Szops.push_back(vec4);
      opindices.resize(4);
      opindices[0].push_back(ind[0]); opindices[1].push_back(ind[1]);
      opindices[2].push_back(ind[2]); opindices[3].push_back(ind[3]);
      rows = 2;
    }
  }

  TensorOp& operator*(const double d)
  {
    for (int i=0; i<Szops.size(); i++)
      for (int j=0; j<Szops[i].size(); j++)
	Szops[i][j] *= d;
    return *this;
  }

    
  TensorOp& product(TensorOp& op1, int pspin, int pirrep, bool identical=false) {

    identical=false;

    if (pspin < abs(Spin - op1.Spin) || pspin > Spin+op1.Spin) {
      perr <<"incorrect spin is chosen"<<endl;
      perr <<"cannot combine "<<Spin<<" and "<<op1.Spin<<" to form "<<pspin<<endl;
    }

    std::vector<IrrepSpace> spinvec = IrrepSpace(irrep)+IrrepSpace(op1.irrep);
    bool allowed=false;
    for (int i=0; i<spinvec.size(); i++) {
      if (spinvec[i].getirrep() == pirrep) {
	allowed = true;
	break;
      }
    }
    if (!allowed) {
      rows=0; Szops.clear(); opindices.clear(); Spin =0; irrep = 0; optypes.clear();
      empty= true;
      return *this;
    }
    
    copy(op1.optypes.begin(), op1.optypes.end(), back_inserter(optypes));

    int newrows = npdm_sizeofIrrep(pirrep);
    
    int length = identical? opindices.size()*(op1.opindices.size()+1)/2 : opindices.size()*(op1.opindices.size());
    std::vector< std::vector<int> > tempopindices(length);
    int index =0;
    for (int i=0; i<opindices.size(); i++) 
      for (int j=0; j<op1.opindices.size(); j++) { 
	if (identical && i>j) continue;
	tempopindices[index] = opindices[i];
	copy(op1.opindices[j].begin(), op1.opindices[j].end(), back_inserter(tempopindices[index]) );
	index++;
      }
    opindices = tempopindices;

    std::vector<double> temp(opindices.size(), 0.0);
    std::vector< std::vector<double> > newSzops(newrows*(pspin+1), temp);
    for (int ilz=0; ilz <newrows; ilz++) 
    for (int sz=pspin; sz> -pspin-1; sz-=2) {
      //int lz0 = templz[ilz];
      std::vector<double>&  iSz = newSzops[ilz*(pspin+1)+(-sz+pspin)/2];

      for (int ilz1=0; ilz1 <rows; ilz1++)
      for (int sz1=Spin; sz1> -Spin-1; sz1-=2)
      for (int ilz2=0; ilz2 <op1.rows; ilz2++)	
      for (int sz2=op1.Spin; sz2> -op1.Spin-1; sz2-=2) {
	//int lz1 = lz[ilz1], lz2 = op1.lz[ilz2];
	std::vector<double>&  iSz1 = Szops[ilz1*(Spin+1)+(-sz1+Spin)/2];
	std::vector<double>&  iSz2 = op1.Szops[ilz2*(op1.Spin+1)+(-sz2+op1.Spin)/2];

	//double cleb = cleb_(Spin, sz1, op1.Spin, sz2, pspin, sz);
	double cleb = clebsch(Spin, sz1, op1.Spin, sz2, pspin, sz);
	double clebdinfh = 1.0;//npdm_spatial_cg(irrep, op1.irrep, pirrep, ilz1, ilz2, ilz);
	if (fabs(cleb) <= 1.0e-14 || fabs(clebdinfh) <= 1.0e-14)
	  continue;

	if (!identical) {
	  for (int i=0; i<iSz1.size(); i++)
	    for (int j=0; j<iSz2.size(); j++) {
	      iSz.at(i*iSz2.size()+j) += cleb*clebdinfh*iSz1.at(i)*iSz2.at(j);
	    }
	}
	else {
	  
	  for (int i=0; i<iSz1.size(); i++)
	    for (int j=0; j<iSz2.size(); j++) {
	      double increment = pspin == 0 ? cleb*clebdinfh*iSz1.at(i)*iSz2.at(j)/sqrt(2.0) : 0.0;
	      if (i <= j) 
		iSz.at( i*(i+1)/2 + (j-i) ) += increment;
	      else
		iSz.at( j*(j+1)/2 + (i-j) ) -= increment;//cleb*clebdinfh*iSz1.at(i)*iSz2.at(j)/sqrt(2.0);
	    }
	} 
      }
    }

    this->Szops = newSzops;
    this->Spin = pspin;
    this->rows = newrows;
    this->irrep = pirrep;
    return *this;
  }    
    
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  friend ostream& operator<<(ostream& os, const TensorOp& op) {
    if (op.Spin%2 == 0) {
      for (int ilz = 0; ilz<op.rows; ilz++){
        //int lz = op.lz[ilz];
        os <<"printing operator with Sz = "<<0<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
        for (int i=0; i<op.Szops[op.Spin/2].size(); i++) {
          if (op.Szops[ilz*(op.Spin+1)+op.Spin/2][i] != 0.0) {
            os <<op.Szops[ilz*(op.Spin+1)+op.Spin/2][i]<<"    ";
            for (int j=0; j<op.opindices[i].size(); j++)
              os<<op.opindices[i][j]<<" ";
            os<<endl;
          }
        }
        os << endl;
      }
    }
    else {
      for (int ilz = 0; ilz<op.rows; ilz++){
        //int lz = op.lz[ilz];
        os <<"printing operator with Sz = "<<-op.Spin<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
        for (int i=0; i<op.Szops[op.Spin].size(); i++) {
          if (op.Szops[ilz*(op.Spin+1)+op.Spin][i] != 0.0) {
            os << op.Szops[ilz*(op.Spin+1)+op.Spin][i]<<"   ";
            for (int j=0; j<op.opindices[i].size(); j++)
              os<<op.opindices[i][j]<<" ";
            os<<endl;
          }
        }
        os << endl;
      }
    }

    return os;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}

#endif
