/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "Operators.h"
#include "csf.h"
#include "spinblock.h"
#include "couplingCoeffs.h"
#include "operatorfunctions.h"
#include "opxop.h"
#include "operatorloops.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "SpinQuantum.h"
#include "pario.h"

//using namespace SpinAdapted::operatorfunctions;

bool SpinAdapted::SparseMatrix::nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb)
{
  if (!dmrginp.spinAdapted()) {
    double clebsp = Symmetry::spatial_cg(ladder.sym_is().getirrep(), opsym.get_symm().getirrep(), c1.sym_is().getirrep(), ladder.row(), 0, c1.row());    
    if(c1.n_is() == ladder.n_is() + opsym.get_n() && c1.S.getirrep() == ladder.S.getirrep()+opsym.get_s().getirrep() &&
       fabs(clebsp) >=1.0e-14) {
      nonzeroindex = 0;
      cleb = 1.0;
      return true;
    }
    else
      return false;
  }
  nonzeroindex = 0;
  cleb = 0.0;
  bool found = false;
  int spin = opsym.get_s().getirrep();

  for (int Lz = 0; Lz< Symmetry::sizeofIrrep(opsym.get_symm().getirrep())&&!found; Lz++)
  for (int sz=spin; sz>-spin-1&&!found; sz-=2) 
  {
    cleb = Symmetry::spatial_cg(ladder.sym_is().getirrep(), opsym.get_symm().getirrep(), c1.sym_is().getirrep(), ladder.row(), Lz, c1.row()); 
    cleb *= cg(ladder.S.getirrep(), spin, c1.S.getirrep(), ladder.Sz, sz, c1.Sz);
    if (fabs(cleb) >= 1.0e-14) 
      found = true;
    else
      nonzeroindex++;
  }
  return found;
}

std::vector<double> SpinAdapted::SparseMatrix::calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2)
{
  vector< vector<double> >& szops = Top.Szops; 
  vector<int>& optypes = Top.optypes;
  std::vector<double> elements(szops.size(), 0.0);

  for (int isz=0; isz<szops.size(); isz++) 
    for (map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) 
      for (map<Slater, double>::iterator it2 = c2.det_rep.begin(); it2!= c2.det_rep.end(); it2++) {
        vector<double>& sz = szops[isz];
    
        for (int j=0; j<sz.size(); j++) {
          if (fabs(sz[j]) < 1.0e-14)
	      continue;
          vector<int>& Ind = Top.opindices[j]; //Indices of c and d operators
          Slater s1 = (*it1).first, s2 = (*it2).first;
          double d1 = (*it1).second, d2 = (*it2).second;

          for (int k=Ind.size()-1; k>=0; k--) { //go from high to low
	        if (optypes[k] == 1) s2.c(Ind[k]);
	        else s2.d(Ind[k]);
	        if (s2.isempty())
	          break;
          }
          elements[isz] += sz[j]*d1*d2*(s1.trace(s2));
        }
      }

  return elements;
}


double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const TwoElectronArray& v_2, int integralIndex)
{
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
  for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;

    int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
    std::vector<double>&  iSz2 = op2.Szops[sz2index];
    
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == CD) {
	  factor += 0.5*(-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]) - v_2(Ind2[0], Ind1[0], Ind1[1], Ind2[1]) 
	    		 + v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1]) + v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == DD) {
	  factor += 0.5*v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
	else if (comp == CCD) {
	  factor += 0.5*(v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]) - v_2(Ind1[1], Ind1[0], Ind2[0], Ind1[2]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == CDD) {
	  factor += 0.5*(v_2(Ind2[0], Ind1[0], Ind1[2], Ind1[1]) - v_2(Ind1[0], Ind2[0], Ind1[2], Ind1[1]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == C) {
	  if (op1.dn() + op2.dn() == 0) {
	    factor += 0.5*v_1[integralIndex](Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  } 
	  else { // D
	    factor += 0.5*(v_cc(Ind2[0], Ind1[0]) - v_cc(Ind1[0], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  }
	}
      }
  }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const PerturbTwoElectronArray& v_2, int integralIndex)
{
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
  for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;

    int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
    std::vector<double>&  iSz2 = op2.Szops[sz2index];
    
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == CDD_CD) {
	  factor += (v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  //factor +=  (v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1])- v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
  else if (comp == CCD_CD){
	  factor += (v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
  }
	else if (comp == DD) {
	  //factor += 0.5*v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
  else if (comp == CDD) {
	  //factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += (v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
  else if (comp == CCD) {
	  //factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += (v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
	else {
          assert(false);
        }
      }
  }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const OneElectronArray& v_1, int integralIndex)
{
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
  for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;

    int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
    std::vector<double>&  iSz2 = op2.Szops[sz2index];
    
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == C) {
	  factor += v_1(Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else {
          assert(false);
        }
      }
  }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const PerturbTwoElectronArray& v_2, int integralIndex)
{
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, v_2, integralIndex);
  double factor = 0.0;
  vector<double>& iSz2 = op2.Szops[op2index];
  bool found = false;
  for (int ilz1=0; ilz1 <op1.rows; ilz1++)	
  for (int sz1=-op1.Spin; sz1< op1.Spin+1; sz1+=2) {
    if (found) break;
    
    int ilz2 = op2index/(op2.Spin+1);
    int sz2index = (op2index - ilz2*(op2.Spin+1)), sz2 = op2.Spin - 2*sz2index;
    std::vector<double>&  iSz1 = op1.Szops[ilz1*(op1.Spin+1)+(-sz1+op1.Spin)/2];
    
    //double cleb = cleb_(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    double cleb = clebsch(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == CDD_CD) {
	  //factor += 0.5*(-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1])  
	  //  		  + v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  factor +=  (v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1])) *iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
  else if (comp == CCD_CD){
	  factor += (v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
  }
	else if (comp == DD) {
	  //factor += 0.5*v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
  else if (comp == CDD) {
	  factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
  else if (comp == CCD) {
	  //factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += (v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
	else {
          assert(false);
        }
      }
  }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2, int integralIndex)
{
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, v_2, integralIndex);
  double factor = 0.0;
  vector<double>& iSz2 = op2.Szops[op2index];
  bool found = false;
  for (int ilz1=0; ilz1 <op1.rows; ilz1++)	
  for (int sz1=-op1.Spin; sz1< op1.Spin+1; sz1+=2) {
    if (found) break;
    
    int ilz2 = op2index/(op2.Spin+1);
    int sz2index = (op2index - ilz2*(op2.Spin+1)), sz2 = op2.Spin - 2*sz2index;
    std::vector<double>&  iSz1 = op1.Szops[ilz1*(op1.Spin+1)+(-sz1+op1.Spin)/2];
    
    //double cleb = cleb_(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    double cleb = clebsch(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); ++i1)
      for (i2 =0; i2<iSz2.size(); ++i2) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == CD) {
	  factor += (-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]) + 
	    	     v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == DD) {
	  factor += 0.5*(v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
	else if (comp == CCD) {
	  factor += v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == CDD) {
	  factor += 0.5*(v_2(Ind2[0], Ind1[0], Ind1[2], Ind1[1]) - v_2(Ind1[0], Ind2[0], Ind1[2], Ind1[1]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  //factor += 0.5*(v_2(Ind2[0], Ind1[0], Ind1[1], Ind1[2]) - v_2(Ind1[0], Ind2[0], Ind1[1], Ind1[2]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == C) {
	  if (op1.dn() + op2.dn() == 0) {
	    factor += 0.5*v_1[integralIndex](Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  } 
	  else { // D
	    factor += 0.5*(v_cc(Ind2[0], Ind1[0]) - v_cc(Ind1[0], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  }
	}
	else {
          abort();
        }
      }
  }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const CCCDArray& vcccd) {
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
    for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
      if (found) break;
      
      int ilz1 = 0;
      int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
      std::vector<double>&  iSz2 = op2.Szops[sz2index];
      
      double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

      cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
      if (fabs(cleb) <= 1.0e-14)
        continue;
      else 
        found = true;
      int i1, i2;

      for (i1=0; i1<iSz1.size(); ++i1)
        for (i2 =0; i2<iSz2.size(); ++i2) {
          vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	  
          if (comp == CD) {
            if (op2.dn() == 2) {
              factor += 0.5*vcccd(Ind1[0],Ind2[0],Ind2[1],Ind1[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else {
              factor += 0.5*vcccd(Ind2[1],Ind2[0],Ind1[1],Ind1[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else if (comp == DD) {
            if (op1.dn() == 2) {
              factor += 0.5*vcccd(Ind1[0],Ind1[1],Ind2[0],Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else { // CC_comp
              factor += 0.5*vcccd(Ind1[1],Ind1[0],Ind2[1],Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else if (comp == CCD) { // two cases CDD and CCC
            if (op1.dn() == 3) { // CCC
              factor += (1./6) * vcccd(Ind1[0], Ind1[1], Ind1[2], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else { // CDD
              factor += 0.5 * vcccd(Ind2[0], Ind1[2], Ind1[1], Ind1[0]) *iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else if (comp == CDD) {
            if (op1.dn() == -3) {
              factor += (1./6) * vcccd(Ind1[2], Ind1[1], Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else {
              factor += 0.5 * vcccd(Ind2[0], Ind1[0], Ind1[1], Ind1[2])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else {
            abort();
          }
        }
    }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const CCCDArray& vcccd) {
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, vcccd);
  pout << "Sorry, SpinAdapted BCS calculation not implemented" << endl;
  abort();
  return 0.;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const CCCCArray& vcccc) {
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
    for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
      if (found) break;
      
      int ilz1 = 0;
      int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
      std::vector<double>&  iSz2 = op2.Szops[sz2index];
      
      double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

      cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
      if (fabs(cleb) <= 1.0e-14)
        continue;
      else 
        found = true;
      int i1, i2;

      for (i1=0; i1<iSz1.size(); ++i1)
        for (i2 =0; i2<iSz2.size(); ++i2) {
          vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
          if (comp == DD) {
            factor += 0.25 * vcccc(Ind1[0], Ind1[1], Ind2[0], Ind2[1]) *iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == CCD) {
            factor += (1./6) * vcccc(Ind2[0], Ind1[2], Ind1[1], Ind1[0]) * iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == CDD) {
            factor += (1./6) * vcccc(Ind2[0], Ind1[0], Ind1[1], Ind1[2]) * iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else {
            abort();
          }
        }
    }
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const CCCCArray& vcccc) {
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, vcccc);
  pout << "Sorry, SpinAdapted BCS calculation not implemented" << endl;
  abort();
  return 0.;
}

//******************CRE*****************

void SpinAdapted::Cre::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE).has(i))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE, deltaQuantum, i);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix>& op = rightBlock->get_op_rep(CRE, deltaQuantum, i);
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    dmrginp.makeopsT -> stop();
    return;
  }  
  abort();  
}


double SpinAdapted::Cre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int Iirrep = SymmetryOfOrb(get_orbs()[0]).getirrep();;
  IrrepSpace sym = deltaQuantum[0].get_symm();
  bool finish = false;
  bool write = false;
  int Sign = 1;

  TensorOp C(get_orbs()[0], 1);

  for (int j = 0; j < deltaQuantum.size(); ++j) {  
    for (int i=0; i<ladder.size(); i++) {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        std::vector<double> MatElements = calcMatrixElements(c1, C, ladder[i]) ;
        element += MatElements[index]/cleb;
        break;
      }
      else
        continue;
    }
  }

  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::Cre::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<Cre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new Cre);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//******************DES*****************

void SpinAdapted::Des::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES).has(i))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(DES, deltaQuantum, i);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(DES).has(i))
  {
    const boost::shared_ptr<SparseMatrix>& op = rightBlock->get_op_rep(DES, deltaQuantum, i);
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    dmrginp.makeopsT -> stop();
    return;
  }  
  else
    abort();  

}


double SpinAdapted::Des::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int Iirrep = SymmetryOfOrb(get_orbs()[0]).getirrep();;
  IrrepSpace sym = deltaQuantum[0].get_symm();
  bool finish = false;
  bool write = false;
  int Sign = 1;

  TensorOp D(get_orbs()[0], -1);

  for (int iq=0; iq<deltaQuantum.size(); iq++)
  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[iq], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, D, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
    
  }

  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::Des::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<Des>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new Des);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}


//******************CREDES*****************

void SpinAdapted::CreDes::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }

    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE_DES).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    dmrginp.makeopsT -> stop();
    return;
  }  
  if (leftBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    if (rightBlock->has(DES) ) {
      boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(j), j);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0);
    }
    else {
      Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(j), j));
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  }
  else if (rightBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    if (leftBlock->has(DES) ) {
      boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(j), j);
      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
    }
    else {
      Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(j), j));
      double parity = getCommuteParity(op1->get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[0]);
    // getCommuteParity doesn't depend on deltaQuantum.get_n()
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
    }
  }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();

  TensorOp C(I, 1), D(J, -1);
  TensorOp CD = C.product(D, spin, irrep);

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        std::vector<double> MatElements = calcMatrixElements(c1, CD, ladder[i]) ;
        element += MatElements[index]/cleb;
        break;
      }
      else
        continue;
    }
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDes::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreDes);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//******************DESCRE*****************

void SpinAdapted::DesCre::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_CRE).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }

    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(DES_CRE).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    dmrginp.makeopsT -> stop();
    return;
  }  
  if (leftBlock->get_op_array(DES).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op2, *op1, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (rightBlock->get_op_array(DES).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(DES, -getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);    
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::DesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();

  TensorOp D(I, -1), C(J, 1);
  
  TensorOp CD = C.product(D, spin, irrep);


  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        std::vector<double> MatElements = calcMatrixElements(c1, CD, ladder[i]) ;
        element += MatElements[index]/cleb;
        break;
      }
      else
        continue;
    }
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesCre::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<DesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new DesCre);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//******************CRECRE*****************

void SpinAdapted::CreCre::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    dmrginp.makeopsT -> stop();
    return;
  }  
  if (leftBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (rightBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();


  TensorOp C1(I,1), C2(J,1);
  TensorOp CC = C1.product(C2, spin, sym.getirrep(), I==J);

  for (int j = 0; j < deltaQuantum.size(); ++j) {  
  for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        std::vector<double> MatElements = calcMatrixElements(c1, CC, ladder[i]) ;
        element += MatElements[index]/cleb;
        break;
      }
      else
        continue;
    }
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCre::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreCre);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}


//******************DESDES*****************

void SpinAdapted::DesDes::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_DES).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(DES_DES, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(DES_DES).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(DES_DES, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    dmrginp.makeopsT -> stop();
    return;
  }  
  if (leftBlock->get_op_array(DES).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(j), j);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (rightBlock->get_op_array(DES).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(DES, -getSpinQuantum(i), i);
    const boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(j), j);
    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
  }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::DesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();


  TensorOp D1(I,-1), D2(J,-1);
  TensorOp DD = D1.product(D2, spin, sym.getirrep(), I==J);

  for (int j = 0; j < deltaQuantum.size(); ++j) {  
  for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        std::vector<double> MatElements = calcMatrixElements(c1, DD, ladder[i]) ;
        element += MatElements[index]/cleb;
        break;
      }
      else
        continue;
    }
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesDes::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<DesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new DesDes);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}



//******************CREDESCOMP*****************


void SpinAdapted::CreDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp C(i,1), D(j,-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DESCOMP).has(i, j))
  { 
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }

  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DESCOMP).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // build CDcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	    double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      if (rightBlock->has(DES)) {
	        boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	        SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	      } else {
	        Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	        SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
	      }
	    }
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());

      	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      if (rightBlock->has(DES)) {
	        boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	        double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	        SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	      } else {
	        Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	        double parity = getCommuteParity(op1->get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[0]);
	        SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	      }
  	    }
      }

      if (dmrginp.hamiltonian() == BCS) {
        CK = TensorOp(k, 1);
        TensorOp CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l); // k cannot equal to l
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd);
          CL = TensorOp(l, 1);
          CK = TensorOp(k, 1);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, CC2_commute, CD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[1]);
            scaleV += parity * scaleV2;
	        if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        TensorOp DK(k, -1);
        DL = TensorOp(l, -1);
        TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
        if (!DD2.empty) {
          double scaleV = calcCompfactor(CD1, DD2, CD, v_cccd);
          DL = TensorOp(l, -1);
          DK = TensorOp(k, -1);
          TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, DD2_commute, CD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            if (leftBlock->has(DES)) {
              boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
              boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
              double parity =  getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
              scaleV += parity * scaleV2;
              if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
                SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
              }
            } else {
              Transposeview top1 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
              Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
              double parity = getCommuteParity(top1.get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[2]);
              scaleV += parity * scaleV2;
	          if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
                SpinAdapted::operatorfunctions::TensorProduct(leftBlock, top1, top2, &b, &(b.get_stateInfo()), *this, scaleV);            
              }
            }
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is(); // classify whether we calculate CC or CD

  TensorOp C(I,1), D(J,-1);

  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep());

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];
            
            if (dmrginp.hamiltonian() == BCS && dn == 2) {
              TensorOp CK(k,1), CL(l,1);
              TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
              if (!CC2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i]);
                double factor = calcCompfactor(CD1, CC2, CD, v_cccd);
                element += MatElements[index]*factor/cleb;
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == -2) {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
              if (!DD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
                double factor = calcCompfactor(CD1, DD2, CD, v_cccd);
                element += MatElements[index]*factor/cleb;
              }
            } else {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
                double factor = calcCompfactor(CD1, CD2, CD, *(b->get_twoInt()), b->get_integralIndex());
                element += MatElements[index]*factor/cleb;
              }
            }
          }
        break;
      }
      else
        continue;
    }
  }
  return element;
}


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreDesComp);
      *rep = *this;

      rep->build(*block);
      return rep;
    } } //******************DESCRECOMP*****************


void SpinAdapted::DesCreComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp D(i,-1), C(j,1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_CRECOMP).has(i, j))
  { 
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(DES_CRECOMP, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }

  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(DES_CRECOMP).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(DES_CRECOMP, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // build DCcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k, 1), DL(l, -1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	    double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(DES).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	    }
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
      	if (leftBlock->get_op_array(DES).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	    }
      }

      if (dmrginp.hamiltonian() == BCS) {
        CK = TensorOp(k, 1);
        TensorOp CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l); // k cannot equal to l
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd);
          CL = TensorOp(l, 1);
          CK = TensorOp(k, 1);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, CC2_commute, CD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[1]);
            scaleV += parity * scaleV2;
	        if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        TensorOp DK(k, -1);
        DL = TensorOp(l, -1);
        TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
        if (!DD2.empty) {
          double scaleV = calcCompfactor(CD1, DD2, CD, v_cccd);
          DL = TensorOp(l, -1);
          DK = TensorOp(k, -1);
          TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, DD2_commute, CD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
            double parity =  getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
            scaleV += parity * scaleV2;
            if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::DesCreComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is(); // classify whether we calculate CC or CD

  TensorOp D(I,-1), C(J, 1);

  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep());

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];
            
            if (dmrginp.hamiltonian() == BCS && dn == 2) {
              TensorOp CK(k,1), CL(l,1);
              TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
              if (!CC2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i]);
                double factor = calcCompfactor(CD1, CC2, CD, v_cccd);
                element += MatElements[index]*factor/cleb;
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == -2) {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
              if (!DD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
                double factor = calcCompfactor(CD1, DD2, CD, v_cccd);
                element += MatElements[index]*factor/cleb;
              }
            } else {
	          TensorOp CK(k, 1), DL(l, -1);
	          TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
	          if (!CD2.empty) {
	            std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
	            double factor = calcCompfactor(CD1, CD2, CD, *(b->get_twoInt()), b->get_integralIndex());
	            element += MatElements[index]*factor/cleb;  // FIXME a factor of half?
	          }
            }
          }
        break;
      }
      else
        continue;
    }
  }
  return element;
}


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesCreComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<DesCreComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new DesCreComp);
      *rep = *this;

      rep->build(*block);
      return rep;
    }

}

//******************DESDESCOMP*****************

void SpinAdapted::DesDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp C(i,1), C2(j,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_DESCOMP).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(DES_DESCOMP).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // explicitly build DD_comp
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp DK(k,-1), DL(l,-1);
      TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
      if (!DD2.empty) {
        double scaleV = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());

        DK=TensorOp(k,-1); DL=TensorOp(l,-1);
        DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
        double scaleV2 = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());
        
        if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	  if (leftBlock->has(DES)) {
	    boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	    
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    scaleV += parity*scaleV2;
	    
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	  }
	  else {
	    Transposeview top1 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	    Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	    
	    double parity = getCommuteParity(top1.get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    scaleV += parity*scaleV2;
	    
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, top1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
	  }
        }
      }

      if (dmrginp.hamiltonian() == BCS) {
        TensorOp CK(k, 1), CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CC1, CC2, DD, v_cccc);
          CK = TensorOp(k, 1);
          CL = TensorOp(l, 1);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CC1, CC2_commute, DD, v_cccc);

          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(2));
            scaleV += parity*scaleV2;

            if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	          SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        // Ck*Dl
        CK = TensorOp(k, 1);
        DL = TensorOp(l, -1);
        TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
        if (!CD2.empty) {
          double scaleV = calcCompfactor(CC1, CD2, DD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            if (rightBlock->has(DES)) {
              boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	          SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            } else {
	          Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	          SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        // Cl*Dk
        CL = TensorOp(l,1);
        DK = TensorOp(k,-1);
        CD2 = CL.product(DK, spin, sym.getirrep());
        if (!CD2.empty) {
          double scaleV = calcCompfactor(CC1, CD2, DD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            if (leftBlock->has(DES)) {
               boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	           double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(1));
	           SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
            } else {
	          Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	          double parity = getCommuteParity(op1->get_deltaQuantum(0), top2.get_deltaQuantum(0), get_deltaQuantum(1));
	          SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
            }
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();

}


double SpinAdapted::DesDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp C(I,1), C2(J,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), I==J);
 
  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {

      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];	
           
            if (dmrginp.hamiltonian() == BCS && dn == 2) {
              TensorOp CK(k,1), CL(l,1);
              TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
              if (!CC2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i]);
                double scale = calcCompfactor(CC1, CC2, DD, v_cccc);
                element += MatElements[index]*scale/cleb;                
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == 0) {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
                double scale = calcCompfactor(CC1, CD2, DD, v_cccd);
                element += MatElements[index]*scale/cleb;
              }
            } else {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);

              if (!DD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
                double scale = calcCompfactor(CC1, DD2, DD, index, *(b->get_twoInt()), b->get_integralIndex());
                element += MatElements[index]*scale/cleb;
              }
            } 
          }
        break;
      }
      else
        continue;

    }
  }
  return element;

}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::DesDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<DesDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new DesDesComp);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}


//******************CRECRECOMP*****************

void SpinAdapted::CreCreComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp D(i, -1), D2(j, -1);
  TensorOp DD1 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRECOMP).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRECOMP, deltaQuantum, i,j);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRECOMP).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_CRECOMP, deltaQuantum, i,j);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  // explicitly build CC_comp
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k, 1), CL(l, 1);
      TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
      if (!CC2.empty) {
        double scaleV = calcCompfactor(CC2, DD1, DD, *(b.get_twoInt()), b.get_integralIndex());

        CK=TensorOp(k, 1); CL=TensorOp(l, 1);
        CC2 = CL.product(CK, spin, sym.getirrep(), k==l);
        double scaleV2 = calcCompfactor(CC2, DD1, DD, *(b.get_twoInt()), b.get_integralIndex());
        
        if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	  boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	  
	  double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	  scaleV += parity*scaleV2;
	  
	  if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
        }
      }

      if (dmrginp.hamiltonian() == BCS) {
        TensorOp DK(k,-1), DL(l,-1);
        TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
        if (!DD2.empty) {
          double scaleV = calcCompfactor(DD1, DD2, DD, v_cccc);
          DK = TensorOp(k,-1);
          DL = TensorOp(l,-1);
          TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(DD1, DD2_commute, DD, v_cccc);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(2));
            scaleV += parity*scaleV2;
            if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	          SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        // Cl*Dk
        CL = TensorOp(l,1);
        DK = TensorOp(k,-1);
        TensorOp CD2 = CL.product(DK, spin, sym.getirrep());
        if (!CD2.empty) {
          double scaleV = calcCompfactor(DD1, CD2, DD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	        double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(1));
	        SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
          }
        }
        // Ck*Dl
        CK = TensorOp(k, 1);
        DL = TensorOp(l, -1);
        CD2 = CK.product(DL, spin, sym.getirrep());
        if (!CD2.empty) {
          double scaleV = calcCompfactor(DD1, CD2, DD, v_cccd);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	        SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();

}


double SpinAdapted::CreCreComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp D(I,-1), D2(J,-1);
  TensorOp DD1 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), I==J);
 
  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {

      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];	
           
            if (dmrginp.hamiltonian() == BCS && dn == -2) {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);

              if (!DD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
                double scale = calcCompfactor(DD1, DD2, DD, v_cccc);
                element += MatElements[index]*scale/cleb;                
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == 0) {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
                double scale = calcCompfactor(DD1, CD2, DD, v_cccd);
                element += MatElements[index]*scale/cleb;                
              }
            } else {
              TensorOp CK(k,1), CL(l,1);
	          TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);

	          if (!CC2.empty) {
	            std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i]);
	            double scale = calcCompfactor(CC2, DD1, DD, index, *(b->get_twoInt()), b->get_integralIndex());
	            element += MatElements[index]*scale/cleb;
	          }
            }
          }
        break;
      }
      else
        continue;

    }
  }
  return element;

}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreCreComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreCreComp);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}




//******************CRECREDESCOMP*****************


void SpinAdapted::CreCreDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int k = get_orbs()[0];


  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  SpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);

  if (leftBlock->get_op_array(CRE_CRE_DESCOMP).has(k)) {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      else {
	//const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE_DESCOMP).has(k))
    {
      const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
      //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  

  // explicit build CCD_comp
  if (dmrginp.hamiltonian() != HUBBARD){
    if (loopBlock->has(CRE_DESCOMP))
      {
      Functor f = boost::bind(&opxop::cxcdcomp, otherBlock, _1, &b, k, this, 1.0); 
      for_all_singlethread(loopBlock->get_op_array(CRE), f);

      f = boost::bind(&opxop::dxcccomp, otherBlock, _1, &b, k, this, 2.0); // factor of 2.0 because CCcomp_{ij} = -CCcomp_{ji}
      for_all_singlethread(loopBlock->get_op_array(CRE), f);

      f = boost::bind(&opxop::cxcdcomp, loopBlock, _1, &b, k, this, 1.0); 
      for_all_singlethread(otherBlock->get_op_array(CRE), f);

      f = boost::bind(&opxop::dxcccomp, loopBlock, _1, &b, k, this, 2.0);
      for_all_singlethread(otherBlock->get_op_array(CRE), f);

    } else if (otherBlock->has(CRE_DESCOMP)) {
      pout << "I should not be here"<<endl;exit(0);
    } 
  }

  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CreCreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int K = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();  

  TensorOp D(K, -1);

  for (int j = 0; j < deltaQuantum.size(); ++j)
  for (int i=0; i<ladder.size(); ++i) {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
      for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++) 
      for (int kl =0; kl<b->get_sites().size(); kl++) {
        int _i = b->get_sites()[ki];
        int _j = b->get_sites()[kj];
        int _l = b->get_sites()[kl];
        SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
        if (dmrginp.hamiltonian() == BCS && dn == 3) { // CCC
          std::vector<SpinQuantum> sij = si+sj;
          for (int ij=0; ij<sij.size(); ++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij+sl;
            for (int ijl=0; ijl<sijl.size();++ijl) {
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp CI(_i, 1), CJ(_j, 1), CL(_l, 1);
              TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
              // FIXME maybe there takes specical care when j = l in spinadapted case
              TensorOp CCCIJL = CCIJ.product(CL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (CCCIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, CCCIJL, ladder[i]);
              double scale = calcCompfactor(CCCIJL, D, CCD, v_cccd);
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }
          }
        } else if (dmrginp.hamiltonian() == BCS && dn == -1) { // CDD
          std::vector<SpinQuantum> sij = si-sj;
          for (int ij=0; ij<sij.size(); ++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij-sl;
            for (int ijl=0; ijl<sijl.size(); ++ijl) {
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp CI(_i, 1), DJ(_j, -1), DL(_l, -1);
              TensorOp CDIJ = CI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep());
              TensorOp CDDIJL = CDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (CDDIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, CDDIJL, ladder[i]);
              double scale = calcCompfactor(CDDIJL, D, CCD, v_cccd);
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }
          }
        } else if (dmrginp.hamiltonian() == BCS && dn == -3) { // DDD
          std::vector<SpinQuantum> sij = (-si)-sj;
          for (int ij=0; ij<sij.size(); ++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij-sl;
            for (int ijl=0; ijl<sijl.size(); ++ijl) {            
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp DI(_i, -1), DJ(_j, -1), DL(_l, -1);
              TensorOp DDIJ = DI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
              TensorOp DDDIJL = DDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (DDDIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, DDDIJL, ladder[i]);
              double scale = calcCompfactor(DDDIJL, D, CCD, v_cccc);
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }
          }
        } else {  // CCD
          std::vector<SpinQuantum> sij = si+sj;
          for (int ij=0; ij<sij.size(); ++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij-sl;
            for (int ijl=0; ijl<sijl.size(); ++ijl) {
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
              TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
              TensorOp CCDIJL = CCIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (CCDIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, CCDIJL, ladder[i]);
              double scale = calcCompfactor(CCDIJL, D, CCD, *(b->get_twoInt()), b->get_integralIndex());
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }
          }
        }
      }

      for (int ki =0; ki<b->get_sites().size(); ki++) {  // add C block to CCD block
        int _i = b->get_sites()[ki];
        if (dmrginp.hamiltonian() == BCS && dn == -1) { // D
          TensorOp DI(_i, -1);
          std::vector<double> MatElements = calcMatrixElements(c1, DI, ladder[i]);
          double factor = calcCompfactor(DI, D, C, *(b->get_twoInt()), b->get_integralIndex());
          if (fabs(factor) > dmrginp.oneindex_screen_tol())
            element += factor*MatElements[index]/cleb;
        } else { // C
          TensorOp CI(_i, 1);
          std::vector<double> MatElements = calcMatrixElements(c1, CI, ladder[i]);
          double factor = calcCompfactor(CI, D, C, *(b->get_twoInt()), b->get_integralIndex());
          if (fabs(factor) > dmrginp.oneindex_screen_tol())
            element += factor*MatElements[index]/cleb;
        }
      }
      break;
    } else
      continue;
  }

  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreCreDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreCreDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreCreDesComp);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}


//******************CREDESDESCOMP*****************


void SpinAdapted::CreDesDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int k = get_orbs()[0];


  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  SpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);

  if (leftBlock->get_op_array(CRE_DES_DESCOMP).has(k))
    {      
      const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DESCOMP, deltaQuantum, k);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      else {
	//const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DES_DESCOMP).has(k))
    {
      const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DES_DESCOMP, deltaQuantum, k);
      //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  

  // explicit build CDD_comp
  if (dmrginp.hamiltonian() != HUBBARD){
    if (loopBlock->has(CRE_DESCOMP))
      {
	Functor f = boost::bind(&opxop::dxcdcomp, otherBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(loopBlock->get_op_array(DES), f);
	
	f = boost::bind(&opxop::cxddcomp, otherBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(loopBlock->get_op_array(CRE), f);
        
	f = boost::bind(&opxop::dxcdcomp, loopBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(otherBlock->get_op_array(DES), f);

	f = boost::bind(&opxop::cxddcomp, loopBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(otherBlock->get_op_array(CRE), f);
      }
    else if (otherBlock->has(CRE_DESCOMP))
      {
	pout << "I should not be here"<<endl;exit(0);
     }
  }
  dmrginp.makeopsT -> stop();


}


double SpinAdapted::CreDesDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int K = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp CK(K, 1);
  for (int j = 0; j<deltaQuantum.size(); ++j)
  for (int i=0; i<ladder.size(); ++i) {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
      for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++) 
      for (int kl =0; kl<b->get_sites().size(); kl++) {
	    int _i = b->get_sites()[ki];
	    int _j = b->get_sites()[kj];
	    int _l = b->get_sites()[kl];
	    SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
	    if (dmrginp.hamiltonian() == BCS && dn == -3) {
          std::vector<SpinQuantum> sij = (-si)-sj;
          for (int ij = 0; ij < sij.size(); ++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij-sl;
            for (int ijl=0;ijl<sijl.size();++ijl){
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp DI(_i, -1), DJ(_j, -1), DL(_l, -1);
              TensorOp DDIJ = DI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
              TensorOp DDDIJL = DDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (DDDIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, DDDIJL, ladder[i]);
              double scale = calcCompfactor(DDDIJL, CK, CDD, v_cccd);
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }
          }
        } else if (dmrginp.hamiltonian() == BCS && dn == 1) {
          std::vector<SpinQuantum> sij = si+sj;
          for (int ij=0; ij<sij.size(); ++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij-sl;
            for (int ijl=0; ijl<sijl.size(); ++ijl) {
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
              TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
              TensorOp CCDIJL = CCIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (CCDIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, CCDIJL, ladder[i]);
              double scale = calcCompfactor(CCDIJL, CK, CDD, v_cccd);
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }            
          }
        } else if (dmrginp.hamiltonian() == BCS && dn == 3) {
          std::vector<SpinQuantum> sij = si+sj;
          for (int ij =0; ij<sij.size();++ij) {
            SpinQuantum symij = sij[ij];
            std::vector<SpinQuantum> sijl = symij+sl;
            for (int ijl=0; ijl<sijl.size(); ++ijl) {            
              if (sijl[ijl] != deltaQuantum[j]) continue;
              SpinQuantum symijl = sijl[ijl];
              TensorOp CI(_i, 1), CJ(_j, 1), CL(_l, 1);
              TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
              TensorOp CCCIJL = CCIJ.product(CL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
              if (CCCIJL.empty) continue;
              std::vector<double> MatElements = calcMatrixElements(c1, CCCIJL, ladder[i]);
              double scale = calcCompfactor(CCCIJL, CK, CDD, v_cccc);
              if (fabs(scale) > dmrginp.oneindex_screen_tol())
                element += MatElements[index]*scale/cleb;
            }
          }
        } else { //CDD
          std::vector<SpinQuantum> sij = si-sj;
	      for (int ij=0; ij<sij.size(); ++ij) {
	        SpinQuantum symij = sij[ij];
	        std::vector<SpinQuantum> sijl = symij-sl;
	        for (int ijl=0; ijl<sijl.size(); ijl++) {
	          SpinQuantum symijl = sijl[ijl];
	          if (symijl != deltaQuantum[j]) continue;
	          
	          TensorOp CI(_i, 1), DJ(_j, -1), DL(_l, -1);
	          
	          TensorOp CDIJ = CI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep());

	          TensorOp CDDIJL = CDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
	          if (CDDIJL.empty) continue;
	          std::vector<double> MatElements = calcMatrixElements(c1, CDDIJL, ladder[i]);
	          double scale = calcCompfactor(CDDIJL, CK, CDD, *(b->get_twoInt()), b->get_integralIndex());
	          if (dmrginp.spinAdapted()) scale*=-1; //terrible hack
	          if (fabs(scale) > dmrginp.oneindex_screen_tol()) 
	            element += MatElements[index]*scale/cleb;
	        }
	      }
        }
      }
      for (int ki =0; ki<b->get_sites().size(); ki++) {
	    int _i = b->get_sites()[ki];
        if (dmrginp.hamiltonian() == BCS && dn == 1) {
          TensorOp CI(_i, 1);
          std::vector<double> MatElements = calcMatrixElements(c1, CI, ladder[i]);
          double factor = calcCompfactor(CI, CK, C, *(b->get_twoInt()), b->get_integralIndex());
          if (fabs(factor) > dmrginp.oneindex_screen_tol())
            element += factor*MatElements[index]/cleb;
        } else {
	      TensorOp DI(_i, -1);
	      std::vector<double> MatElements = calcMatrixElements(c1, DI, ladder[i]);
	      double factor = calcCompfactor(CK, DI, C, *(b->get_twoInt()), b->get_integralIndex());
	      if (fabs(factor) > dmrginp.oneindex_screen_tol())
	        element += factor*MatElements[index]/cleb;
        }
      }
      break;
    }
    else
      continue;
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreDesDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreDesDesComp);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}


//******************HAM*****************

void SpinAdapted::Ham::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  SpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);
  //loopBlock = rightBlock; otherBlock = leftBlock; //**********


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif
  Ham *op_array =0;
  Ham *op_distributed=0;
  Ham *op_add=0;

  //initiateMultiThread(this, op_array, op_distributed, maxt);
  initiateMultiThread(this, op_array, op_distributed, MAX_THRD);


  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_rep(HAM, deltaQuantum);

  if (rightBlock->get_sites().size() == 0) 
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  else {
    //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
  }

  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);
    dmrginp.makeopsT -> stop();
    return;
  }

  op = rightBlock->get_op_rep(HAM, deltaQuantum);
  //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
  SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);


  // CCD_A*D_B + CCD_B*D_A + c.c. 
  op_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  op_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  f = boost::bind(&opxop::cxcddcomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  if (dmrginp.hamiltonian() != HUBBARD) {    
    op_add =  (otherBlock->get_op_array(CRE_DESCOMP).is_local() && loopBlock->get_op_array(CRE_DES).is_local())? op_array : op_distributed;
    f = boost::bind(&opxop::cdxcdcomp, otherBlock, _1, &b, op_add);
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);
    
    op_add =  (otherBlock->get_op_array(DES_DESCOMP).is_local() && loopBlock->get_op_array(CRE_CRE).is_local())? op_array : op_distributed;
    f = boost::bind(&opxop::ddxcccomp, otherBlock, _1, &b, op_add);
    for_all_multithread(loopBlock->get_op_array(CRE_CRE), f);
  }
  //accumulateMultiThread(this, op_array, op_distributed, maxt);
  accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);

  dmrginp.makeopsT -> stop();    
}


double SpinAdapted::Ham::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  const TwoElectronArray& v_2 = *(b->get_twoInt());
  double element = 0.0;
  bool finish = false;
  for (int i=0; i<ladder.size(); i++)
  {
    if (finish) break;
    bool isLallowed = Symmetry::spatial_cg(ladder[i].sym_is().getirrep(), 0, c1.sym_is().getirrep(), ladder[i].row(), 0, c1.row())!=0; // symmetry
    if ((c1.Sz != ladder[i].Sz || c1.S != ladder[i].S) || !isLallowed)
      continue;
    else
      finish = true; // only one element from the ladder has none zero matrix element with c1

    double matrixE = 0.0;
    for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
      for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++) {
	    Slater s1 = it1->first, s2 = it2->first; // slater determinants
	    double d1 = it1->second, d2 = it2->second; // weights

	    std::vector<int> cv, dv;
	    s1.connect(s2, cv, dv); // how to generate s1 from s2
	    if ((dv.size() == 2) && (cv.size() == 2)) {
	      int cI = cv[0]; 
	      int cJ = cv[1]; 
	      int dK = dv[0]; 
	      int dL = dv[1]; 
	      int parity = s1.trace(s2.d(dK).d(dL).c(cJ).c(cI));
	      double factor = parity*d1*d2*0.5;
	      matrixE += factor*(v_2(cI, cJ, dK, dL) - v_2(cJ, cI, dK, dL) - v_2(cI, cJ, dL, dK) + v_2(cJ, cI, dL, dK));
	    } else if ((cv.size() == 1) && (dv.size() == 1)) {
          // from v1
	      int cI = cv[0]; 
	      int dK = dv[0]; 
	      int parity = s1.trace(s2.d(dK).c(cI));
	      double factor = parity*d1*d2;
	      matrixE += factor*v_1[b->get_integralIndex()](cI, dK);
          // from v2
	      if(dmrginp.spinAdapted()) {
	        for (int kj=0; kj<b->get_sites().size(); kj++) {
	          int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
	          int num = 2*Symmetry::sizeofIrrep(SymmetryOfOrb(b->get_sites()[kj]).getirrep());
	          for (int J = jindex; J<num+jindex; J++) {	    
	    	    s1 = it1->first; s2 = it2->first;
	    	    parity = s1.trace(s2.d(dK).d(J).c(J).c(cI));
	    	    factor = parity*d1*d2*0.5;
	    	    matrixE += factor*(v_2(cI, J, dK, J) - v_2(J, cI, dK, J) - v_2(cI, J, J, dK) + v_2(J, cI, J, dK));
	          }
	        }
	      } else {
	        for (int kj=0; kj<b->get_sites().size(); kj++) {
	          int J = b->get_sites()[kj];
	          s1 = it1->first; s2 = it2->first;
	          parity = s1.trace(s2.d(dK).d(J).c(J).c(cI));
	          factor = parity*d1*d2*0.5;
	          matrixE += factor*(v_2(cI, J, dK, J) - v_2(J, cI, dK, J) - v_2(cI, J, J, dK) + v_2(J, cI, J, dK));	      
	        }
	      }
	    } else if ((cv.size() == 0) && (dv.size() == 0)) {
	      if(dmrginp.spinAdapted()) { // spin adapted
	        //T
	        for (int kj=0; kj<b->get_sites().size(); kj++) {
	          int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
	          int num = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kj]).getirrep());
	          for (int J = jindex; J<num+jindex; J++) {	    
	    	    s1 = it1->first; s2 = it2->first;
	    	    matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1[b->get_integralIndex()](J,J);
	          }
	        }
	        //V
	        for (int ki=0; ki<b->get_sites().size(); ki++)
	          for (int kk=0; kk<b->get_sites().size(); kk++) {
	    	    int Iindex = dmrginp.spatial_to_spin()[b->get_sites()[ki]];
	    	    int Inum = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[ki]).getirrep());
	    	    for (int I = Iindex; I<Inum+Iindex; I++) {	    
	    	      int Kindex = dmrginp.spatial_to_spin()[b->get_sites()[kk]];
	    	      int Knum = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kk]).getirrep());
	    	      for (int K = Kindex; K<Knum+Kindex; K++) {	    
	    	        double factor = 0.5*d1*d2; //if (ki == kk) factor = 1.0*d1*d2;
	    	        s1 = it1->first; s2 = it2->first;
	    	        matrixE += factor*s1.trace(s2.d(I).d(K).c(K).c(I))*(v_2(I, K, I, K) - v_2(K, I, I, K));
	    	      }
	    	    }
	          }
	      } else { // spin non-adapted
	        for (int kj=0; kj<b->get_sites().size(); kj++) {
	          int J = b->get_sites()[kj];
	          s1 = it1->first; s2 = it2->first;
	          matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1[b->get_integralIndex()](J,J);
	        }
	        //V
	        for (int ki=0; ki<b->get_sites().size(); ki++)
	          for (int kk=0; kk<b->get_sites().size(); kk++) {
	    	    int K = b->get_sites()[kk];
	    	    int I = b->get_sites()[ki];
	    	    double factor = 0.5*d1*d2; //if (ki == kk) factor = 1.0*d1*d2;
	    	    s1 = it1->first; s2 = it2->first;
	    	    matrixE += factor*s1.trace(s2.d(I).d(K).c(K).c(I))*(v_2(I, K, I, K) - v_2(K, I, I, K));
	          }
	      }
	    } else if (dmrginp.hamiltonian() == BCS && cv.size() == 4 && dv.size() == 0) {
          int cI = cv[0];
	      int cJ = cv[1];
	      int cK = cv[2];
	      int cL = cv[3];
	      int parity = s1.trace(s2.c(cL).c(cK).c(cJ).c(cI));
	      double factor = parity*d1*d2;
          matrixE += factor*v_cccc(cI, cJ, cK, cL);
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 0 && dv.size() == 4) {
          int dI = dv[0];
	      int dJ = dv[1];
	      int dK = dv[2];
	      int dL = dv[3];
	      int parity = s1.trace(s2.d(dL).d(dK).d(dJ).d(dI));
	      double factor = parity*d1*d2;
          matrixE += factor*v_cccc(dL, dK, dJ, dI);  
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 3 && dv.size() == 1) {
          int cI = cv[0];
	      int cJ = cv[1];
	      int cK = cv[2];
	      int dL = dv[0];
	      int parity = s1.trace(s2.d(dL).c(cK).c(cJ).c(cI));
	      double factor = parity*d1*d2;
          matrixE += factor*v_cccd(cI,cJ,cK,dL);
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 1 && dv.size() == 3) {
          int cI = cv[0];
	      int dJ = dv[0];
	      int dK = dv[1];
	      int dL = dv[2];
	      int parity = s1.trace(s2.d(dL).d(dK).d(dJ).c(cI));
	      double factor = parity*d1*d2;
          matrixE += factor*v_cccd(dL,dK,dJ,cI);
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 2 && dv.size() == 0) {
          // from v_cc pairing
          int cI = cv[0];
          int cJ = cv[1];
          int parity = s1.trace(s2.c(cJ).c(cI));
          double factor = parity*d1*d2;
          matrixE += factor * (v_cc(cI,cJ)-v_cc(cJ,cI));
          // from v_cccd
          if (dmrginp.spinAdapted()) {
            pout << "Oops... BCS+SpinAdaption not implemented yet!" << endl;
            abort();
          } else {
            for (int kl = 0; kl < b->get_sites().size(); ++kl) {
              int K = b->get_sites()[kl];
              s1 = it1->first; s2 = it2->first;
              parity = s1.trace(s2.d(K).c(K).c(cJ).c(cI));
              factor = parity*d1*d2;
              matrixE += factor*v_cccd(cI,cJ,K,K);
            }
          }
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 0 && dv.size() == 2) {
          // from v_cc pairing
          int dI = dv[0];
          int dJ = dv[1];
          int parity = s1.trace(s2.d(dI).d(dJ));
          double factor = parity*d1*d2;
          matrixE += factor * (v_cc(dI,dJ)-v_cc(dJ,dI));
          // from v_cccd
          if (dmrginp.spinAdapted()) {
            pout << "Oops... BCS+SpinAdaption not implemented yet!" << endl;
            abort();
          } else {
            for (int kl = 0; kl < b->get_sites().size(); ++kl) {
              int K = b->get_sites()[kl];
              s1 = it1->first; s2 = it2->first;
              parity = s1.trace(s2.d(dI).d(dJ).d(K).c(K));
              factor = parity*d1*d2;
              matrixE += factor*v_cccd(dI,dJ,K,K);;
            }
          }
        }
      }
    }
    element += 	matrixE;
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::Ham::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<Ham>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new Ham);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//******************Overlap*****************

void SpinAdapted::Overlap::build(const SpinBlock& b)
{
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  SpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_rep(OVERLAP, deltaQuantum);

  if (rightBlock->get_sites().size() == 0) 
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  else {
    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(OVERLAP, deltaQuantum);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op2, *op, &b, &(b.get_stateInfo()), *this, 1.0);

  }
  return;


}


double SpinAdapted::Overlap::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  bool finish = false;
  for (int i=0; i<ladder.size(); i++)
  {
    if (finish) break;
    bool isLallowed = Symmetry::spatial_cg(ladder[i].sym_is().getirrep(), 0, c1.sym_is().getirrep(), ladder[i].row(), 0, c1.row())!=0; // symmetry
    if ((c1.Sz != ladder[i].Sz || c1.S != ladder[i].S) || !isLallowed)
      continue;
    else
      finish = true; // only one element from the ladder has none zero matrix element with c1

    double matrixE = 0.0;
    for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) 
      for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++) {
	Slater s1 = it1->first, s2 = it2->first; // slater determinants
	double d1 = it1->second, d2 = it2->second; // weights
	
	std::vector<int> cv, dv;
	s1.connect(s2, cv, dv); // how to generate s1 from s2
	matrixE += d1*d2*s1.trace(s2);
      }
    element += 	matrixE;
  }
  return element;
}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::Overlap::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<Overlap>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new Overlap);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}


//******************CDD_sum*****************

void SpinAdapted::CDD_sum::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif
  CDD_sum *op_array =0;
  CDD_sum *op_distributed=0;
  CDD_sum *op_add=0;

  //initiateMultiThread(this, op_array, op_distributed, maxt);
  initiateMultiThread(this, op_array, op_distributed, MAX_THRD);


  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_rep(CDD_SUM, deltaQuantum);

  if (rightBlock->get_sites().size() == 0) {
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);
    dmrginp.makeopsT -> stop();
    return;
  }
  else {
    //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
  }

  op = rightBlock->get_op_rep(CDD_SUM, deltaQuantum);
  //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
  SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);


  // CD_A*D_B + DD_B*C_A + D_A*CD_B + C_A*DD_B 
  op_add = ( leftBlock->get_op_array(CDD_CRE_DESCOMP).is_local() && rightBlock->get_op_array(DES).is_local() )? op_array : op_distributed;
  Functor f = boost::bind(&opxop::cdd_dxcdcomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(DES), f);

  op_add = ( rightBlock->get_op_array(CDD_CRE_DESCOMP).is_local() && leftBlock->get_op_array(DES).is_local() )? op_array : op_distributed;
  f = boost::bind(&opxop::cdd_dxcdcomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(DES), f);  

  op_add = ( leftBlock->get_op_array(CDD_DES_DESCOMP).is_local() && rightBlock->get_op_array(CRE).is_local() )? op_array : op_distributed;
  f = boost::bind(&opxop::cdd_cxddcomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  op_add = ( rightBlock->get_op_array(CDD_DES_DESCOMP).is_local() && leftBlock->get_op_array(CRE).is_local() )? op_array : op_distributed;
  f = boost::bind(&opxop::cdd_cxddcomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  //accumulateMultiThread(this, op_array, op_distributed, maxt);
  accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);

  dmrginp.makeopsT -> stop();    
}


//double SpinAdapted::CDD_sum::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
//{
//  const TwoElectronArray& v_2 = *(b->get_twoInt());
//  double element = 0.0;
//  bool finish = false;
//  assert(b->nonactive_orb().size() == 1);
//  int orb= b->nonactive_orb(0);
//  for (int j = 0; j < deltaQuantum.size(); ++j) {  
//  for (int i=0; i<ladder.size(); i++)
//  {
//    if (finish) break;
//    int index =0; double cleb =0.0;
//    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
//      finish = true; // only one element from the ladder has none zero matrix element with c1
//
//      double matrixE = 0.0;
//      for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
//        for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++) {
//          Slater s1 = it1->first, s2 = it2->first; // slater determinants
//          double d1 = it1->second, d2 = it2->second; // weights
//
//          std::vector<int> cv, dv;
//          s1.connect(s2, cv, dv); // how to generate s1 from s2
//          if ((dv.size() == 2) && (cv.size() == 1)) {
//            int cI = cv[0]; 
//            int dK = dv[0]; 
//            int dL = dv[1]; 
//            int parity = s1.trace(s2.d(dK).d(dL).c(cI));
//            double factor = parity*d1*d2;
//            //TODO
//            //How to get the integral
//            //matrixE += factor*(Va_integral(orb, cI, dK, dL) - Va_integral(orb,cI, dL, dK));
//            if(dmrginp.spinAdapted())
//              matrixE += factor*(vpt_2[Va](2*orb, cI, dK, dL) - vpt_2[Va](2*orb,cI, dL, dK)+vpt_2[Va](2*orb+1, cI, dK, dL) - vpt_2[Va](2*orb+1,cI, dL, dK));
//            else
//              matrixE += factor*(vpt_2[Va](orb, cI, dK, dL) - vpt_2[Va](orb,cI, dL, dK));
//
//          }
//          else if ((cv.size() == 0) && (dv.size() == 1)) {
//            // from v1
//            int dK = dv[0]; 
//            int parity = s1.trace(s2.d(dK));
//            double factor = parity*d1*d2;
//            //TODO
//            if(dmrginp.spinAdapted())
//              matrixE += factor*(vpt_1(2*orb, dK)+vpt_1(2*orb+1,dK));
//            else
//              matrixE += factor*(vpt_1(orb, dK));
//
//            // from v2
//            if(dmrginp.spinAdapted()) {
//              //TODO
//              for (int kj=0; kj<b->get_sites().size(); kj++) {
//                int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
//                int num = 2*Symmetry::sizeofIrrep(SymmetryOfOrb(b->get_sites()[kj]).getirrep());
//                for (int J = jindex; J<num+jindex; J++) {	    
//                  s1 = it1->first; s2 = it2->first;
//                  parity = s1.trace(s2.d(dK).d(J).c(J));
//                  factor = parity*d1*d2;
//                  //matrixE += factor*(Va_integral(orb, J, dK, J) - Va_integral(orb, J, J, dK)) ;
//                  matrixE += factor*(vpt_2[Va](2*orb, J, dK, J) - vpt_2[Va](2*orb, J, J, dK)+vpt_2[Va](2*orb+1, J, dK, J) - vpt_2[Va](2*orb+1, J, J, dK)) ;
//
//                }
//              }
//            }
//            else {
//              //TODO
//              for (int kj=0; kj<b->get_sites().size(); kj++) {
//                int J = b->get_sites()[kj];
//                s1 = it1->first; s2 = it2->first;
//                parity = s1.trace(s2.d(dK).d(J).c(J));
//                factor = parity*d1*d2;
//                //matrixE += factor*(Va_integral(orb, J, dK, J) - Va_integral(orb, J, J, dK)) ;
//                matrixE += factor*(vpt_2[Va](orb, J, dK, J) - vpt_2[Va](orb, J, J, dK)) ;
//
//              }
//            }
//          }
//        }
//        element += 	matrixE;
//      }
//    }
//    }
//  }
//  return element;
//}
//

double SpinAdapted::CDD_sum::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  const TwoElectronArray& v_2 = *(b->get_twoInt());
  double element = 0.0;
  bool finish = false;
  assert(b->nonactive_orb().size() == 1);
  int orb= b->nonactive_orb(0);
  TensorOp CK(orb, 1);
  for (int j = 0; j < deltaQuantum.size(); ++j) {  
  for (int i=0; i<ladder.size(); i++)
  {
    int index =0; double cleb =0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
      for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++) 
      for (int kl =0; kl<b->get_sites().size(); kl++) {
	  int _i = b->get_sites()[ki];
	  int _j = b->get_sites()[kj];
	  int _l = b->get_sites()[kl];
	  SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
          std::vector<SpinQuantum> sij = si-sj;
	  for (int ij=0; ij<sij.size(); ++ij) {
	    SpinQuantum symij = sij[ij];
	    std::vector<SpinQuantum> sijl = symij-sl;
	    for (int ijl=0; ijl<sijl.size(); ijl++) {
	      SpinQuantum symijl = sijl[ijl];
	      if (symijl != deltaQuantum[j]) continue;
	      
	      TensorOp CI(_i, 1), DJ(_j, -1), DL(_l, -1);
	      
	      TensorOp CDIJ = CI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep());

	      TensorOp CDDIJL = CDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
	      if (CDDIJL.empty) continue;
	      std::vector<double> MatElements = calcMatrixElements(c1, CDDIJL, ladder[i]);
	      double scale = calcCompfactor(CK,CDDIJL, CDD, vpt_2[Va], b->get_integralIndex());
	          //if (dmrginp.spinAdapted()) scale*=-1; //terrible hack
	      if (fabs(scale) > dmrginp.oneindex_screen_tol()) 
	        element += MatElements[index]*scale/cleb;
	    }
	  }
      }
      for (int ki =0; ki<b->get_sites().size(); ki++) {
	int _i = b->get_sites()[ki];
	TensorOp DI(_i, -1);
	std::vector<double> MatElements = calcMatrixElements(c1, DI, ladder[i]);
	double factor = calcCompfactor(CK, DI, C, vpt_1, b->get_integralIndex());
	if (fabs(factor) > dmrginp.oneindex_screen_tol())
	  element += factor*MatElements[index]/cleb;
      }

      break;
    }
    else
      continue;
  }
  }
  return element;
}


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CDD_sum::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CDD_sum>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CDD_sum);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}



//******************CDD_DESDESCOMP*****************

void SpinAdapted::CDD_DesDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];

  TensorOp C(b.nonactive_orb(0),1), C2(i,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), false);

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CDD_DES_DESCOMP).has(i))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CDD_DES_DESCOMP, deltaQuantum, i);
    if (rightBlock->get_sites().size() == 0) {
      //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      return;
    }
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  }

  if (rightBlock->get_op_array(CDD_DES_DESCOMP).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CDD_DES_DESCOMP, deltaQuantum, i);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // explicitly build DD_comp

  if (leftBlock->get_op_array(CDD_DES_DESCOMP).has(i) && rightBlock->get_op_array(CDD_DES_DESCOMP).has(i))
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp DK(k,-1), DL(l,-1);
      TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
      if (!DD2.empty) {
        //double scaleV = calcCompfactor(CC1, DD2, DD, Va_integral, b.get_integralIndex());
        double scaleV = calcCompfactor(CC1, DD2, DD, vpt_2[Va], b.get_integralIndex());

        DK=TensorOp(k,-1); DL=TensorOp(l,-1);
        DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
        //double scaleV2 = calcCompfactor(CC1, DD2, DD, Va_integral, b.get_integralIndex());
        double scaleV2 = calcCompfactor(CC1, DD2, DD, vpt_2[Va], b.get_integralIndex());
        
        if (leftBlock->get_op_array(DES).has(k) && rightBlock->get_op_array(DES).has(l) && (fabs(scaleV)+fabs(scaleV2)) > dmrginp.twoindex_screen_tol()) {
	    boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	    
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    scaleV += parity*scaleV2;
	    
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	  }

      }

    }
  dmrginp.makeopsT -> stop();

}


double SpinAdapted::CDD_DesDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0];
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp C(b->nonactive_orb(0),1), C2(I,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), false);
 
  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {

      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];	
            
            TensorOp DK(k,-1), DL(l,-1);
            TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);

            if (!DD2.empty) {
              std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
              //TODO
              //A new kind of calcCompfactor
              //double scale = calcCompfactor(CC1, DD2, DD, index, Va_integral, b->get_integralIndex());
              double scale = calcCompfactor(CC1, DD2, DD, index, vpt_2[Va], b->get_integralIndex());
              element += MatElements[index]*scale/cleb;
            }
          }
        break;
      }
      else
        continue;

    }
  }
  return element;

}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CDD_DesDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CDD_DesDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CDD_DesDesComp);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//******************CDD_CREDESCOMP*****************


void SpinAdapted::CDD_CreDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];

  TensorOp C(b.nonactive_orb(0),1), D(i,-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CDD_CRE_DESCOMP).has(i)) { 
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CDD_CRE_DESCOMP, deltaQuantum, i);
    if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      return;
    }
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }

  }
  if (rightBlock->get_op_array(CDD_CRE_DESCOMP).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CDD_CRE_DESCOMP, deltaQuantum, i);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // build CDcomp explicitely
  if (leftBlock->get_op_array(CDD_CRE_DESCOMP).has(i) && rightBlock->get_op_array(CDD_CRE_DESCOMP).has(i))
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	    //double scaleV = calcCompfactor(CD1, CD2, CD,Va_integral, b.get_integralIndex());
	    double scaleV = calcCompfactor(CD1, CD2, CDD_CD,vpt_2[Va], b.get_integralIndex());
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(DES).has(l) && (fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	    }
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	//double scaleV = calcCompfactor(CD1, CD2, CD,Va_integral, b.get_integralIndex());
      	double scaleV = calcCompfactor(CD1, CD2, CDD_CD,vpt_2[Va], b.get_integralIndex());

      	if (leftBlock->get_op_array(DES).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
  	    }
      }

    }
  //FIXME
  //No need to add the term (I|a) in cdd_cdcomp
  //They are added in the term cdd_sum\times overlap
  //cdd_cdcomp is \sum_{k,l} (Ia|kl)C_kD_l + (I|a). (I|a) is summed twice.
//  if (leftBlock->get_op_array(CDD_CRE_DESCOMP).has(i) && rightBlock->get_op_array(CDD_CRE_DESCOMP).has(i)) { 
//    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
//    const boost::shared_ptr<SparseMatrix>& op1 = leftBlock->get_op_rep(OVERLAP,hq);
//    const boost::shared_ptr<SparseMatrix>& op2 = rightBlock->get_op_rep(OVERLAP,hq);
//    double scale = -vpt_1(b.nonactive_orb(0),i);
//    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scale);
//  }  
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CDD_CreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is(); // classify whether we calculate CC or CD

  TensorOp C(b->nonactive_orb(0),1), D(I,-1);

  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep());

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];
            
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
                //double factor = calcCompfactor(CD1, CD2, CD, Va_integral, b->get_integralIndex());
                double factor = calcCompfactor(CD1, CD2, CDD_CD, vpt_2[Va], b->get_integralIndex());
                element += MatElements[index]*factor/cleb;
              }
          }
//        double matrixE=0;
//        for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) 
//          for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++) {
//            Slater s1 = it1->first, s2 = it2->first; // slater determinants
//            double d1 = it1->second, d2 = it2->second; // weights
//            
//            std::vector<int> cv, dv;
//            s1.connect(s2, cv, dv); // how to generate s1 from s2
//            matrixE += d1*d2*s1.trace(s2);
//          }
//        element += matrixE*vpt_1(b->nonactive_orb(0),I)/cleb;

        break;
      }
      else
        continue;
    }
  }
  return element;
}


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CDD_CreDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CDD_CreDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CDD_CreDesComp);
      *rep = *this;

      rep->build(*block);
      return rep;
    }

}

//******************CCD_sum*****************

void SpinAdapted::CCD_sum::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif
  CCD_sum *op_array =0;
  CCD_sum *op_distributed=0;
  CCD_sum *op_add=0;

  //initiateMultiThread(this, op_array, op_distributed, maxt);
  initiateMultiThread(this, op_array, op_distributed, MAX_THRD);


  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_rep(CCD_SUM, deltaQuantum);

  if (rightBlock->get_sites().size() == 0) {
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);
    dmrginp.makeopsT -> stop();
    return;
  }
  else {
    //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
  }

  op = rightBlock->get_op_rep(CCD_SUM, deltaQuantum);
  //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
  SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);


  // CD_A*D_B + DD_B*C_A + D_A*CD_B + C_A*DD_B 
  op_add =  leftBlock->get_op_array(CCD_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  Functor f = boost::bind(&opxop::ccd_cxcdcomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  op_add =  rightBlock->get_op_array(CCD_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  f = boost::bind(&opxop::ccd_cxcdcomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  op_add =  leftBlock->get_op_array(CCD_CRE_CRECOMP).is_local() ? op_array : op_distributed;
  f = boost::bind(&opxop::ccd_dxcccomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(DES), f);

  op_add =  rightBlock->get_op_array(CCD_CRE_CRECOMP).is_local() ? op_array : op_distributed;
  f = boost::bind(&opxop::ccd_dxcccomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(DES), f);  

  //accumulateMultiThread(this, op_array, op_distributed, maxt);
  accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);

  dmrginp.makeopsT -> stop();    
}


double SpinAdapted::CCD_sum::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  const TwoElectronArray& v_2 = *(b->get_twoInt());
  double element = 0.0;
  bool finish = false;
  assert(b->nonactive_orb().size() == 1);
  int orb= b->nonactive_orb(0);
  TensorOp DK(orb, -1);
  for (int j = 0; j < deltaQuantum.size(); ++j) {  
  for (int i=0; i<ladder.size(); i++)
  {
    int index =0; double cleb =0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
      for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++) 
      for (int kl =0; kl<b->get_sites().size(); kl++) {
	      int _i = b->get_sites()[ki];
	      int _j = b->get_sites()[kj];
	      int _l = b->get_sites()[kl];
	      SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
        std::vector<SpinQuantum> sij = si+sj;
	      for (int ij=0; ij<sij.size(); ++ij) {
	        SpinQuantum symij = sij[ij];
	        std::vector<SpinQuantum> sijl = symij-sl;
	        for (int ijl=0; ijl<sijl.size(); ijl++) {
	          SpinQuantum symijl = sijl[ijl];
	          if (symijl != deltaQuantum[j]) continue;
	          
	          TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
	          
	          TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep());

	          TensorOp CCDIJL = CCIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
	          if (CCDIJL.empty) continue;
	          std::vector<double> MatElements = calcMatrixElements(c1, CCDIJL, ladder[i]);
	          double scale = calcCompfactor(CCDIJL, DK, CCD, vpt_2[Vi], b->get_integralIndex());
	          if (fabs(scale) > dmrginp.oneindex_screen_tol()) 
	            element += MatElements[index]*scale/cleb;
	        }
	      }
      }
      for (int ki =0; ki<b->get_sites().size(); ki++) {
	      int _i = b->get_sites()[ki];
	      TensorOp CI(_i, 1);
	      std::vector<double> MatElements = calcMatrixElements(c1, CI, ladder[i]);
	      double factor = calcCompfactor(CI, DK, C, vpt_1, b->get_integralIndex());
	      if (fabs(factor) > dmrginp.oneindex_screen_tol())
	        element += factor*MatElements[index]/cleb;
      }

      break;
    }
    else
      continue;
  }
  }
  return element;
}


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CCD_sum::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CCD_sum>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CCD_sum);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}



//******************CCD_CRECRECOMP*****************

void SpinAdapted::CCD_CreCreComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];

  TensorOp D2(b.nonactive_orb(0),-1), D(i,-1);
  TensorOp DD2 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), false);

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CCD_CRE_CRECOMP).has(i))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CCD_CRE_CRECOMP, deltaQuantum, i);
    if (rightBlock->get_sites().size() == 0) {
      //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      return;
    }
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  }

  if (rightBlock->get_op_array(CCD_CRE_CRECOMP).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CCD_CRE_CRECOMP, deltaQuantum, i);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // explicitly build DD_comp

  if (leftBlock->get_op_array(CCD_CRE_CRECOMP).has(i) && rightBlock->get_op_array(CCD_CRE_CRECOMP).has(i))
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), CL(l,1);
      TensorOp CC1 = CK.product(CL, spin, sym.getirrep(), k==l);
      if (!CC1.empty) {
        //double scaleV = calcCompfactor(CC1, DD2, DD, Va_integral, b.get_integralIndex());
        double scaleV = calcCompfactor(CC1, DD2, DD, vpt_2[Vi], b.get_integralIndex());

        CK=TensorOp(k,1); CL=TensorOp(l,1);
        CC1 = CL.product(CK, spin, sym.getirrep(), k==l);
        double scaleV2 = calcCompfactor(CC1, DD2, DD, vpt_2[Vi], b.get_integralIndex());
        
        if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV)+fabs(scaleV2)) > dmrginp.twoindex_screen_tol()) {
	    boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	    boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	    
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    scaleV += parity*scaleV2;
	    
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	  }
      }

    }
  dmrginp.makeopsT -> stop();

}


double SpinAdapted::CCD_CreCreComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0];
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp D2(b->nonactive_orb(0),-1), D(I,-1);
  TensorOp DD2 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), false);
 
  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {

      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];	
            
            TensorOp CK(k,1), CL(l,1);
            TensorOp CC1 = CK.product(CL, spin, sym.getirrep(), k==l);

            if (!CC1.empty) {
              std::vector<double> MatElements = calcMatrixElements(c1, CC1, ladder[i]);
              double scale = calcCompfactor(CC1, DD2, DD, index, vpt_2[Vi], b->get_integralIndex());
              element += MatElements[index]*scale/cleb;
            }
          }
        break;
      }
      else
        continue;

    }
  }
  return element;

}

boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CCD_CreCreComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CCD_CreCreComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CCD_CreCreComp);
      *rep = *this;
      rep->build(*block);
      return rep;
    }

}

//******************CCD_CREDESCOMP*****************


void SpinAdapted::CCD_CreDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];

  TensorOp C(i,1), D(b.nonactive_orb(0),-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CCD_CRE_DESCOMP).has(i)) { 
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CCD_CRE_DESCOMP, deltaQuantum, i);
    if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      return;
    }
    else {
      //const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<SparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }

  }
  if (rightBlock->get_op_array(CCD_CRE_DESCOMP).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CCD_CRE_DESCOMP, deltaQuantum, i);
    //const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<SparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  
  // build CDcomp explicitely
  if (leftBlock->get_op_array(CCD_CRE_DESCOMP).has(i) && rightBlock->get_op_array(CCD_CRE_DESCOMP).has(i))
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	    //double scaleV = calcCompfactor(CD1, CD2, CD,Va_integral, b.get_integralIndex());
	    double scaleV = calcCompfactor(CD2, CD1, CCD_CD,vpt_2[Vi], b.get_integralIndex());
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(DES).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	    }
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	//double scaleV = calcCompfactor(CD1, CD2, CD,Va_integral, b.get_integralIndex());
      	double scaleV = calcCompfactor(CD2, CD1, CCD_CD,vpt_2[Vi], b.get_integralIndex());

      	if (leftBlock->get_op_array(DES).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      boost::shared_ptr<SparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
  	    }
      }

    }
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CCD_CreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is(); // classify whether we calculate CC or CD

  TensorOp C(I,1), D(b->nonactive_orb(0),-1);

  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep());

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
    {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];
            
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
                //double factor = calcCompfactor(CD1, CD2, CD, Va_integral, b->get_integralIndex());
                double factor = calcCompfactor(CD2, CD1, CCD_CD, vpt_2[Vi], b->get_integralIndex());
                element += MatElements[index]*factor/cleb;
              }
          }

        break;
      }
      else
        continue;
    }
  }
  return element;
}


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CCD_CreDesComp::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CCD_CreDesComp>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CCD_CreDesComp);
      *rep = *this;

      rep->build(*block);
      return rep;
    }

}

