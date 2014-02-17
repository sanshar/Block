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

//using namespace SpinAdapted::operatorfunctions;

double SpinAdapted::vcccc_4idx_asymm(int i, int j, int k, int l) {
  return v_cccc(i,j,k,l)-v_cccc(i,j,l,k)+v_cccc(i,l,j,k)-v_cccc(l,i,j,k)
    -v_cccc(i,k,j,l)+v_cccc(i,k,l,j)-v_cccc(i,l,k,j)+v_cccc(l,i,k,j)
    +v_cccc(k,i,j,l)-v_cccc(k,i,l,j)+v_cccc(k,l,i,j)-v_cccc(l,k,i,j)
    -v_cccc(j,i,k,l)+v_cccc(j,i,l,k)-v_cccc(j,l,i,k)+v_cccc(l,j,i,k)
    +v_cccc(j,k,i,l)-v_cccc(j,k,l,i)+v_cccc(j,l,k,i)-v_cccc(l,j,k,i)
    -v_cccc(k,j,i,l)+v_cccc(k,j,l,i)-v_cccc(k,l,j,i)+v_cccc(l,k,j,i);
}

bool SpinAdapted::SparseMatrix::nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb)
{
  if (!dmrginp.spinAdapted()) {
    double clebsp = Symmetry::spatial_cg(ladder.sym_is().getirrep(), opsym.get_symm().getirrep(), c1.sym_is().getirrep(), ladder.row(), 0, c1.row());    
    if(c1.S.getirrep() == ladder.S.getirrep()+opsym.get_s().getirrep() &&
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


double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const TwoElectronArray& v_2)
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
	      factor += 0.5*(v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	    }
	    else if (comp == CCD) {
	      factor += 0.5*(v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]) - v_2(Ind1[1], Ind1[0], Ind2[0], Ind1[2]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	    }
	    else if (comp == C) {
          if (op1.dn() == 1) {
	        factor += 0.5*v_1(Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else { // D
            factor += 0.5*(v_cc(Ind2[0], Ind1[0]) 
                - v_cc(Ind1[0], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
          }
	    }
      }
  }
  
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2)
{
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, v_2);
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
	    else if (comp == C) {
	      factor += 0.5*v_1(Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
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
            factor += (1./6)*(vcccd(Ind1[0],Ind2[0],Ind2[1],Ind1[1]) - vcccd(Ind2[0],Ind1[0],Ind2[1],Ind1[1])
                +vcccd(Ind2[0],Ind2[1],Ind1[0],Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == DD) {
            factor += (1./6)*(vcccd(Ind1[0],Ind1[1],Ind2[0],Ind2[1]) - vcccd(Ind1[0],Ind2[0],Ind1[1],Ind2[1])
                +vcccd(Ind2[0],Ind1[0],Ind1[1],Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == CCD) { // two cases CDD and CCC
            if (op1.dn() == 3) { // CCC
              factor += (1./6) * vcccd(Ind1[0], Ind1[1], Ind1[2], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else { // CDD
              //factor += (1./6) * (vcccd(Ind2[0], Ind1[2], Ind1[1], Ind1[0]) - vcccd(Ind1[2], Ind2[0], Ind1[1], Ind1[0]) 
              //+ vcccd(Ind1[2], Ind1[1], Ind2[0], Ind1[0])) *iSz1.at(i1)*iSz2.at(i2)/cleb;
              factor += 0.5 * vcccd(Ind2[0], Ind1[2], Ind1[1], Ind1[0]) *iSz1.at(i1)*iSz2.at(i2)/cleb; // FIXME is this right?
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
  cout << "Sorry, SpinAdapted BCS calculation not implemented" << endl;
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
            factor += 0.25 * (vcccc(Ind1[0], Ind1[1], Ind2[0], Ind2[1]) 
                - vcccc(Ind1[0], Ind2[0], Ind1[1], Ind2[1]) + vcccc(Ind1[0], Ind2[0], Ind2[1], Ind1[1])
                + vcccc(Ind2[0], Ind1[0], Ind1[1], Ind2[1]) - vcccc(Ind2[0], Ind1[0], Ind2[1], Ind1[1])
                + vcccc(Ind2[0], Ind2[1], Ind1[0], Ind1[1]))* iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == CCD) {
            factor += 0.25*(vcccc(Ind2[0], Ind1[2], Ind1[1], Ind1[0])
                - vcccc(Ind1[2], Ind2[0], Ind1[1], Ind1[0]) + vcccc(Ind1[2], Ind1[1], Ind2[0], Ind1[0])
                - vcccc(Ind1[2], Ind1[1], Ind1[0], Ind2[0]))* iSz1.at(i1)*iSz2.at(i2)/cleb;
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
  cout << "Sorry, SpinAdapted BCS calculation not implemented" << endl;
  abort();
  return 0.;
}

//******************CRE*****************

void SpinAdapted::Cre::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());

  const int i = get_orbs()[0];
  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE).has(i))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE, deltaQuantum, i);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix>& op = rightBlock->get_op_rep(CRE, deltaQuantum, i);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
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


//******************CREDES*****************

void SpinAdapted::CreDes::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE_DES).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }  
  if (leftBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(j), j));
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  else if (rightBlock->get_op_array(CRE).has(i))
  {
    const boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
    Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(j), j));
    double parity = getCommuteParity(op1->get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[0]);
    // getCommuteParity doesn't depend on deltaQuantum.get_n()
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
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

//******************CRECRE*****************

void SpinAdapted::CreCre::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE).has(i, j))
  {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
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



//******************CREDESCOMP*****************


void SpinAdapted::CreDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  IrrepSpace sym = deltaQuantum[0].get_symm();  // sym and spin are the same for all dn part
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
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DESCOMP).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
  }  
  // build CDcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	      double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()));
	      if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	        Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	        SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
	      }
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()));

      	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	        Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	        double parity = getCommuteParity(op1->get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[0]);
	        SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
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
            
            TensorOp CK(k,1), DL(l,-1);
            TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
            if (!CD2.empty) {
              std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
              double factor = calcCompfactor(CD1, CD2, CD, *(b->get_twoInt()));
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
    }

}

//******************CREDESCOMP-No-Symm, used for BCS part **********


void SpinAdapted::CreDesComp_No_Symm::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  IrrepSpace sym = deltaQuantum[0].get_symm();  // sym and spin are the same for all dn part
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp C(i,1), D(j,-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DESCOMP_No_Symm).has(i, j))
  { 
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_DESCOMP_No_Symm, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DESCOMP_No_Symm).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_DESCOMP_No_Symm, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
  }  
  // build CDcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), CL(l,1);
      TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l); // k cannot equal to l
      if (!CC2.empty) {
        double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd);
        CL = TensorOp(l, 1);
        CK = TensorOp(k, 1);
        TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
        double scaleV2 = calcCompfactor(CD1, CC2_commute, CD, v_cccd);

        if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
          boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k); // FIXME is getSpinQuantum() proper?
          boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
          double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[1]);
          scaleV += parity * scaleV2;
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::CreDesComp_No_Symm::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
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
            
            TensorOp CK(k,1), CL(l,1);
            TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
            if (!CC2.empty) {
              std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i]);
              double factor = calcCompfactor(CD1, CC2, CD, v_cccd);
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


boost::shared_ptr<SpinAdapted::SparseMatrix> SpinAdapted::CreDesComp_No_Symm::getworkingrepresentation(const SpinBlock* block)
{
  //assert(this->get_initialised());
  if (this->get_built())
    {
      return boost::shared_ptr<CreDesComp_No_Symm>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
    }
  else
    {
      boost::shared_ptr<SparseMatrix> rep(new CreDesComp_No_Symm);
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
  allocate(b.get_stateInfo());
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
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(DES_DESCOMP).has(i, j))
  {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);
  }  
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp DK(k,-1), DL(l,-1);
      TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
      if (!DD2.empty) {
        double scaleV = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()));

        DK=TensorOp(k,-1); DL=TensorOp(l,-1);
        DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
        double scaleV2 = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()));
        
        if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	      Transposeview top1 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	      Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	
	      double parity = getCommuteParity(top1.get_deltaQuantum()[0], top2.get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      scaleV += parity*scaleV2;

	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	        SpinAdapted::operatorfunctions::TensorProduct(leftBlock, top1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
        }
      }

      if (dmrginp.hamiltonian() == BCS) {
        TensorOp CK(k, 1), CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CC1, CC2, DD, v_cccc);
          CK = TensorOp(k, 1);
          CL = TensorOp(k, 2);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CC1, CC2_commute, DD, v_cccc);

          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	        boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<SparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(0));
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
	        Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	        SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
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
	        Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	        double parity = getCommuteParity(op1->get_deltaQuantum(0), top2.get_deltaQuantum(0), get_deltaQuantum(1));
	        SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
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
  int spin = (-deltaQuantum[0].get_s()).getirrep();
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
                double scale = calcCompfactor(CC1, CC2, DD, index, v_cccc);
                element += MatElements[index]*scale/cleb;                
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == 0) {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
                double scale = calcCompfactor(CC1, CD2, DD, index, v_cccd);
                element += MatElements[index]*scale/cleb;
              }
            } else {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);

              if (!DD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
                double scale = calcCompfactor(CC1, DD2, DD, index, *(b->get_twoInt())); // FIXME what does the index do?
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




//******************CRECREDESCOMP*****************


void SpinAdapted::CreCreDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());

  const int k = get_orbs()[0];


  SpinBlock* leftBlock = b.get_leftBlock();
  SpinBlock* rightBlock = b.get_rightBlock();

  SpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);

  if (leftBlock->get_op_array(CRE_CRE_DESCOMP).has(k)) {      
    const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE_DESCOMP).has(k)) {
    const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
    SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
  }  

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY || dmrginp.hamiltonian() == BCS) {
    if (loopBlock->has(CRE_DESCOMP)) {

	  Functor f = boost::bind(&opxop::cxcdcomp, otherBlock, _1, &b, k, this, 1.0); 
	  for_all_singlethread(loopBlock->get_op_array(CRE), f);
	  
      f = boost::bind(&opxop::dxcccomp, otherBlock, _1, &b, k, this, 2.0); // factor of 2.0 because CCcomp_{ij} = -CCcomp_{ji} not neccesarily true for BCS case
	  for_all_singlethread(loopBlock->get_op_array(CRE), f);
          
	  f = boost::bind(&opxop::cxcdcomp, loopBlock, _1, &b, k, this, 1.0); 
	  for_all_singlethread(otherBlock->get_op_array(CRE), f);

	  f = boost::bind(&opxop::dxcccomp, loopBlock, _1, &b, k, this, 2.0);
	  for_all_singlethread(otherBlock->get_op_array(CRE), f);
	
    } else if (otherBlock->has(CRE_DESCOMP)) {
	  cout << "I should not be here"<<endl;exit(0);
    }
  }

  if (dmrginp.hamiltonian() == BCS) {
    Functor f = boost::bind(&opxop::cxcdcomp_no_symm, otherBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(loopBlock->get_op_array(CRE), f);
	f = boost::bind(&opxop::cxcdcomp_no_symm, loopBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(otherBlock->get_op_array(CRE), f);
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
  for (int i=0; i<ladder.size(); i++) {
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
              double scale = calcCompfactor(CCDIJL, D, CCD, *(b->get_twoInt()));
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
          double factor = calcCompfactor(DI, D, C, *(b->get_twoInt()));
          if (fabs(factor) > dmrginp.oneindex_screen_tol())
            element += factor*MatElements[index]/cleb;
        } else { // C
          TensorOp CI(_i, 1);
          std::vector<double> MatElements = calcMatrixElements(c1, CI, ladder[i]);
          double factor = calcCompfactor(CI, D, C, *(b->get_twoInt()));
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


//******************HAM*****************

void SpinAdapted::Ham::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());

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


  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_rep(HAM, deltaQuantum); // H_A
  
  SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);

  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);
    dmrginp.makeopsT -> stop();
    return;
  }

  op = rightBlock->get_op_rep(HAM, deltaQuantum);  // H_B
  SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);  

  // CCD_A*D_B + CCD_B*D_A + c.c. 
  op_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  op_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  f = boost::bind(&opxop::cxcddcomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY || dmrginp.hamiltonian() == BCS) {    
    op_add =  (otherBlock->get_op_array(CRE_DESCOMP).is_local() && loopBlock->get_op_array(CRE_DES).is_local())? op_array : op_distributed;
    f = boost::bind(&opxop::cdxcdcomp, otherBlock, _1, &b, op_add);
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);
    
    op_add =  (otherBlock->get_op_array(DES_DESCOMP).is_local() && loopBlock->get_op_array(CRE_CRE).is_local())? op_array : op_distributed;
    f = boost::bind(&opxop::ddxcccomp, otherBlock, _1, &b, op_add);
    for_all_multithread(loopBlock->get_op_array(CRE_CRE), f);
  }

  if (dmrginp.hamiltonian() == BCS) {
    op_add = (otherBlock->get_op_array(CRE_DESCOMP_No_Symm).is_local() && loopBlock->get_op_array(CRE_DES).is_local())? op_array : op_distributed;
    f = boost::bind(&opxop::cdxcdcomp_no_symm, otherBlock, _1, &b, op_add);
    for_all_multithread(loopBlock->get_op_array(CRE_DES), f);
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
	      matrixE += factor*v_1(cI, dK);
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
	    	    matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1(J,J);
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
	          matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1(J,J);
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
	      double factor = parity*d1*d2*0.25;
          matrixE += factor*vcccc_4idx_asymm(cI, cJ, cK, cL);
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 0 && dv.size() == 4) {
          int dI = dv[0];
	      int dJ = dv[1];
	      int dK = dv[2];
	      int dL = dv[3];
	      int parity = s1.trace(s2.d(dL).d(dK).d(dJ).d(dI));
	      double factor = parity*d1*d2*0.25;
          matrixE += factor*vcccc_4idx_asymm(dL, dK, dJ, dI);  
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 3 && dv.size() == 1) {
          int cI = cv[0];
	      int cJ = cv[1];
	      int cK = cv[2];
	      int dL = dv[0];
	      int parity = s1.trace(s2.d(dL).c(cK).c(cJ).c(cI));
	      double factor = parity*d1*d2*0.5;
          matrixE += factor*(v_cccd(cI,cJ,cK,dL)-v_cccd(cI,cK,cJ,dL)+v_cccd(cK,cI,cJ,dL)
              -v_cccd(cJ,cI,cK,dL)+v_cccd(cJ,cK,cI,dL)-v_cccd(cK,cJ,cI,dL));
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 1 && dv.size() == 3) {
          int cI = cv[0];
	      int dJ = dv[0];
	      int dK = dv[1];
	      int dL = dv[2];
	      int parity = s1.trace(s2.d(dL).d(dK).d(dJ).c(cI));
	      double factor = parity*d1*d2*0.5;
          matrixE += factor*(v_cccd(dL,dK,dJ,cI)-v_cccd(dL,dJ,dK,cI)+v_cccd(dJ,dL,dK,cI)
              -v_cccd(dK,dL,dJ,cI)+v_cccd(dK,dJ,dL,cI)-v_cccd(dJ,dK,dL,cI));
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 2 && dv.size() == 0) {
          // from v_cc pairing
          int cI = cv[0];
          int cJ = cv[1];
          int parity = s1.trace(s2.c(cJ).c(cI));
          double factor = parity*d1*d2*0.5;
          matrixE += factor * (v_cc(cI,cJ)-v_cc(cJ,cI));
          // from v_cccd
          if (dmrginp.spinAdapted()) {
            cout << "Oops... BCS+SpinAdaption not implemented yet!" << endl;
            abort();
          } else {
            for (int kl = 0; kl < b->get_sites().size(); ++kl) {
              int K = b->get_sites()[kl];
              s1 = it1->first; s2 = it2->first;
              parity = s1.trace(s2.d(K).c(K).c(cJ).c(cI));
              factor = parity*d1*d2*0.5;
              matrixE += factor*(v_cccd(cI,cJ,K,K)-v_cccd(cI,K,cJ,K)+v_cccd(K,cI,cJ,K)
                  -v_cccd(cJ,cI,K,K)+v_cccd(cJ,K,cI,K)-v_cccd(K,cJ,cI,K));
            }
          }
        } else if (dmrginp.hamiltonian() == BCS && cv.size() == 0 && dv.size() == 2) {
          // from v_cc pairing
          int dI = dv[0];
          int dJ = dv[1];
          int parity = s1.trace(s2.d(dI).c(dJ));
          double factor = parity*d1*d2*0.5;          
          matrixE += factor * (v_cc(dI,dJ)-v_cc(dJ,dI));
          // from v_cccd
          if (dmrginp.spinAdapted()) {
            cout << "Oops... BCS+SpinAdaption not implemented yet!" << endl;
            abort();
          } else {
            for (int kl = 0; kl < b->get_sites().size(); ++kl) {
              int K = b->get_sites()[kl];
              s1 = it1->first; s2 = it2->first;
              parity = s1.trace(s2.d(dI).d(dJ).d(K).c(K));
              factor = parity*d1*d2*0.5;
              matrixE += factor*(v_cccd(dI,dJ,K,K)-v_cccd(dI,K,dJ,K)+v_cccd(K,dI,dJ,K)
                  -v_cccd(dJ,dI,K,K)+v_cccd(dJ,K,dI,K)-v_cccd(K,dJ,dI,K));
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


