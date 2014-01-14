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

bool SpinAdapted::SparseMatrix::nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb)
{
  nonzeroindex = 0;
  cleb = 0.0;
  bool found = false;
  int spin = opsym.get_s();

  for (int Lz = 0; Lz< Symmetry::sizeofIrrep(opsym.get_symm().getirrep())&&!found; Lz++)
  for (int sz=spin; sz>-spin-1&&!found; sz-=2) 
  {
    cleb = Symmetry::spatial_cg(ladder.sym_is().getirrep(), opsym.get_symm().getirrep(), c1.sym_is().getirrep(), ladder.row(), Lz, c1.row()); 
    cleb *= cg(ladder.S, spin, c1.S, ladder.Sz, sz, c1.Sz);
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
  for (map<Slater, double>::iterator it2 = c2.det_rep.begin(); it2!= c2.det_rep.end(); it2++) 
  {
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
  for (int sz2=-op2.Spin; sz2< op2.Spin+1; sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;
    //int lz1 = op1.lz[0], lz2 = op2.lz[ilz2];
    std::vector<double>&  iSz2 = op2.Szops[ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2];
    
    //double cleb = cleb_(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);
    //pout << "cleb " <<  cleb << " op1.Spin " <<  op1.Spin << " m1 "<< op1.Spin << " op2.Spin " << op2.Spin << " m2 " << sz2 << endl;
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
	  factor += 0.5*v_1(Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
      }
  }
  
  return factor;
}

double SpinAdapted::SparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2)
{
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
  else
    abort();  

}


double SpinAdapted::Cre::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = dmrginp.spatial_to_spin()[get_orbs()[0]]; //convert spatial id to spin id because slaters need that
  int Iirrep = SymmetryOfSpatialOrb(get_orbs()[0]).getirrep();;
  IrrepSpace sym = deltaQuantum.get_symm();
  bool finish = false;
  bool write = false;
  int Sign = 1;

  TensorOp C(get_orbs()[0], 1);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, C, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
    
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
    double parity = getCommuteParity(op1->get_deltaQuantum(), top2.get_deltaQuantum(), get_deltaQuantum());
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
  IrrepSpace sym = deltaQuantum.get_symm();
  int irrep = deltaQuantum.get_symm().getirrep();
  int spin = deltaQuantum.get_s();

  TensorOp C(I, 1), D(J, -1);
  TensorOp CD = C.product(D, spin, irrep);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CD, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
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
    double parity = getCommuteParity(op1->get_deltaQuantum(), op2->get_deltaQuantum(), get_deltaQuantum());
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
  IrrepSpace sym = deltaQuantum.get_symm();
  int irrep = deltaQuantum.get_symm().getirrep();
  int spin = deltaQuantum.get_s();


  TensorOp C1(I,1), C2(J,1);
  TensorOp CC = C1.product(C2, spin, sym.getirrep(), I==J);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CC, ladder[i]) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
    
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
  IrrepSpace sym = deltaQuantum.get_symm();
  int spin = deltaQuantum.get_s();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp C(i,1), D(j,-1);
  TensorOp CD1 = C.product(D, spin, (-sym).getirrep());

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

  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
    {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()));
	
	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l)) {
	  boost::shared_ptr<SparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	  Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));

	  if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
	}
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()));

	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l)) {

	  boost::shared_ptr<SparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	  Transposeview top2 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	  
	  double parity = getCommuteParity(op1->get_deltaQuantum(), top2.get_deltaQuantum(), get_deltaQuantum());
	  if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
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
  IrrepSpace sym = deltaQuantum.get_symm();
  int spin = deltaQuantum.get_s();
  bool finish = false;

  TensorOp C(I,1), D(J,-1);

  TensorOp CD1 = C.product(D, spin, (-sym).getirrep());

  for (int i=0; i<ladder.size(); i++)
  {

    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      for (int kl =0; kl<b->get_sites().size(); kl++) 
      for (int kk =0; kk<b->get_sites().size(); kk++) {

	int k = b->get_sites()[kk];
	int l = b->get_sites()[kl];
	
	TensorOp CK(k,1), DL(l,-1);      
	TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
	if (CD2.empty) continue;
	std::vector<double> MatElements = calcMatrixElements(c1, CD2, ladder[i]);
	double factor = calcCompfactor(CD1, CD2, CD, *(b->get_twoInt()));
	element += MatElements[index]*factor/cleb;
	
      }
      break;
    }
    else
      continue;
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


//******************DESDESCOMP*****************

void SpinAdapted::DesDesComp::build(const SpinBlock& b)
{
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_stateInfo());
  int spin = deltaQuantum.get_s();
  IrrepSpace sym = deltaQuantum.get_symm();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  TensorOp C(i,1), C2(j,1);
  TensorOp CC1 = C.product(C2, spin, (-sym).getirrep(), i==j);

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
      if (DD2.empty) continue;
      double scaleV = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()));

      DK=TensorOp(k,-1); DL=TensorOp(l,-1);
      DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
      double scaleV2 = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()));

      if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l)) {
	Transposeview top1 = Transposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	Transposeview top2 = Transposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	
	double parity = getCommuteParity(top1.get_deltaQuantum(), top2.get_deltaQuantum(), get_deltaQuantum());      
	scaleV += parity*scaleV2;
	
	if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	  SpinAdapted::operatorfunctions::TensorProduct(leftBlock, top1, top2, &b, &(b.get_stateInfo()), *this, scaleV);
      }
    }
  dmrginp.makeopsT -> stop();

}


double SpinAdapted::DesDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum.get_symm();
  int spin = deltaQuantum.get_s();
  bool finish = false;

  TensorOp C(I,1), C2(J,1);

  TensorOp CC1 = C.product(C2, spin, (-sym).getirrep(), I==J);
 

  std::vector<double> values(2,0.0);
  for (int i=0; i<ladder.size(); i++)
  {

    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      for (int kl =0; kl<b->get_sites().size(); kl++) 
      for (int kk =0; kk<b->get_sites().size(); kk++) {

	int k = b->get_sites()[kk];
	int l = b->get_sites()[kl];	
	
	TensorOp DK(k,-1), DL(l,-1);
	TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
	if (DD2.empty) continue;

	std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i]);
	double scale = calcCompfactor(CC1, DD2, DD, index, *(b->get_twoInt()));

	element += MatElements[index]*scale/cleb;
	
      }
      break;
    }
    else
      continue;

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

  if (leftBlock->get_op_array(CRE_CRE_DESCOMP).has(k))
    {      
      const boost::shared_ptr<SparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE_DESCOMP).has(k))
    {
      const boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
      SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY){
    if (loopBlock->has(CRE_DESCOMP))
      {

	Functor f = boost::bind(&opxop::cxcdcomp, otherBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(loopBlock->get_op_array(CRE), f);
	
	f = boost::bind(&opxop::dxcccomp, otherBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(loopBlock->get_op_array(CRE), f);
        
	f = boost::bind(&opxop::cxcdcomp, loopBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(otherBlock->get_op_array(CRE), f);

	f = boost::bind(&opxop::dxcccomp, loopBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(otherBlock->get_op_array(CRE), f);
	
      }
    else if (otherBlock->has(CRE_DESCOMP))
      {
	cout << "I should not be here"<<endl;exit(0);
      }
  }
  dmrginp.makeopsT -> stop();


}


double SpinAdapted::CreCreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b)
{
  double element = 0.0;
  int K = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum.get_symm();
  int spin = deltaQuantum.get_s();
  bool finish = false;

  TensorOp D(K, -1);
  std::vector<double> values(4,0.0);
  for (int i=0; i<ladder.size(); i++)
  {

    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum, ladder[i], index, cleb)) {
      for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++) 
      for (int kl =0; kl<b->get_sites().size(); kl++) {

	int _i = b->get_sites()[ki];
	int _j = b->get_sites()[kj];
	int _l = b->get_sites()[kl];
	
	
	SpinQuantum si(1,1,SymmetryOfSpatialOrb(_i)), sj(1,1,SymmetryOfSpatialOrb(_j)), sl(-1,1,-SymmetryOfSpatialOrb(_l));

	std::vector<SpinQuantum> sij = si+sj;
	for (int ij=0; ij<sij.size(); ij++) {
	  SpinQuantum symij = sij[ij];
	  std::vector<SpinQuantum> sijk = symij+sl;
	  for (int ijk=0; ijk<sijk.size(); ijk++) {
	    SpinQuantum symijk = sijk[ijk];
	    if (symijk != deltaQuantum) continue;
	    
	    TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
	    
	    TensorOp CCIJ = CI.product(CJ, symij.get_s(), symij.get_symm().getirrep(), _i==_j);

	    TensorOp CCDIJL = CCIJ.product(DL, symijk.get_s(), symijk.get_symm().getirrep());
	    if (CCDIJL.empty) continue;
	    
	    std::vector<double> MatElements = calcMatrixElements(c1, CCDIJL, ladder[i]);
	    double scale = calcCompfactor(CCDIJL, D, CCD, *(b->get_twoInt()));



      	    element += MatElements[index]*scale/cleb;


	  }
	}
      }
      for (int ki =0; ki<b->get_sites().size(); ki++) {
	int _i = b->get_sites()[ki];
	TensorOp CI(_i, 1);
	std::vector<double> MatElements = calcMatrixElements(c1, CI, ladder[i]);
	double factor = calcCompfactor(CI, D, C, *(b->get_twoInt()));
	
	element += factor*MatElements[index]/cleb;
      }
      break;
    }
    else
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


  boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_rep(HAM, deltaQuantum);
  
  SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);

  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    accumulateMultiThread(this, op_array, op_distributed, MAX_THRD);
    dmrginp.makeopsT -> stop();    
    return;
  }

  op = rightBlock->get_op_rep(HAM, deltaQuantum);
  SpinAdapted::operatorfunctions::TensorTrace(rightBlock, *op, &b, &(b.get_stateInfo()), *this);  

  op_add =  leftBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  Functor f = boost::bind(&opxop::cxcddcomp, leftBlock, _1, &b, op_add); 
  for_all_multithread(rightBlock->get_op_array(CRE), f);

  op_add =  rightBlock->get_op_array(CRE_CRE_DESCOMP).is_local() ? op_array : op_distributed;
  f = boost::bind(&opxop::cxcddcomp, rightBlock, _1, &b, op_add); 
  for_all_multithread(leftBlock->get_op_array(CRE), f);  

  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {    
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
    bool isLallowed = Symmetry::spatial_cg(ladder[i].sym_is().getirrep(), 0, c1.sym_is().getirrep(), ladder[i].row(), 0, c1.row())!=0;
    if ((c1.Sz != ladder[i].Sz || c1.S != ladder[i].S) || !isLallowed)
      continue;
    else
      finish = true;

    double matrixE = 0.0;
    for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
      for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++)
      {
	
	Slater s1 = (*it1).first, s2 = (*it2).first;
	double d1 = (*it1).second, d2 = (*it2).second;

	std::vector<int> cv, dv;
	s1.connect(s2, cv, dv);
	if ((dv.size() == 2) && (cv.size() == 2))
	{
	  int cI = cv[0]; 
	  int cJ = cv[1]; 
	  int dK = dv[0]; 
	  int dL = dv[1]; 
	  int parity = s1.trace(s2.d(dK).d(dL).c(cJ).c(cI));
	  double factor = parity*d1*d2*0.5;
	  matrixE += factor*(v_2(cI, cJ, dK, dL) - v_2(cJ, cI, dK, dL) - v_2(cI, cJ, dL, dK) + v_2(cJ, cI, dL, dK));
	}
	if ((cv.size() == 1) && (dv.size() == 1))
	{
	  int cI = cv[0]; 
	  int dK = dv[0]; 
	  int parity = s1.trace(s2.d(dK).c(cI));
	  double factor = parity*d1*d2;
	  matrixE += factor*v_1(cI, dK);
	  
	  for (int kj=0; kj<b->get_sites().size(); kj++)
	  {
	    int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
	    int num = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kj]).getirrep());
	    for (int J = jindex; J<num+jindex; J++) {	    
	      s1 = (*it1).first; s2 = (*it2).first;
	      parity = s1.trace(s2.d(dK).d(J).c(J).c(cI));
	      factor = parity*d1*d2*0.5;
	      matrixE += factor*(v_2(cI, J, dK, J) - v_2(J, cI, dK, J) - v_2(cI, J, J, dK) + v_2(J, cI, J, dK));
	    }
	  }
	}
	
	if ((cv.size() == 0) && (dv.size() == 0))
	{
	  //T
	  for (int kj=0; kj<b->get_sites().size(); kj++)
	  {
	    int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
	    int num = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kj]).getirrep());
	    for (int J = jindex; J<num+jindex; J++) {	    
	      s1 = (*it1).first; s2 = (*it2).first;
	      matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1(J,J);
	    }
	  }
	  
	  //V
	  for (int ki=0; ki<b->get_sites().size(); ki++)
	    for (int kk=0; kk<b->get_sites().size(); kk++)
	    {
	      int Iindex = dmrginp.spatial_to_spin()[b->get_sites()[ki]];
	      int Inum = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[ki]).getirrep());
	      for (int I = Iindex; I<Inum+Iindex; I++) {	    
		
		int Kindex = dmrginp.spatial_to_spin()[b->get_sites()[kk]];
		int Knum = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kk]).getirrep());
		for (int K = Kindex; K<Knum+Kindex; K++) {	    
		  
		  double factor = 0.5*d1*d2; //if (ki == kk) factor = 1.0*d1*d2;
		  s1 = (*it1).first; s2 = (*it2).first;
		  matrixE += factor*s1.trace(s2.d(I).d(K).c(K).c(I))*(v_2(I, K, I, K) - v_2(K, I, I, K));

		}
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


