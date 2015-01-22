/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
n
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "csf.h"
#include "SpinQuantum.h"
#include "global.h"
#include "couplingCoeffs.h"
#include "boost/shared_ptr.hpp"
#include "tensor_operator.h"
#include "sweep_params.h"

SpinAdapted::Csf::Csf( const map<Slater, double>& p_dets, const int p_n, const SpinSpace p_S, const int p_Sz, const IrrepVector p_irrep) : det_rep(p_dets), n(p_n), S(p_S), Sz(p_Sz), irrep(p_irrep)
{
  map<Slater, double>::iterator it = det_rep.begin();
  for (; it!= det_rep.end(); it++)
    if ((*it).first.n != n || (*it).first.Sz != Sz )
      {
	pout << it->first<<endl;
	pout << (*it).first.n<<" "<<n<<endl;
	pout << (*it).first.Sz<<" "<<Sz<<endl;
	pout<<" all slaters in Csf do not have the same n, spin or symmetry"<<endl;
	pout << *this<<endl;
	throw 20;
      }
} 

bool SpinAdapted::Csf::operator< (const Csf& s) const
{return (SpinQuantum(n, S, sym_is()) < SpinQuantum(s.n, s.S, s.sym_is())); }


void SpinAdapted::Csf::outerProd(const Csf& csf, double factor, map<Slater, double>& output) const
{
  map<Slater, double>::const_iterator dets1 = det_rep.begin();
  map<Slater, double>::iterator itout;
  map<Slater, double>::const_iterator dets2 = csf.det_rep.begin();
  for (; dets1 != det_rep.end(); dets1++){
    dets2 = csf.det_rep.begin();
    for (; dets2 != csf.det_rep.end(); dets2++)
    {
      Slater s;
      (*dets1).first.outerProd((*dets2).first, s);
      itout = output.find(s);
      if (itout == output.end()) 
	output[s] = (*dets1).second*(*dets2).second*factor;
      else 
	output[s] += (*dets1).second*(*dets2).second*factor;
      
    }
  }
}

double SpinAdapted::csf_energy(const Csf& s, int integralIndex)
{
  map<Slater, double>::const_iterator it = s.det_rep.begin();
  double energy = 0.0;
  for (; it != s.det_rep.end();it++)
    energy += it->second* det_energy(it->first, integralIndex);

  return energy;
}
  
void SpinAdapted::Csf::applySplus(Csf& output)
{
  map<Slater, double>::iterator itout, it = det_rep.begin();
  map<Slater, double> detsout;
  for ( ; it!= det_rep.end(); it++)
    for (int i=0; i<Slater().size(); i+=2) {
      Slater s = (*it).first;
      s.d(i+1).c(i);
      if ( !s.isempty()) {
	itout = detsout.find(s);
	if (itout == detsout.end())
	  detsout[s] = (*it).second;
	else
	  detsout[s] += (*it).second;
      }
    }
  output.set_det_rep(detsout, S, output.irrep);
  output.set_S(S);
  output.set_Sz(Sz+2);
  output.set_n(n);
  output.normalize();
}

void SpinAdapted::Csf::applySminus(Csf& output)
{
  map<Slater, double>::iterator itout, it = det_rep.begin();
  map<Slater, double> detsout;
  for ( ; it!= det_rep.end(); it++)
    for (int i=0; i<Slater().size(); i+=2) {
      Slater s = (*it).first;
      s.d(i).c(i+1);
      if ( !s.isempty()) {
	itout = detsout.find(s);
	if (itout == detsout.end())
	  detsout[s] = (*it).second;
	else
	  detsout[s] += (*it).second;
      }
    }
  output.set_det_rep(detsout, S, irrep);
  output.set_S(S);
  output.set_Sz(Sz-2);
  output.set_n(n);
  output.normalize();
}

void SpinAdapted::Csf::applyRowminus(Csf& output, int irrep)
{

  map<Slater, double>::iterator itout, it = det_rep.begin();
  map<Slater, double> detsout;
  for ( ; it!= det_rep.end(); it++)
    for (int i=0; i<dmrginp.last_site(); i++) {
      if (SymmetryOfSpatialOrb(i).getirrep() != irrep)
	continue;
      int I = dmrginp.spatial_to_spin()[i];
      for (int j=0; j<2; j++) {
	Slater s = (*it).first;
	int sign = s.getSign();
	s.d(I+2+j).c(I+j);
	s.setSign(sign); //This is a hack.
	if ( !s.isempty()) {
	  itout = detsout.find(s);
	  if (itout == detsout.end())
	    detsout[s] = (*it).second;
	  else
	    detsout[s] += (*it).second;
	}
      }
    }
  output.set_det_rep(detsout, S, IrrepVector(this->irrep.getirrep(), 0));
}

vector<SpinAdapted::Csf> SpinAdapted::Csf::spinLadder(int k)
{

  vector< Csf> ladder;
  ladder.push_back(*this);

  for (int j=0; j<2; j++) {
    Csf s = *this;
    if (j==0)
      s = *this;
    else if (j== 1 && Symmetry::sizeofIrrep(irrep.getirrep()) != 1) {
      Csf sLminus;
      vector< vector<int> > partitions = CSFUTIL::generate_partitions( irrep.getirrep() );
      if (partitions.size() == 0) {pout << "big trouble"<<endl;exit(0);}
      for (int vec=0; vec<partitions.size(); vec++) {
	Csf stemp;
	s = *this;
	for (int l=0; l<partitions[vec].size(); l++) {
	  s.applyRowminus(stemp, partitions[vec][l]);
	  s = stemp;
	}
	if (vec == 0)
	  sLminus = s;
	else {
	  map<Slater, double> tempmap = sLminus.det_rep;
	  for (map<Slater, double>::iterator it2 = s.det_rep.begin(); it2!=s.det_rep.end(); it2++)
	    {
	      map<Slater, double>::iterator itout = tempmap.find(it2->first);
	      if (itout == tempmap.end())
		tempmap[it2->first] = it2->second;
	      else
		tempmap[it2->first] += it2->second;
	    }
	  sLminus.set_det_rep(tempmap, S, IrrepVector(this->irrep.getirrep(), 0));
	}
      }
      sLminus.normalize();
      ladder.push_back(sLminus);
    }
    else
      continue;
    for (int i = this->S.getirrep()-2; i>= this->S.getirrep() - 2*k; i-=2) {
      Csf sm;
      s.applySminus(sm);
      if (sm.size() != 0)
	ladder.push_back( sm);
      s = sm;
    }
  }
  return ladder;
}

void SpinAdapted::CSFUTIL::TensorProduct(Csf& lhs, vector<Csf>& lhs_csfs, Csf& rhs, vector<Csf>& rhs_csfs, vector< Csf >& output, vector< vector<Csf> >& outputladder)
{
  //vector< Csf> lhs_csfs = lhs.spinLadder(lhs.S);
  //vector< Csf> rhs_csfs = rhs.spinLadder(rhs.S);

  int lirrep = lhs.sym_is().getirrep();
  int rirrep = rhs.sym_is().getirrep(); 
  vector<IrrepSpace> symvec = lhs.sym_is()+rhs.sym_is();
  vector<SpinSpace> spinvec = lhs.S_is()+rhs.S_is();

  //now generate all the output vectors |J,J> with J = |j1-j2|,...,|j1+j2|
  //the angular momentum is m = m1+m2, |m1-m2|
  for (int l = 0; l<symvec.size(); l++)
  for (int j = 0; j<spinvec.size(); j++)
  {
    vector<Csf> thisladder;
    int irrep = symvec[l].getirrep();
    int J = spinvec[j].getirrep();
    for (int row=0; row<Symmetry::sizeofIrrep(irrep); row++)
    for (int M=-J; M<=J; M+=2)
    {
      //generating |J1 J2 J M>x|lirrep rirrep irrep row> vector

      map<Slater, double> dets;
      for (int j1m1 = 0; j1m1<lhs_csfs.size(); j1m1++) //loops over rows of J1 lirrep
      for (int j2m2 = 0; j2m2<rhs_csfs.size(); j2m2++) //loops over rows of J2 rirrep
      {
	int J1 = lhs.S.getirrep(), M1 = lhs_csfs[j1m1].Sz, J2 = rhs.S.getirrep(), M2 = rhs_csfs[j2m2].Sz;
	int rowl = lhs_csfs[j1m1].irrep.getrow(), rowr = rhs_csfs[j2m2].irrep.getrow();

	double clebsg = cg(J1, J2, J, M1, M2, M);
	double clebspatial = Symmetry::spatial_cg(lirrep, rirrep, irrep, rowl, rowr, row);
	if( abs(clebsg) < 1.0e-12 || abs(clebspatial) < 1.00e-14)
	  continue;
	lhs_csfs[j1m1].outerProd(rhs_csfs[j2m2], clebsg*clebspatial, dets); 
      }
      
      if (dets.size() == 0)
	continue;
      Csf c(dets, rhs.n+lhs.n, spinvec[j], M, IrrepVector(symvec[l].getirrep(), row)); c.normalize();
      if (M == J && row == 0)
	output.push_back(c);
      thisladder.push_back(c);
    }
    outputladder.push_back(thisladder);
  }
}

void SpinAdapted::CSFUTIL::TensorProduct(Csf& lhs, Csf& rhs, vector< Csf >& output)
{
  //vector< Csf> lhs_csfs = lhs.spinLadder(lhs.S);
  //vector< Csf> rhs_csfs = rhs.spinLadder(rhs.S);

  int lirrep = lhs.sym_is().getirrep();
  int rirrep = rhs.sym_is().getirrep(); 
  vector<IrrepSpace> symvec = lhs.sym_is()+rhs.sym_is();
  vector<SpinSpace> spinvec = lhs.S_is()+rhs.S_is();

  //now generate all the output vectors |J,J> with J = |j1-j2|,...,|j1+j2|
  //the angular momentum is m = m1+m2, |m1-m2|
  for (int l = 0; l<symvec.size(); l++)
  for (int j = 0; j<spinvec.size(); j++)
  {
    int irrep = symvec[l].getirrep();
    int J = spinvec[j].getirrep();
    for (int row=0; row<Symmetry::sizeofIrrep(irrep); row++)
    {
      map<Slater, double> dets;
      lhs.outerProd(rhs, 1.0, dets); 

      Csf c(dets, rhs.n+lhs.n, spinvec[j], J, IrrepVector(symvec[l].getirrep(), row)); c.normalize();

      output.push_back(c);
    }
  }
}


Csf SpinAdapted::CSFUTIL::applyTensorOp(const TensorOp& newop, int spinL)
{
  //for (int k=0; k<newop.Szops[newop.Szops.size()-1-newop.Spin].size(); k++) {
  map<Slater, double> m;
  SpinQuantum sq(newop.optypes.size(), SpinSpace(newop.Spin), IrrepSpace(newop.irrep)); 
  for (int k=0; k<newop.Szops[spinL].size(); k++) {
    if (fabs(newop.Szops[spinL][k]) > 1e-10) {
      std::vector<bool> occ_rep(Slater().size(), 0);
      Slater s(occ_rep,1);
      for (int k2=newop.opindices[k].size()-1; k2>=0; k2--) {
	s.c(newop.opindices[k][k2]);
      }
      if (s.isempty()) {
	continue;
      }
      map<Slater, double>::iterator itout = m.find(s);

      if (itout == m.end()) {
	int sign = s.alpha.getSign();
	s.setSign(1);
	m[s] = sign*newop.Szops[spinL][k];
      }	
      else
	m[s] += s.alpha.getSign()*newop.Szops[spinL][k];
    }
  }
  int irreprow = spinL/(newop.Spin+1); 
  int sz = -2*(spinL%(newop.Spin+1)) + newop.Spin;

  

  Csf csf = Csf(m,sq.particleNumber, sq.totalSpin, sz, IrrepVector(sq.get_symm().getirrep(),irreprow)) ; 
  if (!csf.isempty() && fabs(csf.norm()) > 1e-14)
    csf.normalize();
  
  return csf;
}

template<class T> class sorter {
  const std::vector<T> &values;
public:
  sorter(const std::vector<T> &v) : values(v) {}
  bool operator()(int a, int b) { return values[a] < values[b]; }
};

std::vector<SpinAdapted::Csf > SpinAdapted::CSFUTIL::spinfockstrings(const std::vector<int>& orbs, std::vector< std::vector<Csf> >& ladders)
{

  std::vector<Csf > singleSiteCsf;
  std::vector<int> numcsfs;
  std::vector< std::vector<Csf> > singleSiteLadder;
  int numCsfSoFar = 0;
  
  for (int i=0; i<orbs.size(); i++) {
    std::vector< Csf> thisSiteCsf;
    int I = orbs[i];  
    std::vector<TensorOp>  tensorops(1, TensorOp(I, 1));
    IrrepSpace Irrep = SymmetryOfSpatialOrb(orbs[i]); 
    SpinQuantum sQ(1, SpinSpace(1), Irrep);
    
    int irrepsize = Symmetry::sizeofIrrep(Irrep.getirrep());
    
    std::vector<Csf> ladderentry;
    std::map<Csf, std::vector<Csf> > laddermap;
    
    std::vector<bool> occ_rep1(Slater().size(),0), occ_rep2(Slater().size(),0);
    occ_rep2[dmrginp.spatial_to_spin()[I]+2*irrepsize-2] = 1;
    Slater s1(occ_rep1, 1), s2(occ_rep2, 1); map<Slater, double > m1, m2;
    m1[s1]= 1.0; m2[s2] = 1.0;

    if (find(dmrginp.get_openorbs().begin(), dmrginp.get_openorbs().end(), I) != dmrginp.get_openorbs().end() ) {
      thisSiteCsf.push_back( Csf(m1, 0, SpinSpace(0), 0, IrrepVector(0,0))); //0,0,0
      ladderentry.push_back(Csf(m1, 0, SpinSpace(0), 0, IrrepVector(0,0))); singleSiteLadder.push_back(ladderentry);
      numcsfs.push_back(thisSiteCsf.size());

      for (int i=0; i<thisSiteCsf.size(); i++)
	singleSiteCsf.push_back( thisSiteCsf[i]);
      continue;
    }
    else if (find(dmrginp.get_closedorbs().begin(), dmrginp.get_closedorbs().end(), I) != dmrginp.get_closedorbs().end()) {
      std::vector<bool> occ_rep(Slater().size(),0);
      occ_rep[dmrginp.spatial_to_spin()[I]+2*irrepsize-2] = 1;
      occ_rep[dmrginp.spatial_to_spin()[I]+2*irrepsize-1] = 1;
      Slater s(occ_rep, 1); map<Slater, double > m;
      m[s]= 1.0;
      thisSiteCsf.push_back( Csf(m, 2, SpinSpace(0), 0, IrrepVector(0,0))); //2,0,0
      ladderentry.push_back(Csf(m, 2, SpinSpace(0), 0, IrrepVector(0,0))); singleSiteLadder.push_back(ladderentry);
      numcsfs.push_back(thisSiteCsf.size());

      for (int i=0; i<thisSiteCsf.size(); i++)
	singleSiteCsf.push_back( thisSiteCsf[i]);
      continue;
    }
    else if(dmrginp.hamiltonian() == HEISENBERG) {
      thisSiteCsf.push_back( Csf(m2, 1, SpinSpace(1), 1, IrrepVector(Irrep.getirrep(), irrepsize-1))); //1,1,L    
      for (int i=tensorops[0].Szops.size(); i> 0; i--)
	ladderentry.push_back(applyTensorOp(tensorops[0], i-1));
      singleSiteLadder.push_back(ladderentry);

      numcsfs.push_back(thisSiteCsf.size());    
      for (int i=0; i<thisSiteCsf.size(); i++)
	singleSiteCsf.push_back( thisSiteCsf[i]);
      continue;
    }

    thisSiteCsf.push_back( Csf(m1, 0, SpinSpace(0), 0, IrrepVector(0,0))); //0,0,0
    thisSiteCsf.push_back( Csf(m2, 1, SpinSpace(1), 1, IrrepVector(Irrep.getirrep(), irrepsize-1))); //1,1,L    
    ladderentry.push_back(Csf(m1, 0, SpinSpace(0), 0, IrrepVector(0,0))); singleSiteLadder.push_back(ladderentry);
    ladderentry.clear();
    
    for (int i=tensorops[0].Szops.size(); i> 0; i--)
      ladderentry.push_back(applyTensorOp(tensorops[0], i-1));
    singleSiteLadder.push_back(ladderentry);
    //laddermap[singleSiteCsf.back()] = ladderentry;
    
    for (int nele = 2; nele < 2*irrepsize+1; nele++) {
      
      std::vector<SpinQuantum> quanta;
      std::vector<TensorOp> newtensorops;
      for (int i=0; i<tensorops.size(); i++) {
	quanta = SpinQuantum(nele-1, SpinSpace(tensorops[i].Spin), IrrepSpace(tensorops[i].irrep)) + sQ;
	
	for (int j=0; j<quanta.size(); j++) {
	  TensorOp newop = TensorOp(I,1).product(tensorops[i], quanta[j].totalSpin.getirrep(), quanta[j].get_symm().getirrep());
	  
	  Csf csf = applyTensorOp(newop, newop.Szops.size()-1-newop.Spin);
	  
	  if (!csf.isempty() && csf.norm() >1e-10) {
	    csf.normalize();

	    if (find(thisSiteCsf.begin(), thisSiteCsf.end(), csf) == thisSiteCsf.end()) {
	      thisSiteCsf.push_back( csf);
	      newtensorops.push_back(newop);
	      
	      std::vector<Csf> ladderentry;
	      for (int k=newop.Szops.size(); k>0 ; k--) 
		ladderentry.push_back(applyTensorOp(newop, k-1));
	      
	      singleSiteLadder.push_back(ladderentry);
	      //laddermap[csf] = ladderentry;
	    }
	  }
	}
      }
      
      tensorops = newtensorops;
    }
    
    numcsfs.push_back(thisSiteCsf.size());

    
    for (int i=0; i<thisSiteCsf.size(); i++)
      singleSiteCsf.push_back( thisSiteCsf[i]);
    
  }
  
  std::vector<Csf> prevoutput, output;
  std::vector< std::vector<Csf> > prevladder, ladder;

  for (int i=0; i<numcsfs[0]; i++) {
    output.push_back(singleSiteCsf[i]);
    ladder.push_back(singleSiteLadder[i]);
  }
  int csfindex = numcsfs[0];

  for (int i2=1; i2<orbs.size(); i2++) {


    for (int i=0; i<output.size(); i++) {
      prevoutput.push_back(output[i]);
      prevladder.push_back(ladder[i]);
    }
    output.clear();
    ladder.clear();

    for (int j=0; j<prevoutput.size(); j++) {
      for (int k=0; k<numcsfs[i2]; k++) {
	CSFUTIL::TensorProduct(prevoutput[j], prevladder[j], singleSiteCsf[csfindex+k], singleSiteLadder[csfindex+k], output, ladder);
      }
    }    
    

    prevoutput.clear();
    prevladder.clear();
    csfindex += numcsfs[i2];
  }

  std::vector<int> sortvec(output.size());
  for (int i=0; i<output.size(); i++)
    sortvec[i] = i;

  std::sort(sortvec.begin(), sortvec.end(), sorter<Csf>(output));
  std::sort(output.begin(), output.end());

  for (int i=0; i<output.size(); i++)
    ladders.push_back(ladder[sortvec[i]]);

  return output;
   
}


//this version of the code is used when DMRG is no spin adapted
std::vector<SpinAdapted::Csf > SpinAdapted::CSFUTIL::spinfockstrings(const std::vector<int>& orbs)
{

  std::vector<Csf > singleSiteCsf;
  std::vector<int> numcsfs;
  std::vector< std::vector<Csf> > singleSiteLadder;
  int numCsfSoFar = 0;
  
  for (int i=0; i<orbs.size(); i++) {
    std::vector< Csf> thisSiteCsf;
    int I = orbs[i];  

    IrrepSpace Irrep = SymmetryOfSpatialOrb(orbs[i]); 
    SpinQuantum sQ(1, SpinSpace(1), Irrep);
    
    int irrepsize = Symmetry::sizeofIrrep(Irrep.getirrep());
    
    std::vector<Csf> ladderentry;
    std::map<Csf, std::vector<Csf> > laddermap;
    
    std::vector<bool> occ_rep1(Slater().size(),0), occ_rep2(Slater().size(),0), occ_rep3(Slater().size(), 0), occ_rep4(Slater().size(), 0);
    occ_rep2[dmrginp.spatial_to_spin()[I]] = 1;
    occ_rep3[dmrginp.spatial_to_spin()[I]+1] = 1;
    occ_rep4[dmrginp.spatial_to_spin()[I]] = 1;occ_rep4[dmrginp.spatial_to_spin()[I]+1] = 1;
    
    Slater s1(occ_rep1, 1), s2(occ_rep2, 1), s3(occ_rep3, 1), s4(occ_rep4, 1); map<Slater, double > m1, m2, m3, m4;
    m1[s1]= 1.0; m2[s2] = 1.0;  m3[s3]= 1.0; m4[s4] = 1.0;
    thisSiteCsf.push_back( Csf(m1, 0, SpinSpace(0), 0, IrrepVector(0,0))); //0,0,0
    thisSiteCsf.push_back( Csf(m2, 1, SpinSpace(1), 1, IrrepVector(Irrep.getirrep(), irrepsize-1))); //1,1,L
    thisSiteCsf.push_back( Csf(m3, 1, SpinSpace(-1), -1, IrrepVector(Irrep.getirrep(), irrepsize-1))); //1,1,L
    thisSiteCsf.push_back( Csf(m4, 2, SpinSpace(0), 0, IrrepVector(0, 0))); //2,0,0

    
    numcsfs.push_back(thisSiteCsf.size());

    
    for (int i=0; i<thisSiteCsf.size(); i++)
      singleSiteCsf.push_back( thisSiteCsf[i]);
    
  }
  
  std::vector<Csf> prevoutput, output;

  for (int i=0; i<numcsfs[0]; i++) 
    output.push_back(singleSiteCsf[i]);

  int csfindex = numcsfs[0];

  for (int i2=1; i2<orbs.size(); i2++) {


    for (int i=0; i<output.size(); i++) {
      prevoutput.push_back(output[i]);
    }
    output.clear();

    for (int j=0; j<prevoutput.size(); j++) {
      for (int k=0; k<numcsfs[i2]; k++) {
	CSFUTIL::TensorProduct(prevoutput[j], singleSiteCsf[csfindex+k], output);
      }
    }    
    

    prevoutput.clear();
    csfindex += numcsfs[i2];
  }

  std::vector<int> sortvec(output.size());
  for (int i=0; i<output.size(); i++)
    sortvec[i] = i;

  std::sort(sortvec.begin(), sortvec.end(), sorter<Csf>(output));
  std::sort(output.begin(), output.end());

  return output;
   
}




std::vector< SpinAdapted::Csf > SpinAdapted::Csf::distribute (const int n, const int sp, const IrrepVector &sym, const int left, const int right, const int edge, int integralIndex)
{
  std::vector< Csf > s;

  int na = 0;
  int nb = 0;
  na = (n + sp) / 2;
  nb = (n - sp) / 2;

  if(dmrginp.hamiltonian() == HEISENBERG && nb != 0) return s;
  // let's count how many left and right orbitals there actually are                                                                      
  int actualOrbs = dmrginp.spatial_to_spin()[right] - dmrginp.spatial_to_spin()[left];
  int actualNA, actualNB;
  actualNA = actualOrbs/2;
  actualNB = actualOrbs/2;

  // cannot form this combination 

  if (na > actualNA || nb > actualNB || na < 0 || nb < 0)
    {
      return s;
    }

  // now make orbital occupancies                                                                                                         
  std::vector<int> abSeed (right-left);
  for (int i = 0; i < max(na, nb); ++i)
    abSeed [i] = -1;

  // okay now make some permutations (at most GUESS_PERMUTATIONS)  
  std::vector< std::vector<int> > abList; abList.push_back (abSeed);

  int numberOfGuesses = dmrginp.guess_permutations();
  bool doAgain = true;
  int numtries = 0;
  while(doAgain && numtries <500)
  {
    while (next_permutation (abSeed.begin (), abSeed.end ()) && numberOfGuesses--) 
      abList.push_back (abSeed);

    // make list of available alpha and beta orbitals in energy ordering.
    std::vector<int> orbitalList;
    
    multimap <double, int> orbMap;
    
    int alphaIndex = 0;
    
    // scan orbitals from left to right, pick out orbitals which have smallest v_1(i, i) integrals
    
    for (int i=left; i<right; ++i)
      {
	if (dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[i]]) {
	  orbitalList.push_back(alphaIndex);
	  alphaIndex++;
	  continue;
	}
	orbMap.insert(pair<double, int>(v_1[integralIndex](2*i, 2*i), alphaIndex));
	++alphaIndex;
      }
    
    for (multimap<double, int>::iterator m=orbMap.begin (); m!=orbMap.end (); ++m)
      orbitalList.push_back(m->second);
    

    // Now convert permutation lists to energy ordered form
    for (int i=0; i<abList.size (); ++i)
      ConvertList (abList[i], orbitalList);
    
    
    Slater last_det;
    // Now catenate to make Slaters
    for (int i=0; i<abList.size(); ++i)
      {
	std::vector<bool> lbuffer(dmrginp.spatial_to_spin()[left], 0);
	std::vector<bool> rbuffer(dmrginp.spatial_to_spin()[edge] - dmrginp.spatial_to_spin()[right], 0);
	
	std::vector<bool> orbs(dmrginp.spatial_to_spin()[right] - dmrginp.spatial_to_spin()[left], 0);
	int alphaI = 0;
	int betaI = 0;
	
	for (int orbI=0; orbI<(right-left); orbI++)
	{
	  if (abList[i][alphaI]) {
	    orbs[dmrginp.spatial_to_spin()[orbI]] = 1;
	    if (betaI < min(na, nb)) {
	      orbs[ dmrginp.spatial_to_spin()[orbI]+1] = 1;
	      ++betaI;
	    }
	  }
	  ++alphaI;
	}
	
	std::vector<bool> tmp = lbuffer;
	for (std::vector<bool>::iterator it = orbs.begin(); it!=orbs.end(); ++it) tmp.push_back(*it);
	for (std::vector<bool>::iterator it = rbuffer.begin(); it!=rbuffer.end(); ++it) tmp.push_back(*it);
	Slater new_det = Slater (Orbstring (tmp));
	map<Slater, double> m;
	m[new_det] = 1.0;
	last_det = new_det;

	if(sym.getirrep() == AbelianSymmetryOf(new_det).getirrep())
	  s.push_back ( Csf(m, n, SpinSpace(sp), sp, IrrepVector(sym.getirrep(), 0)) );
      }
    
    if (s.size() == 0) {
      numtries++;
      doAgain = true;
      abList.clear();
      numberOfGuesses = dmrginp.guess_permutations();
    }
    else
      doAgain = false;
  }
  
  return s;
}

std::vector<Csf> Csf::distributeNonSpinAdapted (const int n, const int sp, const IrrepVector &sym, const int left, const int right, const int edge, int integralIndex)
{
  std::vector<Csf> s;
  
  int na = 0;
  int nb = 0;
  na = (n + sp) / 2;
  nb = (n - sp) / 2;

  // let's count how many left and right orbitals there actually are
  int actualOrbs = right - left;
  int actualNA, actualNB;
  actualNA = actualNB = 0;
  for (int i=left; i<right; ++i)
    {
      if((SpinOf(i) == 1))
	++actualNA;
      else
	++actualNB;
    }

  // cannot form this combination
  if (na > actualNA || nb > actualNB || na < 0 || nb < 0) 
    {
      return s;
    }

  // now make orbital occupancies
  std::vector<int> alphaSeed (actualNA);
  std::vector<int> betaSeed (actualNB);
  for (int i = 0; i < na; ++i)
    alphaSeed [i] = -1;
  for (int i = 0; i < nb; ++i)
    betaSeed [i] = -1;

  // okay now make some permutations (at most GUESS_PERMUTATIONS)
  std::vector< std::vector<int> > alphaList; alphaList.push_back (alphaSeed);
  std::vector< std::vector<int> > betaList; betaList.push_back (betaSeed);
  int numberOfGuesses = dmrginp.guess_permutations();
  while (next_permutation (alphaSeed.begin (), alphaSeed.end ()) && numberOfGuesses--)      
    alphaList.push_back (alphaSeed);
  numberOfGuesses = dmrginp.guess_permutations();;
  while (next_permutation (betaSeed.begin (), betaSeed.end ()) && numberOfGuesses--)
    betaList.push_back (betaSeed);
  

  // okay - we have all the permutations. 
  // make list of available alpha and beta orbitals in energy ordering.
  std::vector<int> orbitalAlphaList;
  std::vector<int> orbitalBetaList;

  multimap <double, int> alphaMap;
  multimap <double, int> betaMap;

  int alphaIndex = 0;
  int betaIndex = 0;

  // scan orbitals from left to right, pick out orbitals which are occupied in the hartree-fock reference
  for (int i=left; i<right; ++i)
  {
    if (dmrginp.hf_occupancy()[i])
    {
      if ((SpinOf(i) == 1)) // 0 to take into account case of nospin
      {
	orbitalAlphaList.push_back(alphaIndex);
	++alphaIndex;
      }
      else
      {
	orbitalBetaList.push_back(betaIndex);
	++betaIndex;
      }
    }
    else if ((SpinOf(i) == 1)) // not HF orb, but alpha orb
    {
      alphaMap.insert(pair<double, int>(v_1[integralIndex](i, i), alphaIndex));
      ++alphaIndex;
    }
    else // beta orb
    {
      betaMap.insert(pair<double, int>(v_1[integralIndex](i, i), betaIndex));
      ++betaIndex;
    }
  }

  for (multimap<double, int>::iterator m=alphaMap.begin (); m!=alphaMap.end (); ++m)
    orbitalAlphaList.push_back(m->second);
  for (multimap<double, int>::iterator m=betaMap.begin (); m!=betaMap.end (); ++m)
    orbitalBetaList.push_back(m->second);

  // Now convert permutation lists to energy ordered form
  for (int i=0; i<alphaList.size (); ++i)
    ConvertList (alphaList[i], orbitalAlphaList);
  for (int i=0; i<betaList.size (); ++i)
    ConvertList (betaList[i], orbitalBetaList);
  // Now catenate to make Slaters
  for (int i=0; i<alphaList.size(); ++i)
    for (int j=0; j<betaList.size(); ++j)
    {
      std::vector<bool> lbuffer(left, 0);
      std::vector<bool> rbuffer(edge - right, 0);
      
      std::vector<bool> orbs(right - left);
      int alphaI = 0;
      int betaI = 0;
      for (int orbI=0; orbI<right-left; ++orbI)
      {
	if (SpinOf(orbI + left) == 1)
	{
	  if (alphaList[i][alphaI]) orbs[orbI] = 1;
	  ++alphaI;
	}
	else
	{
	  if (betaList[j][betaI]) orbs[orbI] = 1;
	  ++betaI;
	}
      }
      std::vector<bool> tmp = lbuffer;
      for (std::vector<bool>::iterator it = orbs.begin(); it!=orbs.end(); ++it) tmp.push_back(*it);
      for (std::vector<bool>::iterator it = rbuffer.begin(); it!=rbuffer.end(); ++it) tmp.push_back(*it);
      Slater new_det = Slater (Orbstring (tmp));
      map<Slater, double> m;
      m[new_det] = 1.0;
      if(sym.getirrep() == AbelianSymmetryOf(new_det).getirrep())
	s.push_back (Csf(m,n,SpinSpace(sp), sp, IrrepVector(sym.getirrep(),0)));
    }
  return s;
}


vector< vector<int> > SpinAdapted::CSFUTIL::generate_partitions(int n) 
{
  //n is the irrep of the symmetry element.
  if (Symmetry::sizeofIrrep(n) == 1) {
    pout << "cannot generate partition of irrep which a single row "<<endl;
    exit(0);
  }

  vector<int> thispartition(1, n);
  vector< vector<int> > result(1, thispartition);
  if (n == 4 || n==5) {
    return result;
  }
  else {

    for (int irrep = 4; irrep <n ; irrep++) {
      int remainder = (IrrepSpace(n)+IrrepSpace(irrep))[0].getirrep();
      vector< vector<int> > partitions = generate_partitions(irrep);

      for (int i=0; i<partitions.size(); i ++) {
	vector<int> ithpartition = partitions[i];
	if (ithpartition[0] <= remainder) {
	  vector<int>::iterator it = ithpartition.begin();
	  ithpartition.insert(it, remainder);
	  result.push_back(ithpartition);
	}
      }
    }
    return result;
  }

}

