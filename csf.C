/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is free software: you can redistribute it and/or modify         
it under the terms of the GNU General Public License as published by         
the Free Software Foundation, either version 3 of the License, or            
(at your option) any later version.                                          
                                                                             
This program is distributed in the hope that it will be useful,              
but WITHOUT ANY WARRANTY; without even the implied warranty of               
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
GNU General Public License for more details.                                 
                                                                             
You should have received a copy of the GNU General Public License            
along with this program.  If not, see <http://www.gnu.org/licenses/>.        

This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "csf.h"
#include "SpinQuantum.h"
#include "global.h"
#include "couplingCoeffs.h"
#include "boost/shared_ptr.hpp"

SpinAdapted::Csf::Csf( const map<Slater, double>& p_dets, const int p_n, const int p_S, const int p_Sz, const IrrepVector p_irrep) : det_rep(p_dets), n(p_n), S(p_S), Sz(p_Sz), irrep(p_irrep)
{
  map<Slater, double>::iterator it = det_rep.begin();
  for (; it!= det_rep.end(); it++)
    if ((*it).first.n != n || (*it).first.Sz != Sz )
      {
	cout << it->first<<endl;
	cout << (*it).first.n<<" "<<n<<endl;
	cout << (*it).first.Sz<<" "<<Sz<<endl;
	cout<<" all slaters in Csf dont have the same n, spin or symmetry"<<endl;
	cout << *this<<endl;
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

double SpinAdapted::csf_energy(const Csf& s)
{
  map<Slater, double>::const_iterator it = s.det_rep.begin();
  double energy = 0.0;
  for (; it != s.det_rep.end();it++)
    energy += it->second* det_energy(it->first);

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
  if (sym != "dinfh" || irrep <= 3) return;
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
      if (partitions.size() == 0) {cout << "big trouble"<<endl;exit(0);}
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
    for (int i = this->S-2; i>= this->S - 2*k; i-=2) {
      Csf sm;
      s.applySminus(sm);
      if (sm.size() != 0)
	ladder.push_back( sm);
      s = sm;
    }
  }
  return ladder;
}

void SpinAdapted::CSFUTIL::TensorProduct(Csf& lhs, Csf& rhs, vector< Csf >& output)
{
  vector< Csf> lhs_csfs = lhs.spinLadder(lhs.S);
  vector< Csf> rhs_csfs = rhs.spinLadder(rhs.S);

  int lirrep = lhs.sym_is().getirrep();
  int rirrep = rhs.sym_is().getirrep(); 
  vector<IrrepSpace> symvec = lhs.sym_is()+rhs.sym_is();
  //now generate all the output vectors |J,J> with J = |j1-j2|,...,|j1+j2|
  //the angular momentum is m = m1+m2, |m1-m2|
  for (int l = 0; l<symvec.size(); l++)
  {
    int irrep = symvec[l].getirrep();
    int row = Symmetry::sizeofIrrep(irrep)-1;//***rows go from 0 to 1

    for (int J = abs(lhs.S-rhs.S); J<= lhs.S+rhs.S; J+=2)
    {

      //generating |J1 J2 J J> vector
      map<Slater, double> dets;
      for (int j1m1 = 0; j1m1<lhs_csfs.size(); j1m1++)
	for (int j2m2 = 0; j2m2<rhs_csfs.size(); j2m2++)
	{
	  int J1 = lhs.S, M1 = lhs_csfs[j1m1].Sz, J2 = rhs.S, M2 = rhs_csfs[j2m2].Sz;
	
	  double clebsg = cg(J1, J2, J, M1, M2, J);
	  double clebdinfh = Symmetry::spatial_cg(lirrep, rirrep, irrep, lhs_csfs[j1m1].irrep.getrow(), rhs_csfs[j2m2].irrep.getrow(), row);
	  if( abs(clebsg) < 1.0e-12 || abs(clebdinfh) < 1.00e-14)
	    continue;
	  lhs_csfs[j1m1].outerProd(rhs_csfs[j2m2], clebsg*clebdinfh, dets); 
	}
      
      if (dets.size() == 0)
	continue;
      Csf c(dets, rhs.n+lhs.n, J, J, IrrepVector(symvec[l].getirrep(), 0)); c.normalize();
      output.push_back(c);
    }
  }
}

std::vector<SpinAdapted::Csf > SpinAdapted::CSFUTIL::spinfockstrings(const std::vector<int>& orbs)
{
  int latticelen = Slater().size();
  /* first extract out alpha orbs and beta orbs from orbs */
  std::vector<bool> aorbs, borbs;

  std::vector<Csf > singleSiteCsf;
  std::vector<int> numcsfs;
  
  for (int i=0; i<orbs.size(); i++)
  {
    int I = orbs[i];
    IrrepSpace Irrep = SymmetryOfSpatialOrb(orbs[i]); 
    IrrepSpace Irrep2 = (Irrep+Irrep).back(); 
    IrrepSpace Irrep3 = *(Irrep2+Irrep).begin(); 
    IrrepSpace Irrep4 = *(Irrep3+Irrep).begin(); 

    int index = Symmetry::sizeofIrrep(Irrep.getirrep())==1 ? 0 : 2;

    std::vector<bool> occ_rep1(Slater().size(),0), 
                      occ_rep2(Slater().size(),0), 
                      occ_rep3(Slater().size(),0);

    occ_rep2[dmrginp.spatial_to_spin()[I]+index] = 1;
    occ_rep3[dmrginp.spatial_to_spin()[I]+index] = 1; occ_rep3[dmrginp.spatial_to_spin()[I]+1+index] = 1;

    Slater s1(occ_rep1, 1), s2(occ_rep2, 1), s3(occ_rep3, 1);
    map<Slater, double> m1, m2, m3, m4, m5, m6, m7;
    m1[s1]= 1.0; m2[s2]= 1.0; m3[s3]= 1.0;
    singleSiteCsf.push_back( Csf(m1, 0, 0, 0, IrrepVector(0,0))); //0,0,0
    singleSiteCsf.push_back( Csf(m2, 1, 1, 1, IrrepVector(Irrep.getirrep(), index/2))); //1,1,L
    singleSiteCsf.push_back( Csf(m3, 2, 0, 0, IrrepVector(Irrep2.getirrep(), index/2))); //2,0,2L
    if (Symmetry::sizeofIrrep(Irrep.getirrep()) == 1) {
      numcsfs.push_back(3);
      continue;
    }

    std::vector<bool> occ_rep4(Slater().size(),0), 
                      occ_rep5(Slater().size(),0), 
                      occ_rep6(Slater().size(),0), 
                      occ_rep7(Slater().size(),0),
                      occ_rep8(Slater().size(),0);

    occ_rep4[dmrginp.spatial_to_spin()[I]] = 1; occ_rep4[dmrginp.spatial_to_spin()[I]+3] = 1;
    occ_rep5[dmrginp.spatial_to_spin()[I]+1] = 1; occ_rep5[dmrginp.spatial_to_spin()[I]+2] = 1;
    m4[Slater(occ_rep4, 1)] = -1.0/sqrt(2.0); m4[Slater(occ_rep5, 1)] = 1.0/sqrt(2.0);
    singleSiteCsf.push_back( Csf(m4, 2, 0, 0, IrrepVector((Irrep+Irrep)[0].getirrep(), 0))); //2,0,0

    occ_rep6[dmrginp.spatial_to_spin()[I]] = 1; occ_rep6[dmrginp.spatial_to_spin()[I] +2] = 1;
    m5[Slater(occ_rep6, 1)] = 1.0;
    singleSiteCsf.push_back( Csf(m5, 2, 2, 2, IrrepVector((Irrep+Irrep)[1].getirrep(), 0))); //2,2,0

    occ_rep7[dmrginp.spatial_to_spin()[I]+2] = 1; occ_rep7[dmrginp.spatial_to_spin()[I]+3] = 1; occ_rep7[dmrginp.spatial_to_spin()[I]] = 1;
    m6[Slater(occ_rep7, 1)] = 1.0;
    singleSiteCsf.push_back( Csf(m6, 3, 1, 1, IrrepVector(Irrep3.getirrep(), 1))); //3,1,L

    occ_rep8[dmrginp.spatial_to_spin()[I]] = 1; occ_rep8[dmrginp.spatial_to_spin()[I]+1] = 1;
    occ_rep8[dmrginp.spatial_to_spin()[I] +2] = 1; occ_rep8[dmrginp.spatial_to_spin()[I] +3] = 1;
    m7[Slater(occ_rep8, 1)] = 1.0;
    singleSiteCsf.push_back( Csf(m7, 4, 0, 0, IrrepVector(Irrep4.getirrep(), 0))); //4,0,0

    //change the order of csf2 and csf3 so they are ordered 200 and then 202L
    int len = singleSiteCsf.size();
    Csf temp = singleSiteCsf[len-4];
    singleSiteCsf[len-4] = singleSiteCsf[len-5];
    singleSiteCsf[len-5] = temp;
    numcsfs.push_back(7);
  }

  std::vector< Csf > output, prevoutput;
  for (int i=0;i<numcsfs[0];i++) {
    output.push_back(singleSiteCsf[i]);
  }
  int csfindex = numcsfs[0];

  for (int i2=1; i2<orbs.size(); i2++) {

    int outputsize = output.size();
    for (int i=0;i<outputsize;i++) 
      prevoutput.push_back(output[i]);
    output.clear();


    for (int j=0; j<prevoutput.size(); j++)
      for (int k=0; k<numcsfs[i2]; k++) {
	CSFUTIL::TensorProduct(prevoutput[j], singleSiteCsf[csfindex+k], output);
      }
      

    prevoutput.clear();
    csfindex += numcsfs[i2];
  }	
  const Csfcompare csfcompare;
  sort(output.begin(), output.end());

  return output;
	    
}

std::vector< SpinAdapted::Csf > SpinAdapted::Csf::distribute (const int n, const int sp, const IrrepVector &sym, const int left, const int right, const int edge)
{
  std::vector< Csf > s;

  int na = 0;
  int nb = 0;
  na = (n + sp) / 2;
  nb = (n - sp) / 2;

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
  while(doAgain && numtries <5000)
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
	orbMap.insert(pair<double, int>(v_1(2*i, 2*i), alphaIndex));
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
	copy (orbs.begin(), orbs.end(), back_inserter(tmp));
	copy (rbuffer.begin(), rbuffer.end(), back_inserter(tmp));
	Slater new_det = Slater (Orbstring (tmp));
	map<Slater, double> m;
	m[new_det] = 1.0;
	last_det = new_det;

	if(sym.getirrep() == AbelianSymmetryOf(new_det).getirrep())
	  s.push_back ( Csf(m, n, sp, sp, IrrepVector(sym.getirrep(), 0)) );
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


vector< vector<int> > SpinAdapted::CSFUTIL::generate_partitions(int n) 
{
  //n is the irrep of the symmetry element.
  if (Symmetry::sizeofIrrep(n) == 1) {
    cerr << "cannot generate parition of irrep which a single row "<<endl;
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

