/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include <IntegralMatrix.h>
#include "pario.h"
#include "screen.h"
#include "global.h"
using namespace std;

namespace SpinAdapted{

vector<int, std::allocator<int> > screened_d_indices(const vector<int, std::allocator<int> >& indices,
			       const vector<int, std::allocator<int> >& interactingix,
			       const OneElectronArray& onee, const TwoElectronArray& twoe, double thresh) {
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
    if (dmrginp.use_partial_two_integrals() || screen_d_interaction(indices[i], interactingix, onee, twoe, thresh))
      screened_indices.push_back(indices[i]);
  //pout << "\t\t\tnumber of significant d and d_comp indices: " << screened_indices.size() << endl;
  return screened_indices;
}

bool screen_d_interaction(int index, const vector<int, std::allocator<int> >& interactingix,
			  const OneElectronArray& onee, const TwoElectronArray& twoe, double thresh) {
  if(dmrginp.spinAdapted()) {
    for (int i = 0; i < interactingix.size(); ++i) {
      const int ix = interactingix[i];
      int xl = index;
      for (int ixx = dmrginp.spatial_to_spin(ix); ixx <dmrginp.spatial_to_spin(ix+1); ixx++)
	for (int lxx = dmrginp.spatial_to_spin(xl); lxx <dmrginp.spatial_to_spin(xl+1); lxx++)
	  if (fabs(onee(lxx, ixx)) >= thresh)
	    return true;
    }
    
    for (int i = 0; i < interactingix.size(); ++i)
      for (int j = 0; j < interactingix.size(); ++j)
	for (int k = 0; k < interactingix.size(); ++k)
	  {
	    const int ix = interactingix[i];
	    const int jx = interactingix[j];
	    const int kx = interactingix[k];
	    int xl = index;
	    for (int ixx = dmrginp.spatial_to_spin(ix); ixx <dmrginp.spatial_to_spin(ix+1); ixx++)
	      for (int jxx = dmrginp.spatial_to_spin(jx); jxx <dmrginp.spatial_to_spin(jx+1); jxx++)
		for (int kxx = dmrginp.spatial_to_spin(kx); kxx <dmrginp.spatial_to_spin(kx+1); kxx++)
		  for (int lxx = dmrginp.spatial_to_spin(xl); lxx <dmrginp.spatial_to_spin(xl+1); lxx++)
		    if (fabs(twoe(lxx,ixx,jxx,kxx)) >= thresh)
		      return true;
	  }
    
    if(interactingix.size() == 0)
      return true;
    else
      return false;
  } else {
    for (int i = 0; i < interactingix.size(); ++i){
      const int ix = interactingix[i];
      int xl = index;
      if (fabs(onee(xl, ix)) >= thresh)
	return true;
    }
    
    for (int i = 0; i < interactingix.size(); ++i)
      for (int j = 0; j < interactingix.size(); ++j)
	for (int k = 0; k < interactingix.size(); ++k)
	  {
	    const int ix = interactingix[i];
	    const int jx = interactingix[j];
	    const int kx = interactingix[k];
	    int xl = index;
	    if (fabs(twoe(xl,ix,jx,kx)) >= thresh)
	      return true;
	  }
    
    return (interactingix.size() == 0);
  }
}

vector<int, std::allocator<int> > screened_cddcomp_indices(const vector<int, std::allocator<int> >& otherindices,
				     const vector<int, std::allocator<int> >& selfindices,
				     const OneElectronArray& onee, 
				     const TwoElectronArray& twoe, double thresh)
{
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < otherindices.size(); ++i)
    if (dmrginp.use_partial_two_integrals() || screen_cddcomp_interaction(otherindices[i], selfindices, onee, twoe, thresh))
      screened_indices.push_back(otherindices[i]);
  //pout << "\t\t\tnumber of significant cdd and cdd_comp indices: " << screened_indices.size() << endl;
  return screened_indices;
}

bool screen_cddcomp_interaction(int otherindex, const vector<int, std::allocator<int> >& selfindices,
				const OneElectronArray& onee, 
				const TwoElectronArray& twoe, double thresh)
{
  if(!dmrginp.spinAdapted()) {
    for (int i = 0; i < selfindices.size(); ++i)
      {
	const int ix = selfindices[i];
	int lx = otherindex;
	if (fabs(onee(lx, ix)) >= thresh)
	  return true;      
      }
    
    for (int i = 0; i < selfindices.size(); ++i)
      for (int j = 0; j < selfindices.size(); ++j)
	for (int k = 0; k < selfindices.size(); ++k)
	  {
	    const int ix = selfindices[i];
	    const int jx = selfindices[j];
	    const int kx = selfindices[k];
	    int lx = otherindex;
	    if (fabs(twoe(lx,ix,jx,kx)) >= thresh)
	      return true;
	  }
    return (selfindices.size() == 0);
  }
  else {
    for (int i = 0; i < selfindices.size(); ++i)
      {
	const int ix = selfindices[i];
	int xl = otherindex;
	for (int ixx = dmrginp.spatial_to_spin(ix); ixx <dmrginp.spatial_to_spin(ix+1); ixx++)
	  for (int lxx = dmrginp.spatial_to_spin(xl); lxx <dmrginp.spatial_to_spin(xl+1); lxx++) {
	    if (fabs(onee(lxx, ixx)) >= thresh)
	      return true;
	  }
      }
    
    for (int i = 0; i < selfindices.size(); ++i)
      for (int j = 0; j < selfindices.size(); ++j)
	for (int k = 0; k < selfindices.size(); ++k)
	  {
	    const int ix = selfindices[i];
	    const int jx = selfindices[j];
	    const int kx = selfindices[k];
	    int xl = otherindex;
	    for (int ixx = dmrginp.spatial_to_spin(ix); ixx <dmrginp.spatial_to_spin(ix+1); ixx++)
	      for (int jxx = dmrginp.spatial_to_spin(jx); jxx <dmrginp.spatial_to_spin(jx+1); jxx++)
		for (int kxx = dmrginp.spatial_to_spin(kx); kxx <dmrginp.spatial_to_spin(kx+1); kxx++)
		  for (int lxx = dmrginp.spatial_to_spin(xl); lxx <dmrginp.spatial_to_spin(xl+1); lxx++)
		    if (fabs(twoe(lxx,ixx,jxx,kxx)) >= thresh)
		      return true;
	  }
    return (selfindices.size() == 0);
  }

    
}

/**
 * from a list of (sorted) indices in a given block
 * returns the screened pair indices for the operators
 * cd and cdcomp that
 * interact
 * with the other block indices (interactingix)
 */
vector<pair<int, int> > screened_cd_indices(const vector<int, std::allocator<int> >& indices,
					    const vector<int, std::allocator<int> >& interactingix,
					    const TwoElectronArray& twoe, double thresh)
{
  vector<pair<int, int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
    for (int j = 0; j <= i; ++j)
      if (dmrginp.use_partial_two_integrals() || screen_cd_interaction(indices[i], indices[j], interactingix, twoe, thresh))
	    screened_indices.push_back(make_pair(indices[i], indices[j]));
  return screened_indices;
}

/**
 * from a list of (sorted) indices in a given block
 * returns the screened pair indices for the operators
 * dd and ddcomp that
 * interact
 * with the other block indices (interactingix)
 */
vector<pair<int, int> > screened_dd_indices(const vector<int, std::allocator<int> >& indices,
					    const vector<int, std::allocator<int> >& interactingix,
					    const TwoElectronArray& twoe, double thresh)
{
  vector<pair<int, int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
  for (int j = 0; j <= i; ++j)
    if (dmrginp.use_partial_two_integrals() || screen_dd_interaction(indices[i], indices[j], interactingix, twoe, thresh))
	screened_indices.push_back(make_pair(indices[i], indices[j]));
  return screened_indices;
}

/**
 * given two indices i and j, determine
 * whether we should build c+i dj
 * or the complementary operator for c+i dj
 * by looking at the integrals of the complementary
 * operator
 * interactingix are the indices that we are summing over
 * (i.e. the indices in the block that we are interacting with)
 */
bool screen_cd_interaction(int ci, int dj, const vector<int, std::allocator<int> >& interactingix,
			   const TwoElectronArray& twoe, double thresh)
{
  if (dmrginp.spinAdapted()) {
    int ninter = interactingix.size();
    
    double twoeterm = 0.;
    for (int k = 0; k < ninter; ++k)
      for (int l = 0; l < ninter; ++l)
	{
	  int xk = interactingix[k];
	  int xl = interactingix[l];
	  for (int cix = dmrginp.spatial_to_spin(ci); cix <dmrginp.spatial_to_spin(ci+1); cix++)
	    for (int djx = dmrginp.spatial_to_spin(dj); djx <dmrginp.spatial_to_spin(dj+1); djx++)
	      for (int kxx = dmrginp.spatial_to_spin(xk); kxx <dmrginp.spatial_to_spin(xk+1); kxx++)
		for (int lxx = dmrginp.spatial_to_spin(xl); lxx <dmrginp.spatial_to_spin(xl+1); lxx++)
		  if (fabs(twoe(cix, kxx, lxx, djx))>=thresh || fabs(twoe(kxx, cix, lxx, djx)) >= thresh)
		    return true; // there is a significant integral joining the two regions
	}
    return (ninter == 0);
  }
  else {
    int ninter = interactingix.size();
    
    for (int k = 0; k < ninter; ++k)
      for (int l = 0; l < ninter; ++l)
	{
	  int kx = interactingix[k];
	  int lx = interactingix[l];
	  if (fabs(twoe(ci, kx, lx, dj))>=thresh || fabs(twoe(kx, ci, lx, dj)) >= thresh)
	    return true; // there is a significant integral joining the two regions
	}
    return (ninter == 0);
  }
}

bool screen_dd_interaction(int ci, int cj, const vector<int, std::allocator<int> >& interactingix,
			   const TwoElectronArray& twoe, double thresh)
{
  if(dmrginp.spinAdapted()) {
    int ninter = interactingix.size();
    
    for (int k = 0; k < ninter; ++k)
      for (int l = 0; l < ninter; ++l)
	{
	  int xk = interactingix[k];
	  int xl = interactingix[l];
	  for (int cix = dmrginp.spatial_to_spin(ci); cix <dmrginp.spatial_to_spin(ci+1); cix++)
	    for (int cjx = dmrginp.spatial_to_spin(cj); cjx <dmrginp.spatial_to_spin(cj+1); cjx++)
	      for (int kxx = dmrginp.spatial_to_spin(xk); kxx <dmrginp.spatial_to_spin(xk+1); kxx++)
		for (int lxx = dmrginp.spatial_to_spin(xl); lxx <dmrginp.spatial_to_spin(xl+1); lxx++)
		  if (fabs(twoe(cix, cjx, kxx, lxx))>=thresh)
		    return true; // there is a significant integral joining the two regions
	}
    return (ninter == 0);
  }
  else  {
    int ninter = interactingix.size();
    
    for (int k = 0; k < ninter; ++k)
      for (int l = 0; l < ninter; ++l)
	{
	  int kx = interactingix[k];
	  int lx = interactingix[l];
	  if (fabs(twoe(ci, cj, kx, lx))>=thresh)
	    return true; // there is a significant integral joining the two regions
	}
    return (ninter == 0);
  }
}

// these are for BCS type calculations

std::vector<int, std::allocator<int> > screened_d_indices(const std::vector<int, std::allocator<int> >& indices, const std::vector<int, std::allocator<int> >& interactingix, const OneElectronArray& onee, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh) {
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
    if (dmrginp.use_partial_two_integrals() || screen_d_interaction(indices[i], interactingix, onee, twoe, vcc, vcccc, vcccd, thresh))
      screened_indices.push_back(indices[i]);
  return screened_indices;
}

bool screen_d_interaction(int index, const std::vector<int, std::allocator<int> >& interactingix, const OneElectronArray& onee, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh) {
  if (dmrginp.spinAdapted()) {
    cout << "BCS with spin adaption not implemented!" << endl;
    abort();
  } else {
    for (int i = 0; i < interactingix.size(); ++i) {
      const int ix = interactingix[i];
      int xl = index;
      if (fabs(onee(xl, ix)) >= thresh || fabs(vcc(xl, ix)) >= thresh || fabs(vcc(ix, xl)) >= thresh)
        return true;
    }

    for (int i = 0; i < interactingix.size(); ++i)
    for (int j = 0; j < interactingix.size(); ++j)
	for (int k = 0; k < interactingix.size(); ++k) {
      const int ix = interactingix[i];
	  const int jx = interactingix[j];
	  const int kx = interactingix[k];
      int xl = index;
      if (fabs(twoe(xl, ix, jx, kx)) >= thresh || fabs(vcccd(xl, ix, jx, kx)) >= thresh || fabs(vcccd(ix, jx, kx, xl)) > thresh || fabs(vcccc(xl, ix, jx, kx)) > thresh)
        return true;
    }
    return (interactingix.size() == 0);
  }
}

std::vector<std::pair<int, int> > screened_cd_indices(const std::vector<int, std::allocator<int> >& indices, const std::vector<int, std::allocator<int> >& interactingix, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh, bool symm, bool comp) {
  vector<pair<int, int> > screened_indices;

  for (int i = 0; i < indices.size(); ++i) {
    int max_j = symm ? i: indices.size()-1;
    for (int j = 0; j <= max_j; ++j)
      if (dmrginp.use_partial_two_integrals() || screen_cd_interaction(indices[i], indices[j], interactingix, twoe, vcc, vcccc, vcccd, thresh, symm, comp))
	    screened_indices.push_back(make_pair(indices[i], indices[j]));
  }
  return screened_indices;
}

bool screen_cd_interaction(int ci, int dj, const std::vector<int, std::allocator<int> >& interactingix, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh, bool symm, bool comp) {
  if (dmrginp.spinAdapted()) {
    cout << "BCS with spin adaption not implemented!" << endl;
    abort();
  } else {
    int ninter = interactingix.size();
    for (int k = 0; k < ninter; ++k)
    for (int l = 0; l < ninter; ++l) {
      int kx = interactingix[k];
	  int lx = interactingix[l];
      if (symm && !comp) {
	    if (fabs(twoe(ci, kx, lx, dj))>= thresh || fabs(twoe(kx, ci, lx, dj)) >= thresh || fabs(v_cccd(ci, kx, lx, dj)) >= thresh || fabs(v_cccd(dj, kx, lx, ci)) >= thresh)
          return true;
      } else if (symm && comp) {
        if (fabs(twoe(ci, kx, lx, dj))>=thresh || fabs(twoe(kx, ci, lx, dj)) >= thresh)
          return true;
      } else { // !symm
        if (fabs(v_cccd(ci, kx, lx, dj)) >= thresh || fabs(v_cccd(dj, kx, lx, ci)) >= thresh)
          return true;
      }
    }
    return (ninter == 0);
  }
}

std::vector<std::pair<int, int> > screened_dd_indices(const std::vector<int, std::allocator<int> >& indices, const std::vector<int, std::allocator<int> >& interactingix, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh) {
  vector<pair<int, int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
  for (int j = 0; j <= i; ++j)
    if (dmrginp.use_partial_two_integrals() || screen_dd_interaction(indices[i], indices[j], interactingix, twoe, vcc, vcccc, vcccd, thresh))
	  screened_indices.push_back(make_pair(indices[i], indices[j]));
  return screened_indices;
}

bool screen_dd_interaction(int ci, int cj, const std::vector<int, std::allocator<int> >& interactingix, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh) {
  if (dmrginp.spinAdapted()) {
    cout << "BCS with spin adaption not implemented!" << endl;
    abort();
  } else {
    int ninter = interactingix.size();
    
    for (int k = 0; k < ninter; ++k)
    for (int l = 0; l < ninter; ++l) {
      int kx = interactingix[k];
	  int lx = interactingix[l];
	  if (fabs(twoe(ci, cj, kx, lx))>=thresh || fabs(vcccd(ci, cj, k, l)) >= thresh || fabs(vcccc(ci, cj, kx, lx)) >= thresh)
        return true;
    }
    return (ninter == 0);
  }
}

std::vector<int, std::allocator<int> > screened_cddcomp_indices(const std::vector<int, std::allocator<int> >& otherindices, const std::vector<int, std::allocator<int> >& selfindices, const OneElectronArray& onee, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh) {
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < otherindices.size(); ++i)
    if (dmrginp.use_partial_two_integrals() || screen_cddcomp_interaction(otherindices[i], selfindices, onee, twoe, vcc, vcccc, vcccd, thresh))
      screened_indices.push_back(otherindices[i]);
  //pout << "\t\t\tnumber of significant cdd and cdd_comp indices: " << screened_indices.size() << endl;
  return screened_indices;
}

bool screen_cddcomp_interaction(int otherindex, const std::vector<int, std::allocator<int> >& selfindices, const OneElectronArray& onee, const TwoElectronArray& twoe, const PairArray& vcc, const CCCCArray& vcccc, const CCCDArray& vcccd, double thresh) {
  if (dmrginp.spinAdapted()) {
    cout << "BCS with spin adaption not implemented!" << endl;
    abort();
  } else {
    for (int i = 0; i < selfindices.size(); ++i) {
      const int ix = selfindices[i];
      int lx = otherindex;
      if (fabs(onee(lx, ix)) >= thresh || fabs(vcc(lx, ix)) >= thresh || fabs(vcc(ix, lx)) >= thresh)
        return true;
    }
    for (int i = 0; i < selfindices.size(); ++i)
    for (int j = 0; j < selfindices.size(); ++j)
	for (int k = 0; k < selfindices.size(); ++k) {
      const int ix = selfindices[i];
	  const int jx = selfindices[j];
	  const int kx = selfindices[k];
      int lx = otherindex;
      if (fabs(twoe(lx, ix, jx, kx)) >= thresh || fabs(vcccd(lx, ix, jx, kx)) >= thresh || fabs(vcccd(ix, jx, kx, lx)) > thresh || fabs(vcccc(lx, ix, jx, kx)) > thresh)
        return true;
    }
    return (selfindices.size() == 0);
  }
}

} // namespace SpinAdapted
