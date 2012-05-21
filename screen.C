#include <IntegralMatrix.h>
#include "pario.h"
#include "screen.h"
#include "global.h"
using namespace std;

namespace SpinAdapted{

vector<int, std::allocator<int> > screened_d_indices(const vector<int, std::allocator<int> >& indices,
			       const vector<int, std::allocator<int> >& interactingix,
			       const OneElectronArray& onee, const TwoElectronArray& twoe, double thresh)
{
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
    if (screen_d_interaction(indices[i], interactingix, onee, twoe, thresh))
      screened_indices.push_back(indices[i]);
  //pout << "\t\t\tnumber of significant d and d_comp indices: " << screened_indices.size() << endl;
  return screened_indices;
}

bool screen_d_interaction(int index, const vector<int, std::allocator<int> >& interactingix,
			  const OneElectronArray& onee, const TwoElectronArray& twoe, double thresh)
{
  for (int i = 0; i < interactingix.size(); ++i){
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
}

vector<int, std::allocator<int> > screened_cddcomp_indices(const vector<int, std::allocator<int> >& otherindices,
				     const vector<int, std::allocator<int> >& selfindices,
				     const OneElectronArray& onee, 
				     const TwoElectronArray& twoe, double thresh)
{
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < otherindices.size(); ++i) {
    if (screen_cddcomp_interaction(otherindices[i], selfindices, onee, twoe, thresh))
      screened_indices.push_back(otherindices[i]);
  }
  //pout << "\t\t\tnumber of significant cdd and cdd_comp indices: " << screened_indices.size() << endl;
  return screened_indices;
}

bool screen_cddcomp_interaction(int otherindex, const vector<int, std::allocator<int> >& selfindices,
				const OneElectronArray& onee, 
				const TwoElectronArray& twoe, double thresh)
{
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
  if(selfindices.size() == 0)
    return true;
  else
    return false;
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
      if (screen_cd_interaction(indices[i], indices[j], interactingix, twoe, thresh))
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
      if (screen_dd_interaction(indices[i], indices[j], interactingix, twoe, thresh))
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
  if(ninter == 0)
    return true;
  else
    return false;
}

bool screen_dd_interaction(int ci, int cj, const vector<int, std::allocator<int> >& interactingix,
			   const TwoElectronArray& twoe, double thresh)
{
  int ninter = interactingix.size();

  double twoeterm = 0.;
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
  if(ninter == 0)
    return true;
  else
    return false;
}
}
