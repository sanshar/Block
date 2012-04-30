#include <IntegralMatrix.h>
#include "pario.h"
#include "screen.h"
using namespace std;

namespace SpinAdapted{

vector<int, std::allocator<int> > screened_d_indices(const vector<int, std::allocator<int> >& indices,
			       const vector<int, std::allocator<int> >& interactingix,
			       const OneElectronArray& onee, double thresh)
{
  vector<int, std::allocator<int> > screened_indices;
  for (int i = 0; i < indices.size(); ++i)
    if (screen_d_interaction(indices[i], interactingix, onee, thresh))
      screened_indices.push_back(indices[i]);
  //pout << "\t\t\tnumber of significant d and d_comp indices: " << screened_indices.size() << endl;
  return screened_indices;
}

bool screen_d_interaction(int index, const vector<int, std::allocator<int> >& interactingix,
			  const OneElectronArray& onee, double thresh)
{
  for (int i = 0; i < interactingix.size(); ++i)
    if (abs(onee(2*index, 2*interactingix[i])) >= thresh)
      return true;
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
  for (int i = 0; i < otherindices.size(); ++i)
    if (screen_cddcomp_interaction(otherindices[i], selfindices, onee, twoe, thresh))
      screened_indices.push_back(otherindices[i]);
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
    if (abs(onee(2*otherindex, 2*ix)) >= thresh)
	    return true;
  }

  for (int i = 0; i < selfindices.size(); ++i)
    for (int j = 0; j < selfindices.size(); ++j)
      for (int k = 0; k < selfindices.size(); ++k)
	    {
	      const int ix = selfindices[i];
	      const int jx = selfindices[j];
	      const int kx = selfindices[k];
	      if (abs(twoe(2*otherindex,2*ix,2*jx,2*kx)) >= thresh)
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
      if (screen_cd_interaction(indices[i], indices[j],
				    interactingix, twoe, thresh))
	        screened_indices.push_back(make_pair(indices[i], indices[j]));
  //pout << "\t\t\tnumber of significant cd and cd_comp indices: " << screened_indices.size() << endl;
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
      if (screen_dd_interaction(indices[i], indices[j],
				  interactingix, twoe, thresh))
	      screened_indices.push_back(make_pair(indices[i], indices[j]));
  //pout << "\t\t\tnumber of significant dd and dd_comp indices: " << screened_indices.size() << endl;
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
	    if (fabs(-twoe(2*ci, 2*xk, 2*xl, 2*dj) + 2*twoe(2*xk, 2*ci, 2*xl, 2*dj)) >= thresh ||
		fabs(twoe(2*ci, 2*xk, 2*xl, 2*dj))  >= thresh)
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
	    if (fabs(twoe(2*ci, 2*cj, 2*xk, 2*xl)) >= thresh)
	      return true;
    }
  if(ninter == 0)
    return true;
  else
    return false;
}
}
