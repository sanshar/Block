#ifndef SPIN_SCREEN_HEADER
#define SPIN_SCREEN_HEADER 
#include <vector>
#include <utility>
#include <IntegralMatrix.h>

namespace SpinAdapted{
/**
 * from a list of (sorted) indices in a given block
 * returns the screened pair indices for the operators
 * cd and cdcomp that 
 * interact
 * with the other block indices (interactingix)
 */
std::vector<std::pair<int, int> > screened_cd_indices(const std::vector<int, std::allocator<int> >& indices,
						      const std::vector<int, std::allocator<int> >& interactingix,
						      const TwoElectronArray& twoe, double thresh);
/**
 * from a list of (sorted) indices in a given block
 * returns the screened pair indices for the operators
 * dd and ddcomp that 
 * interact
 * with the other block indices (interactingix)
 */
std::vector<std::pair<int, int> > screened_dd_indices(const std::vector<int, std::allocator<int> >& indices,
						      const std::vector<int, std::allocator<int> >& interactingix,
						      const TwoElectronArray& twoe, double thresh);


/**
 * given two indices i and j, determine
 * whether we should build c+i dj
 * or the complementary operator for c+i dj
 * by looking at the integrals of the complementary
 * operator
 * interactingix are the indices that we are summing over
 * (i.e. the indices in the block that we are interacting with)
 */
bool screen_cd_interaction(int ci, int dj, const std::vector<int, std::allocator<int> >& interactingix,
			   const TwoElectronArray& twoe, double thresh);

/* see comment for screen_cd_interaction */
bool screen_dd_interaction(int di, int dj, const std::vector<int, std::allocator<int> >& interactingix,
			   const TwoElectronArray& twoe, double thresh);

std::vector<int, std::allocator<int> > screened_d_indices(const std::vector<int, std::allocator<int> >& indices,
				    const std::vector<int, std::allocator<int> >& interactingix,
				    const OneElectronArray& onee, double thresh);

bool screen_d_interaction(int index, const std::vector<int, std::allocator<int> >& interactingix,
			  const OneElectronArray& onee, double thresh);

std::vector<int, std::allocator<int> > screened_cddcomp_indices(const std::vector<int, std::allocator<int> >& otherindices,
					  const std::vector<int, std::allocator<int> >& selfindices,
					  const OneElectronArray& onee, 
					  const TwoElectronArray& twoe, double thresh);

bool screen_cddcomp_interaction(int otherindex, const std::vector<int, std::allocator<int> >& selfindices,
				const OneElectronArray& onee, 
				const TwoElectronArray& twoe, double thresh);
}

#endif
