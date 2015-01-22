// Guess configuration generator as a member function of Input
// "hf_occ = integral" with symmetry

// Written by N.Nakatani, 2014.12.25

#include <cassert>
#include <vector>
#include <map>

#include "global.h"
#include "input.h"

// from Symmetry.C
extern array_2d<int> groupTable;

std::vector<int> SpinAdapted::Input::hfOccGenerator_ ()
{
  typedef std::multimap<double,int>::iterator mapIter;

  int irep = m_total_symmetry_number.getirrep();
  if (NonabelianSym) irep = 0; //with nonabelian symmetry this code will only work with Ag symmetry
  int spac_norbs = m_norbs/2;

  // Sorting orbitals by 1-el. Hamiltonian, i.e. h(i, i)
  std::vector<int> reorder(spac_norbs,0);
  {
    // NOTE: v_1 in lattice order
    std::multimap<double,int> orbMap;
    for(int i = 0; i < spac_norbs; ++i)
      orbMap.insert(std::make_pair(v_1[0](2*i, 2*i), i));

    mapIter it = orbMap.begin();
    for(int i = 0; i < spac_norbs; ++i, ++it)
      reorder[i] = it->second;
  }

  int refIrep = 0;

  // NOTE: m_spin_orbs_symmetry in lattice order
  std::vector<int> refOcc(m_norbs,0);
  for(int i = 0; i < m_alpha; ++i) {
    int ir = reorder[i];
    refOcc[2*ir] = 1;
    refIrep = Symmetry::add(refIrep, m_spin_orbs_symmetry[2*ir])[0];
    if(i < m_beta) {
      refOcc[2*ir+1] = 1;
      refIrep = Symmetry::add(refIrep, m_spin_orbs_symmetry[2*ir+1])[0];
    }
  }

  if(refIrep != irep) {
    std::vector<std::vector<int>> exOccs;
    // beta el. excitation
    // --- ---     --- ---
    // --- ---     --- ---
    // -o- ---     -o- ---
    // -o- --- ==> -o- -o-
    // -o- ---     -o- ---
    // -o- ---     -o- ---
    // -o- -o- ==> -o- ---
    // -o- -o-     -o- -o-
    // -o- -o-     -o- -o-
    //  a   b       a   b
    for(int i = 0; i < m_beta; ++i) {
      int ir = reorder[i];
      for(int j = m_beta; j < m_alpha; ++j) {
        int jr = reorder[j];
        int exIrep = Symmetry::add(m_spin_orbs_symmetry[2*ir+1],m_spin_orbs_symmetry[2*jr+1])[0];
        if(Symmetry::add(refIrep,exIrep)[0] == irep) {
          std::vector<int> exOcc = refOcc;
          std::swap(exOcc[2*ir+1],exOcc[2*jr+1]);
          exOccs.push_back(exOcc);
        }
      }
    }
    // alpha el. excitation
    // --- --- ==> -o- ---
    // --- ---     --- ---
    // -o- ---     -o- ---
    // -o- ---     -o- ---
    // -o- ---     -o- ---
    // -o- ---     -o- ---
    // -o- -o-     -o- -o-
    // -o- -o- ==> --- -o-
    // -o- -o-     -o- -o-
    //  a   b       a   b
    for(int i = 0; i < m_alpha; ++i) {
      int ir = reorder[i];
      for(int j = m_alpha; j < spac_norbs; ++j) {
        int jr = reorder[j];
        int exIrep = Symmetry::add(m_spin_orbs_symmetry[2*ir],m_spin_orbs_symmetry[2*jr])[0];
        if(Symmetry::add(refIrep,exIrep)[0] == irep) {
          std::vector<int> exOcc = refOcc;
          std::swap(exOcc[2*ir],exOcc[2*jr]);
          exOccs.push_back(exOcc);
        }
      }
    }

    // No 1-el. configuration was found for the requested state symmetry, irep
    assert(exOccs.size() > 0);

    double confE = 1.0e8;
    for(int i = 0; i < exOccs.size(); ++i) {
      double confE_tmp = 0.0;
      for(int j = 0; j < m_norbs; ++j) {
        if(exOccs[i][j] != 0) confE_tmp += v_1[0](j, j);
      }
      if(confE_tmp < confE) {
        confE = confE_tmp;
        refOcc = exOccs[i];
      }
    }
  }

  // NOTE: hfOcc in original order
  std::vector<int> hfOcc(m_norbs);
  for(int i = 0; i < spac_norbs; ++i) {
    int ir = m_reorder[i];
    hfOcc[2*ir]   = refOcc[2*i];
    hfOcc[2*ir+1] = refOcc[2*i+1];
  }

  return hfOcc;
}
