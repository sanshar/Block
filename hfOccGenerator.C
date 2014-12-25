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
  int norbs_spac = m_norbs/2;

  // Sorting orbitals by 1-el. Hamiltonian, i.e. h(i, i)
  std::vector<int> reorder(norbs_spac,0);
  {
    std::multimap<double,int> orbMap;
    for(int i = 0; i < norbs_spac; ++i)
      orbMap.insert(std::make_pair(v_1[0](2*i, 2*i), i));

    mapIter it = orbMap.begin();
    for(int i = 0; i < norbs_spac; ++i, ++it)
      reorder[i] = it->second;
  }

  int refIrep = 0;

  std::vector<int> refOcc(m_norbs,0);
  for(int i = 0; i < m_alpha; ++i) {
    int ir = reorder[i];
    refOcc[2*ir] = 1;
    refIrep = groupTable(refIrep, m_spin_orbs_symmetry[2*ir]);
    if(i < m_beta) {
      refOcc[2*ir+1] = 1;
      refIrep = groupTable(refIrep, m_spin_orbs_symmetry[2*ir+1]);
    }
  }

  std::vector<int> hfOcc;

  if(refIrep == irep) {
    hfOcc = refOcc;
  }
  else {
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
      for(int j = m_beta; j < m_alpha; ++j) {
        int exIrep = groupTable(m_spin_orbs_symmetry[2*i+1],m_spin_orbs_symmetry[2*j+1]);
        if(groupTable(refIrep,exIrep) == irep) {
          std::vector<int> exOcc = refOcc;
          std::swap(exOcc[2*i+1],exOcc[2*j+1]);
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
      for(int j = m_alpha; j < norbs_spac; ++j) {
        int exIrep = groupTable(m_spin_orbs_symmetry[2*i],m_spin_orbs_symmetry[2*j]);
        if(groupTable(refIrep,exIrep) == irep) {
          std::vector<int> exOcc = refOcc;
          std::swap(exOcc[2*i],exOcc[2*j]);
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
        hfOcc = exOccs[i];
      }
    }
  }

  return hfOcc;
}
