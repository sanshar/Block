/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_npdm_driver.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Nevpt2_npdm_driver::Nevpt2_npdm_driver( int sites ) :
#ifdef DEBUG_NEVPT2NPDM
  // Build NPDMS first, construct A-matrices later
  onepdm_driver( Onepdm_driver(sites) ),
  twopdm_driver( Twopdm_driver(sites) ),
  threepdm_driver( Threepdm_driver(sites) ),
  fourpdm_driver( Fourpdm_driver(sites) )
#else
  // Build A-matrices on the fly
  nevpt2_container( Nevpt2_container( sites ) ), 
  onepdm_driver( Npdm_driver(NPDM_ONEPDM, nevpt2_container ) ),
  twopdm_driver( Npdm_driver(NPDM_TWOPDM, nevpt2_container ) ),
  threepdm_driver( Npdm_driver(NPDM_THREEPDM, nevpt2_container ) ),
  fourpdm_driver( Npdm_driver(NPDM_FOURPDM, nevpt2_container ) )
#endif
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::save_data( const int i, const int j )
{
#ifdef DEBUG_NEVPT2NPDM
  // Build NPDMS first, construct A-matrices later
  onepdm_driver.save_data(i,j);
  twopdm_driver.save_data(i,j);
  threepdm_driver.save_data(i,j);
  fourpdm_driver.save_data(i,j);
  compute_matrices(i,j);
#else
  // Build A-matrices on the fly
  nevpt2_container.save_npdms(i,j);
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::clear()
{
    onepdm_driver.clear();
    twopdm_driver.clear();
    threepdm_driver.clear();
    fourpdm_driver.clear();
#ifndef DEBUG_NEVPT2NPDM
    nevpt2_container.clear();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos) 
{
  // Compute NPDM elements at this sweep position
  onepdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  twopdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  threepdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  fourpdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::compute_matrices( const int i, const int j )
{
  // Note not parallel
  if( mpigetrank() > 0 ) return;

  // Load 1PDM
  char file[5000];
  array_2d<double> onepdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_onepdm.", i, j,".bin");
  std::ifstream ifs1(file, std::ios::binary);
  boost::archive::binary_iarchive load1(ifs1);
  load1 >> onepdm;
  ifs1.close();
  // Load 2PDM
  array_4d<double> twopdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_twopdm.", i, j,".bin");
  std::ifstream ifs2(file, std::ios::binary);
  boost::archive::binary_iarchive load2(ifs2);
  load2 >> twopdm;
  ifs2.close();
  // Load 3PDM
  array_6d<double> threepdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm.", i, j,".bin");
  std::ifstream ifs3(file, std::ios::binary);
  boost::archive::binary_iarchive load3(ifs3);
  load3 >> threepdm;
  ifs3.close();
  // Load 4PDM
  array_8d<double> fourpdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_fourpdm.", i, j,".bin");
  std::ifstream ifs4(file, std::ios::binary);
  boost::archive::binary_iarchive load4(ifs4);
  load4 >> fourpdm;
  ifs4.close();

}

//===========================================================================================================================================================

}
}

