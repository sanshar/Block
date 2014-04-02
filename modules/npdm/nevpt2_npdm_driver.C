/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_npdm_driver.h"
#include "nevpt2_npdm_matrices.h"

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
  nevpt2_A16_matrix( Nevpt2_A16_matrix( sites ) ), 
  onepdm_driver( Npdm_driver(1, nevpt2_A16_matrix ) ),
  twopdm_driver( Npdm_driver(2, nevpt2_A16_matrix ) ),
  threepdm_driver( Npdm_driver(3, nevpt2_A16_matrix ) ),
  fourpdm_driver( Npdm_driver(4, nevpt2_A16_matrix ) )
#endif
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::save_data()
{
#ifdef DEBUG_NEVPT2NPDM
  // Build NPDMS first, construct A-matrices later
  onepdm_driver.save_data();
  twopdm_driver.save_data();
  threepdm_driver.save_data();
  fourpdm_driver.save_data();
  compute_matrices();
#else
  // Build A-matrices on the fly
  nevpt2_A16_matrix.save_npdms(0,0);
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

void Nevpt2_npdm_driver::compute_matrices()
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if ( world.size() > 1 ) abort();  // not parallelized
#endif
  // 1PDM
  char file[5000];
  array_2d<double> onepdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_onepdm.", 0, 0,".bin");
  std::ifstream ifs1(file, std::ios::binary);
  boost::archive::binary_iarchive load1(ifs1);
  load1 >> onepdm;
  ifs1.close();
  // 2PDM
  array_4d<double> twopdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_twopdm.", 0, 0,".bin");
  std::ifstream ifs2(file, std::ios::binary);
  boost::archive::binary_iarchive load2(ifs2);
  load2 >> twopdm;
  ifs2.close();
  // 3PDM
  array_6d<double> threepdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm.", 0, 0,".bin");
  std::ifstream ifs3(file, std::ios::binary);
  boost::archive::binary_iarchive load3(ifs3);
  load3 >> threepdm;
  ifs3.close();
  // 4PDM
  array_8d<double> fourpdm;
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/spatial_fourpdm.", 0, 0,".bin");
  std::ifstream ifs4(file, std::ios::binary);
  boost::archive::binary_iarchive load4(ifs4);
  load4 >> fourpdm;
  ifs4.close();

  // Compute and store NEVPT2 matrices
  Nevpt2_npdm nevpt2;
  array_8d<double> eeee_matrix = nevpt2.compute_EEEE_matrix( onepdm, twopdm, threepdm, fourpdm );
  array_6d<double> eee_matrix = nevpt2.compute_EEE_matrix( onepdm, twopdm, threepdm );
  nevpt2.compute_A16_matrix( eeee_matrix );
  nevpt2.compute_A22_matrix( eee_matrix, eeee_matrix );

}

//===========================================================================================================================================================

}
}

