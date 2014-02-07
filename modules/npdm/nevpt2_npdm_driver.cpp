/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm.h"
#include "nevpt2.h"
#include "nevpt2_npdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Nevpt2_npdm_driver::Nevpt2_npdm_driver( int sites ) :
  twopdm_container( Twopdm_container(sites) ), 
  twopdm_driver( Npdm_driver(2, twopdm_container) ),
  threepdm_container( Threepdm_container(sites) ), 
  threepdm_driver( Npdm_driver(3, threepdm_container) ),
  fourpdm_container( Fourpdm_container(sites) ), 
  fourpdm_driver( Npdm_driver(4, fourpdm_container) )
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::save_data()
{
  twopdm_driver.save_data();
  threepdm_driver.save_data();
  fourpdm_driver.save_data();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos) 
{
  cout << "Computing all relevent NPDM matrix elements for NEVPT2 at this sweep position\n";
  // Compute NPDM elements at this sweep position
  twopdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  threepdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  fourpdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);

  // Increment NEVPT2 matrices with information from this sweep position.
//  compute_A16_matrix( twopdm_container, threepdm_container, fourpdm_container );

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::compute_matrices()
{
  Matrix onepdm; 
  int i=0; int j=0;
  load_onepdm_spatial_binary(onepdm,i,j);
  array_4d<double>& twopdm = twopdm_container.get_spatial_twopdm();
  array_6d<double>& threepdm = threepdm_container.get_spatial_threepdm();
  array_8d<double>& fourpdm = fourpdm_container.get_spatial_fourpdm();

  array_8d<double> eeee_matrix = compute_EEEE_matrix( onepdm, twopdm, threepdm, fourpdm );
  compute_A16_matrix( fourpdm.dim1(), eeee_matrix );
}

//===========================================================================================================================================================

}

