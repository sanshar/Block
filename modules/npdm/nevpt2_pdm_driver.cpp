/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_pdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Nevpt2_pdm_driver::Nevpt2_pdm_driver( int sites )
{
  twopdm_container = boost::shared_ptr<Twopdm_container>( new Twopdm_container(sites) );
  threepdm_container = boost::shared_ptr<Threepdm_container>( new Threepdm_container(sites) );
  fourpdm_container = boost::shared_ptr<Fourpdm_container>( new Fourpdm_container(sites) );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_pdm_driver::save_npdms(const int& i, const int& j)
{
  twopdm_container->save_npdms(i,j);
  threepdm_container->save_npdms(i,j);
  fourpdm_container->save_npdms(i,j);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void update_A16_matrix( Twopdm_container& twopdm_container, Threepdm_container& Threepdm_container, Fourpdm_container& fourpdm_container ) 
{
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_pdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos) 
{

  cout << "Computing all relevent NPDM matrix elements for NEVPT2 at this sweep position\n";
  // Get NPDM elements from this sweep position
  

  twopdm_container->compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  threepdm_container->compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  fourpdm_container->compute_npdm_elements(wavefunctions, big, sweepPos, endPos);

  // Increment NEVPT2 matrices with information from this sweep position.
//  update_A16_matrix( *twopdm_container, *threepdm_container, *fourpdm_container );
//  update_A22_matrix( twopdm_container, threepdm_container, fourpdm_container );

}


//===========================================================================================================================================================

}

