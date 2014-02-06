/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_pdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Nevpt2_pdm_driver::Nevpt2_pdm_driver( int sites ) : Npdm_driver(0) {

  twopdm_driver = boost::shared_ptr<Twopdm_driver>( new Twopdm_driver(sites) );
  threepdm_driver = boost::shared_ptr<Threepdm_driver>( new Threepdm_driver(sites) );
  fourpdm_driver = boost::shared_ptr<Fourpdm_driver>( new Fourpdm_driver(sites) );

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_pdm_driver::save_npdms(const int& i, const int& j)
{
  twopdm_driver->save_npdms(i,j);
  threepdm_driver->save_npdms(i,j);
  fourpdm_driver->save_npdms(i,j);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void update_A16_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& Threepdm_driver, Fourpdm_driver& fourpdm_driver ) {
  assert(false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_pdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos) {

  cout << "Computing all relevent NPDM matrix elements for NEVPT2 at this sweep position\n";
  // Get NPDM elements from this sweep position
  twopdm_driver->compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  threepdm_driver->compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  fourpdm_driver->compute_npdm_elements(wavefunctions, big, sweepPos, endPos);

  // Increment NEVPT2 matrices with information from this sweep position.
//  update_A16_matrix( *twopdm_driver, *threepdm_driver, *fourpdm_driver );
//  update_A22_matrix( twopdm_driver, threepdm_driver, fourpdm_driver );

}


//===========================================================================================================================================================

}

