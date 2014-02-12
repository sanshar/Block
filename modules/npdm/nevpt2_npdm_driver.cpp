/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "npdm.h"
#include "nevpt2_npdm.h"
#include "nevpt2_npdm_driver.h"

namespace SpinAdapted{

//===========================================================================================================================================================

Nevpt2_npdm_driver::Nevpt2_npdm_driver( int sites ) :
//  twopdm_container( Twopdm_container(sites) ), 
//  twopdm_driver( Npdm_driver(2, twopdm_container) ),
//  threepdm_container( Threepdm_container(sites) ), 
//  threepdm_driver( Npdm_driver(3, threepdm_container) ),
//  fourpdm_container( Fourpdm_container(sites) ), 
//  fourpdm_driver( Npdm_driver(4, fourpdm_container) )

  nevpt2_3pdm_container( Nevpt2_3pdm_container( a16_matrix) ), 
  threepdm_driver( Npdm_driver(3, nevpt2_3pdm_container) ),
  nevpt2_4pdm_container( Nevpt2_4pdm_container( a16_matrix ) ), 
  fourpdm_driver( Npdm_driver(4, nevpt2_4pdm_container) )
{ 
  a16_matrix.resize( sites,sites,sites,sites,sites,sites );
  a16_matrix.Clear();

  nevpt2.A16_matrix.resize( sites,sites,sites,sites,sites,sites );
  nevpt2.A16_matrix.Clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::save_data()
{
//  twopdm_driver.save_data();
//  threepdm_driver.save_data();
//  fourpdm_driver.save_data();

  // Finalize NEVPT2 matrices
//  Matrix onepdm; 
//  int i=0; int j=0;
//  load_onepdm_spatial_binary(onepdm,i,j);
//  nevpt2.add_A16_1pdm_part( onepdm );


//compute_matrices();
//  nevpt2.save_A16_matrix_text();

 if( mpigetrank() == 0)
  {
    int dim = a16_matrix.dim1();
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A16_matrix.", 0, 0,".txt");
    ofstream ofs(file);
    ofs << dim << endl;

    double trace = 0.0;
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        for(int k=0; k<dim; ++k)
          for(int l=0; l<dim; ++l)
            for(int m=0; m<dim; ++m)
              for(int n=0; n<dim; ++n) {
                if ( abs(a16_matrix(i,j,k,l,m,n)) > 1e-14 ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % a16_matrix(i,j,k,l,m,n);
                }
              }
    ofs.close();
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::compute_npdm_elements(std::vector<Wavefunction> & wavefunctions, const SpinBlock & big, int sweepPos, int endPos) 
{
  cout << "Computing all relevent NPDM matrix elements for NEVPT2 at this sweep position\n";
  // Compute NPDM elements at this sweep position
//  twopdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  threepdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);
  fourpdm_driver.compute_npdm_elements(wavefunctions, big, sweepPos, endPos);

  // Increment NEVPT2 matrices with information from this sweep position.
//  std::map< std::vector<int>, double > threepdm = threepdm_container.get_sparse_spatial_pdm();
//  std::map< std::vector<int>, double > fourpdm = fourpdm_container.get_sparse_spatial_pdm();
//  nevpt2.update_A16_matrix( threepdm, fourpdm );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm_driver::compute_matrices()
{
  Matrix onepdm; 
  int i=0; int j=0;
  load_onepdm_spatial_binary(onepdm,i,j);
//  array_4d<double>& twopdm = twopdm_container.get_spatial_twopdm();
//  array_6d<double>& threepdm = threepdm_container.get_spatial_threepdm();
//  array_8d<double>& fourpdm = fourpdm_container.get_spatial_fourpdm();

//  array_6d<double> eee_matrix = nevpt2.compute_EEE_matrix( onepdm, twopdm, threepdm );
//  array_8d<double> eeee_matrix = nevpt2.compute_EEEE_matrix( onepdm, twopdm, threepdm, fourpdm );

//  nevpt2.compute_A16_matrix( eeee_matrix );
//  nevpt2.compute_A22_matrix( eee_matrix, eeee_matrix );
}

//===========================================================================================================================================================

}

