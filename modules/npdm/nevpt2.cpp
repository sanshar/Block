/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2.h"

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void compute_EEEE_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& threepdm_driver, Fourpdm_driver& fourpdm_driver ) {

if( mpigetrank() == 0 ) {

std::cout << "Building spatial <0|EEEE|0>\n";

  int dim = threepdm_driver.spatial_threepdm.dim1(); 
  array_8d<double> eeee_matrix(dim,dim,dim,dim,dim,dim,dim,dim);
  eeee_matrix.Clear();

  // Get 1PDM, 2PDM, 3PDM, 4PDM
  Matrix onepdm; 
  int i=0; int j=0;
  load_onepdm_spatial_binary(onepdm,i,j);
  array_4d<double>& twopdm = twopdm_driver.spatial_twopdm;
  array_6d<double>& threepdm = threepdm_driver.spatial_threepdm;
  array_8d<double>& fourpdm = fourpdm_driver.spatial_fourpdm;
  assert( onepdm.Nrows() == twopdm.dim1() );

  // Output text file
  double factor = 1.0;
  char file[5000];
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/EEEE_matrix.", i, j,".txt");
  ofstream ofs(file);
  ofs << fourpdm.dim1() << endl;

  for(int i=0; i<fourpdm.dim1(); ++i) {
    for(int j=0; j<fourpdm.dim2(); ++j) {

      for(int k=0; k<fourpdm.dim3(); ++k) {
        int d_jk = ( j == k );
        for(int l=0; l<fourpdm.dim4(); ++l) {
          for(int m=0; m<fourpdm.dim5(); ++m) {
            int d_lm = ( l == m );
            int d_jm = ( j == m );
            for(int n=0; n<fourpdm.dim6(); ++n) {
              for(int p=0; p<fourpdm.dim7(); ++p) {
                int d_np = ( n == p );
                int d_lp = ( l == p );
                int d_jp = ( j == p );
                for(int q=0; q<fourpdm.dim8(); ++q) {

                  // 1PDM term
                  double val = d_jk * d_lm * d_np * onepdm(i+1,q+1); // Note Matrix indices start at 1 not 0
                  // 2PDM terms 
                  val += d_jk * d_lm * 2.0*twopdm(i,p,q,n); // Note factor of two difference between Block 2PDMs and other codes
                  val += d_jk * d_np * 2.0*twopdm(i,m,q,l);
                  val += d_jk * d_lp * 2.0*twopdm(i,m,n,q);
                  val += d_jm * d_np * 2.0*twopdm(i,k,l,q);
                  val += d_jm * d_lp * 2.0*twopdm(i,k,q,n);
                  val += d_lm * d_np * 2.0*twopdm(i,k,q,j);
                  val += d_lm * d_jp * 2.0*twopdm(i,k,n,q);
                  // 3PDM terms
                  val += d_jk * threepdm(i,m,p,q,n,l);
                  val += d_jm * threepdm(i,k,p,q,l,n);
                  val += d_lm * threepdm(i,k,p,q,n,j);
                  val += d_np * threepdm(i,k,m,q,l,j);
                  val += d_lp * threepdm(i,k,m,n,q,j);
                  val += d_jp * threepdm(i,k,m,n,l,q);
                  // 4PDM term
                  val += fourpdm(i,k,m,p,j,l,n,q);

                  if ( abs(val) > 1e-14 ) {
                    ofs << boost::format("%d %d %d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % p % q % val;
                    eeee_matrix(i,j,k,l,m,n,p,q) = val;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  ofs.close();

  // Save binary matrix
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/EEEE_matrix.", i, j,".bin");
  std::ofstream ofs2(file, std::ios::binary);
  boost::archive::binary_oarchive save(ofs2);
  save << eeee_matrix;
  ofs2.close();

  std::cout << "done spatial <0|EEEE|0>\n";

}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void compute_A16_matrix( Twopdm_driver& twopdm_driver, Threepdm_driver& threepdm_driver, Fourpdm_driver& fourpdm_driver ) {

if( mpigetrank() == 0 ) {

std::cout << "Building NEVPT2 A16 matrix\n";

  const TwoElectronArray& twoInt = v_2;

  int dim = threepdm_driver.spatial_threepdm.dim1(); 

  // Get EEEE matrix from disk
  array_8d<double> eeee(dim,dim,dim,dim,dim,dim,dim,dim);
  char file[5000];
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(), "/EEEE_matrix.", 0, 0,".bin");
  std::ifstream ifs(file, std::ios::binary);
  boost::archive::binary_iarchive load(ifs);
  load >> eeee;
  ifs.close();

  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A16_matrix.", 0, 0,".txt");
  ofstream ofs(file);
  ofs << dim << endl;

  for(int ap=0; ap<dim; ++ap) {
    for(int bp=0; bp<dim; ++bp) {
      for(int cp=0; cp<dim; ++cp) {
        for(int a=0; a<dim; ++a) {
          for(int b=0; b<dim; ++b) {
            for(int c=0; c<dim; ++c) {
          
              double val = 0.0;
              for(int d=0; d<dim; ++d) {
                for(int e=0; e<dim; ++e) {
                  for(int f=0; f<dim; ++f) {
                    // Factor of 2 on indices to recover spatial two-electron integrals
                    val += twoInt(2*d,2*e,2*f,2*a) * eeee(cp,ap,bp,b,d,f,e,c);
                    val -= twoInt(2*d,2*c,2*f,2*e) * eeee(cp,ap,bp,b,d,f,a,e);
                    val -= twoInt(2*d,2*b,2*f,2*e) * eeee(cp,ap,bp,e,d,f,a,c);
                  }
                }
              }
              if ( abs(val) > 1e-14 ) ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % ap % bp % cp % a % b % c % val;
            }
          }
        }
      }
    }
  }
  ofs.close();
  std::cout << "done NEVPT2 A16 matrix\n";

}
}

//===========================================================================================================================================================

}
}

