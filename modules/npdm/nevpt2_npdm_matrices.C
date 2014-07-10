/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <boost/format.hpp>
#include "nevpt2_npdm_matrices.h"

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

array_6d<double> Nevpt2_npdm::compute_EEE_matrix( array_2d<double>& onepdm, array_4d<double>& twopdm, array_6d<double>& threepdm )
{
if( mpigetrank() == 0 ) {

  std::cout << "Building spatial <0|EEE|0>\n";

  int dim = threepdm.dim1(); 
  assert( onepdm.dim1() == twopdm.dim1() );
  array_6d<double> eee_matrix(dim,dim,dim,dim,dim,dim);
  eee_matrix.Clear();

  // Output text file
  double factor = 1.0;
  char file[5000];
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/EEE_matrix.", 0, 0,".txt");
  ofstream ofs(file);
  ofs << dim << endl;

  for(int i=0; i<dim; ++i) {
    for(int j=0; j<dim; ++j) {

      for(int k=0; k<dim; ++k) {
        int d_jk = ( j == k );
        for(int l=0; l<dim; ++l) {
          for(int m=0; m<dim; ++m) {
            int d_lm = ( l == m );
            int d_jm = ( j == m );
            for(int n=0; n<dim; ++n) {

              double val = 0.0;
              // 1PDM terms
              val += d_jk * d_lm * onepdm(i,n);
              // 2PDM terms 
              val += d_jk * 2.0*twopdm(i,m,n,l); // Note factor of two difference between Block 2PDMs and other codes
              val += d_lm * 2.0*twopdm(i,k,n,j);
              val += d_jm * 2.0*twopdm(i,k,l,n);
              // 3PDM terms
              val += threepdm(i,k,m,n,l,j);

              if ( abs(val) > 1e-14 ) {
                ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % val;
                eee_matrix(i,j,k,l,m,n) = val;
              }
            }
          }
        }
      }
    }
  }
  ofs.close();

  // Save binary matrix
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/EEE_matrix.", 0, 0,".bin");
  std::ofstream ofs2(file, std::ios::binary);
  boost::archive::binary_oarchive save(ofs2);
  save << eee_matrix;
  ofs2.close();

  return eee_matrix;

}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

array_8d<double> Nevpt2_npdm::compute_EEEE_matrix(array_2d<double>& onepdm, array_4d<double>& twopdm, array_6d<double>& threepdm, array_8d<double>& fourpdm)
{
if( mpigetrank() == 0 ) {

  std::cout << "Building spatial <0|EEEE|0>\n";

  int dim = threepdm.dim1(); 
  assert( onepdm.dim1() == twopdm.dim1() );
  array_8d<double> eeee_matrix(dim,dim,dim,dim,dim,dim,dim,dim);
  eeee_matrix.Clear();

  // Output text file
  double factor = 1.0;
  char file[5000];
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/EEEE_matrix.", 0, 0,".txt");
  ofstream ofs(file);
  ofs << dim << endl;

  for(int i=0; i<dim; ++i) {
    for(int j=0; j<dim; ++j) {

      for(int k=0; k<dim; ++k) {
        int d_jk = ( j == k );
        for(int l=0; l<dim; ++l) {
          for(int m=0; m<dim; ++m) {
            int d_lm = ( l == m );
            int d_jm = ( j == m );
            for(int n=0; n<dim; ++n) {
              for(int p=0; p<dim; ++p) {
                int d_np = ( n == p );
                int d_lp = ( l == p );
                int d_jp = ( j == p );
                for(int q=0; q<dim; ++q) {

                  double val = 0.0;
                  // 1PDM terms
                  val += d_jk * d_lm * d_np * onepdm(i,q);
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
                  // 4PDM terms
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
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/EEEE_matrix.", 0, 0,".bin");
  std::ofstream ofs2(file, std::ios::binary);
  boost::archive::binary_oarchive save(ofs2);
  save << eeee_matrix;
  ofs2.close();

  return eeee_matrix;

}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm::compute_A16_matrix( array_8d<double>& eeee ) 
{ 
if( mpigetrank() == 0 ) {

  std::cout << "Building NEVPT2 A16 matrix\n";

  int dim = eeee.dim1();
  const TwoElectronArray& twoInt = v_2[0];

  // Output text file
  char file[5000];
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

}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm::compute_A22_matrix( array_6d<double>& eee, array_8d<double>& eeee ) 
{ 
if( mpigetrank() == 0 ) {

  std::cout << "Building NEVPT2 A22 matrix\n";

  int dim = eee.dim1();
  const TwoElectronArray& twoInt = v_2[0];

  // Output text file
  char file[5000];
  sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A22_matrix.", 0, 0,".txt");
  ofstream ofs(file);
  ofs << dim << endl;

  for(int ap=0; ap<dim; ++ap) {
    for(int bp=0; bp<dim; ++bp) {
      for(int cp=0; cp<dim; ++cp) {
        for(int a=0; a<dim; ++a) {
          for(int b=0; b<dim; ++b) {
            int d_bpb = ( bp == b );
            for(int c=0; c<dim; ++c) {
          
              double val = 0.0;
              for(int d=0; d<dim; ++d) {
                for(int e=0; e<dim; ++e) {
                  int d_bpe = ( bp == e );
                  for(int f=0; f<dim; ++f) {
                    // Factor of 2 on indices to recover spatial two-electron integrals
                    val += twoInt(2*d,2*e,2*f,2*a) * ( 2.0 * d_bpb * eee(cp,ap,d,f,e,c) - eeee(cp,ap,b,bp,d,f,e,c) );
                    val -= twoInt(2*d,2*c,2*f,2*e) * ( 2.0 * d_bpb * eee(cp,ap,d,f,a,e) - eeee(cp,ap,b,bp,d,f,a,e) );
                    val += twoInt(2*d,2*e,2*f,2*b) * ( 2.0 * d_bpe * eee(cp,ap,d,f,a,c) - eeee(cp,ap,e,bp,d,f,a,c) );
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

}
}

//===========================================================================================================================================================

}
}

