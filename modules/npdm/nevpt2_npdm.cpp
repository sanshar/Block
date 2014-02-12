/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "nevpt2_npdm.h"

namespace SpinAdapted {

// FIXME parallelization
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm::save_A16_matrix_text()
{
  if( mpigetrank() == 0)
  {
    int dim = A16_matrix.dim1();
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/A16_matrix_direct.", 0, 0,".txt");
    ofstream ofs(file);
    ofs << dim << endl;

    double trace = 0.0;
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        for(int k=0; k<dim; ++k)
          for(int l=0; l<dim; ++l)
            for(int m=0; m<dim; ++m)
              for(int n=0; n<dim; ++n) {
                if ( abs(A16_matrix(i,j,k,l,m,n)) > 1e-14 ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % A16_matrix(i,j,k,l,m,n);
                }
              }
    ofs.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

array_6d<double> Nevpt2_npdm::compute_EEE_matrix( Matrix& onepdm, array_4d<double>& twopdm, array_6d<double>& threepdm )
{
if( mpigetrank() == 0 ) {

  std::cout << "Building spatial <0|EEE|0>\n";

  int dim = threepdm.dim1(); 
  assert( onepdm.Nrows() == twopdm.dim1() );
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
              val += d_jk * d_lm * onepdm(i+1,n+1); // Note Matrix indices start at 1 not 0
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

array_8d<double> Nevpt2_npdm::compute_EEEE_matrix( Matrix& onepdm, array_4d<double>& twopdm, array_6d<double>& threepdm, array_8d<double>& fourpdm )
{
if( mpigetrank() == 0 ) {

  std::cout << "Building spatial <0|EEEE|0>\n";

  int dim = threepdm.dim1(); 
  assert( onepdm.Nrows() == twopdm.dim1() );
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
//                  // 1PDM terms
//                  val += d_jk * d_lm * d_np * onepdm(i+1,q+1); // Note Matrix indices start at 1 not 0
//                  // 2PDM terms 
//                  val += d_jk * d_lm * 2.0*twopdm(i,p,q,n); // Note factor of two difference between Block 2PDMs and other codes
//                  val += d_jk * d_np * 2.0*twopdm(i,m,q,l);
//                  val += d_jk * d_lp * 2.0*twopdm(i,m,n,q);
//                  val += d_jm * d_np * 2.0*twopdm(i,k,l,q);
//                  val += d_jm * d_lp * 2.0*twopdm(i,k,q,n);
//                  val += d_lm * d_np * 2.0*twopdm(i,k,q,j);
//                  val += d_lm * d_jp * 2.0*twopdm(i,k,n,q);
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
  const TwoElectronArray& twoInt = v_2;

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
  const TwoElectronArray& twoInt = v_2;

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm::update_A16_matrix( std::map< std::vector<int>, double >& threepdm, std::map< std::vector<int>, double >& fourpdm ) 
{
cout << "updating A16 matrix using npdm elements from this sweep position\n";

  const TwoElectronArray& twoInt = v_2;
  int dim = A16_matrix.dim1();

  // 2PDM and 1PDM terms are added on using full matrices at end of sweep
  int i,j,k,l,m,n,p,q;

  //--------------
  // 4PDM terms
  //--------------
  for (auto it = fourpdm.begin(); it != fourpdm.end(); ++it) {
    i = (it->first)[0];
    k = (it->first)[1];
    m = (it->first)[2];
    p = (it->first)[3];
    j = (it->first)[4];
    l = (it->first)[5];
    n = (it->first)[6];
    q = (it->first)[7];
    for ( int a=0; a<dim; ++a )
      A16_matrix(j,k,i,a,l,q) += twoInt(2*m,2*p,2*n,2*a) * it->second;
    for ( int c=0; c<dim; ++c )
      A16_matrix(j,k,i,p,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
    for ( int b=0; b<dim; ++b )
      A16_matrix(j,k,i,p,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
  }

  //--------------
  // 3PDM terms
  //--------------
  for (auto it = threepdm.begin(); it != threepdm.end(); ++it) {
    // delta_jk terms
    i = (it->first)[0]; m = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; n = (it->first)[4]; l = (it->first)[5];
    for ( int jk=0; jk<dim; ++jk ) {
      for ( int a=0; a<dim; ++a )
        A16_matrix(jk,jk,i,a,l,q) += twoInt(2*m,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        A16_matrix(jk,jk,i,p,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        A16_matrix(jk,jk,i,p,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
    }
    // delta_jm terms
    i = (it->first)[0]; k = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; l = (it->first)[4]; n = (it->first)[5];
    for ( int jm=0; jm<dim; ++jm ) {
      for ( int a=0; a<dim; ++a )
        A16_matrix(jm,k,i,a,l,q) += twoInt(2*jm,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        A16_matrix(jm,k,i,p,l,c) -= twoInt(2*jm,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        A16_matrix(jm,k,i,p,b,q) -= twoInt(2*jm,2*b,2*n,2*l) * it->second;
    }
    // delta_lm terms
    i = (it->first)[0]; k = (it->first)[1]; p = (it->first)[2]; q = (it->first)[3]; n = (it->first)[4]; j = (it->first)[5];
    for ( int lm=0; lm<dim; ++lm ) {
      for ( int a=0; a<dim; ++a )
        A16_matrix(j,k,i,a,lm,q) += twoInt(2*lm,2*p,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        A16_matrix(j,k,i,p,lm,c) -= twoInt(2*lm,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        A16_matrix(j,k,i,p,b,q) -= twoInt(2*lm,2*b,2*n,2*lm) * it->second;
    }
    // delta_np terms
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; q = (it->first)[3]; l = (it->first)[4]; j = (it->first)[5];
    for ( int np=0; np<dim; ++np ) {
      for ( int a=0; a<dim; ++a )
        A16_matrix(j,k,i,a,l,q) += twoInt(2*m,2*np,2*np,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        A16_matrix(j,k,i,np,l,c) -= twoInt(2*m,2*c,2*np,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        A16_matrix(j,k,i,np,b,q) -= twoInt(2*m,2*b,2*np,2*l) * it->second;
    }
    // delta_lp terms
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; n = (it->first)[3]; q = (it->first)[4]; j = (it->first)[5];
    for ( int lp=0; lp<dim; ++lp ) {
      for ( int a=0; a<dim; ++a )
        A16_matrix(j,k,i,a,lp,q) += twoInt(2*m,2*lp,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        A16_matrix(j,k,i,lp,lp,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        A16_matrix(j,k,i,lp,b,q) -= twoInt(2*m,2*b,2*n,2*lp) * it->second;
    }
    // delta_jp terms
    i = (it->first)[0]; k = (it->first)[1]; m = (it->first)[2]; n = (it->first)[3]; l = (it->first)[4]; q = (it->first)[5];
    for ( int jp=0; jp<dim; ++jp ) {
      for ( int a=0; a<dim; ++a )
        A16_matrix(jp,k,i,a,l,q) += twoInt(2*m,2*jp,2*n,2*a) * it->second;
      for ( int c=0; c<dim; ++c )
        A16_matrix(jp,k,i,jp,l,c) -= twoInt(2*m,2*c,2*n,2*q) * it->second;
      for ( int b=0; b<dim; ++b )
        A16_matrix(jp,k,i,jp,b,q) -= twoInt(2*m,2*b,2*n,2*l) * it->second;
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Nevpt2_npdm::finish_A16_matrix( Matrix& onepdm, array_2d<double> twopdm )
{
  int dim = twopdm.dim1(); 
  assert( onepdm.Nrows() == twopdm.dim1() );

}

//===========================================================================================================================================================

}
