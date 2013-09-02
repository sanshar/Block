#include <algorithm>
#include <string>
#include <set>
#include <Eigen/Dense>
#include "pario.h"
#include "npdm_spin_adaptation.h"
#include "MatrixBLAS.h"

namespace SpinAdapted {
namespace Npdm {

const expression::grammar< std::vector<TensorOp>, std::string::const_iterator > eg;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::vector<int>, double > get_matrix_row ( const TensorOp& op, bool & is_singlet )
{
  std::map< std::vector<int>, double > matrix_row;
  double coeff;
  std::vector<int> indices;

  is_singlet = false; 
  if (op.Spin == 0) is_singlet = true;

  if (op.Spin%2 == 0) {
    for (int ilz = 0; ilz<op.rows; ilz++){
      //int lz = op.lz[ilz];
//      pout <<"printing operator with Sz = "<<0<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
      for (int i=0; i<op.Szops[op.Spin/2].size(); i++) {
        if (op.Szops[ilz*(op.Spin+1)+op.Spin/2][i] != 0.0) {
          coeff = op.Szops[ilz*(op.Spin+1)+op.Spin/2][i];
          indices.clear();
          for (int j=0; j<op.opindices[i].size(); j++) {
            indices.push_back(op.opindices[i][j]);
          }
          // Add non-zero coeff to matrix row
          matrix_row[indices] = coeff;
        }
      }
    }
  }
  else {
    for (int ilz = 0; ilz<op.rows; ilz++){
      //int lz = op.lz[ilz];
//      pout <<"printing operator with Sz = "<<-op.Spin<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
      for (int i=0; i<op.Szops[op.Spin].size(); i++) {
        if (op.Szops[ilz*(op.Spin+1)+op.Spin][i] != 0.0) {
          coeff = op.Szops[ilz*(op.Spin+1)+op.Spin][i];
          indices.clear();
          for (int j=0; j<op.opindices[i].size(); j++) {
            indices.push_back(op.opindices[i][j]);
          }
          // Add non-zero coeff to matrix row
          matrix_row[indices] = coeff;
        }
      }
    }
  }

//  pout << "matrix_row:\n";
//  for ( auto it = matrix_row.begin(); it != matrix_row.end(); ++it ) {
//    assert( it->first.size() == 4 );
//   pout << (it->first)[0] << (it->first)[1] << (it->first)[2] << (it->first)[3] << std::endl;
//  }
//  pout << "-----------------------------\n";

  return matrix_row;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::vector<int>, int > get_map_to_int( std::vector< std::map< std::vector<int>, double > > & matrix_rows )
{
  std::map< std::vector<int>, int > map_to_int;
  std::set< std::vector<int> > indices_set;

  // Make set of unique spin-orbital index tuples
  for ( int i=0; i < matrix_rows.size(); ++i) {
    for ( auto it = matrix_rows[i].begin(); it != matrix_rows[i].end(); ++it ) {
      indices_set.insert(it->first);
    }
  }
//  pout << "set size = " << indices_set.size() << std::endl;

  // Map group of indices to single integer in arbitrary, but well-defined, order
  int k=0;
  for ( auto it = indices_set.begin(); it != indices_set.end(); ++it ) {
    map_to_int[*it] = k;
    k++; 
//    pout << "spin indices = ";
//    for (auto itt = it->begin(); itt != it->end(); ++itt) { pout << *itt << " "; } pout << std::endl;
  }

  return map_to_int;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
   
void parse_result_into_matrix( const std::vector<TensorOp>& tensor_ops, 
                               Matrix& matrix, std::vector< std::vector<int> >& so_indices, std::vector<int>& singlet_rows)
{

  // Note the indices in the matrix/vector notation start at 1 not 0
  std::vector< std::map< std::vector<int>, double > > matrix_rows;
  std::map< std::vector<int>, int > map_to_int;

  // Extract matrix rows from tensor operators
  singlet_rows.clear();
  for (int i=0; i < tensor_ops.size(); ++i) {
    bool is_singlet;
    matrix_rows.push_back( get_matrix_row(tensor_ops[i], is_singlet) );
    if (is_singlet) singlet_rows.push_back(i+1);
  }

  // Make a map of spin-orbital indices to single integers
  map_to_int = get_map_to_int( matrix_rows );

  // Format matrix_rows into an actual matrix
  int dim = matrix_rows.size();
  assert( dim == map_to_int.size() );
  assert( dim == matrix.Nrows() );
  assert( dim == matrix.Ncols() );
  // Don't forget to initialize!
  matrix=0.0;
  for (int i=0; i<dim; ++i) {
    for (auto it = matrix_rows[i].begin(); it != matrix_rows[i].end(); ++it) {
       int col = map_to_int.at( it->first );
       matrix(i+1,col+1) = it->second;
    }
  }

  // Format indices into an ordered vector consistent with matrix A
  so_indices.resize(dim);
  for (auto it = map_to_int.begin(); it != map_to_int.end(); ++it) {
    so_indices.at(it->second) = it->first;
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void apply_permutation( Eigen::MatrixXi & perm_mat, std::vector< std::vector<int> >& so_indices)
{
  int cols = so_indices.size();
  int rows = so_indices[0].size();
  assert( perm_mat.rows() == rows );
  assert( perm_mat.cols() == rows );

  // Convert so_indices to explicit matrix type
  Eigen::MatrixXi so_mat(rows,cols);
  for (int j=0; j<cols; j++) {
    for (int i=0; i<rows; i++) {
      so_mat(i,j) = so_indices[j][i];
    }
  }

  // Transform by matrix multiplication
  Eigen::MatrixXi new_so_mat(rows,cols);
  new_so_mat = perm_mat * so_mat;

  // Over-write so_indices with new ordering
  for (int j=0; j<cols; j++) {
    for (int i=0; i<rows; i++) {
      so_indices[j][i] = new_so_mat(i,j);
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

int commute_so_indices_to_pdm_order( std::string& s, std::vector< std::vector<int> >& so_indices)
{
  // Trim op string to just get ordering of creations and destructions (cre<0, des>0)
  std::vector<int> cd_order;
  int k=0;
  for ( auto it = s.begin(); it != s.end(); ++it ) {
    if (*it=='C') cd_order.push_back((k++));
    if (*it=='D') cd_order.push_back(1000+(k++));
  }
  assert ( cd_order.size() == so_indices.at(0).size() );
  int dim = cd_order.size();

  // Initialize permutation matrix for spin-orbital indices
  std::vector< std::pair< int, std::vector<int> > >  perm_mat_pair;
  int i = 0;
  for (auto cd = cd_order.begin(); cd != cd_order.end(); ++cd) {
    std::vector<int> row(dim,0);
    row[i] = 1;
    perm_mat_pair.push_back( std::make_pair(*cd, row) );
    i++;
  }

//FIXME check that we don't commute two C and D ops with identical indices 

  // Sort into CCC..DDD.. order, noting that original CC and DD sub-orders are preserved to avoid commuting identical indices;
  // hence implicitly build permutation matrix.
  // Note this doesn't preclude commutation of CD with same index, so must manually avoid that by not using such operators! (e.g. CDC 3-index)
  std::sort( perm_mat_pair.begin(), perm_mat_pair.end() );

  // Back out permutation matrix as explicit matrix type
  Eigen::MatrixXi perm_mat(dim,dim);
//pout << "perm_mat(i,j)\n";
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      perm_mat(i,j) = perm_mat_pair[i].second[j];
//pout << perm_mat(i,j) << " ";
    }
//pout << std::endl;
  }

  // Re-order so_indices
  apply_permutation( perm_mat, so_indices );

  // Parity is just the determinant of permutation matrix
  int parity = perm_mat.determinant();
  assert( abs(parity) == 1 );
//pout << "parity = " << parity << std::endl;

  return parity;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_set_up_linear_equations(std::string& s, std::vector<double>& b0, Matrix& A, ColumnVector& b, std::vector< std::vector<int> >& so_indices)
{
  // Parse string of form (((C1C2)(D3))(D4)) etc
  std::string::const_iterator iter = s.begin();        
  std::string::const_iterator end  = s.end();          
  std::vector<TensorOp> result;                                       
  bool success = parse(iter, end, eg, result) ;		 
  //FIXME
  // Edge case: ()((CxCx)... (Arises when LHS and Dot blocks are empty)
  if ( not success ) {
    cout << "WARNING: something wasn't quite perfect in parsing the operator string!\n";
    s.erase( s.begin() );   
    s.erase( s.begin() );   
    success = parse(iter, end, eg, result) ;		 
  }
  assert(success);
//  cout << "Setting up linear equations for spin-adaptation transformation...\n";
//  cout << "Number of compounded tensor operators = " << result.size() << std::endl;
  assert( result.size() == so_indices.size() );
  assert( result.size() == b.Nrows() );

  // RHS indices (i.e. rows of A) corresponding to singlet expectations
  std::vector<int> singlet_rows;

  // Recover transformation matrix A and relevant spin-orbital indices from tensor operators
  parse_result_into_matrix( result, A, so_indices, singlet_rows ); 

  // Get permutation parity and commute so_indices into NPDM order (i.e. cre,cre..,des,des..)
  int parity = commute_so_indices_to_pdm_order( s, so_indices );  

  // Set up RHS of linear equations (note we assume the ordering of the singlets is consistent with earlier)
  assert( b0.size() == singlet_rows.size() );
  b=0.0;
  for (int i=0; i < b0.size(); ++i) {
    assert( singlet_rows[i] <= b.Nrows() );
    b( singlet_rows[i] ) = parity*b0[i];
//pout << singlet_rows[i] << "\t\t" << parity*b0[i] << std::endl;
  }

}

/////////////////////////
//DEBUG
//
//  int dim = 20;
//  Eigen::MatrixXd Amat(dim,dim);
//  Eigen::VectorXd xvec(dim);
//  Eigen::VectorXd bvec(dim);
//
//  for (int i=0; i < dim; ++i) {
//    for (int j=0; j < dim; ++j) {
//      Amat(i,j) = A(i+1,j+1);
//    }
//  }
//////  pout << "Amat\n";
//////  pout << Amat << std::endl;
//
//xvec(0  )=   1.89959055084e-21 ;
//xvec(1  )=   2.75934014869e-22 ;
//xvec(2  )=   3.11219777415e-21 ;
//xvec(3  )=   1.58313329369e-05 ;
//xvec(4  )=   -1.58313329369e-05 ;
//xvec(5  )=   -6.2735506193e-22 ;
//xvec(6  )=   -1.58313329369e-05 ;
//xvec(7  )=   1.58313329369e-05 ;
//xvec(8  )=   2.14466425391e-22 ;
//xvec(9 )=    -1.18199994563e-21 ;
//xvec(10 )=   1.72761681302e-37 ;
//xvec(11 )=   3.11219777415e-21 ;
//xvec(12 )=   1.58313329369e-05 ;
//xvec(13 )=   -1.58313329369e-05 ;
//xvec(14 )=   -4.45652529234e-21 ;
//xvec(15 )=   -1.58313329369e-05 ;
//xvec(16 )=   1.58313329369e-05 ;
//xvec(17 )=   -1.99286234306e-21 ;
//xvec(18 )=   7.73535678409e-22 ;
//xvec(19 )=   -1.64090499098e-21 ;
//
//  bvec = Amat*xvec;
//  pout << "model b vector =\n";
//  pout << bvec << std::endl;
//
//}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}
