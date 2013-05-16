#include <algorithm>
#include <string>
#include <Eigen/Dense>
#include "npdm_spin_adaptation.h"
#include "MatrixBLAS.h"

namespace SpinAdapted {

namespace Npdm {

const expression::grammar<std::vector<TensorOp> ,std::string::const_iterator> eg;

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
//      std::cout <<"printing operator with Sz = "<<0<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
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
//      std::cout <<"printing operator with Sz = "<<-op.Spin<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
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

//  std::cout << "matrix_row:\n";
//  for ( auto it = matrix_row.begin(); it != matrix_row.end(); ++it ) {
//    assert( it->first.size() == 4 );
//   std::cout << (it->first)[0] << (it->first)[1] << (it->first)[2] << (it->first)[3] << std::endl;
//  }
//  std::cout << "-----------------------------\n";

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
  std::cout << "set size = " << indices_set.size() << std::endl;

  // Map group of indices to single integer in arbitrary, but well-defined, order
  int k=0;
  for ( auto it = indices_set.begin(); it != indices_set.end(); ++it ) {
    map_to_int[*it] = k;
    k++; 
    std::cout << (*it)[0] << (*it)[1] << (*it)[2] << (*it)[3] << std::endl;
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
  // Trim op string to just get ordering of creations and destructions (cre=0, des=1)
  std::vector<int> cd_order;
  for ( auto it = s.begin(); it != s.end(); ++it ) {
    if (*it=='C') cd_order.push_back(0);
    if (*it=='D') cd_order.push_back(1);
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

//FIXME check that we don't commute two identical SO indices

  // Sort into CCC..DDD.. order, hence implicitly build permutation matrix
  std::sort( perm_mat_pair.begin(), perm_mat_pair.end() );

  // Back out permutation matrix as explicit matrix type
  Eigen::MatrixXi perm_mat(dim,dim);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      perm_mat(i,j) = perm_mat_pair[i].second[j];
    }
  }

  // Re-order so_indices
  apply_permutation( perm_mat, so_indices );

  // Parity is just the determinant of permutation matrix
  int parity = perm_mat.determinant();
  assert( abs(parity) == 1 );

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
    std::cout << "WARNING: something wasn't quite perfect in parsing the operator string!\n";
    s.erase( s.begin() );   
    s.erase( s.begin() );   
    success = parse(iter, end, eg, result) ;		 
  }
  assert(success);
  std::cout << "Setting up linear equations for spin-adaptation transformation...\n";
  std::cout << "Number of compounded tensor operators = " << result.size() << std::endl;

  // RHS indices (i.e. rows of A) corresponding to singlet expectations
  std::vector<int> singlet_rows;

  // Recover transformation matrix A and relevant spin-orbital indices from tensor operators
  parse_result_into_matrix( result, A, so_indices, singlet_rows ); 

  // Get permutation parity and commute so_indices into NPDM order (i.e. cre,cre..,des,des..)
  int parity = commute_so_indices_to_pdm_order( s, so_indices );  
//std::cout << "parity = " << parity << std::endl;

  // Set up RHS of linear equations (note we assume the ordering of the singlets is consistent with earlier)
  assert( b0.size() == singlet_rows.size() );
  b=0.0;
  for (int i=0; i < b0.size(); ++i) {
    assert( singlet_rows[i] <= b.Nrows() );
    b( singlet_rows[i] ) = parity*b0[i];
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}
