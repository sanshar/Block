#include <set>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "npdm_spin_adaptation.h"
#include "MatrixBLAS.h"
#include "npdm_spin_transformation.h"
#include "pario.h"

namespace SpinAdapted {
namespace Npdm {

//===========================================================================================================================================================

std::map< std::vector<int>, double > Npdm_spin_adaptation::get_matrix_row ( const TensorOp& op, bool & is_singlet )
{
  std::map< std::vector<int>, double > matrix_row;
  double coeff;
  std::vector<int> indices;

  is_singlet = false; 
  if (op.Spin == 0) is_singlet = true;

  if (op.Spin%2 == 0) {
    for (int ilz = 0; ilz<op.rows; ilz++){
      //int lz = op.lz[ilz];
      //pout <<"printing operator with Sz = "<<0<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
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
      //pout <<"printing operator with Sz = "<<-op.Spin<<" and row = "<< ilz<<" and total spin "<<op.Spin<<" and irrep "<<op.irrep<<endl;;
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

  return matrix_row;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::map< std::vector<int>, int > Npdm_spin_adaptation::get_map_to_int( std::vector< std::map< std::vector<int>, double > > & matrix_rows )
{
  std::map< std::vector<int>, int > map_to_int;
  std::set< std::vector<int> > indices_set;

  // Make set of unique spin-orbital index tuples
  for ( int i=0; i < matrix_rows.size(); ++i) {
    for ( auto it = matrix_rows[i].begin(); it != matrix_rows[i].end(); ++it ) {
      indices_set.insert(it->first);
    }
  }

  // Map group of indices to single integer in arbitrary, but well-defined, order
  int k=0;
  for ( auto it = indices_set.begin(); it != indices_set.end(); ++it ) {
    map_to_int[*it] = k;
    k++; 
    //pout << "spin indices = ";
    //for (auto itt = it->begin(); itt != it->end(); ++itt) { pout << *itt << " "; } pout << std::endl;
  }

  return map_to_int;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
   
void Npdm_spin_adaptation::parse_result_into_matrix( const std::vector<TensorOp>& tensor_ops, 
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

void Npdm_spin_adaptation::apply_permutation( Matrix & perm_mat, std::vector< std::vector<int> >& so_indices)
{
  int cols = so_indices.size();
  int rows = so_indices[0].size();
  assert( perm_mat.Nrows() == rows );
  assert( perm_mat.Ncols() == rows );

  // Convert so_indices to explicit matrix type
  Matrix so_mat(rows,cols);
  for (int j=0; j<cols; j++) {
    for (int i=0; i<rows; i++) {
      so_mat(i+1,j+1) = so_indices[j][i];
    }
  }

  // Transform by matrix multiplication
  Matrix new_so_mat(rows,cols);
  new_so_mat = perm_mat * so_mat;

  // Over-write so_indices with new ordering
  for (int j=0; j<cols; j++) {
    for (int i=0; i<rows; i++) {
      so_indices[j][i] = new_so_mat(i+1,j+1);
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

int Npdm_spin_adaptation::commute_so_indices_to_pdm_order( const std::string& s, std::vector< std::vector<int> >& so_indices)
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

  Matrix perm_mat(dim, dim);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      perm_mat(i+1,j+1) = perm_mat_pair[i].second[j];
    }
  }

  apply_permutation( perm_mat, so_indices);

  
  // Parity is just the determinant of permutation matrix
  int parity = perm_mat.Determinant() > 0 ? 1 : -1;
  assert( abs(parity) == 1 );

  return parity;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_spin_adaptation::store_new_A_mat( const int so_dim, const int order, const std::string& cd_string )
{
  // Build an operator string with consecutive spatial indices
  std::vector<int> indices(order);
  for (int i=0; i<order; ++i) {
    indices[i] = i+1;
  }

  // Combine indices and build_pattern into one string
  std::string op_string;
  for (auto it = cd_string.begin(); it != cd_string.end(); ++it) {
    op_string.push_back(*it);
    if ( (*it == 'C') || (*it == 'D') ) {
      std::string index = boost::lexical_cast<std::string>( indices.at(0) );
      op_string.append(index);
      indices.erase( indices.begin() );
    }
  }

  // Parse op_string to build everything we want to store
  Matrix A(so_dim,so_dim);
  std::vector<int> singlet_rows; 
  std::vector<std::vector<int> > so_indices(so_dim);
  build_new_A_mat( op_string, A, singlet_rows, so_indices );

  // Store for later reuse
  stored_A_mats_[ cd_string ] = A;
  stored_singlet_rows_[ cd_string ] = singlet_rows;
  stored_so_indices_[ cd_string ] = so_indices;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_spin_adaptation::build_new_A_mat( const std::string& op_string, Matrix& A, std::vector<int>& singlet_rows, 
                                            std::vector<std::vector<int> >& so_indices)
{
  const expression::grammar< std::vector<TensorOp>, std::string::const_iterator > eg;

  // Parse string of form (((C1C2)(D3))(D4)) etc
  std::string op = op_string;
  std::string::const_iterator iter = op.begin();        
  std::string::const_iterator end  = op.end();          
  std::vector<TensorOp> result;                                       
  bool success = parse(iter, end, eg, result) ;		 
  // Edge case: ()((CxCx)... Arises when LHS and Dot blocks are empty.
  if ( not success ) {
    op.erase( op.begin() );   
    op.erase( op.begin() );   
    success = parse(iter, end, eg, result) ;		 
  }
  assert( success );
  assert( result.size() == so_indices.size() );

  // Get transformation matrix A and relevant spin-orbital indices from tensor operators
  parse_result_into_matrix( result, A, so_indices, singlet_rows ); 

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_spin_adaptation::get_so_indices( std::string& cd_string, const std::vector<int>& indices, std::vector< std::vector<int> >& so_indices )
{
  // Assume all the stored spin-indices were generated from spatial indices of the form 1,2,3,4,5,6...
  so_indices.resize(0);
  std::vector< std::vector<int> > stored_so_indices = stored_so_indices_.at( cd_string );
  // Loop over stored spin-orbital elements
  for (auto vec = stored_so_indices.begin(); vec != stored_so_indices.end(); ++vec) {
    std::vector<int> new_element;
    // Loop over stored spin-indices for one element
    for (auto it = vec->begin(); it != vec->end(); ++it) {
      int i = (*it)/2;
      int j = (*it)%2;
      // Generate new spin indices from spatial indices
      new_element.push_back( 2*indices.at(i-1) + j ); 
    }
    so_indices.push_back( new_element );
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_spin_adaptation::npdm_set_up_linear_equations( const int dim, const std::string& op_string, const std::vector<int>& indices,
                                                         const std::vector<double>& b0, Matrix& A, 
                                                         ColumnVector& b, std::vector< std::vector<int> >& so_indices )
{
  // Separate cre-des pattern and numerical indices
  std::string cd_string;
  for (auto it = op_string.begin(); it != op_string.end(); ++it) {
    if ( (*it == '(') || (*it == ')') || (*it == 'C') || (*it == 'D') ) cd_string.push_back(*it); 
  }
  
  // RHS indices (i.e. rows of A) corresponding to singlet expectations
  std::vector<int> singlet_rows;

  if ( stored_A_mats_.find(cd_string) != stored_A_mats_.end() ) {
    A = stored_A_mats_.at( cd_string );
    singlet_rows = stored_singlet_rows_.at( cd_string );
    get_so_indices( cd_string, indices, so_indices );
  }
  else {
    build_new_A_mat( op_string, A, singlet_rows, so_indices );
    store_new_A_mat( dim, indices.size(), cd_string );
  }

  // Get permutation parity and commute so_indices into NPDM order (i.e. cre,cre..,des,des..)
  int parity = commute_so_indices_to_pdm_order( op_string, so_indices );  

  // Set up RHS of linear equations (note we assume the ordering of the singlets is consistent with earlier)
  assert( b0.size() == singlet_rows.size() );
  b=0.0;
  for (int i=0; i < b0.size(); ++i) {
    assert( singlet_rows[i] <= b.Nrows() );
    b( singlet_rows[i] ) = parity*b0[i];
  }

}
  
//===========================================================================================================================================================

}
}
