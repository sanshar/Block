#include "npdm_tensorop_arithematic.h"
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

  // Map group of indices to single integer in arbitrary order
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
                               Matrix& matrix, std::vector< std::vector<int> >& so_indices, std::vector<int>& singlet_indices)
{

  // Note the indices in the matrix/vector notation start at 1 not 0
  std::vector< std::map< std::vector<int>, double > > matrix_rows;
  std::map< std::vector<int>, int > map_to_int;

  // Extract matrix rows from tensor operators
  singlet_indices.clear();
  for (int i=0; i < tensor_ops.size(); ++i) {
    bool is_singlet;
    matrix_rows.push_back( get_matrix_row(tensor_ops[i], is_singlet) );
    if (is_singlet) singlet_indices.push_back(i+1);
  }

  // Make a map of spin-orbital indices to single integers
  map_to_int = get_map_to_int( matrix_rows );

  // Format matrix_rows into an actual matrix
  int dim = matrix_rows.size();
  assert( dim == map_to_int.size() );
  assert( dim == matrix.Nrows() );
  assert( dim == matrix.Ncols() );
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

void npdm_get_spin_adapt_matrix(std::string& s, Matrix& A, std::vector< std::vector<int> >& so_indices, std::vector<int> & singlet_indices)
{
  // Parse string of form (((C1C2)(D3))(D4)) etc
  std::string::const_iterator iter = s.begin();        
  std::string::const_iterator end  = s.end();          
  std::vector<TensorOp> result;                                       
  bool success = parse(iter, end, eg, result) ;		 
  assert(success);
  std::cout << "Number of compounded tensor operators = " << result.size() << std::endl;

  // Recover transformation matrix A from tensor operators
  parse_result_into_matrix( result, A, so_indices, singlet_indices ); 

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}
