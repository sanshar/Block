#ifndef NPDM_SPIN_ADAPTATION_H
#define NPDM_SPIN_ADAPTATION_H

#include "MatrixBLAS.h"
#include "npdm_tensor_operator.h"

namespace SpinAdapted {
namespace Npdm {

//===========================================================================================================================================================

class Npdm_spin_adaptation {

  public:
    void npdm_set_up_linear_equations( const int dim, const std::string& op_string, const std::vector<int>& indices, const std::vector<double>& b0, 
                                       Matrix& A, ColumnVector& b, std::vector< std::vector<int> >& so_indices );

  private:
    std::map< std::string, Matrix > stored_A_mats_;
    std::map< std::string, std::vector<int> > stored_singlet_rows_;
    std::map< std::string, std::vector< std::vector<int> > > stored_so_indices_;

    void parse_result_into_matrix( const std::vector<TensorOp>& tensor_ops, 
                                   Matrix& matrix, std::vector< std::vector<int> >& so_indices, std::vector<int>& singlet_rows);
    
    void apply_permutation( Matrix & perm_mat, std::vector< std::vector<int> >& so_indices);
    int commute_so_indices_to_pdm_order( const std::string& s, std::vector< std::vector<int> >& so_indices);
    std::map< std::vector<int>, double > get_matrix_row ( const TensorOp& op, bool & is_singlet );
    std::map< std::vector<int>, int > get_map_to_int( std::vector< std::map< std::vector<int>, double > > & matrix_rows );
    void get_so_indices( std::string& cd_string, const std::vector<int>& indices, std::vector< std::vector<int> >& so_indices );

    void build_new_A_mat( const std::string& op_string, Matrix& A, std::vector<int>& singlet_rows, std::vector<std::vector<int> >& so_indices);
    void store_new_A_mat( const int so_dim, const int order, const std::string& cd_string );

};

//===========================================================================================================================================================

}
}

#endif
