#ifndef NPDM_PATTERNS_H
#define NPDM_PATTERNS_H

#include <set>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include "npdm.h"

namespace SpinAdapted {
namespace Npdm {

// Note this definition is not arbitrary as we sort against these integers
enum CD { CREATION=0, DESTRUCTION=1, BREAK=2 };

//===========================================================================================================================================================

class Npdm_patterns
{

  public:
    Npdm_patterns() { pdm_order_=NPDM_EMPTY; };
    Npdm_patterns( NpdmOrder pdm_order, int sweep_pos, int end_pos );

    int size() { return ldr_cd_types_.size(); };
    std::set< std::map< char, std::vector<CD> > >::const_iterator ldr_cd_begin() { return ldr_cd_types_.begin(); };
    std::set< std::map< char, std::vector<CD> > >::const_iterator ldr_cd_end() { return ldr_cd_types_.end(); };

    std::set< std::vector<CD> >::const_iterator lhs_cd_begin() { return lhs_cd_types_.begin(); };
    std::set< std::vector<CD> >::const_iterator dot_cd_begin() { return dot_cd_types_.begin(); };
    std::set< std::vector<CD> >::const_iterator rhs_cd_begin() { return rhs_cd_types_.begin(); };
    std::set< std::vector<CD> >::const_iterator lhs_cd_end() { return lhs_cd_types_.end(); };
    std::set< std::vector<CD> >::const_iterator dot_cd_end() { return dot_cd_types_.end(); };
    std::set< std::vector<CD> >::const_iterator rhs_cd_end() { return rhs_cd_types_.end(); };
    void print_cd_string( const std::vector<CD> & );

    // Screening of duplications
    bool screen_2pdm_strings( const std::vector<int>& indices, const std::string& CD );
    bool screen_3pdm_strings( const std::vector<int>& indices, const std::string& CD );
    bool screen_4pdm_strings( const std::vector<int>& indices, const std::string& CD );

    //FIXME
//  private:
    NpdmOrder pdm_order_;
    // Operator dimensions on LHS, RHS and Dot (add up to 2*order of PDM)
    std::set< std::tuple<int,int,int> > lhs_dot_rhs_types_;
    // Creation/destruction patterns for a given operator string to map to all irreducible permutations
    std::set< std::vector<CD> > cre_des_types_;
    // Aggregation of lhs_dot_rhs_types and cre_des_types (char => LHS, RHS or Dot)
    std::set< std::vector<CD> > lhs_cd_types_;
    std::set< std::vector<CD> > dot_cd_types_;
    std::set< std::vector<CD> > rhs_cd_types_;
    std::set< std::map< char, std::vector<CD> > > ldr_cd_types_;

    // Build all operator partitionings for the NPDM
    void build_lhs_dot_rhs_types( int sweep_pos, int end_pos );
    // Build all creation-destruction patterns
    void build_cre_des_types();
    void add_operator( int cre_ops, int des_ops, std::vector<CD> cd_type1 );
    // Combine ldr and cd patterns together
    void build_ldr_cd_types( int sweep_pos, int end_pos );

    // Test if valid creation-destruction string on dot assuming fermions
    bool is_valid_dot_type( std::vector<CD> );
    // Test if unique creation-destruction string over all blocks
    bool is_valid_ldr_type( std::map< char, std::vector<CD> > & cd_pattern );
    bool is_rhs_gte_lhs( std::vector<int> &, std::vector<int> & );

    // Print functions
    void print_int_string( const std::vector<int> & vec );
//FIXME
// make non-copyable!

};

//===========================================================================================================================================================

}
}

#endif

