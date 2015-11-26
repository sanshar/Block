#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <string>
#include "global.h"
#include "pario.h"
#include "npdm_patterns.h"

namespace SpinAdapted {
namespace Npdm {

//===========================================================================================================================================================

Npdm_patterns::Npdm_patterns( NpdmOrder pdm_order, int sweep_pos, int end_pos )
: pdm_order_(pdm_order)
{
  build_lhs_dot_rhs_types( sweep_pos, end_pos );
  build_cre_des_types( );
  build_ldr_cd_types( sweep_pos, end_pos );
//  pout << "number of cre_des_pattern: " << cre_des_types_.size()<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// The pattern generator for the non-redundant NPDM elements leads to duplicates if indices are repeated, so some can be skipped explicitly.
// If the normal-ordered string is not of non-redundant form, then the original string produces duplicates when permutations are applied.
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_patterns::screen_2pdm_strings( const std::vector<int>& indices, const std::string& CD )
{
  if(!dmrginp.doimplicitTranspose()) return false;
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C' };
    if ( CD == foo ) return true;
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_patterns::screen_3pdm_strings( const std::vector<int>& indices, const std::string& CD )
{
  if(!dmrginp.doimplicitTranspose()) return false;
  if ( (indices[0] == indices[2]) && (indices[1] == indices[3]) ) {
    std::string foo = { 'C', 'C', 'D', 'D', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[1]) && (indices[2] == indices[3]) ) {
    std::string foo = { 'C', 'D', 'C', 'D', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'C', 'D' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'D', 'C', 'C' }; if ( CD == foo ) return true;
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_patterns::screen_4pdm_strings( const std::vector<int>& indices, const std::string& CD )
{
  if(!dmrginp.doimplicitTranspose()) return false;
  if ( (indices[0] == indices[1]) 
    && (indices[2] == indices[4])
    && (indices[3] == indices[5]) ) {
    std::string foo = { 'C', 'D', 'C', 'C', 'D', 'D', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[2]) 
    && (indices[1] == indices[3])
    && (indices[4] == indices[5]) ) {
    std::string foo = { 'C', 'C', 'D', 'D', 'C', 'D', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[1]) 
    && (indices[2] == indices[3])
    && (indices[4] == indices[5]) ) {
    std::string foo = { 'C', 'D', 'C', 'D', 'C', 'D', 'D', 'C' }; if ( CD == foo ) return true;
  }
  // New type
  if ( (indices[0] == indices[1]) 
    && (indices[2] == indices[3]) ) {
    std::string foo = { 'C', 'D', 'C', 'D', 'D', 'C', 'C', 'D' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[1]) 
    && (indices[2] == indices[3]) ) {
    std::string foo = { 'C', 'D', 'C', 'D', 'D', 'C', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[1]) 
    && (indices[2] == indices[3]) ) {
    std::string foo = { 'C', 'D', 'C', 'D', 'D', 'D', 'C', 'C' }; if ( CD == foo ) return true;
  }
  // Full 3 starting with CCDDD
  if ( (indices[0] == indices[2]) 
    && (indices[1] == indices[3]) ) {
    std::string foo = { 'C', 'C', 'D', 'D', 'D', 'C', 'C', 'D' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[2]) 
    && (indices[1] == indices[3]) ) {
    std::string foo = { 'C', 'C', 'D', 'D', 'D', 'C', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( (indices[0] == indices[2]) 
    && (indices[1] == indices[3]) ) {
    std::string foo = { 'C', 'C', 'D', 'D', 'D', 'D', 'C', 'C' }; if ( CD == foo ) return true;
  }
  // Full 10 starting with CDD
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'C', 'D', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'D', 'C', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'D', 'C', 'C', 'D' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'C', 'D', 'C', 'D' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'D', 'D', 'C', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'C', 'C', 'C', 'D', 'D' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'D', 'D', 'C', 'C', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'D', 'C', 'D', 'C', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'D', 'C', 'C', 'D', 'C' }; if ( CD == foo ) return true;
  }
  if ( indices[0] == indices[1] ) {
    std::string foo = { 'C', 'D', 'D', 'D', 'C', 'C', 'C', 'D' }; if ( CD == foo ) return true;
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Generates all required block partitionings for an NPDM

void Npdm_patterns::build_lhs_dot_rhs_types( int sweep_pos, int end_pos )
{

  //---------------
  // General case
  //---------------
  int lhs, rhs, dot, dotmax;
  int order_;
  if (pdm_order_ == NPDM_PAIRMATRIX) {
    order_ = 1;
  } else if (pdm_order_ >= NPDM_ONEPDM || pdm_order_ <= NPDM_FOURPDM) {
    order_ = pdm_order_ - NPDM_ONEPDM + 1;
  } else {
    abort();
  }
  for (lhs = order_; lhs >= 0; lhs--) {
    dotmax = 2*order_ - lhs;
    for (dot = dotmax; dot >= 1; dot--) {
      // Can have no more than 4 on the dot block
      if (dot > 4) continue;
      rhs = 2*order_ - dot - lhs;
      if ( rhs < order_ ) {
        //pout << lhs << " " << dot << " " << rhs << endl;
        lhs_dot_rhs_types_.insert( std::make_tuple(lhs,dot,rhs) );
      }
    }
  }

  //---------------
  // Edge cases 
  //---------------
  // 1PDM and pair matrix
  if (pdm_order_ == NPDM_ONEPDM || pdm_order_ == NPDM_PAIRMATRIX) {
    if ( sweep_pos == 0 ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(2,0,0) );
    }
    else if ( sweep_pos == end_pos ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(0,0,2) );
      lhs_dot_rhs_types_.insert( std::make_tuple(0,1,1) );
      lhs_dot_rhs_types_.insert( std::make_tuple(1,0,1) );
    }
  }
  // 2PDM
  else if (pdm_order_ == NPDM_TWOPDM) {
    if ( sweep_pos == 0 ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(4,0,0) );
      lhs_dot_rhs_types_.insert( std::make_tuple(3,1,0) );
      lhs_dot_rhs_types_.insert( std::make_tuple(3,0,1) );
    }
    else if ( sweep_pos == end_pos ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(0,2,2) );
      lhs_dot_rhs_types_.insert( std::make_tuple(2,0,2) );
      lhs_dot_rhs_types_.insert( std::make_tuple(1,1,2) );
      //
      lhs_dot_rhs_types_.insert( std::make_tuple(0,1,3) );
      lhs_dot_rhs_types_.insert( std::make_tuple(1,0,3) );
      lhs_dot_rhs_types_.insert( std::make_tuple(0,0,4) );
    }
  }
  // 3PDM
  else if (pdm_order_ == NPDM_THREEPDM) {
    if ( sweep_pos == 0 ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(4,2,0) );
      lhs_dot_rhs_types_.insert( std::make_tuple(4,0,2) );
      lhs_dot_rhs_types_.insert( std::make_tuple(4,1,1) );
    }
    else if ( sweep_pos == end_pos ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(0,3,3) );
      lhs_dot_rhs_types_.insert( std::make_tuple(3,0,3) );
      lhs_dot_rhs_types_.insert( std::make_tuple(1,2,3) );
      lhs_dot_rhs_types_.insert( std::make_tuple(2,1,3) );
      //
      lhs_dot_rhs_types_.insert( std::make_tuple(0,2,4) );
      lhs_dot_rhs_types_.insert( std::make_tuple(2,0,4) );
      lhs_dot_rhs_types_.insert( std::make_tuple(1,1,4) );
    }
  }
  // 4PDM
  else if (pdm_order_ == NPDM_FOURPDM) {
    if ( sweep_pos == 0 ) {
      // Nothing extra needed
    }
    else if ( sweep_pos == end_pos ) {
      lhs_dot_rhs_types_.insert( std::make_tuple(0,4,4) );
      lhs_dot_rhs_types_.insert( std::make_tuple(4,0,4) );
      lhs_dot_rhs_types_.insert( std::make_tuple(1,3,4) );
      lhs_dot_rhs_types_.insert( std::make_tuple(3,1,4) );
      lhs_dot_rhs_types_.insert( std::make_tuple(2,2,4) );
    }
  }
  // I don't think higher PDM's have any edge cases since we have 4-index ops max on a 1-site block??
  else abort();

  // Print out
  //pout << "=================================================================\n";
  //pout << "Possible block partitions:\n";
  //for ( auto it = lhs_dot_rhs_types_.begin(); it != lhs_dot_rhs_types_.end(); ++it ) {
  //  pout << get<0>(*it) << "," << get<1>(*it) << "," << get<2>(*it) << endl;
  //}
  //if ( sweep_pos == 0 )
  //  pout << "Added extra partitions for initial sweep position\n";
  //else if ( sweep_pos == end_pos )
  //  pout << "Added extra partitions for final sweep position\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::add_operator( int cre_ops1, int des_ops1, std::vector<CD> cd_type1 )
{

  std::vector<CD> cd_type2 (cd_type1.begin(), cd_type1.end());
  int cre_ops2 = cre_ops1;
  int des_ops2 = des_ops1;
  int order_;
  if (pdm_order_ == NPDM_PAIRMATRIX) {
    order_ = 1;
  } else if (pdm_order_ >= NPDM_ONEPDM || pdm_order_ <= NPDM_FOURPDM) {
    order_ = pdm_order_ - NPDM_ONEPDM + 1;
  } else {
    abort();
  }

  // Add creation operator as first tree branch
  if (cre_ops1 > 0) {
    cre_ops1--;
    cd_type1.push_back( CREATION );
    if ( cd_type1.size() < 2*order_ ) add_operator( cre_ops1, des_ops1, cd_type1 ) ;
  }

  // Add destruction operator as second tree branch
  if (des_ops2 > 0) {
    des_ops2--;
    cd_type2.push_back( DESTRUCTION );
    if ( cd_type2.size() < 2*order_ ) add_operator( cre_ops2, des_ops2, cd_type2 );
  }

  // Add only leaves of tree to final possiblities
  if ( cd_type1.size() == 2*order_ ){
    cre_des_types_.insert( cd_type1 );
  //pout << " cre_des_types1 \n";
  //for(int i=0; i < cd_type1.size();i++)
  //  pout <<cd_type1[i];
  //pout <<endl;
  }
  if ( cd_type2.size() == 2*order_ ){
  cre_des_types_.insert( cd_type2 );
  //pout << " cre_des_types2 \n";
  //for(int i=0; i < cd_type2.size();i++)
  //  pout <<cd_type2[i];
  //pout <<endl;

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::build_cre_des_types()
{
  // Empty vector
  std::vector<CD> cd_type;

//  // Populate the operators we have to work with
//  std::vector<CD> cre_tank(pdm_order_, CREATION);
//  std::vector<CD> des_tank(pdm_order_, DESTRUCTION);

  // Build up tree of valid creation-destruction strings by recursion
  if (pdm_order_ == NPDM_PAIRMATRIX) {
    add_operator(0, 2, cd_type);
  } else if (pdm_order_ >= NPDM_ONEPDM || pdm_order_ <= NPDM_FOURPDM){
    if(dmrginp.doimplicitTranspose()){
      cd_type.push_back( CREATION );
      add_operator( pdm_order_-NPDM_ONEPDM, pdm_order_-NPDM_ONEPDM+1, cd_type );
    }
    else{
      add_operator(pdm_order_-NPDM_ONEPDM+1,pdm_order_-NPDM_ONEPDM+1, cd_type);
    }
  } else {
    abort();
  }

  // Print out
  //pout << "=================================================================\n";
  //pout << "Creation/destruction patterns:\n";
  //for (auto iter = cre_des_types_.begin(); iter != cre_des_types_.end(); iter++) {
  //  print_cd_string(*iter); pout << endl;
  //}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//FIXME explain why (1) holds below???
// Ensure following properties of creation-destruction string on dot
// (1) All creation should be to left of destruction
// (2) No more than 2 creation or 2 destruction, since gives zero in fermion spin-half cases

bool Npdm_patterns::is_valid_dot_type( std::vector<CD> ops )
{
  // No more than 2 creation or destrucution
  if(!dmrginp.spinAdapted()){
    if(ops.size()>4) return false;
    int c = 0;
    int d = 0;
    for (auto op = ops.begin(); op != ops.end(); op++ ) {
      if ( *op == CREATION ) c++;
      if ( *op == DESTRUCTION ) d++;
    }
    if ( c > 2 ) return false;
    if ( d > 2 ) return false;
    //FIXME
    //return true;
    std::vector<int> ops_int;
    for( int i=0;i<ops.size();i++)
      ops_int.push_back(int(ops[i]));

    std::vector<int> order={0,1,0,1};
    if(ops_int.size()==4){
      if(ops_int==order) return true;
      return false;
    }
    if(ops_int.size()==3){
      order={0,0,1};
      if(ops_int==order) return true;
      order={0,1,0};
      if(ops_int==order) return true;
      order={0,1,1};
      if(ops_int==order) return true;
      order={1,0,1};
      if(ops_int==order) return true;
      return false;
    }
    return true;
    //std::vector<CD> sorted_ops( ops.begin(), ops.end() );
    //bool valid = std::is_sorted( sorted_ops.begin(), sorted_ops.end() );
    //if ( not valid ) return false;


  }
  if (ops.size() > 4) return false;
  int c = 0;
  int d = 0;
  for (auto op = ops.begin(); op != ops.end(); op++ ) {
    if ( *op == CREATION ) c++;
    if ( *op == DESTRUCTION ) d++;
    if ( c > 2 ) return false;
    if ( d > 2 ) return false;
  }

  // All creation must be on left
  std::vector<CD> sorted_ops( ops.begin(), ops.end() );
  bool valid = std::is_sorted( sorted_ops.begin(), sorted_ops.end() );
  if ( not valid ) return false;

  return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_patterns::is_rhs_gte_lhs( std::vector<int> & vec1, std::vector<int> & vec2 )
{
  assert(vec1.size() == vec2.size());
  //print_int_string(vec1);
  //print_int_string(vec2);  pout << std::endl;
  for (int i=0; i != vec1.size(); i++) {
    if (vec2[i] > vec1[i]) return true;
    if (vec2[i] < vec1[i]) return false;
  }
  // All indices equal
  return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Not all the creation-destruction patterns map to an irreducible index combination.
// Here we screen away those special cases.

bool Npdm_patterns::is_valid_ldr_type( std::map< char, std::vector<CD> > & cd_pattern )
{
  // To do the test, associate some arbitrary site indices consistent with the block partitioning
  std::vector< std::pair<CD,int> > opstring;
  int i = 1;
  // Add in LHS
  for (auto cd = cd_pattern.at('l').begin(); cd != cd_pattern.at('l').end(); cd++ ) {
    opstring.push_back( std::make_pair( *cd, i++ ) );
  }
  // Add in Dot
  for (auto cd = cd_pattern.at('d').begin(); cd != cd_pattern.at('d').end(); cd++ ) {
    opstring.push_back( std::make_pair( *cd, i ) );
  }
  // Add in RHS
  i++;
  for (auto cd = cd_pattern.at('r').begin(); cd != cd_pattern.at('r').end(); cd++ ) {
    opstring.push_back( std::make_pair( *cd, i++ ) );
  }

  // Sort cre-des operator string as close to irreducible form as possible (creation all on left)
  std::sort( opstring.begin(), opstring.end() );

  // Split into creation and destruction halves
  std::vector<int> cvec, dvec;
  if (pdm_order_ == NPDM_PAIRMATRIX) {
    for (auto it = opstring.begin() + 2; it != opstring.end(); it++ ) {
      dvec.push_back( it->second );
    }
  } else if (pdm_order_ >= NPDM_ONEPDM || pdm_order_ <= NPDM_FOURPDM) {
    for (auto it = opstring.begin(); it != opstring.begin() + pdm_order_-NPDM_ONEPDM+1; it++) {
      cvec.push_back( it->second );
    }
    for (auto it = opstring.begin() + pdm_order_; it != opstring.end(); it++ ) {
      dvec.push_back( it->second );
    }
  } else {
    abort();
  }

  // Test if dvec >= cvec in sense of irreducible operator string generation (triangular loop)
  bool valid = is_rhs_gte_lhs( cvec, dvec );
  if(dmrginp.doimplicitTranspose() && dmrginp.spinAdapted())
    // if implicit Transpose is not used, lhs can be greater than rhs
    // when spin-adpated is closed, dot site has more than one orbitals. Therefore, each orbital number in opstring is different.
    // Since the first one must be creator. rhs is must large than lhs. 
  {
    if ( not valid ) return false;
  }

  return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::build_ldr_cd_types( int sweep_pos, int end_pos )
{

  ldr_cd_types_.clear();
  //pout << "=================================================================\n";
  //pout << "Spin-1/2 fermionic operator patterns for DMRG blocks:\n";

  // Loop over LHS, Dot, RHS patterns
  for (auto ldr_iter = lhs_dot_rhs_types_.begin(); ldr_iter != lhs_dot_rhs_types_.end(); ldr_iter++) {
    int ilhs = std::get<0>(*ldr_iter);
    int idot = std::get<1>(*ldr_iter);
    int irhs = std::get<2>(*ldr_iter);
    //pout << ilhs << "," << idot << "," << irhs << "\n";
    // Loop over creation-destruction patterns
    for (auto cd_iter = cre_des_types_.begin(); cd_iter != cre_des_types_.end(); cd_iter++) {

      // Split CD pattern into LHS, Dot and RHS
      std::vector<CD> lhs_cd( cd_iter->begin() , (cd_iter->begin() + ilhs) );
      std::vector<CD> dot_cd( cd_iter->begin() + ilhs, cd_iter->begin() + ilhs + idot);
      std::vector<CD> rhs_cd( cd_iter->begin() + ilhs + idot, cd_iter->end() );

      // Only allow if it's a valid dot pattern 
      if ( not is_valid_dot_type( dot_cd ) ) continue;
      // Edge case: lhs == dot:
      if ( ( sweep_pos == 0 ) && ( not is_valid_dot_type( lhs_cd ) ) ) continue;
      // Edge case: rhs == dot:
      if ( ( sweep_pos == end_pos ) && ( not is_valid_dot_type( rhs_cd ) ) ) continue;

      // Combine together
      std::map< char, std::vector<CD> > cd_pattern;
      cd_pattern['l'] = lhs_cd;
      cd_pattern['d'] = dot_cd;
      cd_pattern['r'] = rhs_cd;

      // Only allow if it's a valid full pattern
      if ( not is_valid_ldr_type( cd_pattern ) ) continue;

      ldr_cd_types_.insert( cd_pattern );
      //if ( lhs_cd.size() != 0 ) lhs_cd_types_.insert( lhs_cd );
      //if ( dot_cd.size() != 0 ) dot_cd_types_.insert( dot_cd );
      //if ( rhs_cd.size() != 0 ) rhs_cd_types_.insert( rhs_cd );
      lhs_cd_types_.insert( lhs_cd );
      dot_cd_types_.insert( dot_cd );
      rhs_cd_types_.insert( rhs_cd );


    }
  }
// DEBUG add extra patterns
//std::map< char, std::vector<CD> > cd_pattern;
//cd_pattern['l'] = { CREATION, DESTRUCTION, DESTRUCTION, DESTRUCTION };
//cd_pattern['d'] = { };
//cd_pattern['r'] = { CREATION, CREATION };
////cd_pattern['r'] = { };
//ldr_cd_types_.insert( cd_pattern );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//
//  // Print out
//  pout << "----------------------------------\n";
//  for (auto iter = ldr_cd_types_.begin(); iter != ldr_cd_types_.end(); iter++) {
//    print_cd_string( iter->at('l') );
//    print_cd_string( iter->at('d') );
//    print_cd_string( iter->at('r') );
//    pout << std::endl;
//  }
//  pout << "Number of unique patterns = " << ldr_cd_types_.size() << std::endl;
//
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::print_int_string( const std::vector<int> & vec )
{
  pout << "(";
  for (auto op = vec.begin(); op != vec.end(); op++) {
     pout << *op;
  }
  pout << ")";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::print_cd_string( const std::vector<CD> & cdvec )
{
  char cd;
  pout << "(";
  for (auto op = cdvec.begin(); op != cdvec.end(); op++) {
     cd = '.';
     if (*op == CREATION) cd = '+';
     pout << cd;
  }
  pout << ")";

}

//===========================================================================================================================================================

}
}
