#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include "pario.h"
#include "npdm_patterns.h"

using namespace std;  // for pout only

namespace SpinAdapted {
namespace Npdm {

//===========================================================================================================================================================

Npdm_patterns::Npdm_patterns( int pdm_order, int sweep_pos, int end_pos )
: pdm_order_(pdm_order)
{
  build_lhs_dot_rhs_types( sweep_pos, end_pos );
  build_cre_des_types( );
  build_ldr_cd_types( sweep_pos, end_pos );

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Generates all required block partitionings for an NPDM

void Npdm_patterns::build_lhs_dot_rhs_types( int sweep_pos, int end_pos )
{

  //---------------
  // General case
  //---------------
  int lhs, rhs, dot, dotmax;
  for (lhs = pdm_order_; lhs >= 0; lhs--) {
    dotmax = 2*pdm_order_ - lhs;
    for (dot = dotmax; dot >= 1; dot--) {
      // Can have no more than 4 on the dot block
      if (dot > 4) continue;
      rhs = 2*pdm_order_ - dot - lhs;
      if ( rhs < pdm_order_ ) {
        lhs_dot_rhs_types_.insert( std::make_tuple(lhs,dot,rhs) );
      }
    }
  }

  //FIXME make sure no unneccesary patterns
  //---------------
  // Edge cases 
  //---------------
  // 2PDM
  if (pdm_order_ == 2) {
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
  else if (pdm_order_ == 3) {
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
  else if (pdm_order_ == 4) {
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
  else assert(false);

  // Print out
  pout << "=================================================================\n";
  pout << "Possible block partitions:\n";
  for ( auto it = lhs_dot_rhs_types_.begin(); it != lhs_dot_rhs_types_.end(); ++it ) {
    pout << std::get<0>(*it) << "," << std::get<1>(*it) << "," << std::get<2>(*it) << std::endl;
  }
  if ( sweep_pos == 0 )
    pout << "Added extra partitions for initial sweep position\n";
  else if ( sweep_pos == end_pos )
    pout << "Added extra partitions for final sweep position\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::add_operator( int cre_ops1, int des_ops1, std::vector<CD> cd_type1 )
{

  std::vector<CD> cd_type2 (cd_type1.begin(), cd_type1.end());
  int cre_ops2 = cre_ops1;
  int des_ops2 = des_ops1;

  // Add creation operator as first tree branch
  if (cre_ops1 > 0) {
    cre_ops1--;
    cd_type1.push_back( CREATION );
    if ( cd_type1.size() < 2*pdm_order_ ) add_operator( cre_ops1, des_ops1, cd_type1 ) ;
  }

  // Add destruction operator as second tree branch
  if (des_ops2 > 0) {
    des_ops2--;
    cd_type2.push_back( DESTRUCTION );
    if ( cd_type2.size() < 2*pdm_order_ ) add_operator( cre_ops2, des_ops2, cd_type2 );
  }

  // Add only leaves of tree to final possiblities
  if ( cd_type1.size() == 2*pdm_order_ ) cre_des_types_.insert( cd_type1 );
  if ( cd_type2.size() == 2*pdm_order_ ) cre_des_types_.insert( cd_type2 );

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
  cd_type.push_back( CREATION );
  add_operator( pdm_order_-1, pdm_order_, cd_type );

  // Print out
  pout << "=================================================================\n";
  pout << "Creation/destruction patterns:\n";
  for (auto iter = cre_des_types_.begin(); iter != cre_des_types_.end(); iter++) {
    print_cd_string(*iter); pout << std::endl;
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//FIXME explain why (1) holds below???
// Ensure following properties of creation-destruction string on dot
// (1) All creation should be to left of destruction
// (2) No more than 2 creation or 2 destruction, since gives zero in fermion spin-half cases

bool Npdm_patterns::is_valid_dot_type( std::vector<CD> ops )
{
  // No more than 2 creation or destrucution
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
  for (auto it = opstring.begin(); it != opstring.begin() + pdm_order_; it++) {
    cvec.push_back( it->second );
  }
  for (auto it = opstring.begin() + pdm_order_; it != opstring.end(); it++ ) {
    dvec.push_back( it->second );
  }

  // Test if dvec >= cvec in sense of irreducible operator string generation (triangular loop)
  bool valid = is_rhs_gte_lhs( cvec, dvec );
  if ( not valid ) return false;

  return true;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_patterns::build_ldr_cd_types( int sweep_pos, int end_pos )
{

  ldr_cd_types_.clear();
  pout << "=================================================================\n";
  pout << "Spin-1/2 fermionic operator patterns for DMRG blocks:\n";

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
      print_cd_string( lhs_cd );
      print_cd_string( dot_cd );
      print_cd_string( rhs_cd ); pout << "\n";

      ldr_cd_types_.insert( cd_pattern );
      if ( lhs_cd.size() != 0 ) lhs_cd_types_.insert( lhs_cd );
      if ( dot_cd.size() != 0 ) dot_cd_types_.insert( dot_cd );
      if ( rhs_cd.size() != 0 ) rhs_cd_types_.insert( rhs_cd );

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
