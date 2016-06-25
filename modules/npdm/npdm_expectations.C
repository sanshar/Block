/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/serialization/map.hpp>
#include "MatrixBLAS.h"
#include "pario.h"
#include "npdm_expectations_engine.h"
#include "npdm_expectations.h"
#include "pario.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Npdm_expectations::Npdm_expectations( Npdm_spin_adaptation& spin_adaptation, Npdm_patterns& npdm_patterns, 
                                      const NpdmOrder order, Wavefunction & wavefunction0, Wavefunction & wavefunction1, const SpinBlock & big )
: spin_adaptation_(spin_adaptation),
  npdm_patterns_(npdm_patterns),
  npdm_order_(order),
  wavefunction_0(wavefunction0),
  wavefunction_1(wavefunction1),
  big_(big)
{ }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// The pattern generator for the non-redundant NPDM elements leads to duplicates if indices are repeated, so some can be skipped explicitly.
// If the normal-ordered string is not of non-redundant form, then the original string produces duplicates when permutations are applied.

bool Npdm_expectations::screen_op_string_for_duplicates( const std::string& op, const std::vector<int>& indices )
{
  std::string CD;
  for (auto it = op.begin(); it != op.end(); ++it) {
    if ( (*it == 'C') || (*it == 'D') ) CD.push_back(*it);
  }

  if ( indices.size() == 2 ) {
    // 1PDM case
    return false;
  }
  else if ( indices.size() == 4 ) {
    // 2PDM case
    return npdm_patterns_.screen_2pdm_strings( indices, CD ); 
  }
  else if ( indices.size() == 6 ) {
    // 3PDM case
    return npdm_patterns_.screen_3pdm_strings( indices, CD ); 
  }
  else if ( indices.size() == 8 ) {
    // 4PDM case
    return npdm_patterns_.screen_4pdm_strings( indices, CD ); 
  }
  else abort();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::get_full_op_string( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps, 
                                            std::string& op_string, std::vector<int>& indices )
{

  // Set up npdm element indices
  indices.clear();
  indices.insert( indices.end(), lhsOps.indices_.begin(), lhsOps.indices_.end() );
  indices.insert( indices.end(), dotOps.indices_.begin(), dotOps.indices_.end() );
  indices.insert( indices.end(), rhsOps.indices_.begin(), rhsOps.indices_.end() );
  //pout << "lhs indices = "; for (auto it = lhsOps.indices_.begin(); it != lhsOps.indices_.end(); ++it) { pout << *it << " "; } pout << std::endl;
  //pout << "dot indices = "; for (auto it = dotOps.indices_.begin(); it != dotOps.indices_.end(); ++it) { pout << *it << " "; } pout << std::endl;
  //pout << "rhs indices = "; for (auto it = rhsOps.indices_.begin(); it != rhsOps.indices_.end(); ++it) { pout << *it << " "; } pout << std::endl;
  //pout << "spatial indices = "; for (auto it = indices.begin(); it != indices.end(); ++it) { pout << *it << " "; } pout << std::endl;

  // Record how tensor operator is constructed from (compound) block operators
  std::string build_pattern = "(";
  build_pattern.reserve( lhsOps.build_pattern_.size() + dotOps.build_pattern_.size() + rhsOps.build_pattern_.size() + 2 );
  build_pattern.insert( build_pattern.end(), lhsOps.build_pattern_.begin(), lhsOps.build_pattern_.end() );
  build_pattern.insert( build_pattern.end(), dotOps.build_pattern_.begin(), dotOps.build_pattern_.end() );
  build_pattern.push_back( ')' );
  build_pattern.insert( build_pattern.end(), rhsOps.build_pattern_.begin(), rhsOps.build_pattern_.end() );

  // Combine indices and build_pattern into one string
  op_string.clear();
  std::vector<int> idx = indices;
  for (auto it = build_pattern.begin(); it != build_pattern.end(); ++it) {
    op_string.push_back(*it);
    if ( (*it == 'C') || (*it == 'D') ) {
      std::string index = boost::lexical_cast<string>( idx.at(0) );
      op_string.append(index);
      idx.erase( idx.begin() );  
    }
  }
//pout << op_string << std::endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double Npdm_expectations::contract_spin_adapted_operators( int ilhs, int idot, int irhs, 
                                                           NpdmSpinOps_base& lhsOps, NpdmSpinOps_base& rhsOps, NpdmSpinOps_base& dotOps )
{
//FIXME is null_deleter() necessary below?
  SparseMatrix* null = 0; 
  double expectation = 0;

  boost::shared_ptr<SparseMatrix> lhsOp, dotOp, rhsOp;
  if ( lhsOps.opReps_.size() > 0 ) lhsOp = lhsOps.opReps_.at(ilhs); if ( dotOps.opReps_.size() > 0 ) dotOp = dotOps.opReps_.at(idot); if ( rhsOps.opReps_.size() > 0 ) rhsOp = rhsOps.opReps_.at(irhs);

  if(SpinAdapted::dmrginp.doimplicitTranspose()==false){
    if(lhsOps.transpose_ || dotOps.transpose_ || rhsOps.transpose_)
    {
      pout << "Transposeview could not be used if m_implicitTranspose is false" <<endl;
      abort();
    }
  }
    
  // We need to distinguish cases where one or more blocks has an empty operator string
  // X_X_X
  if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *dotOp, *rhsOp, big_);
  }
  // X_X_0
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *dotOp, *null, big_);
  }
  // X_0_X
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *null, *rhsOp, big_);
  }
  // 0_X_X
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *null, *dotOp, *rhsOp, big_);
  }
  // X_0_0
  else if ( (lhsOps.opReps_.size() > 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview lhsOpTr = Transposeview(*lhsOp);
    if ( lhsOps.transpose_ ) lhsOp = boost::shared_ptr<SparseMatrix>( &lhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *lhsOp, *null, *null, big_);
  }
  // 0_X_0
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() > 0) && (rhsOps.opReps_.size() == 0) ) {
    Transposeview dotOpTr = Transposeview(*dotOp);
    if ( dotOps.transpose_ ) dotOp = boost::shared_ptr<SparseMatrix>( &dotOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *null, *dotOp, *null, big_);
  }
  // 0_0_X
  else if ( (lhsOps.opReps_.size() == 0) && (dotOps.opReps_.size() == 0) && (rhsOps.opReps_.size() > 0) ) {
    Transposeview rhsOpTr = Transposeview(*rhsOp);
    if ( rhsOps.transpose_ ) rhsOp = boost::shared_ptr<SparseMatrix>( &rhsOpTr, boostutils::null_deleter() );
    expectation = spinExpectation(wavefunction_0, wavefunction_1, *null, *null, *rhsOp, big_);
  }
  else abort();

  // Modify new element with sign factors and return
  double factor = lhsOps.factor_ * dotOps.factor_ * rhsOps.factor_;
//pout << "expectation, factor = " << expectation << ", " << factor  << std::endl;
  return expectation*factor;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_expectations::test_for_singlet( int ilhs, int idot, int irhs, NpdmSpinOps_base& lhsOps, NpdmSpinOps_base& rhsOps, NpdmSpinOps_base& dotOps )
{
  int lhs2S, rhs2S, dot2S;
  // LHS
  if ( lhsOps.opReps_.size() == 0 ) 
    lhs2S = 0;
  else
    lhs2S = lhsOps.opReps_.at(ilhs)->get_deltaQuantum(0).get_s().getirrep();
  // RHS
  if ( rhsOps.opReps_.size() == 0 ) 
    rhs2S = 0;
  else
    rhs2S = rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s().getirrep();
  // DOT
  if ( dotOps.opReps_.size() == 0 ) 
    dot2S = 0;
  else
    dot2S = dotOps.opReps_.at(idot)->get_deltaQuantum(0).get_s().getirrep();

  // Couple LHS and Dot spin angular momenta and see if any equal RHS  
  for (int twoS = std::abs(lhs2S - dot2S); twoS <= ( lhs2S + dot2S ); twoS += 2 ) {
    if ( twoS == rhs2S ) {
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double DotProduct_spincorrection(const Wavefunction& w1, const Wavefunction& w2, const SpinBlock& big)
{
  // After multipling ket by cd operator, it has the same basis with ket
  int leftOpSz = big.get_leftBlock()->get_ketStateInfo().quanta.size();
  int rightOpSz = big.get_rightBlock()->get_braStateInfo().quanta.size();
  const StateInfo* rS = big.get_braStateInfo().rightStateInfo, *lS = big.get_ketStateInfo().leftStateInfo;

  double output = 0.0;
  SpinQuantum Q= w2.get_deltaQuantum(0);
  for (int lQ =0; lQ < leftOpSz; lQ++)
    for (int rQ = 0; rQ < rightOpSz; rQ++) {
      if (w1.allowed(lQ, rQ) && w2.allowed(lQ, rQ))
      {
	      double b1b2 = MatrixDotProduct(w1(lQ, rQ), w2(lQ, rQ));

              if(abs(b1b2) > NUMERICAL_ZERO )
              {

	        b1b2 *= dmrginp.get_ninej()(lS->quanta[lQ].get_s().getirrep(), (-lS->quanta[lQ].get_s()).getirrep() , 0, 
	          				(-Q.get_s()).getirrep(), Q.get_s().getirrep(), 0,
	          				(-rS->quanta[rQ].get_s()).getirrep(), rS->quanta[rQ].get_s().getirrep() , 0);
	        b1b2 *= Symmetry::spatial_ninej(lS->quanta[lQ].get_symm().getirrep(), (-lS->quanta[lQ].get_symm()).getirrep() , 0, 
	          				(-Q.get_symm()).getirrep(), Q.get_symm().getirrep(), 0,
	          				(-rS->quanta[rQ].get_symm()).getirrep(), rS->quanta[rQ].get_symm().getirrep() , 0);
	        b1b2 /= dmrginp.get_ninej()((-rS->quanta[rQ].get_s()).getirrep(), rS->quanta[rQ].get_s().getirrep() , 0, 
	          				(-Q.get_s()).getirrep(), 0, (-Q.get_s()).getirrep(),
	          				lS->quanta[lQ].get_s().getirrep(), rS->quanta[rQ].get_s().getirrep() , Q.get_s().getirrep());
	        b1b2 /= Symmetry::spatial_ninej((-rS->quanta[rQ].get_symm()).getirrep(), rS->quanta[rQ].get_symm().getirrep() , 0, 
	          				(-Q.get_symm()).getirrep(), 0, (-Q.get_symm()).getirrep(),
	          				lS->quanta[lQ].get_symm().getirrep(), rS->quanta[rQ].get_symm().getirrep() , Q.get_symm().getirrep());
	        b1b2 /= dmrginp.get_ninej()(lS->quanta[lQ].get_s().getirrep(), (-lS->quanta[lQ].get_s()).getirrep() , 0, 
	          				0, Q.get_s().getirrep(), Q.get_s().getirrep(),
	          				lS->quanta[lQ].get_s().getirrep(), rS->quanta[rQ].get_s().getirrep() , Q.get_s().getirrep());
	        b1b2 /= Symmetry::spatial_ninej(lS->quanta[lQ].get_symm().getirrep(), (-lS->quanta[lQ].get_symm()).getirrep() , 0, 
						0, Q.get_symm().getirrep(), Q.get_symm().getirrep(),
						lS->quanta[lQ].get_symm().getirrep(), rS->quanta[rQ].get_symm().getirrep() , Q.get_symm().getirrep());

              //  b1b2 *=Aop.get_scaling(rS->quanta[rQ],lS->quanta[lQ]);
              //  b1b2 /=Transposeview(Aop).get_scaling(lS->quanta[lQ],rS->quanta[rQ]);
              }
	      output += b1b2;
       
      }	
    }
return output;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double Npdm_expectations::build_nonspin_adapted_singlet_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  // IMPORTANT: generate spin-components in the same order as RHS of linear equation solver in npdm_set_up_linear_equations routine
  // i.e. in accordance with the operator string build_pattern
  if(dmrginp.npdm_intermediate() && (npdm_order_== NPDM_NEVPT2 || npdm_order_== NPDM_THREEPDM || npdm_order_== NPDM_FOURPDM))
  {
    assert(false);
    return 0;
  }
  else
    {
      return contract_spin_adapted_operators( 0, 0, 0, lhsOps, rhsOps, dotOps) ;
    }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::build_spin_adapted_singlet_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  expectations_.clear();

  // IMPORTANT: generate spin-components in the same order as RHS of linear equation solver in npdm_set_up_linear_equations routine
  // i.e. in accordance with the operator string build_pattern
  int hilhs = lhsOps.opReps_.size();
  int hirhs = rhsOps.opReps_.size();
  int hidot = dotOps.opReps_.size();
  std::map<std::vector<int>, Wavefunction> leftwaves;
  std::map<std::vector<int>, Wavefunction> rightwaves;
  if(dmrginp.npdm_intermediate() && (npdm_order_== NPDM_NEVPT2 || npdm_order_== NPDM_THREEPDM || npdm_order_== NPDM_FOURPDM))
    assert(false);
  for (int irhs = 0; irhs < std::max(1,hirhs); ++irhs) {
    for (int idot = 0; idot < std::max(1,hidot); ++idot) {
      for (int ilhs = 0; ilhs < std::max(1,hilhs); ++ilhs) {
        // Screen operator combinations that do not combine to give a singlet
        bool singlet = test_for_singlet( ilhs, idot, irhs, lhsOps, rhsOps, dotOps );
        // Contract operators to produce spin-adapted expectation value
        if ( singlet ) 
        {
          if(dmrginp.npdm_intermediate() && (npdm_order_== NPDM_NEVPT2 || npdm_order_== NPDM_THREEPDM || npdm_order_== NPDM_FOURPDM))
          {
            std::vector<int> spin;
            //spin.push_back(hilhs==0? 0 : lhsOps.opReps_.at(ilhs)->get_deltaQuantum(0).get_s().getirrep());
            //spin.push_back(hidot==0? 0 : dotOps.opReps_.at(idot)->get_deltaQuantum(0).get_s().getirrep());
            //spin.push_back(hirhs==0? 0 : (-(rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s())).getirrep());
            spin.push_back(ilhs);
            spin.push_back(idot);
            spin.push_back((-(rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s())).getirrep());
            Wavefunction& lw= leftwaves.at(spin);
            {
              assert(hirhs>0);
              Wavefunction& rw= rightwaves.at(std::vector<int>(1,irhs));
              expectations_.push_back(DotProduct_spincorrection(lw,rw,big_));
            }
          }
          else
            expectations_.push_back( contract_spin_adapted_operators( ilhs, idot, irhs, lhsOps, rhsOps, dotOps) ); }
      }
    }
  }

  assert (expectations_.size() > 0);
//pout << "---------------------------------\n";
//pout << "spin-adapted expectations =\n";
////pout << "mpirank = " << mpigetrank() << endl;
//for (auto it = expectations_.begin(); it != expectations_.end(); ++it) {
//  pout << "p" << mpigetrank() << ":  " << *it << std::endl;
//}
//pout << "---------------------------------\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::build_spin_adapted_singlet_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps, std::map<std::vector<int>, Wavefunction>& leftwaves, std::map<std::vector<int>, Wavefunction>& rightwaves)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  expectations_.clear();

  // IMPORTANT: generate spin-components in the same order as RHS of linear equation solver in npdm_set_up_linear_equations routine
  // i.e. in accordance with the operator string build_pattern
  int hilhs = lhsOps.opReps_.size();
  int hirhs = rhsOps.opReps_.size();
  int hidot = dotOps.opReps_.size();
  for (int irhs = 0; irhs < std::max(1,hirhs); ++irhs) {
    for (int idot = 0; idot < std::max(1,hidot); ++idot) {
      for (int ilhs = 0; ilhs < std::max(1,hilhs); ++ilhs) {
        // Screen operator combinations that do not combine to give a singlet
        bool singlet = test_for_singlet( ilhs, idot, irhs, lhsOps, rhsOps, dotOps );
        // Contract operators to produce spin-adapted expectation value
        if ( singlet ) 
        {
            std::vector<int> spin;
            spin.push_back(ilhs);
            spin.push_back(idot);
            spin.push_back((-(rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s())).getirrep());
            Wavefunction& lw= leftwaves.at(spin);
            Wavefunction& rw= rightwaves.at(std::vector<int>(1,irhs));
            expectations_.push_back(DotProduct_spincorrection(lw,rw,big_));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector< std::pair< std::vector<int>, double > > 
Npdm_expectations::get_nonspin_adapted_expectations( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps )
{
  // Initialize dimension of spin-adapted to non-spin-adapted transformation
  int dim = 0;
  if(dmrginp.spinAdapted()){
  if ( npdm_order_ == NPDM_ONEPDM ) dim = 2;
  else if ( npdm_order_ == NPDM_TWOPDM ) dim = 6;
  else if ( npdm_order_ == NPDM_THREEPDM ) dim = 20;
  else if ( npdm_order_ == NPDM_FOURPDM ) dim = 70;
  else abort();
  }
  else{
    //if ( npdm_order_ == 1 ) dim = 1;
    //else if ( npdm_order_ == 2 ) dim = 1;
    //else if ( npdm_order_ == 3 ) dim = 20;
   // else if ( npdm_order_ == 4 ) dim = 70;
    //else abort();

  }
  std::vector< std::pair< std::vector<int>, double > > new_pdm_elements;

  // Get operator build string. e.g. (C2C4)(D5D6)
  std::string op_string;
  std::vector<int> indices;
  get_full_op_string( lhsOps, rhsOps, dotOps, op_string, indices );

  // Screen away unwanted strings (e.g. those that produce duplicate NPDM elements)
  //for( int i=0;i<indices.size();i++)
  //  pout << indices[i]<<',';
  //pout <<endl;
  if ( screen_op_string_for_duplicates( op_string, indices ) ) return new_pdm_elements;


  if(!dmrginp.spinAdapted()){
    std::vector<int> cd_order;
    int k=0;
    for ( auto it = op_string.begin(); it != op_string.end(); ++it ) {
      if (*it=='C') cd_order.push_back((indices[k++]));
      if (*it=='D') cd_order.push_back(1000+indices[k++]);
    }
    //sort cd_order; bubble sort
    int parity=1;
    double tmp;
    for( int i=cd_order.size()-1;i >0;i--)
      for(int j=0; j <i;j++){
        if(cd_order[j]>cd_order[j+1]){
          tmp=cd_order[j+1];
          cd_order[j+1]=cd_order[j];
          cd_order[j]=tmp;
          parity*=-1;
          
        }
      }

    //change cd_order back to spinorbital number
    //And check whether total S_z of operator is zero
    int sz=0;
    for(int i=0;i<cd_order.size();i++){
      if (cd_order[i] < 1000) {
        if(cd_order[i]%2) sz+=1;
        else sz+=-1;
      } else {
        cd_order[i]-=1000;
        if(cd_order[i]%2) sz+=-1;
        else sz+=1;
      }
    }
    //pout << "nonspinadapted spinorbital ordered: ";
    //for( int i=0;i<cd_order.size();i++)
    //  pout << cd_order[i]<<',';
    //pout <<endl;
    //pout <<"sz: " <<sz<<endl;
    if(sz!=0) return new_pdm_elements;
    double nonspinvalue;

   // pout <<"size of operator: " <<lhsOps.opReps_.size()<<','<<rhsOps.opReps_.size()<<','<<dotOps.opReps_.size()<<endl;
      nonspinvalue=build_nonspin_adapted_singlet_expectations( lhsOps, rhsOps, dotOps);
      {

        boost::shared_ptr<SparseMatrix> lhsOp, dotOp, rhsOp;
          if(lhsOps.transpose_) lhsOp=boost::shared_ptr<SparseMatrix>(new Transposeview(*lhsOps.opReps_.at(0)));
          else lhsOp= lhsOps.opReps_.at(0);

          if(dotOps.transpose_) dotOp=boost::shared_ptr<SparseMatrix>(new Transposeview(*dotOps.opReps_.at(0)));
          else dotOp= dotOps.opReps_.at(0);

          if(rhsOps.transpose_) rhsOp=boost::shared_ptr<SparseMatrix>(new Transposeview(*rhsOps.opReps_.at(0)));
          else rhsOp= rhsOps.opReps_.at(0);

      int ls= !lhsOps.transpose_? lhsOps.opReps_.at(0)->get_deltaQuantum(0).get_s().getirrep(): (-lhsOps.opReps_.at(0)->get_deltaQuantum(0).get_s()).getirrep();
      int ds= !dotOps.transpose_? dotOps.opReps_.at(0)->get_deltaQuantum(0).get_s().getirrep(): (-dotOps.opReps_.at(0)->get_deltaQuantum(0).get_s()).getirrep();
          int total_spin= ls+ds;
          Cre AOp;
          FormLeftOp(big_.get_leftBlock(), *lhsOp, *dotOp, AOp, total_spin);
          Wavefunction opw2, opw_0, opw3;
          vector<SpinQuantum> dQ = wavefunction_0.get_deltaQuantum();
          opw_0.initialise(dQ[0], &big_, true);
          operatorfunctions::TensorMultiply(big_.get_leftBlock(), AOp, *rhsOp, &big_, wavefunction_0, opw_0,dQ[0], lhsOps.factor_*dotOps.factor_*rhsOps.factor_);



          opw2.initialise(dQ[0]-AOp.get_deltaQuantum(0), &big_, true);
          opw3.initialise(dQ[0]+rhsOp->get_deltaQuantum(0), &big_, true);
          //Left part of intermediate wavefuntion should multiply transpose of left ops.
          operatorfunctions::TensorMultiply(big_.get_leftBlock(), Transposeview(AOp), &big_, wavefunction_1, opw2,-AOp.get_deltaQuantum(0), lhsOps.factor_*dotOps.factor_);
          operatorfunctions::TensorMultiply(big_.get_rightBlock(), *rhsOp, &big_, wavefunction_1, opw3,-AOp.get_deltaQuantum(0), lhsOps.factor_*dotOps.factor_);
      }

    new_pdm_elements.push_back(std::make_pair(cd_order,nonspinvalue*parity));
    return new_pdm_elements;
  }
  
  // Contract spin-adapted spatial operators and build singlet expectation values
  build_spin_adapted_singlet_expectations( lhsOps, rhsOps, dotOps);

  // Now transform to non-spin-adapted spin-orbital representation
  // b holds the spin-adapted expectation values (we only care about the singlets)
  // Note the spin-order of elements in b follows a convention
  ColumnVector x(dim), b(dim);
  // x holds the non-spin-adapted expectation values
  x=0.0;

  // Transformation matrix
  Matrix A(dim,dim);
  // Vector of spin-orbital indices ordered according to A
  std::vector< std::vector<int> > so_indices(dim);

  // Parse operator string and set up linear equations
  spin_adaptation_.npdm_set_up_linear_equations(dim, op_string, indices, expectations_, A, b, so_indices );

  // Solve A.x = b to get non-spin-adapted expectations in x
  xsolve_AxeqB(A, b, x);

  // Package transformed elements into container and return
  for (int i=0; i < so_indices.size(); ++i) {
    new_pdm_elements.push_back( std::make_pair(so_indices[i], x(i+1)) );
  } 

//pout << "x vector:\n";
//for (int i=1; i<(dim+1); ++i) { 
//    pout << i << "\t\t" << x(i) << std::endl;
//}
  
  return new_pdm_elements;
}

std::vector< std::pair< std::vector<int>, double > > 
Npdm_expectations::get_nonspin_adapted_expectations(NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & rhsOps, NpdmSpinOps_base & dotOps, std::map<std::vector<int>, Wavefunction>& leftwaves, std::map<std::vector<int>, Wavefunction>& rightwaves )
{
  // Initialize dimension of spin-adapted to non-spin-adapted transformation
  int dim = 0;
  assert(dmrginp.spinAdapted());
  if(dmrginp.spinAdapted()){
  if ( npdm_order_ == NPDM_ONEPDM ) dim = 2;
  else if ( npdm_order_ == NPDM_TWOPDM ) dim = 6;
  else if ( npdm_order_ == NPDM_THREEPDM ) dim = 20;
  else if ( npdm_order_ == NPDM_FOURPDM ) dim = 70;
  else abort();
  }
  std::vector< std::pair< std::vector<int>, double > > new_pdm_elements;

  // Get operator build string. e.g. (C2C4)(D5D6)
  std::string op_string;
  std::vector<int> indices;
  get_full_op_string( lhsOps, rhsOps, dotOps, op_string, indices );

  // Screen away unwanted strings (e.g. those that produce duplicate NPDM elements)
  //for( int i=0;i<indices.size();i++)
  //  pout << indices[i]<<',';
  //pout <<endl;
  if ( screen_op_string_for_duplicates( op_string, indices ) ) return new_pdm_elements;

  // Contract spin-adapted spatial operators and build singlet expectation values
  build_spin_adapted_singlet_expectations(lhsOps, rhsOps, dotOps, leftwaves, rightwaves);

  // Now transform to non-spin-adapted spin-orbital representation
  // b holds the spin-adapted expectation values (we only care about the singlets)
  // Note the spin-order of elements in b follows a convention
  ColumnVector x(dim), b(dim);
  // x holds the non-spin-adapted expectation values
  x=0.0;

  // Transformation matrix
  Matrix A(dim,dim);
  // Vector of spin-orbital indices ordered according to A
  std::vector< std::vector<int> > so_indices(dim);

  // Parse operator string and set up linear equations
  spin_adaptation_.npdm_set_up_linear_equations(dim, op_string, indices, expectations_, A, b, so_indices );

  // Solve A.x = b to get non-spin-adapted expectations in x
  xsolve_AxeqB(A, b, x);

  // Package transformed elements into container and return
  for (int i=0; i < so_indices.size(); ++i) {
    new_pdm_elements.push_back( std::make_pair(so_indices[i], x(i+1)) );
  } 

  return new_pdm_elements;
}

void Npdm_expectations::get_op_string( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & dotOps, 
                                            std::string& op_string )
{

  // Set up npdm element indices
  std::vector<int> indices;
  //if(lhsOps.indices_.size()>0)
  indices.insert( indices.end(), lhsOps.indices_.begin(), lhsOps.indices_.end() );
  //if(dotOps.indices_.size()>0)
  indices.insert( indices.end(), dotOps.indices_.begin(), dotOps.indices_.end() );

  // Record how tensor operator is constructed from (compound) block operators
  std::string build_pattern;
  build_pattern.reserve( lhsOps.build_pattern_.size() + dotOps.build_pattern_.size() );
  //if(lhsOps.build_pattern_.size()>0)
  build_pattern.insert( build_pattern.end(), lhsOps.build_pattern_.begin(), lhsOps.build_pattern_.end() );
  //if(dotOps.build_pattern_.size()>0)
  build_pattern.insert( build_pattern.end(), dotOps.build_pattern_.begin(), dotOps.build_pattern_.end() );

  // Combine indices and build_pattern into one string
  op_string.clear();
  std::vector<int> idx = indices;
  for (auto it = build_pattern.begin(); it != build_pattern.end(); ++it) {
    op_string.push_back(*it);
    if ( (*it == 'C') || (*it == 'D') ) {
      std::string index = boost::lexical_cast<string>( idx.at(0) );
      op_string.append(index);
      idx.erase( idx.begin() );  
    }
  }
//pout << op_string << std::endl;

}

void Npdm_expectations::get_op_string( NpdmSpinOps_base & rhsOps, 
                                            std::string& op_string)
{

  // Set up npdm element indices
  std::vector<int> indices;
  indices.insert( indices.end(), rhsOps.indices_.begin(), rhsOps.indices_.end() );


  // Combine indices and build_pattern into one string
  op_string.clear();
  std::vector<int> idx = indices;
  for (auto it = rhsOps.build_pattern_.begin(); it != rhsOps.build_pattern_.end(); ++it) {
    op_string.push_back(*it);
    if ( (*it == 'C') || (*it == 'D') ) {
      std::string index = boost::lexical_cast<string>( idx.at(0) );
      op_string.append(index);
      idx.erase( idx.begin() );  
    }
  }
//pout << op_string << std::endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_expectations::compute_intermediate( NpdmSpinOps_base & lhsOps, NpdmSpinOps_base & dotOps, std::map<std::vector<int>, Wavefunction> & waves)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  SparseMatrix* null = 0; 
  vector<SpinQuantum> dQ = wavefunction_0.get_deltaQuantum();
  assert(dQ.size()==1 );
  assert(dQ[0].totalSpin.getirrep()== 0);
  int lindices= &lhsOps ? lhsOps.indices_.size(): 0;
  int dindices= &dotOps ? dotOps.indices_.size(): 0;
  int rindices= 2*npdm_order_-lindices- dindices;
  int hilhs = lhsOps.opReps_.size();
  int hidot = dotOps.opReps_.size();
  for (int idot = 0; idot < std::max(1,hidot); ++idot) {
    for (int ilhs = 0; ilhs < std::max(1,hilhs); ++ilhs) {
     // int ls= hilhs? lhsOps.opReps_.at(ilhs)->get_deltaQuantum(0).get_s().getirrep(): 0;
     // int ds= hidot? dotOps.opReps_.at(idot)->get_deltaQuantum(0).get_s().getirrep(): 0;
      int ls= !lhsOps.transpose_? lhsOps.opReps_.at(ilhs)->get_deltaQuantum(0).get_s().getirrep(): (-lhsOps.opReps_.at(ilhs)->get_deltaQuantum(0).get_s()).getirrep();
      int ds= !dotOps.transpose_? dotOps.opReps_.at(idot)->get_deltaQuantum(0).get_s().getirrep(): (-dotOps.opReps_.at(idot)->get_deltaQuantum(0).get_s()).getirrep();
      if((!dmrginp.spinAdapted()? std::abs(ls+ds): std::abs(ls-ds)) <= rindices )
      {
        boost::shared_ptr<SparseMatrix> lhsOp, dotOp;
        {
          assert(lhsOps.opReps_.size() > 0);
          if(lhsOps.transpose_) lhsOp=boost::shared_ptr<SparseMatrix>(new Transposeview(*lhsOps.opReps_.at(ilhs)));
          else lhsOp= lhsOps.opReps_.at(ilhs);
        }

        {
          assert(dotOps.opReps_.size() > 0);
          if(dotOps.transpose_) dotOp=boost::shared_ptr<SparseMatrix>(new Transposeview(*dotOps.opReps_.at(idot)));
          else dotOp= dotOps.opReps_.at(idot);
        }
        for(int total_spin= !dmrginp.spinAdapted()? ls+ds: std::abs(ls-ds);total_spin<= std::min(ls+ds,rindices); total_spin+=2)
        {
          Cre AOp;
          FormLeftOp(big_.get_leftBlock(), *lhsOp, *dotOp, AOp, total_spin);
          Wavefunction opw2;
          opw2.AllowQuantaFor(big_.get_leftBlock()->get_ketStateInfo(),big_.get_rightBlock()->get_braStateInfo(),dQ[0]-AOp.get_deltaQuantum(0));
          //Left part of intermediate wavefuntion should multiply transpose of left ops.
          operatorfunctions::braTensorMultiply(big_.get_leftBlock(), AOp, &big_, wavefunction_0, opw2, lhsOps.factor_*dotOps.factor_);
          std::vector<int> spin;
          //spin.push_back(lhsOp->get_deltaQuantum(0).get_s().getirrep()); spin.push_back(dotOp->get_deltaQuantum(0).get_s().getirrep());spin.push_back(total_spin);
          spin.push_back(ilhs); spin.push_back(idot);spin.push_back(total_spin);
//        waves.emplace(spin,opw2);
          waves.insert(std::make_pair(spin,opw2));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void Npdm_expectations::compute_intermediate( NpdmSpinOps_base & rhsOps, std::map<std::vector<int>, Wavefunction> &  waves)
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  SparseMatrix* null = 0; 
  vector<SpinQuantum> dQ = wavefunction_1.get_deltaQuantum();
  assert(dQ.size()==1 );
  assert(dQ[0].totalSpin.getirrep()== 0);
  int rindices= &rhsOps ? rhsOps.indices_.size(): 0;
  int hirhs = rhsOps.opReps_.size();
  for (int irhs = 0; irhs < hirhs; ++irhs) {
    int rs= !rhsOps.transpose_? rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s().getirrep(): (-rhsOps.opReps_.at(irhs)->get_deltaQuantum(0).get_s()).getirrep();
    if(std::abs(rs)<= 2*npdm_order_-rindices )
    {
      boost::shared_ptr<SparseMatrix> rhsOp;
      if(rhsOps.transpose_) rhsOp=boost::shared_ptr<SparseMatrix>(new Transposeview(*rhsOps.opReps_.at(irhs)));
      else rhsOp= rhsOps.opReps_.at(irhs);
      rs = rhsOp->get_deltaQuantum(0).get_s().getirrep();
      Wavefunction opw2;
      opw2.AllowQuantaFor(big_.get_leftBlock()->get_ketStateInfo(),big_.get_rightBlock()->get_braStateInfo(),dQ[0]+rhsOp->get_deltaQuantum(0));
          //Left part of intermediate wavefuntion should multiply transpose of left ops.
      operatorfunctions::TensorMultiply(big_.get_rightBlock(), *rhsOp, &big_, wavefunction_1, opw2, rhsOp->get_deltaQuantum(0), rhsOps.factor_);
      std::vector<int> spin;
      spin.push_back(irhs);
//    waves.emplace(spin,opw2);
      waves.insert(std::make_pair(spin,opw2));
    }

  }
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//===========================================================================================================================================================

}
}
