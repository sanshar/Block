/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "op_components.h"
#include "BaseOperator.h"
#include "spinblock.h"
#include "operatorfunctions.h"
#include "pario.h"

namespace SpinAdapted{
namespace Three_index_ops{

//FIXME take out common parts with 3- and 4- index functions

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<boost::shared_ptr<SparseMatrix> > get_ops_from_disk( std::ifstream& ifs, int size )
{
  // Only 3-index ops read here
  assert( size == 3 );
  std::vector<boost::shared_ptr<SparseMatrix> > opReps;

  // Read in all spin components for this set of spatial indices
  for ( int i=0; i<size; i++) {
    boost::shared_ptr<SparseMatrix> op (new Cre);
    boost::archive::binary_iarchive load_op(ifs);
    load_op >> *op;
    assert( op->get_built_on_disk() );
    opReps.push_back(op);
  } 
  return opReps;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void store_ops_on_disk( std::ofstream& ofs, std::vector<boost::shared_ptr<SparseMatrix> > spin_ops )
{
  // Only 3-index ops saved here
  assert( spin_ops.size() == 3 );

  // Store all spin components for this set of spatial indices, preserving order
  for ( int i = 0; i<spin_ops.size(); ++i ) {
    boost::shared_ptr<SparseMatrix>& op = spin_ops[i];
    op->set_built_on_disk() = true;
    boost::archive::binary_oarchive save_op(ofs);
    save_op << *op;
    // Deallocate memory for operator representation
    op->set_built() = false;
    op->CleanUp();
  } 
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void finish_tensor_trace( SpinBlock& b, SpinBlock* sysdot, SparseMatrix& sysdot_op, SparseMatrix& op, std::string& build_pattern )
{
  // Build and store new operator
  assert( ! op.get_built() );
  op.set_built() = true;
  op.set_build_pattern() = build_pattern;
//FIXME magic number 1
  op.set_deltaQuantum(1, op.get_quantum_ladder().at( build_pattern ).at(1) );
  op.allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  SpinAdapted::operatorfunctions::TensorTrace(sysdot, sysdot_op, &b, &(b.get_stateInfo()), op);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void finish_tensor_product( SpinBlock& b, SpinBlock* sysdot, 
                            const SparseMatrix& sysdot_op1, const SparseMatrix& sysdot_op2, SparseMatrix& op, 
                            bool include_parity, std::string& build_pattern )
{
  // Build and store new operator
  assert( ! op.get_built() );
  op.set_built() = true;
  op.set_build_pattern() = build_pattern;
//FIXME magic number 1
  op.set_deltaQuantum(1, op.get_quantum_ladder().at( build_pattern ).at(1) );
  op.allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  // Do tensor product
  double parity = 1.0;
  if ( include_parity ) parity = getCommuteParity( sysdot_op1.get_deltaQuantum(0), sysdot_op2.get_deltaQuantum(0), op.get_deltaQuantum(0) );
  SpinAdapted::operatorfunctions::TensorProduct(sysdot, sysdot_op1, sysdot_op2, &b, &(b.get_stateInfo()), op, parity);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_tensor_trace( const opTypes& optype, SpinBlock& big, SpinBlock* sysdot, std::ofstream& ofs,
                             const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo )
{
  // Get pointer to sparse operator array
  Op_component_base& sysdot_array = sysdot->get_op_array(optype);
  // Open filesystem if necessary
  std::ifstream ifs;
  if ( (! dmrginp.do_npdm_in_core()) && sysdot->size() > 1 ) ifs.open( sysdot_array.get_filename().c_str(), std::ios::binary );

//FIXME need reference?  don't want to copy?
  // Loop over all operator indices
  std::vector<boost::shared_ptr<SparseMatrix> > sysdot_ops;
//pout << "trace array.size = " << sysdot_array.get_size() << endl;
  for (int idx = 0; idx < sysdot_array.get_size(); ++idx) {
    if ( dmrginp.do_npdm_in_core() || sysdot->size() <= 1) 
      sysdot_ops = sysdot_array.get_local_element(idx);
    else
//FIXME size
      sysdot_ops = get_ops_from_disk( ifs, sysdot_array.get_local_element(0).size() );

    // Loop over spin-op components
    int i = sysdot_ops[0]->get_orbs()[0]; int j = sysdot_ops[0]->get_orbs()[1]; int k = sysdot_ops[0]->get_orbs()[2];
//pout << "i,j,k = " << i << "," << j << "," << k << endl;
    // In parallel calculations not all operators are built on each proc
    if ( ! big.get_op_array(optype).has_local_index(i,j,k) ) continue;
    std::vector<boost::shared_ptr<SparseMatrix> > new_ops = big.get_op_array(optype).get_element(i,j,k);
    for (int jdx=0; jdx < sysdot_ops.size(); jdx++) {
      boost::shared_ptr<SparseMatrix>& sysdot_op = sysdot_ops[jdx];
      assert( sysdot_op->get_built() );
      std::string build_pattern = sysdot_op->get_build_pattern();

      // Allocate and build new operator
      for (int sx=0; sx < new_ops.size(); sx++) {
        boost::shared_ptr<SparseMatrix>& op = new_ops[sx];
        std::vector<SpinQuantum> s1 = sysdot_op->get_quantum_ladder().at(build_pattern);
        std::vector<SpinQuantum> s2 = op->get_quantum_ladder().at(build_pattern);
        // Store spin component in correct location
        if ( s1 == s2 ) {
          finish_tensor_trace( big, sysdot, *sysdot_op, *op, build_pattern );
          // Renormalise operator
          op->renormalise_transform( rotateMatrix, stateinfo );
        }
      }
    }
    // Store spin-batch on disk 
    if ( ! dmrginp.do_npdm_in_core() ) store_ops_on_disk( ofs, new_ops );
  }

  // Close filesystem if necessary
  if ( ifs.is_open() ) ifs.close();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_tensor_trace( const opTypes& optype, SpinBlock& big, SpinBlock* sysdot, std::ofstream& ofs,
                             const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
{

  const SpinBlock* overlap_block = (big.get_leftBlock() == sysdot) ? big.get_rightBlock() : big.get_leftBlock();
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const boost::shared_ptr<SparseMatrix> Overlap = overlap_block->get_op_rep(OVERLAP, hq);
  // Get pointer to sparse operator array
  Op_component_base& sysdot_array = sysdot->get_op_array(optype);
  // Open filesystem if necessary
  std::ifstream ifs;
  if ( (! dmrginp.do_npdm_in_core()) && sysdot->size() > 1 ) ifs.open( sysdot_array.get_filename().c_str(), std::ios::binary );

//FIXME need reference?  don't want to copy?
  // Loop over all operator indices
  std::vector<boost::shared_ptr<SparseMatrix> > sysdot_ops;
//pout << "trace array.size = " << sysdot_array.get_size() << endl;
  for (int idx = 0; idx < sysdot_array.get_size(); ++idx) {
    if ( dmrginp.do_npdm_in_core() || sysdot->size() <= 1) 
      sysdot_ops = sysdot_array.get_local_element(idx);
    else
//FIXME size
      sysdot_ops = get_ops_from_disk( ifs, sysdot_array.get_local_element(0).size() );

    // Loop over spin-op components
    int i = sysdot_ops[0]->get_orbs()[0]; int j = sysdot_ops[0]->get_orbs()[1]; int k = sysdot_ops[0]->get_orbs()[2];
//pout << "i,j,k = " << i << "," << j << "," << k << endl;
    // In parallel calculations not all operators are built on each proc
    if ( ! big.get_op_array(optype).has_local_index(i,j,k) ) continue;
    std::vector<boost::shared_ptr<SparseMatrix> > new_ops = big.get_op_array(optype).get_element(i,j,k);
    for (int jdx=0; jdx < sysdot_ops.size(); jdx++) {
      boost::shared_ptr<SparseMatrix>& sysdot_op = sysdot_ops[jdx];
      assert( sysdot_op->get_built() );
      std::string build_pattern = sysdot_op->get_build_pattern();

      // Allocate and build new operator
      for (int sx=0; sx < new_ops.size(); sx++) {
        boost::shared_ptr<SparseMatrix>& op = new_ops[sx];
        std::vector<SpinQuantum> s1 = sysdot_op->get_quantum_ladder().at(build_pattern);
        std::vector<SpinQuantum> s2 = op->get_quantum_ladder().at(build_pattern);
        // Store spin component in correct location
        if ( s1 == s2 ) {

          // forwards control whether calculate the commute parity for two operators. Overlap's can commute with any operators.
          bool forwards = false;
          finish_tensor_product( big, sysdot, *sysdot_op, *Overlap, *op, forwards, build_pattern );
          // Renormalise operator
          op->renormalise_transform( leftMat, bra, rightMat, ket);
        }
      }
    }
    // Store spin-batch on disk 
    if ( ! dmrginp.do_npdm_in_core() ) store_ops_on_disk( ofs, new_ops );
  }

  // Close filesystem if necessary
  if ( ifs.is_open() ) ifs.close();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_1_2_tensor_products( bool forwards, const opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& big, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs,
                                    const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo )
{
  // (i | j,k ) partition
  //-------------------------
  Op_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
  assert ( (rhs_array.get_size() == 1) || (lhs_array.get_size() == 1) );
//pout << "tensor_1_2: rhs_size = " << rhs_array.get_size() << "; op = " << rhs_array.get_op_string() << endl;
//pout << "tensor_1_2: lhs_size = " << lhs_array.get_size() << "; op = " << lhs_array.get_op_string() << endl;

  // Initialize filesystem
  std::ifstream lhsifs;
  if (( ! dmrginp.do_npdm_in_core()) && lhsBlock->size() > 1) lhsifs.open( lhs_array.get_filename().c_str(), std::ios::binary );

  // Loop over all lhs operator indices
  for (int idx = 0; idx < lhs_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops;
    // Assume 2-index operators are available on this processor in core
    lhs_ops = lhs_array.get_local_element(idx);

    // Loop over all rhs operator indices
    for (int iidx = 0; iidx < rhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops = rhs_array.get_local_element(iidx);
      int i = rhs_ops[0]->get_orbs()[0];
      int j = lhs_ops[0]->get_orbs()[0]; int k = lhs_ops[0]->get_orbs()[1];
//pout << "i = " << i << endl;
//pout << "j,k = " << j << "," << k << endl;
      // In parallel calculations not all operators are built on each proc
      if ( ! big.get_op_array(optype).has_local_index(i,j,k) ) continue;
//pout << "building i,j,k = " << i << "," << j << "," << k << endl;
      std::vector<boost::shared_ptr<SparseMatrix> > vec = big.get_op_array(optype).get_element(i,j,k);
//pout << "got i,j,k\n";

      // Loop over lhs spin-op components
      for (int jdx=0; jdx < lhs_ops.size(); jdx++) {
        boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jdx];
        assert( lhs_op->get_built() );
        std::string build_23 = lhs_op->get_build_pattern();
//pout << build_23 << endl;
//pout << "getting spin_23\n";
        std::vector<SpinQuantum> spin_23 = lhs_op->get_quantum_ladder().at(build_23);
//pout << "got spin_23\n";

        // Loop over rhs spin-op components
        for (int jjdx=0; jjdx < rhs_ops.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jjdx];
          assert( rhs_op->get_built() );
          std::string build_1 = rhs_op->get_build_pattern();
          std::string build_pattern = "(" + build_1 + build_23 + ")";

          // Allocate and build new operator
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            // Select relevant spin component
            std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0) };
            if ( s == spin_23 ) {
              finish_tensor_product( big, rhsBlock, *rhs_op, *lhs_op, *op, forwards, build_pattern );
              // Renormalise operator
              op->renormalise_transform( rotateMatrix, stateinfo );
            }
          }
        }
      }
      // Store spin-batch on disk 
      if ( ! dmrginp.do_npdm_in_core() ) store_ops_on_disk( ofs, vec );
    }
  }
  if ( lhsifs.is_open() ) lhsifs.close();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_1_2_tensor_products( bool forwards, const opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& big, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs,
                                    const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  // (i | j,k ) partition
  //-------------------------
  Op_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
  assert ( (rhs_array.get_size() == 1) || (lhs_array.get_size() == 1) );
//pout << "tensor_1_2: rhs_size = " << rhs_array.get_size() << "; op = " << rhs_array.get_op_string() << endl;
//pout << "tensor_1_2: lhs_size = " << lhs_array.get_size() << "; op = " << lhs_array.get_op_string() << endl;

  // Initialize filesystem
  std::ifstream lhsifs;
  if (( ! dmrginp.do_npdm_in_core()) && lhsBlock->size() > 1) lhsifs.open( lhs_array.get_filename().c_str(), std::ios::binary );

  // Loop over all lhs operator indices
  for (int idx = 0; idx < lhs_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops;
    // Assume 2-index operators are available on this processor in core
    lhs_ops = lhs_array.get_local_element(idx);

    // Loop over all rhs operator indices
    for (int iidx = 0; iidx < rhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops = rhs_array.get_local_element(iidx);
      int i = rhs_ops[0]->get_orbs()[0];
      int j = lhs_ops[0]->get_orbs()[0]; int k = lhs_ops[0]->get_orbs()[1];
//pout << "i = " << i << endl;
//pout << "j,k = " << j << "," << k << endl;
      // In parallel calculations not all operators are built on each proc
      if ( ! big.get_op_array(optype).has_local_index(i,j,k) ) continue;
//pout << "building i,j,k = " << i << "," << j << "," << k << endl;
      std::vector<boost::shared_ptr<SparseMatrix> > vec = big.get_op_array(optype).get_element(i,j,k);
//pout << "got i,j,k\n";

      // Loop over lhs spin-op components
      for (int jdx=0; jdx < lhs_ops.size(); jdx++) {
        boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jdx];
        assert( lhs_op->get_built() );
        std::string build_23 = lhs_op->get_build_pattern();
//pout << build_23 << endl;
//pout << "getting spin_23\n";
        std::vector<SpinQuantum> spin_23 = lhs_op->get_quantum_ladder().at(build_23);
//pout << "got spin_23\n";

        // Loop over rhs spin-op components
        for (int jjdx=0; jjdx < rhs_ops.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jjdx];
          assert( rhs_op->get_built() );
          std::string build_1 = rhs_op->get_build_pattern();
          std::string build_pattern = "(" + build_1 + build_23 + ")";

          // Allocate and build new operator
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            // Select relevant spin component
            std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0) };
            if ( s == spin_23 ) {
              finish_tensor_product( big, rhsBlock, *rhs_op, *lhs_op, *op, forwards, build_pattern );
              // Renormalise operator
              op->renormalise_transform( leftMat, bra, rightMat, ket );
            }
          }
        }
      }
      // Store spin-batch on disk 
      if ( ! dmrginp.do_npdm_in_core() ) store_ops_on_disk( ofs, vec );
    }
  }
  if ( lhsifs.is_open() ) lhsifs.close();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_2_1_tensor_products( bool forwards, const opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& big, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs,
                                    const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo )
{
  // (i,j | k ) partition
  //-------------------------
  Op_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
  assert ( (rhs_array.get_size() == 1) || (lhs_array.get_size() == 1) );

  // Initialize filesystem
  std::ifstream rhsifs;
  if ( (! dmrginp.do_npdm_in_core()) && rhsBlock->size() > 1 ) rhsifs.open( rhs_array.get_filename().c_str(), std::ios::binary );
//pout << "tensor_2_1: rhs_size = " << rhs_array.get_size() << "; op = " << rhs_array.get_op_string() << endl;
//pout << "tensor_2_1: lhs_size = " << lhs_array.get_size() << "; op = " << lhs_array.get_op_string() << endl;

  // Loop over all rhs operator indices
  for (int idx = 0; idx < rhs_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops;
    // Assume 2-index operators are available on this processor in core
    rhs_ops = rhs_array.get_local_element(idx);

    // Loop over all lhs operator indices
    for (int iidx = 0; iidx < lhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops = lhs_array.get_local_element(iidx);
      int i = rhs_ops[0]->get_orbs()[0]; int j = rhs_ops[0]->get_orbs()[1];
      int k = lhs_ops[0]->get_orbs()[0];
//pout << "i,j = " << i << "," << j << endl;
//pout << "k = " << k << endl;
      // In parallel calculations not all operators are built on each proc
      if ( ! big.get_op_array(optype).has_local_index(i,j,k) ) continue;
//pout << "building i,j,k = " << i << "," << j << "," << k << endl;
      std::vector<boost::shared_ptr<SparseMatrix> > vec = big.get_op_array(optype).get_element(i,j,k);
//pout << "got i,j,k\n";

      // Loop over rhs spin-op components
      for (int jdx=0; jdx < rhs_ops.size(); jdx++) {
        boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jdx];
        assert( rhs_op->get_built() );
//pout << "getting i,j build_pattern\n";
        std::string build_12 = rhs_op->get_build_pattern();
//pout << "getting i,j\n";
        std::vector<SpinQuantum> spin_12 = rhs_op->get_quantum_ladder().at(build_12);

        // Loop over lhs spin-op components //FIXME
        for (int jjdx=0; jjdx < lhs_ops.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jjdx];
          assert( lhs_op->get_built() );
          std::string build_3 = lhs_op->get_build_pattern();
          std::string build_pattern = "(" + build_12 + build_3 + ")";

          // Allocate and build new operator
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            // Select relevant spin component
            std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0) };
//            if ( s == spin_12 ) finish_tensor_product( big, rhsBlock, *rhs_op, Transposeview(lhs_op), *op, forwards, build_pattern );
            if ( s == spin_12 ) {
              finish_tensor_product( big, rhsBlock, *rhs_op, *lhs_op, *op, forwards, build_pattern );
              // Renormalise operator
              op->renormalise_transform( rotateMatrix, stateinfo );
            }
          }
        }
      }
      // Store spin-batch on disk 
      if ( ! dmrginp.do_npdm_in_core() ) store_ops_on_disk( ofs, vec );
    }
  }
  if ( rhsifs.is_open() ) rhsifs.close();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_2_1_tensor_products( bool forwards, const opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& big, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs,
                                    const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket )
{
  // (i,j | k ) partition
  //-------------------------
  Op_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
  assert ( (rhs_array.get_size() == 1) || (lhs_array.get_size() == 1) );

  // Initialize filesystem
  std::ifstream rhsifs;
  if ( (! dmrginp.do_npdm_in_core()) && rhsBlock->size() > 1 ) rhsifs.open( rhs_array.get_filename().c_str(), std::ios::binary );
//pout << "tensor_2_1: rhs_size = " << rhs_array.get_size() << "; op = " << rhs_array.get_op_string() << endl;
//pout << "tensor_2_1: lhs_size = " << lhs_array.get_size() << "; op = " << lhs_array.get_op_string() << endl;

  // Loop over all rhs operator indices
  for (int idx = 0; idx < rhs_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops;
    // Assume 2-index operators are available on this processor in core
    rhs_ops = rhs_array.get_local_element(idx);

    // Loop over all lhs operator indices
    for (int iidx = 0; iidx < lhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops = lhs_array.get_local_element(iidx);
      int i = rhs_ops[0]->get_orbs()[0]; int j = rhs_ops[0]->get_orbs()[1];
      int k = lhs_ops[0]->get_orbs()[0];
//pout << "i,j = " << i << "," << j << endl;
//pout << "k = " << k << endl;
      // In parallel calculations not all operators are built on each proc
      if ( ! big.get_op_array(optype).has_local_index(i,j,k) ) continue;
//pout << "building i,j,k = " << i << "," << j << "," << k << endl;
      std::vector<boost::shared_ptr<SparseMatrix> > vec = big.get_op_array(optype).get_element(i,j,k);
//pout << "got i,j,k\n";

      // Loop over rhs spin-op components
      for (int jdx=0; jdx < rhs_ops.size(); jdx++) {
        boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jdx];
        assert( rhs_op->get_built() );
//pout << "getting i,j build_pattern\n";
        std::string build_12 = rhs_op->get_build_pattern();
//pout << "getting i,j\n";
        std::vector<SpinQuantum> spin_12 = rhs_op->get_quantum_ladder().at(build_12);

        // Loop over lhs spin-op components //FIXME
        for (int jjdx=0; jjdx < lhs_ops.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jjdx];
          assert( lhs_op->get_built() );
          std::string build_3 = lhs_op->get_build_pattern();
          std::string build_pattern = "(" + build_12 + build_3 + ")";

          // Allocate and build new operator
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            // Select relevant spin component
            std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0) };
//            if ( s == spin_12 ) finish_tensor_product( big, rhsBlock, *rhs_op, Transposeview(lhs_op), *op, forwards, build_pattern );
            if ( s == spin_12 ) {
              finish_tensor_product( big, rhsBlock, *rhs_op, *lhs_op, *op, forwards, build_pattern );
              // Renormalise operator
              op->renormalise_transform( leftMat, bra, rightMat, ket);
            }
          }
        }
      }
      // Store spin-batch on disk 
      if ( ! dmrginp.do_npdm_in_core() ) store_ops_on_disk( ofs, vec );
    }
  }
  if ( rhsifs.is_open() ) rhsifs.close();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void build_3index_ops( const opTypes& optype, SpinBlock& big, 
                       const opTypes& lhsType1, const opTypes& lhsType2,
                       const opTypes& rhsType1, const opTypes& rhsType2,
                       const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo )
{
  // 3-index output file
//pout << "build_3index_op, ofs =" <<  big.get_op_array(optype).get_filename() << endl;
  std::ofstream ofs;
  if ( ! dmrginp.do_npdm_in_core() ) ofs.open( big.get_op_array(optype).get_filename().c_str(), std::ios::binary );

  SpinBlock* sysBlock = big.get_leftBlock();
  SpinBlock* dotBlock = big.get_rightBlock();

  // All 3 orbitals on sys or dot block
  do_3index_tensor_trace( optype, big, sysBlock, ofs, rotateMatrix, stateinfo );
  do_3index_tensor_trace( optype, big, dotBlock, ofs, rotateMatrix, stateinfo );

  bool forwards = ! ( sysBlock->get_sites().at(0) > dotBlock->get_sites().at(0) );

  // 2,1 partitioning
  if ( forwards ) {
    do_3index_1_2_tensor_products( forwards, optype, lhsType1, rhsType2, big, dotBlock, sysBlock, ofs, rotateMatrix, stateinfo );
    do_3index_2_1_tensor_products( forwards, optype, lhsType2, rhsType1, big, dotBlock, sysBlock, ofs, rotateMatrix, stateinfo );
  } else {
    do_3index_1_2_tensor_products( forwards, optype, lhsType1, rhsType2, big, sysBlock, dotBlock, ofs, rotateMatrix, stateinfo );
    do_3index_2_1_tensor_products( forwards, optype, lhsType2, rhsType1, big, sysBlock, dotBlock, ofs, rotateMatrix, stateinfo );
  }

  if ( ofs.is_open() ) ofs.close();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void build_3index_ops( const opTypes& optype, SpinBlock& big, 
                       const opTypes& lhsType1, const opTypes& lhsType2,
                       const opTypes& rhsType1, const opTypes& rhsType2,
                       const std::vector<Matrix>& leftMat,  const StateInfo *bra, 
                       const std::vector<Matrix>& rightMat, const StateInfo *ket )
{
  // 3-index output file
//pout << "build_3index_op, ofs =" <<  big.get_op_array(optype).get_filename() << endl;
  std::ofstream ofs;
  if ( ! dmrginp.do_npdm_in_core() ) ofs.open( big.get_op_array(optype).get_filename().c_str(), std::ios::binary );

  SpinBlock* sysBlock = big.get_leftBlock();
  SpinBlock* dotBlock = big.get_rightBlock();

  // All 3 orbitals on sys or dot block
  do_3index_tensor_trace( optype, big, sysBlock, ofs, leftMat, bra, rightMat, ket);
  do_3index_tensor_trace( optype, big, dotBlock, ofs, leftMat, bra, rightMat, ket);

  bool forwards = ! ( sysBlock->get_sites().at(0) > dotBlock->get_sites().at(0) );

  // 2,1 partitioning
  if ( forwards ) {
    do_3index_1_2_tensor_products( forwards, optype, lhsType1, rhsType2, big, dotBlock, sysBlock, ofs, leftMat, bra, rightMat, ket);
    do_3index_2_1_tensor_products( forwards, optype, lhsType2, rhsType1, big, dotBlock, sysBlock, ofs, leftMat, bra, rightMat, ket);
  } else {
    do_3index_1_2_tensor_products( forwards, optype, lhsType1, rhsType2, big, sysBlock, dotBlock, ofs, leftMat, bra, rightMat, ket);
    do_3index_2_1_tensor_products( forwards, optype, lhsType2, rhsType1, big, sysBlock, dotBlock, ofs, leftMat, bra, rightMat, ket);
  }

  if ( ofs.is_open() ) ofs.close();

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
}
}

