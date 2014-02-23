/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


//#include <vector>
//#include <iostream>
//#include <communicate.h>
#include "op_components.h"
#include "BaseOperator.h"
#include "spinblock.h"
#include "operatorfunctions.h"
//#include "tensor_operator.h"
//#include <boost/shared_ptr.hpp>

#include "build_4index_ops.h"

namespace SpinAdapted{

//cout << "indices = " << i << "," << j << "," << k << "," << l << endl;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

std::vector<boost::shared_ptr<SparseMatrix> > get_ops_from_disk( int idx, std::ifstream& ifs, Op_component_base& sysdot_array )
{
//FIXME reference or copy??
  std::vector<boost::shared_ptr<SparseMatrix> > opReps;

  // Read in first spin component
//FIXME optype
  boost::shared_ptr<SparseMatrix> op (new Cre);
  boost::archive::binary_iarchive load_op(ifs);
  load_op >> *op;
  std::string build_pattern = op->get_build_pattern();
  // Get order of spin components
  int ix = op->get_orbs(0); 
  int jx = op->get_orbs(1); 
  int kx = op->get_orbs(2); 
  int lx = op->get_orbs(3);
  opReps = sysdot_array.get_element(ix,jx,kx,lx);
  // Store spin components in correct order
  for (int j = 0; j < opReps.size(); j++) {
    if ( opReps[j]->get_quantum_ladder().at(build_pattern) == op->get_quantum_ladder().at(build_pattern) ) {
      opReps[j] = op; 
      break;
    }
  }

  // Read in remaining spin components
  for (int i = 1; i < opReps.size(); i++) {
//FIXME optype
    boost::shared_ptr<SparseMatrix> op (new Cre);
    boost::archive::binary_iarchive load_op(ifs);
    load_op >> *op;
    for (int j = 0; j < opReps.size(); j++) {
      // Store spin components in correct order
      if ( opReps[j]->get_quantum_ladder().at(build_pattern) == op->get_quantum_ladder().at(build_pattern) ) {
        opReps[j] = op; 
        break;
      }
    }
  }

  return opReps;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void finish_tensor_trace( SpinBlock& b, SpinBlock* sysdot, SparseMatrix& sysdot_op, SparseMatrix& op, std::string& build_pattern, std::ofstream& ofs ) 
{
  // Build and store new operator
  assert( ! op.get_built() );
  op.set_built() = true;
  op.set_build_pattern() = build_pattern;
  op.set_deltaQuantum() = op.get_quantum_ladder().at( build_pattern ).at(2);
  op.allocate(b.get_stateInfo());
  SpinAdapted::operatorfunctions::TensorTrace(sysdot, sysdot_op, &b, &(b.get_stateInfo()), op);

  // Move to disk if necessary
//FIXME
//          if ( ! dmrginp.do_npdm_in_core() ) {
  if ( true ) {
    op.set_built_on_disk() = true;
    boost::archive::binary_oarchive save_op(ofs);
    save_op << op;
    // Deallocate memory for operator representation
    op.set_built() = false;
    op.deallocate(b);
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void finish_tensor_product( SpinBlock& b, SpinBlock* sysdot, 
                            const SparseMatrix& sysdot_op1, const SparseMatrix& sysdot_op2, SparseMatrix& op, 
                            bool include_parity, std::string& build_pattern, std::ofstream& ofs ) 
{
  // Build and store new operator
  assert( ! op.get_built() );
  op.set_built() = true;
  op.set_build_pattern() = build_pattern;
  op.set_deltaQuantum() = op.get_quantum_ladder().at( build_pattern ).at(2);
  op.allocate(b.get_stateInfo());

  // Do tensor product
  double parity = 1.0;
  if ( include_parity ) parity = getCommuteParity( sysdot_op1.get_deltaQuantum(), sysdot_op2.get_deltaQuantum(), op.get_deltaQuantum() );
  SpinAdapted::operatorfunctions::TensorProduct(sysdot, sysdot_op1, sysdot_op2, &b, &(b.get_stateInfo()), op, parity);

  // Move to disk if necessary
//FIXME
  if ( true ) {
    op.set_built_on_disk() = true;
    boost::archive::binary_oarchive save_op(ofs);
    save_op << op;
    // Deallocate memory for operator representation
    op.set_built() = false;
    op.deallocate(b);
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_tensor_trace( opTypes& optype, SpinBlock& b, SpinBlock* sysdot, std::ofstream& ofs ) 
{

  // Get pointer to sparse operator array
  Op_component_base& sysdot_array = sysdot->get_op_array(optype);
//FIXME
  // Open filesystem if necessary
  std::ifstream ifs( sysdot_array.get_filename().c_str(), std::ios::binary );
  if ( true ) {
//  if ( ! dmrginp.do_npdm_in_core() ) {
    ifs.open( sysdot_array.get_filename().c_str(), std::ios::binary );
  }

//FIXME need reference?  don't want to copy?
  std::vector<boost::shared_ptr<SparseMatrix> > sysdot_ops;
  // Loop over all operator indices
  for (int idx = 0; idx < sysdot_array.get_size(); ++idx) {
//    if ( dmrginp.do_npdm_in_core() ) {
    if ( false )
      sysdot_ops = sysdot_array.get_local_element(idx);
    else
      sysdot_ops = get_ops_from_disk( idx, ifs, sysdot_array );

    // Loop over spin-op components
    for (int jdx=0; jdx < sysdot_ops.size(); jdx++) {
      boost::shared_ptr<SparseMatrix>& sysdot_op = sysdot_ops[jdx];
      assert( sysdot_op->get_built() );
      int i = sysdot_op->get_orbs()[0]; int j = sysdot_op->get_orbs()[1]; int k = sysdot_op->get_orbs()[2]; int l = sysdot_op->get_orbs()[3];
      std::string build_pattern = sysdot_op->get_build_pattern();

      // Allocate and build new operator
      std::vector<boost::shared_ptr<SparseMatrix> > new_ops = b.get_op_array(optype).get_element(i,j,k,l);
      for (int sx=0; sx < new_ops.size(); sx++) {
        boost::shared_ptr<SparseMatrix>& op = new_ops[sx];
        std::vector<SpinQuantum> s1 = sysdot_op->get_quantum_ladder().at(build_pattern);
        std::vector<SpinQuantum> s2 = op->get_quantum_ladder().at(build_pattern);
        // Store spin component in correct location
        if ( s1 == s2 ) finish_tensor_trace( b, sysdot, *sysdot_op, *op, build_pattern, ofs );
      }
    }
  }

  // Close filesystem if necessary
//FIXME
  if ( ifs.is_open() ) {
    ifs.close();
//  assert( ! dmrginp.do_npdm_in_core() );
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_2_2_tensor_products( bool forwards, opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& b, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs ) 
{
  // Loop over all rhs operator indices
  Op_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  for (int iidx = 0; iidx < rhs_array.get_size(); ++iidx) {
    std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops = rhs_array.get_local_element(iidx);

    // Loop over all lhs operator indices
    Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
    for (int idx = 0; idx < lhs_array.get_size(); ++idx) {
      std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops = lhs_array.get_local_element(idx);

      // Loop over rhs spin-op components //FIXME
      for (int jjdx=0; jjdx < rhs_ops.size(); jjdx++) {
        boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jjdx];
        int i = rhs_op->get_orbs()[0]; int j = rhs_op->get_orbs()[1];
        const SpinQuantum& spin_12 = rhs_op->get_deltaQuantum();

        // Loop over lhs spin-op components
        for (int jdx=0; jdx < lhs_ops.size(); jdx++) {
          boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jdx];
          int k = lhs_op->get_orbs()[0]; int l = lhs_op->get_orbs()[1];
          std::string build_12 = rhs_op->get_build_pattern();
          std::string build_34 = lhs_op->get_build_pattern();
          std::string build_pattern = "(" + build_12 + build_34 + ")";
          const SpinQuantum& spin_34 = lhs_op->get_deltaQuantum();

          // Allocate and build new operator
          std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(optype).get_element(i,j,k,l);
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            std::vector<SpinQuantum> s1 = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
            std::vector<SpinQuantum> s2 = { spin_12, spin_34 };
            // Select relevant spin component
            if ( s1 == s2 ) finish_tensor_product( b, rhsBlock, *rhs_op, *lhs_op, *op, forwards, build_pattern, ofs );
          }

        }
      }
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_1_3_tensor_products( bool forwards, opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& b, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs ) 
{
  // (i | j,k,l ) partition
  //-------------------------
  // Loop over all rhs operator indices
  Op_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  for (int iidx = 0; iidx < rhs_array.get_size(); ++iidx) {
    std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops = rhs_array.get_local_element(iidx);

    // Loop over all lhs operator indices
    Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
    for (int idx = 0; idx < lhs_array.get_size(); ++idx) {
      std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops = lhs_array.get_local_element(idx);

      // Loop over rhs spin-op components
      for (int jjdx=0; jjdx < rhs_ops.size(); jjdx++) {
        boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jjdx];
        assert( rhs_op->get_built() );
        int i = rhs_op->get_orbs()[0];

        // Loop over lhs spin-op components
        for (int jdx=0; jdx < lhs_ops.size(); jdx++) {
          boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jdx];
          assert( lhs_op->get_built() );
          int j = lhs_op->get_orbs()[0]; int k = lhs_op->get_orbs()[1]; int l = lhs_op->get_orbs()[2];
          std::string build_1 = rhs_op->get_build_pattern();
          std::string build_234 = lhs_op->get_build_pattern();
          std::string build_pattern = "(" + build_1 + build_234 + ")";
          std::vector<SpinQuantum> spin_234 = lhs_op->get_quantum_ladder().at(build_234);

          // Allocate and build new operator
          std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(optype).get_element(i,j,k,l);
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            //FIXME Each 3-index spin component should contribute to two 4-index spin_ops
            // Select relevant spin component
            std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
            if ( s == spin_234 ) finish_tensor_product( b, rhsBlock, *rhs_op, *lhs_op, *op, forwards, build_pattern, ofs );
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_3_1_tensor_products( bool forwards, opTypes& optype, const opTypes& rhsType, const opTypes& lhsType,
                                    SpinBlock& b, SpinBlock* rhsBlock, SpinBlock* lhsBlock, std::ofstream& ofs ) 
{
  // (i,j,k | l ) partition
  //-------------------------
  // Loop over all rhs operator indices
  Op_component_base& rhs_array2 = rhsBlock->get_op_array(rhsType);
  for (int idx = 0; idx < rhs_array2.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > rhs_ops = rhs_array2.get_local_element(idx);

    // Loop over all lhs operator indices
    Op_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
    for (int iidx = 0; iidx < lhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > lhs_ops = lhs_array.get_local_element(iidx);

      // Loop over rhs spin-op components
      for (int jdx=0; jdx < rhs_ops.size(); jdx++) {
        boost::shared_ptr<SparseMatrix>& rhs_op = rhs_ops[jdx];
        assert( rhs_op->get_built() );
        int i = rhs_op->get_orbs()[0]; int j = rhs_op->get_orbs()[1]; int k = rhs_op->get_orbs()[2];
        std::string build_123 = rhs_op->get_build_pattern();
//FIXME only works with "D" as 4th operator!!
assert( optype == CRE_CRE_DES_DES );
        std::string build_pattern = "(" + build_123 + "(D))";
        std::vector<SpinQuantum> spin_123 = rhs_op->get_quantum_ladder().at(build_123);

        // Loop over lhs spin-op components //FIXME
        for (int jjdx=0; jjdx < lhs_ops.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& lhs_op = lhs_ops[jjdx];
          assert( lhs_op->get_built() );
          int l = lhs_op->get_orbs()[0];

          // Allocate and build new operator
          std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(optype).get_element(i,j,k,l);
          for (int sx=0; sx < vec.size(); sx++) {
            boost::shared_ptr<SparseMatrix>& op = vec[sx];
            // Select relevant spin component
            std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
//FIXME Transposeview assumed for lhs_op
            if ( s == spin_123 ) finish_tensor_product( b, rhsBlock, *rhs_op, Transposeview(lhs_op), *op, forwards, build_pattern, ofs );
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void build_4index_ops( opTypes& optype, SpinBlock& big )
{
  // 4-index output file
  std::ofstream ofs(big.get_op_array(optype).get_filename().c_str(), std::ios::binary);

  SpinBlock* sysBlock = big.get_leftBlock();
  SpinBlock* dotBlock = big.get_rightBlock();

  // All 4 orbitals on sys or dot block
  do_4index_tensor_trace( optype, big, sysBlock, ofs );
  do_4index_tensor_trace( optype, big, dotBlock, ofs );

  // 2,2 partitioning
  bool forwards = ! ( sysBlock->get_sites().at(0) > dotBlock->get_sites().at(0) );
  if ( forwards ) {
    do_4index_2_2_tensor_products( forwards, optype, CRE_CRE, DES_DES, big, dotBlock, sysBlock, ofs );
  } else {
    do_4index_2_2_tensor_products( forwards, optype, CRE_CRE, DES_DES, big, sysBlock, dotBlock, ofs );
  }

  // 3,1 partitioning
  if ( forwards ) {
    do_4index_1_3_tensor_products( forwards, optype, CRE, CRE_DES_DES, big, dotBlock, sysBlock, ofs );
    do_4index_3_1_tensor_products( forwards, optype, CRE_CRE_DES, CRE, big, dotBlock, sysBlock, ofs );
  } else {
    do_4index_1_3_tensor_products( forwards, optype, CRE, CRE_DES_DES, big, sysBlock, dotBlock, ofs );
    do_4index_3_1_tensor_products( forwards, optype, CRE_CRE_DES, CRE, big, sysBlock, dotBlock, ofs );
  }

  ofs.close();


}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}

