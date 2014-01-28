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

//
//FIXME clean up these routines to make shorter and simpler and notso much duplicate code!!!
//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_tensor_trace( opTypes& ot, SpinBlock& b, std::ofstream& ofs ) {

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  // Input files for sys and dot 4-index ops
  std::ifstream sys_ifs(sysBlock->get_op_array(ot).get_filename().c_str(), std::ios::binary);
  std::ifstream dot_ifs(dotBlock->get_op_array(ot).get_filename().c_str(), std::ios::binary);

cout << "doing CCDD tensor trace with dot block\n";
  // Tensor trace with dot block
  Op_component_base& dot_array = dotBlock->get_op_array(CRE_CRE_DES_DES);
  // Loop over all operator indices
  assert( dot_array.get_size() == 1 );
  for (int idx = 0; idx < dot_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array.get_local_element(idx);
    // Loop over spin-op components
    for (int jdx=0; jdx < dot_vec.size(); jdx++) {
      boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jdx];
      assert( dot_op->get_built() );
      int i = dot_op->get_orbs()[0]; int j = dot_op->get_orbs()[1]; int k = dot_op->get_orbs()[2]; int l = dot_op->get_orbs()[3];
      assert( (i==j) && (j==k) && (k==l) );
      std::string build_pattern = dot_op->get_build_pattern();
//FIXME are there any repetitions???
      // Allocate and build new operator
      std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
      for (int sx=0; sx < vec.size(); sx++) {
        boost::shared_ptr<SparseMatrix>& op = vec[sx];
        std::vector<SpinQuantum> s1 = dot_op->get_quantum_ladder().at(build_pattern);
        std::vector<SpinQuantum> s2 = op->get_quantum_ladder().at(build_pattern);
        if ( s1 == s2 ) {
cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
          assert( ! op->get_built() );
          op->set_built() = true;
          op->set_build_pattern() = build_pattern;
          op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
          op->allocate(b.get_stateInfo());
          SpinAdapted::operatorfunctions::TensorTrace(dotBlock, *dot_op, &b, &(b.get_stateInfo()), *op);
//FIXME
          // Store on disk
          op->set_built_on_disk() = true;
          boost::archive::binary_oarchive save_op(ofs);
          save_op << *op;
//          // Deallocate memory for operator representation
//          op->set_built() = false;
//          op->deallocate(b);
        }
      }
    }
  }

cout << "doing CCDD tensor trace with sys block\n";
  // Tensor trace with sys block
  Op_component_base& sys_array = sysBlock->get_op_array(CRE_CRE_DES_DES);
  // Loop over all operator indices
  for (int idx = 0; idx < sys_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array.get_local_element(idx);
    // Loop over spin-op components
    for (int jdx=0; jdx < sys_vec.size(); jdx++) {
      boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jdx];
      assert( sys_op->get_built() );
      int i = sys_op->get_orbs()[0]; int j = sys_op->get_orbs()[1]; int k = sys_op->get_orbs()[2]; int l = sys_op->get_orbs()[3];
      std::string build_pattern = sys_op->get_build_pattern();
      // Allocate and build new operator
      std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
      for (int sx=0; sx < vec.size(); sx++) {
        boost::shared_ptr<SparseMatrix>& op = vec[sx];
        std::vector<SpinQuantum> s1 = sys_op->get_quantum_ladder().at(build_pattern);
        std::vector<SpinQuantum> s2 = op->get_quantum_ladder().at(build_pattern);
        if ( s1 == s2 ) {
cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
          assert( ! op->get_built() );
          op->set_built() = true;
          op->set_build_pattern() = build_pattern;
          op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
          op->allocate(b.get_stateInfo());
          SpinAdapted::operatorfunctions::TensorTrace(sysBlock, *sys_op, &b, &(b.get_stateInfo()), *op);
//FIXME
          // Store on disk
          op->set_built_on_disk() = true;
          boost::archive::binary_oarchive save_op(ofs);
          save_op << *op;
//          // Deallocate memory for operator representation
//          op->set_built() = false;
//          op->deallocate(b);
        }
      }
    }
  }
cout << "done all tensor trace ops\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_3_1_tensor_products( opTypes& ot, SpinBlock& b, std::ofstream& ofs ) {

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  bool forward = ! ( sysBlock->get_sites()[0] > dotBlock->get_sites()[0] );

  // For backwards sweep
  //-------------------------
  if ( ! forward ) {
cout << "doing CCDD tensor product with (i) on sys; (j,k,l) on dot\n";
//FIXME note order of these loops

    // Loop over all sys operator indices
    Op_component_base& sys_array = sysBlock->get_op_array(CRE);
    for (int iidx = 0; iidx < sys_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array.get_local_element(iidx);

      // Loop over all dot operator indices
      Op_component_base& dot_array = dotBlock->get_op_array(CRE_DES_DES);
      assert( dot_array.get_size() == 1 );
      for (int idx = 0; idx < dot_array.get_size(); ++idx) {
        std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array.get_local_element(idx);

        // Loop over sys spin-op components
        assert( sys_vec.size() == 1 );
        for (int jjdx=0; jjdx < sys_vec.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jjdx];
          assert( sys_op->get_built() );
          int i = sys_op->get_orbs()[0];

          // Loop over dot spin-op components
          for (int jdx=0; jdx < dot_vec.size(); jdx++) {
            boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jdx];
            assert( dot_op->get_built() );
            int j = dot_op->get_orbs()[0]; int k = dot_op->get_orbs()[1]; int l = dot_op->get_orbs()[2];
            std::string build_234 = dot_op->get_build_pattern();
            std::string build_pattern = "((C)" + build_234 + ")";
            std::vector<SpinQuantum> spin_234 = dot_op->get_quantum_ladder().at(build_234);

            // Allocate and build new operator
            std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
            for (int sx=0; sx < vec.size(); sx++) {
              boost::shared_ptr<SparseMatrix>& op = vec[sx];
              //FIXME Each 3-index spin component should contribute to two 4-index spin_ops
              // Select relevant spin component
              std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
              if ( s == spin_234 ) {
  cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
                assert( ! op->get_built() );
                op->set_built() = true;
                op->set_build_pattern() = build_pattern;
                op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
                op->allocate(b.get_stateInfo());
                SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *sys_op, *dot_op, &b, &(b.get_stateInfo()), *op, 1.0);
//FIXME
                // Store on disk
                op->set_built_on_disk() = true;
                boost::archive::binary_oarchive save_op(ofs);
                save_op << *op;
//cout << "op:\n";
//cout << *op << endl;
      //          // Deallocate memory for operator representation
      //          op->set_built() = false;
      //          op->deallocate(b);
              }
            }
          }
        }
      }
    }

cout << "doing CCDD tensor product with (i,j,k) on sys; (l) on dot\n";
    // Loop over all sys operator indices
    Op_component_base& sys_array2 = sysBlock->get_op_array(CRE_CRE_DES);
    for (int idx = 0; idx < sys_array2.get_size(); ++idx) {
      std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array2.get_local_element(idx);

      // Loop over all dot operator indices
      Op_component_base& dot_array = dotBlock->get_op_array(CRE);
      assert( dot_array.get_size() == 1 );
      for (int iidx = 0; iidx < dot_array.get_size(); ++iidx) {
        std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array.get_local_element(iidx);

        // Loop over sys spin-op components
        for (int jdx=0; jdx < sys_vec.size(); jdx++) {
          boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jdx];
          assert( sys_op->get_built() );
          int i = sys_op->get_orbs()[0]; int j = sys_op->get_orbs()[1]; int k = sys_op->get_orbs()[2];
          std::string build_123 = sys_op->get_build_pattern();
          std::string build_pattern = "(" + build_123 + "(D))";
          std::vector<SpinQuantum> spin_123 = sys_op->get_quantum_ladder().at(build_123);

          // Loop over dot spin-op components //FIXME
          assert( dot_vec.size() == 1 );
          for (int jjdx=0; jjdx < dot_vec.size(); jjdx++) {
            boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jjdx];
            assert( dot_op->get_built() );
            int l = dot_op->get_orbs()[0];

            // Allocate and build new operator
            std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
            for (int sx=0; sx < vec.size(); sx++) {
              boost::shared_ptr<SparseMatrix>& op = vec[sx];
              // Select relevant spin component
              std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
              if ( s == spin_123 ) {
cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
                assert( ! op->get_built() );
                op->set_built() = true;
                op->set_build_pattern() = build_pattern;
                op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
                op->allocate(b.get_stateInfo());
                SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *sys_op, Transposeview(dot_op), &b, &(b.get_stateInfo()), *op, 1.0);
//FIXME
                // Store on disk
                op->set_built_on_disk() = true;
                boost::archive::binary_oarchive save_op(ofs);
                save_op << *op;
//cout << "op:\n";
//cout << *op << endl;
      //          // Deallocate memory for operator representation
      //          op->set_built() = false;
      //          op->deallocate(b);
              }
            }
          }
        }
      }
    }
  }

  // For forwards sweep
  //-------------------------
  else {
cout << "doing CCDD tensor product with (i,j,k) on dot; (l) on sys\n";

    // Loop over all dot operator indices
    Op_component_base& dot_array = dotBlock->get_op_array(CRE_CRE_DES);
    assert( dot_array.get_size() == 1 );
    for (int idx = 0; idx < dot_array.get_size(); ++idx) {
      std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array.get_local_element(idx);

      // Loop over all sys operator indices
      Op_component_base& sys_array = sysBlock->get_op_array(CRE);
      for (int iidx = 0; iidx < sys_array.get_size(); ++iidx) {
        std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array.get_local_element(iidx);

        // Loop over dot spin-op components
        for (int jdx=0; jdx < dot_vec.size(); jdx++) {
          boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jdx];
          assert( dot_op->get_built() );
          int i = dot_op->get_orbs()[0]; int j = dot_op->get_orbs()[1]; int k = dot_op->get_orbs()[2];
          std::string build_123 = dot_op->get_build_pattern();
          std::string build_pattern = "(" + build_123 + "(D))";
          std::vector<SpinQuantum> spin_123 = dot_op->get_quantum_ladder().at(build_123);
  
          // Loop over sys spin-op components
          assert( sys_vec.size() == 1 );
          for (int jjdx=0; jjdx < sys_vec.size(); jjdx++) {
            boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jjdx];
            assert( sys_op->get_built() );
            int l = sys_op->get_orbs()[0];

            // Allocate and build new operator
            std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
            for (int sx=0; sx < vec.size(); sx++) {
              boost::shared_ptr<SparseMatrix>& op = vec[sx];
              // Select relevant spin component
              std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
              if ( s == spin_123 ) {
cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
                assert( ! op->get_built() );
                op->set_built() = true;
                op->set_build_pattern() = build_pattern;
                op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
                op->allocate(b.get_stateInfo());
                // Tensor product of dot*sys so need to take into account parity factors
//FIXME minus sign in sys_op here for DES instead of CRE...
                double parity = getCommuteParity( dot_op->get_deltaQuantum(), -sys_op->get_deltaQuantum(), op->get_deltaQuantum() );
                SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *dot_op, Transposeview(sys_op), &b, &(b.get_stateInfo()), *op, 1.0*parity);
//FIXME
                // Store on disk
                op->set_built_on_disk() = true;
                boost::archive::binary_oarchive save_op(ofs);
                save_op << *op;
//cout << "op:\n";
//cout << *op << endl;
      //          // Deallocate memory for operator representation
      //          op->set_built() = false;
      //          op->deallocate(b);
              }
            }
          }
        }
      }
    }
cout << "doing CCDD tensor product with (i) on dot; (j,k,l) on sys\n";

    // Loop over all dot operator indices
    Op_component_base& dot_array2 = dotBlock->get_op_array(CRE);
    assert( dot_array2.get_size() == 1 );
    for (int iidx = 0; iidx < dot_array2.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array2.get_local_element(iidx);

      // Loop over all sys operator indices
      Op_component_base& sys_array = sysBlock->get_op_array(CRE_DES_DES);
      for (int idx = 0; idx < sys_array.get_size(); ++idx) {
        std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array.get_local_element(idx);

        // Loop over dot spin-op components //FIXME
        assert( dot_vec.size() == 1 );
        for (int jjdx=0; jjdx < dot_vec.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jjdx];
          assert( dot_op->get_built() );
          int i = dot_op->get_orbs()[0];

          // Loop over sys spin-op components
          for (int jdx=0; jdx < sys_vec.size(); jdx++) {
            boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jdx];
            assert( sys_op->get_built() );
            int j = sys_op->get_orbs()[0]; int k = sys_op->get_orbs()[1]; int l = sys_op->get_orbs()[2];
            std::string build_234 = sys_op->get_build_pattern();
            std::string build_pattern = "((C)" + build_234 + ")";
            std::vector<SpinQuantum> spin_234 = sys_op->get_quantum_ladder().at(build_234);

            // Allocate and build new operator
            std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
            for (int sx=0; sx < vec.size(); sx++) {
              boost::shared_ptr<SparseMatrix>& op = vec[sx];
              // Select relevant spin component
              std::vector<SpinQuantum> s = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
              if ( s == spin_234 ) {
cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
                assert( ! op->get_built() );
                op->set_built() = true;
                op->set_build_pattern() = build_pattern;
                op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
                op->allocate(b.get_stateInfo());
                // Tensor product of dot*sys so need to take into account parity factors
                double parity = getCommuteParity( dot_op->get_deltaQuantum(), sys_op->get_deltaQuantum(), op->get_deltaQuantum() );
                SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *dot_op, *sys_op, &b, &(b.get_stateInfo()), *op, 1.0*parity);
//cout << "op:\n";
//cout << *op << endl;
//FIXME
                // Store on disk
                op->set_built_on_disk() = true;
                boost::archive::binary_oarchive save_op(ofs);
                save_op << *op;
      //          // Deallocate memory for operator representation
      //          op->set_built() = false;
      //          op->deallocate(b);
              }
            }
          }
        }
      }
    }
  }
cout << "done all 3,1 tensor product ops\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_4index_2_2_tensor_products( opTypes& ot, SpinBlock& b, std::ofstream& ofs ) {

  SpinBlock* sysBlock = b.get_leftBlock();
  SpinBlock* dotBlock = b.get_rightBlock();

  bool forward = ! ( sysBlock->get_sites()[0] > dotBlock->get_sites()[0] );
  std::string build_pattern = "((CC)(DD))";

  // For backwards sweep
  //-------------------------
  if ( ! forward ) {
cout << "doing CCDD tensor product with (i,j) on sys; (k,l) on dot\n";

    // Loop over all sys operator indices
    Op_component_base& sys_array = sysBlock->get_op_array(CRE_CRE);
    for (int iidx = 0; iidx < sys_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array.get_local_element(iidx);

      // Loop over all dot operator indices
      Op_component_base& dot_array = dotBlock->get_op_array(DES_DES);
      assert( dot_array.get_size() == 1 );
      for (int idx = 0; idx < dot_array.get_size(); ++idx) {
        std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array.get_local_element(idx);

        // Loop over sys spin-op components //FIXME
        for (int jjdx=0; jjdx < sys_vec.size(); jjdx++) {
          boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jjdx];
          int i = sys_op->get_orbs()[0]; int j = sys_op->get_orbs()[1];
          const SpinQuantum& spin_12 = sys_op->get_deltaQuantum();

          // Loop over dot spin-op components
          for (int jdx=0; jdx < dot_vec.size(); jdx++) {
            boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jdx];
            int k = dot_op->get_orbs()[0]; int l = dot_op->get_orbs()[1];
            const SpinQuantum& spin_34 = dot_op->get_deltaQuantum();

            // Allocate and build new operator
            std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
            for (int sx=0; sx < vec.size(); sx++) {
              boost::shared_ptr<SparseMatrix>& op = vec[sx];
              // Select relevant spin component
              std::vector<SpinQuantum> s1 = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
              std::vector<SpinQuantum> s2 = { spin_12, spin_34 };
              if ( s1 == s2 ) {
  cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
                assert( ! op->get_built() );
                op->set_built() = true;
                op->set_build_pattern() = build_pattern;
                op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
                op->allocate(b.get_stateInfo());
                SpinAdapted::operatorfunctions::TensorProduct(sysBlock, *sys_op, *dot_op, &b, &(b.get_stateInfo()), *op, 1.0);
//FIXME
                // Store on disk
                op->set_built_on_disk() = true;
                boost::archive::binary_oarchive save_op(ofs);
                save_op << *op;
      //          // Deallocate memory for operator representation
      //          op->set_built() = false;
      //          op->deallocate(b);
              }
            }
          }
        }
      }
    }
  }
  // For forwards sweep
  //-------------------------
  else {
cout << "doing CCDD tensor product with (i,j) on dot; (k,l) on sys\n";

    // Loop over all dot operator indices
    Op_component_base& dot_array = dotBlock->get_op_array(CRE_CRE);
    assert( dot_array.get_size() == 1 );
    for (int idx = 0; idx < dot_array.get_size(); ++idx) {
      std::vector<boost::shared_ptr<SparseMatrix> > dot_vec = dot_array.get_local_element(idx);

      // Loop over all sys operator indices
      Op_component_base& sys_array = sysBlock->get_op_array(DES_DES);
      for (int iidx = 0; iidx < sys_array.get_size(); ++iidx) {
        std::vector<boost::shared_ptr<SparseMatrix> > sys_vec = sys_array.get_local_element(iidx);

        // Loop over dot spin-op components
        for (int jdx=0; jdx < dot_vec.size(); jdx++) {
          boost::shared_ptr<SparseMatrix>& dot_op = dot_vec[jdx];
          int i = dot_op->get_orbs()[0]; int j = dot_op->get_orbs()[1];
          const SpinQuantum& spin_12 = dot_op->get_deltaQuantum();

          // Loop over sys spin-op components //FIXME
          for (int jjdx=0; jjdx < sys_vec.size(); jjdx++) {
            boost::shared_ptr<SparseMatrix>& sys_op = sys_vec[jjdx];
            int k = sys_op->get_orbs()[0]; int l = sys_op->get_orbs()[1];
            const SpinQuantum& spin_34 = sys_op->get_deltaQuantum();

            // Allocate and build new operator
            std::vector<boost::shared_ptr<SparseMatrix> > vec = b.get_op_array(CRE_CRE_DES_DES).get_element(i,j,k,l);
            for (int sx=0; sx < vec.size(); sx++) {
              boost::shared_ptr<SparseMatrix>& op = vec[sx];
              std::vector<SpinQuantum> s1 = { op->get_quantum_ladder().at(build_pattern).at(0), op->get_quantum_ladder().at(build_pattern).at(1) };
              std::vector<SpinQuantum> s2 = { spin_12, spin_34 };
              if ( s1 == s2 ) {
  cout << "indices = " << i << "," << j << "," << k << "," << l << endl;
                assert( ! op->get_built() );
                op->set_built() = true;
                op->set_build_pattern() = build_pattern;
                op->set_deltaQuantum() = op->get_quantum_ladder().at( build_pattern ).at(2);
                op->allocate(b.get_stateInfo());
                // FIXME is parity even needed???
                // Tensor product of dot*sys so need to take into account parity factors
                double parity = getCommuteParity( dot_op->get_deltaQuantum(), sys_op->get_deltaQuantum(), op->get_deltaQuantum() );
                SpinAdapted::operatorfunctions::TensorProduct(dotBlock, *dot_op, *sys_op, &b, &(b.get_stateInfo()), *op, 1.0*parity);
//FIXME
                // Store on disk
                op->set_built_on_disk() = true;
                boost::archive::binary_oarchive save_op(ofs);
                save_op << *op;
      //          // Deallocate memory for operator representation
      //          op->set_built() = false;
      //          op->deallocate(b);
              }
            }
          }
        }
      }
    }
  }
cout << "done all 2,2 tensor product ops\n";

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//      if(! it->second->is_core()) {
//        // Output file for operators written to disk
//        std::string ofile = it->second->get_filename();
//        // Input file for operators written to disk on sysblock
//        std::string sysfile = get_leftBlock()->ops[ot]->get_filename();
//        // Input file for operators written to disk on dotblock
//        std::string dotfile = get_rightBlock()->ops[ot]->get_filename();
//        // Build operators
//
//        if ( ot == CRE_CRE_DES_DES ) {
//
//      std::ofstream ofs(ofile.c_str(), std::ios::binary);
//      std::ifstream sysfs(sysfile.c_str(), std::ios::binary);
//      std::ifstream dotfs(dotfile.c_str(), std::ios::binary);
//
//
//  bool check_file_open( int idx )
//    {
//      if ( idx == 0 ) {
//        assert( ! ifs_.is_open() );
//        ifs_.open(ifile_.c_str(), ios::binary);
//      }
//      return ifs_.good();
//    }
//    bool check_file_close( int idx ) 
//    { 
//      if ( idx == (size_-1) ) {
//        assert( ifs_.is_open() );
//        ifs_.close(); 
//      }
//      return true;
//    }



//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void build_4index_ops( opTypes& ot, SpinBlock& b )
{
  // Hard code for CCDD to start with
  assert( ot == CRE_CRE_DES_DES );

  // 4-index output file
  std::ofstream ofs(b.get_op_array(ot).get_filename().c_str(), std::ios::binary);

  do_4index_tensor_trace( ot, b, ofs );
  do_4index_3_1_tensor_products( ot, b, ofs );
  do_4index_2_2_tensor_products( ot, b, ofs );

  ofs.close();

cout << "CCDD all done!\n";

}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}

