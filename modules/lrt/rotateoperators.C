 
#include "operatorfunctions.h"
#include "BaseOperator.h"
#include "rotateoperators.h"

// Written by N.N.; rotate operatos by Ritz vectors for DMRG-LRT
// currently, for < 0 | O | I > and < I | O | 0 > only, i.e. < 0 | O | I > = sum_{J} { < 0 | O | J > * alpha(J, I) }

// sort operators by state index and call rotate_operatos
void SpinAdapted::RotateOps::rotate_operator_arrays(const Matrix& alpha, SpinBlock& b, const opTypes& ot, int nroots, bool doTrans)
{
  const int mroots = alpha.Nrows()+1;

  // checking operators
  bool has_all = 1; // true if found all state indices
  bool has_one = 0; // true if found at least one state index
  for(int i = 1; i < mroots; ++i) {
    bool has_i;
    if(!doTrans)
      has_i = b.has(ot | make_state_index(0, i));
    else
      has_i = b.has(ot | make_state_index(i, 0));
    has_all &= has_i;
    has_one |= has_i;
  }
  if(!has_all && has_one) {
    pout << "\t\t\t ERROR @ RotateOps::rotate_operator_arrays: insufficient states were stored" << endl;
    abort();
  }
  else if(!has_one) {
    pout << "\t\t\t state rotation is not performed for dot or starting block" << endl;
    return;
  }

  const int opvec_size = b.get_op_array(ot).get_size(); // 0-th opvec size

  std::vector<Op_component_base*> opvec_array(mroots);

  // for ops[0I]
  for (int i = 0; i < mroots; ++i) {
    Op_component_base* opvec;
    if(!doTrans)
      opvec = &(b.get_op_array(ot | make_state_index(0, i)));
    else
      opvec = &(b.get_op_array(ot | make_state_index(i, 0)));
    assert(opvec->get_size() == opvec_size);
    opvec_array[i] = opvec;
  }

  // omp flags might be inserted here

  for (int m = 0; m < opvec_size; ++m) {
    std::vector<std::vector<boost::shared_ptr<SparseMatrix> > > ops(mroots);
    for (int i = 0; i < mroots; ++i) {
      ops[i] = opvec_array[i]->get_local_element(m);
    }
    rotate_operators(alpha, &b, ops, nroots);
  }

//pout << "DEBUG @ RotateOps::rotate_operator_arrays: check point 3" << endl;
//  // for ops[I0]
//  for (int i = 1; i < mroots; ++i) {
//pout << "DEBUG @ RotateOps::rotate_operator_arrays: check point 1 - " << i << endl;
//    Op_component_base* opvec = &(b.get_op_array(ot | make_state_index(i, 0)));
//    assert(opvec->get_size() == opvec_size);
//    opvec_array[i] = opvec;
//  }

//  // omp flags might be inserted here

//pout << "DEBUG @ RotateOps::rotate_operator_arrays: check point 4" << endl;
//  for (int m = 0; m < opvec_size; ++m) {
//    std::vector<std::vector<boost::shared_ptr<SparseMatrix> > > ops(mroots);
//    for (int i = 0; i < mroots; ++i) {
//      ops[i] = opvec_array[i]->get_local_element(m);
//    }
//    rotate_operators(alpha, &b, ops, nroots);
//  }
//pout << "DEBUG @ RotateOps::rotate_operator_arrays: return" << endl;
}

void SpinAdapted::RotateOps::rotate_operators(const Matrix& alpha, SpinBlock* b, std::vector<std::vector<boost::shared_ptr<SparseMatrix> > >& ops, int nroots)
{
  const int mroots = alpha.Nrows()+1;

  int opsize = ops[0].size();
  for (int opind = 0; opind < opsize; ++opind) {
    std::vector<boost::shared_ptr<SparseMatrix> > opx     (mroots);
    std::vector<boost::shared_ptr<SparseMatrix> > opx_save(mroots);
    for (int i = 1; i < mroots; ++i) {
      opx[i]      = ops[i].at(opind)->getworkingrepresentation(b);
      opx_save[i] = opx[i]->deepCopy(); // returns deep copied shared_ptr
    }
    for(int i = 1; i < nroots; ++i) {
      Scale(alpha(i, i), *opx[i]);
    }
    for(int i = 1; i < nroots; ++i) {
      for(int j = 1; j < mroots; ++j) {
        if(i != j) {
          ScaleAdd(alpha(j, i), *opx_save[j], *opx[i]);
        }
      }
    }
  }
}

