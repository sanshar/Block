/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "spinblock.h"

namespace SpinAdapted{
void SpinBlock::setstoragetype(Storagetype st)
{
  if (st == LOCAL_STORAGE)
  {
    localstorage = true;
//  if (has(CRE))
//    set_op_array(CRE).set_local() = true;
//  if (has(CRE_DES))
//    set_op_array(CRE_DES).set_local() = true;
//  if (has(CRE_CRE))
//    set_op_array(CRE_CRE).set_local() = true;
//  if (has(DES_DESCOMP))
//    set_op_array(DES_DESCOMP).set_local() = true;
//  if (has(CRE_DESCOMP))
//    set_op_array(CRE_DESCOMP).set_local() = true;
//  if (has(CRE_CRE_DESCOMP))
//    set_op_array(CRE_CRE_DESCOMP).set_local() = true;
    for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
      if((it->first & OP_TYPE_MASK) == HAM) continue;
      it->second->set_local() = true;
//    set_op_array(it->first).set_local() = true;
    }
  }
  else if (st == DISTRIBUTED_STORAGE)
  {
    localstorage = false;
//  if (has(CRE))
//    set_op_array(CRE).set_local() = false;
//  if (has(CRE_DES))
//    set_op_array(CRE_DES).set_local() = false;
//  if (has(CRE_CRE))
//    set_op_array(CRE_CRE).set_local() = false;
//  if (has(DES_DESCOMP))
//    set_op_array(DES_DESCOMP).set_local() = false;
//  if (has(CRE_DESCOMP))
//    set_op_array(CRE_DESCOMP).set_local() = false;
//  if (has(CRE_CRE_DESCOMP))
//    set_op_array(CRE_CRE_DESCOMP).set_local() = false;
    for (std::map<opTypes, boost::shared_ptr< Op_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
      if((it->first & OP_TYPE_MASK) == HAM) continue;
      it->second->set_local() = false;
//    set_op_array(it->first).set_local() = false;
    }
  }

}

boost::shared_ptr<Op_component_base> make_new_op(const opTypes &optype, const bool &is_core)
{
  boost::shared_ptr<Op_component_base> ret;
//switch(optype)
  switch(optype & OP_TYPE_MASK)
  {
    case CRE:
      ret = boost::shared_ptr<Op_component<Cre> >(new Op_component<Cre>(is_core, optype & GENERIC_MASK));
      break;
    case CRE_DES:
      ret = boost::shared_ptr<Op_component<CreDes> >(new Op_component<CreDes>(is_core, optype & GENERIC_MASK));
      break;
    case CRE_CRE:
      ret = boost::shared_ptr<Op_component<CreCre> >(new Op_component<CreCre>(is_core, optype & GENERIC_MASK));
      break;
    case CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreDesComp> >(new Op_component<CreDesComp>(is_core, optype & GENERIC_MASK));
      break;
    case DES_DESCOMP:
      ret = boost::shared_ptr<Op_component<DesDesComp> >(new Op_component<DesDesComp>(is_core, optype & GENERIC_MASK));
      break;
    case CRE_CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreCreDesComp> >(new Op_component<CreCreDesComp>(is_core, optype & GENERIC_MASK));
      break;
    case HAM:
      ret = boost::shared_ptr<Op_component<Ham> >(new Op_component<Ham>(is_core, optype & GENERIC_MASK));
      break;
  }
  return ret;
}


void SpinBlock::default_op_components(bool complementary_, int nroots)
{
  if (complementary_)
  {
    this->complementary = true;
    this->normal = false;
  }
  else
  {
    this->complementary = false;
    this->normal = true;
  }
  this->direct = false;

  ops[CRE] = make_new_op(CRE, true);
  ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
  ops[HAM] = make_new_op(HAM, true);
  if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
    ops[CRE_DES] = make_new_op(CRE_DES, true);
    ops[CRE_CRE] = make_new_op(CRE_CRE, true);
    ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, true);
    ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, true);
  }

  // for DMRG-LRT calculation
  if(dmrginp.calc_type() == DMRG_LRT) {
    for(int i = 1; i < nroots; ++i) {
      opTypes state_index_0i = make_state_index(0, i);
      opTypes state_index_i0 = make_state_index(i, 0);

      ops[CRE | state_index_0i] = make_new_op(CRE | state_index_0i, true);
      ops[CRE | state_index_i0] = make_new_op(CRE | state_index_i0, true);

      ops[CRE_CRE_DESCOMP | state_index_0i] = make_new_op(CRE_CRE_DESCOMP | state_index_0i, true);
      ops[CRE_CRE_DESCOMP | state_index_i0] = make_new_op(CRE_CRE_DESCOMP | state_index_i0, true);

      ops[HAM | state_index_0i] = make_new_op(HAM | state_index_0i, true);

      if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
        ops[CRE_DES | state_index_0i] = make_new_op(CRE_DES | state_index_0i, true);
        ops[CRE_DES | state_index_i0] = make_new_op(CRE_DES | state_index_i0, true);

        ops[CRE_CRE | state_index_0i] = make_new_op(CRE_CRE | state_index_0i, true);
        ops[CRE_CRE | state_index_i0] = make_new_op(CRE_CRE | state_index_i0, true);

        ops[CRE_DESCOMP | state_index_0i] = make_new_op(CRE_DESCOMP | state_index_0i, true);
        ops[CRE_DESCOMP | state_index_i0] = make_new_op(CRE_DESCOMP | state_index_i0, true);

        ops[DES_DESCOMP | state_index_0i] = make_new_op(DES_DESCOMP | state_index_0i, true);
        ops[DES_DESCOMP | state_index_i0] = make_new_op(DES_DESCOMP | state_index_i0, true);
      }
    }
  }

  this->loopblock = true;

}


void SpinBlock::set_big_components()
{
  setstoragetype(DISTRIBUTED_STORAGE);

  ops[HAM] = make_new_op(HAM, false);
}

//this is used for the dot block
void SpinBlock::default_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, bool haveNormops, bool haveCompops, int nroots)
{
  this->direct = direct;
  if (lBlock.is_complementary() || rBlock.is_complementary())
  {
    this->complementary = true;
    this->normal = false;
  }
  else
  {
    this->complementary = false;
    this->normal = true;
  }

  if (!is_direct() )
  {
    ops[CRE] = make_new_op(CRE, true);
    ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
    ops[HAM] = make_new_op(HAM, true);
    if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
      if (haveNormops) {
	ops[CRE_DES] = make_new_op(CRE_DES, true);
	ops[CRE_CRE] = make_new_op(CRE_CRE, true);
      }
      if (haveCompops) {
	ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, true);
	ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, true);
      }
    }

    // for DMRG-LRT calculation
    if(dmrginp.calc_type() == DMRG_LRT) {
      for(int i = 1; i < nroots; ++i) {
        opTypes state_index_0i = make_state_index(0, i);
        opTypes state_index_i0 = make_state_index(i, 0);

        ops[CRE | state_index_0i] = make_new_op(CRE | state_index_0i, true);
        ops[CRE | state_index_i0] = make_new_op(CRE | state_index_i0, true);

        ops[CRE_CRE_DESCOMP | state_index_0i] = make_new_op(CRE_CRE_DESCOMP | state_index_0i, true);
        ops[CRE_CRE_DESCOMP | state_index_i0] = make_new_op(CRE_CRE_DESCOMP | state_index_i0, true);

        ops[HAM | state_index_0i] = make_new_op(HAM | state_index_0i, true);

        if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
          if (haveNormops) {
            ops[CRE_DES | state_index_0i] = make_new_op(CRE_DES | state_index_0i, true);
            ops[CRE_DES | state_index_i0] = make_new_op(CRE_DES | state_index_i0, true);

            ops[CRE_CRE | state_index_0i] = make_new_op(CRE_CRE | state_index_0i, true);
            ops[CRE_CRE | state_index_i0] = make_new_op(CRE_CRE | state_index_i0, true);
          }

          if (haveCompops) {
            ops[CRE_DESCOMP | state_index_0i] = make_new_op(CRE_DESCOMP | state_index_0i, true);
            ops[CRE_DESCOMP | state_index_i0] = make_new_op(CRE_DESCOMP | state_index_i0, true);

            ops[DES_DESCOMP | state_index_0i] = make_new_op(DES_DESCOMP | state_index_0i, true);
            ops[DES_DESCOMP | state_index_i0] = make_new_op(DES_DESCOMP | state_index_i0, true);
          }
        }
      }
    }

    if (haveNormops)
      this->loopblock = true;
    else
      this->loopblock = false;
  }
  else
  {
    // op_components for a single dot block
    ops[CRE] = make_new_op(CRE, false); //this should definitely be false, we not have copies of CRE is all the procs
    ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
    ops[HAM] = make_new_op(HAM, true);
    if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
      if (haveNormops || dmrginp.do_cd()) {
	ops[CRE_DES] = make_new_op(CRE_DES, false);
	ops[CRE_CRE] = make_new_op(CRE_CRE, false);
      }
      if (haveCompops) {
	ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, false);
	ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, false);
      }
    }

    // for DMRG-LRT calculation
    if(dmrginp.calc_type() == DMRG_LRT) {
      for(int i = 1; i < nroots; ++i) {
        opTypes state_index_0i = make_state_index(0, i);
        opTypes state_index_i0 = make_state_index(i, 0);

        ops[CRE | state_index_0i] = make_new_op(CRE | state_index_0i, false);
        ops[CRE | state_index_i0] = make_new_op(CRE | state_index_i0, false);

        ops[CRE_CRE_DESCOMP | state_index_0i] = make_new_op(CRE_CRE_DESCOMP | state_index_0i, true);
        ops[CRE_CRE_DESCOMP | state_index_i0] = make_new_op(CRE_CRE_DESCOMP | state_index_i0, true);

        ops[HAM | state_index_0i] = make_new_op(HAM | state_index_0i, true);

        if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) {
          if (haveNormops || dmrginp.do_cd()) {
            ops[CRE_DES | state_index_0i] = make_new_op(CRE_DES | state_index_0i, false);
            ops[CRE_DES | state_index_i0] = make_new_op(CRE_DES | state_index_i0, false);

            ops[CRE_CRE | state_index_0i] = make_new_op(CRE_CRE | state_index_0i, false);
            ops[CRE_CRE | state_index_i0] = make_new_op(CRE_CRE | state_index_i0, false);
          }

          if (haveCompops) {
            ops[CRE_DESCOMP | state_index_0i] = make_new_op(CRE_DESCOMP | state_index_0i, false);
            ops[CRE_DESCOMP | state_index_i0] = make_new_op(CRE_DESCOMP | state_index_i0, false);

            ops[DES_DESCOMP | state_index_0i] = make_new_op(DES_DESCOMP | state_index_0i, false);
            ops[DES_DESCOMP | state_index_i0] = make_new_op(DES_DESCOMP | state_index_i0, false);
          }
        }
      }
    }

    if (haveNormops)
      this->loopblock = true;
    else
      this->loopblock = false;
  }

}
}
