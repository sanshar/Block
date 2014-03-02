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
    if (has(CRE))
      set_op_array(CRE).set_local() = true;
    if (has(CRE_DES))
      set_op_array(CRE_DES).set_local() = true;
    if (has(CRE_CRE))
      set_op_array(CRE_CRE).set_local() = true;
    if (has(DES_DESCOMP))
      set_op_array(DES_DESCOMP).set_local() = true;
    if (has(CRE_DESCOMP))
      set_op_array(CRE_DESCOMP).set_local() = true;
    if (has(CRE_CRE_DESCOMP))
      set_op_array(CRE_CRE_DESCOMP).set_local() = true;

  }
  else if (st == DISTRIBUTED_STORAGE)
  {
    localstorage = false;
    if (has(CRE))
      set_op_array(CRE).set_local() = false;
    if (has(CRE_DES))
      set_op_array(CRE_DES).set_local() = false;
    if (has(CRE_CRE))
      set_op_array(CRE_CRE).set_local() = false;
    if (has(DES_DESCOMP))
      set_op_array(DES_DESCOMP).set_local() = false;
    if (has(CRE_DESCOMP))
      set_op_array(CRE_DESCOMP).set_local() = false;
    if (has(CRE_CRE_DESCOMP))
      set_op_array(CRE_CRE_DESCOMP).set_local() = false;
  }
  //this is needed for onepdm generation, the system block all the cre are local
  //and on the environment block all the cre are distributed, this way in multiple
  //processor runs, we can generate all O_{ij} elements of onepdm where i is on 
  //the system and j is on the environment
  else if (st == DISTRIBUTED_STORAGE_FOR_ONEPDM)
  {
    localstorage = false;
    if (has(CRE))
      set_op_array(CRE).set_local() = true;
    if (has(CRE_DES))
      set_op_array(CRE_DES).set_local() = false;
    if (has(CRE_CRE))
      set_op_array(CRE_CRE).set_local() = false;
    if (has(DES_DESCOMP))
      set_op_array(DES_DESCOMP).set_local() = false;
    if (has(CRE_DESCOMP))
      set_op_array(CRE_DESCOMP).set_local() = false;
    if (has(CRE_CRE_DESCOMP))
      set_op_array(CRE_CRE_DESCOMP).set_local() = false;
  }


}

boost::shared_ptr<Op_component_base> make_new_op(const opTypes &optype, const bool &is_core)
{
  boost::shared_ptr<Op_component_base> ret;
  switch(optype)
  {
    case CRE:
      ret = boost::shared_ptr<Op_component<Cre> >(new Op_component<Cre>(is_core));
      break;
    case CRE_DES:
      ret = boost::shared_ptr<Op_component<CreDes> >(new Op_component<CreDes>(is_core));
      break;
    case CRE_CRE:
      ret = boost::shared_ptr<Op_component<CreCre> >(new Op_component<CreCre>(is_core));
      break;
    case CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreDesComp> >(new Op_component<CreDesComp>(is_core));
      break;
    case DES_DESCOMP:
      ret = boost::shared_ptr<Op_component<DesDesComp> >(new Op_component<DesDesComp>(is_core));
      break;
    case CRE_CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreCreDesComp> >(new Op_component<CreCreDesComp>(is_core));
      break;
    case HAM:
      ret = boost::shared_ptr<Op_component<Ham> >(new Op_component<Ham>(is_core));
      break;
  }
  return ret;
}

//this is used for the dot block
void SpinBlock::default_op_components(bool complementary_)
{
  complementary = complementary_;
  normal = !complementary_;

  this->direct = false;

  //for a dot operator generate all possible operators
  //they are not rigorously needed in all possible scenarios, e.g. not needed
  //for hubbard model. But they are so cheap that there is no need to have special
  //cases
  ops[CRE] = make_new_op(CRE, true);
  ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
  ops[HAM] = make_new_op(HAM, true);

  ops[CRE_DES] = make_new_op(CRE_DES, true);
  ops[CRE_CRE] = make_new_op(CRE_CRE, true);
  ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, true);
  ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, true);

  this->loopblock = true;

}


void SpinBlock::set_big_components()
{
  setstoragetype(DISTRIBUTED_STORAGE);

  ops[HAM] = make_new_op(HAM, false);
}

void SpinBlock::default_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, bool haveNormops, bool haveCompops)
{
  this->direct = direct;
  if (lBlock.is_complementary() || rBlock.is_complementary()) {
    this->complementary = true;
    this->normal = false;
  } else {
    this->complementary = false;
    this->normal = true;
  }

  if (!is_direct()) {
    ops[CRE] = make_new_op(CRE, true);
    ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
    ops[HAM] = make_new_op(HAM, true);
    //for hubbard model if we want to calculate twopdm we still need cd operators
    if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY || dmrginp.hamiltonian() == BCS || dmrginp.do_cd()) {
      if (haveNormops) {
	    ops[CRE_DES] = make_new_op(CRE_DES, true);
	    ops[CRE_CRE] = make_new_op(CRE_CRE, true);
      }
      if (haveCompops) {
	    ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, true);
	    ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, true);
      }
    }
    if (haveNormops)
      this->loopblock = true;
    else
      this->loopblock = false;
  } else {
    //we need CCDcomp to be on core, the rest of them can be generated very quickly
    //and dont really required incore storage
    ops[CRE] = make_new_op(CRE, false); 
    ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
    ops[HAM] = make_new_op(HAM, true);
    
    //for hubbard model if we want to calculate twopdm we still need cd operators
    if (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY || dmrginp.hamiltonian() == BCS || dmrginp.do_cd()) {
      if (haveNormops || dmrginp.do_cd()) {
	    ops[CRE_DES] = make_new_op(CRE_DES, false);
	    ops[CRE_CRE] = make_new_op(CRE_CRE, false);
      }
      if (haveCompops) {
	    ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, false);
	    ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, false);
      }
    }
    if (haveNormops)
      this->loopblock = true;
    else
      this->loopblock = false;
  }

}
}
