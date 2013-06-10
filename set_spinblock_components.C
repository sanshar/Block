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
//FIXME MAW 3PDM
    if (has(DES_CRE))
      set_op_array(DES_CRE).set_local() = true;
    if (has(CRE_CRE_CRE))
      set_op_array(CRE_CRE_CRE).set_local() = true;
    if (has(CRE_CRE_DES))
      set_op_array(CRE_CRE_DES).set_local() = true;
    if (has(CRE_DES_DES))
      set_op_array(CRE_DES_DES).set_local() = true;
    if (has(CRE_DES_CRE))
      set_op_array(CRE_DES_CRE).set_local() = true;

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
//FIXME MAW 3PDM
    if (has(DES_CRE))
      set_op_array(DES_CRE).set_local() = false;
    if (has(CRE_CRE_CRE))
      set_op_array(CRE_CRE_CRE).set_local() = false;
    if (has(CRE_CRE_DES))
      set_op_array(CRE_CRE_DES).set_local() = false;
    if (has(CRE_DES_DES))
      set_op_array(CRE_DES_DES).set_local() = false;
    if (has(CRE_DES_CRE))
      set_op_array(CRE_DES_CRE).set_local() = false;
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
//FIXME MAW 3PDM
    case DES_CRE:
pout << "allocating new DES_CRE\n";
      ret = boost::shared_ptr<Op_component<DesCre> >(new Op_component<DesCre>(is_core));
      break;
    case CRE_CRE_CRE:
pout << "allocating new CreCreCre\n";
      ret = boost::shared_ptr<Op_component<CreCreCre> >(new Op_component<CreCreCre>(is_core));
      break;
    case CRE_CRE_DES:
pout << "allocating new CreCreDes\n";
      ret = boost::shared_ptr<Op_component<CreCreDes> >(new Op_component<CreCreDes>(is_core));
      break;
    case CRE_DES_DES:
pout << "allocating new CreDesDes\n";
      ret = boost::shared_ptr<Op_component<CreDesDes> >(new Op_component<CreDesDes>(is_core));
      break;
    case CRE_DES_CRE:
pout << "allocating new CreDesCre\n";
      ret = boost::shared_ptr<Op_component<CreDesCre> >(new Op_component<CreDesCre>(is_core));
      break;
  }
  return ret;
}


void SpinBlock::default_op_components(bool complementary_)
{
pout << "SpinBlock::default_op_components(bool complementary_)\n";
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
//FIXME MAW 3PDM
    if (dmrginp.do_3ops()) {
//pout << "maw setting 3ops\n";
      ops[DES_CRE] = make_new_op(DES_CRE, true);
      ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, true);
      ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, true);
      ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, true);
      ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, true);
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
void SpinBlock::default_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, bool haveNormops, bool haveCompops)
{
pout << "SpinBlock::default_op_components(..........) for dot block\n";
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
//FIXME MAW 3PDM
//pout << "maw setting 3ops?\n";
        if (dmrginp.do_3ops()) {
//pout << "maw setting 3ops\n";
          ops[DES_CRE] = make_new_op(DES_CRE, true);
          ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, true);
          ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, true);
          ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, true);
          ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, true);
        }
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
//FIXME MAW 3PDM
//pout << "maw setting 3ops?\n";
        if (dmrginp.do_3ops()) {
//pout << "maw setting 3ops\n";
          ops[DES_CRE] = make_new_op(DES_CRE, false);
          ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, false);
          ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, false);
          ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, false);
          ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, false);
        }
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
