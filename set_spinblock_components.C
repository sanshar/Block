/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "perturb.h"
#include "spinblock.h"

namespace SpinAdapted{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::setstoragetype(Storagetype st)
{
  if (st == LOCAL_STORAGE)
  {
    localstorage = true;
    if (has(CRE))
      set_op_array(CRE)->set_local() = true;
    if (has(DES))
      set_op_array(DES)->set_local() = true;
    if (has(CRE_DES))
      set_op_array(CRE_DES)->set_local() = true;
    if (has(DES_CRE))
      set_op_array(DES_CRE)->set_local() = true;
    if (has(CRE_CRE))
      set_op_array(CRE_CRE)->set_local() = true;
    if (has(DES_DES))
      set_op_array(DES_DES)->set_local() = true;
    if (has(DES_DESCOMP))
      set_op_array(DES_DESCOMP)->set_local() = true;
    if (has(CRE_CRECOMP))
      set_op_array(CRE_CRECOMP)->set_local() = true;
    if (has(CRE_DESCOMP))
      set_op_array(CRE_DESCOMP)->set_local() = true;
    if (has(DES_CRECOMP))
      set_op_array(DES_CRECOMP)->set_local() = true;
    if (has(CRE_CRE_DESCOMP))
      set_op_array(CRE_CRE_DESCOMP)->set_local() = true;
    if (has(CRE_DES_DESCOMP))
      set_op_array(CRE_DES_DESCOMP)->set_local() = true;
    // NPDM
    if (has(RI_3INDEX))
      set_op_array(RI_3INDEX)->set_local() = true;
    if (has(RI_4INDEX))
      set_op_array(RI_4INDEX)->set_local() = true;
    if (has(CRE_CRE_DES))
      set_op_array(CRE_CRE_DES)->set_local() = true;
    if (has(CRE_DES_DES))
      set_op_array(CRE_DES_DES)->set_local() = true;
    if (has(CRE_DES_CRE))
      set_op_array(CRE_DES_CRE)->set_local() = true;
    if (has(CRE_CRE_CRE))
      set_op_array(CRE_CRE_CRE)->set_local() = true;
    // 4PDM
    if (has(DES_CRE_DES))
      set_op_array(DES_CRE_DES)->set_local() = true;
    if (has(DES_DES_CRE))
      set_op_array(DES_DES_CRE)->set_local() = true;
    if (has(DES_CRE_CRE))
      set_op_array(DES_CRE_CRE)->set_local() = true;
    if (has(DES_DES_DES))
      set_op_array(DES_DES_DES)->set_local() = true;
    if (has(CRE_CRE_DES_DES))
      set_op_array(CRE_CRE_DES_DES)->set_local() = true;
    if (has(CRE_DES_CRE_DES))
      set_op_array(CRE_DES_CRE_DES)->set_local() = true;
    if (has(CRE_DES_DES_CRE))
      set_op_array(CRE_DES_DES_CRE)->set_local() = true;
    if (has(CRE_DES_DES_DES))
      set_op_array(CRE_DES_DES_DES)->set_local() = true;
    if (has(CRE_CRE_CRE_DES))
      set_op_array(CRE_CRE_CRE_DES)->set_local() = true;
    if (has(CRE_CRE_DES_CRE))
      set_op_array(CRE_CRE_DES_CRE)->set_local() = true;
    if (has(CRE_DES_CRE_CRE))
      set_op_array(CRE_DES_CRE_CRE)->set_local() = true;
    if (has(CRE_CRE_CRE_CRE))
      set_op_array(CRE_CRE_CRE_CRE)->set_local() = true;
    //mps_nevpt2
    if (has(CDD_CRE_DESCOMP))
      set_op_array(CDD_CRE_DESCOMP)->set_local() = true;
    if (has(CDD_DES_DESCOMP))
      set_op_array(CDD_DES_DESCOMP)->set_local() = true;
    if (has(CCD_CRE_DESCOMP))
      set_op_array(CCD_CRE_DESCOMP)->set_local() = true;
    if (has(CCD_CRE_CRECOMP))
      set_op_array(CCD_CRE_CRECOMP)->set_local() = true;

  }
  else if (st == DISTRIBUTED_STORAGE)
  {
    localstorage = false;
    if (has(CRE))
      set_op_array(CRE)->set_local() = false;
    if (has(DES))
      set_op_array(DES)->set_local() = false;
    if (has(CRE_DES))
      set_op_array(CRE_DES)->set_local() = false;
    if (has(DES_CRE))
      set_op_array(DES_CRE)->set_local() = false;
    if (has(CRE_CRE))
      set_op_array(CRE_CRE)->set_local() = false;
    if (has(DES_DES))
      set_op_array(DES_DES)->set_local() = false;
    if (has(DES_DESCOMP))
      set_op_array(DES_DESCOMP)->set_local() = false;
    if (has(CRE_CRECOMP))
      set_op_array(CRE_CRECOMP)->set_local() = false;
    if (has(CRE_DESCOMP))
      set_op_array(CRE_DESCOMP)->set_local() = false;
    if (has(DES_CRECOMP))
      set_op_array(DES_CRECOMP)->set_local() = false;
    if (has(CRE_CRE_DESCOMP))
      set_op_array(CRE_CRE_DESCOMP)->set_local() = false;
    if (has(CRE_DES_DESCOMP))
      set_op_array(CRE_DES_DESCOMP)->set_local() = false;
    // NPDM
    if (has(RI_3INDEX))
      set_op_array(RI_3INDEX)->set_local() = false;
    if (has(RI_4INDEX))
      set_op_array(RI_4INDEX)->set_local() = false;
    if (has(CRE_CRE_DES))
      set_op_array(CRE_CRE_DES)->set_local() = false;
    if (has(CRE_DES_DES))
      set_op_array(CRE_DES_DES)->set_local() = false;
    if (has(CRE_DES_CRE))
      set_op_array(CRE_DES_CRE)->set_local() = false;
    if (has(CRE_CRE_CRE))
      set_op_array(CRE_CRE_CRE)->set_local() = false;
    // 4PDM
    if (has(DES_CRE_DES))
      set_op_array(DES_CRE_DES)->set_local() = false;
    if (has(DES_DES_CRE))
      set_op_array(DES_DES_CRE)->set_local() = false;
    if (has(DES_CRE_CRE))
      set_op_array(DES_CRE_CRE)->set_local() = false;
    if (has(DES_DES_DES))
      set_op_array(DES_DES_DES)->set_local() = false;
    if (has(CRE_CRE_DES_DES))
      set_op_array(CRE_CRE_DES_DES)->set_local() = false;
    if (has(CRE_DES_CRE_DES))
      set_op_array(CRE_DES_CRE_DES)->set_local() = false;
    if (has(CRE_DES_DES_CRE))
      set_op_array(CRE_DES_DES_CRE)->set_local() = false;
    if (has(CRE_DES_DES_DES))
      set_op_array(CRE_DES_DES_DES)->set_local() = false;
    if (has(CRE_CRE_CRE_DES))
      set_op_array(CRE_CRE_CRE_DES)->set_local() = false;
    if (has(CRE_CRE_DES_CRE))
      set_op_array(CRE_CRE_DES_CRE)->set_local() = false;
    if (has(CRE_DES_CRE_CRE))
      set_op_array(CRE_DES_CRE_CRE)->set_local() = false;
    if (has(CRE_CRE_CRE_CRE))
      set_op_array(CRE_CRE_CRE_CRE)->set_local() = false;
    //mps_nevpt2
    if (has(CDD_CRE_DESCOMP))
      set_op_array(CDD_CRE_DESCOMP)->set_local() = false;
    if (has(CDD_DES_DESCOMP))
      set_op_array(CDD_DES_DESCOMP)->set_local() = false;
    if (has(CCD_CRE_DESCOMP))
      set_op_array(CCD_CRE_DESCOMP)->set_local() = false;
    if (has(CCD_CRE_CRECOMP))
      set_op_array(CCD_CRE_CRECOMP)->set_local() = false;
  }

  //this is needed for onepdm generation, the system block all the cre are local
  //and on the environment block all the cre are distributed, this way in multiple
  //processor runs, we can generate all O_{ij} elements of onepdm where i is on 
  //the system and j is on the environment
  else if (st == DISTRIBUTED_STORAGE_FOR_ONEPDM)
  {
    if ( dmrginp.new_npdm_code() && !dmrginp.nevpt2() ) assert(false);
    localstorage = false;
    if (has(CRE))
      set_op_array(CRE)->set_local() = true;
    if (has(DES))
      set_op_array(DES)->set_local() = true;
    if (has(CRE_DES))
      set_op_array(CRE_DES)->set_local() = false;
    if (has(CRE_CRE))
      set_op_array(CRE_CRE)->set_local() = false;
    if (has(DES_DESCOMP))
      set_op_array(DES_DESCOMP)->set_local() = false;
    if (has(CRE_DESCOMP))
      set_op_array(CRE_DESCOMP)->set_local() = false;
    if (has(CRE_CRE_DESCOMP))
      set_op_array(CRE_CRE_DESCOMP)->set_local() = false;
  }


}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<Op_component_base> make_new_op(const opTypes &optype, const bool &is_core)
{
  boost::shared_ptr<Op_component_base> ret;
  switch(optype)
  {
    case CRE:
      ret = boost::shared_ptr<Op_component<Cre> >(new Op_component<Cre>(is_core));
      break;
    case DES:
      ret = boost::shared_ptr<Op_component<Des> >(new Op_component<Des>(is_core));
      break;
    case CRE_DES:
      ret = boost::shared_ptr<Op_component<CreDes> >(new Op_component<CreDes>(is_core));
      break;
    case DES_CRE:
      ret = boost::shared_ptr<Op_component<DesCre> >(new Op_component<DesCre>(is_core));
      break;
    case CRE_CRE:
      ret = boost::shared_ptr<Op_component<CreCre> >(new Op_component<CreCre>(is_core));
      break;
    case DES_DES:
      ret = boost::shared_ptr<Op_component<DesDes> >(new Op_component<DesDes>(is_core));
      break;
    case CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreDesComp> >(new Op_component<CreDesComp>(is_core));
      break;
    case DES_CRECOMP:
      ret = boost::shared_ptr<Op_component<DesCreComp> >(new Op_component<DesCreComp>(is_core));
      break;
    case DES_DESCOMP:
      ret = boost::shared_ptr<Op_component<DesDesComp> >(new Op_component<DesDesComp>(is_core));
      break;
    case CRE_CRECOMP:
      ret = boost::shared_ptr<Op_component<CreCreComp> >(new Op_component<CreCreComp>(is_core));
      break;
    case CRE_CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreCreDesComp> >(new Op_component<CreCreDesComp>(is_core));
      break;
    case CRE_DES_DESCOMP:
      ret = boost::shared_ptr<Op_component<CreDesDesComp> >(new Op_component<CreDesDesComp>(is_core));
      break;
    case HAM:
      ret = boost::shared_ptr<Op_component<Ham> >(new Op_component<Ham>(is_core));
      break;
    case OVERLAP:
      ret = boost::shared_ptr<Op_component<Overlap> >(new Op_component<Overlap>(is_core));
      break;
    // NPDM
    case RI_3INDEX:
      ret = boost::shared_ptr<Op_component<RI3index> >(new Op_component<RI3index>(is_core));
      break;
    case RI_4INDEX:
      ret = boost::shared_ptr<Op_component<RI4index> >(new Op_component<RI4index>(is_core));
      break;
    case CRE_CRE_DES:
      ret = boost::shared_ptr<Op_component<CreCreDes> >(new Op_component<CreCreDes>(is_core));
      break;
    case CRE_DES_DES:
      ret = boost::shared_ptr<Op_component<CreDesDes> >(new Op_component<CreDesDes>(is_core));
      break;
    case CRE_DES_CRE:
      ret = boost::shared_ptr<Op_component<CreDesCre> >(new Op_component<CreDesCre>(is_core));
      break;
    case CRE_CRE_CRE:
      ret = boost::shared_ptr<Op_component<CreCreCre> >(new Op_component<CreCreCre>(is_core));
      break;
    // 4PDM
    case DES_CRE_DES:
      ret = boost::shared_ptr<Op_component<DesCreDes> >(new Op_component<DesCreDes>(is_core));
      break;
    case DES_DES_CRE:
      ret = boost::shared_ptr<Op_component<DesDesCre> >(new Op_component<DesDesCre>(is_core));
      break;
    case DES_CRE_CRE:
      ret = boost::shared_ptr<Op_component<DesCreCre> >(new Op_component<DesCreCre>(is_core));
      break;
    case DES_DES_DES:
      ret = boost::shared_ptr<Op_component<DesDesDes> >(new Op_component<DesDesDes>(is_core));
      break;
    case CRE_CRE_DES_DES:
      ret = boost::shared_ptr<Op_component<CreCreDesDes> >(new Op_component<CreCreDesDes>(is_core));
      break;
    case CRE_DES_CRE_DES:
      ret = boost::shared_ptr<Op_component<CreDesCreDes> >(new Op_component<CreDesCreDes>(is_core));
      break;
    case CRE_DES_DES_CRE:
      ret = boost::shared_ptr<Op_component<CreDesDesCre> >(new Op_component<CreDesDesCre>(is_core));
      break;
    case CRE_DES_DES_DES:
      ret = boost::shared_ptr<Op_component<CreDesDesDes> >(new Op_component<CreDesDesDes>(is_core));
      break;
    case CRE_CRE_CRE_DES:
      ret = boost::shared_ptr<Op_component<CreCreCreDes> >(new Op_component<CreCreCreDes>(is_core));
      break;
    case CRE_CRE_DES_CRE:
      ret = boost::shared_ptr<Op_component<CreCreDesCre> >(new Op_component<CreCreDesCre>(is_core));
      break;
    case CRE_DES_CRE_CRE:
      ret = boost::shared_ptr<Op_component<CreDesCreCre> >(new Op_component<CreDesCreCre>(is_core));
      break;
    case CRE_CRE_CRE_CRE:
      ret = boost::shared_ptr<Op_component<CreCreCreCre> >(new Op_component<CreCreCreCre>(is_core));
      break;
    case CDD_SUM:
      ret = boost::shared_ptr<Op_component<CDD_sum> >(new Op_component<CDD_sum>(is_core));
      break;
    case CDD_CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CDD_CreDesComp> >(new Op_component<CDD_CreDesComp>(is_core));
      break;
    case CDD_DES_DESCOMP:
      ret = boost::shared_ptr<Op_component<CDD_DesDesComp> >(new Op_component<CDD_DesDesComp>(is_core));
      break;
    case CCD_SUM:
      ret = boost::shared_ptr<Op_component<CCD_sum> >(new Op_component<CCD_sum>(is_core));
      break;
    case CCD_CRE_DESCOMP:
      ret = boost::shared_ptr<Op_component<CCD_CreDesComp> >(new Op_component<CCD_CreDesComp>(is_core));
      break;
    case CCD_CRE_CRECOMP:
      ret = boost::shared_ptr<Op_component<CCD_CreCreComp> >(new Op_component<CCD_CreCreComp>(is_core));
      break;
    default:
      assert(false);
      break;
  }
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//this is used for the dot block
void SpinBlock::default_op_components(bool complementary_, bool implicitTranspose)
{
  // New version of NPDM code not yet working with implicit transposes
  if ( dmrginp.new_npdm_code() ) implicitTranspose = false;
  if ( !dmrginp.doimplicitTranspose() ) implicitTranspose = false; //this is usually used for testing

  complementary = complementary_;
  normal = !complementary_;

  this->direct = false;
  this->loopblock = true;

  //TODO
  if(dmrginp.calc_type() == MPS_NEVPT)
  {
     ops[CRE] = make_new_op(CRE, true);
     ops[DES] = make_new_op(DES, true);
     ops[OVERLAP] = make_new_op(OVERLAP, true);
     if(this->nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2)){
     //  ops[CRE_DES] = make_new_op(CRE_DES, true);
     //  ops[DES_DES] = make_new_op(DES_DES, true);
       ops[CDD_CRE_DESCOMP] = make_new_op(CDD_CRE_DESCOMP, true);
       ops[CDD_DES_DESCOMP] = make_new_op(CDD_DES_DESCOMP, true);
       ops[CDD_SUM] = make_new_op(CDD_SUM, true);
     }
     else{
     //  ops[CRE_DES] = make_new_op(CRE_DES, true);
     //  ops[DES_DES] = make_new_op(DES_DES, true);
       ops[CCD_CRE_DESCOMP] = make_new_op(CCD_CRE_DESCOMP, true);
       ops[CCD_CRE_CRECOMP] = make_new_op(CCD_CRE_CRECOMP, true);
       ops[CCD_SUM] = make_new_op(CCD_SUM, true);
     }
  return; 
  }
  //for a dot operator generate all possible operators
  //they are not rigorously needed in all possible scenarios, e.g. not needed
  //for hubbard model. But they are so cheap that there is no need to have special
  //cases
  ops[CRE] = make_new_op(CRE, true);
  ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);

  if (!implicitTranspose) {
    ops[DES] = make_new_op(DES, true);
    ops[CRE_DES_DESCOMP] = make_new_op(CRE_DES_DESCOMP, true);
    ops[DES_CRE] = make_new_op(DES_CRE, true);
    ops[DES_CRECOMP] = make_new_op(DES_CRECOMP, true);
    ops[DES_DES] = make_new_op(DES_DES, true);
    ops[CRE_CRECOMP] = make_new_op(CRE_CRECOMP, true);
  }

  ops[HAM] = make_new_op(HAM, true);
  ops[OVERLAP] = make_new_op(OVERLAP, true);

  ops[CRE_DES] = make_new_op(CRE_DES, true);
  ops[CRE_CRE] = make_new_op(CRE_CRE, true);
  ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, true);
  ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, true);

	if ( dmrginp.new_npdm_code() ) {
	  if (dmrginp.npdm_generate() == true)
	  {
	  	//Now it is generating environment block for the npdm sweep.
	  	//On environment block, the number of indices of operators are less than 
	  	//the order of pdm.
      ops[RI_3INDEX] = make_new_op(RI_3INDEX, true);
      ops[RI_4INDEX] = make_new_op(RI_4INDEX, true);
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)  ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE] = make_new_op(DES_CRE, true);
      }
      if ( (dmrginp.calc_type() == FOURPDM)   ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, true);
        ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, true);
        ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, true);
        ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, true);
        ops[DES_CRE_DES] = make_new_op(DES_CRE_DES, true);
        ops[DES_DES_CRE] = make_new_op(DES_DES_CRE, true);
        ops[DES_CRE_CRE] = make_new_op(DES_CRE_CRE, true);
        ops[DES_DES_DES] = make_new_op(DES_DES_DES, true);
      }
    }
    else{
      ops[RI_3INDEX] = make_new_op(RI_3INDEX, true);
      ops[RI_4INDEX] = make_new_op(RI_4INDEX, true);
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)  ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE] = make_new_op(DES_CRE, true);
        ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, true);
        ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, true);
        ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, true);
        ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, true);
      }
      if ( (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)   ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE_DES] = make_new_op(DES_CRE_DES, true);
        ops[DES_DES_CRE] = make_new_op(DES_DES_CRE, true);
        ops[DES_CRE_CRE] = make_new_op(DES_CRE_CRE, true);
        ops[DES_DES_DES] = make_new_op(DES_DES_DES, true);
      }
      if ( (dmrginp.calc_type() == FOURPDM)   ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[CRE_CRE_DES_DES] = make_new_op(CRE_CRE_DES_DES, true);
        ops[CRE_DES_CRE_DES] = make_new_op(CRE_DES_CRE_DES, true);
        ops[CRE_DES_DES_CRE] = make_new_op(CRE_DES_DES_CRE, true);
        ops[CRE_DES_DES_DES] = make_new_op(CRE_DES_DES_DES, true);
        ops[CRE_CRE_CRE_DES] = make_new_op(CRE_CRE_CRE_DES, true);
        ops[CRE_CRE_DES_CRE] = make_new_op(CRE_CRE_DES_CRE, true);
        ops[CRE_DES_CRE_CRE] = make_new_op(CRE_DES_CRE_CRE, true);
        ops[CRE_CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE_CRE, true);
      }
    }
  }


}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::set_big_components()
{
  setstoragetype(DISTRIBUTED_STORAGE);
   if(dmrginp.calc_type() == MPS_NEVPT)
   {
      if(this->nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2)){
        ops[CDD_SUM] = make_new_op(CDD_SUM, false);
      }
      else{
        ops[CCD_SUM] = make_new_op(CCD_SUM, false);
      }
      return; 
   }

   else
     ops[HAM] = make_new_op(HAM, false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::default_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, bool haveNormops, bool haveCompops, bool implicitTranspose)
{
  // New version of NPDM code not yet working with implicit transposes
  if ( dmrginp.new_npdm_code() ) implicitTranspose = false;
  if ( !dmrginp.doimplicitTranspose() ) implicitTranspose = false; //this is usually used for testing

  this->direct = direct;
  if (lBlock.is_complementary() || rBlock.is_complementary()) {
    this->complementary = true;
    this->normal = false;
  } else {
    this->complementary = false;
    this->normal = true;
  }

  if (haveNormops)
    this->loopblock = true;
  else
    this->loopblock = false;

  // Not direct
  //------------------
  if (!is_direct()) {
   if(dmrginp.calc_type() == MPS_NEVPT)
   {
      ops[CRE] = make_new_op(CRE, true);
      ops[DES] = make_new_op(DES, true);
      ops[OVERLAP] = make_new_op(OVERLAP, true);
      if(lBlock.nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2)){
      //  ops[CRE_DES] = make_new_op(CRE_DES, true);
      //  ops[DES_DES] = make_new_op(DES_DES, true);
        ops[CDD_CRE_DESCOMP] = make_new_op(CDD_CRE_DESCOMP, true);
        ops[CDD_DES_DESCOMP] = make_new_op(CDD_DES_DESCOMP, true);
        ops[CDD_SUM] = make_new_op(CDD_SUM, true);
      }
      else{
      //  ops[CRE_DES] = make_new_op(CRE_DES, true);
      //  ops[DES_DES] = make_new_op(DES_DES, true);
        ops[CCD_CRE_DESCOMP] = make_new_op(CCD_CRE_DESCOMP, true);
        ops[CCD_CRE_CRECOMP] = make_new_op(CCD_CRE_CRECOMP, true);
        ops[CCD_SUM] = make_new_op(CCD_SUM, true);
      }
      return; 
   }
    if ( dmrginp.new_npdm_code() && sites.size() > 1) assert(false);

    ops[CRE] = make_new_op(CRE, true);
    ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, true);
    ops[HAM] = make_new_op(HAM, true);
    ops[OVERLAP] = make_new_op(OVERLAP, true);

    //this option is used when bra and ket states are different
    if (!implicitTranspose) {
      ops[DES] = make_new_op(DES, true);
      ops[CRE_DES_DESCOMP] = make_new_op(CRE_DES_DESCOMP, true);
    }

    //for hubbard model if we want to calculate twopdm we still need cd operators
    if (dmrginp.hamiltonian() != HUBBARD || dmrginp.do_npdm_ops()) {
      if (haveNormops) {
        ops[CRE_DES] = make_new_op(CRE_DES, true);
        ops[CRE_CRE] = make_new_op(CRE_CRE, true);
        if (!implicitTranspose) {
          ops[DES_CRE] = make_new_op(DES_CRE, true);
          ops[DES_DES] = make_new_op(DES_DES, true);
        }
      }
      if (haveCompops) {
        ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, true);
        ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, true);
        if (!implicitTranspose) {
          ops[DES_CRECOMP] = make_new_op(DES_CRECOMP, true);
          ops[CRE_CRECOMP] = make_new_op(CRE_CRECOMP, true);
        }
      }
    }


	if ( dmrginp.new_npdm_code() ) {
	  if (dmrginp.npdm_generate() == true)
	  {
	  	//Now it is generating environment block for the npdm sweep.
	  	//On environment block, the number of indices of operators are less than 
	  	//the order of pdm.
      ops[RI_3INDEX] = make_new_op(RI_3INDEX, true);
      ops[RI_4INDEX] = make_new_op(RI_4INDEX, true);
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)  ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE] = make_new_op(DES_CRE, true);
        if ( (dmrginp.calc_type() == FOURPDM)   ||
             (dmrginp.calc_type() == RESTART_FOURPDM)  ||
             (dmrginp.calc_type() == NEVPT2PDM) ||
             (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
          ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, true);
          ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, true);
          ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, true);
          ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, true);
          ops[DES_CRE_DES] = make_new_op(DES_CRE_DES, true);
          ops[DES_DES_CRE] = make_new_op(DES_DES_CRE, true);
          ops[DES_CRE_CRE] = make_new_op(DES_CRE_CRE, true);
          ops[DES_DES_DES] = make_new_op(DES_DES_DES, true);
        }
      }

	  }
    else{
      ops[RI_3INDEX] = make_new_op(RI_3INDEX, true);
      ops[RI_4INDEX] = make_new_op(RI_4INDEX, true);
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)  ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE] = make_new_op(DES_CRE, true);
        ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, true);
        ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, true);
        ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, true);
        ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, true);
      }
      if ( 
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)   ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE_DES] = make_new_op(DES_CRE_DES, true);
        ops[DES_DES_CRE] = make_new_op(DES_DES_CRE, true);
        ops[DES_CRE_CRE] = make_new_op(DES_CRE_CRE, true);
        ops[DES_DES_DES] = make_new_op(DES_DES_DES, true);
      }
      if ( (dmrginp.calc_type() == FOURPDM)   ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[CRE_CRE_DES_DES] = make_new_op(CRE_CRE_DES_DES, true);
        ops[CRE_DES_CRE_DES] = make_new_op(CRE_DES_CRE_DES, true);
        ops[CRE_DES_DES_CRE] = make_new_op(CRE_DES_DES_CRE, true);
        ops[CRE_DES_DES_DES] = make_new_op(CRE_DES_DES_DES, true);
        ops[CRE_CRE_CRE_DES] = make_new_op(CRE_CRE_CRE_DES, true);
        ops[CRE_CRE_DES_CRE] = make_new_op(CRE_CRE_DES_CRE, true);
        ops[CRE_DES_CRE_CRE] = make_new_op(CRE_DES_CRE_CRE, true);
        ops[CRE_CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE_CRE, true);
      }
    }
  }

  } 

  // Is direct
  //------------------
  else {
   if(dmrginp.calc_type() == MPS_NEVPT)
   {
      ops[CRE] = make_new_op(CRE, false);
      ops[DES] = make_new_op(DES, false);
      ops[OVERLAP] = make_new_op(OVERLAP, false);
      if(lBlock.nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2)){
      //  ops[CRE_DES] = make_new_op(CRE_DES, false);
      //  ops[DES_DES] = make_new_op(DES_DES, false);
        ops[CDD_CRE_DESCOMP] = make_new_op(CDD_CRE_DESCOMP, false);
        ops[CDD_DES_DESCOMP] = make_new_op(CDD_DES_DESCOMP, false);
        ops[CDD_SUM] = make_new_op(CDD_SUM, false);
      }
      else{
      //  ops[CRE_DES] = make_new_op(CRE_DES, false);
      //  ops[DES_DES] = make_new_op(DES_DES, false);
        ops[CCD_CRE_DESCOMP] = make_new_op(CCD_CRE_DESCOMP, false);
        ops[CCD_CRE_CRECOMP] = make_new_op(CCD_CRE_CRECOMP, false);
        ops[CCD_SUM] = make_new_op(CCD_SUM, false);
      }
      return; 
   }
    //we need CCDcomp to be on core, the rest of them can be generated very quickly
    //and dont really required incore storage
    ops[CRE] = make_new_op(CRE, false); 
    ops[CRE_CRE_DESCOMP] = make_new_op(CRE_CRE_DESCOMP, false);
    ops[HAM] = make_new_op(HAM, false);
    ops[OVERLAP] = make_new_op(OVERLAP, false);

    //this option is used when bra and ket states are different
    if (!implicitTranspose) {
      ops[DES] = make_new_op(DES, false);
      ops[CRE_DES_DESCOMP] = make_new_op(CRE_DES_DESCOMP, false);
    }
    
    //for hubbard model if we want to calculate twopdm we still need cd operators
    if (dmrginp.hamiltonian() != HUBBARD || dmrginp.do_npdm_ops()) {
      if (haveNormops || dmrginp.do_npdm_ops()) {
        ops[CRE_DES] = make_new_op(CRE_DES, false);
        ops[CRE_CRE] = make_new_op(CRE_CRE, false);
        if (!implicitTranspose) {
          ops[DES_CRE] = make_new_op(DES_CRE, false);
          ops[DES_DES] = make_new_op(DES_DES, false);
        }
      }
      if (haveCompops) {
        ops[CRE_DESCOMP] = make_new_op(CRE_DESCOMP, false);
        ops[DES_DESCOMP] = make_new_op(DES_DESCOMP, false);
        if (!implicitTranspose) {
          ops[DES_CRECOMP] = make_new_op(DES_CRECOMP, false);
          ops[CRE_CRECOMP] = make_new_op(CRE_CRECOMP, false);
        }
      }
    }

	if ( dmrginp.new_npdm_code() ) {
	  if (dmrginp.npdm_generate() == true)
	  {
	  	//Now it is generating environment block for the npdm sweep.
	  	//On environment block, the number of indices of operators are less than 
	  	//the order of pdm.
      ops[RI_3INDEX] = make_new_op(RI_3INDEX, false);
      ops[RI_4INDEX] = make_new_op(RI_4INDEX, false);
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)  ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE] = make_new_op(DES_CRE, false);
        if ( (dmrginp.calc_type() == FOURPDM)   ||
             (dmrginp.calc_type() == RESTART_FOURPDM)  ||
             (dmrginp.calc_type() == NEVPT2PDM) ||
             (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
          ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, false);
          ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, false);
          ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, false);
          ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, false);
          ops[DES_CRE_DES] = make_new_op(DES_CRE_DES, false);
          ops[DES_DES_CRE] = make_new_op(DES_DES_CRE, false);
          ops[DES_CRE_CRE] = make_new_op(DES_CRE_CRE, false);
          ops[DES_DES_DES] = make_new_op(DES_DES_DES, false);
        }
      }

	  }
	  else
    {
      ops[RI_3INDEX] = make_new_op(RI_3INDEX, false);
      ops[RI_4INDEX] = make_new_op(RI_4INDEX, false);
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
           (dmrginp.calc_type() == FOURPDM)  ||
           (dmrginp.calc_type() == RESTART_FOURPDM)  ||
           (dmrginp.calc_type() == NEVPT2PDM) ||
           (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
        ops[DES_CRE] = make_new_op(DES_CRE, false);
        ops[CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE, false);
        ops[CRE_DES_DES] = make_new_op(CRE_DES_DES, false);
        ops[CRE_CRE_DES] = make_new_op(CRE_CRE_DES, false);
        ops[CRE_DES_CRE] = make_new_op(CRE_DES_CRE, false);
        }
        if ( (dmrginp.calc_type() == TRANSITION_THREEPDM)  ||
             (dmrginp.calc_type() == RESTART_T_THREEPDM)  ||
             (dmrginp.calc_type() == FOURPDM)   ||
             (dmrginp.calc_type() == RESTART_FOURPDM)  ||
             (dmrginp.calc_type() == NEVPT2PDM) ||
             (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
          ops[DES_CRE_DES] = make_new_op(DES_CRE_DES, false);
          ops[DES_DES_CRE] = make_new_op(DES_DES_CRE, false);
          ops[DES_CRE_CRE] = make_new_op(DES_CRE_CRE, false);
          ops[DES_DES_DES] = make_new_op(DES_DES_DES, false);
        }
        if ( (dmrginp.calc_type() == FOURPDM)   ||
             (dmrginp.calc_type() == RESTART_FOURPDM)  ||
             (dmrginp.calc_type() == NEVPT2PDM) ||
             (dmrginp.calc_type() == RESTART_NEVPT2PDM) ) {
          ops[CRE_CRE_DES_DES] = make_new_op(CRE_CRE_DES_DES, false);
          ops[CRE_DES_CRE_DES] = make_new_op(CRE_DES_CRE_DES, false);
          ops[CRE_DES_DES_CRE] = make_new_op(CRE_DES_DES_CRE, false);
          ops[CRE_DES_DES_DES] = make_new_op(CRE_DES_DES_DES, false);
          ops[CRE_CRE_CRE_DES] = make_new_op(CRE_CRE_CRE_DES, false);
          ops[CRE_CRE_DES_CRE] = make_new_op(CRE_CRE_DES_CRE, false);
          ops[CRE_DES_CRE_CRE] = make_new_op(CRE_DES_CRE_CRE, false);
          ops[CRE_CRE_CRE_CRE] = make_new_op(CRE_CRE_CRE_CRE, false);
        }
      }
  }

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void SpinBlock::perturb_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, const perturber& pb)
{
  this->direct = direct;
  if (lBlock.is_complementary() || rBlock.is_complementary()) {
    this->complementary = true;
    this->normal = false;
  } else {
    this->complementary = false;
    this->normal = true;
  }
  if (is_direct()) {
    ops[CRE] = make_new_op(CRE, false);
    ops[DES] = make_new_op(DES, false);
    ops[OVERLAP] = make_new_op(OVERLAP, false);

    if(pb.type() == Va){
    //  ops[CRE_DES] = make_new_op(CRE_DES, false);
    //  ops[DES_DES] = make_new_op(DES_DES, false);
      ops[CDD_CRE_DESCOMP] = make_new_op(CDD_CRE_DESCOMP, false);
      ops[CDD_DES_DESCOMP] = make_new_op(CDD_DES_DESCOMP, false);
      ops[CDD_SUM] = make_new_op(CDD_SUM, false);
    }
    else if(pb.type() == TwoPerturbType::Vi){
    //  ops[CRE_DES] = make_new_op(CRE_DES, false);
    //  ops[DES_DES] = make_new_op(DES_DES, false);
      ops[CCD_CRE_DESCOMP] = make_new_op(CCD_CRE_DESCOMP, false);
      ops[CCD_CRE_CRECOMP] = make_new_op(CCD_CRE_CRECOMP, false);
      ops[CCD_SUM] = make_new_op(CCD_SUM, false);
    }
  } 
  else
  {
    ops[CRE] = make_new_op(CRE, true);
    ops[DES] = make_new_op(DES, true);
    ops[OVERLAP] = make_new_op(OVERLAP, true);

    if(pb.type() == TwoPerturbType::Va){
      ops[CDD_CRE_DESCOMP] = make_new_op(CDD_CRE_DESCOMP, true);
      ops[CDD_DES_DESCOMP] = make_new_op(CDD_DES_DESCOMP, true);
      ops[CDD_SUM] = make_new_op(CDD_SUM, true);
    }
    else if(pb.type() == TwoPerturbType::Vi){
      ops[CCD_CRE_DESCOMP] = make_new_op(CCD_CRE_DESCOMP, true);
      ops[CCD_CRE_CRECOMP] = make_new_op(CCD_CRE_CRECOMP, true);
      ops[CCD_SUM] = make_new_op(CCD_SUM, true);
    }

  } 


}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
