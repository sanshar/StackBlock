/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "Stackspinblock.h"

namespace SpinAdapted{

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void StackSpinBlock::setstoragetype(Storagetype st)
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


    if (has(DES_CRE_DES))
      set_op_array(DES_CRE_DES)->set_local() = true;
    if (has(DES_DES_CRE))
      set_op_array(DES_DES_CRE)->set_local() = true;
    if (has(DES_CRE_CRE))
      set_op_array(DES_CRE_CRE)->set_local() = true;
    if (has(DES_DES_DES))
      set_op_array(DES_DES_DES)->set_local() = true;

  }
  else if (st == DISTRIBUTED_STORAGE)
  {
    localstorage = false;
    //ONLY CRE IS LOCAL
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

    if (has(DES_CRE_DES))
      set_op_array(DES_CRE_DES)->set_local() = false;
    if (has(DES_DES_CRE))
      set_op_array(DES_DES_CRE)->set_local() = false;
    if (has(DES_CRE_CRE))
      set_op_array(DES_CRE_CRE)->set_local() = false;
    if (has(DES_DES_DES))
      set_op_array(DES_DES_DES)->set_local() = false;

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

boost::shared_ptr<StackOp_component_base> make_new_stackop(const opTypes &optype, const bool &is_core)
{
  boost::shared_ptr<StackOp_component_base> ret;
  switch(optype)
  {
    case CRE:
      ret = boost::shared_ptr<StackOp_component<StackCre> >(new StackOp_component<StackCre>(is_core));
      break;
    case DES:
      ret = boost::shared_ptr<StackOp_component<StackDes> >(new StackOp_component<StackDes>(is_core));
      break;
    case CRE_DES:
      ret = boost::shared_ptr<StackOp_component<StackCreDes> >(new StackOp_component<StackCreDes>(is_core));
      break;
    case DES_CRE:
      ret = boost::shared_ptr<StackOp_component<StackDesCre> >(new StackOp_component<StackDesCre>(is_core));
      break;
    case CRE_CRE:
      ret = boost::shared_ptr<StackOp_component<StackCreCre> >(new StackOp_component<StackCreCre>(is_core));
      break;
    case DES_DES:
      ret = boost::shared_ptr<StackOp_component<StackDesDes> >(new StackOp_component<StackDesDes>(is_core));
      break;
    case CRE_DESCOMP:
      ret = boost::shared_ptr<StackOp_component<StackCreDesComp> >(new StackOp_component<StackCreDesComp>(is_core));
      break;
    case DES_CRECOMP:
      ret = boost::shared_ptr<StackOp_component<StackDesCreComp> >(new StackOp_component<StackDesCreComp>(is_core));
      break;
    case DES_DESCOMP:
      ret = boost::shared_ptr<StackOp_component<StackDesDesComp> >(new StackOp_component<StackDesDesComp>(is_core));
      break;
    case CRE_CRECOMP:
      ret = boost::shared_ptr<StackOp_component<StackCreCreComp> >(new StackOp_component<StackCreCreComp>(is_core));
      break;
    case CRE_CRE_DESCOMP:
      ret = boost::shared_ptr<StackOp_component<StackCreCreDesComp> >(new StackOp_component<StackCreCreDesComp>(is_core));
      break;
    case CRE_DES_DESCOMP:
      ret = boost::shared_ptr<StackOp_component<StackCreDesDesComp> >(new StackOp_component<StackCreDesDesComp>(is_core));
      break;
    case HAM:
      ret = boost::shared_ptr<StackOp_component<StackHam> >(new StackOp_component<StackHam>(is_core));
      break;
    case OVERLAP:
      ret = boost::shared_ptr<StackOp_component<StackOverlap> >(new StackOp_component<StackOverlap>(is_core));
      break;
    case RI_3INDEX:
      ret = boost::shared_ptr<StackOp_component<RI3index> >(new StackOp_component<RI3index>(is_core));
      break;
    case RI_4INDEX:
      ret = boost::shared_ptr<StackOp_component<RI4index> >(new StackOp_component<RI4index>(is_core));
      break;

    case CRE_CRE_DES:
      ret = boost::shared_ptr<StackOp_component<StackCreCreDes> >(new StackOp_component<StackCreCreDes>(is_core));
      break;
    case CRE_DES_DES:
      ret = boost::shared_ptr<StackOp_component<StackCreDesDes> >(new StackOp_component<StackCreDesDes>(is_core));
      break;
    case CRE_DES_CRE:
      ret = boost::shared_ptr<StackOp_component<StackCreDesCre> >(new StackOp_component<StackCreDesCre>(is_core));
      break;
    case CRE_CRE_CRE:
      ret = boost::shared_ptr<StackOp_component<StackCreCreCre> >(new StackOp_component<StackCreCreCre>(is_core));
      break;
    // 4PDM and tran_4PDM
    case DES_CRE_DES:
      ret = boost::shared_ptr<StackOp_component<StackDesCreDes> >(new StackOp_component<StackDesCreDes>(is_core));
      break;
    case DES_DES_CRE:
      ret = boost::shared_ptr<StackOp_component<StackDesDesCre> >(new StackOp_component<StackDesDesCre>(is_core));
      break;
    case DES_CRE_CRE:
      ret = boost::shared_ptr<StackOp_component<StackDesCreCre> >(new StackOp_component<StackDesCreCre>(is_core));
      break;
    case DES_DES_DES:
      ret = boost::shared_ptr<StackOp_component<StackDesDesDes> >(new StackOp_component<StackDesDesDes>(is_core));
      break;
    // 4PDM and tran_4PDM
     /*
    case CRE_DES_CRE_DES:
      ret = boost::shared_ptr<StackOp_component<StackCreDesCreDes> >(new StackOp_component<StackCreDesCreDes>(is_core));
      break;
    case CRE_DES_DES_CRE:
      ret = boost::shared_ptr<StackOp_component<StackCreDesDesCre> >(new StackOp_component<StackCreDesDesCre>(is_core));
      break;
    case CRE_CRE_DES_DES:
      ret = boost::shared_ptr<StackOp_component<StackCreCreDesDes> >(new StackOp_component<StackCreCreDesDes>(is_core));
      break;
      */
    default:
      assert(false);
      break;
  }
  return ret;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//this is used for the dot block
void StackSpinBlock::default_op_components(bool complementary_, bool implicitTranspose)
{
  // New version of NPDM code not yet working with implicit transposes
  if ( dmrginp.new_npdm_code() ) implicitTranspose = false;
  if ( !dmrginp.doimplicitTranspose() ) implicitTranspose = false; //this is usually used for testing

  //************
  //implicitTranspose = false;
  //*******

  complementary = complementary_;
  normal = !complementary_;

  this->direct = false;

  //for a dot operator generate all possible operators
  //they are not rigorously needed in all possible scenarios, e.g. not needed
  //for hubbard model. But they are so cheap that there is no need to have special
  //cases
  ops[CRE] = make_new_stackop(CRE, true);
  if (!dmrginp.do_npdm_ops()) {
    ops[CRE_CRE_DESCOMP] = make_new_stackop(CRE_CRE_DESCOMP, true);
    ops[HAM] = make_new_stackop(HAM, true);
  }

  if (!implicitTranspose ) {
    ops[DES] = make_new_stackop(DES, true);
    ops[DES_CRE] = make_new_stackop(DES_CRE, true);
    ops[DES_DES] = make_new_stackop(DES_DES, true);
    if (!dmrginp.do_npdm_ops()) {
      ops[CRE_DES_DESCOMP] = make_new_stackop(CRE_DES_DESCOMP, true);
      ops[DES_CRECOMP] = make_new_stackop(DES_CRECOMP, true);
      ops[CRE_CRECOMP] = make_new_stackop(CRE_CRECOMP, true);
    }
  }

  ops[OVERLAP] = make_new_stackop(OVERLAP, true);

  ops[CRE_DES] = make_new_stackop(CRE_DES, true);
  ops[CRE_CRE] = make_new_stackop(CRE_CRE, true);
  if (!dmrginp.do_npdm_ops()) {
    ops[CRE_DESCOMP] = make_new_stackop(CRE_DESCOMP, true);
    ops[DES_DESCOMP] = make_new_stackop(DES_DESCOMP, true);
  }

  if (dmrginp.do_npdm_ops()) {
    ops[RI_3INDEX] = make_new_stackop(RI_3INDEX, true);
    ops[RI_4INDEX] = make_new_stackop(RI_4INDEX, true);

    if (!dmrginp.set_fullrestart()) {
      if ( (dmrginp.calc_type() == THREEPDM) ||
           (dmrginp.calc_type() == RESTART_THREEPDM)  ||
           (dmrginp.calc_type() == TRANSITION_THREEPDM) ||
           (dmrginp.calc_type() == RESTART_T_THREEPDM)
         ){
      ops[CRE_CRE_CRE] = make_new_stackop(CRE_CRE_CRE, true);
      ops[CRE_DES_DES] = make_new_stackop(CRE_DES_DES, true);
      ops[CRE_CRE_DES] = make_new_stackop(CRE_CRE_DES, true);
      ops[CRE_DES_CRE] = make_new_stackop(CRE_DES_CRE, true);
      if ( !dmrginp.doimplicitTranspose())
      {
          ops[DES_CRE_DES] = make_new_stackop(DES_CRE_DES, true);
          ops[DES_DES_CRE] = make_new_stackop(DES_DES_CRE, true);
          ops[DES_CRE_CRE] = make_new_stackop(DES_CRE_CRE, true);
          ops[DES_DES_DES] = make_new_stackop(DES_DES_DES, true);
      }
      }
    }

  }

  this->loopblock = true;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void StackSpinBlock::set_big_components()
{
  setstoragetype(DISTRIBUTED_STORAGE);

  ops[HAM] = make_new_stackop(HAM, false);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void StackSpinBlock::default_op_components(bool direct, bool haveNormops, bool haveCompops, bool implicitTranspose)
{
  // New version of NPDM code not yet working with implicit transposes
  if ( dmrginp.new_npdm_code() ) implicitTranspose = false;
  if ( !dmrginp.doimplicitTranspose() ) implicitTranspose = false; //this is usually used for testing

  this->direct = direct;

  //************
  //implicitTranspose = false;
  //*******

  // Not direct
  //------------------
  if (!is_direct()) {
    if ( dmrginp.new_npdm_code() && sites.size() > 1) assert(false);

    ops[CRE] = make_new_stackop(CRE, true);
    if (!dmrginp.do_npdm_ops()) {
      ops[CRE_CRE_DESCOMP] = make_new_stackop(CRE_CRE_DESCOMP, true);
      ops[HAM] = make_new_stackop(HAM, true);
    }
    ops[OVERLAP] = make_new_stackop(OVERLAP, true);

    //this option is used when bra and ket states are different
    if (!implicitTranspose ) {
      ops[DES] = make_new_stackop(DES, true);
      if (!dmrginp.do_npdm_ops())
	ops[CRE_DES_DESCOMP] = make_new_stackop(CRE_DES_DESCOMP, true);
    }

    if(dmrginp.do_npdm_ops() && (dmrginp.calc_type() == RESTART_ONEPDM || dmrginp.calc_type() == ONEPDM))
      return;

    //for hubbard model if we want to calculate twopdm we still need cd operators
    if (dmrginp.hamiltonian() != HUBBARD || dmrginp.do_npdm_ops()) {
      if (haveNormops || dmrginp.do_npdm_ops()) {
        ops[CRE_DES] = make_new_stackop(CRE_DES, true);
        ops[CRE_CRE] = make_new_stackop(CRE_CRE, true);
        if (!implicitTranspose) {
          ops[DES_CRE] = make_new_stackop(DES_CRE, true);
          ops[DES_DES] = make_new_stackop(DES_DES, true);
        }
      }
      if (haveCompops && !dmrginp.do_npdm_ops()) {
        ops[CRE_DESCOMP] = make_new_stackop(CRE_DESCOMP, true);
        ops[DES_DESCOMP] = make_new_stackop(DES_DESCOMP, true);
        if (!implicitTranspose) {
          ops[DES_CRECOMP] = make_new_stackop(DES_CRECOMP, true);
          ops[CRE_CRECOMP] = make_new_stackop(CRE_CRECOMP, true);
        }
      }
      if (dmrginp.do_npdm_ops()) {
	ops[RI_3INDEX] = make_new_stackop(RI_3INDEX, true);
	ops[RI_4INDEX] = make_new_stackop(RI_4INDEX, true);

	if (!dmrginp.set_fullrestart()) {
	  if ( (dmrginp.calc_type() == THREEPDM) ||
	       (dmrginp.calc_type() == RESTART_THREEPDM)  ||
	       (dmrginp.calc_type() == TRANSITION_THREEPDM) ||
	       (dmrginp.calc_type() == RESTART_T_THREEPDM)
	       ){
	    ops[CRE_CRE_CRE] = make_new_stackop(CRE_CRE_CRE, true);
	    ops[CRE_DES_DES] = make_new_stackop(CRE_DES_DES, true);
	    ops[CRE_CRE_DES] = make_new_stackop(CRE_CRE_DES, true);
	    ops[CRE_DES_CRE] = make_new_stackop(CRE_DES_CRE, true);
	    if (!dmrginp.doimplicitTranspose()) {
	      ops[DES_CRE_DES] = make_new_stackop(DES_CRE_DES, true);
	      ops[DES_DES_CRE] = make_new_stackop(DES_DES_CRE, true);
	      ops[DES_CRE_CRE] = make_new_stackop(DES_CRE_CRE, true);
	      ops[DES_DES_DES] = make_new_stackop(DES_DES_DES, true);
	    }
	  }
	}
      }
    }


    if (haveNormops)
      this->loopblock = true;
    else
      this->loopblock = false;
  } 

  // Is direct
  //------------------
  else {
    //we need CCDcomp to be on core, the rest of them can be generated very quickly
    //and dont really required incore storage
    ops[CRE] = make_new_stackop(CRE, false); 
    if (!dmrginp.do_npdm_ops()) {
      ops[CRE_CRE_DESCOMP] = make_new_stackop(CRE_CRE_DESCOMP, false);
      ops[HAM] = make_new_stackop(HAM, false);
    }
    ops[OVERLAP] = make_new_stackop(OVERLAP, false);

    //this option is used when bra and ket states are different
    if (!implicitTranspose) {
      ops[DES] = make_new_stackop(DES, false);
      if (!dmrginp.do_npdm_ops())
	ops[CRE_DES_DESCOMP] = make_new_stackop(CRE_DES_DESCOMP, false);
    }

    //Only need one index operators for one pdm calculation except for single dot.
    if(dmrginp.do_npdm_ops() && (dmrginp.calc_type() == RESTART_ONEPDM || dmrginp.calc_type() == ONEPDM))
      return;
    
    //for hubbard model if we want to calculate twopdm we still need cd operators
    if (dmrginp.hamiltonian() != HUBBARD || dmrginp.do_npdm_ops()) {
      if (haveNormops || dmrginp.do_npdm_ops()) {
        ops[CRE_DES] = make_new_stackop(CRE_DES, false);
        ops[CRE_CRE] = make_new_stackop(CRE_CRE, false);
        if (!implicitTranspose ) {
          ops[DES_CRE] = make_new_stackop(DES_CRE, false);
          ops[DES_DES] = make_new_stackop(DES_DES, false);
        }
      }
      if (haveCompops && !dmrginp.do_npdm_ops()) {
        ops[CRE_DESCOMP] = make_new_stackop(CRE_DESCOMP, false);
        ops[DES_DESCOMP] = make_new_stackop(DES_DESCOMP, false);
        if (!implicitTranspose) {
          ops[DES_CRECOMP] = make_new_stackop(DES_CRECOMP, false);
          ops[CRE_CRECOMP] = make_new_stackop(CRE_CRECOMP, false);
        }
      }

      if (dmrginp.do_npdm_ops()) {
	ops[RI_3INDEX] = make_new_stackop(RI_3INDEX, false);
	ops[RI_4INDEX] = make_new_stackop(RI_4INDEX, false);


	if (!dmrginp.set_fullrestart()) {
	  if ( (dmrginp.calc_type() == THREEPDM) ||
	       (dmrginp.calc_type() == RESTART_THREEPDM)  ||
	       (dmrginp.calc_type() == TRANSITION_THREEPDM) ||
	       (dmrginp.calc_type() == RESTART_T_THREEPDM)
	       ){
	    ops[CRE_CRE_CRE] = make_new_stackop(CRE_CRE_CRE, false);
	    ops[CRE_DES_DES] = make_new_stackop(CRE_DES_DES, false);
	    ops[CRE_CRE_DES] = make_new_stackop(CRE_CRE_DES, false);
	    ops[CRE_DES_CRE] = make_new_stackop(CRE_DES_CRE, false);
	    if (!dmrginp.doimplicitTranspose()) {
	      ops[DES_CRE_DES] = make_new_stackop(DES_CRE_DES, false);
	      ops[DES_DES_CRE] = make_new_stackop(DES_DES_CRE, false);
	      ops[DES_CRE_CRE] = make_new_stackop(DES_CRE_CRE, false);
	      ops[DES_DES_DES] = make_new_stackop(DES_DES_DES, false);
	    }
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
