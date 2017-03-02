/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "Stackspinblock.h"
#include "Stack_op_components.h"
#include "screen.h"

namespace SpinAdapted {
  // -------------------- C_S1 ---------------------------  
  template<> string StackOp_component<StackCre>::get_op_string() const {
    return "STACKCRE";
  }

  template<> long StackOp_component<StackCre>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      std::vector<int> screened_ix;
      if(dmrginp.calc_type() == MPS_NEVPT)
      {
        if(b.nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2))
          screened_ix =  screened_cdd_c_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
        else
          screened_ix =  screened_ccd_c_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);

      }
      else{

        int integralIndex = b.get_integralIndex();
        screened_ix = (dmrginp.hamiltonian() == BCS) ? 
        screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) : 
        screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      }
      m_op.set_indices(screened_ix, dmrginp.last_site());  
      std::vector<int> orbs(1);

      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackCre>(new StackCre);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  op.set_deltaQuantum(1, getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
	  //op.set_deltaQuantum() = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOf(orbs[0]));      
	  op.set_quantum_ladder()["(C)"] = { op.get_deltaQuantum(0) };
	  if (calcMem)
	    requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	}
      return requiredMemory;
    }

  template<> void StackOp_component<StackCre>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_indices(screened_c_ix, dmrginp.last_site());  
      std::vector<int> orbs(1);

      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackCre>(new StackCre);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  op.set_deltaQuantum(1, getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
	  //op.set_deltaQuantum() = SpinQuantum(1, SpinOf(orbs[0]), SymmetryOf(orbs[0]));      
	  op.set_quantum_ladder()["(C)"] = { op.get_deltaQuantum(0) };
	}
    }
  
  
  
  template<> std::vector< std::vector<int> > StackOp_component<StackCre>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }
  

  template<> void StackOp_component<StackCre>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<StackCre> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<StackCre>(new StackCre);
    }


  //usually not needed, because it can be calculated as a transpose of C, but
  //when the bra and ket state in the block are different than transpose cannot be used
  // -------------------- D_S1 ---------------------------  
  template<> string StackOp_component<StackDes>::get_op_string() const {
    return "STACKDES";
  }

  template<> long StackOp_component<StackDes>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      double screen_tol = dmrginp.oneindex_screen_tol();

      std::vector<int> screened_ix;
      if(dmrginp.calc_type() == MPS_NEVPT)
      {
        if(b.nonactive_orb()[0] >=  (dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2))
          screened_ix =  screened_cdd_d_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
        else
          screened_ix =  screened_ccd_d_indices(b.get_sites(),b.get_complementary_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);
      }
      else
      {
      int integralIndex = b.get_integralIndex();
      screened_ix = (dmrginp.hamiltonian() == BCS) ? 
        screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) : 
        screened_d_indices(b.get_sites(), b.get_complementary_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      }
      m_op.set_indices(screened_ix, dmrginp.last_site());  
      std::vector<int> orbs(1);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackDes>(new StackDes);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  op.set_deltaQuantum(1, -getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
	  op.set_quantum_ladder()["(D)"] = { op.get_deltaQuantum(0) };
	  if (calcMem)
	    requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackDes>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_d_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_indices(screened_d_ix, dmrginp.last_site());  
      std::vector<int> orbs(1);

      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackDes>(new StackDes);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  op.set_deltaQuantum(1, -getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));      
	  op.set_quantum_ladder()["(D)"] = { op.get_deltaQuantum(0) };
	}
    }
  
  
  
  template<> std::vector< std::vector<int> > StackOp_component<StackDes>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }
  

  template<> void StackOp_component<StackDes>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<StackDes> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<StackDes>(new StackDes);
    }


  
  // -------------------- Cd_ ---------------------------  
  template<> string StackOp_component<StackCreDes>::get_op_string() const {
    return "STACKCREDES";
  }

  template<> long StackOp_component<StackCreDes>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      int integralIndex = b.get_integralIndex();
      const double screen_tol = dmrginp.twoindex_screen_tol();
      vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ? 
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreDes> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1-spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreDes>(new StackCreDes);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(CD)"] = { op.get_deltaQuantum(0) };
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	  }
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackCreDes>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_cd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreDes> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1-spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreDes>(new StackCreDes);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(CD)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  

  // -------------------- dC_ ---------------------------  
  template<> string StackOp_component<StackDesCre>::get_op_string() const {
    return "STACKDESCRE";
  }

  template<> long StackOp_component<StackDesCre>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      int integralIndex = b.get_integralIndex();
      const double screen_tol = dmrginp.twoindex_screen_tol();
      vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ? 
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
        screened_cd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackDesCre> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackDesCre>(new StackDesCre);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(DC)"] = { op.get_deltaQuantum(0) };
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	  }
	}
      return requiredMemory;
    }
  
  
  template<> void StackOp_component<StackDesCre>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_cd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackDesCre> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackDesCre>(new StackDesCre);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(DC)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  
  // -------------------- Cc_ ---------------------------  
  template<> string StackOp_component<StackCreCre>::get_op_string() const {
    return "STACKCRECRE";
  }
  template<> long StackOp_component<StackCreCre>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      int integralIndex = b.get_integralIndex();
      const double screen_tol = dmrginp.twoindex_screen_tol();
      
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :        
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      long requiredMemory =0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreCre> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreCre>(new StackCreCre);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(CC)"] = { op.get_deltaQuantum(0) };
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	  }
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackCreCre>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_dd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);

      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreCre> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreCre>(new StackCreCre);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(CC)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  
  // -------------------- Dd_ ---------------------------  
  template<> string StackOp_component<StackDesDes>::get_op_string() const {
    return "STACKDESDES";
  }
  template<> long StackOp_component<StackDesDes>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :        
        screened_dd_indices(b.get_sites(), b.get_complementary_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);
      long requiredMemory=0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackDesDes> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = -getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = -getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackDesDes>(new StackDesDes);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(DD)"] = { op.get_deltaQuantum(0) };
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	  }
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackDesDes>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_dd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      std::vector<int> orbs(2);

      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackDesDes> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = -getSpinQuantum(orbs[0]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]));
	  SpinQuantum spin2 = -getSpinQuantum(orbs[1]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[1]));
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackDesDes>(new StackDesDes);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    op.set_deltaQuantum(1, spinvec[j]);      
	    op.set_quantum_ladder()["(DD)"] = { op.get_deltaQuantum(0) };
	  }
	}
    }
  
  
  // -------------------- Cdcomp_ ---------------------------  
  template<> string StackOp_component<StackCreDesComp>::get_op_string() const {
    return "STACKCREDES_COMP";
  }
  template<> long StackOp_component<StackCreDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
  {
    if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.twoindex_screen_tol();
      int integralIndex = b.get_integralIndex();
    vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ?
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
    m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
    
    std::vector<int> orbs(2);
    long requiredMemory = 0;
    for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreDesComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreDesComp>(new StackCreDesComp);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    if (dmrginp.hamiltonian() == BCS) {
	      op.resize_deltaQuantum(3);
	      op.set_deltaQuantum(0) = spinvec[j];
	      op.set_deltaQuantum(1) = SpinQuantum(2, spinvec[j].get_s(), spinvec[j].get_symm());
	      op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	      if (calcMem)
		requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    } else {
	      op.set_deltaQuantum(1, spinvec[j]);
	      if (calcMem)
		requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    }
	  }
	}
    return requiredMemory;
  }

  template<> void StackOp_component<StackCreDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_cd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
    
    std::vector<int> orbs(2);
    for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreDesComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin2-spin1;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreDesComp>(new StackCreDesComp);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    if (dmrginp.hamiltonian() == BCS) {
	      op.resize_deltaQuantum(3);
	      op.set_deltaQuantum(0) = spinvec[j];
	      op.set_deltaQuantum(1) = SpinQuantum(2, spinvec[j].get_s(), spinvec[j].get_symm());
	      op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	    } else {
	      op.set_deltaQuantum(1, spinvec[j]);
	    }
	  }
	}
  }
  
  template<> void StackOp_component<StackCreDesComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<StackCreDesComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin2-spin1;
      vec.resize(spinvec.size());
      vector<int> orbs(2); orbs[0] = i; orbs[1] = j;
      for (int j=0; j<spinvec.size(); j++) {
	vec[j]=boost::shared_ptr<StackCreDesComp>(new StackCreDesComp);
	vec[j]->set_orbs() = orbs;
	vec[j]->set_initialised() = true;
	vec[j]->set_fermion() = false;
	vec[j]->set_deltaQuantum(1, spinvec[j]);
      }
    }
  
  // -------------------- dCcomp_ ---------------------------  
  template<> string StackOp_component<StackDesCreComp>::get_op_string() const {
    return "STACKDESCRE_COMP";
  }
  template<> long StackOp_component<StackDesCreComp>::build_iterators(StackSpinBlock& b, bool calcMem)
  {
    if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.twoindex_screen_tol();
    int integralIndex = b.get_integralIndex();
    vector< pair<int, int> > screened_cd_ix = (dmrginp.hamiltonian() == BCS) ?
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
      screened_cd_indices( b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
    m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
    
    std::vector<int> orbs(2);
    long requiredMemory = 0;
    for (int i = 0; i < m_op.local_nnz(); ++i)
    {
      orbs = m_op.unmap_local_index(i);
      std::vector<boost::shared_ptr<StackDesCreComp> >& vec = m_op.get_local_element(i);
      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
      SpinQuantum spin2 = getSpinQuantum(orbs[1]);
      std::vector<SpinQuantum> spinvec = spin1-spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) {
	vec[j]=boost::shared_ptr<StackDesCreComp>(new StackDesCreComp);
	StackSparseMatrix& op = *vec[j];
	op.set_orbs() = orbs;
	op.set_initialised() = true;
	op.set_fermion() = false;
        if (dmrginp.hamiltonian() == BCS) {
          op.resize_deltaQuantum(3);
          op.set_deltaQuantum(0) = spinvec[j];
          op.set_deltaQuantum(1) = SpinQuantum(2, spinvec[j].get_s(), spinvec[j].get_symm());
          op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	  if (calcMem)
	    requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
        } else {
          op.set_deltaQuantum(1, spinvec[j]);
	  if (calcMem)
	    requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
        }
      }
    }
    return requiredMemory;
  }
  
  template<> void StackOp_component<StackDesCreComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_cd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    m_op.set_pair_indices(screened_cd_ix, dmrginp.last_site());      
    
    std::vector<int> orbs(2);
    for (int i = 0; i < m_op.local_nnz(); ++i)
    {
      orbs = m_op.unmap_local_index(i);
      std::vector<boost::shared_ptr<StackDesCreComp> >& vec = m_op.get_local_element(i);
      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
      SpinQuantum spin2 = getSpinQuantum(orbs[1]);
      std::vector<SpinQuantum> spinvec = spin1-spin2;
      vec.resize(spinvec.size());
      for (int j=0; j<spinvec.size(); j++) {
	vec[j]=boost::shared_ptr<StackDesCreComp>(new StackDesCreComp);
	StackSparseMatrix& op = *vec[j];
	op.set_orbs() = orbs;
	op.set_initialised() = true;
	op.set_fermion() = false;
        if (dmrginp.hamiltonian() == BCS) {
          op.resize_deltaQuantum(3);
          op.set_deltaQuantum(0) = spinvec[j];
          op.set_deltaQuantum(1) = SpinQuantum(2, spinvec[j].get_s(), spinvec[j].get_symm());
          op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
        } else {
          op.set_deltaQuantum(1, spinvec[j]);
        }
      }
    }

  }
  
  
  template<> void StackOp_component<StackDesCreComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<StackDesCreComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1-spin2;
      vec.resize(spinvec.size());
      vector<int> orbs(2); orbs[0] = i; orbs[1] = j;
      for (int j=0; j<spinvec.size(); j++) {
	vec[j]=boost::shared_ptr<StackDesCreComp>(new StackDesCreComp);
	vec[j]->set_orbs() = orbs;
	vec[j]->set_initialised() = true;
	vec[j]->set_fermion() = false;
	vec[j]->set_deltaQuantum(1, spinvec[j]);
      }
    }
  
  // -------------------- Ddcomp_ ---------------------------  
  template<> string StackOp_component<StackDesDesComp>::get_op_string() const {
    return "STACKDESDES_COMP";
  }
  template<> long StackOp_component<StackDesDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      long requiredMemory=0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackDesDesComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackDesDesComp>(new StackDesDesComp);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    
	    if (dmrginp.hamiltonian() == BCS) {
	      op.resize_deltaQuantum(3);          
	      op.set_deltaQuantum(0) = -spinvec[j];
	      op.set_deltaQuantum(1) = -SpinQuantum(0, spinvec[j].get_s(), spinvec[j].get_symm());
	      op.set_deltaQuantum(2) = -SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	      if (calcMem)
		requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    } else {
	      op.set_deltaQuantum(1, -spinvec[j]);
	      if (calcMem)
		requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    }  
	  }
	}
      return requiredMemory;
    }

  template<> void StackOp_component<StackDesDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_dd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackDesDesComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackDesDesComp>(new StackDesDesComp);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    
	    if (dmrginp.hamiltonian() == BCS) {
	      op.resize_deltaQuantum(3);          
	      op.set_deltaQuantum(0) = -spinvec[j];
	      op.set_deltaQuantum(1) = -SpinQuantum(0, spinvec[j].get_s(), spinvec[j].get_symm());
	      op.set_deltaQuantum(2) = -SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	    } else {
	      op.set_deltaQuantum(1, -spinvec[j]);
	    }  
	  }
	}
    }
  
  template<> void StackOp_component<StackDesDesComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<StackDesDesComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1+spin2;
      vec.resize(spinvec.size());
      vector<int> orbs(2); orbs[0] = i; orbs[1] = j;
      for (int j=0; j<spinvec.size(); j++) {
	vec[j]=boost::shared_ptr<StackDesDesComp>(new StackDesDesComp);
	vec[j]->set_orbs() = orbs;
	vec[j]->set_initialised() = true;
	vec[j]->set_fermion() = false;
	vec[j]->set_deltaQuantum(1, -spinvec[j]);
      }
    }
  
  
  // -------------------- CCcomp_ ---------------------------  
  template<> string StackOp_component<StackCreCreComp>::get_op_string() const {
    return "STACKCRECRE_COMP";
  }
  template<> long StackOp_component<StackCreCreComp>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.twoindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      vector< pair<int, int> > screened_dd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
        screened_dd_indices(b.get_complementary_sites(), b.get_sites(), *b.get_twoInt(), screen_tol);
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreCreComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreCreComp>(new StackCreCreComp);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    if (dmrginp.hamiltonian() == BCS) {
	      op.resize_deltaQuantum(3);          
	      op.set_deltaQuantum(0) = spinvec[j];
	      op.set_deltaQuantum(1) = SpinQuantum(0, spinvec[j].get_s(), spinvec[j].get_symm());
	      op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	      if (calcMem)
		requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    } else {
	      op.set_deltaQuantum(1, spinvec[j]);
	      if (calcMem)
		requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    }
	  }
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackCreCreComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_dd_ix, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_pair_indices(screened_dd_ix, dmrginp.last_site());      
      
      std::vector<int> orbs(2);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs = m_op.unmap_local_index(i);
	  std::vector<boost::shared_ptr<StackCreCreComp> >& vec = m_op.get_local_element(i);
	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	  SpinQuantum spin2 = getSpinQuantum(orbs[1]);
	  std::vector<SpinQuantum> spinvec = spin1+spin2;
	  vec.resize(spinvec.size());
	  for (int j=0; j<spinvec.size(); j++) {
	    vec[j]=boost::shared_ptr<StackCreCreComp>(new StackCreCreComp);
	    StackSparseMatrix& op = *vec[j];
	    op.set_orbs() = orbs;
	    op.set_initialised() = true;
	    op.set_fermion() = false;
	    if (dmrginp.hamiltonian() == BCS) {
	      op.resize_deltaQuantum(3);          
	      op.set_deltaQuantum(0) = spinvec[j];
	      op.set_deltaQuantum(1) = SpinQuantum(0, spinvec[j].get_s(), spinvec[j].get_symm());
	      op.set_deltaQuantum(2) = SpinQuantum(-2, spinvec[j].get_s(), spinvec[j].get_symm());
	    } else {
	      op.set_deltaQuantum(1, spinvec[j]);
	    }
	  }
	}
    }
  

  template<> void StackOp_component<StackCreCreComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_indices(i,j);
      
      std::vector<boost::shared_ptr<StackCreCreComp> >& vec = m_op(i,j);
      SpinQuantum spin1 = getSpinQuantum(i);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(i));
      SpinQuantum spin2 = getSpinQuantum(j);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(j));
      std::vector<SpinQuantum> spinvec = spin1+spin2;
      vec.resize(spinvec.size());
      vector<int> orbs(2); orbs[0] = i; orbs[1] = j;
      for (int j=0; j<spinvec.size(); j++) {
	vec[j]=boost::shared_ptr<StackCreCreComp>(new StackCreCreComp);
	vec[j]->set_orbs() = orbs;
	vec[j]->set_initialised() = true;
	vec[j]->set_fermion() = false;
	vec[j]->set_deltaQuantum(1, spinvec[j]);
      }
    }
  
  
  // -------------------- Ccdcomp_ ---------------------------  
  template<> string StackOp_component<StackCreCreDesComp>::get_op_string() const {
    return "STACKCRECREDES_COMP";
  }
  template<> long StackOp_component<StackCreCreDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      vector< int > screened_cdd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      m_op.set_indices(screened_cdd_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackCreCreDesComp>(new StackCreCreDesComp);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  if (dmrginp.hamiltonian() == BCS) {
	    op.resize_deltaQuantum(4);
	    SpinQuantum qorb = getSpinQuantum(orbs[0]);
	    op.set_deltaQuantum(0) = qorb;
	    op.set_deltaQuantum(1) = SpinQuantum(3, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(2) = SpinQuantum(-1, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(3) = SpinQuantum(-3, qorb.get_s(), qorb.get_symm());
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	  } else {
	    op.set_deltaQuantum(1, getSpinQuantum(orbs[0]));
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	  }
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackCreCreDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_cdd_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_indices(screened_cdd_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);

      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackCreCreDesComp>(new StackCreCreDesComp);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  if (dmrginp.hamiltonian() == BCS) {
	    op.resize_deltaQuantum(4);
	    SpinQuantum qorb = getSpinQuantum(orbs[0]);
	    op.set_deltaQuantum(0) = qorb;
	    op.set_deltaQuantum(1) = SpinQuantum(3, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(2) = SpinQuantum(-1, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(3) = SpinQuantum(-3, qorb.get_s(), qorb.get_symm());
	  } else {
	    op.set_deltaQuantum(1, getSpinQuantum(orbs[0]));
	  }
	}
    }
  
  
  template<> std::vector<std::vector<int> > StackOp_component<StackCreCreDesComp>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }

  template<> void StackOp_component<StackCreCreDesComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<StackCreCreDesComp> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<StackCreCreDesComp>(new StackCreCreDesComp);
    }


  //usually not needed, because it can be calculated as a transpose of CCDcomp, but
  //when the bra and ket state in the block are different than transpose cannot be used
  // -------------------- Cddcomp_ ---------------------------  
  template<> string StackOp_component<StackCreDesDesComp>::get_op_string() const {
    return "STACKCREDESDES_COMP";
  }
  template<> long StackOp_component<StackCreDesDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      int integralIndex = b.get_integralIndex();
      vector< int > screened_cdd_ix = (dmrginp.hamiltonian() == BCS) ?
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), v_cc[integralIndex], v_cccc[integralIndex], v_cccd[integralIndex], screen_tol) :
        screened_cddcomp_indices(b.get_complementary_sites(), b.get_sites(), v_1[integralIndex], *b.get_twoInt(), screen_tol);
      m_op.set_indices(screened_cdd_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackCreDesDesComp>(new StackCreDesDesComp);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  if (dmrginp.hamiltonian() == BCS) {
	    op.resize_deltaQuantum(4);
	    SpinQuantum qorb = getSpinQuantum(orbs[0]);
	    op.set_deltaQuantum(0) = -qorb;
	    op.set_deltaQuantum(1) = -SpinQuantum(3, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(2) = -SpinQuantum(-1, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(3) = -SpinQuantum(-3, qorb.get_s(), qorb.get_symm());
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    else
	      requiredMemory = 0;
	  } else {
	    op.set_deltaQuantum(1, -getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
	    if (calcMem)
	      requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	    else
	      requiredMemory = 0;
	  }     
	}
      return requiredMemory;
    }
  
  template<> void StackOp_component<StackCreDesDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_cdd_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      m_op.set_indices(screened_cdd_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);

      for (int i = 0; i < m_op.local_nnz(); ++i)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  m_op.get_local_element(i).resize(1);
	  m_op.get_local_element(i)[0]=boost::shared_ptr<StackCreDesDesComp>(new StackCreDesDesComp);
	  StackSparseMatrix& op = *m_op.get_local_element(i)[0];
	  op.set_orbs() = orbs;
	  op.set_initialised() = true;
	  op.set_fermion() = true;
	  if (dmrginp.hamiltonian() == BCS) {
	    op.resize_deltaQuantum(4);
	    SpinQuantum qorb = getSpinQuantum(orbs[0]);
	    op.set_deltaQuantum(0) = -qorb;
	    op.set_deltaQuantum(1) = -SpinQuantum(3, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(2) = -SpinQuantum(-1, qorb.get_s(), qorb.get_symm());
	    op.set_deltaQuantum(3) = -SpinQuantum(-3, qorb.get_s(), qorb.get_symm());
	  } else {
	    op.set_deltaQuantum(1, -getSpinQuantum(orbs[0]));//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
	  }     
	}
    }
  
  
  template<> std::vector<std::vector<int> > StackOp_component<StackCreDesDesComp>::get_array() const 
    {
      std::vector<int> orbs(1);
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      for (int i=0; i<m_op.local_nnz(); i++)
	{
	  orbs[0] = m_op.get_local_indices()[i];
	  ret_val[i] = orbs;
	}
      return ret_val;
    }

  template<> void StackOp_component<StackCreDesDesComp>::add_local_indices(int i, int j , int k)
    {
      m_op.add_local_index(i);
      
      std::vector<boost::shared_ptr<StackCreDesDesComp> >& vec = m_op(i);
      vec.resize(1);
      vec[0]=boost::shared_ptr<StackCreDesDesComp>(new StackCreDesDesComp);
    }

  
  // -------------------- HAM ---------------------------  
  template<> string StackOp_component<StackHam>::get_op_string() const {
    return "STACKHAM";
  }
  template<> long StackOp_component<StackHam>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackHam>(new StackHam);
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = false;
      long requiredMemory = 0;
      if (dmrginp.hamiltonian() == BCS) {
        m_op(0)[0]->resize_deltaQuantum(5);
        for (int i = 0; i <5; ++i) {
          m_op(0)[0]->set_deltaQuantum(i) = SpinQuantum(2*(i-2), SpinSpace(0), IrrepSpace(0) );
	  if (calcMem)
	    requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), m_op(0)[0]->get_deltaQuantum());
	  else
	    requiredMemory = 0;
	}    
      } 
      else {
        m_op(0)[0]->set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
	if (calcMem)
	  requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), m_op(0)[0]->get_deltaQuantum());
	else
	  requiredMemory = 0;
      }      
      return requiredMemory;
    }

  template<> void StackOp_component<StackHam>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackHam>(new StackHam);
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = false;
      if (dmrginp.hamiltonian() == BCS) {
        m_op(0)[0]->resize_deltaQuantum(5);
        for (int i = 0; i <5; ++i) {
          m_op(0)[0]->set_deltaQuantum(i) = SpinQuantum(2*(i-2), SpinSpace(0), IrrepSpace(0) );
	}    
      } 
      else {
        m_op(0)[0]->set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
      }      
    }
  
  template<> std::vector<std::vector<int> > StackOp_component<StackHam>::get_array() const 
    {
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      return ret_val;
    }
  
  // -------------------- Overlap ---------------------------  
  template<> string StackOp_component<StackOverlap>::get_op_string() const {
    return "STACKOVERLAP";
  }
  template<> long StackOp_component<StackOverlap>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackOverlap>(new StackOverlap);
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = false;

      m_op(0)[0]->set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
      if (calcMem)
	return  SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), m_op(0)[0]->get_deltaQuantum());
      else
	return 0;
    }

  template<> void StackOp_component<StackOverlap>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackOverlap>(new StackOverlap);
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = false;

      m_op(0)[0]->set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
    }
  
  template<> std::vector<std::vector<int> > StackOp_component<StackOverlap>::get_array() const 
    {
      std::vector< std::vector<int> > ret_val(m_op.local_nnz());
      return ret_val;
    }
  


}
