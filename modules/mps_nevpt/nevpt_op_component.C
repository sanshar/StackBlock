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

  
  // -------------------- CDD_sum ---------------------------  
  template<> string StackOp_component<StackCDD_sum>::get_op_string() const {
    return "CDD_SUM";
  }
  template<> long StackOp_component<StackCDD_sum>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackCDD_sum>(new StackCDD_sum);
      //m_op(0)[0]->set_orbs() = std::vector<int>(1,b.nonactive_orb(0));
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = true;
      m_op(0)[0]->set_deltaQuantum(1, -getSpinQuantum(b.nonactive_orb(0)));
	    if (calcMem)
	      return SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), m_op(0)[0]->get_deltaQuantum());
	    else
	      return 0;
    }

  template<> void StackOp_component<StackCDD_sum>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackCDD_sum>(new StackCDD_sum);
      //m_op(0)[0]->set_orbs() = std::vector<int>(1,b.nonactive_orb(0));
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = true;
      m_op(0)[0]->set_deltaQuantum(1, -getSpinQuantum(b.nonactive_orb(0)));
    }
  
  // -------------------- CDD_CreDescomp_ ---------------------------  
  template<> string StackOp_component<StackCDD_CreDesComp>::get_op_string() const {
    return "CDD_CREDES_COMP";
  }
  template<> long StackOp_component<StackCDD_CreDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0 ; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      std::vector<int> screened_d_ix = screened_cdd_d_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
      m_op.set_indices(screened_d_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	    {
	      orbs[0] = m_op.get_local_indices()[i];
	      std::vector<boost::shared_ptr<StackCDD_CreDesComp> >& vec = m_op.get_local_element(i);
	      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	      SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
	      std::vector<SpinQuantum> spinvec = spin2-spin1;
        vec.resize(spinvec.size());
        for(int j=0; j < spinvec.size(); j++)
        {
	        vec[j]=boost::shared_ptr<StackCDD_CreDesComp>(new StackCDD_CreDesComp);
	        StackSparseMatrix& op = *m_op.get_local_element(i)[j];
	        op.set_orbs() = orbs;
	        op.set_initialised() = true;
	        op.set_fermion() = false;
	        op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
	        if (calcMem)
	          requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	        else
	          requiredMemory = 0;
        }
	    }
    return requiredMemory;
    }

  template<> void StackOp_component<StackCDD_CreDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    //screen_tol = 0.0;
      m_op.set_indices(screened_c_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      for (int i = 0; i < m_op.local_nnz(); ++i)
	    {
	      orbs[0] = m_op.get_local_indices()[i];
	      std::vector<boost::shared_ptr<StackCDD_CreDesComp> >& vec = m_op.get_local_element(i);
	      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	      SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
	      std::vector<SpinQuantum> spinvec = spin2-spin1;
        vec.resize(spinvec.size());
        for(int j=0; j < spinvec.size(); j++)
        {
	        vec[j]=boost::shared_ptr<StackCDD_CreDesComp>(new StackCDD_CreDesComp);
	        StackSparseMatrix& op = *m_op.get_local_element(i)[j];
	        op.set_orbs() = orbs;
	        op.set_initialised() = true;
	        op.set_fermion() = false;
	        op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
        }
	    }
    }
  
  // -------------------- CDD_DesDescomp_ ---------------------------  
  template<> string StackOp_component<StackCDD_DesDesComp>::get_op_string() const {
    return "CDD_DESDES_COMP";
  }

  template<> long StackOp_component<StackCDD_DesDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
  {
    if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.oneindex_screen_tol();
    //screen_tol = 0.0;
    std::vector<int> screened_c_ix = screened_cdd_c_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Va], screen_tol);
    m_op.set_indices(screened_c_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    long requiredMemory = 0;
    for (int i = 0; i < m_op.local_nnz(); ++i)
  	{
  	  orbs[0] = m_op.get_local_indices()[i];
  	  std::vector<boost::shared_ptr<StackCDD_DesDesComp> >& vec = m_op.get_local_element(i);
  	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
  	  SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
  	  std::vector<SpinQuantum> spinvec = spin2+spin1;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
  	    vec[j]=boost::shared_ptr<StackCDD_DesDesComp>(new StackCDD_DesDesComp);
  	    StackSparseMatrix& op = *m_op.get_local_element(i)[j];
  	    op.set_orbs() = orbs;
  	    op.set_initialised() = true;
  	    op.set_fermion() = false;
  	    op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
	        if (calcMem)
	          requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	        else
	          requiredMemory = 0;
      }
  	}
    return requiredMemory;
  }
  
  template<> void StackOp_component<StackCDD_DesDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    m_op.set_indices(screened_c_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    for (int i = 0; i < m_op.local_nnz(); ++i)
  	{
  	  orbs[0] = m_op.get_local_indices()[i];
  	  std::vector<boost::shared_ptr<StackCDD_DesDesComp> >& vec = m_op.get_local_element(i);
  	  SpinQuantum spin1 = getSpinQuantum(orbs[0]);
  	  SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
  	  std::vector<SpinQuantum> spinvec = spin2+spin1;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
  	    vec[j]=boost::shared_ptr<StackCDD_DesDesComp>(new StackCDD_DesDesComp);
  	    StackSparseMatrix& op = *m_op.get_local_element(i)[j];
  	    op.set_orbs() = orbs;
  	    op.set_initialised() = true;
  	    op.set_fermion() = false;
  	    op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
      }
  	}
  }
  
  // -------------------- CCD_sum ---------------------------  
  template<> string StackOp_component<StackCCD_sum>::get_op_string() const {
    return "CCD_SUM";
  }
  template<> long StackOp_component<StackCCD_sum>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackCCD_sum>(new StackCCD_sum);
      //m_op(0)[0]->set_orbs() = std::vector<int>(1,b.nonactive_orb(0));
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = true;
      m_op(0)[0]->set_deltaQuantum(1, getSpinQuantum(b.nonactive_orb(0)));
	    if (calcMem)
	      return SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), m_op(0)[0]->get_deltaQuantum());
	    else
	      return 0;
    }
  
  template<> void StackOp_component<StackCCD_sum>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
    {
      m_op.set_indices();
      m_op(0).resize(1);
      m_op(0)[0]=boost::shared_ptr<StackCCD_sum>(new StackCCD_sum);
      //m_op(0)[0]->set_orbs() = std::vector<int>(1,b.nonactive_orb(0));
      m_op(0)[0]->set_orbs() = std::vector<int>();
      m_op(0)[0]->set_initialised() = true;
      m_op(0)[0]->set_fermion() = true;
      m_op(0)[0]->set_deltaQuantum(1, getSpinQuantum(b.nonactive_orb(0)));
    }
  
  // -------------------- CCD_CreDescomp_ ---------------------------  
  template<> string StackOp_component<StackCCD_CreDesComp>::get_op_string() const {
    return "CCD_CREDES_COMP";
  }
  template<> long StackOp_component<StackCCD_CreDesComp>::build_iterators(StackSpinBlock& b, bool calcMem)
    {
      if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
      const double screen_tol = dmrginp.oneindex_screen_tol();
      std::vector<int> screened_c_ix = screened_ccd_c_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);

      m_op.set_indices(screened_c_ix, dmrginp.last_site());      
      std::vector<int> orbs(1);
      long requiredMemory = 0;
      for (int i = 0; i < m_op.local_nnz(); ++i)
	    {
	      orbs[0] = m_op.get_local_indices()[i];
	      std::vector<boost::shared_ptr<StackCCD_CreDesComp> >& vec = m_op.get_local_element(i);
	      SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	      SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
	      std::vector<SpinQuantum> spinvec = spin1-spin2;
        vec.resize(spinvec.size());
        for(int j=0; j < spinvec.size(); j++)
        {
	        vec[j]=boost::shared_ptr<StackCCD_CreDesComp>(new StackCCD_CreDesComp);
	        StackSparseMatrix& op = *m_op.get_local_element(i)[j];
	        op.set_orbs() = orbs;
	        op.set_initialised() = true;
	        op.set_fermion() = false;
	        op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
	        if (calcMem)
	          requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	        else
	          requiredMemory = 0;
        }
	    }
    return requiredMemory;
    }
  
  template<> void StackOp_component<StackCCD_CreDesComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    m_op.set_indices(screened_c_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    for (int i = 0; i < m_op.local_nnz(); ++i)
	  {
	    orbs[0] = m_op.get_local_indices()[i];
	    std::vector<boost::shared_ptr<StackCCD_CreDesComp> >& vec = m_op.get_local_element(i);
	    SpinQuantum spin1 = getSpinQuantum(orbs[0]);
	    SpinQuantum spin2 = getSpinQuantum(b.nonactive_orb(0));
	    std::vector<SpinQuantum> spinvec = spin1-spin2;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
	      vec[j]=boost::shared_ptr<StackCCD_CreDesComp>(new StackCCD_CreDesComp);
	      StackSparseMatrix& op = *m_op.get_local_element(i)[j];
	      op.set_orbs() = orbs;
	      op.set_initialised() = true;
	      op.set_fermion() = false;
	      op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
      }
	  }
  }
  

  // -------------------- CCD_CreCrecomp_ ---------------------------  
  template<> string StackOp_component<StackCCD_CreCreComp>::get_op_string() const {
    return "CCD_CRECRE_COMP";
  }
  template<> long StackOp_component<StackCCD_CreCreComp>::build_iterators(StackSpinBlock& b, bool calcMem)
  {
    if (b.get_sites().size () == 0) return 0; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    const double screen_tol = dmrginp.oneindex_screen_tol();
    int integralIndex = b.get_integralIndex();
    std::vector<int> screened_d_ix = screened_ccd_d_indices(b.get_complementary_sites(), b.get_sites(), b.nonactive_orb()[0], vpt_1, vpt_2[Vi], screen_tol);
    m_op.set_indices(screened_d_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    long requiredMemory = 0;
    for (int i = 0; i < m_op.local_nnz(); ++i)
  	{
  	  orbs[0] = m_op.get_local_indices()[i];
  	  std::vector<boost::shared_ptr<StackCCD_CreCreComp> >& vec = m_op.get_local_element(i);
  	  SpinQuantum spin1 = -getSpinQuantum(orbs[0]);
  	  SpinQuantum spin2 = -getSpinQuantum(b.nonactive_orb(0));
  	  std::vector<SpinQuantum> spinvec = spin2+spin1;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
  	    vec[j]=boost::shared_ptr<StackCCD_CreCreComp>(new StackCCD_CreCreComp);
  	    StackSparseMatrix& op = *m_op.get_local_element(i)[j];
  	    op.set_orbs() = orbs;
  	    op.set_initialised() = true;
  	    op.set_fermion() = false;
  	    op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
	      if (calcMem)
	        requiredMemory += SpinAdapted::getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), op.get_deltaQuantum());
	      else
	        requiredMemory = 0;
      }
  	}
    return requiredMemory;
  }
  
  template<> void StackOp_component<StackCCD_CreCreComp>::build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)
  {
    if (b.get_sites().size () == 0) return; // blank construction (used in unset_initialised() Block copy construction, for use with STL)
    m_op.set_indices(screened_c_ix, dmrginp.last_site());      
    std::vector<int> orbs(1);
    for (int i = 0; i < m_op.local_nnz(); ++i)
  	{
  	  orbs[0] = m_op.get_local_indices()[i];
  	  std::vector<boost::shared_ptr<StackCCD_CreCreComp> >& vec = m_op.get_local_element(i);
  	  SpinQuantum spin1 = -getSpinQuantum(orbs[0]);
  	  SpinQuantum spin2 = -getSpinQuantum(b.nonactive_orb(0));
  	  std::vector<SpinQuantum> spinvec = spin2+spin1;
      vec.resize(spinvec.size());
      for(int j=0; j < spinvec.size(); j++)
      {
  	    vec[j]=boost::shared_ptr<StackCCD_CreCreComp>(new StackCCD_CreCreComp);
  	    StackSparseMatrix& op = *m_op.get_local_element(i)[j];
  	    op.set_orbs() = orbs;
  	    op.set_initialised() = true;
  	    op.set_fermion() = false;
  	    op.set_deltaQuantum(1, -spinvec[j]);//SpinQuantum(1, 1, SymmetryOfSpatialOrb(orbs[0]) );
      }
  	}
  }
  
}
