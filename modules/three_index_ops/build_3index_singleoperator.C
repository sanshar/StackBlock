/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "Stack_op_components.h"
#include "StackBaseOperator.h"
#include "Stackspinblock.h"
#include "operatorfunctions.h"
#include "pario.h"
#include <boost/algorithm/string.hpp>

namespace SpinAdapted{
namespace Three_index_ops{

//FIXME take out common parts with 3- and 4- index functions


//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void finish_tensor_trace( const StackSpinBlock& b, StackSpinBlock* sysdot, StackSparseMatrix& sysdot_op, StackSparseMatrix& op, std::string& build_pattern )
{
  op.set_build_pattern() = build_pattern;
  op.set_deltaQuantum(1, op.get_quantum_ladder().at( build_pattern ).at(1) );
  SpinAdapted::operatorfunctions::TensorTrace(sysdot, sysdot_op, &b, &(b.get_stateInfo()), op);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void finish_tensor_product( const StackSpinBlock& b, StackSpinBlock* sysdot, 
                            const StackSparseMatrix& sysdot_op1, const StackSparseMatrix& sysdot_op2, StackSparseMatrix& op, 
                            bool include_parity, std::string& build_pattern )
{
  op.set_build_pattern() = build_pattern;
  op.set_deltaQuantum(1, op.get_quantum_ladder().at( build_pattern ).at(1) );
  // Do tensor product
  double parity = 1.0;
  if ( include_parity ) parity = getCommuteParity( sysdot_op1.get_deltaQuantum(0), sysdot_op2.get_deltaQuantum(0), op.get_deltaQuantum(0) );
  SpinAdapted::operatorfunctions::TensorProduct(sysdot, sysdot_op1, sysdot_op2, &b, &(b.get_stateInfo()), op, parity);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_single_op_tensor_trace( const opTypes& optype, StackSparseMatrix& op, const StackSpinBlock& big, StackSpinBlock* sysdot)
{
  int i = op.get_orbs()[0]; int j = op.get_orbs()[1]; int k = op.get_orbs()[2];
  StackOp_component_base& sysdot_array = sysdot->get_op_array(optype);
  if (!sysdot_array.has(i,j,k)) return;
  
  if ( (! dmrginp.do_npdm_in_core()) && sysdot->size() > 1 ) {
    std::vector<boost::shared_ptr<StackSparseMatrix> > sysdot_ops = sysdot_array.get_element(i,j,k);;
    
    //loop over the sysdotops and pick the one corresponding to op
    for (int jdx=0; jdx < sysdot_ops.size(); jdx++) {
      boost::shared_ptr<StackSparseMatrix>& sysdot_op = sysdot_ops[jdx];

      int len1 = sysdot_op->get_filename().length(),
	        len2 = op.get_filename().length();
      char char1 = sysdot_op->get_filename()[len1-1], char2 = op.get_filename()[len2-1];
      if ( char1 == char2) {
	std::string build_pattern = sysdot_op->get_build_pattern();
        if ( dmrginp.doimplicitTranspose() )
	  finish_tensor_trace( big, sysdot, *sysdot_op, op, build_pattern );
        else{
	  
          const StackSpinBlock* overlap_block = (big.get_leftBlock() == sysdot) ? big.get_rightBlock() : big.get_leftBlock();
          SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
          const boost::shared_ptr<StackSparseMatrix> Overlap = overlap_block->get_op_rep(OVERLAP, hq);
          bool forwards = false;
          finish_tensor_product( big, sysdot, *sysdot_op, *Overlap, op, forwards, build_pattern );
        }
	
	return;
      }
    }
  }
  else {  
    // Loop over all operator indices
    std::vector<boost::shared_ptr<StackSparseMatrix> > sysdot_ops = sysdot_array.get_element(i,j,k);;
    
    //loop over the sysdotops and pick the one corresponding to op
    for (int jdx=0; jdx < sysdot_ops.size(); jdx++) {
      boost::shared_ptr<StackSparseMatrix>& sysdot_op = sysdot_ops[jdx];
      std::string build_pattern = sysdot_op->get_build_pattern();
      
      std::vector<SpinQuantum> s1 = sysdot_op->get_quantum_ladder().at(build_pattern);
      std::vector<SpinQuantum> s2 = op.get_quantum_ladder().at(build_pattern);
      // Store spin component in correct location
      if ( s1 == s2 ) 
	{
	  if ( dmrginp.doimplicitTranspose() )
	    finish_tensor_trace( big, sysdot, *sysdot_op, op, build_pattern );
	  else{
	    const StackSpinBlock* overlap_block = (big.get_leftBlock() == sysdot) ? big.get_rightBlock() : big.get_leftBlock();
	    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	    const boost::shared_ptr<StackSparseMatrix> Overlap = overlap_block->get_op_rep(OVERLAP, hq);
	    bool forwards = false;
	    finish_tensor_product( big, sysdot, *sysdot_op, *Overlap, op, forwards, build_pattern );
	  }
	  
	}
    }
  }
  
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_1_2_single_op_tensor_products( bool forwards, const opTypes& optype, StackSparseMatrix& op, const opTypes& rhsType, 
					      const opTypes& lhsType, const StackSpinBlock& big, StackSpinBlock* rhsBlock, StackSpinBlock* lhsBlock)
{
  // (i | j,k ) partition
  //-------------------------
  StackOp_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  StackOp_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
  assert ( (rhs_array.get_size() == 1) || (lhs_array.get_size() == 1) );

  int io = op.get_orbs()[0]; int jo = op.get_orbs()[1]; int ko = op.get_orbs()[2];

  // Loop over all lhs operator indices
  for (int idx = 0; idx < lhs_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<StackSparseMatrix> > lhs_ops;
    // Assume 2-index operators are available on this processor in core
    lhs_ops = lhs_array.get_local_element(idx);

    // Loop over all rhs operator indices
    for (int iidx = 0; iidx < rhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<StackSparseMatrix> > rhs_ops = rhs_array.get_local_element(iidx);
      int i = rhs_ops[0]->get_orbs()[0];
      int j = lhs_ops[0]->get_orbs()[0]; int k = lhs_ops[0]->get_orbs()[1];

      // In parallel calculations not all operators are built on each proc
      if ( io==i && jo==j && ko==k ) {

	// Loop over lhs spin-op components
	for (int jdx=0; jdx < lhs_ops.size(); jdx++) {
	  boost::shared_ptr<StackSparseMatrix>& lhs_op = lhs_ops[jdx];
	  //assert( lhs_op.get_built() );
	  std::string build_23 = lhs_op->get_build_pattern();
	  std::vector<SpinQuantum> spin_23 = lhs_op->get_quantum_ladder().at(build_23);

	  // Loop over rhs spin-op components
	  for (int jjdx=0; jjdx < rhs_ops.size(); jjdx++) {
	    boost::shared_ptr<StackSparseMatrix>& rhs_op = rhs_ops[jjdx];
	    //assert( rhs_op->get_built() );
	    std::string build_1 = rhs_op->get_build_pattern();
	    std::string build_pattern = "(" + build_1 + build_23 + ")";

            std::vector<SpinQuantum> s = { op.get_quantum_ladder().at(build_pattern).at(0) };
            if ( s == spin_23 ) {
              finish_tensor_product( big, rhsBlock, *rhs_op, *lhs_op, op, forwards, build_pattern );
            }
          }
        }
	break;
      }

    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void do_3index_2_1_single_op_tensor_products( bool forwards, const opTypes& optype, StackSparseMatrix& op, const opTypes& rhsType, const opTypes& lhsType,
					      const StackSpinBlock& big, StackSpinBlock* rhsBlock, StackSpinBlock* lhsBlock)
{
  // (i,j | k ) partition
  //-------------------------
  StackOp_component_base& rhs_array = rhsBlock->get_op_array(rhsType);
  StackOp_component_base& lhs_array = lhsBlock->get_op_array(lhsType);
  assert ( (rhs_array.get_size() == 1) || (lhs_array.get_size() == 1) );

  int io = op.get_orbs()[0]; int jo = op.get_orbs()[1]; int ko = op.get_orbs()[2];

  // Loop over all rhs operator indices
  for (int idx = 0; idx < rhs_array.get_size(); ++idx) {
    std::vector<boost::shared_ptr<StackSparseMatrix> > rhs_ops;
    // Assume 2-index operators are available on this processor in core
    rhs_ops = rhs_array.get_local_element(idx);

    // Loop over all lhs operator indices
    for (int iidx = 0; iidx < lhs_array.get_size(); ++iidx) {
      std::vector<boost::shared_ptr<StackSparseMatrix> > lhs_ops = lhs_array.get_local_element(iidx);
      int i = rhs_ops[0]->get_orbs()[0]; int j = rhs_ops[0]->get_orbs()[1];
      int k = lhs_ops[0]->get_orbs()[0];

      if (io==i && jo==j && ko==k ) {

	// Loop over rhs spin-op components
	for (int jdx=0; jdx < rhs_ops.size(); jdx++) {
	  boost::shared_ptr<StackSparseMatrix>& rhs_op = rhs_ops[jdx];
	  //assert( rhs_op->get_built() );

	  std::string build_12 = rhs_op->get_build_pattern();
	  std::vector<SpinQuantum> spin_12 = rhs_op->get_quantum_ladder().at(build_12);

	  // Loop over lhs spin-op components //FIXME
	  for (int jjdx=0; jjdx < lhs_ops.size(); jjdx++) {
	    boost::shared_ptr<StackSparseMatrix>& lhs_op = lhs_ops[jjdx];
	    //assert( lhs_op->get_built() );
	    std::string build_3 = lhs_op->get_build_pattern();
	    std::string build_pattern = "(" + build_12 + build_3 + ")";

            std::vector<SpinQuantum> s = { op.get_quantum_ladder().at(build_pattern).at(0) };

            if ( s == spin_12 ) {
              finish_tensor_product( big, rhsBlock, *rhs_op, *lhs_op, op, forwards, build_pattern );
            }
          }
        }
	break;
      }
    }
  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void build_3index_single_op( const opTypes& optype, const StackSpinBlock& big, 
			     const opTypes& lhsType1, const opTypes& lhsType2,
			     const opTypes& rhsType1, const opTypes& rhsType2,
			     StackSparseMatrix& op)
{
  //std::vector<SpinQuantum> sq = op->get_quantum_ladder()[op->get_build_pattern()];
  //string ladderop = to_string(sq[0].get_s().getirrep())+ to_string(sq[1].get_s().getirrep()) ;

  // Open filesystem if necessary
  //std::ofstream ofs;
  //string filename = (big.get_op_array(optype).get_filename()+ladderop);
  //if ( (! dmrginp.do_npdm_in_core()) && big.size() > 1 ) ofs.open( filename.c_str(), std::ios::binary );

  StackSpinBlock* sysBlock = big.get_leftBlock();
  StackSpinBlock* dotBlock = big.get_rightBlock();

  // All 3 orbitals on sys or dot block
  do_3index_single_op_tensor_trace( optype, op, big, sysBlock);
  do_3index_single_op_tensor_trace( optype, op, big, dotBlock);

  bool forwards = ! ( sysBlock->get_sites().at(0) > dotBlock->get_sites().at(0) );

  // 2,1 partitioning
  if ( forwards ) {
    do_3index_1_2_single_op_tensor_products( forwards, optype, op, lhsType1, rhsType2, big, dotBlock, sysBlock);
    do_3index_2_1_single_op_tensor_products( forwards, optype, op, lhsType2, rhsType1, big, dotBlock, sysBlock);
  } else {
    do_3index_1_2_single_op_tensor_products( forwards, optype, op, lhsType1, rhsType2, big, sysBlock, dotBlock);
    do_3index_2_1_single_op_tensor_products( forwards, optype, op, lhsType2, rhsType1, big, sysBlock, dotBlock);
  }

  //if ( ! dmrginp.do_npdm_in_core() ) store_single_op_on_disk( ofs, op );
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
}
