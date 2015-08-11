#include "Stackspinblock.h"
#include "StackOperators.h"
#include "Stackwavefunction.h"
#include "stackopxop.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pario.h"

//using namespace operatorfunctions;


/********************************************
Formulas for making hamiltonian matrix while blocking a block with a dot block
********************************************/


void SpinAdapted::stackopxop::cdxcdcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o) {
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();    

  for (int opind=0; opind<opvec1.size(); opind++) { // this is CreDes_{ij}
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor, numthreads);
    if (i != j) {
      //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
      if (otherblock->has(DES_CRECOMP)) {
	op1 = loopblock->get_op_array(DES_CRE).get_element(i,j).at(opind);
	double parity = 1.0;
	if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	  parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
	op2 = otherblock->get_op_array(DES_CRECOMP).get_element(i,j).at(opind);
	SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity, numthreads);
      } else {
	//op2->set_conjugacy('t');
	//op1->set_conjugacy('t');
	SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], factor, numthreads);
	//op2->set_conjugacy('n');op1->set_conjugacy('n');
      }
    }
  }
}

void SpinAdapted::stackopxop::ddxcccomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o)
{
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind); // CC_{ij}
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
      return;
    double factor = 2.0; if (i==j) factor = 1.0;
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(opind);

    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));

    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthreads);

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (otherblock->has(CRE_CRECOMP)) {
      op2 = otherblock->get_op_array(CRE_CRECOMP).get_element(i, j).at(opind);
      op1 = loopblock->get_op_array(DES_DES).get_element(i, j).at(opind);
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor, numthreads);
    } else {
      //StackTransposeview top1 = StackTransposeview(*op1);
      //StackTransposeview top2 = StackTransposeview(*op2);
      
      SpinQuantum sq1 = op1->get_deltaQuantum(0);
      SpinQuantum sq2 = op2->get_deltaQuantum(0);
      double parity2 =TensorOp::getTransposeFactorDD(i, j, sq1.get_s().getirrep(), sq1.get_symm().getirrep());
      parity2*=TensorOp::getTransposeFactorDD(i, j, sq2.get_s().getirrep(), sq2.get_symm().getirrep());
      
      parity *= parity2;
      //op2->set_conjugacy('t');op1->set_conjugacy('t');      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], parity*factor, numthreads);
      //op2->set_conjugacy('n');op1->set_conjugacy('n');
    }
  }
}



void SpinAdapted::stackopxop::cxcddcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o)
{
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
      return;

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES)) {
      boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(CRE_DES_DESCOMP).get_element(i).at(opind);
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthreads);	    

      if (otherblock == b->get_rightBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;

      op1 = loopblock->get_op_array(DES).get_element(i).at(opind);
      op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind);
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthreads);      
    } else {
      //StackTransposeview top1 = StackTransposeview(op1);  // DES_i
      boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(opind);
      
      double scale = 1.0;
      double parity = 1.0;
      
      if (otherblock == b->get_rightBlock()) parity = getCommuteParity(-op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      //op1->set_conjugacy('t');
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], scale*parity, numthreads);
      //op1->set_conjugacy('n');

      // complex conjugate
      //StackTransposeview top2 = StackTransposeview(op2); // CDD_i
      if (otherblock == b->get_leftBlock()) parity = getCommuteParity(-op2->get_deltaQuantum(0), op1->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;
      //op2->set_conjugacy('t');//op1->set_conjugacy('n');
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transpose(*op2), *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity, numthreads);	    
      //op2->set_conjugacy('n');//op1->set_conjugacy('n');
    }
  }
}

//***************************************************************************************





/********************************************
Formulas for multiplying hamiltonian with wavefunction without ever making the hamiltonian explicitly
********************************************/

void SpinAdapted::stackopxop::cdxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  if (otherblock->has(DES_CRECOMP)) {
    SpinQuantum opq = op1->get_deltaQuantum()[0];
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    {
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      
      if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
	return;
      boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
	
      double factor = 1.0;
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor);

      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    if (i != j) {
      boost::shared_ptr<StackSparseMatrix> op1t = loopblock->get_op_rep(DES_CRE, -opq, i, j);
      double parity = 1.0;
      if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1t->get_deltaQuantum()[0]);
      boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(DES_CRECOMP, opq, i, j);
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, parity);
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  }
  else {
    SpinQuantum opq = op1->get_deltaQuantum()[0];
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
    
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor);
    if (i != j) {
      //op2->set_conjugacy('t');op1->set_conjugacy('t');
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, Transpose(*op2), Transpose(*op1), b, c, v, hq, factor);
      //op2->set_conjugacy('n');op1->set_conjugacy('n');
    }
    
    if (deallocate2) op2->deallocate();
    if (deallocate1) op1->deallocate();
    
  }
}

void SpinAdapted::stackopxop::ddxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  int j = op1->get_orbs(1);
  if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
    return;
  double factor = 2.0; if (i==j) factor = 1.0;
  boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(DES_DESCOMP, -opq, i, j);
  
  
  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  
  if (otherblock->has(CRE_CRECOMP)) {
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
    {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor*parity);
      
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    {
      op2 = otherblock->get_op_rep(CRE_CRECOMP, opq, i, j);
      op1 = loopblock->get_op_rep(DES_DES, -opq, i, j);
      
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor*parity);
      
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  }
  else {
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
    
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
    
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor*parity);
    
    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    //StackTransposeview top1 = StackTransposeview(*op1);
    //StackTransposeview top2 = StackTransposeview(*op2);
    
    SpinQuantum sq1 = op1->get_deltaQuantum(0);
    SpinQuantum sq2 = op2->get_deltaQuantum(0);
    double parity2 =TensorOp::getTransposeFactorDD(i, j, sq1.get_s().getirrep(), sq1.get_symm().getirrep());
    parity2*=TensorOp::getTransposeFactorDD(i, j, sq2.get_s().getirrep(), sq2.get_symm().getirrep());
    
    parity *= parity2;
    //op2->set_conjugacy('t');op1->set_conjugacy('t');
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, Transpose(*op2), Transpose(*op1), b, c, v, hq, factor*parity);
    //op2->set_conjugacy('n');op1->set_conjugacy('n');
    
    if (deallocate2) op2->deallocate();
    if (deallocate1) op1->deallocate();
  }
}




void SpinAdapted::stackopxop::cxcddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));  // in get_parity, number part is not used
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
    return;
  
  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
  if (loopblock->has(DES) ) {
    double scale = 1.0;
    double parity = 1.0;
    {
      boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(CRE_DES_DESCOMP, -opq, i);
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      
      if (otherblock == b->get_leftBlock())
	parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, scale*parity);	    
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    {
      boost::shared_ptr<StackSparseMatrix> op1 = loopblock->get_op_rep(DES, -opq, i);
      boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      
      if (otherblock == b->get_rightBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
      else parity = 1.0;
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, scale*parity);
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  } else {
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
    //StackTransposeview top1 = StackTransposeview(op1);  // DES_i
    
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(-op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
    
    //op1->set_conjugacy('t');
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, Transpose(*op1), b, c, v, hq, scale*parity);
    //op1->set_conjugacy('n');
    
    // complex conjugate
    //StackTransposeview top2 = StackTransposeview(op2); // CDD_i
    if (otherblock == b->get_leftBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), hq);
    else parity = 1.0;
    
    //op2->set_conjugacy('t');//op1->set_conjugacy('n');
    SpinAdapted::operatorfunctions::TensorMultiply(otherblock, Transpose(*op2), *op1, b, c, v, hq, scale*parity);	    
    //op2->set_conjugacy('n');//op1->set_conjugacy('n');
    if (deallocate2) op2->deallocate();
    if (deallocate1) op1->deallocate();
    
  }
}



void SpinAdapted::stackopxop::hamandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q, double scale)
{
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));  // in get_parity, number part is not used
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
  op1->build(*loopblock);
      
  boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(OVERLAP, hq);
  bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
  op2->build(*otherblock);
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, scale);	    
  if (deallocate2) op2->deallocate();


  op2 = otherblock->get_op_rep(HAM, hq);
  deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
  op2->build(*otherblock);
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, 1.0);	    
  if (deallocate2) op2->deallocate();
  if (deallocate1) op1->deallocate();

  op1 = loopblock->get_op_rep(HAM, hq);
  deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
  op1->build(*loopblock);

  op2 = otherblock->get_op_rep(OVERLAP, hq);
  deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
  op2->build(*otherblock);
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, 1.0);	    
  if (deallocate2) op2->deallocate();
  if (deallocate1) op1->deallocate();
}


//***************************************************************************************************

/********************************************
Formulas for making diagonal hamiltonian matrix while blocking system and environment blocks
********************************************/


void SpinAdapted::stackopxop::cdxcdcomp_d(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, DiagonalMatrix* e)
{
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    //boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind)->getworkingrepresentation(loopblock);
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind);
    double factor = 1.0;
    {
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);

      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), e[ilock], factor);
      if (i != j) {
	//op2->set_conjugacy('t');op1->set_conjugacy('t');
	SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), e[ilock], factor);
	//op2->set_conjugacy('n');op1->set_conjugacy('n'); 
     }

      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  }
}


//************************************************************************************



/********************************************
Formulas for making CCdcomp operators while blocking a block with a dot block
********************************************/

void SpinAdapted::stackopxop::cxcdcomp(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double scale)
{
  int ilock = 0;//omp_get_thread_num();
  int numthreads = 1;
  const StackSpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  if (opvec1[0]->get_orbs(0) >= I) // opvec1 is CRE
    {
      for (int opind=0; opind<opvec1.size(); opind++) {    
	boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
	if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(op1->get_orbs(0), I)) // we have c_J d_I *c_K d_L
	  return;
	
	const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherBlock->get_op_array(CRE_DESCOMP).get_element(op1->get_orbs(0), I); // CD_comp(j,i) have multiple matrices because of spin adaption
	for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	  boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2); // CD
	  vector<SpinQuantum> op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum(); // o is the resulted CCD
	  int j2 = op2q[0].get_s().getirrep(), j1 = op1q[0].get_s().getirrep(), j21 = oq[0].get_s().getirrep();
	  int l2 = op2q[0].get_symm().getirrep(), l1 = op1q[0].get_symm().getirrep(), l21 = oq[0].get_symm().getirrep(), l3 = (-SymmetryOfOrb(I)).getirrep();
	  double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	  if (NonabelianSym)
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	  
	  double parity = 1.0;
	  if (otherBlock == b->get_rightBlock())
	    parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
	  factor*= parity;
	  
	  SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale, numthreads); // CD*C
	}
      }
    } 
  else {
    for (int opind=0; opind<opvec1.size(); opind++) {    
      boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
      if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(I, op1->get_orbs(0)))
	return;
      const std::vector<boost::shared_ptr<StackSparseMatrix>>& opvec2 = otherBlock->has(DES_CRECOMP) ? \
	otherBlock->get_op_array(DES_CRECOMP).get_element(I, op1->get_orbs(0)) : \
	otherBlock->get_op_array(CRE_DESCOMP).get_element(I, op1->get_orbs(0));
      
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
	vector<SpinQuantum> op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum();
	int j2 = op2q[0].get_s().getirrep(), j1 = op1q[0].get_s().getirrep(), j21 = oq[0].get_s().getirrep();
	int l2 = (-op2q[0].get_symm()).getirrep(), l1 = op1q[0].get_symm().getirrep(), l21 = oq[0].get_symm().getirrep(), l3 = (-SymmetryOfOrb(I)).getirrep();
	double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((1+1+0+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	if (NonabelianSym)
	  factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	
	double parity = 1.0;
	if (otherBlock == b->get_rightBlock() && !otherBlock->has(DES_CRECOMP))
	  parity *= getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
	else if (otherBlock == b->get_rightBlock())
	  parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
	
	parity *= TensorOp::getTransposeFactorCD(I, op1->get_orbs(0), j2, l2);
	
	//If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
	if (!otherBlock->has(DES_CRECOMP)) {
	  //op2->set_conjugacy('t');
	  SpinAdapted::operatorfunctions::TensorProduct(otherBlock, Transpose(*op2), *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale*parity, numthreads);
	  //op2->set_conjugacy('n');
        } 
	else {
	  double parity = 1.0;
	  if (otherBlock == b->get_rightBlock())
	    parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
	  SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale*parity, numthreads);
	}
      }
    }
  }
}

void SpinAdapted::stackopxop::dxcccomp(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int K, StackSparseMatrix* o, double scale)
{ 
  int ilock = 0;//omp_get_thread_num();
  int numthreads = 1;
  //int numthreads = dmrginp.thrds_per_node()[mpigetrank()];
  const StackSpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES) ) {
      boost::shared_ptr<StackSparseMatrix> op1 = loopblock->get_op_array(DES).get_element(opvec1.at(opind)->get_orbs(0)).at(opind); // DES_j
      
      bool transpose = false;
      int k = K, i = op1->get_orbs(0); // P_{ij}=-P_{ji} so only one of them is stored --- P_{ij} where i>j
      if (k < i) { 
        k=i; i=K; transpose = true;
      }
      SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
      
      if (!otherBlock->get_op_array(CRE_CRECOMP).has_local_index(k,i))
	    return;
      
      const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherBlock->get_op_array(CRE_CRECOMP).get_element(k,i); // P_{ki}
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	    boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);  // P_{ki}^\dagger
	    
	    SpinQuantum op2q = op2->get_deltaQuantum(0), op1q = op1->get_deltaQuantum(0), oq = o->get_deltaQuantum(0);
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfOrb(K)).getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    if (NonabelianSym)
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	    
	    double parity = 1.0;

	    if (!transpose)  //this is because when you go from CC_{ij} to CC_{ji} there is a phase factor
	      parity *= getCommuteParity(iq, kq, op2->get_deltaQuantum(0)); 
	    if (loopblock == b->get_leftBlock()) //this is because you have CC_{ji} d_j 
	      parity*= getCommuteParity(-iq, op2->get_deltaQuantum(0), kq); 

	    SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthreads);
      }
    } else {
      boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind); // CRE_j
      
      bool transpose = false;
      int k = K, i = op1->get_orbs(0); // P_{ij}=-P_{ji} so only one of them is stored --- P_{ij} where i>j
      if (k < i) { 
        k=i; i=K; transpose = true;
      }
      SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
      
      if (!otherBlock->get_op_array(DES_DESCOMP).has_local_index(k,i))
	    return;
      
      const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherBlock->get_op_array(DES_DESCOMP).get_element(k,i); // P_{ki}
      for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
	//StackTransposeview top = StackTransposeview(opvec2.at(opind2));  // P_{ki}^\dagger
	boost::shared_ptr<StackSparseMatrix> op = opvec2.at(opind2);
	//op->set_conjugacy('t');
	    
	SpinQuantum op2q = -op->get_deltaQuantum(0), op1q = -op1->get_deltaQuantum(0), oq = o->get_deltaQuantum(0);
	int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = (-SymmetryOfOrb(K)).getirrep();
	double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	if (NonabelianSym)
	  factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
	
	double parity = 1.0;
	//this is because DD_ij^dag = - CC_ij for spin 0 and certain spatial irreps
	parity*=TensorOp::getTransposeFactorDD(K, op1->get_orbs(0), j2, l2);
	if (transpose)  //this is because when you go from CC_{ij} to CC_{ji} there is a phase factor
	  parity *= getCommuteParity(iq, kq, op->get_deltaQuantum(0)); 
	if (loopblock == b->get_leftBlock()) //this is because you have CC_{ji} d_j 
	  parity*= getCommuteParity(-iq, op->get_deltaQuantum(0), kq); 
	
	    //pout << k<<"  "<<i<<"  "<<factor<<"  "<<scale<<"  "<<parity<<endl;
	//op1->set_conjugacy('t');
	SpinAdapted::operatorfunctions::TensorProduct(otherBlock, Transpose(*op), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthreads);
	//op->set_conjugacy('n');op1->set_conjugacy('n');
      }
    }
  }
}


/********************************************
Formulas for making Cddcomp operators while blocking a block with a dot block
********************************************/

void SpinAdapted::stackopxop::dxcdcomp(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double scale)
{
  int ilock = 0;//omp_get_thread_num();
  int numthreads = 1;
  const StackSpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {    
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
    if (!otherBlock->get_op_array(CRE_DESCOMP).has_local_index(max(I, op1->get_orbs(0)), min(I, op1->get_orbs(0)))) // we have c_J d_I *c_K d_L
      return;
    
    std::vector<boost::shared_ptr<StackSparseMatrix> > opvec2 ; 
    if (I > opvec1[0]->get_orbs(0) ) // opvec1 is DES
      opvec2 = otherBlock->get_op_array(CRE_DESCOMP).get_element(I, op1->get_orbs(0)); 
    else
      opvec2 = otherBlock->get_op_array(DES_CRECOMP).get_element(op1->get_orbs(0), I); 
    
    for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2); // CD
      vector<SpinQuantum> op2q = op2->get_deltaQuantum(), op1q = op1->get_deltaQuantum(), oq = o->get_deltaQuantum(); // o is the resulted CCD
      int j2 = op2q[0].get_s().getirrep(), j1 = op1q[0].get_s().getirrep(), j21 = oq[0].get_s().getirrep();
      int l2 = op2q[0].get_symm().getirrep(), l1 = op1q[0].get_symm().getirrep(), l21 = oq[0].get_symm().getirrep(), l3 = (-SymmetryOfOrb(I)).getirrep();
      double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
      if (NonabelianSym)
      factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
      
      double parity = 1.0;
      if (otherBlock == b->get_leftBlock())
	parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
      
      SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*scale*parity, numthreads); // CD*D
    }
  
  } 
}

void SpinAdapted::stackopxop::cxddcomp(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int K, StackSparseMatrix* o, double scale)
{ 
  int ilock = 0;//omp_get_thread_num();
  int numthreads = 1;
  //int numthreads = dmrginp.thrds_per_node()[mpigetrank()];
  const StackSpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {    
    
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind); // CRE_j
    
    bool transpose = false;
    int k = K, i = op1->get_orbs(0); // P_{ij}=-P_{ji} so only one of them is stored --- P_{ij} where i>j
    if (k < i) 
      { k=i; i=K; transpose = true;}
    SpinQuantum iq = getSpinQuantum(i), kq = getSpinQuantum(k);
    
    if (!otherBlock->get_op_array(DES_DESCOMP).has_local_index(k,i))
      return;
    
    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherBlock->get_op_array(DES_DESCOMP).get_element(k,i); // P_{ki}
    for (int opind2 = 0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);  // P_{ki}^\dagger
      
      SpinQuantum op2q = op2->get_deltaQuantum(0), op1q = op1->get_deltaQuantum(0), oq = o->get_deltaQuantum(0);
      int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
      int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = (-oq.get_symm()).getirrep(), l3 = (-SymmetryOfOrb(K)).getirrep();
      double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
      if (NonabelianSym)
	factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
      
      double parity = 1.0;
      
      if (transpose)  //this is because when you go from CC_{ij} to CC_{ji} there is a phase factor
	parity *= getCommuteParity(iq, kq, -op2->get_deltaQuantum(0)); 
      if (loopblock == b->get_rightBlock()) //this is because you have CC_{ij} d_j 
	parity*= getCommuteParity(iq, op2->get_deltaQuantum(0), -kq); 

      SpinAdapted::operatorfunctions::TensorProduct(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], parity*factor*scale, numthreads);

    }
  }
}


//**********************************************************************************************************

