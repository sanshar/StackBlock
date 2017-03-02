#include <newmat.h>
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
#include "MatrixBLAS.h"

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
	      if (dmrginp.spinAdapted()  && dmrginp.hamiltonian() != BCS)
	        parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
	      op2 = otherblock->get_op_array(DES_CRECOMP).get_element(i,j).at(opind);
	      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity, numthreads);
      } else {
	      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], factor, numthreads);
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
    //pout << "building ham dd "<<i<<"  "<<j<<endl;
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

    //pout << "building ham ccd "<<i<<endl;
    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES) && otherblock->has(CRE_DES_DESCOMP)) {
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

void SpinAdapted::stackopxop::cdxcdcomp_Element(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o, StackMatrix& m,int row, int col) {
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();    
  for (int opind=0; opind<opvec1.size(); opind++) { // this is CreDes_{ij}
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    //pout << "building ham cd "<<i<<"  "<<j<<endl;
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(opind);
    double factor = 1.0;
    SpinAdapted::operatorfunctions::TensorProductElement(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, factor);
    if (i != j) {
      //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
      if (otherblock->has(DES_CRECOMP)) {
	pout << "NOT YET IMPLEMENTED "<<endl;exit(0);
	op1 = loopblock->get_op_array(DES_CRE).get_element(i,j).at(opind);
	double parity = 1.0;
	if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	  parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
	op2 = otherblock->get_op_array(DES_CRECOMP).get_element(i,j).at(opind);
	SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity, numthreads);
      } else {
	SpinAdapted::operatorfunctions::TensorProductElement(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], m, row, col, factor);
      }
    }
  }
}

void SpinAdapted::stackopxop::ddxcccomp_Element(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o, StackMatrix& m, int row, int col)
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
    //pout << "building ham dd "<<i<<"  "<<j<<endl;
    double factor = 2.0; if (i==j) factor = 1.0;
    boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(opind);

    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));

    SpinAdapted::operatorfunctions::TensorProductElement(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, parity*factor);

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (otherblock->has(CRE_CRECOMP)) {
      pout << "NOT YET IMPLEMENTED "<<endl;exit(0);
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
      SpinAdapted::operatorfunctions::TensorProductElement(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], m, row, col, parity*factor);
      //op2->set_conjugacy('n');op1->set_conjugacy('n');
    }
  }
}

void SpinAdapted::stackopxop::cxcddcomp_Element(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o, StackMatrix& m, int row, int col)
{
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind);
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CRE_CRE_DESCOMP).has_local_index(i))
      return;

    //pout << "building ham ccd "<<i<<endl;
    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES) && otherblock->has(CRE_DES_DESCOMP)) {
      pout << "NOT YET IMPLEMENTED "<<endl;exit(0);
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
      SpinAdapted::operatorfunctions::TensorProductElement(otherblock, *op2, Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], m, row, col, scale*parity);
      //op1->set_conjugacy('n');

      // complex conjugate
      //StackTransposeview top2 = StackTransposeview(op2); // CDD_i
      if (otherblock == b->get_leftBlock()) parity = getCommuteParity(-op2->get_deltaQuantum(0), op1->get_deltaQuantum(0), o->get_deltaQuantum(0));
      else parity = 1.0;
      //op2->set_conjugacy('t');//op1->set_conjugacy('n');
      SpinAdapted::operatorfunctions::TensorProductElement(otherblock, Transpose(*op2), *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, scale*parity);	    
      //op2->set_conjugacy('n');//op1->set_conjugacy('n');
    }
  }
}


//***************************************************************************************





/********************************************
Formulas for multiplying hamiltonian with wavefunction without ever making the hamiltonian explicitly
********************************************/

void SpinAdapted::stackopxop::ddxcccomp_3index(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  dmrginp.cctime->start();
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  int j = op1->get_orbs(1);
  if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
    return;
  double factor = 2.0; if (i==j) factor = 1.0;
  boost::shared_ptr<StackSparseMatrix> op2;
  if (dmrginp.spinAdapted() && dmrginp.hamiltonian() != BCS)
    op2 = otherblock->get_op_rep(DES_DESCOMP, -opq, i, j);
  else
    op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(0);
  
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
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) {
	      op1 = loopblock->get_op_array(DES_DES).get_element(i, j).at(0);
	      op2 = otherblock->get_op_array(CRE_CRECOMP).get_element(i, j).at(0);
      }
      else {
	      op1 = loopblock->get_op_rep(DES_DES, -opq, i, j);
	      op2 = otherblock->get_op_rep(CRE_CRECOMP, opq, i, j);
      }
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      
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
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
    SpinQuantum sq1 = op1->get_deltaQuantum(0);
    SpinQuantum sq2 = op2->get_deltaQuantum(0);
    double parity2 =TensorOp::getTransposeFactorDD(i, j, sq1.get_s().getirrep(), sq1.get_symm().getirrep());
    parity2*=TensorOp::getTransposeFactorDD(i, j, sq2.get_s().getirrep(), sq2.get_symm().getirrep());    
    parity2 *= parity;

    
    //SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor*parity);        
    //SpinAdapted::operatorfunctions::TensorMultiply(otherblock, Transpose(*op2), Transpose(*op1), b, c, v, hq, factor*parity2);

    const std::vector<int>& lsites = loopblock->get_leftBlock()->get_sites();
    bool iinleft = find(lsites.begin(), lsites.end(), i) != lsites.end();
    bool jinleft = find(lsites.begin(), lsites.end(), j) != lsites.end();

    boost::shared_ptr<StackSparseMatrix> dotop, leftop;
    SpinQuantum opq = sq1;
    if (loopblock == b->get_leftBlock()) {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	//double parity3 = getCommuteParity(leftop->get_deltaQuantum()[0], dotop->get_deltaQuantum()[0], opq);
	double parity3 = getCommuteParity(dotop->get_deltaQuantum()[0], leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity*parity3);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2*parity3);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2);
      }
    }
    else {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	//double parity3 = getCommuteParity(leftop->get_deltaQuantum()[0], dotop->get_deltaQuantum()[0], opq);
	double parity3 = getCommuteParity(dotop->get_deltaQuantum()[0], leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity*parity3);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2*parity3);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, factor*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, factor*parity2);
      }
    }


    if (deallocate2) op2->deallocate();
  }
  dmrginp.cctime->stop();
}

void SpinAdapted::stackopxop::ddxcccomp_3indexElement(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, int index, const SpinQuantum& q)
{
  dmrginp.cctime->start();
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  int j = op1->get_orbs(1);
  if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
    return;
  double factor = 2.0; if (i==j) factor = 1.0;
  boost::shared_ptr<StackSparseMatrix> op2;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
    op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i, j).at(0);
  else
    op2 = otherblock->get_op_rep(DES_DESCOMP, -opq, i, j);
  
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
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) {
	op1 = loopblock->get_op_array(DES_DES).get_element(i, j).at(0);
	op2 = otherblock->get_op_array(CRE_CRECOMP).get_element(i, j).at(0);
      }
      else {
	op1 = loopblock->get_op_rep(DES_DES, -opq, i, j);
	op2 = otherblock->get_op_rep(CRE_CRECOMP, opq, i, j);
      }
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      
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
    //op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    //op2->build(*otherblock);
    
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
    SpinQuantum sq1 = op1->get_deltaQuantum(0);
    SpinQuantum sq2 = op2->get_deltaQuantum(0);
    double parity2 =TensorOp::getTransposeFactorDD(i, j, sq1.get_s().getirrep(), sq1.get_symm().getirrep());
    parity2*=TensorOp::getTransposeFactorDD(i, j, sq2.get_s().getirrep(), sq2.get_symm().getirrep());    
    parity2 *= parity;

    

    const std::vector<int>& lsites = loopblock->get_leftBlock()->get_sites();
    bool iinleft = find(lsites.begin(), lsites.end(), i) != lsites.end();
    bool jinleft = find(lsites.begin(), lsites.end(), j) != lsites.end();

    boost::shared_ptr<StackSparseMatrix> dotop, leftop;
    SpinQuantum opq = sq1;
    if (loopblock == b->get_leftBlock()) {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	//double parity3 = getCommuteParity(leftop->get_deltaQuantum()[0], dotop->get_deltaQuantum()[0], opq);
	double parity3 = getCommuteParity(dotop->get_deltaQuantum()[0], leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity*parity3);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity);
      }
    }
    else {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity, true);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity, true);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	//double parity3 = getCommuteParity(leftop->get_deltaQuantum()[0], dotop->get_deltaQuantum()[0], opq);
	double parity3 = getCommuteParity(dotop->get_deltaQuantum()[0], leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity*parity3, true);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_CRE, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, factor*parity, true);
      }
    }


    //if (deallocate2) op2->deallocate();
  }
  dmrginp.cctime->stop();
}


void SpinAdapted::stackopxop::cdxcdcomp_3index(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  dmrginp.cdtime->start();
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
      
      if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j)) {
	dmrginp.cdtime->stop();
	return;
      }
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
	op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(0);
      else
	op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
	
      double factor = 1.0;
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor);

      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    if (i != j) {
      boost::shared_ptr<StackSparseMatrix> op1;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
	op1 = loopblock->get_op_array(DES_CRE).get_element( i, j).at(0);
      else
	op1 = loopblock->get_op_rep(DES_CRE, -opq, i, j);
      double parity = 1.0;
      if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
	op2 = otherblock->get_op_array(DES_CRECOMP).get_element( i, j).at(0);
      else
	op2 = otherblock->get_op_rep(DES_CRECOMP, opq, i, j);

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
    
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j)) {
      dmrginp.cdtime->stop();
      return;
    }
    boost::shared_ptr<StackSparseMatrix> op2;
    if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
      op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(0);
    else
      op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    if (deallocate2) op2->build(*otherblock);
    
    boost::shared_ptr<StackSparseMatrix> dotop, leftop;

    const std::vector<int>& lsites = loopblock->get_leftBlock()->get_sites();
    bool iinleft = find(lsites.begin(), lsites.end(), i) != lsites.end();
    bool jinleft = find(lsites.begin(), lsites.end(), j) != lsites.end();

    if (loopblock == b->get_leftBlock()) {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_DES, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, *op1, b, c, v, hq, 1.0);
	if (i!=j)
	  SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, 1.0);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, Transpose(*dotop), *op1, b, c, v, hq, 1.0);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, Transpose(*dotop), Transpose(*op1), b, c, v, hq, 1.0);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	double parity = getCommuteParity(dotop->get_deltaQuantum()[0], -leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, Transpose(*leftop), *dotop, *op1, b, c, v, hq, parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), Transpose(*leftop), *dotop, Transpose(*op1), b, c, v, hq, parity);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_DES, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, *op1, b, c, v, hq, 1.0);
	if (i!=j)
	  SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, 1.0);
      }
    }
    else {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_DES, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, 1.0);
	if (i!=j)
	  SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, 1.0);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, Transpose(*dotop), *op1, b, c, v, hq, 1.0);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, Transpose(*dotop), Transpose(*op1), b, c, v, hq, 1.0);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	double parity = getCommuteParity(dotop->get_deltaQuantum()[0], -leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, Transpose(*leftop), *dotop, *op1, b, c, v, hq, 1.0*parity);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), Transpose(*leftop), *dotop, Transpose(*op1), b, c, v, hq, 1.0*parity);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_DES, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, 1.0);
	if (i!=j)
	  SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, 1.0);
      }

    }

    
    if (deallocate2) op2->deallocate();
    
  }
  dmrginp.cdtime->stop();

}

void SpinAdapted::stackopxop::cdxcdcomp_3indexElement(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, int index, const SpinQuantum& q)
{
  dmrginp.cdtime->start();
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
      
      if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j)) {
	dmrginp.cdtime->stop();
	return;
      }
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
	op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(0);
      else
	op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
	
      double factor = 1.0;
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor);

      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    if (i != j) {
      boost::shared_ptr<StackSparseMatrix> op1;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
	op1 = loopblock->get_op_array(DES_CRE).get_element( i, j).at(0);
      else
	op1 = loopblock->get_op_rep(DES_CRE, -opq, i, j);
      double parity = 1.0;
      if (dmrginp.spinAdapted() == true && dmrginp.hamiltonian() != BCS)
	parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
	op2 = otherblock->get_op_array(DES_CRECOMP).get_element( i, j).at(0);
      else
	op2 = otherblock->get_op_rep(DES_CRECOMP, opq, i, j);

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
    
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j)) {
      dmrginp.cdtime->stop();
      return;
    }
    boost::shared_ptr<StackSparseMatrix> op2;
    if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
      op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(0);
    else
      op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
    
    boost::shared_ptr<StackSparseMatrix> dotop, leftop;

    const std::vector<int>& lsites = loopblock->get_leftBlock()->get_sites();
    bool iinleft = find(lsites.begin(), lsites.end(), i) != lsites.end();
    bool jinleft = find(lsites.begin(), lsites.end(), j) != lsites.end();

    if (loopblock == b->get_leftBlock()) {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_DES, opq, i, j);
	if (i!=j) {
	  SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index,  1.0);
	}
	else
	  SpinAdapted::operatorfunctions::TensorMultiplysplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index,  1.0);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, Transpose(*dotop), *op1, b, c, v, hq, index,  1.0);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitLeftElement(Transpose(*op2), *leftop, Transpose(*dotop), Transpose(*op1), b, c, v, hq, index,  1.0);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	double parity = getCommuteParity(dotop->get_deltaQuantum()[0], -leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, Transpose(*leftop), *dotop, *op1, b, c, v, hq, index,  parity);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitLeftElement(Transpose(*op2), Transpose(*leftop), *dotop, Transpose(*op1), b, c, v, hq, index,  parity);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_DES, opq, i, j);
	if (i!=j)
	  SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index,  1.0);
	else
	  SpinAdapted::operatorfunctions::TensorMultiplysplitLeftElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index,  1.0);
      }
    }
    else {
      if (iinleft && jinleft) {
	dotop = loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	leftop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_leftBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_leftBlock()->get_op_rep(CRE_DES, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, 1.0, i!=j);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, 1.0);
	//if (i!=j)
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, 1.0);
      }
      else if (iinleft && !jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(j).at(0);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, Transpose(*dotop), *op1, b, c, v, hq, index, 1.0, i!=j);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, Transpose(*dotop), *op1, b, c, v, hq, 1.0);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, Transpose(*dotop), Transpose(*op1), b, c, v, hq, 1.0);
      }
      else if (!iinleft && jinleft) {
	leftop = loopblock->get_leftBlock()->get_op_array(CRE).get_element(j).at(0);
	dotop = loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0);
	double parity = getCommuteParity(dotop->get_deltaQuantum()[0], -leftop->get_deltaQuantum()[0], opq);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, Transpose(*leftop), *dotop, *op1, b, c, v, hq, index, 1.0*parity, i!=j);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, Transpose(*leftop), *dotop, *op1, b, c, v, hq, 1.0*parity);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), Transpose(*leftop), *dotop, Transpose(*op1), b, c, v, hq, 1.0*parity);
      }
      else {
	leftop = loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0);
	dotop = (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) ? loopblock->get_rightBlock()->get_op_array(CRE_DES).get_element(i, j).at(0) : loopblock->get_rightBlock()->get_op_rep(CRE_DES, opq, i, j);
	SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, 1.0, i!=j);
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, *op1, b, c, v, hq, 1.0);
	//if (i!=j)
	//SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, Transpose(*op1), b, c, v, hq, 1.0);
      }

    }

    
    //if (deallocate2) op2->deallocate();
    
  }
  dmrginp.cdtime->stop();

}

void SpinAdapted::stackopxop::cxcddcomp_3index(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
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
      boost::shared_ptr<StackSparseMatrix> op2;// 
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
	op2 = otherblock->get_op_array(CRE_DES_DESCOMP).get_element(i).at(0);
      else
	op2 = otherblock->get_op_rep(CRE_DES_DESCOMP, -opq, i);
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
      boost::shared_ptr<StackSparseMatrix> op1;
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) {
	op1 = loopblock->get_op_array(DES).get_element(i).at(0);
	op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0);
      }
      else {
	op1 = loopblock->get_op_rep(DES, -opq, i);
	op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
      }
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
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
  } 
  else {
    
    boost::shared_ptr<StackSparseMatrix> op2;
    if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
      op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0);
    else
      op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(-op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);

    int parity2 = 1.0;
    if (otherblock == b->get_leftBlock()) parity2 = getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), hq);

    const std::vector<int>& lsites = loopblock->get_leftBlock()->get_sites();
    bool iinleft = find(lsites.begin(), lsites.end(), i) != lsites.end();
    boost::shared_ptr<StackSparseMatrix> dotop = iinleft ? 
      loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0) :
      loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0) ;

    boost::shared_ptr<StackSparseMatrix> leftop = iinleft ? 
      loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0) :
      loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0) ;

    if (loopblock == b->get_leftBlock()) {
      SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(*op2, *leftop, *dotop, Transpose(*op1), b, c, v, hq, scale*parity);
      SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(Transpose(*op2), *leftop, *dotop, *op1, b, c, v, hq, scale*parity2);
    }
    else {
      SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, Transpose(*op1), b, c, v, hq, scale*parity);
      SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, *op1, b, c, v, hq, scale*parity2);
    }

    if (deallocate2) op2->deallocate();

  }
}


void SpinAdapted::stackopxop::cxcddcomp_3indexElement(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, int index, const SpinQuantum& q)
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
      boost::shared_ptr<StackSparseMatrix> op2;// 
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
	op2 = otherblock->get_op_array(CRE_DES_DESCOMP).get_element(i).at(0);
      else
	op2 = otherblock->get_op_rep(CRE_DES_DESCOMP, -opq, i);
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
      boost::shared_ptr<StackSparseMatrix> op1;
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) {
	op1 = loopblock->get_op_array(DES).get_element(i).at(0);
	op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0);
      }
      else {
	op1 = loopblock->get_op_rep(DES, -opq, i);
	op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
      }
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
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
  } 
  else {
    
    boost::shared_ptr<StackSparseMatrix> op2;
    if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) 
      op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0);
    else
      op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
    //bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    //op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    //op2->build(*otherblock);
    
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_rightBlock())
      parity = getCommuteParity(-op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);

    int parity2 = 1.0;
    if (otherblock == b->get_leftBlock()) parity2 = getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), hq);

    const std::vector<int>& lsites = loopblock->get_leftBlock()->get_sites();
    bool iinleft = find(lsites.begin(), lsites.end(), i) != lsites.end();
    boost::shared_ptr<StackSparseMatrix> dotop = iinleft ? 
      loopblock->get_rightBlock()->get_op_array(OVERLAP).get_element(0).at(0) :
      loopblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0) ;

    boost::shared_ptr<StackSparseMatrix> leftop = iinleft ? 
      loopblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0) :
      loopblock->get_leftBlock()->get_op_array(OVERLAP).get_element(0).at(0) ;

    if (loopblock == b->get_leftBlock()) {
      //pout << parity <<"  "<<parity2<<endl;
      SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(*op2, *leftop, *dotop, Transpose(*op1), b, c, v, hq, index, scale*parity);
    }
    else {
      //SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, *op1, b, c, v, hq, index, scale*parity, true);
      SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(*op2, *leftop, *dotop, Transpose(*op1), b, c, v, hq, index, scale*parity, true);
      //SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(Transpose(*op2), *leftop, *dotop, *op1, b, c, v, hq, index, scale*parity, false);
      //SpinAdapted::operatorfunctions::TensorMultiplysplitRight(*op2, *leftop, *dotop, Transpose(*op1), b, c, v, hq, scale*parity);
      //SpinAdapted::operatorfunctions::TensorMultiplysplitRight(Transpose(*op2), *leftop, *dotop, *op1, b, c, v, hq, scale*parity2);
    }

    //if (deallocate2) op2->deallocate();

  }
}


void SpinAdapted::stackopxop::cdxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  dmrginp.cdtime->start();
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
  if (otherblock->has(DES_CRECOMP)) {
    SpinQuantum opq = op1->get_deltaQuantum()[0];
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    {
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false;
      if (deallocate1) {
        op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
        op1->build(*loopblock);
      }
      
      if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
	return;
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
        op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i,j).at(0);
      else
        op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false;
      if (deallocate2) { 
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
      double factor = 1.0;
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor);

      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    if (i != j) {
      boost::shared_ptr<StackSparseMatrix> op1;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
        op1 = loopblock->get_op_array(DES_CRE).get_element(i,j).at(0);
      else      
        op1 = loopblock->get_op_rep(DES_CRE, -opq, i, j);
      double parity = 1.0;
      if (dmrginp.spinAdapted() && dmrginp.hamiltonian() != BCS)
	      parity = getCommuteParity(-getSpinQuantum(i), getSpinQuantum(j), op1->get_deltaQuantum()[0]);
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
        op2 = otherblock->get_op_array(DES_CRECOMP).get_element(i,j).at(0);
      else
        op2 = otherblock->get_op_rep(DES_CRECOMP, opq, i, j);
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false;
     if (deallocate1) { 
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
     }
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
     if (deallocate2) {
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
     }
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, parity);
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  }
  else {
    SpinQuantum opq = op1->get_deltaQuantum()[0];
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false;
   if (deallocate1) { 
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
   }
    int i = op1->get_orbs(0);
    int j = op1->get_orbs(1);
    if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
      return;
    boost::shared_ptr<StackSparseMatrix> op2;
    if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
      op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i, j).at(0);
    else
      op2 = otherblock->get_op_rep(CRE_DESCOMP, -opq, i, j);
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    if (deallocate2) {
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    }
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
  dmrginp.cdtime->stop();
}

void SpinAdapted::stackopxop::ddxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  dmrginp.cctime->start();
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
  
  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  int j = op1->get_orbs(1);
  if (!otherblock->get_op_array(DES_DESCOMP).has_local_index(i,j))
    return;
  double factor = 2.0; if (i==j) factor = 1.0;
  boost::shared_ptr<StackSparseMatrix> op2;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op2 = otherblock->get_op_array(DES_DESCOMP).get_element(i,j).at(0);
  else
    op2 = otherblock->get_op_rep(DES_DESCOMP, -opq, i, j);
  
  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  
  if (otherblock->has(CRE_CRECOMP)) {
    double scale = 1.0;
    double parity = 1.0;
    if (otherblock == b->get_leftBlock())
      parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
    {
      if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      }

      if (deallocate2) {
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      }
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor*parity);
      
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    {
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) {
        op1 = loopblock->get_op_array(DES_DES).get_element(i, j).at(0);
        op2 = otherblock->get_op_array(CRE_CRECOMP).get_element(i, j).at(0);
      } else {
        op1 = loopblock->get_op_rep(DES_DES, -opq, i, j);
        op2 = otherblock->get_op_rep(CRE_CRECOMP, opq, i, j);
      }
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      
      if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
      }
      if (deallocate2) {
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
      }
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, factor*parity);
      
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  }
  else {
    if (deallocate1) {
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
    }
    if (deallocate2) {
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
    }
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
  dmrginp.cctime->stop();
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
      boost::shared_ptr<StackSparseMatrix> op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
        op2 = otherblock->get_op_array(CRE_DES_DESCOMP).get_element(i).at(0);
      else
        op2 = otherblock->get_op_rep(CRE_DES_DESCOMP, -opq, i);

      if (deallocate1) {
        op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
        op1->build(*loopblock);
      }
      
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false;
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
      
      if (otherblock == b->get_leftBlock())
	parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, scale*parity);	    
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
    {
      boost::shared_ptr<StackSparseMatrix> op1, op2;
      if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS) {
        op1 = loopblock->get_op_array(DES).get_element(i).at(0);
        op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0);
      } else {
        op1 = loopblock->get_op_rep(DES, -opq, i);
        op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
      }
      bool deallocate1 = op1->memoryUsed() == 0 ? true : false;
      if (deallocate1) { 
        op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
        op1->build(*loopblock);
      }
      
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false;
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
      
      if (otherblock == b->get_rightBlock()) parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), hq);
      else parity = 1.0;
      
      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, scale*parity);
      if (deallocate2) op2->deallocate();
      if (deallocate1) op1->deallocate();
    }
  } else {
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }
    //StackTransposeview top1 = StackTransposeview(op1);  // DES_i
    
    boost::shared_ptr<StackSparseMatrix> op2;
    if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
      op2 = otherblock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0);
    else
      op2 = otherblock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false;
    if (deallocate2) {
      op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
      op2->build(*otherblock);
    }
    
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


void SpinAdapted::stackopxop::hamandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q, double scale, int proc)
{

  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));  // in get_parity, number part is not used
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();


  boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(OVERLAP, hq);
  bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  if (deallocate2) {
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
  }

  boost::shared_ptr<StackSparseMatrix> op1ham;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op1ham = loopblock->get_op_array(HAM).get_element(0).at(0);
  else
    op1ham = loopblock->get_op_rep(HAM, hq);
  bool deallocate1ham = op1ham->memoryUsed() == 0 ? true : false;
  if (deallocate1ham) {
    op1ham->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1ham->build(*loopblock);
  }

  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1ham, b, c, v, hq, 1.0);	    
  if (deallocate1ham) op1ham->deallocate();


  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  op1 = loopblock->get_op_rep(OVERLAP, hq);
  if (deallocate1) {
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
  }

  boost::shared_ptr<StackSparseMatrix> op2ham;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op2ham = otherblock->get_op_array(HAM).get_element(0).at(0);
  else
    op2ham = otherblock->get_op_rep(HAM, hq);

  bool deallocate2ham = op2ham->memoryUsed() == 0 ? true : false;
  if (deallocate2ham) {
    op2ham->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2ham->build(*otherblock);
  }
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2ham, *op1, b, c, v, hq, 1.0);	    
  if (deallocate2ham) op2ham->deallocate();
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, hq, scale);	    

  if (deallocate1) op1->deallocate();
  if (deallocate2) op2->deallocate();


}


//***************************************************************************************************

/********************************************
Formulas for making diagonal hamiltonian matrix while blocking system and environment blocks
********************************************/


void SpinAdapted::stackopxop::cdxcdcomp_d(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, DiagonalMatrix* e)
{
  int ilock = 0;
  int numthreads = 1;//MAX_THRD;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  SpinQuantum sq = op1->get_deltaQuantum()[0];
  int i = op1->get_orbs(0);
  int j = op1->get_orbs(1);
  
  if (!otherblock->get_op_array(CRE_DESCOMP).has_local_index(i,j))
    return;
  boost::shared_ptr<StackSparseMatrix> op2;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op2 = otherblock->get_op_array(CRE_DESCOMP).get_element(i,j).at(0);
  else
    op2 = otherblock->get_op_rep(CRE_DESCOMP, -sq, i, j);

  double factor = 1.0;
  {
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    if (deallocate1) op1->build(*loopblock);
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    if (deallocate1) op2->build(*otherblock);
    SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), e, factor);
    if (i != j) {
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, Transpose(*op2), Transpose(*op1), b, &(b->get_stateInfo()), e, factor);
    }
    
    if (deallocate2) op2->deallocate();
    if (deallocate1) op1->deallocate();
  }

}

void SpinAdapted::stackopxop::ham_d(const StackSpinBlock* thisBlock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, DiagonalMatrix* e, int proc)
{
  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  op1->allocate(thisBlock->get_braStateInfo(), thisBlock->get_ketStateInfo());
  if (deallocate1) op1->build(*thisBlock);
  SpinAdapted::operatorfunctions::TensorTrace(thisBlock, *op1, b, &(b->get_stateInfo()), e, 1.0);
  if (deallocate1) op1->deallocate();
}


//************************************************************************************

/********************************************
Formulas for making CCdcomp operators while blocking a block with a dot block
********************************************/

void SpinAdapted::stackopxop::cxcdcompElement(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int I, StackSparseMatrix* o, StackMatrix& m, int row, int col, double scale)
{
  int ilock = 0;//omp_get_thread_num();
  int numthreads = 1;
  const StackSpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  if (opvec1[0]->get_orbs(0) >= I) {// opvec1 is CRE
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
	      SpinAdapted::operatorfunctions::TensorProductElement(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, factor*scale); // CD*C
	    }
    }
  } else {
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
	  SpinAdapted::operatorfunctions::TensorProductElement(otherBlock, Transpose(*op2), *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, factor*scale*parity);
	  //op2->set_conjugacy('n');
        } 
	else {
	  double parity = 1.0;
	  if (otherBlock == b->get_rightBlock())
	    parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0)); // doesn't depend on nelec
	  SpinAdapted::operatorfunctions::TensorProductElement(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, factor*scale*parity);
	}
      }
    }
  }
}


void SpinAdapted::stackopxop::dxcccompElement(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int K, StackSparseMatrix* o, StackMatrix& m, int row, int col, double scale)
{ 
  int ilock = 0;//omp_get_thread_num();
  int numthreads = 1;
  //int numthreads = dmrginp.thrds_per_node()[mpigetrank()];
  const StackSpinBlock* loopblock = (otherBlock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind=0; opind<opvec1.size(); opind++) {

    //If we have all the operators we dont have to take transposes, useful for <bra|H|ket> evaluation
    if (loopblock->has(DES) && otherBlock->has(CRE_CRECOMP)) {
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
	
	SpinAdapted::operatorfunctions::TensorProductElement(otherBlock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], m, row, col, parity*factor*scale);
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
	SpinAdapted::operatorfunctions::TensorProductElement(otherBlock, Transpose(*op), Transpose(*op1), b, &(b->get_stateInfo()), o[ilock], m, row, col, parity*factor*scale);
	//op->set_conjugacy('n');op1->set_conjugacy('n');
      }
    }
  }
}



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
    if (loopblock->has(DES) && otherBlock->has(CRE_CRECOMP)) {
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



void SpinAdapted::stackopxop::CreDesonLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;


  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int I = op1->get_orbs(0);
  int J = op1->get_orbs(1);
  StackSpinBlock *leftBlock = cblock->get_leftBlock()->get_leftBlock(), *dotBlock = cblock->get_leftBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock();

  StackSparseMatrix& LEFTOP = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
      *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I,J).at(0) :
      *cblock->get_leftBlock()->get_op_rep(CRE_DES, opq, I, J);
  StackSparseMatrix& leftOp = *op1;

  int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

  const std::vector<int>& ccolinds = c.getActiveCols(luncollectedQPrime); //colinds
  assert(ccolinds.size() == 1); //the c wavefunction has spin = 0
  int rQPrime = ccolinds[0];


  {
    const std::vector<int>& rowinds = leftOp.getActiveRows(lQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;    

    for (int r = 0; r < rowinds.size(); r++) {
      int lQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(leftOp.operator_element(lQ, lQPrime).Nrows() * c.operator_element(luncollectedQPrime, rQPrime).Ncols());
      blocks[r].allocate(data, leftOp.operator_element(lQ, lQPrime).Nrows(), c.operator_element(luncollectedQPrime, rQPrime).Ncols());
      memToDeallocate += blocks[r].Storage();
      //L(l, l') c(l',d',r') ->  b(l, d', r')
      MatrixMultiply (leftOp.operator_element(lQ, lQPrime), 'n', c.operator_element(luncollectedQPrime, rQPrime), 'n', 
		      blocks[r], 1.0, 0.);	      
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackSparseMatrix& rightop = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
      *rightBlock->get_op_array(CRE_DESCOMP).get_element(I,J).at(0) :
      *rightBlock->get_op_rep(CRE_DESCOMP, opq, I, J);
    operatorfunctions::multiplyDotRight(LEFTOP, leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);

    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }

  if (I != J)
  {
    const std::vector<int>& rowinds = leftOp.getActiveCols(lQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;
    
    for (int r = 0; r < rowinds.size(); r++) {
      int rQPrime = ccolinds[0];
      int lQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(lketS->quantaStates[lQ] * rketS->quantaStates[rQPrime]);
      blocks[r].allocate(data, lketS->quantaStates[lQ], rketS->quantaStates[rQPrime]);
      memToDeallocate += blocks[r].Storage();
      //L(l', l) c(l',d',r') ->  b(l, d', r')
      MatrixMultiply (leftOp.operator_element(lQPrime, lQ), 't', c.operator_element(luncollectedQPrime, rQPrime), 'n', 
		      blocks[r], 1.0, 0.);	      
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackTransposeview rightop(
        !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ? 
        *rightBlock->get_op_array(CRE_DESCOMP).get_element(I,J).at(0): 
        *rightBlock->get_op_rep(CRE_DESCOMP, opq, I, J));
    operatorfunctions::multiplyDotRight(Transpose(LEFTOP), leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);

    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }
  


}



void SpinAdapted::stackopxop::CreCreonLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;


  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int I = op1->get_orbs(0);
  int J = op1->get_orbs(1);
  double factor = I==J ? 1.0 : 2.0;
  StackSpinBlock *leftBlock = cblock->get_leftBlock()->get_leftBlock(), *dotBlock = cblock->get_leftBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock();

  StackSparseMatrix& LEFTOP = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
    *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I,J).at(0):
    *cblock->get_leftBlock()->get_op_rep(CRE_CRE, opq, I, J);
  StackSparseMatrix& leftOp = *op1;

  int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

  const std::vector<int>& ccolinds = c.getActiveCols(luncollectedQPrime); //colinds
  assert(ccolinds.size() == 1); //the c wavefunction has spin = 0
  int rQPrime = ccolinds[0];


  {
    const std::vector<int>& rowinds = leftOp.getActiveRows(lQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;    

    for (int r = 0; r < rowinds.size(); r++) {
      int lQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(leftOp.operator_element(lQ, lQPrime).Nrows() * c.operator_element(luncollectedQPrime, rQPrime).Ncols());
      blocks[r].allocate(data, leftOp.operator_element(lQ, lQPrime).Nrows(), c.operator_element(luncollectedQPrime, rQPrime).Ncols());
      memToDeallocate += blocks[r].Storage();
      //L(l, l') c(l',d',r') ->  b(l, d', r')
      MatrixMultiply (leftOp.operator_element(lQ, lQPrime), 'n', c.operator_element(luncollectedQPrime, rQPrime), 'n', 
		      blocks[r], 1.0, 0.);	      
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackSparseMatrix& rightop = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
      *rightBlock->get_op_array(DES_DESCOMP).get_element(I,J).at(0):
      *rightBlock->get_op_rep(DES_DESCOMP, -opq, I, J);
    operatorfunctions::multiplyDotRight(LEFTOP, leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, factor);

    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }

  {
    const std::vector<int>& rowinds = leftOp.getActiveCols(lQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;
    
    for (int r = 0; r < rowinds.size(); r++) {
      int rQPrime = ccolinds[0];
      int lQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(lketS->quantaStates[lQ] * rketS->quantaStates[rQPrime]);
      blocks[r].allocate(data, lketS->quantaStates[lQ], rketS->quantaStates[rQPrime]);
      memToDeallocate += blocks[r].Storage();
      //L(l', l) c(l',d',r') ->  b(l, d', r')
      MatrixMultiply (leftOp.operator_element(lQPrime, lQ), 't', c.operator_element(luncollectedQPrime, rQPrime), 'n', 
		      blocks[r], 1.0, 0.);	      
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackTransposeview rightop(!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *rightBlock->get_op_array(DES_DESCOMP).get_element(I,J).at(0):
        *rightBlock->get_op_rep(DES_DESCOMP, -opq, I, J));
    operatorfunctions::multiplyDotRight(Transpose(LEFTOP), leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, factor);

    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }
  


}
  



void SpinAdapted::stackopxop::CreonLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;


  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  StackSpinBlock *leftBlock = cblock->get_leftBlock()->get_leftBlock(), *dotBlock = cblock->get_leftBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock();
  StackSparseMatrix& LEFTOP = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
    *cblock->get_leftBlock()->get_op_array(CRE).get_element(i).at(0):
    *cblock->get_leftBlock()->get_op_rep(CRE, opq, i);
  StackSparseMatrix& leftOp = *op1;
  const std::vector<int>& dotindices = dotBlock->get_sites();

  int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

  const std::vector<int>& ccolinds = c.getActiveCols(luncollectedQPrime); //colinds
  assert(ccolinds.size() == 1); //the c wavefunction has spin = 0
  int rQPrime = ccolinds[0];


  {
    const std::vector<int>& rowinds = leftOp.getActiveRows(lQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;    

    for (int r = 0; r < rowinds.size(); r++) {
      int lQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(leftOp.operator_element(lQ, lQPrime).Nrows() * c.operator_element(luncollectedQPrime, rQPrime).Ncols());
      blocks[r].allocate(data, leftOp.operator_element(lQ, lQPrime).Nrows(), c.operator_element(luncollectedQPrime, rQPrime).Ncols());
      memToDeallocate += blocks[r].Storage();
      //L(l, l') c(l',d',r') ->  b(l, d', r')
      MatrixMultiply (leftOp.operator_element(lQ, lQPrime), 'n', c.operator_element(luncollectedQPrime, rQPrime), 'n', 
		      blocks[r], 1.0, 0.);	      
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackTransposeview rightop(!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS?
        rightBlock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0):
        rightBlock->get_op_rep(CRE_CRE_DESCOMP, opq, i));
    operatorfunctions::multiplyDotRight(LEFTOP, leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);
    

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      int I = dx;
      int J = i;
      if (rightBlock->get_op_array(CRE_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackTransposeview rightop1(rightBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(0));
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I, J).at(0);
	double parity = getCommuteParity(dotop.get_deltaQuantum()[0], -leftOp.get_deltaQuantum()[0], LEFTOP1.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(Transpose(LEFTOP1), Transpose(leftOp), dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, parity);

	StackTransposeview rightop2(*rightBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(1));
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I, J).at(1);
	double parity2 = getCommuteParity(dotop.get_deltaQuantum()[0], -leftOp.get_deltaQuantum()[0], LEFTOP2.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(Transpose(LEFTOP2), Transpose(leftOp), dotop, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, parity2);

      }

      //c c  ddcomp
      if (rightBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackSparseMatrix& rightop1 = *rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	double parity1 = getCommuteParity(dotop.get_deltaQuantum()[0], leftOp.get_deltaQuantum()[0], LEFTOP1.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(LEFTOP1, leftOp, dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, 2.0*parity1);

	StackSparseMatrix& rightop2 = *rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	double parity2 = getCommuteParity(dotop.get_deltaQuantum()[0], leftOp.get_deltaQuantum()[0], LEFTOP2.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(LEFTOP2, leftOp, dotop, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, 2.0*parity2);

      }

    }


    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }


  //CRE IS TRANSPOSED
  {
    const std::vector<int>& rowinds = leftOp.getActiveCols(lQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;
    
    for (int r = 0; r < rowinds.size(); r++) {
      int rQPrime = ccolinds[0];
      int lQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(lketS->quantaStates[lQ] * rketS->quantaStates[rQPrime]);
      blocks[r].allocate(data, lketS->quantaStates[lQ], rketS->quantaStates[rQPrime]);
      memToDeallocate += blocks[r].Storage();
      //L(l', l) c(l',d',r') ->  b(l, d', r')
      MatrixMultiply (leftOp.operator_element(lQPrime, lQ), 't', c.operator_element(luncollectedQPrime, rQPrime), 'n', 
		      blocks[r], 1.0, 0.);	      
    }

    //c I  ccdcomp
    {
      StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
      StackSparseMatrix& rightop = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *rightBlock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0) :
        *rightBlock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
      operatorfunctions::multiplyDotRight(Transpose(LEFTOP), leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);
    }

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      int I = dx;
      int J = i;

      //c d  cdcomp
      if (rightBlock->get_op_array(CRE_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackSparseMatrix& rightop1 = *rightBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I, J).at(0);
	double parity = getCommuteParity(dotop.get_deltaQuantum()[0], -leftOp.get_deltaQuantum()[0], LEFTOP1.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(LEFTOP1, Transpose(leftOp), dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, parity);

	StackSparseMatrix& rightop2 = *rightBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I, J).at(1);
	double parity2 = getCommuteParity(dotop.get_deltaQuantum()[0], -leftOp.get_deltaQuantum()[0], LEFTOP2.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(LEFTOP2, Transpose(leftOp), dotop, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, parity2);

      }

      //c c  ddcomp
      if (rightBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackTransposeview rightop1(*rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0));
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	double parity1 = getCommuteParity(dotop.get_deltaQuantum()[0], leftOp.get_deltaQuantum()[0], LEFTOP1.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(Transpose(LEFTOP1), leftOp, dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, 2.0*parity1);

	StackTransposeview rightop2(*rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1));
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	double parity2 = getCommuteParity(dotop.get_deltaQuantum()[0], leftOp.get_deltaQuantum()[0], LEFTOP2.get_deltaQuantum()[0]);
	operatorfunctions::multiplyDotRight(Transpose(LEFTOP2), leftOp, dotop, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, 2.0*parity2);

      }

    }
    
    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }
  


}
  
  

void SpinAdapted::stackopxop::OverlaponLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;


  StackSpinBlock *leftBlock = cblock->get_leftBlock()->get_leftBlock(), *dotBlock = cblock->get_leftBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock();
  StackSparseMatrix& leftOp = *leftBlock->get_op_array(OVERLAP).get_element(0).at(0);
  const std::vector<int>& dotindices = dotBlock->get_sites();

  int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

  const std::vector<int>& ccolinds = c.getActiveCols(luncollectedQPrime); //colinds
  assert(ccolinds.size() == 1); //the c wavefunction has spin = 0
  int rQPrime = ccolinds[0];

  {
    vector<StackMatrix> blocks(1, c.operator_element(luncollectedQPrime, rQPrime));

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
      SpinQuantum opq = dotop.get_deltaQuantum()[0];
      StackTransposeview rightop(!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
          rightBlock->get_op_array(CRE_CRE_DESCOMP).get_element(dx).at(0) :
          rightBlock->get_op_rep(CRE_CRE_DESCOMP, opq, dx));
      StackSparseMatrix& LEFTOP1 = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *cblock->get_leftBlock()->get_op_array(CRE).get_element(dx).at(0) :
        *cblock->get_leftBlock()->get_op_rep(CRE, opq, dx);
      operatorfunctions::multiplyDotRight(LEFTOP1, leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);
    }

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
    for (int dotx1=0; dotx1<=dotx; dotx1++) {
      int dx = dotindices[dotx];
      int dx1 = dotindices[dotx1];
      int I = dx > dx1 ? dx : dx1;
      int J = dx > dx1 ? dx1 : dx;
      if (rightBlock->get_op_array(CRE_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE_DES).get_element(I, J).at(0);
	StackSparseMatrix& leftop = *leftBlock->get_op_array(OVERLAP).get_element(0).at(0);
	StackSparseMatrix& rightop1 = *rightBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I, J).at(0);
	operatorfunctions::multiplyDotRight(LEFTOP1, leftOp, dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);

	StackSparseMatrix& rightop2 = *rightBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& dotop2 = *dotBlock->get_op_array(CRE_DES).get_element(I, J).at(1);
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_DES).get_element(I, J).at(1);
	operatorfunctions::multiplyDotRight(LEFTOP2, leftOp, dotop2, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);
	
      }

      if (rightBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	double factor = I==J ? 1.0 : 2.0;
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(0);
	StackSparseMatrix& leftop = *leftBlock->get_op_array(OVERLAP).get_element(0).at(0);

	StackSparseMatrix& rightop1 = *rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	operatorfunctions::multiplyDotRight(LEFTOP1, leftOp, dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, factor);

	StackSparseMatrix& rightop2 = *rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& dotop2 = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(1);
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	operatorfunctions::multiplyDotRight(LEFTOP2, leftOp, dotop2, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, factor);
	
      }

    }
    }
  }


  {
    vector<StackMatrix> blocks(1, c.operator_element(luncollectedQPrime, rQPrime));

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
      SpinQuantum opq = dotop.get_deltaQuantum()[0];    
      StackSparseMatrix& rightop = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *rightBlock->get_op_array(CRE_CRE_DESCOMP).get_element(dx).at(0) :
        *rightBlock->get_op_rep(CRE_CRE_DESCOMP, opq, dx);
      StackSparseMatrix& LEFTOP1 = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *cblock->get_leftBlock()->get_op_array(CRE).get_element(dx).at(0) :
        *cblock->get_leftBlock()->get_op_rep(CRE, opq, dx);
      operatorfunctions::multiplyDotRight(Transpose(LEFTOP1), leftOp, dotop, rightop, blocks, v, cblock, luncollectedQPrime, rQPrime, 1.0);
    }

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
    for (int dotx1=0; dotx1<=dotx; dotx1++) {
      int dx = dotindices[dotx];
      int dx1 = dotindices[dotx1];
      int I = dx > dx1 ? dx : dx1;
      int J = dx > dx1 ? dx1 : dx;

      if (rightBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	double factor = I==J ? 1.0 : 2.0;
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(0);
	StackSparseMatrix& leftop = *leftBlock->get_op_array(OVERLAP).get_element(0).at(0);
	StackTransposeview rightop1(*rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0));
	StackSparseMatrix& LEFTOP1 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	operatorfunctions::multiplyDotRight(Transpose(LEFTOP1), leftOp, dotop, rightop1, blocks, v, cblock, luncollectedQPrime, rQPrime, factor);

	StackTransposeview rightop2(*rightBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1));
	StackSparseMatrix& dotop2 = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(1);
	StackSparseMatrix& LEFTOP2 = *cblock->get_leftBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	operatorfunctions::multiplyDotRight(Transpose(LEFTOP2), leftOp, dotop2, rightop2, blocks, v, cblock, luncollectedQPrime, rQPrime, factor);
	
      }

    }
    }

  }


}
  

void SpinAdapted::stackopxop::CreCreonRight(boost::shared_ptr<StackSparseMatrix> op1, int runcollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;


  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int I = op1->get_orbs(0);
  int J = op1->get_orbs(1);
  double factor = I==J ? 1.0 : 2.0;
  StackSpinBlock *leftBlock = cblock->get_leftBlock(), *dotBlock = cblock->get_rightBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock()->get_leftBlock();

  StackSparseMatrix& RIGHTOP = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
    *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I,J).at(0):
    *cblock->get_rightBlock()->get_op_rep(CRE_CRE, opq, I, J);
  StackSparseMatrix& rightop = *op1;

  int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

  const std::vector<int>& crowinds = c.getActiveRows(runcollectedQPrime); //colinds
  assert(crowinds.size() == 1); //the c wavefunction has spin = 0
  int lQPrime = crowinds[0];


  {
    const std::vector<int>& rowinds = rightop.getActiveRows(rQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;    

    for (int r = 0; r < rowinds.size(); r++) {
      int rQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(c.operator_element(lQPrime, runcollectedQPrime).Nrows()* rightop.operator_element(rQ, rQPrime).Nrows());
      blocks[r].allocate(data, c.operator_element(lQPrime, runcollectedQPrime).Nrows(), rightop.operator_element(rQ, rQPrime).Nrows());
      memToDeallocate += blocks[r].Storage();
      //R(r, r') c(l',r',d') ->  b(l', r, d')
      MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightop.operator_element(rQ, rQPrime), 't', blocks[r], 1.0, 0.0);

    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackSparseMatrix& leftOp = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
      *leftBlock->get_op_array(DES_DESCOMP).get_element(I,J).at(0):
      *leftBlock->get_op_rep(DES_DESCOMP, -opq, I, J);
    operatorfunctions::multiplyDotLeft(RIGHTOP, rightop, dotop, leftOp, blocks, v, cblock, lQPrime, runcollectedQPrime, factor);


    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }


  //CRE IS TRANSPOSED
  {
    const std::vector<int>& colinds = rightop.getActiveCols(rQPrime);
    vector<StackMatrix> blocks(colinds.size(), StackMatrix());
    long memToDeallocate = 0;
    
    for (int r = 0; r < colinds.size(); r++) {
      int rQ = colinds[r];
      double* data = Stackmem[omprank].allocate(c.operator_element(lQPrime, runcollectedQPrime).Nrows()* rightop.operator_element(rQPrime, rQ).Ncols());
      blocks[r].allocate(data, c.operator_element(lQPrime, runcollectedQPrime).Nrows(), rightop.operator_element(rQPrime, rQ).Ncols());
      memToDeallocate += blocks[r].Storage();
      //R(r', r) c(l',r',d') ->  b(l', r, d')
      MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightop.operator_element(rQPrime, rQ), 'n', blocks[r], 1.0, 0.0);
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackTransposeview leftOp(!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0) :
        *leftBlock->get_op_rep(DES_DESCOMP, -opq, I, J));
    operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP), rightop, dotop, leftOp, blocks, v, cblock, lQPrime, runcollectedQPrime, factor);


    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }



}
  
void SpinAdapted::stackopxop::CreDesonRight(boost::shared_ptr<StackSparseMatrix> op1, int runcollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;


  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int I = op1->get_orbs(0);
  int J = op1->get_orbs(1);
  StackSpinBlock *leftBlock = cblock->get_leftBlock(), *dotBlock = cblock->get_rightBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock()->get_leftBlock();

  StackSparseMatrix& RIGHTOP = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
    *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I,J).at(0):
    *cblock->get_rightBlock()->get_op_rep(CRE_DES, opq, I, J);
  StackSparseMatrix& rightop = *op1;

  int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

  const std::vector<int>& crowinds = c.getActiveRows(runcollectedQPrime); //colinds
  assert(crowinds.size() == 1); //the c wavefunction has spin = 0
  int lQPrime = crowinds[0];


  {
    const std::vector<int>& rowinds = rightop.getActiveRows(rQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;    

    for (int r = 0; r < rowinds.size(); r++) {
      int rQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(c.operator_element(lQPrime, runcollectedQPrime).Nrows()* rightop.operator_element(rQ, rQPrime).Nrows());
      blocks[r].allocate(data, c.operator_element(lQPrime, runcollectedQPrime).Nrows(), rightop.operator_element(rQ, rQPrime).Nrows());
      memToDeallocate += blocks[r].Storage();
      //R(r, r') c(l',r',d') ->  b(l', r, d')
      MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightop.operator_element(rQ, rQPrime), 't', blocks[r], 1.0, 0.0);

    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackSparseMatrix& leftOp = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
      *leftBlock->get_op_array(CRE_DESCOMP).get_element(I,J).at(0):
      *leftBlock->get_op_rep(CRE_DESCOMP, opq, I, J);
    operatorfunctions::multiplyDotLeft(RIGHTOP, rightop, dotop, leftOp, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);


    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }


  //CRE IS TRANSPOSED
  if (I != J)
  {
    const std::vector<int>& colinds = rightop.getActiveCols(rQPrime);
    vector<StackMatrix> blocks(colinds.size(), StackMatrix());
    long memToDeallocate = 0;
    
    for (int r = 0; r < colinds.size(); r++) {
      int rQ = colinds[r];
      double* data = Stackmem[omprank].allocate(c.operator_element(lQPrime, runcollectedQPrime).Nrows()* rightop.operator_element(rQPrime, rQ).Ncols());
      blocks[r].allocate(data, c.operator_element(lQPrime, runcollectedQPrime).Nrows(), rightop.operator_element(rQPrime, rQ).Ncols());
      memToDeallocate += blocks[r].Storage();
      //R(r', r) c(l',r',d') ->  b(l', r, d')
      MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightop.operator_element(rQPrime, rQ), 'n', blocks[r], 1.0, 0.0);
    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackTransposeview leftOp(!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *leftBlock->get_op_array(CRE_DESCOMP).get_element(I,J).at(0) :
        *leftBlock->get_op_rep(CRE_DESCOMP, opq, I, J));
    operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP), rightop, dotop, leftOp, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);


    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }



}

void SpinAdapted::stackopxop::CreonRight(boost::shared_ptr<StackSparseMatrix> op1, int runcollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;


  SpinQuantum opq = op1->get_deltaQuantum()[0];    
  int i = op1->get_orbs(0);
  StackSpinBlock *leftBlock = cblock->get_leftBlock(), *dotBlock = cblock->get_rightBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock()->get_leftBlock();
  StackSparseMatrix& RIGHTOP = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
    *cblock->get_rightBlock()->get_op_array(CRE).get_element(i).at(0) :
    *cblock->get_rightBlock()->get_op_rep(CRE, opq, i);
  StackSparseMatrix& rightop = *op1;
  const std::vector<int>& dotindices = dotBlock->get_sites();

  int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

  const std::vector<int>& crowinds = c.getActiveRows(runcollectedQPrime); //colinds
  assert(crowinds.size() == 1); //the c wavefunction has spin = 0
  int lQPrime = crowinds[0];


  {
    const std::vector<int>& rowinds = rightop.getActiveRows(rQPrime);
    vector<StackMatrix> blocks(rowinds.size(), StackMatrix());
    long memToDeallocate = 0;    

    for (int r = 0; r < rowinds.size(); r++) {
      int rQ = rowinds[r];
      double* data = Stackmem[omprank].allocate(c.operator_element(lQPrime, runcollectedQPrime).Nrows()* rightop.operator_element(rQ, rQPrime).Nrows());
      blocks[r].allocate(data, c.operator_element(lQPrime, runcollectedQPrime).Nrows(), rightop.operator_element(rQ, rQPrime).Nrows());
      memToDeallocate += blocks[r].Storage();
      //R(r, r') c(l',r',d') ->  b(l', r, d')
      MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightop.operator_element(rQ, rQPrime), 't', blocks[r], 1.0, 0.0);

    }

    StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
    StackTransposeview leftOp(!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        leftBlock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0) :
        leftBlock->get_op_rep(CRE_CRE_DESCOMP, opq, i));
    SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));  // in get_parity, number part is not used
    int parity = getCommuteParity(RIGHTOP.get_deltaQuantum(0), -(leftOp.get_deltaQuantum(0)), hq);
    operatorfunctions::multiplyDotLeft(RIGHTOP, rightop, dotop, leftOp, blocks, v, cblock, lQPrime, runcollectedQPrime, parity);

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      int I = i;
      int J = dx;
      //cdcomp c d  left right dot
      if (leftBlock->get_op_array(CRE_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackSparseMatrix& leftop1 = *leftBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(RIGHTOP1, rightop, Transpose(dotop), leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);

	StackSparseMatrix& leftop2 = *leftBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(RIGHTOP2, rightop, Transpose(dotop), leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);

      }
    }

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      int I = i;
      int J = dx;
      //cdcomp c d  left right dot
      if (leftBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackSparseMatrix& leftop1 = *leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(RIGHTOP1, rightop, dotop, leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, 2.0);

	StackSparseMatrix& leftop2 = *leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(RIGHTOP2, rightop, dotop, leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, 2.0);

      }
    }

    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }


  //CRE IS TRANSPOSED
  {
    const std::vector<int>& colinds = rightop.getActiveCols(rQPrime);
    vector<StackMatrix> blocks(colinds.size(), StackMatrix());
    long memToDeallocate = 0;
    
    for (int r = 0; r < colinds.size(); r++) {
      int rQ = colinds[r];
      double* data = Stackmem[omprank].allocate(c.operator_element(lQPrime, runcollectedQPrime).Nrows()* rightop.operator_element(rQPrime, rQ).Ncols());
      blocks[r].allocate(data, c.operator_element(lQPrime, runcollectedQPrime).Nrows(), rightop.operator_element(rQPrime, rQ).Ncols());
      memToDeallocate += blocks[r].Storage();
      //R(r', r) c(l',r',d') ->  b(l', r, d')
      MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightop.operator_element(rQPrime, rQ), 'n', blocks[r], 1.0, 0.0);
    }

    //c I  ccdcomp
    {
      StackSparseMatrix& dotop = *dotBlock->get_op_array(OVERLAP).get_element(0).at(0);
      StackSparseMatrix& leftOp = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *leftBlock->get_op_array(CRE_CRE_DESCOMP).get_element(i).at(0):
        *leftBlock->get_op_rep(CRE_CRE_DESCOMP, opq, i);
      SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));  // in get_parity, number part is not used
      int parity = getCommuteParity(RIGHTOP.get_deltaQuantum(0), -(leftOp.get_deltaQuantum(0)), hq);
      operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP), rightop, dotop, leftOp, blocks, v, cblock, lQPrime, runcollectedQPrime, parity);
    }

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      int I = i;
      int J = dx;

      //cdcomp c d  left right dot
      if (leftBlock->get_op_array(CRE_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackTransposeview leftop1(*leftBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(0));
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP1), rightop, Transpose(dotop), leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);

	StackTransposeview leftop2(*leftBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(1));
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP2), rightop, Transpose(dotop), leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);

      }

      //cdcomp c d  left right dot
      if (leftBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
	StackTransposeview leftop1(*leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0));
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP1), rightop, dotop, leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, 2.0);

	StackTransposeview leftop2(*leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1));
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP2), rightop, dotop, leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, 2.0);

      }

    }

    if (blocks.size() != 0)
      Stackmem[omprank].deallocate(blocks[0].Store(), memToDeallocate);
  }



}


void SpinAdapted::stackopxop::OverlaponRight(boost::shared_ptr<StackSparseMatrix> op1, int runcollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;


  StackSpinBlock *leftBlock = cblock->get_leftBlock(), *dotBlock = cblock->get_rightBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock()->get_leftBlock();
  StackSparseMatrix& rightop = *rightBlock->get_op_array(OVERLAP).get_element(0).at(0);
  const std::vector<int>& dotindices = dotBlock->get_sites();

  int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

  const std::vector<int>& crowinds = c.getActiveRows(runcollectedQPrime); //colinds
  assert(crowinds.size() == 1); //the c wavefunction has spin = 0
  int lQPrime = crowinds[0];


  {
    vector<StackMatrix> blocks(1, c.operator_element(lQPrime, runcollectedQPrime));

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
      SpinQuantum opq = dotop.get_deltaQuantum()[0];    
      StackTransposeview leftop(*leftBlock->get_op_array(CRE_CRE_DESCOMP).get_element(dx).at(0));
      StackSparseMatrix& RIGHTOP1 = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *cblock->get_rightBlock()->get_op_array(CRE).get_element(dx).at(0) :
        *cblock->get_rightBlock()->get_op_rep(CRE, opq, dx);
      operatorfunctions::multiplyDotLeft(RIGHTOP1, rightop, dotop, leftop, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);
    }

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
    for (int dotx1=0; dotx1<=dotx; dotx1++) {
      int dx = dotindices[dotx];
      int dx1 = dotindices[dotx1];
      int I = dx > dx1 ? dx : dx1;
      int J = dx > dx1 ? dx1 : dx;
      if (leftBlock->get_op_array(CRE_DESCOMP).has_local_index(I, J)) {
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE_DES).get_element(I, J).at(0);
	StackSparseMatrix& leftop1 = *leftBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(RIGHTOP1, rightop, dotop, leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);

	StackSparseMatrix& leftop2 = *leftBlock->get_op_array(CRE_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& dotop2 = *dotBlock->get_op_array(CRE_DES).get_element(I, J).at(1);
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_DES).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(RIGHTOP2, rightop, dotop2, leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);
	
      }

      if (leftBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	double factor = I==J ? 1.0 : 2.0;
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(0);
	StackSparseMatrix& leftop1 = *leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0);
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(RIGHTOP1, rightop, dotop, leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, factor);

	StackSparseMatrix& leftop2 = *leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1);
	StackSparseMatrix& dotop2 = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(1);
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(RIGHTOP2, rightop, dotop2, leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, factor);
	
      }
    }
    }

  }


  //CRE IS TRANSPOSED
  {
    vector<StackMatrix> blocks(1, c.operator_element(lQPrime, runcollectedQPrime));

    for (int dotx=0; dotx<dotindices.size(); dotx++) {
      int dx = dotindices[dotx];
      StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE).get_element(dx).at(0);
      SpinQuantum opq = dotop.get_deltaQuantum()[0];
      StackSparseMatrix& leftop = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *leftBlock->get_op_array(CRE_CRE_DESCOMP).get_element(dx).at(0) :
        *leftBlock->get_op_rep(CRE_CRE_DESCOMP, opq, dx);
      StackSparseMatrix& RIGHTOP1 = !dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS ?
        *cblock->get_rightBlock()->get_op_array(CRE).get_element(dx).at(0) :
        *cblock->get_rightBlock()->get_op_rep(CRE, opq, dx);
      operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP1), rightop, dotop, leftop, blocks, v, cblock, lQPrime, runcollectedQPrime, 1.0);
    }


    for (int dotx=0; dotx<dotindices.size(); dotx++) {
    for (int dotx1=0; dotx1<=dotx; dotx1++) {
      int dx = dotindices[dotx];
      int dx1 = dotindices[dotx1];
      int I = dx > dx1 ? dx : dx1;
      int J = dx > dx1 ? dx1 : dx;

      if (leftBlock->get_op_array(DES_DESCOMP).has_local_index(I, J)) {
	double factor = I==J ? 1.0 : 2.0;
	StackSparseMatrix& dotop = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(0);
	StackTransposeview leftop1(*leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(0));
	StackSparseMatrix& RIGHTOP1 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(0);
	operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP1), rightop, dotop, leftop1, blocks, v, cblock, lQPrime, runcollectedQPrime, factor);

	StackTransposeview leftop2(*leftBlock->get_op_array(DES_DESCOMP).get_element(I, J).at(1));
	StackSparseMatrix& dotop2 = *dotBlock->get_op_array(CRE_CRE).get_element(I, J).at(1);
	StackSparseMatrix& RIGHTOP2 = *cblock->get_rightBlock()->get_op_array(CRE_CRE).get_element(I, J).at(1);
	operatorfunctions::multiplyDotLeft(Transpose(RIGHTOP2), rightop, dotop2, leftop2, blocks, v, cblock, lQPrime, runcollectedQPrime, factor);
	
      }

    }
    }
  }



}
  

//**********************************************************************************************************
  
void SpinAdapted::stackopxop::cdd_cxddcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_DES_DESCOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }

    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_DES_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
      
	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity);	    
      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();

  }
}

void SpinAdapted::stackopxop::cdd_dxcdcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_CRE_DESCOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }

    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

    
      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity);	    
      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();
  }
}

void SpinAdapted::stackopxop::ccd_dxcccomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1); // DES_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_CRECOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }

    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_CRECOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
      
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity);	    
      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();

  }
}

void SpinAdapted::stackopxop::ccd_cxcdcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();

  for (int opind1=0; opind1<opvec1.size(); opind1++) {
    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1); // CRE_i
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_DESCOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }

    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
    
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
      double scale = 1.0;
      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
      
      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity);	    
      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();
  }
}

//void SpinAdapted::stackopxop::cdd_cxddcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const StackSpinQuantum& q)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//    
//  for (int opind1=0; opind1<opvec1.size(); opind1++) {
//    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1);
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CDD_DES_DESCOMP).has_local_index(i))
//      return;
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_DES_DESCOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//
//	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
//	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
//	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
//	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
//	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
//
//      double parity = 1.0;
//      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
//
//      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
//    }
//  }
//}
//
//void SpinAdapted::stackopxop::cdd_dxcdcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//    
//  for (int opind1=0; opind1<opvec1.size(); opind1++) {
//    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1);
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CDD_CRE_DESCOMP).has_local_index(i))
//      return;
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_CRE_DESCOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//
//	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
//	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
//	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
//	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
//	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
//
//
//      double parity = 1.0;
//      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
//      //if (loopblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
//
//      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
//
//    }
//  }
//}
//
////**********************************************************************************************************

//void SpinAdapted::stackopxop::ccd_dxcccomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//    
//  for (int opind1=0; opind1<opvec1.size(); opind1++) {
//    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1);
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CCD_CRE_CRECOMP).has_local_index(i))
//      return;
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_CRECOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//
//	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CCD
//	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
//	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
//	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
//	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
//
//      double parity = 1.0;
//      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
//
//      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
//    }
//  }
//}
//
//void SpinAdapted::stackopxop::ccd_cxcdcomp(const StackSpinBlock* otherblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//    
//  for (int opind1=0; opind1<opvec1.size(); opind1++) {
//    boost::shared_ptr<StackSparseMatrix> op1 = opvec1.at(opind1);
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CCD_CRE_DESCOMP).has_local_index(i))
//      return;
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_DESCOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//
//	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
//	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
//	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
//	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
//	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
//
//
//      double parity = 1.0;
//      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
//      //if (loopblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
//
//      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
//
//    }
//  }
//}
//

//void SpinAdapted::stackopxop::cdd_cxddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CDD_DES_DESCOMP).has_local_index(i))
//      return;
//
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_DES_DESCOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//      
//	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
//	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
//	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
//	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
//	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
//
//      double parity = 1.0;
//      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
//      
//      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity);	    
//    }
//
//}
//
//void SpinAdapted::stackopxop::cdd_dxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CDD_CRE_DESCOMP).has_local_index(i))
//      return;
//
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_CRE_DESCOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//
//	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
//	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
//	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
//	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
//	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());
//
//    
//      double parity = 1.0;
//      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
//      
//      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], factor*parity);	    
//    }
//}
//
//void SpinAdapted::stackopxop::ccd_dxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CCD_CRE_CRECOMP).has_local_index(i))
//      return;
//
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_CRECOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//      
//      double scale = 1.0;
//      double parity = 1.0;
//      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
//      
//      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity);	    
//    }
//
//}
//
//void SpinAdapted::stackopxop::ccd_cxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o)
//{
//  int ilock = 0;
//  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
//
//    int i = op1->get_orbs(0);
//    if (!otherblock->get_op_array(CCD_CRE_DESCOMP).has_local_index(i))
//      return;
//
//    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_DESCOMP).get_element(i); // P_{ki}
//    for (int opind2=0; opind2<opvec2.size(); opind2++) {
//      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
//    
//      double scale = 1.0;
//      double parity = 1.0;
//      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), o->get_deltaQuantum(0));
//      
//      SpinAdapted::operatorfunctions::TensorProduct(otherblock, *op2, *op1, b, &(b->get_stateInfo()), o[ilock], scale*parity);	    
//    }
//}
//
void SpinAdapted::stackopxop::cdd_cxddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_DES_DESCOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }
    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_DES_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);

      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();
}

void SpinAdapted::stackopxop::cdd_dxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CDD_CRE_DESCOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }
    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CDD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);

      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }
	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());


      double parity = 1.0;
      if (otherblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
      //if (loopblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);

      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();
}

void SpinAdapted::stackopxop::ccd_dxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_CRECOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }
    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_CRECOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CCD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());

      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
      if (deallocate2) op2->deallocate();
    }
    if (deallocate1) op1->deallocate();
}

void SpinAdapted::stackopxop::ccd_cxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{
  int ilock = 0;
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();
    
    int i = op1->get_orbs(0);
    if (!otherblock->get_op_array(CCD_CRE_DESCOMP).has_local_index(i))
      return;
    bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
      op1->build(*loopblock);
    }
    const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec2 = otherblock->get_op_array(CCD_CRE_DESCOMP).get_element(i); // P_{ki}
    for (int opind2=0; opind2<opvec2.size(); opind2++) {
      boost::shared_ptr<StackSparseMatrix> op2 = opvec2.at(opind2);
      bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
      if (deallocate2) {
        op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
        op2->build(*otherblock);
      }

	    SpinQuantum op1q = op1->get_deltaQuantum()[0], op2q = op2->get_deltaQuantum()[0], oq = -getSpinQuantum(b->nonactive_orb(0)); // o is the resulted CDD
	    int j2 = op2q.get_s().getirrep(), j1 = op1q.get_s().getirrep(), j21 = oq.get_s().getirrep();
	    int l2 = op2q.get_symm().getirrep(), l1 = op1q.get_symm().getirrep(), l21 = oq.get_symm().getirrep(), l3 = oq.get_symm().getirrep();
	    double factor = dmrginp.spinAdapted() ? pow(-1.0, static_cast<int>((2+j2)/2)) * sixj(j2, j1, j21, 1, 0, j2) * sqrt((j21+1)*(j2+1)) : 1.0;
	    factor *= Symmetry::spatial_sixj(l2, l1, l21, l3, 0, (-IrrepSpace(l2)).getirrep());


      double parity = 1.0;
      if (otherblock == b->get_rightBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);
      //if (loopblock == b->get_leftBlock()) parity *= getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), q);

      SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1, b, c, v, q, factor*parity);
      if (deallocate2) op2->deallocate();

    }
    if (deallocate1) op1->deallocate();
}

void SpinAdapted::stackopxop::CDDandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{

  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();


  boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(OVERLAP, hq);
  bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  if (deallocate2) {
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
  }

  boost::shared_ptr<StackSparseMatrix> op1ham;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op1ham = loopblock->get_op_array(CDD_SUM).get_element(0).at(0);
  else
    op1ham = loopblock->get_op_rep(CDD_SUM, q);
  bool deallocate1ham = op1ham->memoryUsed() == 0 ? true : false;
  if (deallocate1ham) {
    op1ham->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1ham->build(*loopblock);
  }

  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1ham, b, c, v, q, 1.0);	    
  if (deallocate1ham) op1ham->deallocate();


  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  op1 = loopblock->get_op_rep(OVERLAP, hq);
  if (deallocate1) {
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
  }

  boost::shared_ptr<StackSparseMatrix> op2ham;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op2ham = otherblock->get_op_array(CDD_SUM).get_element(0).at(0);
  else
    op2ham = otherblock->get_op_rep(CDD_SUM, q);

  bool deallocate2ham = op2ham->memoryUsed() == 0 ? true : false;
  if (deallocate2ham) {
    op2ham->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2ham->build(*otherblock);
  }
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2ham, *op1, b, c, v, q, 1.0);	    
  if (deallocate2ham) op2ham->deallocate();

  if (deallocate1) op1->deallocate();
  if (deallocate2) op2->deallocate();


}

void SpinAdapted::stackopxop::CCDandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q)
{

  
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const StackSpinBlock* loopblock = (otherblock==b->get_leftBlock()) ? b->get_rightBlock() : b->get_leftBlock();


  boost::shared_ptr<StackSparseMatrix> op2 = otherblock->get_op_rep(OVERLAP, hq);
  bool deallocate2 = op2->memoryUsed() == 0 ? true : false; 
  if (deallocate2) {
    op2->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2->build(*otherblock);
  }

  boost::shared_ptr<StackSparseMatrix> op1ham;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op1ham = loopblock->get_op_array(CCD_SUM).get_element(0).at(0);
  else
    op1ham = loopblock->get_op_rep(CCD_SUM, q);
  bool deallocate1ham = op1ham->memoryUsed() == 0 ? true : false;
  if (deallocate1ham) {
    op1ham->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1ham->build(*loopblock);
  }

  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2, *op1ham, b, c, v, q, 1.0);	    
  if (deallocate1ham) op1ham->deallocate();


  bool deallocate1 = op1->memoryUsed() == 0 ? true : false; 
  op1 = loopblock->get_op_rep(OVERLAP, hq);
  if (deallocate1) {
    op1->allocate(loopblock->get_braStateInfo(), loopblock->get_ketStateInfo());
    op1->build(*loopblock);
  }

  boost::shared_ptr<StackSparseMatrix> op2ham;
  if (!dmrginp.spinAdapted() || dmrginp.hamiltonian() == BCS)
    op2ham = otherblock->get_op_array(CCD_SUM).get_element(0).at(0);
  else
    op2ham = otherblock->get_op_rep(CCD_SUM, q);

  bool deallocate2ham = op2ham->memoryUsed() == 0 ? true : false;
  if (deallocate2ham) {
    op2ham->allocate(otherblock->get_braStateInfo(), otherblock->get_ketStateInfo());
    op2ham->build(*otherblock);
  }
  SpinAdapted::operatorfunctions::TensorMultiply(otherblock, *op2ham, *op1, b, c, v, q, 1.0);	    
  if (deallocate2ham) op2ham->deallocate();
  if (deallocate1) op1->deallocate();
  if (deallocate2) op2->deallocate();


}

