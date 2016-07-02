/*                                                                           
									     Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
									     Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
									     This program is integrated in Molpro with the permission of 
									     Sandeep Sharma and Garnet K.-L. Chan
*/

#include "MatrixBLAS.h"
#include "StackMatrix.h"
#include "Stackspinblock.h"
#include "StackOperators.h"
#include "csf.h"
#include "couplingCoeffs.h"
#include "operatorfunctions.h"
#include "stackopxop.h"
#include "operatorloops.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "SpinQuantum.h"
#include "pario.h"
#include "stackopxop.h"
#include "IntegralMatrix.h"

//******************CRE*****************

void SpinAdapted::StackCre::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  const int i = get_orbs()[0];


  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE).has(i))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE, deltaQuantum, i);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      else {
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }
      return;
    }
  if (rightBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix>& op = rightBlock->get_op_rep(CRE, deltaQuantum, i);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      return;
    }  
  abort();  
}

void SpinAdapted::StackCre::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();
  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());

  const int i = get_orbs()[0];

  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE).has(i))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE, deltaQuantum, i);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
      dmrginp.makeopsT -> stop();
      return;
    }
  if (rightBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix>& op = rightBlock->get_op_rep(CRE, deltaQuantum, i);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  abort();  
}

double SpinAdapted::StackCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int Iirrep = SymmetryOfOrb(get_orbs()[0]).getirrep();;
  IrrepSpace sym = deltaQuantum[0].get_symm();
  bool finish = false;
  bool write = false;
  int Sign = 1;

  TensorOp C(get_orbs()[0], 1);

  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {  
    for (int i=0; i<ladder.size(); i++) {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	double MatElements = calcMatrixElements(c1, C, ladder[i], backupSlater1, backupSlater2, index) ;
	//std::vector<double> MatElements = calcMatrixElements(c1, C, ladder[i], backupSlater1, backupSlater2) ;
        element += MatElements/cleb;
        break;
      }
      else
        continue;
    }
  }
  return element;
}

//******************DES*****************

void SpinAdapted::StackDes::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  const int i = get_orbs()[0];


  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES).has(i))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES, deltaQuantum, i);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      else {
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }
      return;
    }
  if (rightBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix>& op = rightBlock->get_op_rep(DES, deltaQuantum, i);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
    
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      return;
    }  
  abort();  
}

void SpinAdapted::StackDes::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();
  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());

  const int i = get_orbs()[0];

  memset(data, 0, totalMemory * sizeof(double));
 
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES).has(i))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES, deltaQuantum, i);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
      dmrginp.makeopsT -> stop();
      return;
    }
  if (rightBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix>& op = rightBlock->get_op_rep(DES, deltaQuantum, i);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  else
    abort();  

}

double SpinAdapted::StackDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int Iirrep = SymmetryOfOrb(get_orbs()[0]).getirrep();;
  IrrepSpace sym = deltaQuantum[0].get_symm();
  bool finish = false;
  bool write = false;
  int Sign = 1;

  TensorOp D(get_orbs()[0], -1);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int iq=0; iq<deltaQuantum.size(); iq++)
    for (int i=0; i<ladder.size(); i++)
      {
	int index = 0; double cleb=0.0;
	if (nonZeroTensorComponent(c1, deltaQuantum[iq], ladder[i], index, cleb)) {
	  double MatElements = calcMatrixElements(c1, D, ladder[i], backupSlater1, backupSlater2, index) ;
	  element = MatElements/cleb;
	  break;
	}
	else
	  continue;
    
      }

  return element;
}


//******************CREDES*****************

void SpinAdapted::StackCreDes::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();


  if (leftBlock->get_op_array(CRE_DES).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }

      return;
    }
  if (rightBlock->get_op_array(CRE_DES).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      return;
    }  
  if (leftBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      if (rightBlock->has(DES) ) {
	boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(j), j);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }
      else {
	boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
	//op2->set_conjugacy('t');
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
	//op2->set_conjugacy('n');
      }
    }
  else if (rightBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      if (leftBlock->has(DES) ) {
	boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(j), j);
	double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0*parity);
      }
      else {
	boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
	//op2->set_conjugacy('t');
	double parity = getCommuteParity(op1->get_deltaQuantum()[0], -op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	// getCommuteParity doesn't depend on deltaQuantum.get_n()
	SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, 1.0*parity);
	//op2->set_conjugacy('n');
      }
    }
  else
    abort();  
}

void SpinAdapted::StackCreDes::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();
  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());

  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  memset(data, 0, totalMemory * sizeof(double));

  if (leftBlock->get_op_array(CRE_DES).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }

      dmrginp.makeopsT -> stop();
      return;
    }
  if (rightBlock->get_op_array(CRE_DES).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_DES, deltaQuantum, i,j);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  if (leftBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      if (rightBlock->has(DES) ) {
	boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(j), j);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0);
      }
      else {
	boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
	//op2->set_conjugacy('t');
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, 1.0);
	//op2->set_conjugacy('n');
      }
    }
  else if (rightBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      if (leftBlock->has(DES) ) {
	boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(j), j);
	double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
      }
      else {
	boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
	//op2->set_conjugacy('t');
	double parity = getCommuteParity(op1->get_deltaQuantum()[0], -op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	// getCommuteParity doesn't depend on deltaQuantum.get_n()
	SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, 1.0*parity);
	//op2->set_conjugacy('n');
      }
    }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}

void SpinAdapted::StackCreDes::buildUsingCre(const StackSpinBlock* b) {
  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));
  
  if (b->get_op_array(CRE).has(i) && b->get_op_array(CRE).has(j)) {
    const boost::shared_ptr<StackSparseMatrix> op1 = b->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<StackSparseMatrix> op2 = b->get_op_rep(CRE, getSpinQuantum(j), j);
    SpinAdapted::operatorfunctions::Product(b, *op1, Transpose(*op2), *this, 1.0);
  }
}

double SpinAdapted::StackCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();

  TensorOp C(I, 1), D(J, -1);
  TensorOp CD = C.product(D, spin, irrep);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
      {
	int index = 0; double cleb=0.0;
	if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	  double MatElements = calcMatrixElements(c1, CD, ladder[i], backupSlater1, backupSlater2, index) ;
	  element += MatElements/cleb;
	  break;
	}
	else
	  continue;
      }
  }
  return element;
}


//******************DESCRE*****************

void SpinAdapted::StackDesCre::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();


  if (leftBlock->get_op_array(DES_CRE).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }

      return;
    }
  if (rightBlock->get_op_array(DES_CRE).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      return;
    }  
  if (leftBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(i), i);
      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
    }
  else if (rightBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(DES, -getSpinQuantum(i), i);
      boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      double parity = getCommuteParity(-op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
      SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0*parity);
    }
  else
    abort();  
}

void SpinAdapted::StackDesCre::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_CRE).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }

      dmrginp.makeopsT -> stop();
      return;
    }
  if (rightBlock->get_op_array(DES_CRE).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(DES_CRE, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  if (leftBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op2, *op1, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  else if (rightBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(DES, -getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);    
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
    }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}

double SpinAdapted::StackDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();

  TensorOp D(I, -1), C(J, 1);
  
  TensorOp CD = C.product(D, spin, irrep);

  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
      {
	int index = 0; double cleb=0.0;
	if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	  double MatElements = calcMatrixElements(c1, CD, ladder[i], backupSlater1, backupSlater2, index) ;
	  element += MatElements/cleb;
	  break;
	}
	else
	  continue;
      }
  }
  return element;
}

//******************CRECRE*****************

void SpinAdapted::StackCreCre::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }
      return;
    }
  if (rightBlock->get_op_array(CRE_CRE).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  if (leftBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
    }
  else if (rightBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
      SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0*parity);
    }
  else
    abort();  
}


void SpinAdapted::StackCreCre::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
      dmrginp.makeopsT -> stop();
      return;
    }
  if (rightBlock->get_op_array(CRE_CRE).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_CRE, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  if (leftBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  else if (rightBlock->get_op_array(CRE).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(j), j);
      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
    }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}

void SpinAdapted::StackCreCre::buildUsingCre(const StackSpinBlock* b) {
  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));
  
  if (b->get_op_array(CRE).has(i) && b->get_op_array(CRE).has(j)) {
    const boost::shared_ptr<StackSparseMatrix> op1 = b->get_op_rep(CRE, getSpinQuantum(i), i);
    const boost::shared_ptr<StackSparseMatrix> op2 = b->get_op_rep(CRE, getSpinQuantum(j), j);
    SpinAdapted::operatorfunctions::Product(b, *op1, *op2, *this, 1.0);
  }
}

double SpinAdapted::StackCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();


  TensorOp C1(I,1), C2(J,1);
  TensorOp CC = C1.product(C2, spin, sym.getirrep(), I==J);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {  
    for (int i=0; i<ladder.size(); i++)
      {
	int index = 0; double cleb=0.0;
	if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	  double MatElements = calcMatrixElements(c1, CC, ladder[i], backupSlater1, backupSlater2, index) ;
	  element += MatElements/cleb;
	  break;
	}
	else
	  continue;
      }
  }
  return element;
}

//******************DESDES*****************

void SpinAdapted::StackDesDes::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();
  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());
  Sign = 1;

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_DES).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_DES, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
      dmrginp.makeopsT -> stop();
      return;
    }
  if (rightBlock->get_op_array(DES_DES).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(DES_DES, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      dmrginp.makeopsT -> stop();
      return;
    }  
  if (leftBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(j), j);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  else if (rightBlock->get_op_array(DES).has(i))
    {
      const boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(DES, -getSpinQuantum(i), i);
      const boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(j), j);
      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, 1.0*parity);
    }
  else
    abort();  
  dmrginp.makeopsT -> stop();
}



double SpinAdapted::StackDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int irrep = deltaQuantum[0].get_symm().getirrep();
  int spin = deltaQuantum[0].get_s().getirrep();


  TensorOp D1(I,-1), D2(J,-1);
  TensorOp DD = D1.product(D2, spin, sym.getirrep(), I==J);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {  
    for (int i=0; i<ladder.size(); i++)
      {
	int index = 0; double cleb=0.0;
	if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	  double MatElements = calcMatrixElements(c1, DD, ladder[i], backupSlater1, backupSlater2, index) ;
	  element += MatElements/cleb;
	  break;
	}
	else
	  continue;
      }
  }
  return element;
}

//******************CREDESCOMP*****************


void SpinAdapted::StackCreDesComp::buildUsingCre(const StackSpinBlock* b) {
  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));
 
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
 
  TensorOp C(i,1), D(j,-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // th
  for (int kx = 0; kx < b->get_sites().size(); ++kx)
    for (int lx = 0; lx < b->get_sites().size(); ++lx) {
      int k = b->get_sites()[kx];
      int l = b->get_sites()[lx];

      if (!(b->get_op_array(CRE).has(k) && b->get_op_array(CRE).has(l))) continue;

      if (dmrginp.hamiltonian() == BCS) {
	// CC
	{
	  TensorOp CK(k,1), CL(l,1);
	  TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
	  if (!CC2.empty) {
	    double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd[b->get_integralIndex()]);
	    const boost::shared_ptr<StackSparseMatrix> op1 = b->get_op_rep(CRE, getSpinQuantum(k), k);
	    const boost::shared_ptr<StackSparseMatrix> op2 = b->get_op_rep(CRE, getSpinQuantum(l), l);
	    SpinAdapted::operatorfunctions::Product(b, *op1, *op2, *this, scaleV);
	  }
	}
	// DD
	{
	  TensorOp DK(k,-1), DL(l,-1);
	  TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
	  if (!DD2.empty) {
	    double scaleV = calcCompfactor(CD1, DD2, CD, v_cccd[b->get_integralIndex()]);
	    const boost::shared_ptr<StackSparseMatrix> op1 = b->get_op_rep(CRE, getSpinQuantum(k), k);
	    const boost::shared_ptr<StackSparseMatrix> op2 = b->get_op_rep(CRE, getSpinQuantum(l), l);
	    SpinAdapted::operatorfunctions::Product(b, Transpose(*op1), Transpose(*op2), *this, scaleV);
	  }
	}
      }

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	double scaleV = calcCompfactor(CD1, CD2, CD,*(b->get_twoInt()), b->get_integralIndex());
	const boost::shared_ptr<StackSparseMatrix> op1 = b->get_op_rep(CRE, getSpinQuantum(k), k);
	const boost::shared_ptr<StackSparseMatrix> op2 = b->get_op_rep(CRE, getSpinQuantum(l), l);
	SpinAdapted::operatorfunctions::Product(b, *op1, Transpose(*op2), *this, scaleV);
      }
    }

}


void SpinAdapted::StackCreDesComp::buildfromCreDes(StackSpinBlock& b) 
{
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  memset(data, 0, totalMemory * sizeof(double));

  //StackCreDesComp op_array = this; 
  //initiateMultiThread(this, op_array, numthrds);

  TensorOp C(i,1), D(j,-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops1;
  int numCD = 0, numDC = 0;
  std::vector<double> scaleCD, scaleCC;
  if (dmrginp.hamiltonian() == BCS) scaleCC.resize(b.get_op_array(CRE_CRE).get_size()*2, 0.0);
  //for each CD there is spin0, spin 1 and normal and transpose
  if (dmrginp.spinAdapted()) scaleCD.resize(b.get_op_array(CRE_DES).get_size()*8, 0.0);
  else scaleCD.resize(b.get_op_array(CRE_DES).get_size()*4, 0.0);
  for (int ii=0; ii<b.get_op_array(CRE_DES).get_size(); ii++)
    for (int ji=0; ji<b.get_op_array(CRE_DES).get_local_element(ii).size(); ji++) 
      if (b.get_op_array(CRE_DES).get_local_element(ii)[ji]->get_deltaQuantum(0) == deltaQuantum[0] || 
          b.get_op_array(CRE_DES).get_local_element(ii)[ji]->get_deltaQuantum(0) == -deltaQuantum[0]) {
	allops1.push_back(b.get_op_array(CRE_DES).get_local_element(ii)[ji]);
	const int k = allops1.back()->get_orbs()[0];
	const int l = allops1.back()->get_orbs()[1];
	
	TensorOp CK(k,1), DL(l,-1);
	TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
	if (!CD2.empty) 
	  scaleCD[2*(allops1.size()-1)] = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	
	CK=TensorOp(l,1); DL=TensorOp(k,-1);      
	CD2 = CK.product(DL, spin, sym.getirrep());
	if (!b.has(DES) && l!=k && !CD2.empty) {
	  double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	  if (dmrginp.spinAdapted()) {
	    scaleV *= getCommuteParity(getSpinQuantum(l), getSpinQuantum(k), get_deltaQuantum()[0]);
	  }
	  scaleCD[2*(allops1.size()-1) + 1] = scaleV;
	}
      }
  numCD = allops1.size();

  if(b.has(DES)) {
    for (int ii=0; ii<b.get_op_array(DES_CRE).get_size(); ii++)
      for (int ji=0; ji<b.get_op_array(DES_CRE).get_local_element(ii).size(); ji++) 
	if (b.get_op_array(DES_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0) == -deltaQuantum[0]) {
	  allops1.push_back(b.get_op_array(DES_CRE).get_local_element(ii)[ji]);
	  
	  const int l = allops1.back()->get_orbs()[0];
	  const int k = allops1.back()->get_orbs()[1];
	  
	  if (k!=l) {
	    TensorOp CK(k,1), DL(l,-1);      
	    TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
	    if (!CD2.empty) {
	      double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	      if (dmrginp.spinAdapted()) {
		scaleV *= getCommuteParity(getSpinQuantum(l), getSpinQuantum(k), get_deltaQuantum()[0]);
	      }
	      scaleCD[2*(allops1.size()-1)] = scaleV;
	    }
	  }
	}
  }
  numDC = allops1.size();

  if (dmrginp.hamiltonian() == BCS) {
    for (int ii = 0; ii < b.get_op_array(CRE_CRE).get_size(); ++ii)
      for (int ji = 0; ji < b.get_op_array(CRE_CRE).get_local_element(ii).size(); ++ji)
        if (b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0).get_s() == deltaQuantum[0].get_s() || 
            b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0).get_s() == -deltaQuantum[0].get_s()) {
          allops1.push_back(b.get_op_array(CRE_CRE).get_local_element(ii)[ji]);

          const int k = allops1.back()->get_orbs()[0];
          const int l = allops1.back()->get_orbs()[1];
          TensorOp CK(k,1), CL(l,1);
          TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k == l);

          if (!CC2.empty) {
            // also consider CL*CK
            CL = TensorOp(l,1);
            CK = TensorOp(k,1);
            TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
            int parity = getCommuteParity(getSpinQuantum(l),getSpinQuantum(k),get_deltaQuantum()[0]);              
            scaleCC[2*(allops1.size()-numDC-1)] = calcCompfactor(CD1, CC2, CD, v_cccd[b.get_integralIndex()]);
            if (k != l)
              scaleCC[2*(allops1.size()-numDC-1)] += parity * calcCompfactor(CD1, CC2_commute, CD, v_cccd[b.get_integralIndex()]);
          }

          // now consider DL*DK
          TensorOp DL(l,-1), DK(k,-1);
          TensorOp DD2 = DL.product(DK, spin, sym.getirrep(), k == l);

          if (!b.has(DES) && !DD2.empty) {
            DK = TensorOp(k,-1);
            DL = TensorOp(l,-1);
            TensorOp DD2_commute = DK.product(DL, spin, sym.getirrep(), k == l);
            int parity = getCommuteParity(getSpinQuantum(l),getSpinQuantum(k),get_deltaQuantum()[0]);
            scaleCC[2*(allops1.size()-numCD-1)+1] = calcCompfactor(CD1, DD2, CD, v_cccd[b.get_integralIndex()]);
            if (k != l)
              scaleCC[2*(allops1.size()-numCD-1)+1] += parity * calcCompfactor(CD1, DD2_commute, CD, v_cccd[b.get_integralIndex()]);
          }
        }

    if (b.has(DES)) {
      pout << "buildfromCreDes with DES in BCS not implemented" << endl;
      abort();
    }
  }

  const int quantaSz = b.get_braStateInfo().quanta.size () * b.get_ketStateInfo().quanta.size();
  std::multimap<long, std::pair<int, int> > reorder;
  for (int i=0; i<b.get_braStateInfo().quanta.size(); i++)
    for (int j=0; j<b.get_ketStateInfo().quanta.size(); j++)
      reorder.insert( std::pair<long, pair<int,int> >(b.get_braStateInfo().getquantastates(i)*b.get_ketStateInfo().getquantastates(j), std::pair<int, int>(i, j) ));

  std::vector<pair<int, int> > reorderedVector(quantaSz);
  int index = quantaSz-1;
  for (std::multimap<long, std::pair<int,int> >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
    reorderedVector[index] = it->second;
    index--;
  }  

  //#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii<allops1.size()*reorderedVector.size(); ii++) {
    int opindex = (ii)%allops1.size(), quantaindex = (ii)/allops1.size();
    int lQ = reorderedVector[quantaindex].first, rQ = reorderedVector[quantaindex].second;

    const int k = allops1[opindex]->get_orbs()[0];
    const int l = allops1[opindex]->get_orbs()[1];
    //this is a cd operator
    if (opindex <numCD) {
      if (allops1[opindex]->allowed(lQ, rQ) && allowed(lQ, rQ) && fabs(scaleCD[2*opindex]) >TINY) 
	MatrixScaleAdd(scaleCD[2*opindex], allops1[opindex]->operator_element(lQ, rQ), this->operator_element(lQ, rQ));
      
      if (!b.has(DES) && l!=k && fabs(scaleCD[2*opindex+1]) > TINY) {
	if (allowed(lQ, rQ) && allops1[opindex]->allowed(rQ, lQ)) {
	  double scaling = getStandAlonescaling(-(allops1[opindex]->get_deltaQuantum(0)), b.get_braStateInfo().quanta[lQ], b.get_ketStateInfo().quanta[rQ]);
	  int nrows = operator_element(lQ, rQ).Nrows();
	  int ncols = operator_element(lQ, rQ).Ncols();
	  for (int row=0; row<nrows; row++)
	    DAXPY(ncols, scaling*scaleCD[2*opindex+1], &(allops1[opindex]->operator_element(rQ, lQ)(1, row+1)), nrows, &(this->operator_element(lQ, rQ)(row+1, 1)), 1); 
	}
      }
    } else if (opindex < numDC) {
      if (k!=l && allops1[opindex]->allowed(lQ, rQ) && allowed(lQ, rQ) && fabs(scaleCD[2*opindex]) > TINY)  
	MatrixScaleAdd(scaleCD[2*opindex], allops1[opindex]->operator_element(lQ, rQ), this->operator_element(lQ, rQ));
    } else { // CC
      if (allops1[opindex]->allowed(lQ, rQ) && allowed(lQ, rQ) && fabs(scaleCC[2*(opindex-numDC)]) > TINY)
        MatrixScaleAdd(scaleCC[2*(opindex-numDC)], allops1[opindex]->operator_element(lQ, rQ), this->operator_element(lQ, rQ));
      
      if (!b.has(DES) &&  fabs(scaleCC[2*(opindex-numDC)+1]) > TINY) {
        if (allowed(lQ, rQ) && allops1[opindex]->allowed(rQ, lQ)) {
	  double scaling = getStandAlonescaling(-(allops1[opindex]->get_deltaQuantum(0)), b.get_braStateInfo().quanta[lQ], b.get_ketStateInfo().quanta[rQ]);
          int nrows = operator_element(lQ, rQ).Nrows();
          int ncols = operator_element(lQ, rQ).Ncols();
          for (int row=0; row<nrows; ++row)
            DAXPY(ncols, scaling*scaleCC[2*(opindex-numDC)+1], &(allops1[opindex]->operator_element(rQ, lQ)(1, row+1)), nrows, &(this->operator_element(lQ, rQ)(row+1, 1)), 1);
        }
      }
    }
  }

  //accumulateMultiThread(this, op_array, numthrds);
}


void SpinAdapted::StackCreDesComp::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];


  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DESCOMP).has(i, j))
    { 
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }

    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DESCOMP).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
    }  

  const TwoElectronArray& v2 = *(b.get_twoInt());
  // build CDcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      int spink = dmrginp.spin_orbs_symmetry()[k], spinl = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[l]);
      if (dmrginp.spinAdapted()) {
	spink = dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[k]], \
	spinl = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[l]]);
      }
      if (Symmetry::negativeofAbelian(sym.getirrep()) == Symmetry::addAbelian(spink, spinl)) {
	//if ((-sym).getirrep() == Symmetry::addAbelian(spink, spinl)) {
	double scaleV = calcCompfactor(i,j,k,l, spin, CD,*(b.get_twoInt()), b.get_integralIndex());

	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	  if (rightBlock->has(DES)) {
	    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	    SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
	  } 
	  else {
	    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	    SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
	  }
	}
      }
      
      spinl = dmrginp.spin_orbs_symmetry()[l]; spink = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[k]);
      if (dmrginp.spinAdapted()) {
	spinl = dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[l]], \
	spink = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[k]]);
      }

      if (Symmetry::negativeofAbelian(sym.getirrep()) == Symmetry::addAbelian(spink, spinl)) {
	double scaleV = calcCompfactor(i,j,l,k, spin, CD,*(b.get_twoInt()), b.get_integralIndex());
	
      	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	  if (rightBlock->has(DES) && leftBlock->has(DES)) {
	    boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV*parity);
	  } else {
	    //StackTransposeview top2 = StackTransposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	    boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    //op2->set_conjugacy('t');
	    SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, scaleV*parity);
	    //op2->set_conjugacy('n');
	  }
	}
      }
      
      if (dmrginp.hamiltonian() == BCS) {
	TensorOp C(i,1), D(j,-1);
	TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried
        TensorOp CK = TensorOp(k, 1);
        TensorOp CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l); // k cannot equal to l
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd[b.get_integralIndex()]);
          CL = TensorOp(l, 1);
          CK = TensorOp(k, 1);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, CC2_commute, CD, v_cccd[b.get_integralIndex()]);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[1]);
            scaleV += parity * scaleV2;
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
            }
          }
        }
        TensorOp DK(k, -1);
        TensorOp DL = TensorOp(l, -1);
        TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
        if (!DD2.empty) {
          double scaleV = calcCompfactor(CD1, DD2, CD, v_cccd[b.get_integralIndex()]);
          DL = TensorOp(l, -1);
          DK = TensorOp(k, -1);
          TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, DD2_commute, CD, v_cccd[b.get_integralIndex()]);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            if (leftBlock->has(DES)) {
              boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
              boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
              double parity =  getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
              scaleV += parity * scaleV2;
              if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
                SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
              }
            } else {
              boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
              boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      //op1->set_conjugacy('t'); op2->set_conjugacy('t');
              double parity = getCommuteParity(-op1->get_deltaQuantum()[0], -op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
              scaleV += parity * scaleV2;
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
                SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, Transpose(*op1), Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
              }
	      //op1->set_conjugacy('n'); op2->set_conjugacy('n');
            }
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();
}


void SpinAdapted::StackCreDesComp::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

 
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  memset(data, 0, totalMemory * sizeof(double));


  TensorOp C(i,1), D(j,-1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DESCOMP).has(i, j))
    { 
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }

    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DESCOMP).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_DESCOMP, deltaQuantum, i,j);
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  

  // build CDcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k,1), DL(l,-1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	  if (rightBlock->has(DES)) {
	    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	  } 
	  else {
	    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, scaleV);
	  }
	}
      }
      
      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	
      	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	  if (rightBlock->has(DES) && leftBlock->has(DES)) {
	    boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	  } 
	  else {
	    boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, scaleV*parity);
	  }
	}
      }
      
      if (dmrginp.hamiltonian() == BCS) {
        CK = TensorOp(k, 1);
        TensorOp CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l); // k cannot equal to l
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd[b.get_integralIndex()]);
          CL = TensorOp(l, 1);
          CK = TensorOp(k, 1);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, CC2_commute, CD, v_cccd[b.get_integralIndex()]);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[1]);
            scaleV += parity * scaleV2;
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        TensorOp DK(k, -1);
        DL = TensorOp(l, -1);
        TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
        if (!DD2.empty) {
          double scaleV = calcCompfactor(CD1, DD2, CD, v_cccd[b.get_integralIndex()]);
          DL = TensorOp(l, -1);
          DK = TensorOp(k, -1);
          TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, DD2_commute, CD, v_cccd[b.get_integralIndex()]);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            if (leftBlock->has(DES)) {
              boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
              boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
              double parity =  getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
              scaleV += parity * scaleV2;
              if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
                SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
              }
            } 
	    else {
              boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
              boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      //op1->set_conjugacy('t'); op2->set_conjugacy('t');
              double parity = getCommuteParity(-op1->get_deltaQuantum()[0], -op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
              scaleV += parity * scaleV2;
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
                SpinAdapted::operatorfunctions::TensorProduct(leftBlock, Transpose(*op1), Transpose(*op2), &b, &(b.get_stateInfo()), *this, scaleV);
              }
	      //op1->set_conjugacy('n'); op2->set_conjugacy('n');
            }
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::StackCreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is(); // classify whether we calculate CC or CD

  TensorOp C(I,1), D(J,-1);

  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep());
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  if (dmrginp.hamiltonian() != BCS && c1.n_is() != ladder[0].n_is()) return 0.0;
  if (dmrginp.hamiltonian() == BCS && abs(dn) > 2) return 0.0;
  if (dmrginp.spinAdapted() && (c1.S_is().getirrep() > ladder[0].S_is().getirrep()+spin || c1.S_is().getirrep() <ladder[0].S_is().getirrep() -spin)) return 0.0;

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {
      if (!dmrginp.spinAdapted() && c1.S_is().getirrep() != ladder[i].S_is().getirrep()+spin) continue;
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];

	    
	    bool isZero = true;
            if (dmrginp.hamiltonian() == BCS) { // BCS
              if (ladder[i].det_rep.size() > 1) {
                pout << "StackCreDesComp::redMatrixElement failed" << endl;
                abort();
              }
              if (dn == 0) {
                for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
                isZero = true;
                for (auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[l] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
              } else if (dn == 2) {
                if (k == l) continue;
                for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 1 && s1.get_orbstring().get_occ_rep()[l] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
                isZero = true;
                for (auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 0 && s1.get_orbstring().get_occ_rep()[l] == 0) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
              } else if (dn == -2) {
                if (k == l) continue;
                for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 0 && s1.get_orbstring().get_occ_rep()[l] == 0) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
                isZero = true;
                for (auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 1 && s1.get_orbstring().get_occ_rep()[l] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
              }
            } else if (dmrginp.spinAdapted()) { // spin-adapted
	      for (auto it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
		const Slater &s1 = it1->first;
		if ( s1.get_orbstring().get_occ_rep()[2*k] == 1 || s1.get_orbstring().get_occ_rep()[2*k+1] == 1) {
		  isZero = false;
		  break;
		}
	      }
	      if (isZero) continue;

	      isZero = true;
	      for (auto it1 = ladder[i].det_rep.begin(); it1!= ladder[i].det_rep.end(); it1++) {
		const Slater &s1 = it1->first;
		if ( s1.get_orbstring().get_occ_rep()[2*l] == 1 || s1.get_orbstring().get_occ_rep()[2*l+1] == 1) {
		  isZero = false;
		  break;
		}
	      }
	      if (isZero) continue;
            } else { // non-spinadapted
              if (ladder[i].det_rep.size() > 1) {
                pout << "StackCreDesComp::redMatrixElement failed" << endl;
                abort();
              }
              for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                const Slater &s1 = it1->first;
                if (s1.get_orbstring().get_occ_rep()[k] == 1) {
                  isZero = false;
                  break;
                }
              }
              if (isZero) continue;
              
              isZero = true;
              for(auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                const Slater &s1 = it1->first;
                if (s1.get_orbstring().get_occ_rep()[l] == 1) {
                  isZero = false;
                  break;
                }
              }
              if (isZero) continue;
            }
            
            if (dmrginp.hamiltonian() == BCS && dn == 2) {
              TensorOp CK(k,1), CL(l,1);
              TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
              if (!CC2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i], backupSlater1, backupSlater2) ;
                double factor = calcCompfactor(CD1, CC2, CD, v_cccd[b->get_integralIndex()]);
                element += MatElements[index]*factor/cleb;
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == -2) {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
              if (!DD2.empty) {
                std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i], backupSlater1, backupSlater2) ;
                double factor = calcCompfactor(CD1, DD2, CD, v_cccd[b->get_integralIndex()]);
                element += MatElements[index]*factor/cleb;
              }
            } else {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
		double MatElements = calcMatrixElements(c1, CD2, ladder[i], backupSlater1, backupSlater2, index) ;
                double factor = calcCompfactor(CD1, CD2, CD, *(b->get_twoInt()), b->get_integralIndex());
                element += MatElements*factor/cleb;
              }
            }
          }
        break;
      }
      else
        continue;
    }
  }
  return element;
}



//******************DESCRECOMP*****************

void SpinAdapted::StackDesCreComp::buildfromDesCre(StackSpinBlock& b) 
{
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];

  memset(data, 0, totalMemory * sizeof(double));

  //StackDesCreComp* op_array; 
  //initiateMultiThread(this, op_array, numthrds);


  TensorOp D(i,-1), C(j,1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops1, allops2;
  for (int ii=0; ii<b.get_op_array(CRE_DES).get_size(); ii++)
    for (int ji=0; ji<b.get_op_array(CRE_DES).get_local_element(ii).size(); ji++) 
      if (b.get_op_array(CRE_DES).get_local_element(ii)[ji]->get_deltaQuantum(0) == -deltaQuantum[0])
	allops1.push_back(b.get_op_array(CRE_DES).get_local_element(ii)[ji]);

  for (int ii=0; ii<b.get_op_array(DES_CRE).get_size(); ii++)
    for (int ji=0; ji<b.get_op_array(DES_CRE).get_local_element(ii).size(); ji++) 
      if (b.get_op_array(DES_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0) == deltaQuantum[0])
	allops2.push_back(b.get_op_array(DES_CRE).get_local_element(ii)[ji]);

  //#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii<allops1.size(); ii++) {
    const int k = allops1[ii]->get_orbs()[0];
    const int l = allops1[ii]->get_orbs()[1];

    TensorOp CK(k,1), DL(l,-1);      
    TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
    if (!CD2.empty) {
      double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
      ScaleAdd(scaleV, *allops1[ii], *this);
    }

  }
  
  //#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii<allops2.size(); ii++) {
    const int l = allops2[ii]->get_orbs()[0];
    const int k = allops2[ii]->get_orbs()[1];

    if (k==l) continue;
    TensorOp CK(k,1), DL(l,-1);      
    TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
    if (!CD2.empty) {
      double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
      //scaleV *= getCommuteParity(getSpinQuantum(i), getSpinQuantum(j), get_deltaQuantum()[0]);
      scaleV *= getCommuteParity(getSpinQuantum(l), -getSpinQuantum(k), get_deltaQuantum()[0]);
      ScaleAdd(scaleV, *allops2[ii], *this);
    }
  }

  //accumulateMultiThread(this, op_array, numthrds);
}


void SpinAdapted::StackDesCreComp::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  TensorOp D(i,-1), C(j,1);
  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep()); // the operator to be complimentaried

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_CRECOMP).has(i, j))
    { 
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_CRECOMP, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }

    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(DES_CRECOMP).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(DES_CRECOMP, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  
  // build DCcomp explicitely
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx) {
      int k = leftBlock->get_sites()[kx];
      int l = rightBlock->get_sites()[lx];

      TensorOp CK(k, 1), DL(l, -1);      
      TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
	if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(DES).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	  boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	  SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	}
      }

      CK=TensorOp(l,1); DL=TensorOp(k,-1);      
      CD2 = CK.product(DL, spin, sym.getirrep());
      if (!CD2.empty) {
      	double scaleV = calcCompfactor(CD1, CD2, CD,*(b.get_twoInt()), b.get_integralIndex());
      	if (leftBlock->get_op_array(DES).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	  boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	  boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	  double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	  SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	}
      }

      if (dmrginp.hamiltonian() == BCS) {
        CK = TensorOp(k, 1);
        TensorOp CL(l, 1);
        TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l); // k cannot equal to l
        if (!CC2.empty) {
          double scaleV = calcCompfactor(CD1, CC2, CD, v_cccd[b.get_integralIndex()]);
          CL = TensorOp(l, 1);
          CK = TensorOp(k, 1);
          TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, CC2_commute, CD, v_cccd[b.get_integralIndex()]);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
            boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
            double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[1]);
            scaleV += parity * scaleV2;
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
        TensorOp DK(k, -1);
        DL = TensorOp(l, -1);
        TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
        if (!DD2.empty) {
          double scaleV = calcCompfactor(CD1, DD2, CD, v_cccd[b.get_integralIndex()]);
          DL = TensorOp(l, -1);
          DK = TensorOp(k, -1);
          TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
          double scaleV2 = calcCompfactor(CD1, DD2_commute, CD, v_cccd[b.get_integralIndex()]);
          if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) &&  fabs(scaleV2)+fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
            boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
            boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
            double parity =  getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[2]);
            scaleV += parity * scaleV2;
            if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
              SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
            }
          }
        }
      }
    }
  dmrginp.makeopsT -> stop();
}


double SpinAdapted::StackDesCreComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is(); // classify whether we calculate CC or CD

  TensorOp D(I,-1), C(J, 1);

  TensorOp CD1 = C.product(D, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep());
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++)
      {
	int index = 0; double cleb=0.0;
	if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	  for (int kl =0; kl<b->get_sites().size(); kl++) 
	    for (int kk =0; kk<b->get_sites().size(); kk++) {

	      int k = b->get_sites()[kk];
	      int l = b->get_sites()[kl];
            
	      if (dmrginp.hamiltonian() == BCS && dn == 2) {
		TensorOp CK(k,1), CL(l,1);
		TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
		if (!CC2.empty) {
		  std::vector<double> MatElements = calcMatrixElements(c1, CC2, ladder[i], backupSlater1, backupSlater2) ;
		  double factor = calcCompfactor(CD1, CC2, CD, v_cccd[b->get_integralIndex()]);
		  element += MatElements[index]*factor/cleb;
		}
	      } else if (dmrginp.hamiltonian() == BCS && dn == -2) {
		TensorOp DK(k,-1), DL(l,-1);
		TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
		if (!DD2.empty) {
		  std::vector<double> MatElements = calcMatrixElements(c1, DD2, ladder[i], backupSlater1, backupSlater2) ;
		  double factor = calcCompfactor(CD1, DD2, CD, v_cccd[b->get_integralIndex()]);
		  element += MatElements[index]*factor/cleb;
		}
	      } else {
		TensorOp CK(k, 1), DL(l, -1);
		TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
		if (!CD2.empty) {
		  double MatElements = calcMatrixElements(c1, CD2, ladder[i], backupSlater1, backupSlater2, index) ;
		  double factor = calcCompfactor(CD1, CD2, CD, *(b->get_twoInt()), b->get_integralIndex());
		  element += MatElements*factor/cleb;  // FIXME a factor of half?
		}
	      }
	    }
	  break;
	}
	else
	  continue;
      }
  }
  return element;
}



//******************DESDESCOMP*****************

void SpinAdapted::StackDesDesComp::buildfromDesDes(StackSpinBlock& b) 
{
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();
  
  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  //StackDesDesComp* op_array; 
  //initiateMultiThread(this, op_array, numthrds);
  
  TensorOp C(i,1), C2(j,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);
  
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  int numCC = 0;
  std::vector<double> scaleDD, scaleCD, scaleCC;
  // TODO resize scaleCC and scaleCD
  if (dmrginp.hamiltonian() == BCS) {
    scaleCC.resize(b.get_op_array(CRE_CRE).get_size(), 0.);
    scaleCD.resize(b.get_op_array(CRE_DES).get_size()*2, 0.);
  }
  if (dmrginp.spinAdapted()) scaleDD.resize(b.get_op_array(CRE_CRE).get_size()*2, 0.);
  else scaleDD.resize(b.get_op_array(CRE_CRE).get_size(), 0.);

  if(!b.has(DES))  {
    for (int ii=0; ii<b.get_op_array(CRE_CRE).get_size(); ii++)
      for (int ji=0; ji<b.get_op_array(CRE_CRE).get_local_element(ii).size(); ji++) {
        bool nonzero = false;

        if (b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0) == -deltaQuantum[0]) {
          const int k = b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_orbs()[0];
          const int l = b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_orbs()[1]; 

          TensorOp DK(k,-1), DL(l,-1);
          TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
          if (!DD2.empty) {
            nonzero = true;
            double parity = -1;
            scaleDD[allops.size()] = parity * calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());
            DK = TensorOp(k,-1); DL = TensorOp(l,-1);
            DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
            double parity2 = getCommuteParity(-getSpinQuantum(k), -getSpinQuantum(l), get_deltaQuantum()[0]);
            if (k != l)
              scaleDD[allops.size()] += parity * parity2 * calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());
          }
        }

        if (dmrginp.hamiltonian() == BCS && b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0).get_s() == deltaQuantum[0].get_s()) { // then we need CC
          const int k = b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_orbs()[0];
          const int l = b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_orbs()[1]; 
          TensorOp CK(k,1), CL(l,1);
          TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k == l);
          if (!CC2.empty) {
            nonzero = true;
            scaleCC[allops.size()] = calcCompfactor(CC1, CC2, DD, v_cccc[b.get_integralIndex()]);
            CK = TensorOp(k,1); CL = TensorOp(l,1);
            CC2 = CL.product(CK, spin, sym.getirrep(), k == l);
            double parity = getCommuteParity(getSpinQuantum(k), getSpinQuantum(l), get_deltaQuantum()[0]);
            if (k != l)
              scaleCC[allops.size()] += parity * calcCompfactor(CC1, CC2, DD, v_cccc[b.get_integralIndex()]);
          }
        }

        if (nonzero) allops.push_back(b.get_op_array(CRE_CRE).get_local_element(ii)[ji]);
      }
  } else {
    if (dmrginp.hamiltonian() == BCS) {
      pout << "buildfromDesDes with DES in BCS not implemented" << endl;
      abort();
    }
    for (int ii=0; ii<b.get_op_array(DES_DES).get_size(); ii++)
      for (int ji=0; ji<b.get_op_array(DES_DES).get_local_element(ii).size(); ji++) 
        if (b.get_op_array(DES_DES).get_local_element(ii)[ji]->get_deltaQuantum(0) == deltaQuantum[0]) {
          allops.push_back(b.get_op_array(DES_DES).get_local_element(ii)[ji]);
          const int k = allops.back()->get_orbs()[0];
          const int l = allops.back()->get_orbs()[1];
          TensorOp DK(k,-1), DL(l,-1);
          TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
          if (!DD2.empty) {
            scaleDD[allops.size()-1] = calcCompfactor(CC1, DD2, DD,  *(b.get_twoInt()), b.get_integralIndex());
            DK = TensorOp(k,-1); DL = TensorOp(l, -1);
            TensorOp DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
            double parity = getCommuteParity(-getSpinQuantum(k), -getSpinQuantum(l), get_deltaQuantum()[0]);
            if (k != l)
              scaleDD[allops.size()-1] += parity * calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());
          }
        }
  }
  numCC = allops.size();

  if (dmrginp.hamiltonian() == BCS) {
    for (int ii=0; ii<b.get_op_array(CRE_DES).get_size(); ii++)
      for (int ji=0; ji<b.get_op_array(CRE_DES).get_local_element(ii).size(); ji++)
        if (b.get_op_array(CRE_DES).get_local_element(ii)[ji]->get_deltaQuantum(0).get_s() == deltaQuantum[0].get_s() || 
	    b.get_op_array(CRE_DES).get_local_element(ii)[ji]->get_deltaQuantum(0).get_s() == -deltaQuantum[0].get_s()) {
          allops.push_back(b.get_op_array(CRE_DES).get_local_element(ii)[ji]);
          const int k = allops.back()->get_orbs()[0];
          const int l = allops.back()->get_orbs()[1];

          TensorOp CK(k,1), DL(l,-1);
          TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
          if (!CD2.empty)
            scaleCD[2*(allops.size()-1)] = calcCompfactor(CC1, CD2, DD, v_cccd[b.get_integralIndex()]);

          TensorOp CL(l,1), DK(k,-1);
          CD2 = CL.product(DK, spin, sym.getirrep());
          if (k!=l && !CD2.empty)
            scaleCD[2*(allops.size()-1)+1] = calcCompfactor(CC1, CD2, DD, v_cccd[b.get_integralIndex()]);
        }
  }

  const int quantaSz = b.get_braStateInfo().quanta.size () * b.get_ketStateInfo().quanta.size();
  std::multimap<long, std::pair<int, int> > reorder;
  for (int i=0; i<b.get_braStateInfo().quanta.size(); i++)
    for (int j=0; j<b.get_ketStateInfo().quanta.size(); j++)
      reorder.insert( std::pair<long, pair<int,int> >(b.get_braStateInfo().getquantastates(i)*b.get_ketStateInfo().getquantastates(j), std::pair<int, int>(i, j) ));

  std::vector<pair<int, int> > reorderedVector(quantaSz);
  int index = quantaSz-1;
  for (std::multimap<long, std::pair<int,int> >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
    reorderedVector[index] = it->second;
    index--;
  }

  //#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii < allops.size()*reorderedVector.size(); ++ii) {
    int opindex = (ii)%allops.size(), quantaindex = (ii)/allops.size();
    int lQ = reorderedVector[quantaindex].first, rQ = reorderedVector[quantaindex].second;

    if (opindex < numCC) {
      if (!b.has(DES)) { // operator is CC
        if (allowed(lQ, rQ) && allops[opindex]->allowed(rQ,lQ) && fabs(scaleDD[opindex]) > TINY) {
          double scaling = getStandAlonescaling(-(allops[opindex]->get_deltaQuantum(0)), b.get_braStateInfo().quanta[lQ], b.get_ketStateInfo().quanta[rQ]);
          int nrows = operator_element(lQ, rQ).Nrows();
          int ncols = operator_element(lQ, rQ).Ncols();
          for (int row=0; row < nrows; ++row)
            DAXPY(ncols, scaling*scaleDD[opindex], &(allops[opindex]->operator_element(rQ,lQ)(1, row+1)), nrows, &(this->operator_element(lQ, rQ)(row+1, 1)), 1);
        }
      } else { // operator DD
        if (allowed(lQ, rQ) && allops[opindex]->allowed(lQ,rQ) && fabs(scaleDD[opindex])> TINY) {
          MatrixScaleAdd(scaleDD[opindex], allops[opindex]->operator_element(lQ,rQ), this->operator_element(lQ, rQ));
        }
      }

      if (dmrginp.hamiltonian() == BCS && allowed(lQ, rQ) 
          && allops[opindex]->allowed(lQ,rQ) && fabs(scaleCC[opindex]) > TINY) {
        MatrixScaleAdd(scaleCC[opindex], allops[opindex]->operator_element(lQ,rQ), this->operator_element(lQ, rQ));
      }
    } else { // CkDl and ClDk
      const int k = allops[opindex]->get_orbs()[0];    
      const int l = allops[opindex]->get_orbs()[1];
      if (allowed(lQ,rQ) && allops[opindex]->allowed(lQ,rQ) && fabs(scaleCD[2*opindex]) > TINY)
        MatrixScaleAdd(scaleCD[2*opindex], allops[opindex]->operator_element(lQ,rQ), this->operator_element(lQ, rQ));

      if (!b.has(DES) && k != l && allowed(lQ,rQ) && allops[opindex]->allowed(rQ,lQ) && fabs(scaleCD[2*opindex+1]) > TINY) {
        double scaling = getStandAlonescaling(-(allops[opindex]->get_deltaQuantum(0)), b.get_braStateInfo().quanta[lQ], b.get_ketStateInfo().quanta[rQ]);
        int nrows = operator_element(lQ, rQ).Nrows();
        int ncols = operator_element(lQ, rQ).Ncols();
        for (int row=0; row < nrows; ++row)
          DAXPY(ncols, scaling*scaleCD[2*opindex+1], &(allops[opindex]->operator_element(rQ,lQ)(1, row+1)), nrows, &(this->operator_element(lQ, rQ)(row+1, 1)), 1);
      }
    }
  }

  //accumulateMultiThread(this, op_array, numthrds);

}

void SpinAdapted::StackDesDesComp::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];


  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_DESCOMP).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
      }
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(DES_DESCOMP).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
    } 

  const TwoElectronArray& v2 = *(b.get_twoInt()); 
  // explicitly build DD_comp
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
      {
	int k = leftBlock->get_sites()[kx];
	int l = rightBlock->get_sites()[lx];

      int spink = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[k]),\
	spinl = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[l]);
      if (dmrginp.spinAdapted()) {
	spink = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[k]]), \
	spinl = Symmetry::negativeofAbelian(dmrginp.spin_orbs_symmetry()[dmrginp.spatial_to_spin()[l]]);
      }

      if (Symmetry::negativeofAbelian(sym.getirrep()) == Symmetry::addAbelian(spink, spinl)) {
	//if ((-sym).getirrep() == Symmetry::addAbelian(spink, spinl)) {
	double scaleV = calcCompfactor(i,j,k,l, spin, DD,*(b.get_twoInt()), b.get_integralIndex());

	  double scaleV2 = calcCompfactor(i,j,l,k, spin, DD,*(b.get_twoInt()), b.get_integralIndex());
        
	  if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	    if (leftBlock->has(DES)) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	    
	      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      scaleV += parity*scaleV2;
	    
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
		SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
	    }
	    else {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      //StackTransposeview top1 = StackTransposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	      //StackTransposeview top2 = StackTransposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	      //op1->set_conjugacy('t'); op2->set_conjugacy('t');
	      double parity = getCommuteParity(-op1->get_deltaQuantum()[0], -op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      scaleV += parity*scaleV2;
	    
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
		SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, Transpose(*op1), Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
	      //op1->set_conjugacy('n'); op2->set_conjugacy('n');
	    }
	  }
	}

	if (dmrginp.hamiltonian() == BCS) {
	  TensorOp C(i,1), C2(j,1);
	  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);

	  TensorOp CK(k, 1), CL(l, 1);
	  TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
	  if (!CC2.empty) {
	    double scaleV = calcCompfactor(CC1, CC2, DD, v_cccc[b.get_integralIndex()]);
	    CK = TensorOp(k, 1);
	    CL = TensorOp(l, 1);
	    TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
	    double scaleV2 = calcCompfactor(CC1, CC2_commute, DD, v_cccc[b.get_integralIndex()]);

	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(2));
	      scaleV += parity*scaleV2;

	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
		SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
	      }
	    }
	  }
	  // Ck*Dl
	  CK = TensorOp(k, 1);
	  TensorOp DL = TensorOp(l, -1);
	  TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
	  if (!CD2.empty) {
	    double scaleV = calcCompfactor(CC1, CD2, DD, v_cccd[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      if (rightBlock->has(DES)) {
		boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
		SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
	      } else {
		//StackTransposeview top2 = StackTransposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
		boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
		//op2->set_conjugacy('t');
		SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, scaleV);
		//op2->set_conjugacy('n');
	      }
	    }
	  }
	  // Cl*Dk
	  CL = TensorOp(l,1);
	  TensorOp DK = TensorOp(k,-1);
	  CD2 = CL.product(DK, spin, sym.getirrep());
	  if (!CD2.empty) {
	    double scaleV = calcCompfactor(CC1, CD2, DD, v_cccd[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      if (leftBlock->has(DES)) {
		boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
		double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(1));
		SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, m, row, col, scaleV*parity);
	      } else {
		boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
		//op2->set_conjugacy('t');
		//StackTransposeview top2 = StackTransposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
		double parity = getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), get_deltaQuantum(1));
		SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, m, row, col, scaleV*parity);
		//op2->set_conjugacy('n');
	      }
	    }
	  }
	}
      }
  dmrginp.makeopsT -> stop();

}


void SpinAdapted::StackDesDesComp::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  TensorOp C(i,1), C2(j,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(DES_DESCOMP).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(DES_DESCOMP).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(DES_DESCOMP, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  
  // explicitly build DD_comp
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
      {
	int k = leftBlock->get_sites()[kx];
	int l = rightBlock->get_sites()[lx];

	TensorOp DK(k,-1), DL(l,-1);
	TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
	if (!DD2.empty) {
	  double scaleV = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());

	  DK=TensorOp(k,-1); DL=TensorOp(l,-1);
	  DD2 = DL.product(DK, spin, sym.getirrep(), k==l);
	  double scaleV2 = calcCompfactor(CC1, DD2, DD, *(b.get_twoInt()), b.get_integralIndex());
        
	  if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	    if (leftBlock->has(DES)) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	    
	      double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      scaleV += parity*scaleV2;
	    
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
		SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	    }
	    else {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      //StackTransposeview top1 = StackTransposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
	      //StackTransposeview top2 = StackTransposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
	      //op1->set_conjugacy('t'); op2->set_conjugacy('t');
	      double parity = getCommuteParity(-op1->get_deltaQuantum()[0], -op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	      scaleV += parity*scaleV2;
	    
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
		SpinAdapted::operatorfunctions::TensorProduct(leftBlock, Transpose(*op1), Transpose(*op2), &b, &(b.get_stateInfo()), *this, scaleV);
	      //op1->set_conjugacy('n'); op2->set_conjugacy('n');
	    }
	  }
	}

	if (dmrginp.hamiltonian() == BCS) {
	  TensorOp CK(k, 1), CL(l, 1);
	  TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
	  if (!CC2.empty) {
	    double scaleV = calcCompfactor(CC1, CC2, DD, v_cccc[b.get_integralIndex()]);
	    CK = TensorOp(k, 1);
	    CL = TensorOp(l, 1);
	    TensorOp CC2_commute = CL.product(CK, spin, sym.getirrep(), k==l);
	    double scaleV2 = calcCompfactor(CC1, CC2_commute, DD, v_cccc[b.get_integralIndex()]);

	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(2));
	      scaleV += parity*scaleV2;

	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
		SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	      }
	    }
	  }
	  // Ck*Dl
	  CK = TensorOp(k, 1);
	  DL = TensorOp(l, -1);
	  TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
	  if (!CD2.empty) {
	    double scaleV = calcCompfactor(CC1, CD2, DD, v_cccd[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      if (rightBlock->has(DES)) {
		boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
		SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	      } else {
		//StackTransposeview top2 = StackTransposeview(rightBlock->get_op_rep(CRE, getSpinQuantum(l), l));
		boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
		//op2->set_conjugacy('t');
		SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, scaleV);
		//op2->set_conjugacy('n');
	      }
	    }
	  }
	  // Cl*Dk
	  CL = TensorOp(l,1);
	  DK = TensorOp(k,-1);
	  CD2 = CL.product(DK, spin, sym.getirrep());
	  if (!CD2.empty) {
	    double scaleV = calcCompfactor(CC1, CD2, DD, v_cccd[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      if (leftBlock->has(DES)) {
		boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
		double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(1));
		SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	      } else {
		boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
		//op2->set_conjugacy('t');
		//StackTransposeview top2 = StackTransposeview(leftBlock->get_op_rep(CRE, getSpinQuantum(k), k));
		double parity = getCommuteParity(op1->get_deltaQuantum(0), -op2->get_deltaQuantum(0), get_deltaQuantum(1));
		SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, Transpose(*op2), &b, &(b.get_stateInfo()), *this, scaleV*parity);
		//op2->set_conjugacy('n');
	      }
	    }
	  }
	}
      }
  dmrginp.makeopsT -> stop();

}

double SpinAdapted::StackDesDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp C(I,1), C2(J,1);
  TensorOp CC1 = C.product(C2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), I==J);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  if (dmrginp.hamiltonian() != BCS && c1.n_is() != ladder[0].n_is()-2) return 0.0;
  if (dmrginp.hamiltonian() == BCS && abs(dn) > 2) return 0.0;
  if (dmrginp.spinAdapted() && (c1.S_is().getirrep() > ladder[0].S_is().getirrep()+spin || c1.S_is().getirrep() <ladder[0].S_is().getirrep()-spin)) return 0.0;

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {
      if (!dmrginp.spinAdapted() && c1.S_is().getirrep() != ladder[i].S_is().getirrep()+spin) continue;
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];	

            bool isZero = true;
            if (dmrginp.hamiltonian() == BCS) {
              if (ladder[i].det_rep.size() > 1) {
                pout << "StackDesDesComp::redMatrixElement failed" << endl;
                abort();
              }
              if (dn == 0) {
                for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
                isZero = true;
                for (auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[l] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
              } else if (dn == 2) {
                if (k == l) continue;
                for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 1 && s1.get_orbstring().get_occ_rep()[l] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
                isZero = true;
                for (auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 0 && s1.get_orbstring().get_occ_rep()[l] == 0) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
              } else if (dn == -2) {
                if (k == l) continue;
                for (auto it1 = c1.det_rep.begin(); it1 != c1.det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 0 && s1.get_orbstring().get_occ_rep()[l] == 0) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
                isZero = true;
                for (auto it1 = ladder[i].det_rep.begin(); it1 != ladder[i].det_rep.end(); ++it1) {
                  const Slater &s1 = it1->first;
                  if (s1.get_orbstring().get_occ_rep()[k] == 1 && s1.get_orbstring().get_occ_rep()[l] == 1) {
                    isZero = false;
                    break;
                  }
                }
                if (isZero) continue;
              }
            } else if (dmrginp.spinAdapted()) {
	      for (auto it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
		const Slater &s1 = it1->first;
		if ( (s1.get_orbstring().get_occ_rep()[2*k] == 0 || s1.get_orbstring().get_occ_rep()[2*k+1] == 0) && 
		     (s1.get_orbstring().get_occ_rep()[2*l] == 0 || s1.get_orbstring().get_occ_rep()[2*l+1] == 0 )) {
		  isZero = false;
		  break;
		}
	      }
	      if (isZero) continue;
	      isZero = true;
	      for (auto it1 = ladder[i].det_rep.begin(); it1!= ladder[i].det_rep.end(); it1++) {
		const Slater &s1 = it1->first;
		if ( (s1.get_orbstring().get_occ_rep()[2*k] == 1 || s1.get_orbstring().get_occ_rep()[2*k+1] == 1) && 
		     (s1.get_orbstring().get_occ_rep()[2*l] == 1 || s1.get_orbstring().get_occ_rep()[2*l+1] == 1 ) ){
		  isZero = false;
		  break;
		}
	      }
	      if (isZero) continue;
	    } else {
              if (ladder[i].det_rep.size() > 1) {
                pout << "StackDesDesComp::redMatrixElement failed" << endl;
                abort();
              }
              if (k == l) continue;
	      bool isZero = true;
	      for (auto it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
		const Slater &s1 = it1->first;
		if (s1.get_orbstring().get_occ_rep()[k] == 0 && s1.get_orbstring().get_occ_rep()[l] == 0) {
		  isZero = false;
		  break;
		}
	      }
	      if (isZero) continue;
	      isZero = true;
	      for (auto it1 = ladder[i].det_rep.begin(); it1!= ladder[i].det_rep.end(); it1++) {
		const Slater &s1 = it1->first;
		if (s1.get_orbstring().get_occ_rep()[k] == 1 && s1.get_orbstring().get_occ_rep()[l] == 1){
		  isZero = false;
		  break;
		}
	      }
	      if (isZero) continue;
            }
           
            if (dmrginp.hamiltonian() == BCS && dn == 2) {
              TensorOp CK(k,1), CL(l,1);
              TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
              if (!CC2.empty) {
                double MatElements = calcMatrixElements(c1, CC2, ladder[i], backupSlater1, backupSlater2, index) ;
                double scale = calcCompfactor(CC1, CC2, DD, v_cccc[b->get_integralIndex()]);
                element += MatElements*scale/cleb;                
              }
            } else if (dmrginp.hamiltonian() == BCS && dn == 0) {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                double MatElements = calcMatrixElements(c1, CD2, ladder[i], backupSlater1, backupSlater2, index) ;
                double scale = calcCompfactor(CC1, CD2, DD, v_cccd[b->get_integralIndex()]);
                element += MatElements*scale/cleb;
              }
            } else {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);

              if (!DD2.empty) {
		double MatElements = calcMatrixElements(c1, DD2, ladder[i], backupSlater1, backupSlater2, index) ;
                double scale = calcCompfactor(CC1, DD2, DD, index, *(b->get_twoInt()), b->get_integralIndex());
                element += MatElements*scale/cleb;
              }
            } 
          }
        break;
      }
      else
        continue;

    }
  }
  return element;

}


//******************CRECRECOMP*****************

void SpinAdapted::StackCreCreComp::buildfromCreCre(StackSpinBlock& b) 
{
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();
 
  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  //StackCreCreComp* op_array; 
  //initiateMultiThread(this, op_array, numthrds);
 
  TensorOp D(i, -1), D2(j, -1);
  TensorOp DD1 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);
 
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;

  for (int ii=0; ii<b.get_op_array(CRE_CRE).get_size(); ii++)
    for (int ji=0; ji<b.get_op_array(CRE_CRE).get_local_element(ii).size(); ji++) 
      if (b.get_op_array(CRE_CRE).get_local_element(ii)[ji]->get_deltaQuantum(0) == deltaQuantum[0])
	allops.push_back(b.get_op_array(CRE_CRE).get_local_element(ii)[ji]);

  //#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii<allops.size(); ii++) {
    const int k = allops[ii]->get_orbs()[0];
    const int l = allops[ii]->get_orbs()[1];
   
    TensorOp CK(k, 1), CL(l, 1);      
    TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
    if (!CC2.empty) {
      double scaleV = calcCompfactor(CC2, DD1, DD,*(b.get_twoInt()), b.get_integralIndex());      
     
      CK=TensorOp(k, 1); CL=TensorOp(l, 1);
      CC2 = CL.product(CK, spin, sym.getirrep(), k==l);
      double scaleV2 = calcCompfactor(CC2, DD1, DD, *(b.get_twoInt()), b.get_integralIndex());
     
      double parity = getCommuteParity(getSpinQuantum(k), getSpinQuantum(l), get_deltaQuantum()[0]);
      if (k != l)
	scaleV += parity*scaleV2;

      ScaleAdd(scaleV, *allops[ii], *this);
     
    }
  }

  //accumulateMultiThread(this, op_array, numthrds);

}


void SpinAdapted::StackCreCreComp::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());
  int spin = deltaQuantum[0].get_s().getirrep();
  IrrepSpace sym = deltaQuantum[0].get_symm();

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  memset(data, 0, totalMemory * sizeof(double));

  TensorOp D(i, -1), D2(j, -1);
  TensorOp DD1 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), i==j);

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRECOMP).has(i, j))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRECOMP, deltaQuantum, i,j);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRECOMP).has(i, j))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_CRECOMP, deltaQuantum, i,j);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  // explicitly build CC_comp
  for (int kx = 0; kx < leftBlock->get_sites().size(); ++kx)
    for (int lx = 0; lx < rightBlock->get_sites().size(); ++lx)
      {
	int k = leftBlock->get_sites()[kx];
	int l = rightBlock->get_sites()[lx];

	TensorOp CK(k, 1), CL(l, 1);
	TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
	if (!CC2.empty) {
	  double scaleV = calcCompfactor(CC2, DD1, DD, *(b.get_twoInt()), b.get_integralIndex());

	  CK=TensorOp(k, 1); CL=TensorOp(l, 1);
	  CC2 = CL.product(CK, spin, sym.getirrep(), k==l);
	  double scaleV2 = calcCompfactor(CC2, DD1, DD, *(b.get_twoInt()), b.get_integralIndex());
        
	  if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	    boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	  
	    double parity = getCommuteParity(op1->get_deltaQuantum()[0], op2->get_deltaQuantum()[0], get_deltaQuantum()[0]);
	    scaleV += parity*scaleV2;
	  
	    if (fabs(scaleV) > dmrginp.twoindex_screen_tol())
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	  }
	}

	if (dmrginp.hamiltonian() == BCS) {
	  TensorOp DK(k,-1), DL(l,-1);
	  TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);
	  if (!DD2.empty) {
	    double scaleV = calcCompfactor(DD1, DD2, DD, v_cccc[b.get_integralIndex()]);
	    DK = TensorOp(k,-1);
	    DL = TensorOp(l,-1);
	    TensorOp DD2_commute = DL.product(DK, spin, sym.getirrep(), k==l);
	    double scaleV2 = calcCompfactor(DD1, DD2_commute, DD, v_cccc[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && (fabs(scaleV2)+fabs(scaleV)) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	      double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(2));
	      scaleV += parity*scaleV2;
	      if (fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
		SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	      }
	    }
	  }
	  // Cl*Dk
	  CL = TensorOp(l,1);
	  DK = TensorOp(k,-1);
	  TensorOp CD2 = CL.product(DK, spin, sym.getirrep());
	  if (!CD2.empty) {
	    double scaleV = calcCompfactor(DD1, CD2, DD, v_cccd[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = rightBlock->get_op_rep(CRE, getSpinQuantum(l), l);
	      boost::shared_ptr<StackSparseMatrix> op2 = leftBlock->get_op_rep(DES, -getSpinQuantum(k), k);
	      double parity = getCommuteParity(op1->get_deltaQuantum(0), op2->get_deltaQuantum(0), get_deltaQuantum(1));
	      SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV*parity);
	    }
	  }
	  // Ck*Dl
	  CK = TensorOp(k, 1);
	  DL = TensorOp(l, -1);
	  CD2 = CK.product(DL, spin, sym.getirrep());
	  if (!CD2.empty) {
	    double scaleV = calcCompfactor(DD1, CD2, DD, v_cccd[b.get_integralIndex()]);
	    if (leftBlock->get_op_array(CRE).has(k) && rightBlock->get_op_array(CRE).has(l) && fabs(scaleV) > dmrginp.twoindex_screen_tol()) {
	      boost::shared_ptr<StackSparseMatrix> op1 = leftBlock->get_op_rep(CRE, getSpinQuantum(k), k);
	      boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(DES, -getSpinQuantum(l), l);
	      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op1, *op2, &b, &(b.get_stateInfo()), *this, scaleV);
	    }
	  }
	}
      }
  dmrginp.makeopsT -> stop();

}

double SpinAdapted::StackCreCreComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int I = get_orbs()[0], 
    J = get_orbs()[1]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();

  TensorOp D(I,-1), D2(J,-1);
  TensorOp DD1 = D.product(D2, (-deltaQuantum[0].get_s()).getirrep(), (-sym).getirrep(), I==J);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  for (int j = 0; j < deltaQuantum.size(); ++j) {
    for (int i=0; i<ladder.size(); i++) {

      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
        for (int kl =0; kl<b->get_sites().size(); kl++) 
          for (int kk =0; kk<b->get_sites().size(); kk++) {

            int k = b->get_sites()[kk];
            int l = b->get_sites()[kl];	
           
            if (dmrginp.hamiltonian() == BCS && dn == -2) {
              TensorOp DK(k,-1), DL(l,-1);
              TensorOp DD2 = DK.product(DL, spin, sym.getirrep(), k==l);

              if (!DD2.empty) {
                double MatElements = calcMatrixElements(c1, DD2, ladder[i], backupSlater1, backupSlater2, index) ;
                double scale = calcCompfactor(DD1, DD2, DD, v_cccc[b->get_integralIndex()]);
                element += MatElements*scale/cleb;                
              }
            } 
	    else if (dmrginp.hamiltonian() == BCS && dn == 0) {
              TensorOp CK(k,1), DL(l,-1);
              TensorOp CD2 = CK.product(DL, spin, sym.getirrep());
              if (!CD2.empty) {
                double MatElements = calcMatrixElements(c1, CD2, ladder[i], backupSlater1, backupSlater2, index) ;
                double scale = calcCompfactor(DD1, CD2, DD, v_cccd[b->get_integralIndex()]);
                element += MatElements*scale/cleb;                
              }
            } 
	    else {
              TensorOp CK(k,1), CL(l,1);
	      TensorOp CC2 = CK.product(CL, spin, sym.getirrep(), k==l);
	      
	      if (!CC2.empty) {
		double MatElements = calcMatrixElements(c1, CC2, ladder[i], backupSlater1, backupSlater2, index) ;
		double scale = calcCompfactor(CC2, DD1, DD, index, *(b->get_twoInt()), b->get_integralIndex());
		element += MatElements*scale/cleb;
	      }
            }
          }
        break;
      }
      else
        continue;

    }
  }
  return element;

}


//******************CREDESDESCOMP*****************

void SpinAdapted::StackCreDesDesComp::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  const int k = get_orbs()[0];

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());

  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  StackSpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);

  if (leftBlock->get_op_array(CRE_DES_DESCOMP).has(k))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DESCOMP, deltaQuantum, k);
      if (rightBlock->get_sites().size() == 0) 
	SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
      else {
	//const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
	SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
	const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
	SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
      }
    }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_DES_DESCOMP).has(k))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_DES_DESCOMP, deltaQuantum, k);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  

  // explicit build CDD_comp
  if (dmrginp.hamiltonian() != HUBBARD){
    if (loopBlock->has(CRE_DESCOMP))
      {
	FUNCTOR f = boost::bind(&stackopxop::dxcdcomp, otherBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(loopBlock->get_op_array(DES), f);
	
	f = boost::bind(&stackopxop::cxddcomp, otherBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(loopBlock->get_op_array(CRE), f);
        
	f = boost::bind(&stackopxop::dxcdcomp, loopBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(otherBlock->get_op_array(DES), f);

	f = boost::bind(&stackopxop::cxddcomp, loopBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(otherBlock->get_op_array(CRE), f);
      }
    else if (otherBlock->has(CRE_DESCOMP))
      {
	pout << "I should not be here"<<endl;exit(0);
      }
  }
  dmrginp.makeopsT -> stop();


}

double SpinAdapted::StackCreDesDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int K = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  TensorOp CK(K, 1);
  for (int j = 0; j<deltaQuantum.size(); ++j)
    for (int i=0; i<ladder.size(); ++i) {
      int index = 0; double cleb=0.0;
      if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	for (int ki =0; ki<b->get_sites().size(); ki++) 
	  for (int kj =0; kj<b->get_sites().size(); kj++) 
	    for (int kl =0; kl<b->get_sites().size(); kl++) {
	      int _i = b->get_sites()[ki];
	      int _j = b->get_sites()[kj];
	      int _l = b->get_sites()[kl];
	      SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
	      if (dmrginp.hamiltonian() == BCS && dn == -3) {
		std::vector<SpinQuantum> sij = (-si)-sj;
		for (int ij = 0; ij < sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij-sl;
		  for (int ijl=0;ijl<sijl.size();++ijl){
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp DI(_i, -1), DJ(_j, -1), DL(_l, -1);
		    TensorOp DDIJ = DI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		    TensorOp DDDIJL = DDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (DDDIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, DDDIJL, ladder[i], backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(DDDIJL, CK, CDD, v_cccd[b->get_integralIndex()]);
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }
		}
	      } 
	      else if (dmrginp.hamiltonian() == BCS && dn == 1) {
		std::vector<SpinQuantum> sij = si+sj;
		for (int ij=0; ij<sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij-sl;
		  for (int ijl=0; ijl<sijl.size(); ++ijl) {
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
		    TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		    TensorOp CCDIJL = CCIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (CCDIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, CCDIJL, ladder[i], backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CCDIJL, CK, CDD, v_cccd[b->get_integralIndex()]);
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }            
		}
	      } 
	      else if (dmrginp.hamiltonian() == BCS && dn == 3) {
		std::vector<SpinQuantum> sij = si+sj;
		for (int ij =0; ij<sij.size();++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij+sl;
		  for (int ijl=0; ijl<sijl.size(); ++ijl) {            
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp CI(_i, 1), CJ(_j, 1), CL(_l, 1);
		    TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		    TensorOp CCCIJL = CCIJ.product(CL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (CCCIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, CCCIJL, ladder[i], backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CCCIJL, CK, CDD, v_cccc[b->get_integralIndex()]);
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }
		}
	      } 
	      else { //CDD
		std::vector<SpinQuantum> sij = si-sj;
		for (int ij=0; ij<sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij-sl;
		  for (int ijl=0; ijl<sijl.size(); ijl++) {
		    SpinQuantum symijl = sijl[ijl];
		    if (symijl != deltaQuantum[j]) continue;
	          
		    TensorOp CI(_i, 1), DJ(_j, -1), DL(_l, -1);
	          
		    TensorOp CDIJ = CI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep());

		    TensorOp CDDIJL = CDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (CDDIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, CDDIJL, ladder[i], backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CDDIJL, CK, CDD, *(b->get_twoInt()), b->get_integralIndex());
		    if (dmrginp.spinAdapted()) scale*=-1; //terrible hack
		    if (fabs(scale) > dmrginp.oneindex_screen_tol()) 
		      element += MatElements*scale/cleb;
		  }
		}
	      }
	    }
	for (int ki =0; ki<b->get_sites().size(); ki++) {
	  int _i = b->get_sites()[ki];
	  if (dmrginp.hamiltonian() == BCS && dn == 1) {
	    TensorOp CI(_i, 1);
	    double MatElements = calcMatrixElements(c1, CI, ladder[i], backupSlater1, backupSlater2, index) ;
	    double factor = calcCompfactor(CI, CK, C, *(b->get_twoInt()), b->get_integralIndex());
	    if (fabs(factor) > dmrginp.oneindex_screen_tol())
	      element += factor*MatElements/cleb;
	  } 
	  else {
	    TensorOp DI(_i, -1);
	    double MatElements = calcMatrixElements(c1, DI, ladder[i], backupSlater1, backupSlater2, index) ;
	    double factor = calcCompfactor(CK, DI, C, *(b->get_twoInt()), b->get_integralIndex());
	    if (fabs(factor) > dmrginp.oneindex_screen_tol())
	      element += factor*MatElements/cleb;
	  }
	}
	break;
      }
      else
	continue;
    }
  return element;
}


//******************CRECREDESCOMP*****************


void SpinAdapted::StackCreCreDesComp::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }


  const int k = get_orbs()[0];


  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  StackSpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);

  if (leftBlock->get_op_array(CRE_CRE_DESCOMP).has(k)) {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    else {
      //const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
    }
  }
  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE_DESCOMP).has(k))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
    }  
  // explicit build CCD_comp
  if (dmrginp.hamiltonian() != HUBBARD){
    if (loopBlock->has(CRE_DESCOMP))
      {

	FUNCTOR f = boost::bind(&stackopxop::cxcdcompElement, otherBlock, _1, &b, k, this, m, row, col, 1.0); 
	for_all_singlethread(loopBlock->get_op_array(CRE), f);

	f = boost::bind(&stackopxop::dxcccompElement, otherBlock, _1, &b, k, this, m, row, col, 2.0); // factor of 2.0 because CCcomp_{ij} = -CCcomp_{ji}
	for_all_singlethread(loopBlock->get_op_array(CRE), f);

	f = boost::bind(&stackopxop::cxcdcompElement, loopBlock, _1, &b, k, this, m, row, col, 1.0); 
	for_all_singlethread(otherBlock->get_op_array(CRE), f);

	f = boost::bind(&stackopxop::dxcccompElement, loopBlock, _1, &b, k, this, m, row, col, 2.0);
	for_all_singlethread(otherBlock->get_op_array(CRE), f);

      } else if (otherBlock->has(CRE_DESCOMP)) {
      pout << "I should not be here"<<endl;exit(0);
    } 
  }

  dmrginp.makeopsT -> stop();
}

void SpinAdapted::StackCreCreDesComp::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  const int k = get_orbs()[0];

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());

  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  StackSpinBlock* loopBlock, *otherBlock;
  assignloopblock(loopBlock, otherBlock, leftBlock, rightBlock);

  if (leftBlock->get_op_array(CRE_CRE_DESCOMP).has(k)) {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    else {
      //const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
    }
  }

  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    return;
  }
  if (rightBlock->get_op_array(CRE_CRE_DESCOMP).has(k))
    {
      const boost::shared_ptr<StackSparseMatrix> op = rightBlock->get_op_rep(CRE_CRE_DESCOMP, deltaQuantum, k);
      //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
      SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
      const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
      SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);
    }  

  // explicit build CCD_comp
  if (dmrginp.hamiltonian() != HUBBARD){
    if (loopBlock->has(CRE_DESCOMP))
      {
	FUNCTOR f = boost::bind(&stackopxop::cxcdcomp, otherBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(loopBlock->get_op_array(CRE), f);

	f = boost::bind(&stackopxop::dxcccomp, otherBlock, _1, &b, k, this, 2.0); // factor of 2.0 because CCcomp_{ij} = -CCcomp_{ji}
	for_all_singlethread(loopBlock->get_op_array(CRE), f);

	f = boost::bind(&stackopxop::cxcdcomp, loopBlock, _1, &b, k, this, 1.0); 
	for_all_singlethread(otherBlock->get_op_array(CRE), f);

	f = boost::bind(&stackopxop::dxcccomp, loopBlock, _1, &b, k, this, 2.0);
	for_all_singlethread(otherBlock->get_op_array(CRE), f);

      } else if (otherBlock->has(CRE_DESCOMP)) {
      pout << "I should not be here"<<endl;exit(0);
    } 
  }

  dmrginp.makeopsT -> stop();
}


double SpinAdapted::StackCreCreDesComp::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  int K = get_orbs()[0]; //convert spatial id to spin id because slaters need that
  IrrepSpace sym = deltaQuantum[0].get_symm();
  int spin = deltaQuantum[0].get_s().getirrep();
  bool finish = false;
  int dn = c1.n_is() - ladder[0].n_is();  

  TensorOp D(K, -1);
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  if (dmrginp.hamiltonian() != BCS && c1.n_is() != ladder[0].n_is()+1) return 0.0;
  if (dmrginp.hamiltonian() == BCS && abs(c1.n_is()-ladder[0].n_is()) > 3) return 0.0;
  if (dmrginp.spinAdapted() && (c1.S_is().getirrep() > ladder[0].S_is().getirrep()+1 || c1.S_is().getirrep() <ladder[0].S_is().getirrep() -1)) return 0.0;

  if (dmrginp.hamiltonian() != BCS) {
    int ladidx;
    if (!dmrginp.spinAdapted()) {
      ladidx = -1;
      for (int i = 0; i < ladder.size(); ++i)
        if (ladder[i].S_is().getirrep() + spin == c1.S_is().getirrep()) {
          ladidx = i; break;
        }
      if (ladidx < 0) return 0.;
    } else {
      ladidx = 0;
    }

    Csf& detladder = ladder[ladidx];
    if (!dmrginp.spinAdapted() && detladder.det_rep.size() > 1) {
      pout << detladder << endl;
      pout << "CreCreDesComp::redMatrixElement failed" << endl;
      abort();
    }
    for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++)
	for (int kl =0; kl<b->get_sites().size(); kl++) {
	  int _i = b->get_sites()[ki];
	  int _j = b->get_sites()[kj];
	  int _l = b->get_sites()[kl];

	  if (!dmrginp.spinAdapted()) {
	    for (auto it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
	      const Slater &s1 = it1->first;
	      if (s1.get_orbstring().get_occ_rep()[_i] == 1 && s1.get_orbstring().get_occ_rep()[_j] == 1) {
		if (_l != _i && _l != _j && s1.get_orbstring().get_occ_rep()[_l] == 0) {
		  for (auto it2 = detladder.det_rep.begin(); it2 != detladder.det_rep.end(); ++it2) {
		    const Slater &s2 = it2 -> first;
		    if (s2.get_orbstring().get_occ_rep()[_i] == 0 && s2.get_orbstring().get_occ_rep()[_j] == 0 
			&& s2.get_orbstring().get_occ_rep()[_l] == 1) {
		      goto dontstop;
		    }
		  }
		} else {
		  for (auto it2 = detladder.det_rep.begin(); it2 != detladder.det_rep.end(); ++it2) {
		    const Slater &s2 = it2->first;
		    if (s2.get_orbstring().get_occ_rep()[_l] == 1) {
		      goto dontstop;
		    }
		  }
		}
	      }
	    }
	  } else {
	    for (map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
	      const Slater &s1 = it1->first;
	      for (int ixx = dmrginp.spatial_to_spin(_i); ixx <dmrginp.spatial_to_spin(_i+1); ixx++)
		if ( (s1.get_orbstring().get_occ_rep()[ixx] == 1) ) {
		  for (int jxx = dmrginp.spatial_to_spin(_j); jxx <dmrginp.spatial_to_spin(_j+1); jxx++)
		    if ( (s1.get_orbstring().get_occ_rep()[jxx] == 1) ) {
		      if (_l != _i && _l != _j) {
			for (int lxx1 = dmrginp.spatial_to_spin(_l); lxx1 <dmrginp.spatial_to_spin(_l+1); lxx1++)
			  if ( s1.get_orbstring().get_occ_rep()[lxx1] == 0) 
			    for (map<Slater, double>::iterator it2 = detladder.det_rep.begin(); it2!= detladder.det_rep.end(); it2++) {
			      const Slater &s2 = it2->first;
			      for (int ixx2 = dmrginp.spatial_to_spin(_i); ixx2 <dmrginp.spatial_to_spin(_i+1); ixx2++)
				if ( (s2.get_orbstring().get_occ_rep()[ixx2] == 0) )
				  for (int jxx2 = dmrginp.spatial_to_spin(_j); jxx2 <dmrginp.spatial_to_spin(_j+1); jxx2++)
				    if ( (s2.get_orbstring().get_occ_rep()[jxx2] == 0) )
				      for (int lxx = dmrginp.spatial_to_spin(_l); lxx <dmrginp.spatial_to_spin(_l+1); lxx++)
					if ( s2.get_orbstring().get_occ_rep()[lxx] == 1) { goto dontstop; }
			    }
		      } else {      
			for (map<Slater, double>::iterator it2 = detladder.det_rep.begin(); it2!= detladder.det_rep.end(); it2++) {
			  const Slater &s2 = it2->first;
			  for (int lxx = dmrginp.spatial_to_spin(_l); lxx <dmrginp.spatial_to_spin(_l+1); lxx++)
			    if ( s2.get_orbstring().get_occ_rep()[lxx] == 1) 
			      goto dontstop;
			}
		      }
		    }
		}
	    }
	  }
	  continue;    
    
	dontstop:
      
	  for (int j = 0; j < deltaQuantum.size(); ++j) {
	    SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
	    std::vector<SpinQuantum> sij = si+sj;
	    for (int ij=0; ij<sij.size(); ++ij) {
	      SpinQuantum symij = sij[ij];
	      std::vector<SpinQuantum> sijl = symij-sl;
	      for (int ijl=0; ijl<sijl.size(); ++ijl) {
		if (sijl[ijl] != deltaQuantum[j]) continue;
		SpinQuantum symijl = sijl[ijl];
		TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
		TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		TensorOp CCDIJL = CCIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		if (CCDIJL.empty) continue;

		for (int i=0; i<ladder.size(); ++i) {
		  int index = 0; double cleb=0.0;
		  if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
		    double MatElements = calcMatrixElements(c1, CCDIJL, ladder[i], backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CCDIJL, D, CCD, *(b->get_twoInt()), b->get_integralIndex());
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		    break;
		  }
		}
	      }
	    }
	  }
	}

    for (int ki =0; ki<b->get_sites().size(); ki++) {  // add C block to CCD block
      int _i = b->get_sites()[ki];
      for (int j = 0; j < deltaQuantum.size(); ++j)
	for (int i=0; i<ladder.size(); ++i) {
	  int index = 0; double cleb=0.0;
	  if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
	    TensorOp CI(_i, 1);
	    double MatElements = calcMatrixElements(c1, CI, ladder[i], backupSlater1, backupSlater2, index) ;
	    double factor = calcCompfactor(CI, D, C, *(b->get_twoInt()), b->get_integralIndex());
	    if (fabs(factor) > dmrginp.oneindex_screen_tol())
	      element += factor*MatElements/cleb;
	    break;
	  }
	} 
    }
  } else { // BCS
    if (dmrginp.spinAdapted()) abort();

    int ladidx = -1;
    for (int i = 0; i < ladder.size(); ++i)
      if (ladder[i].S_is().getirrep() + spin == c1.S_is().getirrep()) {
        ladidx = i; break;
      }

    if (ladidx < 0) return 0.;

    Csf& detladder = ladder[ladidx];

    if (detladder.det_rep.size() > 1) {
      pout << detladder << endl;
      pout << "CreCreDesComp redMatrixElement failed" << endl;
      abort();
    }

    for (int ki =0; ki<b->get_sites().size(); ki++) 
      for (int kj =0; kj<b->get_sites().size(); kj++)
	for (int kl =0; kl<b->get_sites().size(); kl++) {
	  int _i = b->get_sites()[ki];
	  int _j = b->get_sites()[kj];
	  int _l = b->get_sites()[kl];
	  // bypass zero elements
	  for (auto it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); ++it1) {
	    const Slater &s1 = it1->first;
	    if (_i == _j && _j == _l) {
	      continue;
	    } else if (_i != _j && _i != _l && _j != _l) {
	      bool occ1_i = s1.get_orbstring().get_occ_rep()[_i]; 
	      bool occ1_j = s1.get_orbstring().get_occ_rep()[_j];
	      bool occ1_l = s1.get_orbstring().get_occ_rep()[_l];
	      if (occ1_i < occ1_j || occ1_j < occ1_l) {
		continue;
	      }
	      for (auto it2 = detladder.det_rep.begin(); it2 != detladder.det_rep.end(); ++it2) {
		const Slater &s2 = it2->first;
		bool occ2_i = s2.get_orbstring().get_occ_rep()[_i];
		bool occ2_j = s2.get_orbstring().get_occ_rep()[_j];
		bool occ2_l = s2.get_orbstring().get_occ_rep()[_l];
		if (occ2_i != occ1_i && occ2_j != occ1_j && occ2_l != occ1_l) goto dontstop1;
	      }
	    } else if (_i == _j) {
	      bool occ1_i = s1.get_orbstring().get_occ_rep()[_i]; 
	      bool occ1_l = s1.get_orbstring().get_occ_rep()[_l];
	      if (!occ1_i || occ1_l) continue;         
	      for (auto it2 = detladder.det_rep.begin(); it2 != detladder.det_rep.end(); ++it2) {
		const Slater &s2 = it2->first;
		bool occ2_i = s2.get_orbstring().get_occ_rep()[_i];
		bool occ2_l = s2.get_orbstring().get_occ_rep()[_l];
		if (occ2_i && occ2_l) goto dontstop1;            
	      }
	    } else if (_i == _l) {
	      bool occ1_i = s1.get_orbstring().get_occ_rep()[_i]; 
	      bool occ1_j = s1.get_orbstring().get_occ_rep()[_j];
	      if (!occ1_i) continue;
	      for (auto it2 = detladder.det_rep.begin(); it2 != detladder.det_rep.end(); ++it2) {
		const Slater &s2 = it2->first;
		bool occ2_i = s2.get_orbstring().get_occ_rep()[_i];
		bool occ2_j = s2.get_orbstring().get_occ_rep()[_j];
		if (occ2_i && (occ1_j != occ2_j)) goto dontstop1;
	      }
	    } else if (_j == _l) {
	      bool occ1_i = s1.get_orbstring().get_occ_rep()[_i]; 
	      bool occ1_j = s1.get_orbstring().get_occ_rep()[_j];
	      if (!occ1_i || !occ1_j) continue;
	      for (auto it2 = detladder.det_rep.begin(); it2 != detladder.det_rep.end(); ++it2) {
		const Slater &s2 = it2->first;
		bool occ2_i = s2.get_orbstring().get_occ_rep()[_i];
		bool occ2_j = s2.get_orbstring().get_occ_rep()[_j];
		if (!occ2_i && occ2_j) goto dontstop1;
	      }
	    }
	  }

	  continue;

	dontstop1:
	  for (int j = 0; j < deltaQuantum.size(); ++j) {
	    int index = 0; double cleb=0.0;
	    if (nonZeroTensorComponent(c1, deltaQuantum[j], detladder, index, cleb)) {
	      if (dmrginp.hamiltonian() == BCS && dn == 3) { // CCC
		SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
		std::vector<SpinQuantum> sij = si+sj;
		for (int ij=0; ij<sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij+sl;
		  for (int ijl=0; ijl<sijl.size();++ijl) {
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp CI(_i, 1), CJ(_j, 1), CL(_l, 1);
		    TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		    // FIXME maybe there takes specical care when j = l in spinadapted case
		    TensorOp CCCIJL = CCIJ.product(CL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (CCCIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, CCCIJL, detladder, backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CCCIJL, D, CCD, v_cccd[b->get_integralIndex()]);
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }
		}
	      } else if (dmrginp.hamiltonian() == BCS && dn == -1) { // CDD
		SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
		std::vector<SpinQuantum> sij = si-sj;
		for (int ij=0; ij<sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij-sl;
		  for (int ijl=0; ijl<sijl.size(); ++ijl) {
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp CI(_i, 1), DJ(_j, -1), DL(_l, -1);
		    TensorOp CDIJ = CI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep());
		    TensorOp CDDIJL = CDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (CDDIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, CDDIJL, detladder, backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CDDIJL, D, CCD, v_cccd[b->get_integralIndex()]);
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }
		}
	      } else if (dmrginp.hamiltonian() == BCS && dn == -3) { // DDD
		SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
		std::vector<SpinQuantum> sij = (-si)-sj;
		for (int ij=0; ij<sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij-sl;
		  for (int ijl=0; ijl<sijl.size(); ++ijl) {            
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp DI(_i, -1), DJ(_j, -1), DL(_l, -1);
		    TensorOp DDIJ = DI.product(DJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		    TensorOp DDDIJL = DDIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (DDDIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, DDDIJL, detladder, backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(DDDIJL, D, CCD, v_cccc[b->get_integralIndex()]);
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }
		}
	      } else {  // CCD
		SpinQuantum si=getSpinQuantum(_i), sj=getSpinQuantum(_j), sl=getSpinQuantum(_l);
		std::vector<SpinQuantum> sij = si+sj;
		for (int ij=0; ij<sij.size(); ++ij) {
		  SpinQuantum symij = sij[ij];
		  std::vector<SpinQuantum> sijl = symij-sl;
		  for (int ijl=0; ijl<sijl.size(); ++ijl) {
		    if (sijl[ijl] != deltaQuantum[j]) continue;
		    SpinQuantum symijl = sijl[ijl];
		    TensorOp CI(_i, 1), CJ(_j, 1), DL(_l, -1);
		    TensorOp CCIJ = CI.product(CJ, symij.get_s().getirrep(), symij.get_symm().getirrep(), _i==_j);
		    TensorOp CCDIJL = CCIJ.product(DL, symijl.get_s().getirrep(), symijl.get_symm().getirrep());
		    if (CCDIJL.empty) continue;
		    double MatElements = calcMatrixElements(c1, CCDIJL, detladder, backupSlater1, backupSlater2, index) ;
		    double scale = calcCompfactor(CCDIJL, D, CCD, *(b->get_twoInt()), b->get_integralIndex());
		    if (fabs(scale) > dmrginp.oneindex_screen_tol())
		      element += MatElements*scale/cleb;
		  }
		}
	      }
	      break;
	    }
	  }
	}

    for (int ki =0; ki<b->get_sites().size(); ki++) {  // add C block to CCD block
      int _i = b->get_sites()[ki];
      for (int j = 0; j < deltaQuantum.size(); ++j) {
        int index = 0; double cleb=0.0;
        if (nonZeroTensorComponent(c1, deltaQuantum[j], detladder, index, cleb)) {
          if (dmrginp.hamiltonian() == BCS && dn == -1) { // D
            TensorOp DI(_i, -1);
            double MatElements = calcMatrixElements(c1, DI, detladder, backupSlater1, backupSlater2, index) ;
            double factor = calcCompfactor(DI, D, C, *(b->get_twoInt()), b->get_integralIndex());
            if (fabs(factor) > dmrginp.oneindex_screen_tol())
              element += factor*MatElements/cleb;
          } else { // C
            TensorOp CI(_i, 1);
	    double MatElements = calcMatrixElements(c1, CI, detladder, backupSlater1, backupSlater2, index) ;
            double factor = calcCompfactor(CI, D, C, *(b->get_twoInt()), b->get_integralIndex());
            if (fabs(factor) > dmrginp.oneindex_screen_tol())
              element += factor*MatElements/cleb;
          }
	  break;
        }
      } 
    }
  }
  return element;
}


//******************HAM*****************

//the memory must already have been allocated
void SpinAdapted::StackHam::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();


  StackSpinBlock* loopBlock=rightBlock, *otherBlock=leftBlock;


  StackHam *op_array = this;

  //initiateMultiThread(this, op_array, 1);

  boost::shared_ptr<StackSparseMatrix> op = leftBlock->get_op_rep(HAM, deltaQuantum);


  //this action is only performed on the 0th process
  if (rightBlock->get_sites().size() == 0) 
    SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
  else {
    //const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
  }
  

  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    //distributedaccumulate(*this);
    dmrginp.makeopsT -> stop();
    return;
  }

  op = rightBlock->get_op_rep(HAM, deltaQuantum);
  //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
  SpinAdapted::operatorfunctions::TensorProductElement(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);

  // CCD_A*D_B + CCD_B*D_A + c.c. 
  FUNCTOR f = boost::bind(&stackopxop::cxcddcomp_Element, leftBlock, _1, &b, op_array, boost::ref(m), row, col); 
  for_all_singlethread(rightBlock->get_op_array(CRE), f);

  f = boost::bind(&stackopxop::cxcddcomp_Element, rightBlock, _1, &b, op_array, boost::ref(m), row, col); 
  for_all_singlethread(leftBlock->get_op_array(CRE), f);  

  if (dmrginp.hamiltonian() != HUBBARD) {    
    f = boost::bind(&stackopxop::cdxcdcomp_Element, otherBlock, _1, &b, op_array, boost::ref(m), row, col );
    for_all_singlethread(loopBlock->get_op_array(CRE_DES), f);
    
    f = boost::bind(&stackopxop::ddxcccomp_Element, otherBlock, _1, &b, op_array, boost::ref(m), row, col);
    for_all_singlethread(loopBlock->get_op_array(CRE_CRE), f);
  }
  
}


//the memory must already have been allocated
void SpinAdapted::StackHam::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //cannot build
  dmrginp.makeopsT -> start();

  //check if we have enough memory
  //assert(totalMemory == getRequiredMemory(b, deltaQuantum));

  //check if the operatorMatrix has been initialized appropriately
  //assert(operatorMatrix.nrows() == b.get_braStateInfo().quanta.size() && operatorMatrix.ncols() == b.get_ketStateInfo().quanta.size());

  memset(data, 0, totalMemory * sizeof(double));

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();


  StackSpinBlock* loopBlock=rightBlock, *otherBlock=leftBlock;



#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  StackHam *op_array = this;

  //initiateMultiThread(this, op_array, 1);

  boost::shared_ptr<StackSparseMatrix> op = leftBlock->get_op_rep(HAM, deltaQuantum);


  //this action is only performed on the 0th process
  if (rightBlock->get_sites().size() == 0) 
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  else {
    //const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->getOverlap();
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    const boost::shared_ptr<StackSparseMatrix> Overlap = rightBlock->get_op_rep(OVERLAP, hq);
    SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *op, *Overlap, &b, &(b.get_stateInfo()), *this, 1.0);
  }
  

  if (rightBlock->get_sites().size() == 0) {
    //this is a special case where the right block is just a dummy block to make the effective wavefunction have spin 0
    //distributedaccumulate(*this);
    dmrginp.makeopsT -> stop();
    return;
  }

  op = rightBlock->get_op_rep(HAM, deltaQuantum);
  //const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->getOverlap();
  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  const boost::shared_ptr<StackSparseMatrix> Overlap = leftBlock->get_op_rep(OVERLAP, hq);
  SpinAdapted::operatorfunctions::TensorProduct(leftBlock, *Overlap, *op, &b, &(b.get_stateInfo()), *this, 1.0);

  // CCD_A*D_B + CCD_B*D_A + c.c. 
  FUNCTOR f = boost::bind(&stackopxop::cxcddcomp, leftBlock, _1, &b, op_array); 
  for_all_singlethread(rightBlock->get_op_array(CRE), f);

  f = boost::bind(&stackopxop::cxcddcomp, rightBlock, _1, &b, op_array); 
  for_all_singlethread(leftBlock->get_op_array(CRE), f);  

  if (dmrginp.hamiltonian() != HUBBARD) {    
    f = boost::bind(&stackopxop::cdxcdcomp, otherBlock, _1, &b, op_array);
    for_all_singlethread(loopBlock->get_op_array(CRE_DES), f);
    
    f = boost::bind(&stackopxop::ddxcccomp, otherBlock, _1, &b, op_array);
    for_all_singlethread(loopBlock->get_op_array(CRE_CRE), f);
  }
  //accumulateSinglethread(this, op_array, op_distributed, maxt);
  //accumulateMultiThread(this, op_array, 1);
  //distributedaccumulate(*this);

  dmrginp.makeopsT -> stop();    
  
}

double SpinAdapted::StackHam::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  const TwoElectronArray& v_2 = *(b->get_twoInt());
  double element = 0.0;
  bool finish = false;
  for (int i=0; i<ladder.size(); i++)
    {
      if (finish) break;
      bool isLallowed = Symmetry::spatial_cg(ladder[i].sym_is().getirrep(), 0, c1.sym_is().getirrep(), ladder[i].row(), 0, c1.row())!=0; // symmetry
      if ((c1.Sz != ladder[i].Sz || c1.S != ladder[i].S) || !isLallowed)
	continue;
      else
	finish = true; // only one element from the ladder has none zero matrix element with c1

      double matrixE = 0.0;
      for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) {
	for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++) {
	  Slater s1 = it1->first, s2 = it2->first; // slater determinants
	  double d1 = it1->second, d2 = it2->second; // weights

	  std::vector<int> cv, dv;
	  s1.connect(s2, cv, dv); // how to generate s1 from s2
	  if ((dv.size() == 2) && (cv.size() == 2)) {
	    int cI = cv[0]; 
	    int cJ = cv[1]; 
	    int dK = dv[0]; 
	    int dL = dv[1]; 
	    int parity = s1.trace(s2.d(dK).d(dL).c(cJ).c(cI));
	    double factor = parity*d1*d2*0.5;
	    matrixE += factor*(v_2(cI, cJ, dK, dL) - v_2(cJ, cI, dK, dL) - v_2(cI, cJ, dL, dK) + v_2(cJ, cI, dL, dK));
	  } else if ((cv.size() == 1) && (dv.size() == 1)) {
	    // from v1
	    int cI = cv[0]; 
	    int dK = dv[0]; 
	    int parity = s1.trace(s2.d(dK).c(cI));
	    double factor = parity*d1*d2;
	    matrixE += factor*v_1[b->get_integralIndex()](cI, dK);
	    // from v2
	    if(dmrginp.spinAdapted()) {
	      for (int kj=0; kj<b->get_sites().size(); kj++) {
		int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
		int num = 2*Symmetry::sizeofIrrep(SymmetryOfOrb(b->get_sites()[kj]).getirrep());
		for (int J = jindex; J<num+jindex; J++) {	    
		  s1 = it1->first; s2 = it2->first;
		  parity = s1.trace(s2.d(dK).d(J).c(J).c(cI));
		  factor = parity*d1*d2*0.5;
		  matrixE += factor*(v_2(cI, J, dK, J) - v_2(J, cI, dK, J) - v_2(cI, J, J, dK) + v_2(J, cI, J, dK));
		}
	      }
	    } else {
	      for (int kj=0; kj<b->get_sites().size(); kj++) {
		int J = b->get_sites()[kj];
		s1 = it1->first; s2 = it2->first;
		parity = s1.trace(s2.d(dK).d(J).c(J).c(cI));
		factor = parity*d1*d2*0.5;
		matrixE += factor*(v_2(cI, J, dK, J) - v_2(J, cI, dK, J) - v_2(cI, J, J, dK) + v_2(J, cI, J, dK));	      
	      }
	    }
	  } else if ((cv.size() == 0) && (dv.size() == 0)) {
	    if(dmrginp.spinAdapted()) { // spin adapted
	      //T
	      for (int kj=0; kj<b->get_sites().size(); kj++) {
		int jindex = dmrginp.spatial_to_spin()[b->get_sites()[kj]];
		int num = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kj]).getirrep());
		for (int J = jindex; J<num+jindex; J++) {	    
		  s1 = it1->first; s2 = it2->first;
		  matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1[b->get_integralIndex()](J,J);
		}
	      }
	      //V
	      for (int ki=0; ki<b->get_sites().size(); ki++)
		for (int kk=0; kk<b->get_sites().size(); kk++) {
		  int Iindex = dmrginp.spatial_to_spin()[b->get_sites()[ki]];
		  int Inum = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[ki]).getirrep());
		  for (int I = Iindex; I<Inum+Iindex; I++) {	    
		    int Kindex = dmrginp.spatial_to_spin()[b->get_sites()[kk]];
		    int Knum = 2*Symmetry::sizeofIrrep(SymmetryOfSpatialOrb(b->get_sites()[kk]).getirrep());
		    for (int K = Kindex; K<Knum+Kindex; K++) {	    
		      double factor = 0.5*d1*d2; //if (ki == kk) factor = 1.0*d1*d2;
		      s1 = it1->first; s2 = it2->first;
		      matrixE += factor*s1.trace(s2.d(I).d(K).c(K).c(I))*(v_2(I, K, I, K) - v_2(K, I, I, K));
		    }
		  }
		}
	    } else { // spin non-adapted
	      for (int kj=0; kj<b->get_sites().size(); kj++) {
		int J = b->get_sites()[kj];
		s1 = it1->first; s2 = it2->first;
		matrixE += d1*d2*s1.trace(s2.d(J).c(J))*v_1[b->get_integralIndex()](J,J);
	      }
	      //V
	      for (int ki=0; ki<b->get_sites().size(); ki++)
		for (int kk=0; kk<b->get_sites().size(); kk++) {
		  int K = b->get_sites()[kk];
		  int I = b->get_sites()[ki];
		  double factor = 0.5*d1*d2; //if (ki == kk) factor = 1.0*d1*d2;
		  s1 = it1->first; s2 = it2->first;
		  matrixE += factor*s1.trace(s2.d(I).d(K).c(K).c(I))*(v_2(I, K, I, K) - v_2(K, I, I, K));
		}
	    }
	  } else if (dmrginp.hamiltonian() == BCS && cv.size() == 4 && dv.size() == 0) {
	    int cI = cv[0];
	    int cJ = cv[1];
	    int cK = cv[2];
	    int cL = cv[3];
	    int parity = s1.trace(s2.c(cL).c(cK).c(cJ).c(cI));
	    double factor = parity*d1*d2;
	    matrixE += factor*v_cccc[b->get_integralIndex()](cI, cJ, cK, cL);
	  } else if (dmrginp.hamiltonian() == BCS && cv.size() == 0 && dv.size() == 4) {
	    int dI = dv[0];
	    int dJ = dv[1];
	    int dK = dv[2];
	    int dL = dv[3];
	    int parity = s1.trace(s2.d(dL).d(dK).d(dJ).d(dI));
	    double factor = parity*d1*d2;
	    matrixE += factor*v_cccc[b->get_integralIndex()](dL, dK, dJ, dI);  
	  } else if (dmrginp.hamiltonian() == BCS && cv.size() == 3 && dv.size() == 1) {
	    int cI = cv[0];
	    int cJ = cv[1];
	    int cK = cv[2];
	    int dL = dv[0];
	    int parity = s1.trace(s2.d(dL).c(cK).c(cJ).c(cI));
	    double factor = parity*d1*d2;
	    matrixE += factor*v_cccd[b->get_integralIndex()](cI,cJ,cK,dL);
	  } else if (dmrginp.hamiltonian() == BCS && cv.size() == 1 && dv.size() == 3) {
	    int cI = cv[0];
	    int dJ = dv[0];
	    int dK = dv[1];
	    int dL = dv[2];
	    int parity = s1.trace(s2.d(dL).d(dK).d(dJ).c(cI));
	    double factor = parity*d1*d2;
	    matrixE += factor*v_cccd[b->get_integralIndex()](dL,dK,dJ,cI);
	  } else if (dmrginp.hamiltonian() == BCS && cv.size() == 2 && dv.size() == 0) {
	    // from v_cc pairing
	    int cI = cv[0];
	    int cJ = cv[1];
	    int parity = s1.trace(s2.c(cJ).c(cI));
	    double factor = parity*d1*d2;
	    matrixE += factor * (v_cc[b->get_integralIndex()](cI,cJ)-v_cc[b->get_integralIndex()](cJ,cI));
	    // from v_cccd
	    if (dmrginp.spinAdapted()) {
	      pout << "Oops... BCS+SpinAdaption not implemented yet!" << endl;
	      abort();
	    } else {
	      for (int kl = 0; kl < b->get_sites().size(); ++kl) {
		int K = b->get_sites()[kl];
		s1 = it1->first; s2 = it2->first;
		parity = s1.trace(s2.d(K).c(K).c(cJ).c(cI));
		factor = parity*d1*d2;
		matrixE += factor*v_cccd[b->get_integralIndex()](cI,cJ,K,K);
	      }
	    }
	  } else if (dmrginp.hamiltonian() == BCS && cv.size() == 0 && dv.size() == 2) {
	    // from v_cc pairing
	    int dI = dv[0];
	    int dJ = dv[1];
	    int parity = s1.trace(s2.d(dI).d(dJ));
	    double factor = parity*d1*d2;
	    matrixE += factor * (v_cc[b->get_integralIndex()](dI,dJ)-v_cc[b->get_integralIndex()](dJ,dI));
	    // from v_cccd
	    if (dmrginp.spinAdapted()) {
	      pout << "Oops... BCS+SpinAdaption not implemented yet!" << endl;
	      abort();
	    } else {
	      for (int kl = 0; kl < b->get_sites().size(); ++kl) {
		int K = b->get_sites()[kl];
		s1 = it1->first; s2 = it2->first;
		parity = s1.trace(s2.d(dI).d(dJ).d(K).c(K));
		factor = parity*d1*d2;
		matrixE += factor*v_cccd[b->get_integralIndex()](dI,dJ,K,K);;
	      }
	    }
	  }
	}
      }
      element += 	matrixE;
    }
  return element;
}


//******************Overlap*****************
void SpinAdapted::StackOverlap::build(StackMatrix& m, int row, int col, const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0 || memoryUsed() != 0) {
    m = operator_element(row, col);
    return;
  }

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  boost::shared_ptr<StackSparseMatrix> op = leftBlock->get_op_rep(OVERLAP, deltaQuantum);

  if (rightBlock->get_sites().size() == 0) 
    SpinAdapted::operatorfunctions::TensorTraceElement(leftBlock, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);
  else {
    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(OVERLAP, deltaQuantum);
    SpinAdapted::operatorfunctions::TensorProductElement(rightBlock, *op2, *op, &b, &(b.get_stateInfo()), *this, m, row, col, 1.0);

  }

  return;


}

void SpinAdapted::StackOverlap::build(const StackSpinBlock& b)
{
  if (b.get_rightBlock() == 0) return; //already built
  dmrginp.makeopsT -> start();

  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  memset(data, 0, totalMemory * sizeof(double));


#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  boost::shared_ptr<StackSparseMatrix> op = leftBlock->get_op_rep(OVERLAP, deltaQuantum);

  if (rightBlock->get_sites().size() == 0) 
    SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
  else {
    boost::shared_ptr<StackSparseMatrix> op2 = rightBlock->get_op_rep(OVERLAP, deltaQuantum);
    SpinAdapted::operatorfunctions::TensorProduct(rightBlock, *op2, *op, &b, &(b.get_stateInfo()), *this, 1.0);

  }
  dmrginp.makeopsT -> stop();

  return;


}

double SpinAdapted::StackOverlap::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  double element = 0.0;
  bool finish = false;
  for (int i=0; i<ladder.size(); i++)
    {
      if (finish) break;
      bool isLallowed = Symmetry::spatial_cg(ladder[i].sym_is().getirrep(), 0, c1.sym_is().getirrep(), ladder[i].row(), 0, c1.row())!=0; // symmetry
      if ((c1.Sz != ladder[i].Sz || c1.S != ladder[i].S) || !isLallowed)
	continue;
      else
	finish = true; // only one element from the ladder has none zero matrix element with c1

      double matrixE = 0.0;
      for(map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) 
	for (map<Slater, double>::iterator it2 = ladder[i].det_rep.begin(); it2 != ladder[i].det_rep.end(); it2++) {
	  Slater s1 = it1->first, s2 = it2->first; // slater determinants
	  double d1 = it1->second, d2 = it2->second; // weights
	
	  std::vector<int> cv, dv;
	  s1.connect(s2, cv, dv); // how to generate s1 from s2
	  matrixE += d1*d2*s1.trace(s2);
	}
      element += 	matrixE;
    }
  return element;
}
