/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Extension of Operators.C for 3-index operators
//FIXME there's a lot of duplication, especially in build_from_disk... Templates??
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "Stack_op_components.h"
#include "StackBaseOperator.h"
#include "Stackspinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"
#include "three_index_ops.h"
#include "pario.h"
#include "build_3index_ops.h"

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackDesDesDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(DES_DES_DES).has(i,j,k))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_DES_DES, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op( DES_DES_DES, b, DES, DES_DES, DES, DES_DES, *this );
  }
}

double SpinAdapted::StackDesDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "((DD)(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((DD)(D))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

  TensorOp D1(I,-1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
  TensorOp DD = D1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp DDD = DD.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DDD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//===========================================================================================================================================================
// 3PDM operators
//===========================================================================================================================================================

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreCreDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(CRE_CRE_DES).has(i,j,k) )
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op( CRE_CRE_DES, b, CRE, CRE_CRE, DES, CRE_DES, *this );
  }

}

double SpinAdapted::StackCreCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
//pout << "building CreCreDes explicitly from CSF..\n";
  assert( build_pattern == "((CC)(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CC)(D))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

//pout << "redMatrixElement indices:\n";
//pout << I << "  " << J << "  " << K << std::endl;
//pout << "spin composition:\n";
//pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
//FIXME  TensorOp CC = C1.product(C2, spin12, irrep12);  I==J has no affect
  TensorOp CC = C1.product(C2, spin12, irrep12, I==J);
  // Combine with third operator
  TensorOp CCD = CC.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreDesDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(CRE_DES_DES).has(i,j,k) )
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DES, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op( CRE_DES_DES, b, CRE, CRE_DES, DES, DES_DES, *this );
  }

  //assert(false && "Only build CREDESDES in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
//pout << "building CreDesDes explicitly from CSF..\n";
  assert( build_pattern == "((CD)(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CD)(D))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

//FIXME
//pout << "redMatrixElement indices:\n";
//pout << I << "  " << J << "  " << K << std::endl;
//pout << "spin composition:\n";
//pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
  TensorOp CD = C1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp CDD = CD.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreDesCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(CRE_DES_CRE).has(i,j,k) )
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_CRE, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op( CRE_DES_CRE, b, CRE, CRE_DES, CRE, DES_CRE, *this );
  }

  //assert(false && "Only build CREDESCRE in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
//pout << "building CreDesCre in CSF space explicitly..\n";
  assert( build_pattern == "((CD)(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CD)(C))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  assert( deltaQuantum[0] == deltaQuantum123 );

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

//pout << "redMatrixElement indices:\n";
//pout << I << "  " << J << "  " << K << std::endl;
//pout << "spin composition:\n";
//pout << spin12/2.0 << "  " << spin123/2.0 << std::endl;
  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 

  // Combine first two operators
  TensorOp CD = C1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp CDC = CD.product(C3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreCreCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(CRE_CRE_CRE).has(i,j,k) )
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_CRE, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op( CRE_CRE_CRE, b, CRE, CRE_CRE, CRE, CRE_CRE, *this );
  }

  //assert(false && "Only build CRECRECRE in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
//pout << "building CreCreCre explicitly from CSF..\n";
//pout << "mpirank = " << mpigetrank() << endl;
  assert( build_pattern == "((CC)(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
//pout << "i,j,k = " << I << " " << J << " " << K << endl;
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((CC)(C))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
//FIXME components != 0
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

  TensorOp C1(I,1); 
  TensorOp C2(J,1); 
  TensorOp C3(K,1); 

  // Combine first two operators
//FIXME I=J argument
  TensorOp CC = C1.product(C2, spin12, irrep12);
  // Combine with third operator
  TensorOp CCC = CC.product(C3, spin123, irrep123);

//FIXME loop over deltaQuantum components
  int j = 0;
  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}


//===========================================================================================================================================================
// 4PDM operators
//===========================================================================================================================================================

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackDesCreDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(DES_CRE_DES).has(i,j,k) && rightBlock->get_sites().size() == 0)
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_CRE_DES, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op(DES_CRE_DES, b, DES, DES_CRE, DES, CRE_DES, *this );
  }
}

double SpinAdapted::StackDesCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "((DC)(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((DC)(D))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

  TensorOp D1(I,-1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 

  // Combine first two operators
  TensorOp DC = D1.product(C2, spin12, irrep12);
  // Combine with third operator
  TensorOp DCD = DC.product(D3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DCD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackDesDesCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0)
  {
    if (leftBlock->get_op_array(DES_DES_CRE).has(i,j,k))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_DES_CRE, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else {
    Three_index_ops::build_3index_single_op(DES_DES_CRE, b, DES, DES_DES, CRE, DES_CRE, *this );
  }
}

double SpinAdapted::StackDesDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "((DD)(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((DD)(C))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

  TensorOp D1(I,-1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 

  // Combine first two operators
//FIXME I=J argument??
  TensorOp DD = D1.product(D2, spin12, irrep12);
  // Combine with third operator
  TensorOp DDC = DD.product(C3, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DDC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Des,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackDesCreCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  //allocate(b.get_braStateInfo(), b.get_ketStateInfo());
  memset(data, 0, totalMemory * sizeof(double));

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (rightBlock->get_sites().size() == 0) {
    if (leftBlock->get_op_array(DES_CRE_CRE).has(i,j,k))
    {      
      const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(DES_CRE_CRE, quantum_ladder, i,j,k);
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
      dmrginp.makeopsT -> stop();
      return;
    }
  }
  else
  {
    Three_index_ops::build_3index_single_op(DES_CRE_CRE, b, DES, DES_CRE, CRE, CRE_CRE, *this );
  }
}

double SpinAdapted::StackDesCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "((D)(CC))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 3-index is built from a combination of 2-index and 1-index
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("((D)(CC))");
  assert( quantum_ladder.size() == 2 );
  SpinQuantum deltaQuantum23 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  deltaQuantum[0] = deltaQuantum123;

  // Spin quantum data for first pair of operators combined
  IrrepSpace sym23 = deltaQuantum23.get_symm();
  int irrep23 = deltaQuantum23.get_symm().getirrep();
  int spin23 = deltaQuantum23.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123 = deltaQuantum123.get_s().getirrep();

  TensorOp D1(I,-1); 
  TensorOp C2(J, 1); 
  TensorOp C3(K, 1); 

  // Combine first two operators
//FIXME I=J argument??
  TensorOp CC = C2.product(C3, spin23, irrep23);
  // Combine with third operator
  TensorOp DCC = D1.product(CC, spin123, irrep123);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, DCC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------

