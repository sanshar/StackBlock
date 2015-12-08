/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Extension of Operators.C for 4-index operators
//FIXME there's a lot of duplication, especially in build_from_disk... Templates??
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "Stack_op_components.h"
#include "StackBaseOperator.h"
#include "Stackspinblock.h"
#include "operatorfunctions.h"
#include "tensor_operator.h"
#include "four_index_ops.h"

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
 void SpinAdapted::StackCreCreDesDes::build(const StackSpinBlock& b) { 
      dmrginp.makeopsT -> start();
      built = true;
      allocate(b.get_braStateInfo(), b.get_ketStateInfo());
    
      const int i = get_orbs()[0];
      const int j = get_orbs()[1];
      const int k = get_orbs()[2];
      const int l = get_orbs()[3];
      StackSpinBlock* leftBlock = b.get_leftBlock();
      StackSpinBlock* rightBlock = b.get_rightBlock();
    
      if (leftBlock->get_op_array(CRE_CRE_DES_DES).has(i,j,k,l))
      {      
        const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES_DES, quantum_ladder, i,j,k,l);
        if (rightBlock->get_sites().size() == 0) 
          SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
        dmrginp.makeopsT -> stop();
        return;
      }
      assert(false && "Only build CRECREDESDES in the starting block when spin-embeding is used");
    }

double SpinAdapted::StackCreCreDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CC)(D))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(D))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
//FIXME components != 0
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 
  TensorOp D4(L,-1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCD = CC.product(D3, spin123, irrep123);
  TensorOp CCDD = CCD.product(D4, spin1234, irrep1234);

//FIXME loop over deltaQuantum components
  int j = 0;
  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[j], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCDD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreCreDesDes::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreCreDesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreCreDesDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreDesCreDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_CRE_DES).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_CRE_DES, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESCREDES in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreDesCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CD)(C))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(C))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 
  TensorOp D4(L,-1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDC = CD.product(C3, spin123, irrep123);
  TensorOp CDCD = CDC.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDCD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreDesCreDes::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreDesCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreDesCreDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreDesDesCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_DES_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DES_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESDESCRE in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreDesDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CD)(D))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(D))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 
  TensorOp C4(L, 1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDD = CD.product(D3, spin123, irrep123);
  TensorOp CDDC = CDD.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDDC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreDesDesCre::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreDesDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreDesDesCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Des,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreDesDesDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_DES_DES).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_DES_DES, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESDESDES in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreDesDesDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CD)(D))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(D))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp D3(K,-1); 
  TensorOp D4(L,-1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDD = CD.product(D3, spin123, irrep123);
  TensorOp CDDD = CDD.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDDD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreDesDesDes::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreDesDesDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreDesDesDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Cre,Des)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreCreCreDes::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE_CRE_DES).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_CRE_DES, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CRECRECREDES in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreCreCreDes::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CC)(C))(D))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(C))(D))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp C3(K, 1); 
  TensorOp D4(L,-1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCC = CC.product(C3, spin123, irrep123);
  TensorOp CCCD = CCC.product(D4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCCD, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreCreCreDes::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreCreCreDes>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreCreCreDes);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Des,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreCreDesCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE_DES_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_DES_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CRECREDESCRE in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreCreDesCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CC)(D))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(D))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)D
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp D3(K,-1); 
  TensorOp C4(L, 1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCD = CC.product(D3, spin123, irrep123);
  TensorOp CCDC = CCD.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCDC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreCreDesCre::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreCreDesCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreCreDesCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Des,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreDesCreCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_DES_CRE_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_DES_CRE_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CREDESCRECRE in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreDesCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CD)(C))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CD)(C))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CD
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CD)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp D2(J,-1); 
  TensorOp C3(K, 1); 
  TensorOp C4(L, 1); 

  TensorOp CD = C1.product(D2, spin12, irrep12);
  TensorOp CDC = CD.product(C3, spin123, irrep123);
  TensorOp CDCC = CDC.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CDCC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreDesCreCre::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreDesCreCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreDesCreCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//  (Cre,Cre,Cre,Cre)
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
void SpinAdapted::StackCreCreCreCre::build(const StackSpinBlock& b) { 
  dmrginp.makeopsT -> start();
  built = true;
  allocate(b.get_braStateInfo(), b.get_ketStateInfo());

  const int i = get_orbs()[0];
  const int j = get_orbs()[1];
  const int k = get_orbs()[2];
  const int l = get_orbs()[3];
  StackSpinBlock* leftBlock = b.get_leftBlock();
  StackSpinBlock* rightBlock = b.get_rightBlock();

  if (leftBlock->get_op_array(CRE_CRE_CRE_CRE).has(i,j,k,l))
  {      
    const boost::shared_ptr<StackSparseMatrix>& op = leftBlock->get_op_rep(CRE_CRE_CRE_CRE, quantum_ladder, i,j,k,l);
    if (rightBlock->get_sites().size() == 0) 
      SpinAdapted::operatorfunctions::TensorTrace(leftBlock, *op, &b, &(b.get_stateInfo()), *this);
    dmrginp.makeopsT -> stop();
    return;
  }
  assert(false && "Only build CRECRECRECRE in the starting block when spin-embeding is used");
}

double SpinAdapted::StackCreCreCreCre::redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b)
{
  assert( build_pattern == "(((CC)(C))(C))" );
  double element = 0.0;
  int I = get_orbs()[0]; 
  int J = get_orbs()[1];
  int K = get_orbs()[2];
  int L = get_orbs()[3];
  int Slaterlength = c1.det_rep.begin()->first.size();
  vector<bool> backupSlater1(Slaterlength,0), backupSlater2(Slaterlength,0);

  // Must take into account how the 4-index is built from a combination of the 2-index ops
  std::vector<SpinQuantum> quantum_ladder = get_quantum_ladder().at("(((CC)(C))(C))");
  assert( quantum_ladder.size() == 3 );

  SpinQuantum deltaQuantum12 = quantum_ladder.at(0);
  SpinQuantum deltaQuantum123 = quantum_ladder.at(1);
  SpinQuantum deltaQuantum1234 = quantum_ladder.at(2);
  deltaQuantum[0] = deltaQuantum1234;

  // Spin quantum data for CC
  IrrepSpace sym12 = deltaQuantum12.get_symm();
  int irrep12 = deltaQuantum12.get_symm().getirrep();
  int spin12 = deltaQuantum12.get_s().getirrep();
  // Spin quantum data for (CC)C
  IrrepSpace sym123 = deltaQuantum123.get_symm();
  int irrep123 = deltaQuantum123.get_symm().getirrep();
  int spin123= deltaQuantum123.get_s().getirrep();
  // Spin quantum data for total operator
  IrrepSpace sym1234 = deltaQuantum1234.get_symm();
  int irrep1234 = deltaQuantum1234.get_symm().getirrep();
  int spin1234 = deltaQuantum1234.get_s().getirrep();

  TensorOp C1(I, 1); 
  TensorOp C2(J, 1); 
  TensorOp C3(K, 1); 
  TensorOp C4(L, 1); 

  TensorOp CC = C1.product(C2, spin12, irrep12);
  TensorOp CCC = CC.product(C3, spin123, irrep123);
  TensorOp CCCC = CCC.product(C4, spin1234, irrep1234);

  for (int i=0; i<ladder.size(); i++)
  {
    int index = 0; double cleb=0.0;
    if (nonZeroTensorComponent(c1, deltaQuantum[0], ladder[i], index, cleb)) {
      std::vector<double> MatElements = calcMatrixElements(c1, CCCC, ladder[i], backupSlater1, backupSlater2) ;
      element = MatElements[index]/cleb;
      break;
    }
    else
      continue;
  }
  return element;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

boost::shared_ptr<SpinAdapted::StackSparseMatrix> SpinAdapted::StackCreCreCreCre::getworkingrepresentation(const StackSpinBlock* block)
{
  assert(this->get_initialised());
  if (this->get_built()) {
    return boost::shared_ptr<StackCreCreCreCre>(this, boostutils::null_deleter()); // boost::shared_ptr does not own op
  }
  else {
    boost::shared_ptr<StackSparseMatrix> rep(new StackCreCreCreCre);
    *rep = *this;
    rep->build(*block);

    return rep;
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
