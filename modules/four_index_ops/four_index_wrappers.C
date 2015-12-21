/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "npdm_patterns.h"
#include "npdm_operators.h"

//FIXME streamline code and reuse??  functions are almost identical here
//template these???
namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Npdm_op_wrapper_CCDD::Npdm_op_wrapper_CCDD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCDD::set_local_ops( int idx )
{
//pout << "getting CCDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME combine 3-index and 4-indexx read from disk.... UNIFY......
    int size = spinBlock_->get_op_array(CRE_CRE_DES_DES).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDCD::Npdm_op_wrapper_CDCD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDCD::set_local_ops( int idx )
{
//pout << "getting CDCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE_DES).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";
//FIXME
    int size = spinBlock_->get_op_array(CRE_DES_CRE_DES).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( jx == kx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDDC::Npdm_op_wrapper_CDDC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDDC::set_local_ops( int idx )
{
//pout << "getting CDDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES_CRE).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME
    int size = spinBlock_->get_op_array(CRE_DES_DES_CRE).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( (kx == lx) || (jx == lx) ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDDD::Npdm_op_wrapper_CDDD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDDD::set_local_ops( int idx )
{
//pout << "getting CDDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES_DES).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME
    int size = spinBlock_->get_op_array(CRE_DES_DES_DES).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  //if ( kx == lx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCCD::Npdm_op_wrapper_CCCD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCCD::set_local_ops( int idx )
{
//pout << "getting CCCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME
    int size = spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  //if ( kx == lx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCDC::Npdm_op_wrapper_CCDC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCDC::set_local_ops( int idx )
{
//pout << "getting CCDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME
    int size = spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( kx == lx ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDCC::Npdm_op_wrapper_CDCC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDCC::set_local_ops( int idx )
{
//pout << "getting CDCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME
    int size = spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);
  // Skip if commutation problematic
  if ( (jx==kx) || (jx==lx) ) return true;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCCC::Npdm_op_wrapper_CCCC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=4;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCCC::set_local_ops( int idx )
{
//pout << "getting CCCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx, lx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() ) {
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_local_element(idx);
  }
  else {
//pout << "...from disk...\n";

//FIXME combine 3-index and 4-indexx read from disk.... UNIFY......
    int size = spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_local_element(idx).size();
    assert( size == 6 );
    opReps_.clear();
    // Read in all spin components for this set of spatial indices (note order matters)
    for ( int i=0; i<size; i++) {
      boost::shared_ptr<StackSparseMatrix> op (new StackCre);
      //boost::archive::binary_iarchive load_op(ifs_);
      //load_op >> *op;
      assert( op->get_built_on_disk() );
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  lx = opReps_.at(0)->get_orbs(3);

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  indices_.push_back( lx );

  return false;
}

//===========================================================================================================================================================

}
}

