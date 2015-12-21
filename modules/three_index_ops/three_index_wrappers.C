/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include "npdm_patterns.h"
#include "npdm_operators.h"
#include "pario.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================
// 3PDM operators
//===========================================================================================================================================================

Npdm_op_wrapper_CCC::Npdm_op_wrapper_CCC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCC::set_local_ops( int idx )
{
//pout << "getting CCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    for (int i = 0; i < opReps_tmp.size(); i++) {
      //string filename = spinBlock_->get_op_array(CRE_CRE_CRE).get_filename()+to_string(i);
      string filename = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx)[i]->get_filename();
      boost::shared_ptr<StackSparseMatrix> op (new StackCreCreCre);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);

      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
//pout << "indices  " << ix << " " << jx << " " << kx << std::endl;
//pout << "build pattern " << opReps_.at(0)->get_build_pattern() << std::endl;
//pout << "2a CCC operator elements:\n";
//pout << *(opReps_[0]);
//pout << "2b CCC operator elements:\n";
//pout << *(opReps_[1]);
//pout << "4  CCC operator elements:\n";
//pout << *(opReps_[2]);
//  //assert( build_pattern_ == opReps_.at(0)->get_build_pattern() );

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( (ix == jx) && (jx == kx) ) {
    //FIXME This operator should not be built; it's zero as we cannot create 3 spin-1/2 particles with different spins
    //pout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CCD::Npdm_op_wrapper_CCD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CCD::set_local_ops( int idx )
{
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core() || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;

    opReps_tmp = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {

      string filename = spinBlock_->get_op_array(CRE_CRE_DES).get_local_element(idx)[i]->get_filename();
      boost::shared_ptr<StackSparseMatrix> op (new StackCreCreDes);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }

  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();
  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDD::Npdm_op_wrapper_CDD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  // For disk-based storage
  opIndices_=3;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDD::set_local_ops( int idx )
{
//pout << "getting CDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(CRE_DES_DES).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_DES_DES).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {


      string filename = spinBlock_->get_op_array(CRE_DES_DES).get_local_element(idx)[i]->get_filename();
      //string filename = spinBlock_->get_op_array(CRE_DES_DES).get_filename()+to_string(i);
      boost::shared_ptr<StackSparseMatrix> op (new StackCreDesDes);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);
  build_pattern_ = opReps_.at(0)->get_build_pattern();


  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  return false;

//FIXME GET AS TRANSPOSE of CCD and compare!!!
}

//===========================================================================================================================================================

Npdm_op_wrapper_CDC::Npdm_op_wrapper_CDC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(CRE_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(CRE_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_CDC::set_local_ops( int idx )
{
//pout << "getting CDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(CRE_DES_CRE).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(CRE_DES_CRE).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
      string filename = spinBlock_->get_op_array(CRE_DES_CRE).get_local_element(idx)[i]->get_filename();
      //string filename = spinBlock_->get_op_array(CRE_DES_CRE).get_filename()+to_string(i);
      boost::shared_ptr<StackSparseMatrix> op (new StackCreDesCre);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

    ix = opReps_.at(0)->get_orbs(0);
    jx = opReps_.at(0)->get_orbs(1);
    kx = opReps_.at(0)->get_orbs(2);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
//pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( jx == kx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================
// 4PDM operators
//===========================================================================================================================================================

Npdm_op_wrapper_DCD::Npdm_op_wrapper_DCD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_CRE_DES).get_size();
  is_local_ = spinBlock_->get_op_array(DES_CRE_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DCD::set_local_ops( int idx )
{
//pout << "getting DCD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(DES_CRE_DES).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_CRE_DES).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
      string filename = spinBlock_->get_op_array(DES_CRE_DES).get_local_element(idx)[i]->get_filename();
      //string filename = spinBlock_->get_op_array(DES_CRE_DES).get_filename()+to_string(i);
      boost::shared_ptr<StackSparseMatrix> op (new StackDesCreDes);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

    ix = opReps_.at(0)->get_orbs(0);
    jx = opReps_.at(0)->get_orbs(1);
    kx = opReps_.at(0)->get_orbs(2);
  build_pattern_ = opReps_.at(0)->get_build_pattern();
//pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( ix == jx ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================

Npdm_op_wrapper_DDC::Npdm_op_wrapper_DDC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_DES_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(DES_DES_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DDC::set_local_ops( int idx )
{
//pout << "getting DDC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    for (int i = 0; i < opReps_tmp.size(); i++) {
      string filename = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx)[i]->get_filename();
      //string filename = spinBlock_->get_op_array(DES_DES_CRE).get_filename()+to_string(i);
      boost::shared_ptr<StackSparseMatrix> op (new StackDesDesCre);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

    ix = opReps_.at(0)->get_orbs(0);
    jx = opReps_.at(0)->get_orbs(1);
    kx = opReps_.at(0)->get_orbs(2);
  build_pattern_ = opReps_.at(0)->get_build_pattern();

//pout << "indices  " << ix << " " << jx << " " << kx << std::endl;

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( (jx == kx) || (ix == kx) ) {
    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }
  return false;
}

//===========================================================================================================================================================
// FIXME This DCC operator wrapper using tranpose of DDC seems to need extra minus sign 
//
////Npdm_op_wrapper_DCC::Npdm_op_wrapper_DCC( StackSpinBlock * spinBlock )
////{
////  opReps_.clear();
////  indices_.clear();
////  spinBlock_ = spinBlock;
////  size_ = spinBlock_->get_op_array(DES_DES_CRE).get_size();
////  is_local_ = spinBlock_->get_op_array(DES_DES_CRE).is_local();
////  factor_ = 1.0;
////  transpose_ = true;
////  build_pattern_ = "0";
////  // For disk-based storage
////  ifile_ = spinBlock_->get_op_array(DES_DES_CRE).get_filename();
////}
////
//////-----------------------------------------------------------------------------------------------------------------------------------------------------------
////
////bool Npdm_op_wrapper_DCC::set_local_ops( int idx )
////{
////pout << "getting DCC operator...\n";
////  // Spatial orbital indices
////  indices_.clear();
////  int ix, jx, kx;
////
////  // Read in operator representations from disk or memory
////  if ( dmrginp.do_npdm_in_core() )
////    opReps_ = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
////  else {
////    //assert( check_file_open( idx ) );
////    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
////    opReps_tmp = spinBlock_->get_op_array(DES_DES_CRE).get_local_element(idx);
////    //assert( opReps_tmp.at(0)->get_built_on_disk() );
////    opReps_.clear();
////    // Read in full spin-set from disk
////    for (int i = 0; i < opReps_tmp.size(); i++) {
////       boost::archive::binary_iarchive load_op(ifs_);
////       boost::shared_ptr<StackSparseMatrix> op (new StackCre);
////       load_op >> *op;
////       opReps_.push_back(op);
////    }
////    //assert( check_file_close( idx ) );
////  }
////
////  std::string tmp = opReps_tmp.at(0)->get_build_pattern();
////  if ( tmp == "((DD)C)" ) build_pattern_ = "(D(CC))"; //  <--------- seem to need minus sign here, don't know why
////  else if ( tmp == "(D(DC))" ) build_pattern_ = "((DC)C)";
////  else abort();
////
////  ix = opReps_tmp.at(0)->get_orbs(0);
////  jx = opReps_tmp.at(0)->get_orbs(1);
////  kx = opReps_tmp.at(0)->get_orbs(2);
////
////  // Note use of transpose means we store this as (k,j,i) not (i,j,k)
////  indices_.push_back( kx );
////  indices_.push_back( jx );
////  indices_.push_back( ix );
////  if ( (jx == kx) || (ix == kx) ) {
////    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
////    //pout << "WARNING: skipping this operator\n";
////    return true;
////  }
////  return false;
////}
////
//===========================================================================================================================================================

Npdm_op_wrapper_DCC::Npdm_op_wrapper_DCC( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_CRE_CRE).get_size();
  is_local_ = spinBlock_->get_op_array(DES_CRE_CRE).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DCC::set_local_ops( int idx )
{
//pout << "getting DCC operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

//FIXME SKIP THEM WHEN SOME INDICES ARE EQUAL ?
//return true;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(DES_CRE_CRE).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_CRE_CRE).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    ix = opReps_tmp.at(0)->get_orbs(0);
    jx = opReps_tmp.at(0)->get_orbs(1);
    kx = opReps_tmp.at(0)->get_orbs(2);
    for (int i = 0; i < opReps_tmp.size(); i++) {
      string filename = spinBlock_->get_op_array(DES_CRE_CRE).get_local_element(idx)[i]->get_filename();
      //string filename = spinBlock_->get_op_array(DES_CRE_CRE).get_filename()+to_string(i);
      boost::shared_ptr<StackSparseMatrix> op (new StackDesCreCre);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

    ix = opReps_.at(0)->get_orbs(0);
    jx = opReps_.at(0)->get_orbs(1);
    kx = opReps_.at(0)->get_orbs(2);
  build_pattern_ = opReps_.at(0)->get_build_pattern();

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );
  if ( (ix == jx) || (ix == kx) ) {
////    //FIXME I think this fails because of potential problems commuting operators with same indices in spin-transformation
    //pout << "WARNING: skipping this operator\n";
    return true;
  }

  return false;
}

//===========================================================================================================================================================
// Build using Transpose(CCC) -- check rationale of minus signs!
////
////Npdm_op_wrapper_DDD::Npdm_op_wrapper_DDD( StackSpinBlock * spinBlock )
////{
////  opReps_.clear();
////  indices_.clear();
////  spinBlock_ = spinBlock;
////  size_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_size();
////  is_local_ = spinBlock_->get_op_array(CRE_CRE_CRE).is_local();
////  //FIXME Minus sign here because of a (DD) transpose ??
////  factor_ = -1.0;
////  transpose_ = true;
////  build_pattern_ = "0";
////  // For disk-based storage
////  ifile_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_filename();
////}
////
//////-----------------------------------------------------------------------------------------------------------------------------------------------------------
////
////bool Npdm_op_wrapper_DDD::set_local_ops( int idx )
////{
//////pout << "getting DDD operator...\n";
////  // Spatial orbital indices
////  indices_.clear();
////  int ix, jx, kx;
////
////  // Read in operator representations from disk or memory
////  if ( dmrginp.do_npdm_in_core() )
////    opReps_ = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
////  else {
////abort();
////    //assert( check_file_open( idx ) );
////    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
////    opReps_tmp = spinBlock_->get_op_array(CRE_CRE_CRE).get_local_element(idx);
////    //assert( opReps_tmp.at(0)->get_built_on_disk() );
////    opReps_.clear();
////    // Read in full spin-set from disk
////    for (int i = 0; i < opReps_tmp.size(); i++) {
////       boost::archive::binary_iarchive load_op(ifs_);
////       boost::shared_ptr<StackSparseMatrix> op (new StackCre);
////       load_op >> *op;
////       opReps_.push_back(op);
////    }
////    //assert( check_file_close( idx ) );
////  }
////
////  std::string tmp = opReps_.at(0)->get_build_pattern();
////  if ( tmp == "((CC)C)" ) build_pattern_ = "(D(DD))";
////  else if ( tmp == "(C(CC))" ) build_pattern_ = "((DD)D)";
////  else abort();
////
////  ix = opReps_tmp.at(0)->get_orbs(0);
////  jx = opReps_tmp.at(0)->get_orbs(1);
////  kx = opReps_tmp.at(0)->get_orbs(2);
////
////  // Note use of transpose means we store this as (k,j,i) not (i,j,k)
////  indices_.push_back( kx );
////  indices_.push_back( jx );
////  indices_.push_back( ix );
////  if ( (ix == jx) && (jx == kx) ) {
////    //FIXME This operator should not be built; it's zero as we cannot destroy 3 spin-1/2 particles with different spins
////    return true;
////  }
////  return false;
////}
////
//===========================================================================================================================================================

Npdm_op_wrapper_DDD::Npdm_op_wrapper_DDD( StackSpinBlock * spinBlock )
{
  opReps_.clear();
  indices_.clear();
  spinBlock_ = spinBlock;
  size_ = spinBlock_->get_op_array(DES_DES_DES).get_size();
  is_local_ = spinBlock_->get_op_array(DES_DES_DES).is_local();
  factor_ = 1.0;
  transpose_ = false;
  build_pattern_ = "0";
  opIndices_=3;
  // For disk-based storage

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

bool Npdm_op_wrapper_DDD::set_local_ops( int idx )
{
//pout << "getting DDD operator...\n";
  // Spatial orbital indices
  indices_.clear();
  int ix, jx, kx;

  // Read in operator representations from disk or memory
  if ( dmrginp.do_npdm_in_core()  || spinBlock_->size() <= 1)
    opReps_ = spinBlock_->get_op_array(DES_DES_DES).get_local_element(idx);
  else {
    //if ( ! check_file_open(idx) ) abort();
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_tmp;
    opReps_tmp = spinBlock_->get_op_array(DES_DES_DES).get_local_element(idx);
    //assert( opReps_tmp.at(0)->get_built_on_disk() );
    opReps_.clear();
    // Read in full spin-set from disk
    for (int i = 0; i < opReps_tmp.size(); i++) {
      string filename = spinBlock_->get_op_array(DES_DES_DES).get_local_element(idx)[i]->get_filename();
      //string filename = spinBlock_->get_op_array(DES_DES_DES).get_filename()+to_string(i);
      boost::shared_ptr<StackSparseMatrix> op (new StackDesDesDes);
      op->set_filename() = filename;
      op->LoadThreadSafe(true);
      
      op->allocateOperatorMatrix(); 
      opReps_.push_back(op);
    }
    //if ( ! check_file_close(idx) ) abort();
  }

  build_pattern_ = opReps_.at(0)->get_build_pattern();

  ix = opReps_.at(0)->get_orbs(0);
  jx = opReps_.at(0)->get_orbs(1);
  kx = opReps_.at(0)->get_orbs(2);

  indices_.push_back( ix );
  indices_.push_back( jx );
  indices_.push_back( kx );

  return false;
}

//===========================================================================================================================================================

}
}

