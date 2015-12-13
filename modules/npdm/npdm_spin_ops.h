/*
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012
Copyright (c) 2012, Garnet K.-L. Chan

This program is integrated in Molpro with the permission of
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_SPIN_OPS_H
#define NPDM_SPIN_OPS_H

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

#include "StackBaseOperator.h"
#include "operatorfunctions.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include "Stackspinblock.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

class NpdmSpinOps_base {

  public:
//FIXME
    NpdmSpinOps_base() = default;
// FIXME Shallow copy constructor
    NpdmSpinOps_base( const NpdmSpinOps_base & obj ) {
      //opReps_ = obj.opReps_;  //dont copy opReps_
      //opReps_.resize(obj.opReps_.size());
      build_pattern_ = obj.build_pattern_;
      transpose_ = obj.transpose_;
      factor_ = obj.factor_;
      indices_ = obj.indices_;
      is_local_ = obj.is_local_;
      opIndices_ = obj.opIndices_;
    }
    virtual bool set_local_ops( int idx ) { abort(); }

    // Numerical representation of the operators for several total spins (e.g. 2-index op has two forms with spin-1/2 particles)
//FIXME should this be a reference?  Don't want to copy!!
    std::vector< boost::shared_ptr<StackSparseMatrix> > opReps_;
    // is_local_ == true mean operators are duplicated on all mpi ranks
    bool is_local_;
    // How the operator is built (e.g. 3-index from product of 2-index cre-cre and 1-index destruction)
    std::string build_pattern_;
    // Do we need to transpose the representation before using it?
    bool transpose_;
    // Do we need to multiply by any constant factors when using it (due to implicit use of commutation relations or such like)?
    double factor_;
    // Effective spatial orbital indices (since due to use of transposition / commutation may not match OpRep.get_orbs() etc)
    std::vector<int> indices_;

    //the number of indices of the operators stored in opReps_
    int opIndices_=0;

#ifndef SERIAL
//FIXME put in implementation file
//FIXME  Do like this since serialization of whole object broken!!
    std::vector< boost::mpi::request > isend_mpi_obj( int rank, unsigned int tag_lo, unsigned int tag_hi )
    {
//FIXME!! This seems brittle and unstable. e.g. moved order of calls and it broke!!  Issues with bools ???
      boost::mpi::communicator world;
      assert( tag_hi < boost::mpi::environment::max_tag() );
      std::vector< boost::mpi::request > reqs;
      int k = tag_lo;
      reqs.push_back( world.isend(rank, k++, build_pattern_) );
      reqs.push_back( world.isend(rank, k++, transpose_) );
      reqs.push_back( world.isend(rank, k++, factor_) );
      reqs.push_back( world.isend(rank, k++, indices_) );
      reqs.push_back( world.isend(rank, k++, is_local_) );
      for ( int i = 0; i < opReps_.size(); ++i) {
        reqs.push_back( world.isend(rank, k++, *(opReps_.at(i))) );
      }
      assert( k < tag_hi );
      return reqs;
    }
      
    std::vector< boost::mpi::request > irecv_mpi_obj( int rank, unsigned int tag_lo, unsigned int tag_hi, int size )
    {
      boost::mpi::communicator world;
      assert( tag_hi < boost::mpi::environment::max_tag() );
      std::vector< boost::mpi::request > reqs;
      assert( opReps_.size() == 0 );
      int k = tag_lo;
      reqs.push_back( world.irecv(rank, k++, build_pattern_) );
      reqs.push_back( world.irecv(rank, k++, transpose_) );
      reqs.push_back( world.irecv(rank, k++, factor_) );
      reqs.push_back( world.irecv(rank, k++, indices_) );
      reqs.push_back( world.isend(rank, k++, is_local_) );
      for ( int i = 0; i < size; ++i) {
        boost::shared_ptr<StackSparseMatrix> op (new StackCre);
        reqs.push_back( world.irecv(rank, k++, *op ) );
        opReps_.push_back(op);
      }
      assert( k < tag_hi );
      return reqs;
    }
#endif

};

//===========================================================================================================================================================
//FIXME constructors / destructors
class NpdmSpinOps : public NpdmSpinOps_base {

  public:
    int size() { return size_; }
    virtual bool set_local_ops( int idx ) { abort(); }

    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<NpdmSpinOps>(new NpdmSpinOps(*this));}
//FIXME public??
    // Input file stream for disk-based operators used to build NPDM
    //std::ifstream ifs_;
    //std::string ifile_;
    // Number of spatial orbital combinations
    int size_;

    NpdmSpinOps() = default;
 NpdmSpinOps( const NpdmSpinOps & obj ) : NpdmSpinOps_base(obj) {
      spinBlock_ = obj.spinBlock_;
    }

    virtual std::vector< std::vector<int> > get_indices() { abort(); }
//    virtual const std::vector< int >& get_1d_indices() { abort(); }

    StackSpinBlock* spinBlock_;
  protected:
    //FIXME
    boost::shared_ptr<StackSparseMatrix> build_compound_operator( bool is_fermion, int sign,
                                                             boost::shared_ptr<StackSparseMatrix> lhsOp,
                                                             boost::shared_ptr<StackSparseMatrix> rhsOp,
                                                             int ispin, std::vector<int> indices, bool transpose );



};

//===========================================================================================================================================================

}
}

#endif

