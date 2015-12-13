#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

#include "Stackwavefunction.h"
#include "npdm_driver.h"
#include "npdm_patterns.h"
#include "npdm_expectations.h"
#include "pario.h"
#include <stdio.h>
#include "distribute.h"

namespace SpinAdapted{
namespace Npdm{

// DEBUG only
extern double DEBUG_COMM_TIME;
extern int DEBUG_CALL_GET_EXPECT;
extern double DEBUG_STORE_ELE_TIME;

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( StackSpinBlock * spinBlock,const std::vector<Npdm::CD> & cd_type );

//===========================================================================================================================================================

#ifndef SERIAL
unsigned int get_mpi_tag( int rank0, int rank1, int lda )
{
  unsigned int tag = rank0 * lda + rank1;
  assert( tag < 42949672 );
  return 100 * tag;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef SERIAL
int Npdm_driver::get_mpi_max_size( int my_size )
{
  int maxsize;
  boost::mpi::communicator world;
  std::vector<int> all_sizes;
  boost::mpi::all_gather(world, my_size, all_sizes);
  maxsize = *std::max_element( all_sizes.begin(), all_sizes.end() );
  assert( my_size <= maxsize );
  return maxsize;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Do we parallelize by broadcasting LHS or RHS operators?  This is a very simple heuristic for now.  With disk-access, it may not be so good.

#ifndef SERIAL
bool Npdm_driver::broadcast_lhs( int lhs_size, int rhs_size )
{
  // Note all ranks have to make the same decision!
  bool do_lhs = true;
  int lhs_maxsize = get_mpi_max_size( lhs_size );
  int rhs_maxsize = get_mpi_max_size( rhs_size );
  if (rhs_maxsize < lhs_maxsize) do_lhs = false;
  return do_lhs;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef SERIAL
bool Npdm_driver::skip_this_mpi_rank( NpdmSpinOps & lhsOps, NpdmSpinOps & rhsOps )
{
  boost::mpi::communicator world;
  bool skip = ( mpigetrank() > 0    && 
                lhsOps.is_local_    && 
                rhsOps.is_local_ );
  return skip;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef SERIAL
bool Npdm_driver::skip_parallel( NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, bool lhsrhsdot )
{
  boost::mpi::communicator world;
  // Don't want to parallelize if any of these conditions hold (either makes no sense or creates duplicated work)
  // Assumes 1-index ops are duplicated on all mpi ranks //FIXME check this?
  // (There might be some overlap in these criteria)
  bool skip = ( world.size() == 1               ||
                outerOps.is_local_  );
  return skip;
}
#endif
  
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Originally wrote this routine for LHS and RHS operators, but also works interchanging "lhs" and "rhs" everywhere when called in reverse

#ifndef SERIAL
void Npdm_driver::do_parallel_lhs_loop( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                                        NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps, bool skip )
{
  boost::mpi::communicator world;

  // Parallelize by broadcasting LHS ops
  NpdmSpinOps_base local_base(outerOps);
  std::vector< boost::mpi::request > reqs;
  std::vector< NpdmSpinOps_base > nonlocal_base( world.size() );

  Timer timer;
  // Communicate basic op info
  std::vector< int > nonlocal_size( world.size() );
  std::vector< int > nonlocal_skip( world.size() );
  int local_skip = skip; // apparently boost::mpi::all_gather fails with bools...??
  boost::mpi::all_gather(world, local_skip, nonlocal_skip);
  int local_size = local_base.opReps_.size();
  boost::mpi::all_gather(world, local_size, nonlocal_size);

  // Communicate operator reps
  for (int rank = 0; rank < world.size(); ++rank) {
    if ( rank != mpigetrank() ) {
      // Get unique tag for send-recv pair (asymmetric)
      unsigned int send_tag = get_mpi_tag(mpigetrank(), rank, world.size()); assert( send_tag%100 == 0 );
      unsigned int recv_tag = get_mpi_tag(rank, mpigetrank(), world.size()); assert( send_tag%100 == 0 );
      if ( ! local_skip ) {
        std::vector< boost::mpi::request > new_reqs = local_base.isend_mpi_obj(rank, send_tag+2, send_tag+50);
        reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );
      }
      if ( ! nonlocal_skip.at(rank) ) {
        std::vector< boost::mpi::request > new_reqs = nonlocal_base.at(rank).irecv_mpi_obj(rank, recv_tag+2, recv_tag+50, nonlocal_size.at(rank));
        reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );
      }
    }
  }
  DEBUG_COMM_TIME+= timer.elapsedwalltime();
  
  // Do loop over RHS with local LHS operator while waiting for all non-local to be communicated 
  // FIXME Can we extend this idea to do batches while other batches are communicating
  if ( ! local_skip ) do_inner_loop( inner, npdm_expectations, local_base, innerOps, dotOps ); 

  // Contract all nonlocal LHS ops with local RHS ops; must wait for communication to be finished first
  Timer timer2;
  boost::mpi::wait_all( reqs.begin(), reqs.end() );
  DEBUG_COMM_TIME += timer2.elapsedwalltime();
  for (int rank = 0; rank < world.size(); ++rank) {
    if ( rank != mpigetrank() ) {
      if ( ! nonlocal_skip.at(rank) ) do_inner_loop( inner, npdm_expectations, nonlocal_base.at(rank), innerOps, dotOps ); 
    }
  }

  // Synchronize all MPI ranks here
  pout.flush();
  world.barrier();

}

void Npdm_driver::do_parallel_intermediate_loop( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                                        NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps, bool skip )
{
  boost::mpi::communicator world;
  std::map<std::vector<int>, StackWavefunction> local_waves;
  if(!skip)
  {

    if( inner =='r')
    {
      npdm_expectations.compute_intermediate(outerOps,dotOps,local_waves);
    }

    else if ( inner =='l')
    {
      npdm_expectations.compute_intermediate(outerOps,local_waves);
    }
    else assert(false);
  }



  if(outerOps.is_local_ && innerOps.is_local_  )
  {
    if(mpigetrank()==0)
    {
      if(!skip) do_inner_loop( inner, npdm_expectations, outerOps, dotOps, local_waves); 
    }
    return;
  }
  else if(outerOps.is_local_ || innerOps.is_local_ )
  {
    if(!skip) do_inner_loop( inner, npdm_expectations, outerOps, dotOps, local_waves); 
    return;
  }



  // Parallelize by broadcasting LHS or RHS intermediates
  NpdmSpinOps_base local_base(outerOps);
  std::vector< NpdmSpinOps_base > nonlocal_base( world.size() );
  std::vector< boost::mpi::request > reqs;
  std::vector<std::map<std::vector<int>, StackWavefunction>> nonlocal_waves( world.size());
  std::vector< int > nonlocal_size( world.size() );
  std::vector< int > nonlocal_skip( world.size() );

  Timer timer;
  // Communicate basic op info
  int local_skip = skip; // apparently boost::mpi::all_gather fails with bools...??
  boost::mpi::all_gather(world, local_skip, nonlocal_skip);
  int local_size = local_base.opReps_.size();
  boost::mpi::all_gather(world, local_size, nonlocal_size);


  // Communicate operator reps
  for (int rank = 0; rank < world.size(); ++rank) {
    if ( rank != mpigetrank() ) {
      // Get unique tag for send-recv pair (asymmetric)
      unsigned int send_tag = get_mpi_tag(mpigetrank(), rank, world.size()); assert( send_tag%100 == 0 );
      unsigned int recv_tag = get_mpi_tag(rank, mpigetrank(), world.size()); assert( send_tag%100 == 0 );
      if ( ! local_skip ) {
        std::vector< boost::mpi::request > new_reqs = local_base.isend_mpi_obj(rank, send_tag+2, send_tag+50);
        reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );


        boost::mpi::request new_req = world.isend(rank,mpigetrank()+1024,local_waves);
        reqs.push_back(new_req);
      }
      if ( ! nonlocal_skip.at(rank) ) {
        std::vector< boost::mpi::request > new_reqs = nonlocal_base.at(rank).irecv_mpi_obj(rank, recv_tag+2, recv_tag+50, nonlocal_size.at(rank));
        reqs.insert( reqs.end(), new_reqs.begin(), new_reqs.end() );

        boost::mpi::request new_req = world.irecv(rank,rank+1024,nonlocal_waves.at(rank));
        reqs.push_back(new_req);
      }
    }
  }
  DEBUG_COMM_TIME += timer.elapsedwalltime();
  
  // Do loop over RHS with local LHS operator while waiting for all non-local to be communicated 
  // FIXME Can we extend this idea to do batches while other batches are communicating

  if ( ! local_skip ) do_inner_loop( inner, npdm_expectations, local_base, dotOps, local_waves); 

  // Contract all nonlocal LHS ops with local RHS ops; must wait for communication to be finished first
  Timer timer2;
  boost::mpi::wait_all( reqs.begin(), reqs.end() );
  DEBUG_COMM_TIME += timer2.elapsedwalltime();
  for (int rank = 0; rank < world.size(); ++rank) {
    if ( rank != mpigetrank() ) {
      if ( ! nonlocal_skip.at(rank) ) do_inner_loop( inner, npdm_expectations, nonlocal_base.at(rank), dotOps, nonlocal_waves.at(rank)); 
    }
  }

  // Synchronize all MPI ranks here
  pout.flush();
  world.barrier();

}
#endif

  

#ifndef SERIAL
void Npdm_driver::par_loop_over_block_operators( const char inner, Npdm::Npdm_expectations & npdm_expectations,
                                                 NpdmSpinOps & outerOps, NpdmSpinOps & innerOps, NpdmSpinOps & dotOps, bool lhsrhsdot ) 
{
  //FIXME
  //Read inner loop operator and intermediate and store them in the memory to resuse them. 
  int lhs_maxsize = get_mpi_max_size( outerOps.size() );

  // Skip parallelization completely if it generates duplicates
  if ( skip_this_mpi_rank( outerOps, innerOps ) ) return;

  // Many spatial combinations on left block
  for ( int ilhs = 0; ilhs < lhs_maxsize; ++ilhs ) {
    size_t mem = Stackmem[0].memused;
    double *ptr = Stackmem[0].data+mem;

    bool skip_op = true;
    if ( ilhs < outerOps.size() ) {
      Timer timer;
      skip_op = outerOps.set_local_ops( ilhs );
      diskread_time += timer.elapsedwalltime();
    }
    
    if(dmrginp.npdm_intermediate() && (npdm_order_== NPDM_NEVPT2 || npdm_order_== NPDM_THREEPDM || npdm_order_== NPDM_FOURPDM) && dmrginp.npdm_multinode())
      do_parallel_intermediate_loop(inner, npdm_expectations, outerOps, innerOps, dotOps, skip_op );
    else{
      
      // Set local operators as dummy if load-balancing isn't perfect
      
      if ( skip_parallel( outerOps, innerOps, lhsrhsdot ) ) {
	if ( ! skip_op ) {
	  do_inner_loop( inner, npdm_expectations, outerOps, innerOps, dotOps );
	}
      }
      else {
	// Parallelize by broadcasting LHS ops
        do_parallel_lhs_loop( inner, npdm_expectations, outerOps, innerOps, dotOps, skip_op );
      }
    }

    //clear all the memory used so far
    Stackmem[0].deallocate(ptr, Stackmem[0].memused-mem);
  }
  
}
#endif
  
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
void Npdm_driver::do_inner_loop( const char inner, Npdm::Npdm_expectations& npdm_expectations, 
				 NpdmSpinOps_base& outerOps, NpdmSpinOps& innerOps, NpdmSpinOps& dotOps ) 
{
  // Many spatial combinations on right block
  if(innerOps.is_local_ && mpigetrank()>0) return;


  //SplitStackmem();
  //#pragma omp parallel for schedule(dynamic)
  for ( int iop = 0; iop < innerOps.size(); ++iop ) {
    Timer timer2;
    //bool skip = innerOps.set_local_ops( iop );
    bool skip = innerOps.set_local_ops( iop );
    diskread_time += timer2.elapsedwalltime();
    if (skip) continue;
    
    // Get non-spin-adapated spin-orbital 3PDM elements after building spin-adapted elements
    std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements;
    DEBUG_CALL_GET_EXPECT += 1;
    // This should always work out as calling in order (lhs,rhs,dot)
    if ( inner == 'r' )
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( outerOps, innerOps, dotOps );
    else if ( inner == 'l' ) 
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations( innerOps, outerOps, dotOps );
    else
      abort();
    
    Timer timer;
    // Store new npdm elements
    #pragma omp critical
    {
     if ( new_spin_orbital_elements.size() > 0 ) container_.store_npdm_elements( new_spin_orbital_elements );
     DEBUG_STORE_ELE_TIME += timer.elapsedwalltime();
    }
  }
  //MergeStackmem();

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::do_inner_loop( const char inner, Npdm::Npdm_expectations& npdm_expectations, 
                                 NpdmSpinOps_base& outerOps, NpdmSpinOps& dotOps, std::map<std::vector<int>, StackWavefunction>& outerwaves) 
{
  //if(innerOps.is_local_ && mpigetrank()>0) return;
  // Many spatial combinations on right block

  npdm_expectations.expectations_.resize(numthrds, std::vector<double>(1,0.0));
  npdm_expectations.spin_adaptation_.stored_A_mats_.resize(numthrds);
  npdm_expectations.spin_adaptation_.stored_singlet_rows_.resize(numthrds);
  npdm_expectations.spin_adaptation_.stored_so_indices_.resize(numthrds);


  //std::vector< boost::shared_ptr<NpdmSpinOps> > inneropsvector;
  //for (int i=0; i<numthrds; i++)
  //inneropsvector.push_back( innerOps.getcopy() ); 

  SplitStackmem();
#pragma omp parallel for schedule(dynamic)
  for ( int i = 0; i < inner_Operators.size(); ++i ) {
    if(inner_Operators[i] == NULL) continue;
    
    // Get non-spin-adapated spin-orbital 3PDM elements after building spin-adapted elements
    std::vector< std::pair< std::vector<int>, double > > new_spin_orbital_elements;
    DEBUG_CALL_GET_EXPECT += 1;
    // This should always work out as calling in order (lhs,rhs,dot)
    if ( inner == 'r' )
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations(outerOps, *inner_Operators[i], dotOps, outerwaves, *inner_intermediate[i] );
    else if ( inner == 'l' ) 
      new_spin_orbital_elements = npdm_expectations.get_nonspin_adapted_expectations(*inner_Operators[i], outerOps, dotOps, *inner_intermediate[i], outerwaves );
    else
      abort();
    
    Timer timer;
#pragma omp critical
    {
      // Store new npdm elements
      if ( new_spin_orbital_elements.size() > 0 ) container_.store_npdm_elements( new_spin_orbital_elements );
      DEBUG_STORE_ELE_TIME += timer.elapsedwalltime();
    }
  }

  MergeStackmem();

}



//-----------------------------------------------------------------------------------------------------------------------------------------------------------


bool Npdm_driver::screen(const std::vector<CD> &lhs_cd_type,const std::vector<CD> &dot_cd_type)
{
  int cre_num=0;
  int des_num=0;
  for(auto i:lhs_cd_type)
  {
    if(i== CREATION) cre_num++;
    else if (i== DESTRUCTION) des_num++;
  }

  for(auto i:dot_cd_type)
  {
    if(i== CREATION) cre_num++;
    else if (i== DESTRUCTION) des_num++;
  }
  if(cre_num> npdm_order_) return true;
  if(des_num> npdm_order_) return true;
  return false;

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void Npdm_driver::clear_npdm_intermediate(Npdm::Npdm_expectations& expectations)
{
  for(std::string filename: expectations.intermediate_filenames)
    remove(filename.c_str());
  expectations.intermediate_filenames.clear();

}

}
}
