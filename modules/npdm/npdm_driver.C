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
double DEBUG_COMM_TIME;
int DEBUG_CALL_GET_EXPECT;
double DEBUG_STORE_ELE_TIME;

// Forward declaration
boost::shared_ptr<NpdmSpinOps> select_op_wrapper( StackSpinBlock * spinBlock,const std::vector<Npdm::CD> & cd_type );


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
  void Npdm_driver::get_inner_Operators( const char inner, Npdm_expectations& npdm_expectations, boost::shared_ptr<NpdmSpinOps> lhsOps, boost::shared_ptr<NpdmSpinOps> dotOps, boost::shared_ptr<NpdmSpinOps> rhsOps, int procrank) 
{
  // Many spatial combinations on right block
  if( inner == 'l')
    {
      cout << "the inner loop should always be the RHS"<<endl;
      exit(0);
    }
  else if( inner == 'r')
    {
      int rhsopsize = rhsOps->size();
#ifndef SERIAL
      mpi::broadcast(calc, rhsopsize, procrank);
#endif

      for (int i=0; i<rhsopsize; i++)
	inner_Operators.push_back(rhsOps->getcopy());

      inner_intermediate.resize(rhsopsize);

      if (mpigetrank() == procrank) {
	//allocate memory for the inner_intermediates
	for ( int i = 0; i < rhsopsize; ++i ) {
	  bool skip = rhsOps->set_local_ops( i );
	  if (!skip) {
	    //boost::shared_ptr<NpdmSpinOps> newOps( new NpdmSpinOps(*rhsOps));
	    inner_Operators[i]->set_local_ops(i);
	    inner_intermediate[i] = boost::shared_ptr<std::map<std::vector<int>, StackWavefunction> >(new std::map<std::vector<int>, StackWavefunction>);
	    npdm_expectations.AllocateInitialiseWavefunctions(*inner_Operators[i], *inner_intermediate[i]);
	  }
	  else {
	    inner_Operators[i] = boost::shared_ptr<NpdmSpinOps>();
	    inner_intermediate[i] = boost::shared_ptr<std::map<std::vector<int>, StackWavefunction> >();
	  }
	}
	
	std::vector<boost::shared_ptr<NpdmSpinOps> > rhsopsvec;
	for (int i=0; i<numthrds; i++)
	  rhsopsvec.push_back(rhsOps->getcopy());

	SplitStackmem();      
#pragma omp parallel for schedule(dynamic)
	for ( int i = 0; i < rhsOps->size(); ++i ) {
	  if(inner_Operators[i] == NULL)
	    inner_intermediate[i] = boost::shared_ptr<std::map<std::vector<int>, StackWavefunction> >();
	  else{
	    npdm_expectations.compute_intermediate( *inner_Operators[i], *inner_intermediate[i]);
	  }
	}
	MergeStackmem();
      }
      
#ifndef SERIAL
      //In this messy code we take the inner_Operators and inner_intermediates from the processor "procrank"
      //and broadcast it to all procs. It is more messy than it needs to be because a lot of data is stored 
      //in boost pointers. To broadcast boost pointers we have to first allocate the memory in boost pointers.
      // So a big part of the following code is devoted to checking the length of vectors containing boost:pointers,
      //making the vectors and assining memory to each vector element and finally broadcasting the data.
      for (int i=0; i<rhsopsize; i++) {
	int isNull = 1;
	if (mpigetrank() == procrank) 
	  isNull = inner_Operators[i] == NULL ? 1 : 0;
	mpi::broadcast(calc, isNull, procrank);
	if (isNull==1) {
	  if (mpigetrank() != procrank) {
	    inner_Operators[i] = boost::shared_ptr<NpdmSpinOps>();
	    inner_intermediate[i] = boost::shared_ptr<std::map<std::vector<int>, StackWavefunction> >();
	  }
	}
	else {
	  //for inner_operators, we dont need to transmit the data, just the shell
	  mpi::broadcast(calc, inner_Operators[i]->build_pattern_ , procrank); 
	  mpi::broadcast(calc, inner_Operators[i]->transpose_ , procrank); 
	  mpi::broadcast(calc, inner_Operators[i]->factor_ , procrank); 
	  mpi::broadcast(calc, inner_Operators[i]->indices_ , procrank); 
	  mpi::broadcast(calc, inner_Operators[i]->is_local_ , procrank); 
	  mpi::broadcast(calc, inner_Operators[i]->opIndices_ , procrank); 

	  //initialize the vector,Wavefunction map
	  if (mpigetrank() != procrank) {
	    inner_intermediate[i] = boost::shared_ptr<std::map<std::vector<int>, StackWavefunction> >(new std::map<std::vector<int>, StackWavefunction>);
	  } 
	  mpi::broadcast(calc, *inner_intermediate[i] , procrank);

	  int index = 0;
	  for (std::map<std::vector<int>, StackWavefunction>::iterator it = inner_intermediate[i]->begin(); it != inner_intermediate[i]->end(); it++){
	    //allocate the memory for the wavefunction
	    if (mpigetrank() != procrank) {
	      it->second.initialise(it->second.get_deltaQuantum(), npdm_expectations.big_.get_leftBlock()->get_ketStateInfo(), npdm_expectations.big_.get_rightBlock()->get_braStateInfo(),true);
	      inner_Operators[i]->opReps_.push_back(boost::shared_ptr<StackSparseMatrix>(new StackSparseMatrix));
	    }

	    mpi::broadcast(calc, *inner_Operators[i]->opReps_[index] , procrank);
	    //broadcast the data
	    MPI_Bcast(it->second.get_data(), it->second.memoryUsed(), MPI_DOUBLE, procrank, Calc);
	    index++;
	  }
	}
      }
#endif

    }
  
}
  
  //-----------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::do_inner_loop( const char inner, Npdm::Npdm_expectations& npdm_expectations, 
                                 NpdmSpinOps_base& outerOps, NpdmSpinOps& dotOps, std::map<std::vector<int>, StackWavefunction>& outerwaves, int i) 
{
  //npdm_expectations.expectations_.resize(numthrds, std::vector<double>(1,0.0));
  //npdm_expectations.spin_adaptation_.stored_A_mats_.resize(numthrds);
  //npdm_expectations.spin_adaptation_.stored_singlet_rows_.resize(numthrds);
  //npdm_expectations.spin_adaptation_.stored_so_indices_.resize(numthrds);

  if(inner_Operators[i] == NULL) return;
    
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



//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::loop_over_block_operators( Npdm::Npdm_expectations& npdm_expectations, NpdmSpinOps& outerOps, NpdmSpinOps& innerOps, NpdmSpinOps& dotOps)
{
  npdm_expectations.expectations_.resize(numthrds, std::vector<double>(1,0.0));
  npdm_expectations.spin_adaptation_.stored_A_mats_.resize(numthrds);
  npdm_expectations.spin_adaptation_.stored_singlet_rows_.resize(numthrds);
  npdm_expectations.spin_adaptation_.stored_so_indices_.resize(numthrds);

  std::vector< boost::shared_ptr<NpdmSpinOps> > outeropsvector;
  std::vector< boost::shared_ptr<NpdmSpinOps> > dotopsvector;
  for (int i=0; i<numthrds; i++) {
    outeropsvector.push_back( outerOps.getcopy() ); 
    dotopsvector.push_back( dotOps.getcopy() ); 
    dotopsvector[i]->set_local_ops(0);
  }

  std::vector<bool> skip_op(numthrds, false);
  //This is so messy because we want only one thread to read from the disk
  //this reading happens in the function set_local_ops
  //I dont know why, but when I let all threads to read their own operators, the program does not work
  if (!dmrginp.do_npdm_in_core() && outerOps.opIndices_ > 2) {
    SplitStackmem();
#pragma omp parallel
    {
      for ( int ilhs = 0; ilhs < outerOps.size()/numthrds+1; ++ilhs ) {
	
	size_t mem = Stackmem[omprank].memused;
	double *ptr = Stackmem[omprank].data+mem;
	
	if (omprank == 0) {
	  for (int i=0; i<numthrds; i++) {
	    if (ilhs*numthrds+i < outerOps.size())
	      skip_op[i] = outeropsvector[i]->set_local_ops(ilhs*numthrds+i);
	  }
	}
#pragma omp barrier
	
	
	if (!skip_op[omprank] && ilhs*numthrds+omprank<outerOps.size()) {
	  std::map<std::vector<int>, StackWavefunction> outerWaves;
	  npdm_expectations.compute_intermediate(*outeropsvector[omprank], *dotopsvector[omprank], outerWaves);
	  for (int k=0; k<inner_intermediate.size(); k++) 
	    do_inner_loop( 'r', npdm_expectations, *outeropsvector[omprank], *dotopsvector[omprank], outerWaves, k);      
	  outerWaves.clear();
	}
	
#pragma omp barrier
	if (omprank == 0) {
	  for (int i=0; i<numthrds; i++) {
	    if (ilhs*numthrds+i < outerOps.size())
	      for (int j=0; j<outeropsvector[i]->opReps_.size(); j++) 
		outeropsvector[i]->opReps_[j]->CleanUp();
	  }
	}

	Stackmem[omprank].deallocate(ptr, Stackmem[omprank].memused-mem);
	
      }    
    }
    MergeStackmem();
  }
  else {

    SplitStackmem();
    // Many spatial combinations on left block
#pragma omp parallel for schedule(dynamic)
    for ( int ilhs = 0; ilhs < outerOps.size(); ++ilhs ) {
      // Set local operators as dummy if load-balancing isn't perfect
      bool skip_op = true;
      //Timer timer2;
      
      size_t mem = Stackmem[omprank].memused;
      double *ptr = Stackmem[omprank].data+mem;
      
      skip_op = outeropsvector[omprank]->set_local_ops( ilhs );
      //diskread_time += timer2.elapsedwalltime();

      if (!skip_op) {
	std::map<std::vector<int>, StackWavefunction> outerWaves;
	npdm_expectations.compute_intermediate(*outeropsvector[omprank], *dotopsvector[omprank], outerWaves);
	for (int k=0; k<inner_intermediate.size(); k++) 
	  do_inner_loop( 'r', npdm_expectations, *outeropsvector[omprank], *dotopsvector[omprank], outerWaves, k);      
	outerWaves.clear();
      }
      
      Stackmem[omprank].deallocate(ptr, Stackmem[omprank].memused-mem);
    }
    MergeStackmem();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::loop_over_operator_patterns( Npdm::Npdm_patterns& patterns, Npdm::Npdm_expectations& expectations, const StackSpinBlock& big )
{
#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  // Get LHS, Dot and RHS spin-blocks
  StackSpinBlock* rhsBlock = big.get_rightBlock();
  StackSpinBlock* lhsdotBlock = big.get_leftBlock();
  StackSpinBlock* lhsBlock = lhsdotBlock->get_leftBlock();
  StackSpinBlock* dotBlock = lhsdotBlock->get_rightBlock();

  //std::cout << *big.get_rightBlock()->get_op_rep(DES_DES, SpinQuantum(-2,SpinSpace(0),SpinAdapted::IrrepSpace(0)), 0,0)<<endl;

  int count = 0;
  for (auto pattern = patterns.ldr_cd_begin(); pattern != patterns.ldr_cd_end(); ++pattern) {
    count++;
    DEBUG_CALL_GET_EXPECT= 0;

    //if (count == 5) exit(0);
    //if (count == 4) {
    //pout << "at 4 "<<endl;
    //}

#ifndef SERIAL
    // MPI threads must be synchronised here so they all work on same operator pattern simultaneously
    pout.flush();
    calc.barrier();
#endif
    //pout << "-------------------------------------------------------------------------------------------\n";
    pout << "Doing pattern " << count << " of " << patterns.size() <<"   ";
    patterns.print_cd_string( pattern->at('l') );
    patterns.print_cd_string( pattern->at('d') );
    patterns.print_cd_string( pattern->at('r') );
    pout << std::endl; 
    pout.flush();

    // Choice of read from disk or not done inside the wrapper
    std::vector<Npdm::CD> lhs_cd_type = pattern->at('l');
    std::vector<Npdm::CD> dot_cd_type = pattern->at('d');
    std::vector<Npdm::CD> rhs_cd_type = pattern->at('r');


    boost::shared_ptr<NpdmSpinOps> rhsOps = select_op_wrapper( rhsBlock, rhs_cd_type );
    boost::shared_ptr<NpdmSpinOps> dotOps = select_op_wrapper( dotBlock, dot_cd_type );
    boost::shared_ptr<NpdmSpinOps> lhsOps = select_op_wrapper( lhsBlock, lhs_cd_type );


    // Only one spatial combination on the dot block (including NULL)
    if(dmrginp.spinAdapted()){
      assert( dotOps->size() == 1 );

      for(int i=0; i< dotOps->size(); i++){
      
	bool skip = dotOps->set_local_ops( i );
	if ( ! skip ) {
	  // <Psi_1| L_i d_j  R_k |Psi_0> 
	  //we form a series of  |Psi_k> = R_k |Psi_0>
	  
	  //then we loop over ij and contract <Psi_1| L_i d_j |Psi_k> 
	  //these expectation values are then converted to usable RDM elements
	  Timer timer2;
	  inner_Operators.clear();
	  inner_intermediate.clear();
	  
#ifndef SERIAL
	  for (int proc=0; proc<calc.size(); proc++) {
#else
	  {
	    int proc = 0;
#endif
	    if (rhsOps->is_local_ && proc > 0) continue;
	    size_t mem = Stackmem[0].memused;
	    double *ptr = Stackmem[0].data+mem;
	    
	    get_inner_Operators( 'r', expectations, lhsOps, dotOps , rhsOps, proc) ;  // this makes the |Psi_k> and stores them as intermediates
	    diskread_time += timer2.elapsedwalltime();
	    loop_over_block_operators( expectations, *lhsOps, *rhsOps, *dotOps ); //contracts <Psi_1|L_i d_j|Psi_k>
	    
	    inner_Operators.clear();
	    inner_intermediate.clear();
	    Stackmem[0].deallocate(ptr, Stackmem[0].memused-mem);
	  }
	}
      }
    }
    else{
      //bool skip=true; 
      // build all kind of dotOps
      for(int i=0; i< dotOps->size(); i++){
        // if it is valid
        if(!dotOps->set_local_ops( i )){

//pout << "p" << mpigetrank() << ": lhs = " << lhsOps->size() << endl;
//pout << "p" << mpigetrank() << ": rhs = " << rhsOps->size() << endl;
      // Compute all irreducible PDM elements generated by this block operator pattern at this sweep position
#ifndef SERIAL
      bool lhs_or_rhs_dot = ( (lhsBlock->size() == 1) || (rhsBlock->size() == 1) );
//pout << "lhs_or_rhs_dot " << lhs_or_rhs_dot << endl;
      if ( broadcast_lhs( lhsOps->size(), rhsOps->size() ) ) {
//pout << "broadcast lhs\n";
//pout.flush();
	
      par_loop_over_block_operators( 'r', expectations, *lhsOps, *rhsOps, *dotOps, lhs_or_rhs_dot );
    }
    else {
//pout << "broadcast rhs\n";
//pout.flush();
        par_loop_over_block_operators( 'l', expectations, *rhsOps, *lhsOps, *dotOps, lhs_or_rhs_dot );
      }
#else
      loop_over_block_operators( expectations, *lhsOps, *rhsOps, *dotOps );
#endif
      }
      }
    }


  }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Npdm_driver::compute_npdm_elements(std::vector<StackWavefunction> & wavefunctions, const StackSpinBlock & big, int sweepPos, int endPos)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  pout.flush();
  calc.barrier();
#endif
  DEBUG_COMM_TIME = 0;
  DEBUG_STORE_ELE_TIME = 0;
  DEBUG_CALL_GET_EXPECT = 0;
	double total_time =0;
  Timer timer;
  pout << "===========================================================================================\n";
  pout << "Current NPDM sweep position = "<< sweepPos+1 << " of " << endPos+1 << "\n";

  pout << *big.get_leftBlock() <<endl;
  big.get_leftBlock()->printOperatorSummary();
  //TODO
  //Store intermidiate of O_l|\Psi> and <\Psi|O_r
  
  // Loop over NPDM operator patterns
  Npdm_patterns npdm_patterns( npdm_order_, sweepPos, endPos );
  StackWavefunction& wave1= wavefunctions.size()==2? wavefunctions.at(1): wavefunctions.at(0);
  Npdm_expectations npdm_expectations( spin_adaptation_, npdm_patterns, npdm_order_, wavefunctions.at(0), wave1, big );

  loop_over_operator_patterns( npdm_patterns, npdm_expectations, big );
#ifndef SERIAL
  calc.barrier();
#endif
  if(dmrginp.npdm_intermediate() && (npdm_order_== NPDM_NEVPT2 || npdm_order_== NPDM_THREEPDM || npdm_order_== NPDM_FOURPDM))
    clear_npdm_intermediate(npdm_expectations);

	diskread_time +=npdm_expectations.diskread_time; 
  // Print outs
#ifndef SERIAL
  if (mpigetrank() == 0) {
    int sum;
    reduce(calc, DEBUG_CALL_GET_EXPECT, sum, std::plus<int>(), 0);
    p3out << "NPDM calls to expectation engine " << sum << endl;
  } else {
    reduce(calc, DEBUG_CALL_GET_EXPECT, std::plus<int>(), 0);
  }

  if (mpigetrank() == 0) {
    double sum;
    reduce(calc, DEBUG_COMM_TIME, sum, std::plus<double>(), 0);
    p3out << "NPDM mpi communications time " << sum << endl;
  } else {
    reduce(calc, DEBUG_COMM_TIME, std::plus<double>(), 0);
  }

  if (mpigetrank() == 0) {
    double sum;
    reduce(calc, DEBUG_STORE_ELE_TIME, sum, std::plus<double>(), 0);
    p3out << "NPDM store elements time " << sum << endl;
  } else {
    reduce(calc, DEBUG_STORE_ELE_TIME, std::plus<double>(), 0);
  }

  if (mpigetrank() == 0) {
    double sum;
    reduce(calc, diskread_time, sum, std::plus<double>(), 0);
    p3out << "NPDM operators reading time " << sum << endl;
  } else {
    reduce(calc, diskread_time, std::plus<double>(), 0);
  }

  double ecpu = timer.elapsedcputime();double ewall=timer.elapsedwalltime();
  p3out << "NPDM compute elements time " << ewall << " "<< ecpu << endl;
#else
  p3out << "NPDM compute elements time " << timer.elapsedwalltime() << " " << timer.elapsedcputime() << endl;
#endif
  pout << "===========================================================================================\n";
  
}

//===========================================================================================================================================================

}
}
