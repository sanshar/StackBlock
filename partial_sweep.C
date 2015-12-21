/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include "blas_calls.h"
#include "IntegralMatrix.h"
#include "stackguess_wavefunction.h"
#include "sweep.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "rotationmat.h"
#include "Stackdensity.h"
#include "pario.h"
#include "Stackwavefunction.h"
#include "Stackspinblock.h"
#include "sweep_params.h"
#include "operatorfunctions.h"

using namespace boost;
using namespace std;



double SpinAdapted::Sweep::do_one_partial(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize)
{
  Timer sweeptimer;
  int integralIndex = 0; //By default we assume that we only have one set of integrals and its index is 0
  StackSpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  std::vector<double> finalEnergy(nroots,1.0e10);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;
  if (restart) {
    finalEnergy = sweepParams.get_lowest_energy();
    finalEnergy_spins = sweepParams.get_lowest_energy();
    finalError = sweepParams.get_lowest_error();
  }

  sweepParams.set_sweep_parameters();
  if (dmrginp.spinAdapted()) {
    sweepParams.set_n_iters() =  sweepParams.set_n_iters()+(dmrginp.getPartialSweep() - dmrginp.last_site())/sweepParams.get_sys_add();
  }
  else {
    sweepParams.set_n_iters() =  sweepParams.set_n_iters()+(dmrginp.getPartialSweep() - dmrginp.last_site())/(2*sweepParams.get_sys_add());
  }
  sweepParams.set_backward_starting_size() = dmrginp.last_site()-dmrginp.getPartialSweep()+1;

  // a new renormalisation sweep routine
  pout << endl;
  if (forward) {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
  }
  else {
    pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
  }
  pout << "\t\t\t ============================================================================ " << endl;

  int len = forward? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();
  vector<int> sites(len);
  if (forward)
    for (int i=0; i<len; i++)
      sites[i] = i;
  else
    for (int i=0; i<len; i++)
      sites[i] = dmrginp.last_site() - len + i;

  StackSpinBlock::restore (forward, sites, system, sweepParams.current_root(), sweepParams.current_root());
  system.set_twoInt(integralIndex);

  if(!restart)
    sweepParams.set_block_iter() = 0;


  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  StackSpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root()); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (restart || dmrginp.get_sweep_type() == PARTIAL)
    set_dot_with_sys(dot_with_sys, system, sweepParams, forward);

  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  bool useRGStartUp = false;

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ============================" << endl;
      if (forward) {
       p1out << "\t\t\t Current direction is :: Forwards " << endl;
      }
      else {
       p1out << "\t\t\t Current direction is :: Backwards " << endl;
      }

      if (sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if (sweepParams.get_block_iter() == 0 && 
                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;

      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp && (dmrginp.warmup() == WILSON || (sym=="dinfh" || NonabelianSym || dmrginp.hamiltonian()==HEISENBERG))) {
	useRGStartUp = true;
	Startup(sweepParams, system, newSystem);
      }
      else {
         if (sweepParams.set_sweep_iter() == 1 && sweepParams.get_block_iter() == 0)
           sweepParams.set_guesstype() = BASIC;
         if(sweepParams.set_sweep_iter() == 1 && sweepParams.get_largest_dw()<=NUMERICAL_ZERO)
           sweepParams.set_additional_noise() = dmrginp.get_twodot_noise();
         BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys);
      }

      //mcheck("In sweep after block and decimate");

      //Need to substitute by?

      //if (!(warmUp && (sym=="trans" || sym == "dinfh_abelian" || NonabelianSym || dmrginp.hamiltonian()==HEISENBERG))){
      if (!useRGStartUp) {
	for(int j=0;j<nroots;++j)
	{
	  int istate = dmrginp.setStateSpecific() ? sweepParams.current_root() : j;
	  
#ifndef MOLPRO
	  pout << "\t\t\t Total block energy for State [ " << istate << 
	    " ] with " << sweepParams.get_keep_states()<<" States :: " << fixed << sweepParams.get_lowest_energy()[j] <<endl;              
#else 
	  //We might want to relax the output restrictions here, so it prints out with outputlevel=0
          p1out << "\t\t\t Total block energy for State [ " << istate << 
	      " ] with " << sweepParams.get_keep_states()<<" States :: " << fixed << setprecision(10) << sweepParams.get_lowest_energy()[j] <<endl;              
#endif
	}
	
	//this criteria should work for state average or state specific because the lowest sweep energy is always the lowest of the average
	finalEnergy_spins = ( (std::accumulate(sweepParams.get_lowest_energy().begin(), sweepParams.get_lowest_energy().end(),0.0) < std::accumulate(finalEnergy.begin(), finalEnergy.end(),0.0)) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
	finalEnergy = ((std::accumulate(sweepParams.get_lowest_energy().begin(), sweepParams.get_lowest_energy().end(),0.0) < std::accumulate(finalEnergy.begin(), finalEnergy.end(),0.0)) ? sweepParams.get_lowest_energy() : finalEnergy);
	finalError = max(sweepParams.get_lowest_error(),finalError);
      }
      
      system = newSystem;
      system.printOperatorSummary();

      

      StackSpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root());	 	
      pout << system<<endl;
      //if (sweepParams.set_block_iter() == 4) exit(0);
      set_dot_with_sys(dot_with_sys, system, sweepParams, forward);

      p1out << "\t\t\t Saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      
#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(calc,finalError,0);
      calc.barrier();
#endif
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }

  system.deallocate();
  system.clear();

  for(int j=0;j<nroots;++j) {
    int istate = dmrginp.setStateSpecific() ? sweepParams.current_root() : j;
    pout << "\n\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << istate 
	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j] << endl;
  }

  pout << "\n\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  sweepParams.set_largest_dw() = finalError;

  for(int j=0;j<nroots;++j){
    int istate = dmrginp.setStateSpecific() ? sweepParams.current_root() : j;
#ifndef MOLPRO
//  printf("\t\t\t M = %6i  state = %4i  Largest Discarded Weight = %8.3e  Sweep Energy = %20.10f \n",sweepParams.get_keep_states(), istate, finalError, finalEnergy[j]+dmrginp.get_coreenergy());
    pout << "\t\t\t M = " << setw(6) << sweepParams.get_keep_states()
         << "  state = " << setw(4) << istate
         << "  Largest Discarded Weight = " << finalError
         << "  Sweep Energy = " << finalEnergy[j]
         << " " << endl;
#else 
    //printf("\t\t\t M = %6i   Largest Discarded Weight = %8.3e  Sweep Energy = %20.10f \n",sweepParams.get_keep_states(), finalError, finalEnergy[j]+dmrginp.get_coreenergy());
    pout << "\t\t\t M = " <<  setw(6) << sweepParams.get_keep_states() ; 
    pout << "\t Largest Discarded Weight = " << scientific << setprecision(3) << finalError ;
    pout << "\t Sweep Energy = " << fixed << setprecision(10) << finalEnergy[j] << endl;
#endif
  }
  pout << "\t\t\t ============================================================================ " << endl;
  double cputime = sweeptimer.elapsedcputime();
  double walltime = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << cputime << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << walltime<< endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();
  //if (!(warmUp && (sym=="trans" || sym == "dinfh_abelian" || NonabelianSym || dmrginp.hamiltonian()==HEISENBERG))){
  if (!useRGStartUp) {
    if (!mpigetrank())
    {
      std::string efile;
      efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );

      //if state specific only write back the current energy to the dmrg.e file and leave the rest unchanged
      if (dmrginp.setStateSpecific() ) {
	sweepParams.set_lowest_energy().resize(dmrginp.nroots());
	sweepParams.set_lowest_energy()[sweepParams.current_root()] = sweepParams.get_lowest_energy()[0];
	FILE* fin = fopen(efile.c_str(), "rb");
	for(int j=0;j<dmrginp.nroots();++j) {
	  double e;
	  fread( &e, 1, sizeof(double), fin);
	  if (j != sweepParams.current_root())
	    sweepParams.set_lowest_energy()[j] = e;
	}
	fclose(fin);
      }

      FILE* f = fopen(efile.c_str(), "wb");      
      for(int j=0;j<dmrginp.nroots();++j) {
	double e = sweepParams.get_lowest_energy()[j]; //instead of the lowest energy of the sweep, we record the last energy of the sweep
	fwrite( &e, 1, sizeof(double), f);
      }
      fclose(f);
    }
  }

  return std::accumulate(finalEnergy.begin(), finalEnergy.end(),0.0)/dmrginp.nroots(sweepParams.get_sweep_iter());
}

