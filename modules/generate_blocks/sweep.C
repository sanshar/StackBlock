/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "Stackspinblock.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "rotationmat.h"
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "operatorfunctions.h"
#include "sweepgenblock.h"
#include "stackguess_wavefunction.h"
#include "davidson.h"
#include "pario.h"
#ifdef USE_BTAS
#include "overlaptensor.h"
#endif
#include "sweep.h"
using namespace boost;
using namespace std;

namespace SpinAdapted{
void SweepGenblock::BlockAndDecimate (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int stateA, int stateB)
{
  if (dmrginp.outputlevel() > 0) 
    mcheck("at the start of block and decimate");
  p1out << "\t\t\t Performing Blocking"<<endl;
  dmrginp.guessgenT -> start();
  // figure out if we are going forward or backwards  
  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot;
  int systemDotStart, systemDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
  }
  vector<int> spindotsites(2); 
  spindotsites[0] = systemDotStart;
  spindotsites[1] = systemDotEnd;

  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  dmrginp.guessgenT -> stop();

  bool doNorms = (dot_with_sys || dmrginp.new_npdm_code() == true);
  bool doComp = dmrginp.get_lowMemoryAlgorithm() ? !dot_with_sys : true;
  //bool doComp = true;
  if (doComp && !system.has(CRE_DESCOMP) && !dmrginp.npdm_generate())
    system.addAllCompOps();
  if (dmrginp.npdm_generate()) {
    dmrginp.datatransfer->start();
    system.addOneIndexNormOps();    
    dmrginp.datatransfer->stop();
  } else {
    system.addAdditionalOps();
  }

  //bool doComp = true;
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, stateA, stateB, sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, doNorms, doComp);


  pout << "\t\t\t System  Block"<<newSystem;
  newSystem.printOperatorSummary();

  std::vector<Matrix> leftrotateMatrix, rightrotateMatrix;

  LoadRotationMatrix (newSystem.get_sites(), leftrotateMatrix, stateA);
  LoadRotationMatrix (newSystem.get_sites(), rightrotateMatrix, stateB);

    
#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, leftrotateMatrix, 0);
  broadcast(calc, rightrotateMatrix, 0);
#endif

  p1out <<"\t\t\t Performing Renormalization "<<endl<<endl;

  dmrginp.operrotT->start();
  if (stateB == stateA && dmrginp.setimplicitTranspose()==true)
    newSystem.transform_operators(leftrotateMatrix);
  else
    newSystem.transform_operators(leftrotateMatrix, rightrotateMatrix);
  dmrginp.operrotT->stop();

  //if (system.get_sites().size() != 1) {
  //if (system.get_sites().size() != 1 || (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted())) {
  {
    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysmem = newSystem.memoryUsed();
    newSystem.moveToNewMemory(system.getdata());
    Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
    system.clear();
  }

  p2out << str(boost::format("%-40s - %-10.4f\n") % "Total walltime" % globaltimer.totalwalltime());
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  -->Blocking" % *(dmrginp.guessgenT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  -->Wavefunction Solution" % *(dmrginp.davidsonT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  -->Add noise" % *(dmrginp.rotmatrixT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  -->Renormalisation" % *(dmrginp.operrotT));
  //p2out << str(boost::format("%-40s - %-10.4f\n") % "diskio" % *(dmrginp.diskio)); 
  p2out << str(boost::format("%-40s - %-10.4f\n") % "mpicommunication" % *(dmrginp.datatransfer)); 

}

double SweepGenblock::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int stateA, int stateB)
{
  Timer sweeptimer;
  int integralIndex = 0;

  StackSpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, stateA, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;

  p2out << "\t\t\t Starting block is :: " << endl << system << endl;
  //if (!restart) 
  StackSpinBlock::store (forward, system.get_sites(), system, stateA, stateB); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  if (restart)
    Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);
  //Sweep::set_dot_with_sys(dot_with_sys, system, const_cast<bool&>(forward));

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if (!warmUp && sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if (!warmUp && sweepParams.get_block_iter() == 0 && 
                ((dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter() != sweepParams.get_sweep_iter()) ||
                  dmrginp.algorithm_method() != TWODOT_TO_ONEDOT))
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, warmUp, dot_with_sys, stateA, stateB);

      
      system = newSystem;

      Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);

      StackSpinBlock::store (forward, system.get_sites(), system, stateA, stateB);	 	
      pout << system<<endl;
      p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
      ++sweepParams.set_block_iter();
      //if (sweepParams.get_onedot())
      //pout << "\t\t\tUsing one dot algorithm!!"<<endl; 
      sweepParams.savestate(forward, system.get_sites().size());
    }

  system.deallocate();
  system.clear();
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();
  double cputime = sweeptimer.elapsedcputime();
  double walltime = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << cputime << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << walltime << endl;

  return finalEnergy[0];
}


void SweepGenblock::do_one(SweepParams &sweepParams, const bool &forward, int stateA, int stateB)
{
  Timer sweeptimer;
  int integralIndex = 0;
  StackSpinBlock system;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  InitBlocks::InitStartingBlock (system,forward, stateA, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 0, false, false, integralIndex);

  sweepParams.set_block_iter() = 0;

  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  bool dot_with_sys = true;

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }

  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if ( sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if ( sweepParams.get_block_iter() == 0 )
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, false, dot_with_sys, stateA, stateB);

      system = newSystem;

      StackSpinBlock::store(forward, system.get_sites(), system, stateA, stateB);

      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	dot_with_sys = false;

      ++sweepParams.set_block_iter();
    }
  system.deallocate();
  system.clear();
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  double cputime = sweeptimer.elapsedcputime();
  double walltime = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << cputime << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << walltime << endl;

}


void SweepGenblock::do_one_partialSweep(SweepParams &sweepParams, const bool &forward, int stateA, int stateB, int integralIndex)
{
  Timer sweeptimer;
  StackSpinBlock system;

  sweepParams.set_sweep_parameters();
  if (dmrginp.get_sweep_type() == PARTIAL) {
    if (dmrginp.spinAdapted()) {
      sweepParams.set_n_iters() =  sweepParams.set_n_iters()+(dmrginp.getPartialSweep() - dmrginp.last_site())/sweepParams.get_sys_add();
    }
    else {
      sweepParams.set_n_iters() =  sweepParams.set_n_iters()+(dmrginp.getPartialSweep() - dmrginp.last_site())/(2*sweepParams.get_sys_add());
    }
  
    sweepParams.set_backward_starting_size() = dmrginp.last_site()-dmrginp.getPartialSweep()+1;
  }

  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  if (dmrginp.get_sweep_type() == PARTIAL && !forward) {
    int len = forward? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();
    std::vector<int> sites(len);  
    if (forward)
      for (int i=0; i<len; i++)
	sites[i] = i;
    else
      for (int i=0; i<len; i++)
	sites[i] = dmrginp.last_site() - len + i;

    StackSpinBlock::restore (forward, sites, system, stateA, stateB);
  }
  else
    InitBlocks::InitStartingBlock (system,forward, stateA, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 0, false, false, integralIndex);

  sweepParams.set_block_iter() = 0;

  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  bool dot_with_sys = true;
  StackSpinBlock::store (forward, system.get_sites(), system, stateA, stateB); 

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }

  
      if (dmrginp.no_transform())
	      sweepParams.set_guesstype() = BASIC;
      else if ( sweepParams.get_block_iter() != 0) 
  	    sweepParams.set_guesstype() = TRANSFORM;
      else if ( sweepParams.get_block_iter() == 0 )
        sweepParams.set_guesstype() = TRANSPOSE;
      else
        sweepParams.set_guesstype() = BASIC;
      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem;

      BlockAndDecimate (sweepParams, system, newSystem, false, dot_with_sys, stateA, stateB);

      system = newSystem;

      StackSpinBlock::store(forward, system.get_sites(), system, stateA, stateB);

      Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);

      ++sweepParams.set_block_iter();
    }
  system.deallocate();
  system.clear();
  pout << "\t\t\t Finished Generate-Blocks Sweep. " << endl;
  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  double cputime = sweeptimer.elapsedcputime();
  double walltime = sweeptimer.elapsedwalltime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << cputime << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << walltime << endl;

}


}
