/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <newmat.h>
#include <newmatio.h>
#include "npdm.h"
#include "sweep.h"
#include "sweepgenblock.h"
#include "Stackdensity.h"
//#include "sweeponepdm.h"  // For legacy version of 1pdm
//#include "sweeptwopdm.h"  // For legacy version of 2pdm
#include "npdm_driver.h"
#include "nevpt2_npdm_driver.h"
#include "pario.h"
#include "initblocks.h"
#include "stackguess_wavefunction.h"
#include "Stackwavefunction.h"
#include "SpinQuantum.h"

void dmrg(double sweep_tol);
void restart(double sweep_tol, bool reset_iter);
void ReadInput(char* conf);
void fullrestartGenblock();

namespace SpinAdapted {
namespace Npdm {

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm_block_and_decimate( Npdm_driver_base& npdm_driver, SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, 
                              const bool &useSlater, const bool& dot_with_sys, const int state, const int stateB)
{
  Timer timer;
  //mcheck("at the start of block and decimate");
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot;
  StackSpinBlock envDot;
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
  //if (useSlater) {
    systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
    //StackSpinBlock::store(true, systemDot.get_sites(), systemDot);
    //}
    //else
    //StackSpinBlock::restore(true, spindotsites, systemDot);
  StackSpinBlock environment, environmentDot, newEnvironment;

  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  system.addOneIndexNormOps();
  if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep()){
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, state, stateB,
                                   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, true, false);
    
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      state, stateB,
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, environment.get_integralIndex(), true, false, true);
  }
  else{
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, sweepParams.current_root(), sweepParams.current_root(),
                                   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, true, false);
    
    InitBlocks::InitNewEnvironmentBlock(environment, systemDot, newEnvironment, system, systemDot,
                                      sweepParams.current_root(), sweepParams.current_root(),
                                      sweepParams.get_sys_add(), sweepParams.get_env_add(), forward, dmrginp.direct(),
                                      sweepParams.get_onedot(), nexact, useSlater, environment.get_integralIndex(), true, false, true);

  }
  StackSpinBlock big;
  newSystem.set_loopblock(true);
  system.set_loopblock(false);
  newEnvironment.set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 


  const int nroots = dmrginp.nroots();
  std::vector<StackWavefunction> solution;
  if(state==stateB){
    solution.resize(1);
    solution[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), true);
    solution[0].Clear();
    DiagonalMatrix e;
    GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0); 
#ifndef SERIAL
    MPI_Bcast(solution[0].get_data(), solution[0].memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

  }
  else{
    solution.resize(2);
    DiagonalMatrix e;
    solution[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(), true);
    solution[1].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_ketStateInfo(), big.get_rightBlock()->get_ketStateInfo(), true);
    solution[0].Clear();
    solution[1].Clear();
    GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.get_guesstype(), true, state, true, 0.0,false); 
    GuessWave::guess_wavefunctions(solution[1], e, big, sweepParams.get_guesstype(), true, stateB, true, 0.0,true); 
#ifndef SERIAL
    MPI_Bcast(solution[0].get_data(), solution[0].memoryUsed(), MPI_DOUBLE, 0, Calc);
    MPI_Bcast(solution[1].get_data(), solution[1].memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif
  }


  std::vector<Matrix> rotateMatrix;
  std::vector<Matrix> rotateMatrixB;

  if(state!=stateB){

    StackDensityMatrix tracedMatrix(newSystem.get_braStateInfo());
    tracedMatrix.allocate(newSystem.get_braStateInfo());
    tracedMatrix.makedensitymatrix(solution[0], big, 1.0);
    rotateMatrix.clear();

    StackDensityMatrix tracedMatrixB(newSystem.get_ketStateInfo());
    tracedMatrixB.allocate(newSystem.get_ketStateInfo());
    tracedMatrixB.makedensitymatrix(solution[1], big, 1.0);
    rotateMatrixB.clear();
    if (!mpigetrank()){
      double error = makeRotateMatrix(tracedMatrixB, rotateMatrixB, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
      if (dmrginp.bra_M() == 0)
        error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
      else
        error = makeRotateMatrix(tracedMatrix, rotateMatrix, dmrginp.bra_M(), sweepParams.get_keep_qstates());
    }
    tracedMatrixB.deallocate();
    tracedMatrix.deallocate();
  }
  else{
    StackDensityMatrix tracedMatrix(newSystem.get_stateInfo());
    tracedMatrix.allocate(newSystem.get_stateInfo());
    tracedMatrix.makedensitymatrix(solution[0], big, 1.0);
    rotateMatrix.clear();
    if (!mpigetrank()){
      double error = makeRotateMatrix(tracedMatrix, rotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    }
    tracedMatrix.deallocate();

  }


  


  int sweepPos = sweepParams.get_block_iter();
  int endPos = sweepParams.get_n_iters()-1;
  //cout <<"before "<< Stackmem[0].memused<<endl;
  size_t mem = Stackmem[0].memused;
  double *ptr = Stackmem[0].data+mem;
  npdm_driver.compute_npdm_elements(solution, big, sweepPos, endPos);
  //clear all the memory used so far
  Stackmem[0].deallocate(ptr, Stackmem[0].memused-mem);

  //cout <<"before "<< Stackmem[0].memused<<endl;

  SaveRotationMatrix (newSystem.get_sites(), rotateMatrix, state);
  solution[0].SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), state);
  if(state!=stateB){
    SaveRotationMatrix (newSystem.get_sites(), rotateMatrixB, stateB);
    solution[1].SaveWavefunctionInfo (big.get_ketStateInfo(), big.get_leftBlock()->get_sites(), stateB);
    solution[1].deallocate();
  }
  solution[0].deallocate();

  //FIXME
  //Maybe, for StateSpecific calculations, we can load rotation matrices, wavefuntions from the disk. 
  //There is no longer need to transform wavefuntions and to make rotation matrices from the density matrices.
  
  //FIXME
  //If in the state-average pdm, different states do not share the same rotation matrices as they do in energy calculations. Making rotation matrices from 
  //density matrices of different states is neccessary. 
  
  //if(newSystem.get_sites().size()>1)
  //if (!mpigetrank()){
  //LoadRotationMatrix (newSystem.get_sites(), rotateMatrix, state);
  //LoadRotationMatrix (newSystem.get_sites(), rotateMatrixB, stateB);
  //}
  #ifndef SERIAL
  mpi::broadcast(calc,rotateMatrix,0);
  if(state!=stateB)
    mpi::broadcast(calc,rotateMatrixB,0);
  #endif

  // Do we need to do this step for NPDM on the last sweep site? (It's not negligible cost...?)
  //It crashes at the last sweep site. 
  //Since it is useless, Just omit it at the last sweep site.
 // if( sweepParams.get_block_iter()  != sweepParams.get_n_iters() - 1)
  {
    if(state!=stateB)
      newSystem.transform_operators(rotateMatrix,rotateMatrixB);
    else {
      newSystem.transform_operators(rotateMatrix, rotateMatrix);
    }
  }

  {
    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysmem = newSystem.memoryUsed();
    newSystem.moveToNewMemory(system.getdata());
    Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
    //system.clear();
  }

  pout << newSystem <<endl;
  newSystem.printOperatorSummary();
  p2out << str(boost::format("%-40s - %-10.4f\n") % "Total memory" % (Stackmem[0].size*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Memory used" % (Stackmem[0].memused*8/1.e9));

  //newSystem.transform_operators(rotateMatrix,rotateMatrixB);
  double cputime = timer.elapsedcputime();
  double walltime = timer.elapsedwalltime();
  p3out << "NPDM block and decimate and compute elements " << walltime << " " << cputime << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

double npdm_do_one_sweep(Npdm_driver_base &npdm_driver, SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
                         const bool &restart, const int &restartSize, const int state, const int stateB)
{
  Timer sweeptimer;
  pout.precision(12);
  StackSpinBlock system;
  const int nroots = dmrginp.nroots();
  std::vector<double> finalEnergy(nroots,0.);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << ((forward) ? "\t\t\t Starting renormalisation sweep in forwards direction" : "\t\t\t Starting renormalisation sweep in backwards direction") << endl;
  pout << "\t\t\t ============================================================================ " << endl;
  
  int integralIndex = 0;
  if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep()) {
    InitBlocks::InitStartingBlock( system, forward, state, stateB, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  }
  else 
    InitBlocks::InitStartingBlock( system, forward, sweepParams.current_root(), sweepParams.current_root(), sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);

  pout << "\t\t\t Starting block is :: " << endl << system << endl;

  if (!restart) sweepParams.set_block_iter() = 0;
  if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep()){
    if (!restart) StackSpinBlock::store (forward, system.get_sites(), system, state, stateB ); // if restart, just restoring an existing block --
  }
  else{
    if (!restart) StackSpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root() ); // if restart, just restoring an existing block --
  }

  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;

  // Loop over all block sites
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) {
    Timer timer;
    pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
    pout << "\t\t\t ----------------------------" << endl;
    if (forward) { p1out << "\t\t\t Current direction is :: Forwards " << endl; }
    else { p1out << "\t\t\t Current direction is :: Backwards " << endl; }

    //if (SHOW_MORE) pout << "system block" << endl << system << endl;

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
    
    //pout << "guess wave funtion type: " << sweepParams.get_guesstype()<<endl;
    p1out << "\t\t\t Blocking and Decimating " << endl;
 
    StackSpinBlock newSystem;

    // Build npdm elements
    npdm_block_and_decimate(npdm_driver, sweepParams, system, newSystem, warmUp, dot_with_sys, state, stateB);

//    for(int j=0;j<nroots;++j)
//      pout << "\t\t\t Total block energy for State [ " << j << 
// " ] with " << sweepParams.get_keep_states()<<" :: " << sweepParams.get_lowest_energy()[j]+dmrginp.get_coreenergy() <<endl;              
//
//    finalEnergy_spins = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy_spins() : finalEnergy_spins);
//    finalEnergy = ((sweepParams.get_lowest_energy()[0] < finalEnergy[0]) ? sweepParams.get_lowest_energy() : finalEnergy);
//    finalError = max(sweepParams.get_lowest_error(),finalError);

    system = newSystem;

    pout << system<<endl;
    
    if(dmrginp.setStateSpecific() || dmrginp.transition_diff_irrep())
      StackSpinBlock::store (forward, system.get_sites(), system, state, stateB);
    else
      StackSpinBlock::store (forward, system.get_sites(), system, sweepParams.current_root(), sweepParams.current_root() );

    p1out << "\t\t\t saving state " << system.get_sites().size() << endl;
    ++sweepParams.set_block_iter();
    //sweepParams.savestate(forward, system.get_sites().size());

    double cputime = timer.elapsedcputime();
    p3out << "NPDM do one site time " << timer.elapsedwalltime() << " " << cputime << endl;
  }
  system.deallocate();
  system.clear();
  //for(int j=0;j<nroots;++j)
//  {int j = state;
//    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
//	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
//  }
//  {int j = stateB;
//    pout << "\t\t\t Finished Sweep with " << sweepParams.get_keep_states() << " states and sweep energy for State [ " << j 
//	 << " ] with Spin [ " << dmrginp.molecule_quantum().get_s()  << " ] :: " << finalEnergy[j]+dmrginp.get_coreenergy() << endl;
//  }
//  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
//  pout << "\t\t\t ============================================================================ " << endl;

  // Dump NPDM to disk if necessary
  npdm_driver.save_data( state, stateB );

  // Update the static number of iterations
  ++sweepParams.set_sweep_iter();

  double cputime = sweeptimer.elapsedcputime();
  pout << "\t\t\t Elapsed Sweep CPU  Time (seconds): " << std::setprecision(3) << cputime << endl;
  pout << "\t\t\t Elapsed Sweep Wall Time (seconds): " << std::setprecision(3) << sweeptimer.elapsedwalltime() << endl;

  return finalEnergy[0];

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void npdm(NpdmOrder npdm_order, bool transitionpdm)
{
  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;

  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }

  if (dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
    abort();
  }

  if(transitionpdm )
    dmrginp.setimplicitTranspose() = false;

  dmrginp.do_pdm() = true;

  // Screening can break things for NPDM (e.g. smaller operators won't be available from which to build larger ones etc...?)
  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  dmrginp.do_npdm_ops() = true;
  sweep_copy.restorestate(direction_copy, restartsize_copy);

  // Initialize npdm_driver
  boost::shared_ptr<Npdm_driver_base> npdm_driver;
  //if ( (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) && dmrginp.spinAdapted() ) {
  if ( (dmrginp.hamiltonian() == QUANTUM_CHEMISTRY) || (dmrginp.hamiltonian() == BCS)) {
    //By default, new_npdm_code is false.
    //For npdm_order 1 or 2. new_npdm_code is determined by default or manual setting.
    //For the other situation, only old or new code is suitable.
    //if(npdm_order == NPDM_PAIRMATRIX || npdm_order == NPDM_THREEPDM || npdm_order == NPDM_FOURPDM || npdm_order == NPDM_NEVPT2 ||  transitionpdm == true  || dmrginp.spinAdapted() == false || dmrginp.setStateSpecific())
    dmrginp.set_new_npdm_code();

    if(dmrginp.new_npdm_code()){
    if      (npdm_order == NPDM_ONEPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Onepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_TWOPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Twopdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_THREEPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Threepdm_driver( dmrginp.last_site() ) );
    else if (npdm_order == NPDM_FOURPDM) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Fourpdm_driver( dmrginp.last_site() ) );
    //else if (npdm_order == NPDM_NEVPT2) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Nevpt2_npdm_driver( dmrginp.last_site() ) );
    //else if (npdm_order == NPDM_PAIRMATRIX) npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Pairpdm_driver( dmrginp.last_site() ) );
    else abort();
    }
  }

  if (dmrginp.specificpdm().size()!=0)
  {
    Timer timer;
    dmrginp.set_fullrestart() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
	  dmrginp.npdm_generate() = true;
    if (dmrginp.specificpdm().size()==1)
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, dmrginp.specificpdm()[0], dmrginp.specificpdm()[0]); //this will generate the cd operators                               
    else if (dmrginp.specificpdm().size()==2)
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, dmrginp.specificpdm()[0], dmrginp.specificpdm()[1]); //this will generate the cd operators                               
    else 
      abort();
		dmrginp.npdm_generate() = false;
    double ecpu = timer.elapsedcputime();
    double ewall=timer.elapsedwalltime();
    p3out << "\t\t\t NPDM SweepGenblock time " << ewall << " " << ecpu << endl;
    dmrginp.set_fullrestart() = false;

    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    Timer timerX;
    npdm_driver->clear();
    if (dmrginp.specificpdm().size()==1)
      npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0,dmrginp.specificpdm()[0],dmrginp.specificpdm()[0]);
    else if (dmrginp.specificpdm().size()==2)
      npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0,dmrginp.specificpdm()[0],dmrginp.specificpdm()[1]);
    else
      abort();
    ecpu = timerX.elapsedcputime();ewall=timerX.elapsedwalltime();
    p3out << "\t\t\t NPDM sweep time " << ewall << " " << ecpu << endl;
    return;
  }



  if(dmrginp.transition_diff_irrep()){
    // It is used when bra and ket has different spatial irrep
    // For now, only the transtion pdm between two wavefuntions( 1 as bra and 0 as ket) are calculation
    // If the spatial irrep information is stored in wavefuntions, transition pdm among  several wavefunctions( i for bra and j for ket, there are n(n-1)/2 kinds of situations.) are possible.
    for (int state=0; state<dmrginp.nroots(); state++) {
      for(int stateB=0; stateB<state; stateB++){
        Timer timer;
        dmrginp.set_fullrestart() = true;
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
				dmrginp.npdm_generate() = true;
        SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, stateB); //this will generate the cd operators                               
				dmrginp.npdm_generate() = false;
        double cputime = timer.elapsedcputime();
        p3out << "\t\t\t NPDM SweepGenblock time " << timer.elapsedwalltime() << " " << cputime << endl;
        dmrginp.set_fullrestart() = false;

        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        Timer timerX;
        npdm_driver->clear();
        npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,stateB);
        cputime = timerX.elapsedcputime(); 
        p3out << "\t\t\t NPDM sweep time " << timerX.elapsedwalltime() << " " << cputime << endl;
       }
    }
  }

  else if( !dmrginp.setStateSpecific()){
    Timer timer;
    dmrginp.set_fullrestart() = true;
    dmrginp.npdm_generate() = true;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;

    dmrginp.setimplicitTranspose() = false;
    if (dmrginp.new_npdm_code())
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, -1, -1); //this will generate the cd operators                               

    if(!transitionpdm )
      dmrginp.setimplicitTranspose() = true;

    dmrginp.npdm_generate() = false;
    
    double cputime = timer.elapsedcputime();
    p3out << "\t\t\t NPDM SweepGenblock time " << timer.elapsedwalltime() << " " << cputime << endl;
    dmrginp.set_fullrestart() = false;
    
    
    if(transitionpdm){
      //  <\Phi_k|a^+_ia_j|\Phi_l> = <\Phi_l|a^+_ja_i|\Phi_k>*
      //  Therefore, only calculate the situations with k >= l.
      //  for(int stateB=0; stateB<= dmrginp.nroots(); stateB++){
      for (int state=0; state<dmrginp.nroots(); state++) {
	for(int stateB=0; stateB<=state; stateB++){
          sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
          Timer timerX;

          npdm_driver->clear();
          npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,stateB);
          double cputime = timerX.elapsedcputime();
          p3out << "\t\t\t NPDM sweep time " << timerX.elapsedwalltime() << " " << cputime << endl;

          if (dmrginp.hamiltonian() == BCS && npdm_order == NPDM_ONEPDM) {
            Timer timerX1;            
            npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Pairpdm_driver( dmrginp.last_site() ) );
            npdm_driver->clear();
            npdm_do_one_sweep(*npdm_driver, sweepParams, false, !direction, false, 0, state,stateB);
            cputime = timerX1.elapsedcputime();
            p3out << "\t\t\t NPDM sweep time " << timerX1.elapsedwalltime() << " " << cputime << endl; 
          }
	}
      }
      
    }
    else {
      for (int state=0; state<dmrginp.nroots(); state++) {
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        if ( dmrginp.new_npdm_code() ) {
          Timer timerX;
          npdm_driver->clear();
          npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,state);
          double cputime = timerX.elapsedcputime();
          p3out << "\t\t\t NPDM sweep time " << timerX.elapsedwalltime() << " " << cputime << endl;
          if (dmrginp.hamiltonian() == BCS && npdm_order == NPDM_ONEPDM) {
            Timer timerX1;            
            npdm_driver = boost::shared_ptr<Npdm_driver_base>( new Pairpdm_driver( dmrginp.last_site() ) );
            npdm_driver->clear();
            npdm_do_one_sweep(*npdm_driver, sweepParams, false, !direction, false, 0, state,state);
            cputime = timerX1.elapsedcputime();
            p3out << "\t\t\t NPDM sweep time " << timerX1.elapsedwalltime() << " " << cputime << endl; 
          }
        } 
        else{
	  pout << "Old code is not available anymore "<<endl;
	  print_trace(2);
	  //SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, state); //this will generate the cd operators                               
          //if (npdm_order == NPDM_ONEPDM) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);     
          //else if (npdm_order == NPDM_TWOPDM) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state, state);
          //else abort();
        }
      }
    }
    
  }

  else {
    // state-specific
    if(transitionpdm){
      for (int state=0; state<dmrginp.nroots(); state++) {
        for (int stateB=0; stateB<=state; stateB++) {
  
        dmrginp.set_fullrestart() = true;
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
  
        if (mpigetrank() == 0) {
          Sweep::InitializeStateInfo(sweepParams, direction, state);
          Sweep::InitializeStateInfo(sweepParams, !direction, state);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
          Sweep::CanonicalizeWavefunction(sweepParams, !direction, state);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
        }
  
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        if (mpigetrank() == 0) {
          Sweep::InitializeStateInfo(sweepParams, direction, stateB);
          Sweep::InitializeStateInfo(sweepParams, !direction, stateB);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, stateB);
          Sweep::CanonicalizeWavefunction(sweepParams, !direction, stateB);
          Sweep::CanonicalizeWavefunction(sweepParams, direction, stateB);
        }
        // Prepare NPDM operators
        Timer timer;
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
				dmrginp.npdm_generate() = true;
        SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, stateB); //this will generate the cd operators
				dmrginp.npdm_generate() = false;
        dmrginp.set_fullrestart() = false;
        double cputime = timer.elapsedcputime();
        p3out << "\t\t\t NPDM SweepGenblock time " << timer.elapsedwalltime() << " " << cputime << endl;
  
        Timer timerX;
        npdm_driver->clear();
        sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
        npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,stateB);
        cputime = timerX.elapsedcputime();
        p3out << "\t\t\t NPDM sweep time " << timerX.elapsedwalltime() << " " << cputime << endl;
      }
    }
  }
  else{
    for (int state=0; state<dmrginp.nroots(); state++) {
      dmrginp.set_fullrestart() = true;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;

      if (mpigetrank() == 0) {
        Sweep::InitializeStateInfo(sweepParams, direction, state);
        Sweep::InitializeStateInfo(sweepParams, !direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, !direction, state);
        Sweep::CanonicalizeWavefunction(sweepParams, direction, state);
      }
      // Prepare NPDM operators
			dmrginp.npdm_generate() = true;
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, state); //this will generate the cd operators
			dmrginp.npdm_generate() = false;
      dmrginp.set_fullrestart() = false;
      // Do NPDM sweep
      npdm_driver->clear();
      npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, state,state);
   }
  }
  
}
  sweep_copy.savestate(direction_copy, restartsize_copy);
}

}
}

