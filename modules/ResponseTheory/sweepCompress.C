/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "stackguess_wavefunction.h"
#include "sweepCompress.h"
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
#include "davidson.h"
#include "sweep.h"
#include "sweepResponse.h"
#include "Stackwavefunction.h"
#include "Stackspinblock.h"
#include "sweep_params.h"

using namespace boost;
using namespace std;


void SpinAdapted::SweepCompress::BlockDecimateAndCompress (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, int baseState)
{
  int sweepiter = sweepParams.get_sweep_iter();

  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;

  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot;
  StackSpinBlock environment, environmentDot, newEnvironment;
  StackSpinBlock envDot, big;
  int systemDotStart, systemDotEnd;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() -1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];

  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), targetState, baseState);


  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot()) { pout << "\t\t\t System  Block"<<system; }   
  else { pout << "\t\t\t System  Block"<<newSystem; }
  pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<StackWavefunction> solution; solution.resize(1);
  solution[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_ketStateInfo(), big.get_rightBlock()->get_ketStateInfo(), sweepParams.get_onedot());
  solution[0].set_onedot(sweepParams.get_onedot());
  solution[0].Clear();
  DiagonalMatrix e;


  //read the 0th wavefunction which we keep on the ket side because by default the ket stateinfo is used to initialize wavefunction
  //also when you use spinblock operators to multiply a state, it does so from the ket side i.e.  H|ket>

  //**********************
  GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 1, 0.0, baseState); 
#ifndef SERIAL
  mpi::communicator world;
  MPI_Bcast(solution[0].get_data(), solution[0].memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif
  
  //*********************
  multiply_h_2Index davidson_f(big, sweepParams.get_onedot());
  vector<StackWavefunction> outputState; outputState.resize(1);
  outputState[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(), sweepParams.get_onedot());
  outputState[0].set_onedot(sweepParams.get_onedot());
  outputState[0].Clear();

  //*************
  davidson_f(solution[0], outputState[0]);
  double overlap = 1.0;

  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
  sweepParams.set_lowest_energy() = std::vector<double>(1,overlap);

  StackSpinBlock newbig;

  if (sweepParams.get_onedot() && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, baseState, targetState, systemDot.size(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newSystem, environment, newbig); 

    StackWavefunction tempwave; tempwave.initialise(solution[0]); tempwave.Clear();
    GuessWave::onedot_shufflesysdot(big.get_ketStateInfo(), newbig.get_ketStateInfo(), solution[0], tempwave);  
    //copy(tempwave.get_operatorMatrix(), solution[0].get_operatorMatrix());
    double* backup = solution[0].get_data();
    solution[0] = tempwave;
    solution[0].set_data(backup);
    DCOPY(solution[0].memoryUsed(), tempwave.get_data(), 1, solution[0].get_data(), 1);
    solution[0].allocateOperatorMatrix();
    tempwave.deallocate();

    tempwave.initialise(outputState[0]); tempwave.Clear();
    GuessWave::onedot_shufflesysdot(big.get_braStateInfo(), newbig.get_braStateInfo(), outputState[0], tempwave);  
    backup = outputState[0].get_data();
    outputState[0] = tempwave;
    outputState[0].set_data(backup);
    DCOPY(outputState[0].memoryUsed(), tempwave.get_data(), 1, outputState[0].get_data(), 1);
    outputState[0].allocateOperatorMatrix();
    tempwave.deallocate();

    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;

  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo()), kettracedMatrix(newSystem.get_ketStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo());

  kettracedMatrix.allocate(newSystem.get_ketStateInfo());
  //bratracedMatrix.allocate(newSystem.get_braStateInfo()); kettracedMatrix.allocate(newSystem.get_ketStateInfo());

//**********************
  bratracedMatrix.makedensitymatrix(outputState, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
  if (sweepParams.get_noise() > NUMERICAL_ZERO) {
    pout << "adding noise  "<<trace(bratracedMatrix)<<"  "<<sweepiter<<"  "<<dmrginp.weights(sweepiter)[0]<<endl;

    double* backupData;
    if (mpigetrank() == 0) {
      backupData = Stackmem[omprank].allocate(bratracedMatrix.memoryUsed());
      DCOPY(bratracedMatrix.memoryUsed(), bratracedMatrix.get_data(), 1, backupData, 1);
    }
    bratracedMatrix.Clear();
    //************************
    bratracedMatrix.add_onedot_noise(solution[0], newbig);

    if (mpigetrank() == 0) {
      double norm = trace(bratracedMatrix);
      if (fabs(norm) > 1e-8) {
	DAXPY(bratracedMatrix.memoryUsed(), sweepParams.get_noise()*max(1.0, trace(bratracedMatrix))/norm, bratracedMatrix.get_data(), 1, backupData, 1);
      }
      DCOPY(bratracedMatrix.memoryUsed(), &backupData[0], 1, bratracedMatrix.get_data(), 1);
      
      if (trace(bratracedMatrix) <1e-14) 
	bratracedMatrix.SymmetricRandomise();
      
      pout << "after noise  "<<trace(bratracedMatrix)<<"  "<<sweepParams.get_noise()<<endl;
      Stackmem[omprank].deallocate(backupData, bratracedMatrix.memoryUsed());
    }
  }

  //****************************
  kettracedMatrix.makedensitymatrix(solution, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
  double braerror, keterror;
  if (!mpigetrank()) {
    keterror = makeRotateMatrix(kettracedMatrix, ketrotateMatrix, newbig.get_rightBlock()->get_ketStateInfo().totalStates, 0);
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }
  kettracedMatrix.deallocate();
  bratracedMatrix.deallocate();

#ifndef SERIAL
  broadcast(calc, ketrotateMatrix, 0);
  broadcast(calc, brarotateMatrix, 0);
#endif

  //assert(keterror < NUMERICAL_ZERO);
  pout << "\t\t\t Total ket discarded weight "<<keterror<<endl<<endl;
  pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;

  sweepParams.set_lowest_error() = braerror;

  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);
  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, targetState);
  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
  outputState[0].SaveWavefunctionInfo (newbig.get_braStateInfo(), newbig.get_leftBlock()->get_sites(), targetState);

  outputState[0].deallocate();
  solution[0].deallocate();

  environment.clear();
  newEnvironment.removeAdditionalOps();
  newEnvironment.deallocate();
  newEnvironment.clear();
  environment.removeAdditionalOps();
  environment.clear();
  environment.deallocate();


  p1out <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix, false, false);

  //if (system.get_sites().size() != 1) {
  //if (system.get_sites().size() != 1 || (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted())) {
  {
    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysmem = newSystem.memoryUsed();
    newSystem.moveToNewMemory(system.getdata());
    Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
    system.clear();
  }


  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;

}

double SpinAdapted::SweepCompress::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int targetState, int baseState)
{
  int integralIndex = 0;
  StackSpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  std::vector<double> finalEnergy(nroots,-1.0e10);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;
  if (restart) {
    finalEnergy = sweepParams.get_lowest_energy();
    finalEnergy_spins = sweepParams.get_lowest_energy();
    finalError = sweepParams.get_lowest_error();
  }

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward)
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl; }
  else
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl; }
  pout << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system,forward, baseState, targetState, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  StackSpinBlock::store (forward, system.get_sites(), system, targetState, baseState); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (restart)
    Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);
   
  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }

      if (sweepParams.get_block_iter() != 0) 
	sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp )
	Startup(sweepParams, system, newSystem, dot_with_sys, targetState, baseState);
      else {
	BlockDecimateAndCompress (sweepParams, system, newSystem, false, dot_with_sys, targetState, baseState);
      }
      
      //Need to substitute by?

      if (!warmUp ){

	//this criteria should work for state average or state specific because the lowest sweep energy is always the lowest of the average
	finalError = max(sweepParams.get_lowest_error(),finalError);
	finalEnergy[0] = max(sweepParams.get_lowest_energy()[0], finalEnergy[0]);
	pout << "final energy "<<finalEnergy[0]<<"  "<<sweepParams.get_lowest_energy()[0]<<endl;
      }
      
      system = newSystem;
      p2out << system<<endl;
      p2out << system.get_braStateInfo()<<endl;
      system.printOperatorSummary();
      
      //system size is going to be less than environment size
      Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);
      
      StackSpinBlock::store (forward, system.get_sites(), system, targetState, baseState);	 	
      syssites = system.get_sites();
      p1out << "\t\t\t saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      
#ifndef SERIAL
      mpi::communicator world;
      calc.barrier();
#endif
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }


  //when we are doing twodot, we still need to do the last sweep to make sure that the
  //correctionVector and base wavefunction are propogated correctly across sweeps
  //especially when we switch from twodot to onedot algorithm
  if (!sweepParams.get_onedot() && !warmUp) {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
    sweepParams.set_onedot() = true;
    sweepParams.set_env_add() = 0;
    bool dot_with_sys = true;
    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, baseState);
    sweepParams.set_onedot() = false;
    sweepParams.set_env_add() = 1;
  }

  //system.deallocate();

  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  pout << "\t\t\t Largest overlap for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalError;
}


void SpinAdapted::SweepCompress::Startup (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool& dot_with_sys, int targetState, int baseState)
{
  bool useSlater = false;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);

  StackSpinBlock systemDot;
  StackSpinBlock environment, environmentDot, newEnvironment;
  StackSpinBlock envDot, big;
  int systemDotStart, systemDotEnd;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() -1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];

  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), baseState, baseState);
  
  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot()) {
    pout << "\t\t\t System  Block"<<system;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  }
  else {
    pout << "\t\t\t System  Block"<<newSystem;
    pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  }
  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<StackWavefunction> solution; solution.resize(1);
  solution[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); solution[0].Clear();

  DiagonalMatrix e;
  e.ReSize(big.get_stateInfo().totalStates); e= 0;

  //read the 0th wavefunction which we keep on the ket side because by default the ket stateinfo is used to initialize wavefunction
  //also when you use spinblock operators to multiply a state, it does so from the ket side i.e.  H|ket>
  //**********************
  GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 1, 0.0, baseState); 

  StackSpinBlock newbig;

  if (sweepParams.get_onedot() && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, targetState, baseState, systemDot.size(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newSystem, environment, newbig); 

    StackWavefunction tempwave; tempwave.initialise(solution[0]);
    tempwave.Clear();
    GuessWave::onedot_shufflesysdot(big.get_ketStateInfo(), newbig.get_ketStateInfo(), solution[0], tempwave);  
    DCOPY(solution[0].memoryUsed(), tempwave.get_data(), 1, solution[0].get_data(), 1);

    //cout << "After shuffle "<<solution[0]<<endl;

    tempwave.deallocate();

    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;


  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));

  std::vector<Matrix> ketrotateMatrix, brarotateMatrix;

  StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo()), kettracedMatrix(newSystem.get_ketStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo());
  kettracedMatrix.allocate(newSystem.get_ketStateInfo());


  //************************
  bratracedMatrix.makedensitymatrix(solution, newbig, dmrginp.weights(0), 0.0, 
				    0.0, true);


  //******************************
  kettracedMatrix.makedensitymatrix(solution, newbig, dmrginp.weights(0), 0.0, 
				    0.0, true);

  double keterror, braerror;
  if (!mpigetrank()) {
    keterror = makeRotateMatrix(kettracedMatrix, ketrotateMatrix, newbig.get_rightBlock()->get_ketStateInfo().totalStates, 0);
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, newbig.get_rightBlock()->get_braStateInfo().totalStates, 0);
  }
#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, ketrotateMatrix, 0);
  broadcast(calc, brarotateMatrix, 0);
#endif
  kettracedMatrix.deallocate();
  bratracedMatrix.deallocate();

  //assert(keterror < NUMERICAL_ZERO);
  p1out <<"\t\t\t Performing Renormalization "<<endl;
  pout << "\t\t\t Total ket discarded weight "<<keterror<<endl<<endl;
  pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;
  sweepParams.set_lowest_error() = keterror;

  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);
  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, targetState);
  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
  solution[0].SaveWavefunctionInfo (newbig.get_braStateInfo(), newbig.get_leftBlock()->get_sites(), targetState);

  solution[0].deallocate();

  environment.clear();
  newEnvironment.deallocate();
  newEnvironment.clear();
  environment.removeAdditionalOps();
  environment.clear();
  environment.deallocate();

  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix, false, false);

  long memoryToFree = newSystem.getdata() - system.getdata();
  long newsysmem = newSystem.memoryUsed();
  newSystem.moveToNewMemory(system.getdata());
  Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
  system.clear();


  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;

}

void SpinAdapted::SweepCompress::WavefunctionCanonicalize (SweepParams &sweepParams, StackSpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int correctionVector, int baseState)
{
  if (dmrginp.outputlevel() > 0)
    mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  
  StackSpinBlock newSystem;
  
  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = 0;
  if (forward)
    {
      systemDotStart = dmrginp.spinAdapted() ? 
	*system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
      systemDotEnd = systemDotStart + systemDotSize;
      environmentDotStart = systemDotEnd + 1;
      environmentDotEnd = environmentDotStart + environmentDotSize;
    }
  else
    {
      systemDotStart = dmrginp.spinAdapted() ? 
	system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
      systemDotEnd = systemDotStart - systemDotSize;
      environmentDotStart = systemDotEnd - 1;
      environmentDotEnd = environmentDotStart - environmentDotSize;
    }
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];
  //StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  vector<int> sitesenvdot(environmentDotSize+1, 0);
  int index = 0;
  for (int i=min(environmentDotStart, environmentDotEnd); i<max(environmentDotStart, environmentDotEnd)+1; i++) {
    sitesenvdot[index] = (i);
    index++;
  }

  //environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];
  StackSpinBlock::restore(!forward, sitesenvdot, environmentDot, correctionVector, baseState); 

  StackSpinBlock environment, newEnvironment;
  
  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  system.addAdditionalOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, correctionVector, baseState, sweepParams.get_sys_add(), dmrginp.direct(), 
				 system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);

  newSystem.set_loopblock(false);  environmentDot.set_loopblock(false); 
  InitBlocks::InitBigBlock(newSystem, environmentDot, big);

  
  
  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;
  
  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot()) { pout << "\t\t\t System  Block"<<system; }   
  else { pout << "\t\t\t System  Block"<<newSystem; }
  pout << "\t\t\t Environment Block"<<environmentDot<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;
  
  
  
  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  //dmrginp.setOutputlevel() = -1;
  
  
  
  std::vector< StackWavefunction> solution(1);
  solution[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); solution[0].Clear();

  //**************************
  if (!mpigetrank())
    GuessWave::transform_previous_twodot_to_onedot_wavefunction(solution[0], big, baseState);
  solution[0].set_onedot(true);

#ifndef SERIAL
  mpi::communicator world;
  MPI_Bcast(solution[0].get_data(), solution[0].memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

  //****************
  multiply_h_2Index davidson_f(big, sweepParams.get_onedot());
  vector<StackWavefunction> outputState; outputState.resize(1);
  outputState[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(), true);
  outputState[0].set_onedot(true);
  outputState[0].Clear();

  //***********
  davidson_f(solution[0], outputState[0]);

  pout << solution[0]<<endl;
  pout << outputState[0]<<endl;
  
  std::vector<Matrix> ketrotateMatrix, brarotateMatrix;
  StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo()), kettracedMatrix(newSystem.get_ketStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo());
  kettracedMatrix.allocate(newSystem.get_ketStateInfo());

  bratracedMatrix.makedensitymatrix(outputState, big, dmrginp.weights(0), 0.0, 0.0, true);
  kettracedMatrix.makedensitymatrix(solution, big, dmrginp.weights(0), 0.0, 0.0, true);
  double braerror, keterror;
  int largeNumber = 1000000;
  if (!mpigetrank()) {
    keterror = makeRotateMatrix(kettracedMatrix, ketrotateMatrix, largeNumber, 0);
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }
  kettracedMatrix.deallocate();
  bratracedMatrix.deallocate();

#ifndef SERIAL
  broadcast(calc, ketrotateMatrix, 0);
  broadcast(calc, brarotateMatrix, 0);
#endif


  pout << "\t\t\t Total ket discarded weight "<<keterror<<endl<<endl;
  pout << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;

  sweepParams.set_lowest_error() = braerror;

  SaveRotationMatrix (big.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);
  SaveRotationMatrix (big.get_leftBlock()->get_sites(), brarotateMatrix, correctionVector);
  solution[0].SaveWavefunctionInfo (big.get_ketStateInfo(), big.get_leftBlock()->get_sites(), baseState);
  outputState[0].SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), correctionVector);

  outputState[0].Clear();
  solution[0].Clear();

  //environmentDot.deallocate();
  p1out <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix, false, false);


  
  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  
}
