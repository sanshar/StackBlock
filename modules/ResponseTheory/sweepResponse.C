/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include "Stackdensity.h"
#include "Stackwavefunction.h"
#include "stackguess_wavefunction.h"
#include "sweepResponse.h"
#include "IntegralMatrix.h"
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
//#include "density.h"
#include "Stackdensity.h"
#include "pario.h"
#include "davidson.h"
#include "initblocks.h"
#include "operatorfunctions.h"
#include "Stackspinblock.h"
#include "sweep_params.h"
#include "sweep.h"

using namespace boost;
using namespace std;


void SpinAdapted::SweepResponse::BlockAndDecimate (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, vector<int>& projectors, vector<int>& baseStates)
{
  if (dmrginp.outputlevel() > 0)
    mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();

  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() - 1;
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
  pout << "**** STACK MEMORY REMAINING before block***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
  StackSpinBlock environment, newEnvironment;

  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)
  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, 
					newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(),  
					targetState, targetState);

  pout << "**** STACK MEMORY REMAINING after block***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;

  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 0)
   mcheck(""); 
  pout << "\t\t\t System  Block"<<*big.get_leftBlock();
  pout << "\t\t\t Environment Block"<<*big.get_rightBlock()<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<StackWavefunction> lowerStates;

  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
  lowerStates.resize(projectors.size()+1); //only make the 
  lowerStates[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot());
  lowerStates[0].Clear();
  if (mpigetrank() == 0) {
    for (int i=1; i<lowerStates.size(); i++) {
      lowerStates[i].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot());
      lowerStates[i].Clear();
    }
  }

  DiagonalMatrix e;
  guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;

  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }

  //<target|O|firstOrderState>
  StackDensityMatrix branoiseMatrix;
  branoiseMatrix.allocate(newSystem.get_braStateInfo()); branoiseMatrix.Clear();

  if (dmrginp.outputlevel() > 0)
    mcheck("before making correction vector");
  for (int l=0; l<baseStates.size(); l++)
  {  
    //now one needs to make |phi_0> = O|psi_0> so that the |phi_0> has the same dimensions as our target state
    int firstOrderState = baseStates[l];

    int perturbationIntegral = l+1;
    StackSpinBlock perturbationBig, perturbationsystemdot, perturbationenvironmentdot;
    StackSpinBlock perturbationsystem, perturbationenvironment, perturbationnewsystem, perturbationnewenvironment;
    perturbationsystemdot = StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, false);
    perturbationenvironmentdot = StackSpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, false);
    perturbationsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, perturbationsystem, targetState, firstOrderState);
    perturbationsystem.set_twoInt(perturbationIntegral);

    Sweep::makeSystemEnvironmentBigBlocks(perturbationsystem, perturbationsystemdot, 
					  perturbationnewsystem, perturbationenvironment, 
					  perturbationenvironmentdot, perturbationnewenvironment, 
					  perturbationBig, sweepParams, dot_with_sys, useSlater,
					  perturbationIntegral, targetState, firstOrderState);
    
    StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), perturbationBig.get_leftBlock()->get_stateInfo(), perturbationBig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();

    //************************************
    GuessWave::guess_wavefunctions(iwave, e, perturbationBig, guesstype, 
				   sweepParams.get_onedot(), firstOrderState, dot_with_sys, 0.0);

#ifndef SERIAL
    mpi::communicator world;
    MPI_Bcast(iwave.get_data(), iwave.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

    //dont add noise in the onedot algorithm
    if (!sweepParams.get_onedot() && sweepParams.get_noise() > NUMERICAL_ZERO && l == 0) { //only add noise using one basestate
      int sweepiter = sweepParams.get_sweep_iter();

      //**********************************
      //if (!(dmrginp.calc_type() == RESPONSEAAAV && systemDotStart >= dmrginp.num_occupied_orbitals()))

      pout << perturbationBig<<endl;
      branoiseMatrix.add_onedot_noise(iwave, perturbationBig);
      DSCAL( branoiseMatrix.memoryUsed(), 1.0/DotProduct(iwave, iwave)/trace(branoiseMatrix), branoiseMatrix.get_data() , 1);
    }


    StackWavefunction temp; 
    temp.set_onedot( sweepParams.get_onedot());
    temp.initialise(dmrginp.effective_molecule_quantum_vec(), *perturbationBig.get_braStateInfo().leftStateInfo,*perturbationBig.get_braStateInfo().rightStateInfo, temp.get_onedot());
    temp.Clear();

    perturbationBig.multiplyH_2index(iwave, &temp, MAX_THRD);

    if (mpigetrank() == 0) {
      if(l==0)
	DCOPY(temp.memoryUsed(), temp.get_data(), 1, lowerStates[0].get_data(), 1);
      else
	ScaleAdd(1.0, temp, lowerStates[0]);
    }

    temp.deallocate();
    iwave.deallocate();

    perturbationnewenvironment.deallocate();
    perturbationenvironment.removeAdditionalOps();
    perturbationenvironment.deallocate();

    perturbationnewsystem.deallocate(); 
    perturbationsystem.removeAdditionalOps();
    perturbationsystem.deallocate(); 
    perturbationenvironmentdot.deallocate();
    perturbationsystemdot.deallocate();
  }

#ifndef SERIAL
    mpi::communicator world;
    MPI_Bcast(lowerStates[0].get_data(), lowerStates[0].memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif
  pout << "**** STACK MEMORY REMAINING before projector***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;

  //<target|O|projectors>
  for (int l=0; l<projectors.size(); l++)
  {  
    
    //now one needs to make |phi_0> = O|psi_0> so that the |phi_0> has the same dimensions as our target state
    int perturbationIntegral = 0;
    StackSpinBlock overlapBig, overlapsystemdot, overlapenvironmentdot;
    StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
    overlapsystemdot = StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, false);
    overlapenvironmentdot = StackSpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, false);
    overlapsystem.set_integralIndex() = perturbationIntegral;
    overlapsystem.set_sites() = systemsites; 

    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, overlapsystemdot, overlapenvironmentdot, 
						 overlapsystem, overlapnewsystem, overlapenvironment, 
						 overlapnewenvironment, overlapBig, sweepParams, dot_with_sys, useSlater,
						 perturbationIntegral, targetState, projectors[l]);
    
    StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();

    //****************************
    GuessWave::guess_wavefunctions(iwave, e, overlapBig, guesstype, 
				   sweepParams.get_onedot(), projectors[l], dot_with_sys, 0.0);
    
#ifndef SERIAL
    mpi::communicator world;
    MPI_Bcast(iwave.get_data(), iwave.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

    if (mpigetrank() == 0)
      overlapBig.multiplyOverlap(iwave, &lowerStates[l+1], MAX_THRD);


    int success = 0;
    if (mpigetrank() == 0) {
      for (int istate=1; istate<l+1; istate++)  {
	double overlap = pow(DotProduct(lowerStates[istate], lowerStates[istate]), 0.5);
	ScaleAdd(-DotProduct(lowerStates[istate], lowerStates[l+1])/overlap, lowerStates[istate], lowerStates[l+1]);
      }
    
      lowerStates[l+1].Normalise(&success);
    }
    iwave.deallocate();

    overlapnewenvironment.deallocate();
    overlapenvironment.deallocate();

    overlapnewsystem.deallocate(); 
    overlapsystem.deallocate(); 

    overlapenvironmentdot.deallocate();
    overlapsystemdot.deallocate();

  }


  dmrginp.setOutputlevel() = originalOutputlevel;

  StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo());
  if (mpigetrank() == 0) {
    bratracedMatrix.allocate(newSystem.get_braStateInfo());
    bratracedMatrix.Clear();
  }
  pout << "**** STACK MEMORY REMAINING before renormalization***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;

  newSystem.RenormaliseFrom (sweepParams.set_lowest_energy(), sweepParams.set_lowest_energy_spins(),
			     sweepParams.set_lowest_error(), rotatematrix, 
			     sweepParams.get_keep_states(), sweepParams.get_keep_qstates(), 
			     sweepParams.get_davidson_tol(), big, sweepParams.get_guesstype(), 
			     sweepParams.get_noise(), sweepParams.get_additional_noise(), //noise 
			     sweepParams.get_onedot(), system, systemDot, environment, 
			     dot_with_sys, useSlater, sweepParams.get_sweep_iter(), targetState, 
			     lowerStates, &bratracedMatrix);



  p1out <<"\t\t\t Performing Renormalization "<<endl;

  //rotatematrix.resize(0);
  pout << "**** STACK MEMORY REMAINING before renormalization***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;

  if(mpigetrank() == 0) {
    //ScaleAdd(sweepParams.get_noise()*(max(1.e-5, trace(branoiseMatrix))), branoiseMatrix, bratracedMatrix);
    //sweepParams.set_lowest_error() = makeRotateMatrix(bratracedMatrix, rotatematrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
    bratracedMatrix.deallocate();
  }
  branoiseMatrix.deallocate();

  if (mpigetrank() == 0) {
    for (int istate=lowerStates.size()-1; istate>0; istate--)
      lowerStates[istate].deallocate();
  }
  lowerStates[0].deallocate();

  //deallocate all the envirnonment blocks
  newEnvironment.deallocate();
  environment.removeAdditionalOps();
  environment.deallocate();


  StackWavefunction targetWave; StateInfo braStateInfo;
  targetWave.initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); targetWave.Clear();
  targetWave.LoadWavefunctionInfo (braStateInfo, newSystem.get_sites(), targetState);
  targetWave.deallocate();

#ifndef SERIAL
  broadcast(calc, rotatematrix, 0);
#endif

  //<target|O|firstOrderState>
  dmrginp.setOutputlevel() = -1; 
  for (int l=0; l<baseStates.size(); l++) 
  {
    int firstOrderState = baseStates[l];
    int originalOutputlevel = dmrginp.outputlevel();
    //dmrginp.setOutputlevel() = -1;

    int perturbationIntegral = l+1;
    StackSpinBlock perturbationBig, perturbationsystemdot, perturbationenvironmentdot;
    StackSpinBlock perturbationsystem, perturbationenvironment, perturbationnewsystem, perturbationnewenvironment;
    perturbationsystemdot = StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, false);
    perturbationenvironmentdot = StackSpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, false);
    perturbationsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, perturbationsystem, targetState, firstOrderState);
    perturbationsystem.set_twoInt(perturbationIntegral);

    Sweep::makeSystemEnvironmentBigBlocks(perturbationsystem, perturbationsystemdot, 
					  perturbationnewsystem, perturbationenvironment, 
					  perturbationenvironmentdot, perturbationnewenvironment, 
					  perturbationBig, sweepParams, dot_with_sys, useSlater,
					  perturbationIntegral, targetState, firstOrderState);

    StackSpinBlock perturbationnewbig;
    if (sweepParams.get_onedot() && !dot_with_sys)
    {
      InitBlocks::InitNewSystemBlock(perturbationsystem, perturbationsystemdot, 
				     perturbationnewsystem, targetState, firstOrderState, 
				     perturbationsystemdot.size(), dmrginp.direct(), 
				     perturbationIntegral, DISTRIBUTED_STORAGE, false, true);
      InitBlocks::InitBigBlock(perturbationnewsystem, perturbationenvironment, perturbationnewbig); 

      perturbationBig.get_rightBlock()->clear();
      perturbationBig.clear();
    }
    else
      perturbationnewbig = perturbationBig;

    std::vector<Matrix> ketrotatematrix;
    if (mpigetrank() == 0) {
      StackWavefunction iwave;
      iwave.initialise(dmrginp.effective_molecule_quantum_vec(), perturbationnewbig.get_leftBlock()->get_stateInfo(), perturbationnewbig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();
      //*************************
      GuessWave::guess_wavefunctions(iwave, e, perturbationnewbig, guesstype, 
				     sweepParams.get_onedot(), firstOrderState, true, 0.0);
      
      StackDensityMatrix tracedMatrix;
      tracedMatrix.allocate(perturbationnewbig.get_leftBlock()->get_ketStateInfo());

      //operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<StackWavefunction&> (iwave)), tracedMatrix, 1.0);
      operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
      int largeNumber = 1000000;
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
      tracedMatrix.deallocate();

      iwave.SaveWavefunctionInfo (perturbationnewbig.get_ketStateInfo(), perturbationnewbig.get_leftBlock()->get_sites(), firstOrderState);
      iwave.deallocate();
      if ( l == 0)
	SaveRotationMatrix (perturbationnewsystem.get_sites(), ketrotatematrix);
      SaveRotationMatrix (perturbationnewsystem.get_sites(), ketrotatematrix, firstOrderState);
      
    }

#ifndef SERIAL
    broadcast(calc, ketrotatematrix, 0);
#endif
        
    perturbationnewenvironment.deallocate();
    perturbationenvironment.removeAdditionalOps();
    perturbationenvironment.deallocate();

    perturbationnewsystem.transform_operators(rotatematrix, ketrotatematrix, false, false);

    StackSpinBlock::store(forward, perturbationnewsystem.get_sites(), perturbationnewsystem, targetState, firstOrderState);

    Stackmem[omprank].deallocate(perturbationsystem.getdata(), (perturbationnewsystem.getdata()-perturbationsystem.getdata())+perturbationnewsystem.memoryUsed());
    perturbationenvironmentdot.deallocate();
    perturbationsystemdot.deallocate();

    //perturbationsystem.clear(); perturbationenvironment.clear(); perturbationnewsystem.clear(); perturbationnewenvironment.clear();

  }


  //<target|O|projectors>
  for (int l=0; l<projectors.size(); l++)
  {
    int originalOutputlevel = dmrginp.outputlevel();
    //dmrginp.setOutputlevel() = -1;

    int perturbationIntegral = 0;
    StackSpinBlock perturbationBig, perturbationsystemdot, perturbationenvironmentdot;
    StackSpinBlock perturbationsystem, perturbationenvironment, perturbationnewsystem, perturbationnewenvironment;
    perturbationsystemdot = StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, false);
    perturbationenvironmentdot = StackSpinBlock(environmentDotStart, environmentDotEnd, perturbationIntegral, false);
    perturbationsystem.set_integralIndex() = perturbationIntegral;
    //StackSpinBlock::restore(forward, systemsites, perturbationsystem, targetState, projectors[l]);
    perturbationsystem.set_sites() = systemsites; 

    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, perturbationsystemdot, perturbationenvironmentdot, 
						 perturbationsystem, perturbationnewsystem, perturbationenvironment, perturbationnewenvironment, 
						 perturbationBig, sweepParams, dot_with_sys, useSlater,
						 perturbationIntegral, targetState, projectors[l]);

    StackSpinBlock perturbationnewbig;
    if (sweepParams.get_onedot() && !dot_with_sys)
    {
      perturbationnewsystem.set_integralIndex() = perturbationsystem.get_integralIndex();
      perturbationnewsystem.initialise_op_array(OVERLAP, false);
      perturbationnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
      perturbationnewsystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, perturbationsystem, perturbationsystemdot);
      InitBlocks::InitBigBlock(perturbationnewsystem, perturbationenvironment, perturbationnewbig); 

      perturbationBig.get_rightBlock()->clear();
      perturbationBig.clear();
    }
    else
      perturbationnewbig = perturbationBig;

    std::vector<Matrix> ketrotatematrix;
    if (mpigetrank() == 0) {
      
      StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), perturbationnewbig.get_leftBlock()->get_stateInfo(), perturbationnewbig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();
      //*********************
      GuessWave::guess_wavefunctions(iwave, e, perturbationnewbig, guesstype, 
				     sweepParams.get_onedot(), projectors[l], true, 0.0);
      
      dmrginp.setOutputlevel() = 10;
      StackDensityMatrix tracedMatrix;
      tracedMatrix.allocate(perturbationnewbig.get_leftBlock()->get_ketStateInfo());

      //operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<StackWavefunction&> (iwave)), tracedMatrix, 1.0);
      operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
      int largeNumber = 1000000;
      if (!mpigetrank())
	double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
      tracedMatrix.deallocate();

      iwave.SaveWavefunctionInfo (perturbationnewbig.get_ketStateInfo(), perturbationnewbig.get_leftBlock()->get_sites(), projectors[l]);
      iwave.deallocate();
      SaveRotationMatrix (perturbationnewsystem.get_sites(), ketrotatematrix, projectors[l]);
    }
    
#ifndef SERIAL
    broadcast(calc, ketrotatematrix, 0);
#endif
    
    perturbationnewenvironment.deallocate();
    perturbationenvironment.removeAdditionalOps();
    perturbationenvironment.deallocate();
    
    perturbationnewsystem.transform_operators(rotatematrix, ketrotatematrix, false, false);

    StackSpinBlock::store(forward, perturbationnewsystem.get_sites(), perturbationnewsystem, targetState, projectors[l]);

    Stackmem[omprank].deallocate(perturbationsystem.getdata(), (perturbationnewsystem.getdata()-perturbationsystem.getdata())+perturbationnewsystem.memoryUsed());
    perturbationenvironmentdot.deallocate();
    perturbationsystemdot.deallocate();


  }

#ifndef SERIAL
    broadcast(calc, rotatematrix, 0);
#endif


  dmrginp.setOutputlevel() = originalOutputlevel;

  pout << "**** STACK MEMORY REMAINING before block renorm***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;
  dmrginp.operrotT -> start();
  newSystem.transform_operators(rotatematrix);
  SaveRotationMatrix (newSystem.get_sites(), rotatematrix, targetState);
  dmrginp.operrotT -> stop();

  if (system.get_sites().size() != 1 || (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted())) {
    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysmem = newSystem.memoryUsed();
    newSystem.moveToNewMemory(system.getdata());
    Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
    //system.clear();
  }

  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  pout << "**** STACK MEMORY REMAINING after block renorm***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << *dmrginp.blockintegrals<<"  "<<*dmrginp.blocksites<<"  "<<*dmrginp.statetensorproduct<<"  "<<*dmrginp.statecollectquanta<<"  "<<*dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build sum block"<<endl;
  

}



double SpinAdapted::SweepResponse::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
					    const bool &restart, const int &restartSize, int targetState, 
					  vector<int>& projectors, vector<int>& baseStates, int correctionVector)
{
  StackSpinBlock system;
  int activeSpaceIntegral = 0;
  std::vector<int> perturbationIntegral(dmrginp.getNumIntegrals()-1,0);
  for (int i=0; i<perturbationIntegral.size(); i++)
    perturbationIntegral[i] = i+1;
  std::vector<double> finalEnergy(1,0.0e100);
  double finalError = 0.;

  if (restart) {
    finalEnergy = sweepParams.get_lowest_energy();
    finalError = sweepParams.get_lowest_error();
  }

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
  pout << endl;
  if (forward)
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl; }
  else
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl; }
  pout << "\t\t\t ============================================================================ " << endl;


  if (dmrginp.get_sweep_type() == PARTIAL) {
    int len = forward? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();
    vector<int> sites(len);  
    if (forward)
      for (int i=0; i<len; i++)
	sites[i] = i;
    else
      for (int i=0; i<len; i++)
	sites[i] = dmrginp.last_site() - len + i;
    
    StackSpinBlock::restore (forward, sites, system, targetState, targetState);

    system.set_twoInt(activeSpaceIntegral);

    for (int l=0; l<projectors.size(); l++)
    {
      StackSpinBlock perturbationSystem;
      perturbationSystem.set_integralIndex() = 0;
      //if (sweepParams.get_sweep_iter() == 0)
      //InitBlocks::InitStartingBlock (perturbationSystem,forward, targetState, projectors[l],
      //			       sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
      //			       restartSize, restart, warmUp, 0);
      //else
      StackSpinBlock::restore (forward, sites, perturbationSystem, targetState, projectors[l]);
      StackSpinBlock::store (forward, system.get_sites(), perturbationSystem, targetState, projectors[l]);
    }
    for (int l=0; l<baseStates.size(); l++)
    {
      StackSpinBlock overlapSystem;
      overlapSystem.set_integralIndex() = l+1;
      //if (sweepParams.get_sweep_iter() == 0)
      //InitBlocks::InitStartingBlock (overlapSystem,forward, targetState, baseStates[l],
      //			       sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
      //			       restartSize, restart, warmUp, perturbationIntegral[l]);
      //else
      StackSpinBlock::restore (forward, sites, overlapSystem, targetState, baseStates[l]);
      StackSpinBlock::store (forward, system.get_sites(), overlapSystem, targetState, baseStates[l]);
    }

  }
  else {
    InitBlocks::InitStartingBlock (system,forward, targetState, targetState,
				   sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				   restartSize, restart, warmUp, activeSpaceIntegral);
  
    for (int l=0; l<projectors.size(); l++)
    {
      StackSpinBlock perturbationSystem;
      perturbationSystem.set_integralIndex() = 0;
      InitBlocks::InitStartingBlock (perturbationSystem,forward, targetState, projectors[l],
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     restartSize, restart, warmUp, 0);
      StackSpinBlock::store (forward, system.get_sites(), perturbationSystem, targetState, projectors[l]);
    }
    
    for (int l=0; l<baseStates.size(); l++)
    {
      StackSpinBlock overlapSystem;
      overlapSystem.set_integralIndex() = l+1;
      InitBlocks::InitStartingBlock (overlapSystem,forward, targetState, baseStates[l],
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     restartSize, restart, warmUp, perturbationIntegral[l]);
      StackSpinBlock::store (forward, system.get_sites(), overlapSystem, targetState, baseStates[l]);
    }
  }

  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  // if restart, just restoring an existing block --
  StackSpinBlock::store (forward, system.get_sites(), system, targetState, targetState);
  sweepParams.savestate(forward, system.get_sites().size());

  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (restart || dmrginp.get_sweep_type() == PARTIAL)
    Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);

 // get_n_iters() returns the number of blocking iterations needed in one sweep
  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); )
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }


      if (sweepParams.get_block_iter() == 0 && sweepParams.get_sweep_iter() == 1)
	sweepParams.set_guesstype() = BASIC;
      else if (sweepParams.get_block_iter() != 0) 
	sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp ) {
	int correctionVector = dmrginp.guessState();
	p2out << "USING state "<<correctionVector<<" as initial guess"<<endl;

	StartUp(sweepParams, system, newSystem, dot_with_sys, 
		targetState, correctionVector, projectors, baseStates);
      }
      else {
	double cE = 0.0;
	if (dmrginp.calc_type() == RESPONSEBW) {
	  //When doing BW perturbation theory in the denominator we have (H0 - E) and not (H0-E0)
	  cE = coreEnergy[system.get_integralIndex()];
	  double e2 = BWPTenergy;
	  coreEnergy[system.get_integralIndex()] = cE - e2;
	  p2out << "\t\t\t BW perturbation theory  "<<cE<<"  "<<e2<<endl; 
	}

	BlockAndDecimate(sweepParams, system, newSystem, warmUp, 
			 dot_with_sys, targetState, projectors, baseStates);

	if (dmrginp.calc_type() == RESPONSEBW) 
	  coreEnergy[system.get_integralIndex()] = cE ;

      }
      
      //Need to substitute by?

      if (!warmUp ){

	//this criteria should work for state average or state specific because the lowest sweep energy is always the lowest of the average
	finalError = max(sweepParams.get_lowest_error(),finalError);
	finalEnergy[0] = min(sweepParams.get_lowest_energy()[0], finalEnergy[0]);
	BWPTenergy = finalEnergy[0];
	pout << "final energy "<<finalEnergy[0]<<"  "<<sweepParams.get_lowest_energy()[0]<<endl;
      }
      
      system = newSystem;
      p2out << system<<endl;
      system.printOperatorSummary();
      
      Sweep::set_dot_with_sys(dot_with_sys, system, sweepParams, forward);

      StackSpinBlock::store (forward, system.get_sites(), system, targetState, targetState);
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
  if (!sweepParams.get_onedot() && !warmUp && dmrginp.get_sweep_type() == FULL) {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
    sweepParams.set_onedot() = true;
    sweepParams.set_env_add() = 0;
    bool dot_with_sys = true;
    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, projectors, baseStates);
    sweepParams.set_onedot() = false;
    sweepParams.set_env_add() = 1;
  }

  //pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  //pout << "\t\t\t Sweep Energy for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  
  if (mpigetrank() == 0)
    printf("\t\t\t M = %6i  state = %4i  Largest Discarded Weight = %8.3e  Sweep Energy = %20.10e \n",sweepParams.get_keep_states(), 0, finalError, finalEnergy[0]);

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  if (mpigetrank()==0)
  {
    pout << "About to write dmrg energy"<<endl;
    std::string efile;
    efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
    
    
    FILE* f = fopen(efile.c_str(), "wb");      
    double e = finalEnergy[0]; //sweepParams.get_lowest_energy()[0]; //instead of the lowest energy of the sweep, we record the last energy of the sweep
    fwrite( &e, 1, sizeof(double), f);
    fclose(f);
  }



  return finalError;
}



void SpinAdapted::SweepResponse::StartUp (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool& dot_with_sys, int targetState, int correctionVector, vector<int>& projectors, vector<int>& baseStates)
{
  //if (dmrginp.outputlevel() > 0)
  mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();

  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() - 1;
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
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];

  StackSpinBlock environment, newEnvironment;

  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)
  bool haveNormOps = dot_with_sys, haveCompOps = !dot_with_sys;

  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  //<target| H |target>
  {
    system.addAdditionalOps();
    //pout <<system<<endl;
    //system.printOperatorSummary();
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, correctionVector, correctionVector, 
				   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);

    //pout << newSystem<<endl;
    //newSystem.printOperatorSummary();
    dmrginp.guessgenT -> stop();
    LoadRotationMatrix (newSystem.get_sites(), brarotateMatrix, correctionVector);
    StackWavefunction targetWave;
    StateInfo s;
    if (mpigetrank() == 0) {
      targetWave.LoadWavefunctionInfo(s, newSystem.get_sites(), correctionVector, true);
      targetWave.SaveWavefunctionInfo(s, newSystem.get_sites(), targetState);
      SaveRotationMatrix (newSystem.get_sites(), brarotateMatrix, targetState);
      targetWave.deallocate();
    }
#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, brarotateMatrix, 0);
#endif


    dmrginp.operrotT -> start();
    newSystem.transform_operators(brarotateMatrix);
    dmrginp.operrotT -> stop();

    //if (system.get_sites().size() != 1 || (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted())) {
    {
      long memoryToFree = newSystem.getdata() - system.getdata();
      long newsysmem = newSystem.memoryUsed();
      newSystem.moveToNewMemory(system.getdata());
      Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
      //system.clear();
    }
  }

  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }
  
  //<target| V |baseStates>
  for(int l=0; l<baseStates.size(); l++) 
  {
    dmrginp.guessgenT -> start();
    int perturbationIntegral = l+1;
    StackSpinBlock perturbationSystemDot, perturbationSystem, perturbationNewSystem;
    perturbationSystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, perturbationSystem, targetState, baseStates[l]);
    perturbationSystem.set_twoInt(perturbationIntegral);
    perturbationSystemDot= StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, false);//singleSiteBlocks[perturbationIntegral][systemDotStart];

    perturbationSystem.addAdditionalOps();
    InitBlocks::InitNewSystemBlock(perturbationSystem, perturbationSystemDot, perturbationNewSystem,
				   targetState, baseStates[l], 
				   sweepParams.get_sys_add(), dmrginp.direct(), perturbationIntegral,
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);

    
    dmrginp.guessgenT -> stop();

    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseStates[l]);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, ketrotateMatrix, 0);
#endif

    SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
    dmrginp.operrotT -> start();
    perturbationNewSystem.transform_operators(brarotateMatrix, ketrotateMatrix, false, false);
    dmrginp.operrotT -> stop();
    StackSpinBlock::store(forward, perturbationNewSystem.get_sites(), perturbationNewSystem, targetState, baseStates[l]);

    if(perturbationSystem.get_sites().size() != 1) 
      Stackmem[omprank].deallocate(perturbationSystem.getdata(), (perturbationNewSystem.getdata()-perturbationSystem.getdata())+perturbationNewSystem.memoryUsed());
  }


  //<target|O|baseState>
  //save the updated overlap spinblock
  for(int l=0; l<projectors.size(); l++) 
  {
    int perturbationIntegral = 0;
    
    StackSpinBlock overlapsystem, overlapsystemDot, overlapnewSystem;
    overlapsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    overlapsystemDot= StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, false);//singleSiteBlocks[perturbationIntegral][systemDotStart];

    overlapnewSystem.set_integralIndex() = perturbationIntegral;
    overlapnewSystem.initialise_op_array(OVERLAP, false);
    overlapnewSystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewSystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);

    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, projectors[l]);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, ketrotateMatrix, 0);
#endif

    dmrginp.operrotT -> start();
    overlapnewSystem.transform_operators(brarotateMatrix, ketrotateMatrix, false, false);
    dmrginp.operrotT -> stop();
    StackSpinBlock::store(forward, overlapnewSystem.get_sites(), overlapnewSystem, targetState, projectors[l]);

    Stackmem[omprank].deallocate(overlapsystem.getdata(), (overlapnewSystem.getdata()-overlapsystem.getdata())+overlapnewSystem.memoryUsed());

  }



  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << *dmrginp.blockintegrals<<"  "<<*dmrginp.blocksites<<"  "<<*dmrginp.statetensorproduct<<"  "<<*dmrginp.statecollectquanta<<"  "<<*dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build sum block"<<endl;

}


void SpinAdapted::SweepResponse::WavefunctionCanonicalize (SweepParams &sweepParams, StackSpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int targetState, vector<int>& projectors, vector<int>& baseStates)
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

  vector<int> sitesenvdot(environmentDotSize+1, 0);
  int index = 0;
  for (int i=min(environmentDotStart, environmentDotEnd); i<max(environmentDotStart, environmentDotEnd)+1; i++) {
    sitesenvdot[index] = (i);
    index++;
  }

  //environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][sitesenvdot[0]];
  StackSpinBlock::restore(!forward, sitesenvdot, environmentDot, targetState, targetState); 

  StackSpinBlock environment, newEnvironment;
  
  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  system.addAdditionalOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, targetState, targetState, 
				 sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				 DISTRIBUTED_STORAGE, false, true);

  newSystem.set_loopblock(false);  environmentDot.set_loopblock(false); 
  InitBlocks::InitBigBlock(newSystem, environmentDot, big);

  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;
  
  
  std::vector<StackWavefunction> lowerStates;
  
  
  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
  DiagonalMatrix e;
  guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
  
  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }
  
  
  StackWavefunction targetWave; targetWave.initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); targetWave.Clear();
  pout << "transform previous"<<endl;

  //**************************
  if (!mpigetrank())
    GuessWave::transform_previous_twodot_to_onedot_wavefunction(targetWave, big, targetState);

  targetWave.set_onedot(true);

  
  StackDensityMatrix bratracedMatrix;
  bratracedMatrix.allocate(newSystem.get_braStateInfo());

  //operatorfunctions::MultiplyProduct(targetWave, Transpose(const_cast<StackWavefunction&> (targetWave)), bratracedMatrix, 1.0);
  operatorfunctions::MultiplyWithOwnTranspose (targetWave, bratracedMatrix, 1.0);  
  int largeNumber = 1000000;
  if (!mpigetrank())
    double error = makeRotateMatrix(bratracedMatrix, rotatematrix, largeNumber, sweepParams.get_keep_qstates());
  bratracedMatrix.deallocate();
  pout << "broadcast rotate"<<endl;
#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, rotatematrix, 0);
#endif

  
  targetWave.SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), targetState);
  targetWave.deallocate();

  SaveRotationMatrix (newSystem.get_sites(), rotatematrix, targetState);
  
  newSystem.transform_operators(rotatematrix);
  long memoryToFree = newSystem.getdata() - system.getdata();
  long newsysmem = newSystem.memoryUsed();
  newSystem.moveToNewMemory(system.getdata());
  Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);


  dmrginp.setOutputlevel() = originalOutputlevel;
  
  
  
  
  
  
  //<target|O|firstOrderState>
  //save the updated overlap spinblock
  for (int l=0; l<baseStates.size(); l++)
  {
    int perturbationIntegral = l+1;
    
    StackSpinBlock overlapBig;
    StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment, overlapenvironmentDot;
    StackSpinBlock overlapsystemDot(systemDotStart, systemDotEnd, perturbationIntegral, false);

    overlapenvironmentDot.set_integralIndex() = perturbationIntegral;

    //overlapenvironmentDot = StackSpinBlock(sitesenvdot[0], sitesenvdot[0], perturbationIntegral, false);//singleSiteBlocks[perturbationIntegral][sitesenvdot[0]];
    StackSpinBlock::restore(!forward, sitesenvdot, overlapenvironmentDot, targetState, baseStates[l]); 

    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, overlapsystem, targetState, baseStates[l]);
    overlapsystem.set_twoInt(perturbationIntegral);
    overlapsystem.addAdditionalOps();
    InitBlocks::InitNewSystemBlock(overlapsystem, overlapsystemDot, 
				   overlapnewsystem, targetState, baseStates[l], 
				   overlapsystemDot.size(), dmrginp.direct(), 
				   perturbationIntegral, DISTRIBUTED_STORAGE, false, true);

    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    pout << "transform iwave"<<endl;
    StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();
    //****************************

    if (!mpigetrank())
      GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, baseStates[l]);
    iwave.set_onedot(true);

    std::vector<Matrix> ketrotatematrix;
    StackDensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    //tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
    //operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<StackWavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    tracedMatrix.deallocate();

#ifndef SERIAL
    broadcast(calc, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), baseStates[l]);
    iwave.deallocate();
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, baseStates[l]);

    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix, false, false);
    StackSpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, baseStates[l]);

    if(overlapsystem.get_sites().size() != 1) 
      Stackmem[omprank].deallocate(overlapsystem.getdata(), (overlapnewsystem.getdata()-overlapsystem.getdata())+overlapnewsystem.memoryUsed());

    overlapenvironmentDot.deallocate();
    overlapsystemDot.deallocate();

    //overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }


  //<target|O|projectors>
  //save the updated overlap spinblock
  for (int l=0; l<projectors.size(); l++)
  {
    int perturbationIntegral = 0;
    
    StackSpinBlock overlapBig;
    StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment, overlapenvironmentDot;
    StackSpinBlock overlapsystemDot(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapenvironmentDot.set_integralIndex() = perturbationIntegral;
    //overlapenvironmentDot = StackSpinBlock(sitesenvdot[0], sitesenvdot[0], perturbationIntegral, false);//singleSiteBlocks[perturbationIntegral][sitesenvdot[0]];
    StackSpinBlock::restore(!forward, sitesenvdot, overlapenvironmentDot, targetState, projectors[l]); 

    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    overlapnewsystem.set_integralIndex() = perturbationIntegral;
    overlapnewsystem.initialise_op_array(OVERLAP, false);
    overlapnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewsystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);

    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    pout << "transform iwave"<<endl;
    StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();
    //******************************
    if (!mpigetrank())
      GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, projectors[l]);
    iwave.set_onedot(true);


    std::vector<Matrix> ketrotatematrix;
    StackDensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());

    //tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
    //operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<StackWavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    tracedMatrix.deallocate();

#ifndef SERIAL
    broadcast(calc, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), projectors[l]);
    iwave.deallocate();

    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, projectors[l]);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix, false, false);
    StackSpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, projectors[l]);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();

    if(overlapsystem.get_sites().size() != 1) 
      Stackmem[omprank].deallocate(overlapsystem.getdata(), (overlapnewsystem.getdata()-overlapsystem.getdata())+overlapnewsystem.memoryUsed());
    overlapenvironmentDot.deallocate();
    overlapsystemDot.deallocate();
  }

  dmrginp.setOutputlevel() = originalOutputlevel;
  
  pout << newSystem<<endl;
  StackSpinBlock::store(forward, newSystem.get_sites(), newSystem, targetState, targetState);
  
  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");
  
  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
}

/*

void SpinAdapted::SweepResponse::StartUp (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool& dot_with_sys, int targetState, int correctionVector, vector<int>& projectors, vector<int>& baseStates)
{
  //if (dmrginp.outputlevel() > 0)
  mcheck("at the start of block and decimate");
  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();

  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot, environmentDot;
  int systemDotStart, systemDotEnd, environmentDotStart, environmentDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() - 1;
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
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);
  StackSpinBlock environment, newEnvironment;

  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)
  bool haveNormOps = dot_with_sys, haveCompOps = true;

  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  //<target| H |target>
  {
    system.addAdditionalOps();
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, correctionVector, correctionVector, 
				   sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);
    
    dmrginp.guessgenT -> stop();
    LoadRotationMatrix (newSystem.get_sites(), brarotateMatrix, correctionVector);
    StackWavefunction targetWave;
    StateInfo s;
    if (mpigetrank() == 0) {
      targetWave.LoadWavefunctionInfo(s, newSystem.get_sites(), correctionVector, true);
      targetWave.SaveWavefunctionInfo(s, newSystem.get_sites(), targetState);
      SaveRotationMatrix (newSystem.get_sites(), brarotateMatrix, targetState);
      targetWave.deallocate();
    }
#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, brarotateMatrix, 0);
#endif
    
    dmrginp.operrotT -> start();
    newSystem.transform_operators(brarotateMatrix);
    dmrginp.operrotT -> stop();
  }

  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }
  
  //<target| V |baseStates>
  for(int l=0; l<baseStates.size(); l++) 
  {
    dmrginp.guessgenT -> start();
    int perturbationIntegral = l+1;
    StackSpinBlock perturbationSystemDot, perturbationSystem, perturbationNewSystem;
    perturbationSystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, perturbationSystem, targetState, baseStates[l]);
    perturbationSystemDot = StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);

    perturbationSystem.addAdditionalOps();
    InitBlocks::InitNewSystemBlock(perturbationSystem, perturbationSystemDot, perturbationNewSystem,
				   targetState, baseStates[l], 
				   sweepParams.get_sys_add(), dmrginp.direct(), perturbationIntegral,
				   DISTRIBUTED_STORAGE, haveNormOps, haveCompOps);

    
    dmrginp.guessgenT -> stop();

    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseStates[l]);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, ketrotateMatrix, 0);
#endif

    dmrginp.operrotT -> start();
    perturbationNewSystem.transform_operators(brarotateMatrix, ketrotateMatrix);
    dmrginp.operrotT -> stop();
    StackSpinBlock::store(forward, perturbationNewSystem.get_sites(), perturbationNewSystem, targetState, baseStates[l]);
  }


  //<target|O|baseState>
  //save the updated overlap spinblock
  for(int l=0; l<projectors.size(); l++) 
  {
    int perturbationIntegral = 0;
    
    StackSpinBlock overlapsystem, overlapsystemDot, overlapnewSystem;
    overlapsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    overlapsystemDot = StackSpinBlock(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapnewSystem.set_integralIndex() = perturbationIntegral;
    overlapnewSystem.initialise_op_array(OVERLAP, false);
    overlapnewSystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewSystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);

    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, projectors[l]);
    
#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, ketrotateMatrix, 0);
#endif

    dmrginp.operrotT -> start();
    overlapnewSystem.transform_operators(brarotateMatrix, ketrotateMatrix);
    dmrginp.operrotT -> stop();
    StackSpinBlock::store(forward, overlapnewSystem.get_sites(), overlapnewSystem, targetState, projectors[l]);

  }



  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << *dmrginp.blockintegrals<<"  "<<*dmrginp.blocksites<<"  "<<*dmrginp.statetensorproduct<<"  "<<*dmrginp.statecollectquanta<<"  "<<*dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build sum block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;

}


void SpinAdapted::SweepResponse::WavefunctionCanonicalize (SweepParams &sweepParams, StackSpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int targetState, vector<int>& projectors, vector<int>& baseStates)
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
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);

  vector<int> sitesenvdot(environmentDotSize+1, 0);
  int index = 0;
  for (int i=min(environmentDotStart, environmentDotEnd); i<max(environmentDotStart, environmentDotEnd)+1; i++) {
    sitesenvdot[index] = (i);
    index++;
  }

  StackSpinBlock::restore(!forward, sitesenvdot, environmentDot, targetState, targetState); 

  StackSpinBlock environment, newEnvironment;
  
  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  system.addAdditionalOps();
  InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, targetState, targetState, sweepParams.get_sys_add(), dmrginp.direct(), system.get_integralIndex(), 
				 DISTRIBUTED_STORAGE, false, true);

  newSystem.set_loopblock(false);  environmentDot.set_loopblock(false); 
  InitBlocks::InitBigBlock(newSystem, environmentDot, big);

  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;
  
  if (dmrginp.outputlevel() > 0)
    mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot()) pout << "\t\t\t System  Block"<<system;    
  else pout << "\t\t\t System  Block"<<newSystem;
  pout << "\t\t\t Environment Block"<<newEnvironment<<endl;
  p1out << "\t\t\t Solving wavefunction "<<endl;
  
  std::vector<StackWavefunction> lowerStates;
  
  
  //make the baseState
  int originalOutputlevel = dmrginp.outputlevel();
  dmrginp.setOutputlevel() = -1;
  DiagonalMatrix e;
  guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
  
  std::vector<int> systemsites;
  if (dmrginp.spinAdapted())
    systemsites = system.get_sites();
  else {
    if (forward) {
      systemsites.push_back(0); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
    else {
      systemsites.push_back(system.get_sites()[0]/2); systemsites.push_back(*system.get_sites().rbegin()/2);
    }
  }
  
  
  StackWavefunction targetWave; targetWave.initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); targetWave.Clear();
  pout << "transform previous"<<endl;


  //if (!mpigetrank())
  //GuessWave::transform_previous_twodot_to_onedot_wavefunction(targetWave, big, targetState);

  targetWave.set_onedot(true);

  
  StackDensityMatrix bratracedMatrix;
  bratracedMatrix.allocate(newSystem.get_braStateInfo());

  //operatorfunctions::MultiplyProduct(targetWave, Transpose(const_cast<StackWavefunction&> (targetWave)), bratracedMatrix, 1.0);
  operatorfunctions::MultiplyWithOwnTranspose (targetWave, bratracedMatrix, 1.0);  
  int largeNumber = 1000000;
  if (!mpigetrank())
    double error = makeRotateMatrix(bratracedMatrix, rotatematrix, largeNumber, sweepParams.get_keep_qstates());
  bratracedMatrix.deallocate();
  pout << "broadcast rotate"<<endl;
#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, rotatematrix, 0);
#endif

  
  targetWave.SaveWavefunctionInfo (big.get_braStateInfo(), big.get_leftBlock()->get_sites(), targetState);
  targetWave.deallocate();
  SaveRotationMatrix (newSystem.get_sites(), rotatematrix, targetState);
  
  newSystem.transform_operators(rotatematrix);
  
  
  dmrginp.setOutputlevel() = originalOutputlevel;
  
  
  
  
  
  
  //<target|O|firstOrderState>
  //save the updated overlap spinblock
  for (int l=0; l<baseStates.size(); l++)
  {
    int perturbationIntegral = l+1;
    
    StackSpinBlock overlapBig;
    StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment, overlapenvironmentDot;
    StackSpinBlock overlapsystemDot(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapenvironmentDot.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(!forward, sitesenvdot, overlapenvironmentDot, targetState, baseStates[l]); 

    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, overlapsystem, targetState, baseStates[l]);
    overlapsystem.addAdditionalOps();
    InitBlocks::InitNewSystemBlock(overlapsystem, overlapsystemDot, 
				   overlapnewsystem, targetState, baseStates[l], 
				   overlapsystemDot.size(), dmrginp.direct(), 
				   perturbationIntegral, DISTRIBUTED_STORAGE, false, true);

    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    pout << "transform iwave"<<endl;
    StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();

    //if (!mpigetrank())
    //GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, baseStates[l]);
    iwave.set_onedot(true);

    std::vector<Matrix> ketrotatematrix;
    StackDensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    //tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
    //operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<StackWavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    tracedMatrix.deallocate();

#ifndef SERIAL
    broadcast(calc, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), baseStates[l]);
    iwave.deallocate();
    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, baseStates[l]);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    StackSpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, baseStates[l]);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }


  //<target|O|projectors>
  //save the updated overlap spinblock
  for (int l=0; l<projectors.size(); l++)
  {
    int perturbationIntegral = 0;
    
    StackSpinBlock overlapBig;
    StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment, overlapenvironmentDot;
    StackSpinBlock overlapsystemDot(systemDotStart, systemDotEnd, perturbationIntegral, true);

    overlapenvironmentDot.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(!forward, sitesenvdot, overlapenvironmentDot, targetState, projectors[l]); 

    guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
    
    DiagonalMatrix e;

    overlapsystem.set_integralIndex() = perturbationIntegral;
    StackSpinBlock::restore(forward, systemsites, overlapsystem, targetState, projectors[l]);
    overlapnewsystem.set_integralIndex() = perturbationIntegral;
    overlapnewsystem.initialise_op_array(OVERLAP, false);
    overlapnewsystem.setstoragetype(DISTRIBUTED_STORAGE);
    overlapnewsystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, overlapsystem, overlapsystemDot);

    overlapnewsystem.set_loopblock(false);  overlapenvironmentDot.set_loopblock(false); 
    InitBlocks::InitBigBlock(overlapnewsystem, overlapenvironmentDot, overlapBig);

    
    pout << "transform iwave"<<endl;
    StackWavefunction iwave; iwave.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot()); iwave.Clear();

    //if (!mpigetrank())
    //GuessWave::transform_previous_twodot_to_onedot_wavefunction(iwave, overlapBig, projectors[l]);
    iwave.set_onedot(true);


    std::vector<Matrix> ketrotatematrix;
    StackDensityMatrix tracedMatrix;
    tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());

    //tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());
    operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
    //operatorfunctions::MultiplyProduct(iwave, Transpose(const_cast<StackWavefunction&> (iwave)), tracedMatrix, 1.0);
    int largeNumber = 1000000;
    if (!mpigetrank())
      double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());
    tracedMatrix.deallocate();

#ifndef SERIAL
    broadcast(calc, ketrotatematrix, 0);
#endif
    
    iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), projectors[l]);
    iwave.deallocate();

    SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, projectors[l]);
    
    overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix);
    StackSpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, targetState, projectors[l]);
    overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
    
  }

  dmrginp.setOutputlevel() = originalOutputlevel;
  
  pout << newSystem<<endl;
  StackSpinBlock::store(forward, newSystem.get_sites(), newSystem, targetState, targetState);
  
  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");
  
  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop   S_2_opxop"<<endl;
  p3out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<" "<<*dmrginp.s2time<<endl;
}
*/
