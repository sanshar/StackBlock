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



void SpinAdapted::Sweep::BlockAndDecimate (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys)
{
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
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];
  StackSpinBlock environment, newEnvironment;
  StackSpinBlock big;  // new_sys = sys+sys_dot; new_env = env+env_dot; big = new_sys + new_env then renormalize to find new_sys(new)

  makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, 
				 big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), 
				 sweepParams.current_root(), sweepParams.current_root());

  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  //if (dmrginp.outputlevel() > 0)
  //mcheck(""); 
  if (!dot_with_sys && sweepParams.get_onedot()) {
    pout << "\t\t\t System  Block"<<system;    
    system.printOperatorSummary();
  }
  else {
    pout << "\t\t\t System  Block"<<newSystem;
    newSystem.printOperatorSummary();
  }
  pout << endl<<"\t\t\t Environment Block"<<newEnvironment<<endl;
  newEnvironment.printOperatorSummary();

  p1out << "\t\t\t Solving wavefunction "<<endl;

  std::vector<StackWavefunction> lowerStates; 

  if(sweepParams.current_root() >= 0 ) {
    if (mpigetrank() == 0) {
      lowerStates.resize(sweepParams.current_root());
      for (int istate=0; istate<sweepParams.current_root(); istate++)
	lowerStates[istate].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), sweepParams.get_onedot());
    }

    int originalOutputlevel = dmrginp.outputlevel();
    dmrginp.setOutputlevel() = -1;

    DiagonalMatrix e;
    for (int istate = 0; istate<sweepParams.current_root(); istate++) {
      guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;

      //now one needs to make |phi> = O|psi> so that the |phi> has the same dimensions as our target state
      StackSpinBlock overlapBig;
      StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
      makeSystemEnvironmentBigOverlapBlocks(system.get_sites(), systemDot, environmentDot,
					    overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment,
					    overlapBig, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), 
					    sweepParams.current_root(), istate);


      if (mpigetrank() == 0) {
	lowerStates[istate].Clear();

	StackWavefunction temp; temp.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), true);
	temp.Clear();

	//************************
	GuessWave::guess_wavefunctions(temp, e, overlapBig, guesstype, sweepParams.get_onedot(), istate, dot_with_sys, 0.0);
	overlapBig.multiplyOverlap(temp, &lowerStates[istate], MAX_THRD);
	temp.deallocate();
      }
      //overlapsystem.clear(); overlapenvironment.clear(); overlapnewsystem.clear(); overlapnewenvironment.clear();
      overlapnewenvironment.deallocate();
      overlapenvironment.deallocate();
      overlapnewsystem.deallocate();
      overlapsystem.deallocate();
      

    }
    dmrginp.setOutputlevel() = originalOutputlevel;
  }

  double Noise = sweepParams.get_noise(), Additionalnoise = sweepParams.get_additional_noise();
  if (find(dmrginp.get_openorbs().begin(), dmrginp.get_openorbs().end(), systemDotStart) != dmrginp.get_openorbs().end()) {
    Noise = 0.0; Additionalnoise = 0.0;
  }

  newSystem.RenormaliseFrom (sweepParams.set_lowest_energy(), sweepParams.set_lowest_energy_spins(), sweepParams.set_lowest_error(), 
                             rotatematrix, sweepParams.get_keep_states(), 
                             sweepParams.get_keep_qstates(), sweepParams.get_davidson_tol(), big, sweepParams.get_guesstype(), Noise, 
                             Additionalnoise, sweepParams.get_onedot(), system, systemDot, environment, 
			     dot_with_sys, useSlater, sweepParams.get_sweep_iter(), sweepParams.current_root(), lowerStates);

  if (mpigetrank() == 0 && sweepParams.current_root() >= 0 ) 
    for (int istate = sweepParams.current_root()-1; istate>-1; istate--) 
      lowerStates[istate].deallocate();


  //newEnvironment.clear();
  newEnvironment.removeAdditionalOps();
  newEnvironment.deallocate();
  environment.removeAdditionalOps();
  //environment.clear();
  environment.deallocate();


  p1out <<"\t\t\t Performing Renormalization "<<endl;
  pout << "\n\t\t\t Total discarded weight "<<sweepParams.get_lowest_error()<<endl<<endl;


  dmrginp.multiplierT -> stop();
  dmrginp.operrotT -> start();
  newSystem.transform_operators(rotatematrix);
  SpinAdapted::SpinQuantum hq(0,SpinAdapted::SpinSpace(0),SpinAdapted::IrrepSpace(0));

  //if (system.get_sites().size() != 1 || (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted())) {
  {
    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysmem = newSystem.memoryUsed();
    newSystem.moveToNewMemory(system.getdata());
    Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
    //system.clear();
  }



  //save the updated overlap spinblock
  if( sweepParams.current_root() >= 0 ) {
    int originalOutputlevel = dmrginp.outputlevel();
    dmrginp.setOutputlevel() = -1;
    for (int istate = 0; istate<sweepParams.current_root(); istate++) {
      StackSpinBlock overlapBig;
      StackSpinBlock overlapsystem, overlapenvironment, overlapnewsystem, overlapnewenvironment;
      StackSpinBlock overlapsystemDot= StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];
      StackSpinBlock overlapenvironmentDot=StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), true);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];
      guessWaveTypes guesstype = sweepParams.get_block_iter() == 0 ? TRANSPOSE : TRANSFORM;
      
      DiagonalMatrix e;
      makeSystemEnvironmentBigOverlapBlocks(system.get_sites(), overlapsystemDot, overlapenvironmentDot,
					    overlapsystem, overlapnewsystem, overlapenvironment, overlapnewenvironment,
					    overlapBig, sweepParams, true, useSlater, newSystem.get_integralIndex(), 
					    sweepParams.current_root(), istate);

      StackWavefunction iwave;
      iwave.initialise(dmrginp.effective_molecule_quantum_vec(), overlapBig.get_leftBlock()->get_stateInfo(), overlapBig.get_rightBlock()->get_stateInfo(), true);
      iwave.Clear();

      //*********************************
      GuessWave::guess_wavefunctions(iwave, e, overlapBig, guesstype, sweepParams.get_onedot(), istate, true, 0.0);
      std::vector<Matrix> ketrotatematrix;
      StackDensityMatrix tracedMatrix;
      tracedMatrix.allocate(overlapnewsystem.get_ketStateInfo());

      operatorfunctions::MultiplyWithOwnTranspose (iwave, tracedMatrix, 1.0);  
      int largeNumber = 1000000;
      if (!mpigetrank())
	double error = makeRotateMatrix(tracedMatrix, ketrotatematrix, largeNumber, sweepParams.get_keep_qstates());

#ifndef SERIAL
      mpi::communicator world;
      broadcast(calc, ketrotatematrix, 0);
#endif
      tracedMatrix.deallocate();
      iwave.SaveWavefunctionInfo (overlapBig.get_ketStateInfo(), overlapBig.get_leftBlock()->get_sites(), istate);
      SaveRotationMatrix (overlapnewsystem.get_sites(), ketrotatematrix, istate);
      iwave.deallocate();

      overlapnewenvironment.deallocate();
      overlapenvironment.deallocate();

      overlapnewsystem.transform_operators(rotatematrix, ketrotatematrix, false, false);
      StackSpinBlock::store(forward, overlapnewsystem.get_sites(), overlapnewsystem, sweepParams.current_root(), istate);

      Stackmem[omprank].deallocate(overlapsystem.getdata(), (overlapnewsystem.getdata()-overlapsystem.getdata())+overlapnewsystem.memoryUsed());
      overlapenvironmentDot.deallocate();
      overlapsystemDot.deallocate();
      
    }
    dmrginp.setOutputlevel() = originalOutputlevel;
  }
  dmrginp.operrotT -> stop();


  p2out << str(boost::format("%-40s - %-10.4f\n") % "Total walltime" % globaltimer.totalwalltime());
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Blocking (includes first sweep)" % *(dmrginp.guessgenT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->diski" % *(dmrginp.diski));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->makeiter" % *(dmrginp.readmakeiter));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->allocop" % *(dmrginp.readallocatemem));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->rawdata" % *(dmrginp.rawdatai));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->mpicomm" % *(dmrginp.datatransfer));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->builditerators" % *(dmrginp.builditeratorsT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Wavefunction Solution" % *(dmrginp.multiplierT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->davidson/guesswf/diagonal" % *(dmrginp.davidsonT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->guesswf" % *(dmrginp.guesswf));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->diagonal" % *(dmrginp.makediagonal));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->davidson" % *(dmrginp.blockdavid));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->makerotation(includes noise)" % *(dmrginp.rotmatrixT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "          |-->Add noise" % *(dmrginp.rotmatrixT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->wave and rotation io" % *(dmrginp.diskwo));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Renormalisation" % *(dmrginp.operrotT));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->in parallel region" % *(dmrginp.parallelrenorm));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Save block" % *(dmrginp.disko)); 
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->rawdata" % *(dmrginp.rawdatao)); 
}


