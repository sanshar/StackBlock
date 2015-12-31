/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "Stackspinblock.h"
#include <boost/bind.hpp>
#include <boost/functional.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include "solver.h"
#include "operatorloops.h"
#include <numeric>
#include "rotationmat.h"
#include "Stackdensity.h"
#include "initblocks.h"
#include "stackguess_wavefunction.h"
#include "linear.h"
#include "davidson.h"
#include <stdlib.h>
//#include "diis.h"
#include "Stackwavefunction.h"
#include "Stackspinblock.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

#include "pario.h"

using namespace boost;
using namespace std;

namespace SpinAdapted{
void StackSpinBlock::RenormaliseFrom(vector<double> &energies, vector<double> &spins, double& error, 
				vector<Matrix>& rotateMatrix, const int keptstates, 
				const int keptqstates, const double tol, StackSpinBlock& big, 
				const guessWaveTypes &guesswavetype, const double noise, 
				const double additional_noise, const bool &onedot, StackSpinBlock& System,
				StackSpinBlock& sysDot, StackSpinBlock& environment, const bool& dot_with_sys,
				const bool& warmUp, int sweepiter, int currentRoot, 
				std::vector<StackWavefunction>& lowerStates, StackDensityMatrix* ReducedDM)
{
  dmrginp.davidsonT -> start();
  int nroots = dmrginp.setStateSpecific() ? 1 : dmrginp.nroots(sweepiter);
  vector<StackWavefunction> wave_solutions(nroots);

  wave_solutions[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), onedot);
  wave_solutions[0].Clear();

  if (mpigetrank() == 0) {
    for (int i=1; i<nroots; i++) {
      wave_solutions[i].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), onedot);
      wave_solutions[i].Clear();
    }
  }

  if (dmrginp.outputlevel() > 0)
    mcheck("before davidson but after all blocks are built");

  dmrginp.solvewf -> start();

  SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));

  Solver::solve_wavefunction(wave_solutions, energies, big, tol, guesswavetype, onedot, 
			     dot_with_sys, warmUp, false, additional_noise, currentRoot, lowerStates);
  dmrginp.solvewf -> stop();

  StackSpinBlock newsystem;
  StackSpinBlock newenvironment;
  StackSpinBlock newbig;
  dmrginp.postwfrearrange -> start();

  if (onedot && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(System, sysDot, newsystem, currentRoot, currentRoot, 
				   sysDot.size(), dmrginp.direct(), System.get_integralIndex(), DISTRIBUTED_STORAGE, false, true);
    InitBlocks::InitBigBlock(newsystem, environment, newbig); 
    for (int i=0; i<nroots&& mpigetrank()==0; i++) 
    {
      StackWavefunction tempwave; tempwave.initialise(wave_solutions[i]);
      GuessWave::onedot_shufflesysdot(big.get_stateInfo(), newbig.get_stateInfo(),wave_solutions[i], 
				      tempwave);  
      DCOPY(wave_solutions[i].memoryUsed(), tempwave.get_data(), 1, wave_solutions[i].get_data(), 1);
      wave_solutions[i].initialise(wave_solutions[i].get_deltaQuantum(), newbig.get_leftBlock()->get_stateInfo(), newbig.get_rightBlock()->get_stateInfo(), wave_solutions[i].get_onedot(), wave_solutions[i].get_data(), wave_solutions[i].memoryUsed());
      tempwave.deallocate();
    }

#ifndef SERIAL
    mpi::communicator world;
    broadcast(calc, wave_solutions[0], 0);
    if (mpigetrank() != 0)
      wave_solutions[0].allocateOperatorMatrix();
#endif

    *this = newsystem;
    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;
  dmrginp.postwfrearrange -> stop();

  if (dmrginp.outputlevel() > 0)
    mcheck("after davidson before noise");

  dmrginp.davidsonT -> stop();

  dmrginp.rotmatrixT -> start();
  StackDensityMatrix tracedMatrix(braStateInfo);
  tracedMatrix.allocate(braStateInfo);

  //mcheck("after allocating density matrix");

  bool normalnoise = warmUp;
  if (newbig.get_rightBlock()->size() < 2)
    normalnoise = true;
  
  dmrginp.addnoise -> start();
  double twodotnoise = 0.0;
  if (dmrginp.noise_type() == RANDOM)
    twodotnoise = additional_noise;

  //************************
  tracedMatrix.makedensitymatrix(wave_solutions, newbig, dmrginp.weights(sweepiter), noise, twodotnoise, normalnoise);

  if (ReducedDM != 0 && mpigetrank() == 0) {
    DCOPY(ReducedDM->memoryUsed(), tracedMatrix.get_data(), 1, ReducedDM->get_data(), 1);
  } 

  dmrginp.addnoise -> stop();

  if (!mpigetrank())
    error = makeRotateMatrix(tracedMatrix, rotateMatrix, keptstates, keptqstates);

  tracedMatrix.deallocate();

#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, rotateMatrix, 0);
#endif

  SaveRotationMatrix (newbig.leftBlock->sites, rotateMatrix);
  for (int i=0; i<nroots && mpigetrank() == 0; i++) {
    int state = dmrginp.setStateSpecific() ? currentRoot : i;
    if (dmrginp.solve_method() == CONJUGATE_GRADIENT) 
      state = currentRoot;

    SaveRotationMatrix (newbig.leftBlock->sites, rotateMatrix, state);
    wave_solutions[i].SaveWavefunctionInfo (newbig.braStateInfo, newbig.leftBlock->sites, state);
  }

  if (mpigetrank() == 0) {
    for (int i=nroots-1; i>0; i--)
      wave_solutions[i].deallocate();
  }
  wave_solutions[0].deallocate();
  dmrginp.rotmatrixT -> stop();
  //if (dmrginp.outputlevel() > 0)
  //mcheck("after noise and calculation of density matrix");
}

  
double makeRotateMatrix(StackDensityMatrix& tracedMatrix, vector<Matrix>& rotateMatrix, const int& keptstates, const int& keptqstates, vector<DiagonalMatrix> *eigs)
{

  std::vector<DiagonalMatrix> eigenMatrix;
  if (dmrginp.hamiltonian() == BCS)
    svd_densitymat(tracedMatrix, eigenMatrix);
  else
    diagonalise_dm(tracedMatrix, eigenMatrix);

  if (eigs != 0)
    *eigs = eigenMatrix;

  vector<pair<int, int> > inorderwts;
  vector<vector<int> > wtsbyquanta;
  

  sort_weights(eigenMatrix, inorderwts, wtsbyquanta);
  
  
  // make transformation matrix by various algorithms
  int totalstatesbydm = min(static_cast<int>(inorderwts.size()), keptstates);
  int totalstatesbyquanta = min(static_cast<int>(inorderwts.size()), keptstates + keptqstates) - totalstatesbydm;
  if (totalstatesbyquanta < 0) totalstatesbyquanta = 0;
  
  p3out << "\t\t\t total states using dm and quanta " << totalstatesbydm << " " << totalstatesbyquanta << endl;
  
  return assign_matrix_by_dm(rotateMatrix, eigenMatrix, tracedMatrix, inorderwts, wtsbyquanta, totalstatesbydm, 
			     totalstatesbyquanta, 0, 0);
}

}
