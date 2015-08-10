/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "StateInfo.h"
#include "solver.h"
#include "linear.h"
#include "davidson.h"
#include "stackguess_wavefunction.h"
#include "blas_calls.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"
#include "StackBaseOperator.h"
#include "Stackspinblock.h"
#include "Stackwavefunction.h"
#include <newmat.h>


void SpinAdapted::Solver::buildDiagonal(StackSpinBlock& big, DiagonalMatrix& e)
{
  p1out << "\t\t\t Building Diagonal Hamiltonian " << endl;
  big.diagonalH(e);
  p1out << "\t\t\t Done building diagonal hamiltonian "<<endl;
  FORTINT m, n=1, nsize=e.Storage();
  p2out << "\t\t\t Number of elements in wavefunction :: " << e.Ncols() << "  "<<endl;
  if (mpigetrank()==0) {
    m = idamax_(nsize,e.Store(), n); 
    p3out << "\t\t\t highest diagonal value "<<m<<" "<<e(m)<<endl;
  }
  else 
    e.ReSize(0);
}

void SpinAdapted::Solver::unCollectSolutions(vector<StackWavefunction>& solution, int nroots, StateInfo& leftState, StateInfo& rightState)
{
  //uncollectquanta in all the solutions
  StateInfo l = leftState, r = rightState;
  if (leftState.rightStateInfo != 0) {
    solution[0].UnCollectQuantaAlongRows(l, r);
    l = *leftState.unCollectedStateInfo;
  }
  if (rightState.rightStateInfo != 0)
    solution[0].UnCollectQuantaAlongColumns(l, r);

  for (int i=1; i<nroots; i++) {
    l = leftState; r = rightState;
    if (mpigetrank() == 0) {
      if (leftState.rightStateInfo != 0) {
	solution[i].UnCollectQuantaAlongRows(l, r);
	l = *leftState.unCollectedStateInfo;
      }
      if (rightState.rightStateInfo != 0)
	solution[i].UnCollectQuantaAlongColumns(l, r);
    }
  }
}

void SpinAdapted::Solver::collectSolutions(vector<StackWavefunction>& solution, int nroots, StateInfo& leftState, StateInfo& rightState)
{
  for (int i=0; i<nroots; i++) {
    StateInfo left = leftState, right = rightState;
    if (mpigetrank() == 0) {
      if (left.rightStateInfo != 0) {
	solution[i].CollectQuantaAlongRows(left, right);
	left.CollectQuanta();
      }
      if (right.rightStateInfo != 0)
	solution[i].CollectQuantaAlongColumns(left, right);
    }
  }
}

void SpinAdapted::Solver::collectBlocks(StackSpinBlock& big)
{
  if (big.get_leftBlock()->get_rightBlock() != 0) {
    big.get_leftBlock()->collectQuanta();
    big.get_leftBlock()->CleanUpOperators();
  } 
  if (big.get_rightBlock()->get_rightBlock() != 0) {
    big.get_rightBlock()->collectQuanta();
    big.get_rightBlock()->CleanUpOperators();
  }
  big.recreateStateInfo(PARTICLE_SPIN_NUMBER_CONSTRAINT);
}

void SpinAdapted::Solver::unCollectBlocks(StackSpinBlock& big)
{
  //uncollectquanta in all the blocks
  big.get_leftBlock()->recreateStateInfo(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT); 
  big.get_rightBlock()->recreateStateInfo(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

  if (big.get_leftBlock()->get_rightBlock() != 0)
    big.get_leftBlock()->CleanUpOperators(); 
  if (big.get_rightBlock()->get_rightBlock() != 0)
    big.get_rightBlock()->CleanUpOperators();

  big.recreateStateInfo(PARTICLE_SPIN_NUMBER_CONSTRAINT);
}

void SpinAdapted::Solver::solve_wavefunction(vector<StackWavefunction>& solution, vector<double>& energies, StackSpinBlock& big, const double tol, 
					     const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp,
					     double additional_noise, int currentRoot, std::vector<StackWavefunction>& lowerStates)
{
  const int nroots = dmrginp.setStateSpecific() ? 1 : dmrginp.nroots();

  DiagonalMatrix e;
  bool useprecond = true;
  e.ReSize(big.get_stateInfo().totalStates); e= 0;


  bool haveEnoughStates = (e.Ncols()< nroots) ? false : true;
#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, haveEnoughStates, 0);
#endif

  if (!haveEnoughStates) {
    //sometimes when you need many roots and at the start of the sweep the hilbert space is not big
    //enough to support all the roots
    
    if (dmrginp.calc_type() != RESPONSE) {
      for (int i=0; i<nroots&&mpigetrank() == 0; i++) {
	solution[i].Randomise();
	Normalise(solution[i]);
      }

    }
    else {

      //****************************
      //GuessWave::guess_wavefunctions(solution[0], e, big, guesswavetype, onedot, currentRoot, 
      //dot_with_sys, 0.0); 
    }
    
  }
  else {
    if(dmrginp.solve_method() == DAVIDSON) {
      //mcheck ("before guess wavefunction");

      //collectedStateInfo
      StateInfo ls=big.get_leftBlock()->get_stateInfo(), rs = big.get_rightBlock()->get_stateInfo();

      if (guesswavetype != BASIC)  {
	GuessWave::guess_wavefunctions(solution, e, big, guesswavetype, onedot, dot_with_sys, nroots, additional_noise, currentRoot); 
      }
      unCollectSolutions(solution, solution.size(), ls, rs);
      unCollectBlocks(big);
      buildDiagonal(big, e);

      if (guesswavetype == BASIC)  
	GuessWave::guess_wavefunctions(solution, e, big, guesswavetype, onedot, dot_with_sys, nroots, additional_noise, currentRoot); 

      multiply_h davidson_f(big, onedot);
      
      solution.resize(dmrginp.deflation_max_size());
      for (int i=nroots; i<solution.size(); i++)
	if (mpigetrank() == 0)
	  solution[i].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), onedot);
      
      for (int istate=0; istate<lowerStates.size(); istate++)  {
	for (int jstate=istate+1; jstate<lowerStates.size(); jstate++) {
	  double overlap = DotProduct(lowerStates[istate], lowerStates[jstate]);
	  ScaleAdd(-overlap/DotProduct(lowerStates[istate], lowerStates[istate]), lowerStates[istate], lowerStates[jstate]);
	}
      }
    

      Linear::block_davidson(solution, e, tol, warmUp, davidson_f, useprecond, currentRoot, lowerStates);

      if (mpigetrank() == 0) {
	for (int i=solution.size()-1; i>=nroots; i--) 
	  solution[i].deallocate();
      }

      //uncollect stateinfo on system and environment
      ls=big.get_leftBlock()->get_stateInfo(); rs = big.get_rightBlock()->get_stateInfo();      
      collectSolutions(solution, nroots, ls, rs);
      collectBlocks(big);

    }
    else if (dmrginp.solve_method() == CONJUGATE_GRADIENT) {

      multiply_h davidson_f(big, onedot);

      if (mpigetrank()!=0) 
	e.ReSize(0);

      //**************************
      //GuessWave::guess_wavefunctions(solution[0], e, big, guesswavetype, onedot, currentRoot, 
      //dot_with_sys, 0.0); 

      if (guesswavetype == BASIC)
	solution[0].Clear();

      double functional = Linear::MinResMethod(solution[0], tol, davidson_f, lowerStates);
      if (mpigetrank() == 0)
	e(1) = functional;

    }
    else {
      pout << "Lanczos is no longer supported"<<endl;
      abort();
      solution.resize(1);
      multiply_h davidson_f(big, onedot);
      //*************************
      //GuessWave::guess_wavefunctions(solution, e, big, guesswavetype, onedot, dot_with_sys, additional_noise, currentRoot); 
    }
  }

  solution.resize(nroots);
  energies.resize(nroots);
  if (haveEnoughStates) {
    for (int i=0; i<nroots&& mpigetrank() == 0;i++) {
      energies[i] = e(i+1);
      //pout << "\t\t\t Energy of wavefunction "<<i<<"  =  "<<e(i+1)<<endl;
    }
  }
  else {
    for (int i=0; i<nroots&& mpigetrank() == 0;i++) {
      if (dmrginp.calc_type() == RESPONSE)
	energies[i] = 1.e10;
      else
	energies[i] = e(1);
    }
  }
#ifndef SERIAL
  broadcast(world, energies, 0);
#endif
  pout<<endl;
}

