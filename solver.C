/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


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


void memorySummary(StackSpinBlock& big, vector<StackWavefunction>& solution) {

  long sysmem, envmem, sysop=0, envop=0;
  StackSpinBlock* sys = big.get_leftBlock()->get_leftBlock() == 0 ? big.get_leftBlock() : big.get_leftBlock()->get_leftBlock();
  StackSpinBlock* env = big.get_rightBlock()->get_leftBlock() == 0 ? big.get_rightBlock() : big.get_rightBlock()->get_leftBlock();
  int opidx = dmrginp.spinAdapted() ? 1 : 0;
  if (!big.get_leftBlock()->has(CRE_DES)) {
    if (big.get_leftBlock()->get_op_array(CRE_DESCOMP).get_size()!=0)
      sysop = big.get_leftBlock()->get_op_array(CRE_DESCOMP).get_local_element(0)[opidx]->memoryUsed();
    if (big.get_rightBlock()->get_op_array(CRE_DES).get_size()!=0)
      envop = big.get_rightBlock()->get_op_array(CRE_DES).get_local_element(0)[opidx]->memoryUsed();
  }
  else {
    if (big.get_leftBlock()->get_op_array(CRE_DES).get_size()!=0)
      sysop = big.get_leftBlock()->get_op_array(CRE_DES).get_local_element(0)[opidx]->memoryUsed();
    if (big.get_rightBlock()->has(CRE_DESCOMP)) {
      if (big.get_rightBlock()->get_op_array(CRE_DESCOMP).get_size()!=0)
	envop = big.get_rightBlock()->get_op_array(CRE_DESCOMP).get_local_element(0)[opidx]->memoryUsed();
    }
  }
  p2out << str(boost::format("%-40s - %-10.4f\n") % "Total memory" % (Stackmem[0].size*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Memory used" % (Stackmem[0].memused*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->System" % ((sys->memoryUsed()+sys->additionalMemoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->Envrionment" % ((env->memoryUsed()+env->additionalMemoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Memory left" % ((Stackmem[0].size-Stackmem[0].memused)*8/1.e9));
  p2out << str(boost::format("%-40s - %3i x %-10.4f = %-10.4f\n") % "      |-->wavefunction" %(3*solution.size()+2*numthrds) %((solution[0].memoryUsed())*8/1.e9) %(((3*solution.size()+2*numthrds)*solution[0].memoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %3i x %-10.4f = %-10.4f\n") % "      |-->sys op" %(numthrds) %((sysop)*8/1.e9) %(numthrds*sysop*8/1.e9));
  p2out << str(boost::format("%-40s - %3i x %-10.4f = %-10.4f\n") % "      |-->env op" %(numthrds) %((envop)*8/1.e9) %(numthrds*envop*8/1.e9));

}

void SpinAdapted::Solver::solve_wavefunction(vector<StackWavefunction>& solution, vector<double>& energies, StackSpinBlock& big, const double tol, 
					     const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp, const bool& twoindex,
					     double additional_noise, int currentRoot, std::vector<StackWavefunction>& lowerStates)
{
  for (int thrd=0; thrd<numthrds; thrd++) 
    dmrginp.matmultFlops[thrd] = 0.0;

  const int nroots = dmrginp.setStateSpecific() ? 1 : dmrginp.nroots();
  dmrginp.makediagonal->start();
  DiagonalMatrix e;
  bool useprecond = true;

  e.ReSize(big.get_stateInfo().totalStates); e= 0;
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

  dmrginp.makediagonal->stop();

  bool haveEnoughStates = (e.Ncols()< nroots) ? false : true;
#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, haveEnoughStates, 0);
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
      dmrginp.guesswf->start();
      //mcheck ("before guess wavefunction");
      solution.resize(dmrginp.deflation_max_size());
      if (mpigetrank()==0) {
        memorySummary(big, solution);
        for (int i=nroots; i<solution.size(); i++) {
	  solution[i].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), onedot);
        }
      }

      //multiply_h davidson_f(big, onedot);
      Davidson_functor* davidson_f;
      if (twoindex) {
        davidson_f = new multiply_h_2Index(big, onedot);
      } else {
        davidson_f = new multiply_h(big, onedot);
      }

      guessWaveTypes guesstype = guesswavetype;
      if (guesswavetype == TRANSPOSE && big.get_leftBlock()->get_rightBlock() == 0)
	guesstype = BASIC;
      GuessWave::guess_wavefunctions(solution, e, big, guesstype, onedot, dot_with_sys, nroots, additional_noise, currentRoot); 
      dmrginp.guesswf->stop();
      
      for (int istate=0; istate<lowerStates.size(); istate++)  {
	for (int jstate=istate+1; jstate<lowerStates.size(); jstate++) {
	  double overlap = DotProduct(lowerStates[istate], lowerStates[jstate]);
	  ScaleAdd(-overlap/DotProduct(lowerStates[istate], lowerStates[istate]), lowerStates[istate], lowerStates[jstate]);
	}
      }
    
      dmrginp.blockdavid->start();
      Linear::block_davidson(solution, e, tol, warmUp, *davidson_f, useprecond, currentRoot, lowerStates);
      dmrginp.blockdavid->stop();

      if (mpigetrank() == 0) {
	for (int i=solution.size()-1; i>=nroots; i--) 
	  solution[i].deallocate();
      }
      delete davidson_f;
    }
    else if (dmrginp.solve_method() == CONJUGATE_GRADIENT) {

      //multiply_h davidson_f(big, onedot);
      Davidson_functor* davidson_f;
      if (twoindex) {
        davidson_f = new multiply_h_2Index(big, onedot);
      } else {
        davidson_f = new multiply_h(big, onedot);
      }

      if (mpigetrank()!=0) 
	e.ReSize(0);

      //**************************
      GuessWave::guess_wavefunctions(solution[0], e, big, guesswavetype, onedot, currentRoot, 
				     dot_with_sys, 0.0); 

      if (guesswavetype == BASIC)
	solution[0].Clear();

      double functional = Linear::MinResMethod(solution[0], tol, *davidson_f, lowerStates);
      //double functional = Linear::ConjugateGradient(solution[0], tol, davidson_f, lowerStates);
      if (mpigetrank() == 0)
	e(1) = functional;
      delete davidson_f;
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
  broadcast(calc, energies, 0);
#endif
  pout<<endl;
}

