/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SOLVER_HEADER_H
#define SPIN_SOLVER_HEADER_H
#include <vector>
#include "enumerator.h"

namespace SpinAdapted{
  class StackSpinBlock;
  class StackWavefunction;

  namespace Solver
  {
    void solve_wavefunction(std::vector<StackWavefunction>& solution, std::vector<double>& energies, StackSpinBlock& big, const double tol, 
			    const guessWaveTypes& guesswavetype, const bool &onedot, const bool& dot_with_sys, const bool& warmUp, const bool& twoindex, 
			    double additional_noise, int currentRoot, std::vector<StackWavefunction>& lowerStates);
  };
}
#endif
