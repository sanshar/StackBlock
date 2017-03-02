/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEP_HEADER
#define SPIN_SWEEP_HEADER

#include <vector>
#include "SpinQuantum.h"

class Matrix;
namespace SpinAdapted{
  class StackSpinBlock;
  class SweepParams;
namespace Sweep
{
  void BlockAndDecimate (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys);
  void Startup (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem);
  double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
  double do_one_partial(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);

  void do_overlap(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
  void fullci(double sweep_tol);
  void tiny(double sweep_tol);

  void CanonicalizeWavefunctionPartialSweep(SweepParams &sweepParams, const bool &forward, int currentstate);
  void CanonicalizeWavefunction(SweepParams &sweepParams, const bool &forward, int currentstate);
  void InitializeStateInfoPartialSweep(SweepParams &sweepParams, const bool &forward, int currentstate);
  void InitializeStateInfo(SweepParams &sweepParams, const bool &forward, int currentstate);
  void InitializeOverlapSpinBlocks(SweepParams &sweepParams, const bool &forward, int stateA, int stateB, int integralIndex);
  void calculateAllOverlap(Matrix& overlap);
  void calculateHMatrixElements(Matrix& H);
  void makeSystemEnvironmentBigBlocks(StackSpinBlock& system, StackSpinBlock& systemDot, StackSpinBlock& newSystem, 
				      StackSpinBlock& environment, StackSpinBlock& environmentDot, StackSpinBlock& newEnvironment,
				      StackSpinBlock& big, SweepParams& sweepParams, const bool& dot_with_sys, const bool& useSlater,
				      int integralIndex, int braState=-1, int ketState=-1, const vector<SpinQuantum>& braquanta= vector<SpinQuantum>(), const vector<SpinQuantum>& ketquanta= vector<SpinQuantum>());
  void makeSystemEnvironmentBigOverlapBlocks(const std::vector<int>& systemSites, StackSpinBlock& systemDot, StackSpinBlock& environmentDot,
					     StackSpinBlock& system, StackSpinBlock& newSystem, StackSpinBlock& environment, StackSpinBlock& newEnvironment,
					     StackSpinBlock& big, SweepParams& sweepParams, const bool& dot_with_sys, const bool& useSlater,
					     int integralIndex, int braState, int ketState);
  
  void set_dot_with_sys(bool& dot_with_sys, const StackSpinBlock& system, const SweepParams& sweepParams, const bool& forward) ;
}
}
#endif

