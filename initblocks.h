/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_INIT_BLOCKS_HEADER
#define SPIN_INIT_BLOCKS_HEADER

#include "enumerator.h"

namespace SpinAdapted{
  class StackSpinBlock;
  class SpinQuantum;

namespace InitBlocks
{
  void InitStartingBlock (StackSpinBlock& startingBlock, const bool &forward, int leftState, int rightState,
                          const int & forward_starting_size, const int &backward_starting_size,
                          const int& restartSize, const bool &restart, const bool& warmUp, int integralIndex, const vector<SpinQuantum>& braquanta = vector<SpinQuantum>(), const vector<SpinQuantum>& ketquanta = vector<SpinQuantum>());

  void InitNewSystemBlock(StackSpinBlock &system, StackSpinBlock &systemDot, StackSpinBlock &newSystem, int leftState, 
			  int rightState, const int &sys_add, const bool &direct, int integralIndex, 
			  const Storagetype &storage, bool haveNormops, bool haveCompops, int constraint=NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

  void InitNewSystemBlock(StackSpinBlock &system, StackSpinBlock &systemDot, StackSpinBlock &newSystem, int leftState, int rightState, const int& sys_add, const bool &direct, int integralIndex, const Storagetype &storage, bool haveNormops, bool haveCompops, int constraint, const std::vector<SpinQuantum>& braquanta, const std::vector<SpinQuantum>& ketquanta);




  void InitBigBlock(StackSpinBlock &leftBlock, StackSpinBlock &rightBlock, StackSpinBlock &big, const vector<SpinQuantum>& braquanta = vector<SpinQuantum>() , const vector<SpinQuantum>& ketquanta = vector<SpinQuantum>() );



  void InitNewEnvironmentBlock(StackSpinBlock &environment, StackSpinBlock& environmentDot, StackSpinBlock &newEnvironment,
                               const StackSpinBlock &system, StackSpinBlock &systemDot, int leftState, int rightState,
                               const int &sys_add, const int &env_add, const bool &forward, const bool &direct, const bool &onedot, 
			       const bool &nexact, const bool &useSlater, int integralIndex, 
			       bool haveNormops, bool haveCompops, const bool& dot_with_sys, int constraint=NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

  void InitNewEnvironmentBlock(StackSpinBlock &environment, StackSpinBlock& environmentDot, StackSpinBlock &newEnvironment, 
						      const StackSpinBlock &system, StackSpinBlock &systemDot, int leftState, int rightState,
						      const int &sys_add, const int &env_add, const bool &forward, const bool &direct, 
						      const bool &onedot, const bool &nexact, const bool &useSlater, int integralIndex, 
						      bool haveNormops, bool haveCompops, const bool& dot_with_sys, int constraint, const std::vector<SpinQuantum>& braquanta, const std::vector<SpinQuantum>& ketquanta);

  void InitNewOverlapEnvironmentBlock(StackSpinBlock &environment, StackSpinBlock& environmentDot, StackSpinBlock &newEnvironment, 
				      const StackSpinBlock &system, StackSpinBlock &systemDot, int leftState, int rightState,
				      const int &sys_add, const int &env_add, const bool &forward, int integralIndex,
				      const bool &onedot, const bool& dot_with_sys, int constraint=NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);


}
}
#endif

