/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SWEEPGENBLOCK_HEADER
#define SWEEPGENBLOCK_HEADER


namespace SpinAdapted{
  class SweepParams;
  class StackSpinBlock;

namespace SweepGenblock
{
  void BlockAndDecimate (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int stateA, int stateB);
  double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int stateA, int stateB);
  void do_one(SweepParams &sweepParams, const bool &forward, int stateA, int stateB);
  void do_one_partialSweep(SweepParams &sweepParams, const bool &forward, int stateA, int stateB, int integralIndex);

};
}

#endif
