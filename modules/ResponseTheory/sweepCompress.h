/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEPCOMPRESS_HEADER
#define SPIN_SWEEPCOMPRESS_HEADER

namespace SpinAdapted{
  class SweepParams;
  class StackSpinBlock;
namespace SweepCompress
{
  void BlockDecimateAndCompress (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, int baseState);
  double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int targetState, int baseState);
  void Startup (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool& dot_with_sys, int targetState, int baseState);
  void WavefunctionCanonicalize (SweepParams &sweepParams, StackSpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int correctionVector, int baseState);
};
}
#endif

