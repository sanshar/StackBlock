/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_SWEEPRESPONSE_HEADER
#define SPIN_SWEEPRESPONSE_HEADER

#include <vector>

namespace SpinAdapted{
  class SweepParams;
  class StackSpinBlock;
namespace SweepResponse
{
  void BlockAndDecimate (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, int targetState, std::vector<int>& projectors, std::vector<int>& baseState);
  double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int targetState, std::vector<int>& projectors, std::vector<int>& baseState, int correctionVector=1);
  void StartUp (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool& dot_with_sys, int targetState, int correctionVector, std::vector<int>& baseState, std::vector<int>& firstorderstate);
  void WavefunctionCanonicalize (SweepParams &sweepParams, StackSpinBlock& system, const bool &useSlater, const bool& dot_with_sys, int targetState, std::vector<int>& correctionVector, std::vector<int>& baseState);


};


};
#endif

