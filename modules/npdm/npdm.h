/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_HEADER
#define NPDM_HEADER

namespace SpinAdapted{
  class SweepParams;
  enum NpdmOrder{NPDM_NEVPT2, NPDM_ONEPDM, NPDM_TWOPDM, NPDM_THREEPDM, NPDM_FOURPDM, NPDM_PAIRMATRIX, NPDM_OVERLAP, NPDM_EMPTY};
namespace Npdm{
  class Npdm_driver_base;

  void npdm(NpdmOrder npdm_order, bool transitionpdm=false);
  double npdm_do_one_sweep(Npdm_driver_base &npdm_driver, SweepParams &sweepParams, const bool &warmUp, const bool &forward, 
			   const bool &restart, const int &restartSize, const int state, const int stateB);

}
}

#endif
