/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef NPDM_EXPECTATIONS_ENGINE_H
#define NPDM_EXPECTATIONS_ENGINE_H

#include "Stackspinblock.h"
#include "Stackwavefunction.h"
#include "StackBaseOperator.h"

namespace SpinAdapted{
namespace Npdm{

  void FormLeftOp(const StackSpinBlock* leftBlock, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, StackSparseMatrix& Aop, int totalspin);
  double DotProduct(const StackWavefunction& w1, const StackWavefunction& w2, const StackSpinBlock& big);
  double spinExpectation(StackWavefunction& wave1, StackWavefunction& wave2, StackSparseMatrix &leftOp, StackSparseMatrix& dotOp, StackSparseMatrix& rightOp, const StackSpinBlock& big);
}
}

#endif

