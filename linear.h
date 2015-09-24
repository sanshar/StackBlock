/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_LINEAR_HEADER_H
#define SPIN_LINEAR_HEADER_H
#include <vector>

class DiagonalMatrix;

namespace SpinAdapted{
  class Davidson_functor;
  class StackWavefunction;
  namespace Linear
  {
    void precondition(StackWavefunction& op, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
    void olsenPrecondition(StackWavefunction& op, StackWavefunction& C0, double e, DiagonalMatrix& diagonal, double levelshift=0.0);
    void block_davidson(std::vector<StackWavefunction>& b, DiagonalMatrix& e, double normtol, const bool &warmUp, Davidson_functor& h_mult, bool& useprecond, int currentRoot, std::vector<StackWavefunction>& lowerStates);
    double MinResMethod(StackWavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<StackWavefunction> &lowerStates);
    double ConjugateGradient(StackWavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<StackWavefunction> &lowerStates);
  };
}
#endif

