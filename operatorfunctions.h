/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_OPERATORFUNCTIONS_HEADER_H
#define SPIN_OPERATORFUNCTIONS_HEADER_H
#include "StateInfo.h"
#include "StackBaseOperator.h"
#include "timer.h"
#include "Stackspinblock.h"
#include "MatrixBLAS.h"
#include <math.h>
#include "global.h"
#include "StackMatrix.h"
#include <omp.h>
#include <iostream>
#include <map>
#include <vector>


#define TINY 1.e-20

namespace SpinAdapted{
namespace operatorfunctions
{
  //TENSOR TRACE A x I  ->  C
  void TensorTraceElement(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, StackMatrix& cel, int cq, int cqprime, double scale);


//TENSOR TRACE A x I  ->  C  
void TensorTrace(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, double scale= 1.0, int num_thrds = 1) ;

//***********************************************************


//TENSOR TRACE A x I -> cD  
 void TensorTrace (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix& cDiagonal, Real scale);


void TensorProduct (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix& cDiagonal, double scale);


//*****************************************************



//TENSOR PRODUCT A x B -> C
 void TensorProductElement(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, StackMatrix& cel, int cq, int cqprime, double scale);

//TENSOR PRODUCT A x B -> C
  void TensorProduct (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, double scale, int num_thrds=1);

//************************************************



 void TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale, int num_thrds=1);

void TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
    
    

 void Product (const StackSpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, double scale);



 void OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const Baseoperator<Matrix>& op1, Baseoperator<Matrix>& op2);

void MultiplyWithOwnTranspose(const StackSparseMatrix& a, StackSparseMatrix& c, Real scale);

}
}
#endif
