/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_OPERATORFUNCTIONS_HEADER_H
#define SPIN_OPERATORFUNCTIONS_HEADER_H


#define TINY 1.e-20

class DiagonalMatrix;
namespace SpinAdapted{
  class StackSpinBlock;
  class StackSparseMatrix;
  class StateInfo;
  class StackMatrix;
  class StackWavefunction;
  class SpinQuantum;

namespace operatorfunctions
{
  //TENSOR TRACE A x I  ->  C
  void TensorTraceElement(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, StackMatrix& cel, int cq, int cqprime, double scale);


//TENSOR TRACE A x I  ->  C  
void TensorTrace(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, double scale= 1.0, int num_thrds = 1) ;

//***********************************************************


//TENSOR TRACE A x I -> cD  
 void TensorTrace (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix* cDiagonal, double scale);


void TensorProduct (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix* cDiagonal, double scale);


//*****************************************************



//TENSOR PRODUCT A x B -> C
 void TensorProductElement(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, StackMatrix& cel, int cq, int cqprime, double scale);

//TENSOR PRODUCT A x B -> C
  void TensorProduct (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, double scale, int num_thrds=1);

//************************************************



 void TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale, int num_thrds=1);

void TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
 void TensorMultiplysplitLeft(const StackSparseMatrix& RightO, const StackSparseMatrix& LeftO, const StackSparseMatrix& DotO, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
 void TensorMultiplysplitRight(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, const StackSparseMatrix& DotO, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
    
    

 void Product (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, StackSparseMatrix& c, double scale);



 void OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const StackSparseMatrix& op1, StackSparseMatrix& op2);

void MultiplyWithOwnTranspose(const StackSparseMatrix& a, StackSparseMatrix& c, double scale);

}
}
#endif
