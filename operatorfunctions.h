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
 //***************************************


 //*****************WHEN LOOP BLOCK IS SPLIT*************
 double getScaling(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& leftOp, const StackSparseMatrix& ldotOp, 
		   const StackSparseMatrix& RIGHTOP, const StackSparseMatrix& rightOp, const StackSparseMatrix& rdotOp, 
		   const SpinQuantum& luncollectedQ, const SpinQuantum& lQ, const SpinQuantum& ldotQ,
		   const SpinQuantum& luncollectedQPrime, const SpinQuantum& lQPrime, const SpinQuantum& ldotQPrime,
		   const SpinQuantum& runcollectedQ, const SpinQuantum& rQ, const SpinQuantum& rdotQ,
		   const SpinQuantum& runcollectedQPrime, const SpinQuantum& rQPrime, const SpinQuantum& rdotQPrime);

 void TensorMultiplysplitLeftsplitRight(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, 
					const StackSparseMatrix& leftOp, const StackSparseMatrix& ldotOp, 
					const StackSparseMatrix& rightOp, const StackSparseMatrix& rdotOp, 
					const StackSpinBlock *cblock, StackWavefunction& ropCmat, StackWavefunction* v, 
					double scale);
 void TensorMultiplysplitLeftsplitRight00(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, 
					  const StackSparseMatrix& leftOp, const StackSparseMatrix& ldotOp, 
					  const StackSparseMatrix& rightOp, const StackSparseMatrix& rdotOp, 
					  const StackSpinBlock *cblock, StackWavefunction& ropCmat, StackWavefunction* v, 
					  double scale);

 void TensorMultiplyleftdot(const StackSparseMatrix& leftOp, StateInfo *cstate, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
 void TensorMultiplydotop(const StackSparseMatrix& dotOp, StateInfo *cstate, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
 void TensorMultiplysplitLeft(const StackSparseMatrix& RightO, const StackSparseMatrix& LeftO, const StackSparseMatrix& DotO, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
 void TensorMultiplysplitLeftElement(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int index, double scale);
 void TensorMultiplyCDxCDsplitLeftElement(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int rQPrime, double scale);
 void TensorMultiplyCDxCDsplitLeftElementcopy(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int rQPrime, double scale);
 void TensorMultiplysplitRight(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, const StackSparseMatrix& DotO, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
 void TensorMultiplysplitRightElement(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, const StackSparseMatrix& DotO, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int lQPrime, double scale);
 void TensorMultiplyCDxCDsplitRightElement(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, const StackSparseMatrix& DotO, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int lQPrime, double scale, bool doTranspose);
 void TensorMultiplyCDxCDsplitRightElementcopy(const StackSparseMatrix& LeftO, const StackSparseMatrix& RightO, const StackSparseMatrix& DotO, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int lQPrime, double scale, bool doTranspose);
 //****************************************************

 //***************the order of contraction is different, but the loopblock is still split******************** 
 void multiplyDotRightElement(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& rightOp,
			      const StackMatrix& cMat, const StackMatrix& rightOpmat, StackMatrix& v,  
			      const SpinQuantum& luncollectedQ, const SpinQuantum& lQ, const SpinQuantum& dotQ, const SpinQuantum& rightQ,
			      const SpinQuantum& luncollectedQPrime, const SpinQuantum& lQPrime, const SpinQuantum& dotQPrime, const SpinQuantum& rightQPrime,		      
			      double scale);
 
 void multiplyDotRight(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotop, 
		       StackSparseMatrix& rightop, std::vector<StackMatrix>& lopCmat, 
		       StackWavefunction* v,  const StackSpinBlock* cblock, int luncollectedQPrime, int rQPrime, double scale);


 void multiplyDotLeftElement(const StackSparseMatrix& RIGHTOP, const StackSparseMatrix& rightOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& leftOp,
			      const StackMatrix& cMat, const StackMatrix& leftOpmat, StackMatrix& v,  
			      const SpinQuantum& lQ, const SpinQuantum& runcollectedQ, const SpinQuantum& dotQ, const SpinQuantum& leftQ,
			      const SpinQuantum& lQPrime, const SpinQuantum& runcollectedQPrime, const SpinQuantum& dotQPrime, const SpinQuantum& leftQPrime,		      
			      double scale); 
 void multiplyDotLeft(const StackSparseMatrix& RIGHTOP, const StackSparseMatrix& rightop, const StackSparseMatrix& dotop, 
		      StackSparseMatrix& leftop, std::vector<StackMatrix>& ropCmat, 
		      StackWavefunction* v,  const StackSpinBlock* cblock, int lQPrime, int luncollectedQPrime, double scale);
 //*********************************************** 

 void TensorMultiplyLeftLeft(const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale);
 void TensorMultiplyRight(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);

 void TensorMultiplyRightLeft(const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale);
 void TensorMultiplyLeft(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale);
    
    

 void Product (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, StackSparseMatrix& c, double scale);



 void OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const StackSparseMatrix& op1, StackSparseMatrix& op2);

 void OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const StackSparseMatrix& op1, StackSparseMatrix& op2, const std::vector<int>& rows, const std::vector<int>& cols);


void MultiplyWithOwnTranspose(const StackSparseMatrix& a, StackSparseMatrix& c, double scale);

 void braTensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix &a, const StackSpinBlock *cblock, StackWavefunction &c, StackWavefunction &v, double scale, int num_thrds=1);

  void TensorMultiply(const StackSparseMatrix& a, const StackSparseMatrix& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, bool aIsLeft, double scale);

  void TensorMultiply(const StackSparseMatrix a, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, bool left, double scale);
}
}
#endif
