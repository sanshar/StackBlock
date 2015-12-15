/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "davidson.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Stackwavefunction.h"
#include "Stackspinblock.h"
#include "StackBaseOperator.h"
#include "MatrixBLAS.h"


SpinAdapted::multiply_h::multiply_h(const StackSpinBlock& b, const bool &onedot_) : block(b){}


void SpinAdapted::multiply_h::operator()(StackWavefunction& c, StackWavefunction& v)
{
  block.multiplyH( c, &v, MAX_THRD);
  //block.multiplyH_2index( c, &v, MAX_THRD);
  //pout << "c" << endl;
  //pout << c << endl;
  //pout << "v" << endl;
  //pout << v << endl;
}

SpinAdapted::multiply_h_2Index::multiply_h_2Index(const StackSpinBlock& b, const bool &onedot_) : block(b){}


void SpinAdapted::multiply_h_2Index::operator()(StackWavefunction& c, StackWavefunction& v)
{
  block.multiplyH_2index( c, &v, MAX_THRD);
  //block.multiplyH3index( c, &v, MAX_THRD);
}





