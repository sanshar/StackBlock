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
#include "MatrixBLAS.h"

SpinAdapted::multiply_h::multiply_h(const StackSpinBlock& b, const bool &onedot_) : block(b){
  if (dmrginp.prebuild() && dmrginp.doimplicitTranspose())
    twoidxcomps = b.prebuild(MAX_THRD);
}

void SpinAdapted::multiply_h::operator()(StackWavefunction& c, StackWavefunction& v)
{
  block.multiplyH( c, &v, MAX_THRD);
}

SpinAdapted::multiply_h::~multiply_h() {
  if (dmrginp.prebuild() && dmrginp.doimplicitTranspose()) {
    for (int i = twoidxcomps.size()-1; i >= 0; --i)
      twoidxcomps.at(i)->deallocate();
  }
}

SpinAdapted::multiply_h_2Index::multiply_h_2Index(const StackSpinBlock& b, const bool &onedot_) : block(b){
  if (dmrginp.prebuild() && dmrginp.doimplicitTranspose())
    twoidxcomps = b.prebuild(MAX_THRD);
}


void SpinAdapted::multiply_h_2Index::operator()(StackWavefunction& c, StackWavefunction& v)
{
  block.multiplyH_2index( c, &v, MAX_THRD);
}


SpinAdapted::multiply_h_2Index::~multiply_h_2Index() {
  if (dmrginp.prebuild() && dmrginp.doimplicitTranspose()) {
    for (int i = twoidxcomps.size()-1; i >= 0; --i)
      twoidxcomps.at(i)->deallocate();
  }
}
