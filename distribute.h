/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_DISTRIBUTE_HEADER_H
#define SPIN_DISTRIBUTE_HEADER_H

#include <iostream>
#include <communicate.h>
#include "timer.h"
#include "pario.h"
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include <vector>

class DiagonalMatrix;
namespace SpinAdapted
{
  class StackSparseMatrix;
  
  void SplitStackmem();
  void MergeStackmem();
  
  template<class T> void initiateMultiThread(T* op, T* &op_array, int MAX_THRD)
  {
    if (MAX_THRD == 1) {
      op_array = op;
    }
    else if (MAX_THRD > 1) {
      op_array = new T[MAX_THRD];
      op_array[0] = *op;
      for (int i=1; i<MAX_THRD; i++)
	op_array[i].deepClearCopy(*op);
    }
  }
  
  
  template <class T> void accumulateMultiThread(T* op, T* &op_array, int MAX_THRD)
  {
    if ( MAX_THRD == 1)
      return;
    else {  //only multithreaded
      for (int i=MAX_THRD-1; i>0; i--) {
	ScaleAdd(1.0, op_array[i], op_array[0]);
	op_array[i].deallocate();
      }
      delete [] op_array;
    }
  }

#ifndef SERIAL
  void distributedaccumulate(DiagonalMatrix& component);
  void distributedaccumulate(SpinAdapted::StackSparseMatrix& component);
#else
  void distributedaccumulate(DiagonalMatrix& component) ;
  void distributedaccumulate(SpinAdapted::StackSparseMatrix& component);
#endif

}

#endif
