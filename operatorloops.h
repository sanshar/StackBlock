/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_OPERATOR_LOOPS_HEADER_H
#define SPIN_OPERATOR_LOOPS_HEADER_H
#include <vector>
#include <iostream>
#include <communicate.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <boost/shared_ptr.hpp>
#include "pario.h"
#include "StackBaseOperator.h"

/**
 * Distributed loops to be used functors on OperatorArrays
 * 
 */

/**
 * Loop over all local (i.e. stored on current processor) in array
 * 
 */


namespace SpinAdapted{
class StackSpinBlock;



//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//execute a function on all elements of an array
//single thread and multithread versions of the code
template<typename T2, class A> void for_all_singlethread_hmult(A& array, const T2& func)
{
  int i;
  {
    for (i = 0; i < array.get_size(); ++i) {
      for (int j=0; j<array.get_local_element(i).size(); j++)
	func(array.get_local_element(i)[j]);
    }
  }
}

//execute a function on all elements of an array
//single thread and multithread versions of the code
template<typename T2, class A> void for_all_singlethread(A& array, const T2& func)
{
  int i;
  {
    for (i = 0; i < array.get_size(); ++i) {
      func(array.get_local_element(i));
    }
  }
}


template<typename T2, class A> void for_all_operators_singlethread(A& array, const T2& func)
{
  int i;
    for (i = 0; i < array.get_size(); ++i) {
      int vecsize = array.get_local_element(i).size();
      for (int j=0; j<vecsize; j++){
	func(*(array.get_local_element(i)[j]));
      }
    }

}




template<class A> void singlethread_build(A& array, StackSpinBlock& b, std::vector< Csf >& s, vector< vector<Csf> >& ladders)
{
  for (int i = 0; i < array.get_size(); ++i) {
    //typedef typename A::OpType Op;
    int vecsize = array.get_local_element(i).size();
    for (int j=0; j<vecsize; j++) {
      array.get_local_element(i)[j]->buildUsingCsf(b, ladders, s);
    }
  }
}

template<class A> void singlethread_build(A& array, StackSpinBlock& b)
{
  for (int i = 0; i < array.get_size(); ++i) {
    //typedef typename A::OpType Op;
    //std::vector<boost::shared_ptr<Op> >& vec = array.get_local_element(i);
    int vecsize = array.get_local_element(i).size();
    for (int j=0; j<vecsize; j++) {
      array.get_local_element(i)[j]->build(b);
      //pout << *vec[j]<<endl;
    }
  }
}




}
#endif


