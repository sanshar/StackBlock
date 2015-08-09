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


/**
 * Distributed loops to be used functors on OperatorArrays
 * 
 */

/**
 * Loop over all local (i.e. stored on current processor) in array
 * 
 */


namespace SpinAdapted{
class SpinBlock;
class StackSpinBlock;


template<class A> void singlethread_build(A& array, SpinBlock& b, std::vector< Csf >& s, vector< vector<Csf> >& ladders)
{
  for (int i = 0; i < array.get_size(); ++i) {
    //typedef typename A::OpType Op;
    std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
    for (int j=0; j<vec.size(); j++) {
      vec[j]->buildUsingCsf(b, ladders, s);
      //pout << *vec[j]<<endl;
    }
  }
}

template<class A> void singlethread_build(A& array, SpinBlock& b)
{
  for (int i = 0; i < array.get_size(); ++i) {
    //typedef typename A::OpType Op;
    //std::vector<boost::shared_ptr<Op> >& vec = array.get_local_element(i);
    std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
    for (int j=0; j<vec.size(); j++) {
      vec[j]->build(b);
      //pout << *vec[j]<<endl;
    }
  }
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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



//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Used with functions designed to build operators in core, but here we actually write to disk instead
template<typename T2, class A> void for_all_operators_to_disk(A& array, const SpinBlock& b, std::ofstream& ofs, const T2& func)
{     
  for (int i = 0; i < array.get_size(); ++i) {
    std::vector<boost::shared_ptr<SparseMatrix> > vec = array.get_local_element(i);
    for (int j=0; j<vec.size(); j++){
      // Don't build if already built!
      assert( ! vec.at(j)->get_built() );
      assert( ! vec.at(j)->get_built_on_disk() );
        
      // Apply function to operator
      func( *(vec.at(j)) );
      
      // Store on disk
      vec.at(j)->set_built_on_disk() = true;
      boost::archive::binary_oarchive save_op(ofs);
      save_op << *(vec.at(j));
     
      // Deallocate memory for operator representation
      vec.at(j)->set_built() = false;
//FIXME is Clear() sufficient to deallocate?
      vec.at(j)->Clear();
    }
  }
}   

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

}
#endif


