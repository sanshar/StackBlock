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

namespace SpinAdapted
{
  class SparseMatrix;
  class Ham;
  class DensityMatrix;
  class Wavefunction;
  class StackSparseMatrix;
  class StackWavefunction;
  class StackHam;

  // tree mpi::broadcast algorithm ...
  void makesendlist(std::vector<int>& tolist, int offsetproc = 0);

  int receivefrom(int offsetproc = 0);

  namespace Broadcastsettings
  {
    const int npyramid = 2;
  };
#ifdef SERIAL
  template<class T> void distributedaccumulate(T& component){;}
#else
  //void broadcast(SpinAdapted::Ham& object, int proc);
  //void broadcast(SpinAdapted::DensityMatrix& object, int proc);
  //void broadcast(SpinAdapted::Wavefunction& object, int proc);

  template<class T> void plus_equals(T& lhs, const T& rhs)
    {
      lhs += rhs;
    }

  template<class T> void distributedaccumulate(T& component)
    {
      Timer distributetimer;
      boost::mpi::communicator world;
      int size = world.size();
      int rank = world.rank();
      if (size > 1)
	{
	  // i.  look at our sendlist
	  // ii. receive stuff from all the nodes in our sendlist
	  // iii. accumulate
	  // iv. send to receivefromnode
	  // v. mpi::broadcast
      
	  std::vector<int> sendlist;
	  makesendlist(sendlist);
      
	  // receive and accumulate
	  int listsize = sendlist.size();
	  for (int i = 0; i < listsize; ++i)
	    {
	      T part;
	      receiveobject(part, sendlist[i]);
	      plus_equals(component, part);
	    }
      
	  // send to root (if not already root)
	  if (rank != 0) {
	    sendobject(component, receivefrom());
	  }

	  boost::mpi::broadcast(world,component,0);
	}
      //     pout << "distribution time " << distributetimer.elapsedwalltime() << " " << distributetimer.elapsedcputime() << endl;
    }
  template<> void distributedaccumulate<DiagonalMatrix>(DiagonalMatrix& component);
  template<> void distributedaccumulate<SpinAdapted::StackSparseMatrix>(SpinAdapted::StackSparseMatrix& component);
  template<> void distributedaccumulate<SpinAdapted::StackHam>(SpinAdapted::StackHam& component);
  template<> void distributedaccumulate<SpinAdapted::StackWavefunction>(SpinAdapted::StackWavefunction& component);
#endif
#ifdef SERIAL
  template<class T> void accumulateto(T& component, int proc){;}
#else
  template<class T> void accumulateto(T& component, int proc)
    {
      Timer distributetimer;
      boost::mpi::communicator world;
      int size = world.size();
      int rank = world.rank();
      if (size > 1)
	{
	  // i.  look at our sendlist
	  // ii. receive stuff from all the nodes in our sendlist
	  // iii. accumulate
	  // iv. send to receivefromnode
	  // v. mpi::broadcast
      
	  std::vector<int> sendlist;
	  makesendlist(sendlist);
      
	  // receive and accumulate
	  int listsize = sendlist.size();
	  for (int i = 0; i < listsize; ++i)
	    {
	      T part;
	      receiveobject(part, sendlist[i]);
	      component += part;	  
	    }
      
	  // send to root (if not already root)
	  if (rank != 0)
	    sendobject(component, receivefrom());

	  // send back out to proc from root
	  if (rank == 0)
	    sendobject(component, proc);
	  else if (rank == proc)
	    receiveobject(component, 0);
	}
      //  pout << "distribution time " << distributetimer.elapsedwalltime() << " " << distributetimer.elapsedcputime() << endl;
    }



#endif


  template<class T> void initiateMultiThread( T* op, T* &op_array, int MAX_THRD)
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
  
  
  
  
  template<class T> void accumulateMultiThread(T* op, T* &op_array, int MAX_THRD)
    {
#ifndef SERIAL
      boost::mpi::communicator world;
      int size = world.size();
      
      if (size == 1 && MAX_THRD == 1)
	return;
      else if (size > 1 && MAX_THRD == 1) { //only mpi
	distributedaccumulate(*op);
      }
      else if (size == 1 && MAX_THRD > 1) {  //only multithreaded
	for (int i=MAX_THRD-1; i>0; i--) {
	  ScaleAdd(1.0, op_array[i], op_array[0]);
	  op_array[i].deallocate();
	}
	delete [] op_array;
      }
      else {  //multithreaded and mpi
	for (int i=MAX_THRD-1; i>0; i--) {
	  ScaleAdd(1.0, op_array[i], op_array[0]);
	  op_array[i].deallocate();
	}
	delete [] op_array;
	
	distributedaccumulate(*op);
	
      }
#else
      if ( MAX_THRD == 1)
	return;
      else {  //only multithreaded
	for (int i=MAX_THRD-1; i>0; i--) {
	  ScaleAdd(1.0, op_array[i], op_array[0]);
	  op_array[i].deallocate();
	}
	delete [] op_array;
      }
#endif
      
    }

  //only accumulate mpi
  template<class T> void accumulateMultiProc(T* op)
    {
#ifndef SERIAL
      boost::mpi::communicator world;
      int size = world.size();
      
      if (size == 1)
	return;
      else if (size > 1) { //only mpi
	distributedaccumulate(*op);
      }
#else
      return;
#endif
    }      

}

#endif
