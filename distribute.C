/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include "communicate.h"
#include "distribute.h"
#include "StackOperators.h"
#include "Stackwavefunction.h"
#ifndef SERIAL
#include "mpi.h"
#endif

namespace SpinAdapted{

void SplitStackmem()
{
  //now we have to distribute remaining memory equally among different threads
  long originalSize = Stackmem[0].size;
  long remainingMem = Stackmem[0].size - Stackmem[0].memused;
  long memPerThrd = remainingMem/numthrds;
  Stackmem[0].size = Stackmem[0].memused+memPerThrd;
  for (int i=1; i<numthrds; i++) {
    Stackmem[i].data = Stackmem[i-1].data+Stackmem[i-1].size;
    Stackmem[i].memused = 0;
    Stackmem[i].size = memPerThrd;
  }
  Stackmem[numthrds-1].size += remainingMem%numthrds; 
}

void MergeStackmem()
{
  //put all the memory again in the zeroth thrd
  for (int i=1; i<numthrds; i++) {
    Stackmem[0].size += Stackmem[i].size;
    Stackmem[i].data = 0;
    Stackmem[i].memused = 0;
    Stackmem[i].size = 0;
  }
}

  
#ifndef SERIAL

#include <boost/mpi/communicator.hpp>

void distributedaccumulate(StackSparseMatrix& component)
{
  dmrginp.datatransfer->start();
  Timer distributetimer;
  boost::mpi::communicator world;
  int size = world.size();
  int rank = world.rank();
  if (size > 1)
  {
    MPI_Allreduce(MPI_IN_PLACE, component.get_data(), component.memoryUsed(), MPI_DOUBLE, MPI_SUM, Calc);
    //MPI::COMM_WORLD.Allreduce(component.get_data(), &tempArray[0], component.memoryUsed(), MPI_DOUBLE, MPI_SUM);
  }
  dmrginp.datatransfer->stop();
}

void distributedaccumulate(DiagonalMatrix& component)
{
  dmrginp.datatransfer->start();
  Timer distributetimer;
  boost::mpi::communicator world;
  int size = world.size();
  int rank = world.rank();
  if (size > 1)
  {
    MPI_Allreduce(MPI_IN_PLACE, component.Store(), component.Ncols(), MPI_DOUBLE, MPI_SUM, Calc);
  }
  dmrginp.datatransfer->stop();
}

#else

void distributedaccumulate(DiagonalMatrix& component) {;}
void distributedaccumulate(SpinAdapted::StackSparseMatrix& component) {;}

#endif

}

