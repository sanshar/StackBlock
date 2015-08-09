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
#ifdef SERIAL
int receivefrom(int offsetproc)
{
  return 0;
}
void makesendlist(vector<int>& tolist, int offsetproc) {;}
#else
#include <boost/mpi/communicator.hpp>


int receivefrom(int offsetproc)
{
  boost::mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  rank -= offsetproc;

  
  if (rank < 0) rank += size; // modulo size arithmetic

  
  return ((rank - 1) / Broadcastsettings::npyramid + offsetproc) % size;
}



template<> void distributedaccumulate<SpinAdapted::StackSparseMatrix>(SpinAdapted::StackSparseMatrix& component)
{
  Timer distributetimer;
  boost::mpi::communicator world;
  int size = world.size();
  int rank = world.rank();
  if (size > 1)
  {
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, component.get_data(), component.memoryUsed(), MPI_DOUBLE, MPI_SUM);
    //MPI::COMM_WORLD.Allreduce(component.get_data(), &tempArray[0], component.memoryUsed(), MPI_DOUBLE, MPI_SUM);
  }
}

template<> void distributedaccumulate<DiagonalMatrix>(DiagonalMatrix& component)
{
  Timer distributetimer;
  boost::mpi::communicator world;
  int size = world.size();
  int rank = world.rank();
  if (size > 1)
  {
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, component.Store(), component.Ncols(), MPI_DOUBLE, MPI_SUM);
  }
}

template<> void distributedaccumulate<SpinAdapted::StackHam>(SpinAdapted::StackHam& component)
{
  Timer distributetimer;
  boost::mpi::communicator world;
  int size = world.size();
  int rank = world.rank();
  if (size > 1)
  {
    
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, component.get_data(), component.memoryUsed(), MPI_DOUBLE, MPI_SUM);
  }

}

template<> void distributedaccumulate<SpinAdapted::StackWavefunction>(SpinAdapted::StackWavefunction& component)
{
  Timer distributetimer;
  boost::mpi::communicator world;
  int size = world.size();
  int rank = world.rank();
  if (size > 1)
  {
    
    MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, component.get_data(), component.memoryUsed(), MPI_DOUBLE, MPI_SUM);
  }

}



void makesendlist(vector<int>& tolist, int offsetproc)
{
  boost::mpi::communicator world;
  int rank = world.rank();
  int size = world.size();
  rank -= offsetproc; // if offsetproc == rank, then treat current proc as the root node
  if (rank < 0) rank += size;
  for (int i = 1; i <= Broadcastsettings::npyramid; ++i)
    {
      int tosend = rank * Broadcastsettings::npyramid +i;      
      if (tosend < size)
	tolist.push_back((tosend + offsetproc) % size); // if offsetproc is not zero, add it back on, and then wraparound
    }
}


#endif

}

