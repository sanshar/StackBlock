#ifndef SERIAL
#include "communicate.h"
#include "BaseOperator.h"
#include "Wavefunction.h"

template<> void receiveobject<SparseMatrix>(SparseMatrix& object, int from)
{
  boost::mpi::communicator world;
  int tag = 0;
  world.recv(from, tag, object.set_orbs());
  world.recv(from, tag, object.set_deltaQuantum());
  world.recv(from, tag, object.set_quantum_ladder());
  world.recv(from, tag, object.set_build_pattern());
  world.recv(from, tag, object.set_fermion());
  world.recv(from, tag, object.set_initialized());
  world.recv(from, tag, object.set_built());
  world.recv(from, tag, object.set_built_on_disk());
  world.recv(from, tag, object.set_allowedQuantaMatrix());
  world.recv(from, tag, object.set_build_pattern());
  world.recv(from, tag, object.set_sign());

  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	world.recv(from, tag, object.operator_element(i,j));
}

template<> void receiveobject<Wavefunction>(Wavefunction& object, int from)
{
  boost::mpi::communicator world;
  int tag = 0;
  world.recv(from, tag, object.set_orbs());
  world.recv(from, tag, object.set_deltaQuantum());
  world.recv(from, tag, object.set_quantum_ladder());
  world.recv(from, tag, object.set_build_pattern());
  world.recv(from, tag, object.set_fermion());
  world.recv(from, tag, object.set_initialized());
  world.recv(from, tag, object.set_built());
  world.recv(from, tag, object.set_built_on_disk());
  world.recv(from, tag, object.set_allowedQuantaMatrix());
  world.recv(from, tag, object.set_build_pattern());
  world.recv(from, tag, object.set_sign());

  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	world.recv(from, tag, object.operator_element(i,j));
  world.recv(from, tag, object.set_onedot());
}

template<> void sendobject<SparseMatrix>(const SparseMatrix& object, int to)  
{
  boost::mpi::communicator world;
  int tag = 0;

  world.send(to, tag, object.get_orbs());
  world.send(to, tag, object.get_deltaQuantum());
  world.send(to, tag, object.get_quantum_ladder());
  world.send(to, tag, object.get_build_pattern());
  world.send(to, tag, object.get_fermion());
  world.send(to, tag, object.get_initialized());
  world.send(to, tag, object.get_built());
  world.send(to, tag, object.get_built_on_disk());
  world.send(to, tag, object.get_allowedQuantaMatrix());
  world.send(to, tag, object.get_build_pattern());
  world.send(to, tag, object.get_sign());

  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	world.send(to, tag, object.operator_element(i,j));
}

template<> void sendobject<Wavefunction>(const Wavefunction& object, int to)  
{
  boost::mpi::communicator world;
  int tag = 0;
  
  world.send(to, tag, object.get_orbs());
  world.send(to, tag, object.get_deltaQuantum());
  world.send(to, tag, object.get_quantum_ladder());
  world.send(to, tag, object.get_build_pattern());
  world.send(to, tag, object.get_fermion());
  world.send(to, tag, object.get_initialized());
  world.send(to, tag, object.get_built());
  world.send(to, tag, object.get_built_on_disk());
  world.send(to, tag, object.get_allowedQuantaMatrix());
  world.send(to, tag, object.get_build_pattern());
  world.send(to, tag, object.get_sign());
  
  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	world.send(to, tag, object.operator_element(i,j));
  world.send(to, tag, object.get_onedot());
}

#endif
