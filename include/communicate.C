#ifndef SERIAL
#include "communicate.h"
#include "BaseOperator.h"
#include "Wavefunction.h"

template<> void receiveobject<SparseMatrix>(SparseMatrix& object, int from)
{
  boost::mpi::communicator world;
  int tag = 0;
  calc.recv(from, tag, object.set_orbs());
  calc.recv(from, tag, object.set_deltaQuantum());
  calc.recv(from, tag, object.set_quantum_ladder());
  calc.recv(from, tag, object.set_build_pattern());
  calc.recv(from, tag, object.set_fermion());
  calc.recv(from, tag, object.set_initialized());
  calc.recv(from, tag, object.set_built());
  calc.recv(from, tag, object.set_built_on_disk());
  calc.recv(from, tag, object.set_allowedQuantaMatrix());
  calc.recv(from, tag, object.set_build_pattern());
  calc.recv(from, tag, object.set_sign());

  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	calc.recv(from, tag, object.operator_element(i,j));
}

template<> void receiveobject<Wavefunction>(Wavefunction& object, int from)
{
  boost::mpi::communicator world;
  int tag = 0;
  calc.recv(from, tag, object.set_orbs());
  calc.recv(from, tag, object.set_deltaQuantum());
  calc.recv(from, tag, object.set_quantum_ladder());
  calc.recv(from, tag, object.set_build_pattern());
  calc.recv(from, tag, object.set_fermion());
  calc.recv(from, tag, object.set_initialized());
  calc.recv(from, tag, object.set_built());
  calc.recv(from, tag, object.set_built_on_disk());
  calc.recv(from, tag, object.set_allowedQuantaMatrix());
  calc.recv(from, tag, object.set_build_pattern());
  calc.recv(from, tag, object.set_sign());

  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	calc.recv(from, tag, object.operator_element(i,j));
  calc.recv(from, tag, object.set_onedot());
}

template<> void sendobject<SparseMatrix>(const SparseMatrix& object, int to)  
{
  boost::mpi::communicator world;
  int tag = 0;

  calc.send(to, tag, object.get_orbs());
  calc.send(to, tag, object.get_deltaQuantum());
  calc.send(to, tag, object.get_quantum_ladder());
  calc.send(to, tag, object.get_build_pattern());
  calc.send(to, tag, object.get_fermion());
  calc.send(to, tag, object.get_initialized());
  calc.send(to, tag, object.get_built());
  calc.send(to, tag, object.get_built_on_disk());
  calc.send(to, tag, object.get_allowedQuantaMatrix());
  calc.send(to, tag, object.get_build_pattern());
  calc.send(to, tag, object.get_sign());

  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	calc.send(to, tag, object.operator_element(i,j));
}

template<> void sendobject<Wavefunction>(const Wavefunction& object, int to)  
{
  boost::mpi::communicator world;
  int tag = 0;
  
  calc.send(to, tag, object.get_orbs());
  calc.send(to, tag, object.get_deltaQuantum());
  calc.send(to, tag, object.get_quantum_ladder());
  calc.send(to, tag, object.get_build_pattern());
  calc.send(to, tag, object.get_fermion());
  calc.send(to, tag, object.get_initialized());
  calc.send(to, tag, object.get_built());
  calc.send(to, tag, object.get_built_on_disk());
  calc.send(to, tag, object.get_allowedQuantaMatrix());
  calc.send(to, tag, object.get_build_pattern());
  calc.send(to, tag, object.get_sign());
  
  int r = object.get_allowedQuantaMatrix.nrows(), c = object.get_allowedQuantaMatrix.ncols();
  for (int i=0; i<r; i++)
    for (int j=0; j<c; j++)
      if (object.allowed(i,j))
	calc.send(to, tag, object.operator_element(i,j));
  calc.send(to, tag, object.get_onedot());
}

#endif
