/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef COMMUNICATE_HEADER_H
#define COMMUNICATE_HEADER_H
#ifndef SERIAL
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
namespace SpinAdapted{
class Ham;
class Wavefunction;
extern boost::mpi::communicator calc;
}

inline int mpigetrank() { return SpinAdapted::calc.rank(); }
inline int mpigetsize() { return SpinAdapted::calc.size(); }
#else
inline int mpigetrank() { return 0; }
inline int mpigetsize() { return 1; }
#endif

#ifdef SERIAL 
template<class T> void sendobject(const T& object, int to) {}
#else
template<class T> void sendobject(const T& object, int to)  
{
  int tag = 0;
  SpinAdapted::calc.send(to,tag,object);
}
#endif
// default argument + template specialization bug workaround
//template<class T> void receiveobject(T& object, int from, int tag = 0)
#ifdef SERIAL
template<class T> void receiveobject(T& object, int from) {}
#else
template<class T> void receiveobject(T& object, int from)
{
  int tag = 0;
  SpinAdapted::calc.recv(from,tag,object);
}

//These are implemented because sometimes the Operator size is very big (maybe more than 2Gb) and this causes 
//mpi implementation to crash. Atleast for now, the parameters in the mpi function calls are 32 big integers.
//which means that there is overflow when the size of buffer being passed in mpi calls is bigger than about 2Gb.

template<> void receiveobject<SpinAdapted::Ham>(SpinAdapted::Ham& object, int from);

template<> void receiveobject<SpinAdapted::Wavefunction>(SpinAdapted::Wavefunction& object, int from);

template<> void sendobject<SpinAdapted::Ham>(const SpinAdapted::Ham& object, int to);  

template<> void sendobject<SpinAdapted::Wavefunction>(const SpinAdapted::Wavefunction& object, int to);  

#endif

#endif
