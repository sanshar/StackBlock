/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_GLOBAL_HEADER
#define SPIN_GLOBAL_HEADER
 
#define WANT_MATH
#define WANT_STREAM
// Matrix library header files


#include "timer.h"
//#include <malloc.h>
#include <newmat.h>
#include <newmatap.h>
#include <newmatio.h>
//#include <IntegralMatrix.h>
// STL headers
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
#include <numeric>
// C headers
//#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Symmetry.h"
#include "input.h"
#include "pario.h"
#include <ctime>
#include <iostream>
// Global variables
#include <boost/interprocess/managed_shared_memory.hpp>
#include "alloc.h"

using namespace std;

namespace SpinAdapted{
class OneElectronArray;
class TwoElectronArray;
class PairArray;
class CCCCArray;
class CCCDArray;
class StackSpinBlock;
class PerturbTwoElectronArray;
enum OnePerturbType;
enum TwoPerturbType;

extern boost::interprocess::shared_memory_object segment;
extern boost::interprocess::mapped_region region;

extern int CACHEBUFFER;
extern int MAX_THRD;

extern Timer globaltimer;

extern std::vector<OneElectronArray> v_1;
extern std::vector<TwoElectronArray> v_2;
extern std::map<TwoPerturbType,PerturbTwoElectronArray> vpt_2;
extern OneElectronArray vpt_1;
extern std::vector<PairArray> v_cc;
extern std::vector<CCCCArray> v_cccc;
extern std::vector<CCCDArray> v_cccd;
extern std::vector<double> coreEnergy;
extern std::vector< std::vector<StackSpinBlock> > singleSiteBlocks;
extern double BWPTenergy;

extern Input dmrginp;

extern bool SHOW_MORE;
extern bool DEBUG_MEMORY;
extern bool RESTART;
extern bool FULLRESTART;
extern bool BACKWARD;
extern bool reset_iter;
extern bool restartwarm;
extern string sym;
extern bool NonabelianSym;
extern std::vector<int> NPROP;
extern int PROPBITLEN;
extern double NUMERICAL_ZERO;
extern std::vector<StackAllocator<double> > Stackmem; 
#ifndef SERIAL
extern boost::mpi::communicator calc;
extern MPI_Comm Calc;
#endif
}
#endif
