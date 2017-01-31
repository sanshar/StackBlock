/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include <sched.h>
#include <unistd.h>
#include "IntegralMatrix.h"
#include <fstream>
#include "input.h"
#include "pario.h"
#include "global.h"
#include "orbstring.h"
#include "least_squares.h"
#include <include/communicate.h>
#include "sweepgenblock.h"
#include "stdlib.h"
#include "npdm.h"

#ifdef _OPENMP
#include <omp.h>
#ifdef _HAS_INTEL_MKL 
#include "mkl.h"
#include "mkl_cblas.h"
#endif
#endif

//the following can be removed later
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/format.hpp>
#include "Stackspinblock.h"
#include "StateInfo.h"
#include "operatorfunctions.h"
#include "solver.h"
#include "davidson.h"
#include "rotationmat.h"
#include "sweep.h"
#include "sweepCompress.h"
#include "sweepResponse.h"
#include "dmrg_wrapper.h"
#include "sweeponepdm.h"
#include "screen.h"
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "pario.h"
#include "Stackwavefunction.h"
#include "StateInfo.h"
#include "Stackdensity.h"
#include "initblocks.h"
#include <boost/filesystem.hpp>

#ifdef USE_BTAS
void calculateOverlap();
#endif
void dmrg(double sweep_tol);
void partialsweepDMRG(double sweep_tol);
void compress(double sweep_tol, int targetState, int baseState);
void responseSweep(double sweep_tol, int targetState, vector<int>& correctionVector, vector<int>& baseState);
void responsepartialSweep(double sweep_tol, int targetState, vector<int>& correctionVector, vector<int>& baseState);
void restart(double sweep_tol, bool reset_iter);
void dmrg_stateSpecific(double sweep_tol, int targetState);
void ReadInput(char* conf);
void fullrestartGenblock();
void Npdm(int pdm, bool restartpdm, bool transitionpdm);

void license() {
#ifndef MOLPRO
  pout << "Block  Copyright (C) 2012  Garnet K.-L. Chan"<<endl;
  pout << "This program comes with ABSOLUTELY NO WARRANTY; for details see license file."<<endl;
  pout << "This is free software, and you are welcome to redistribute it"<<endl;
  pout << "under certain conditions; see license file for details."<<endl;
#endif
}


namespace SpinAdapted{
  Timer globaltimer(false);
  bool DEBUGWAIT = false;
  bool DEBUG_MEMORY = false;
  bool restartwarm = false;
  double NUMERICAL_ZERO = 1e-15;
  double BWPTenergy = 0.0;
  std::vector<OneElectronArray> v_1;
  std::vector<TwoElectronArray> v_2;
  std::vector<PairArray> v_cc;
  std::vector<CCCCArray> v_cccc;
  std::vector<CCCDArray> v_cccd;
  std::vector<std::vector<StackSpinBlock> > singleSiteBlocks;
  std::vector<double> coreEnergy;
  std::vector<std::vector<std::vector<int> > > screened_cre;
  std::vector<std::vector<std::vector<int> > > screened_credesdes;
  std::vector< std::vector<std::vector<pair<int,int> > > > screened_credes;
  std::vector< std::vector<std::vector<pair<int,int> > > > screened_credescomp;
  std::vector< std::vector<std::vector<pair<int,int> > > > screened_desdes;
  std::vector< std::vector<std::vector<pair<int,int> > > > screened_desdescomp;

  Input dmrginp;
  int CACHEBUFFER = 0;
  int MAX_THRD = 1;
  bool FULLRESTART;
  bool RESTART;
  bool BACKWARD;
  bool reset_iter;
  std::vector<int> NPROP;
  int PROPBITLEN=1;

  boost::interprocess::shared_memory_object segment(boost::interprocess::open_or_create, ("Integrals" + to_string(time(NULL) % 1000000)).c_str(), boost::interprocess::read_write);
  boost::interprocess::mapped_region region;

  std::vector<StackAllocator<double> > Stackmem;
#ifndef SERIAL
  boost::mpi::communicator calc;
  MPI_Comm Calc;
#endif

}

using namespace SpinAdapted;

#ifndef SERIAL
void sleepBarrier(boost::mpi::communicator world, int tag, double psleep)
{
  int size = world.size(), rank = world.rank();

  if (size == 1)
    return;

  int mask = 1;
  while (mask < size) {
    int dst = (rank + mask) % size;
    int src = (rank - mask + size) % size;
    string s = "h";
    boost::mpi::request r = world.isend(dst, tag, s);
    //req = comm.isend(None, dst, tag);
    while (! world.iprobe(src, tag)) 
      sleep(psleep);

    world.recv(src, tag, s);
    r.wait();
    mask *= 2;
  }

}
#endif

int calldmrg(char* input, char* output)
{
  srand(1000);
  streambuf *backup;
  backup = cout.rdbuf();
  ofstream file;
  if (output != 0) {
    file.open(output);
    pout.rdbuf(file.rdbuf());
  }
  license();


#ifndef SERIAL
  boost::mpi::communicator world;
#endif

  ReadInput(input);
  dmrginp.matmultFlops.resize(numthrds, 0.);

  int size = 1, rank =0;

#ifndef SERIAL
  size = world.size(); rank = world.rank();
  MPI_Group orig_group, calc_group;
  MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
  std::vector<int>& m_calc_procs = dmrginp.calc_procs();


  //Make the new mpi comm which only has the calc procs 
  if (find(m_calc_procs.begin(), m_calc_procs.end(), rank)!=m_calc_procs.end()) {
    string& lprefix = dmrginp.load_prefix();
    string& sprefix = dmrginp.save_prefix();
    string oldrank = str(boost::format("%s%i")  %"/node" % rank);
    string newrank = str(boost::format("%s%i")  %"/node" % (find(m_calc_procs.begin(), m_calc_procs.end(), rank) - m_calc_procs.begin()));
    int index=0; index = lprefix.find(oldrank, index);
    lprefix.replace(index, oldrank.length(), newrank);
    int index2=0; index2 = sprefix.find(oldrank, index2);
    sprefix.replace(index2, oldrank.length(), newrank);
    MPI_Comm_split(MPI_COMM_WORLD, 0, find(m_calc_procs.begin(), m_calc_procs.end(), rank) - m_calc_procs.begin(), &Calc);
  }
  else
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &Calc);
  calc = boost::mpi::communicator(Calc, boost::mpi::comm_attach);


  cpu_set_t  *oldmask, newmask;
  int oldset, newset;
  //if only a subset of processors are being used, then unset processor affinity
  if (m_calc_procs.size() < world.size()) {
    sched_getaffinity(0, oldset, oldmask);
    CPU_ZERO( &newmask);
    int maxcpus = 2500;
    for (int i=0; i<2500; i++)
      CPU_SET(i, &newmask);
    int success = sched_setaffinity(0, sizeof(&newmask), &newmask);
    if (success == -1) {
      cout << world.rank()<<endl;
    }
  }
#endif

#ifndef SERIAL
  if (find(m_calc_procs.begin(), m_calc_procs.end(), rank)!=m_calc_procs.end()) {
#endif

  MAX_THRD = dmrginp.thrds_per_node()[mpigetrank()];
  int mkl_thrd = dmrginp.mkl_thrds();

#ifndef SERIAL
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  for (int i=0; i<calc.size(); i++) {
    if (mpigetrank() == i) {
      cout << "processor: "<<rank<<"  numthrds: "<<MAX_THRD<<"  node name: "<<processor_name<<endl;
    }
    calc.barrier();
  }
#endif

#ifdef _OPENMP
  omp_set_num_threads(MAX_THRD);
#ifdef _HAS_INTEL_MKL 
  mkl_set_num_threads(mkl_thrd);
  mkl_set_dynamic(0);
#endif
  omp_set_nested(1);
#endif
  cout.precision (12);
  cout << std::fixed;

  cout <<"allocating " <<dmrginp.getMemory()<<" doubles "<<endl;
  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();
 //************
  //memset(stackmemory, 0, dmrginp.getMemory()*sizeof(double));
  dmrginp.initCumulTimer();



  pout << "**** STACK MEMORY REMAINING ***** "<<1.0*(Stackmem[0].size-Stackmem[0].memused)*sizeof(double)/1.e9<<" GB"<<endl;
  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  SweepParams sweepParams;

  SweepParams sweep_copy;
  bool direction_copy; int restartsize_copy;
  Matrix O, H;
  
  //switch(dmrginp.calc_type()) {
    
  if (dmrginp.calc_type() == COMPRESS)
  {
    bool direction; int restartsize;
    //sweepParams.restorestate(direction, restartsize);
    //sweepParams.set_sweep_iter() = 0;
    restartsize = 0;

    int targetState, baseState, correctionVector, firstorderstate;
    {
      direction = true;

      //base state is always defined
      baseState = dmrginp.baseStates()[0];

      //if targetstate is given use it otherwise use basestate+1
      if(dmrginp.targetState() == -1)
	targetState = dmrginp.baseStates()[0]+1;
      else
	targetState = dmrginp.targetState();

      algorithmTypes atype = dmrginp.algorithm_method();
      dmrginp.set_algorithm_method() = ONEDOT;
      //initialize state info and canonicalize wavefunction is always done using onedot algorithm
      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, baseState);
	Sweep::InitializeStateInfo(sweepParams, !direction, baseState);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, baseState);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, baseState);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, baseState);
      }
      dmrginp.set_algorithm_method() = atype;
    }

    //this genblock is required to generate all the nontranspose operators
    dmrginp.setimplicitTranspose() = false;
    SweepGenblock::do_one(sweepParams, false, false, false, restartsize, baseState, baseState);


    compress(sweep_tol, targetState, baseState);
  }
  else if (dmrginp.calc_type() == RESPONSEBW)
  {
    //compressing the V|\Psi_0>, here \Psi_0 is the basestate and 
    //its product with V will have a larger bond dimension and is being compressed
    //it is called the target state

    //dmrginp.setimplicitTranspose() = false;


    sweepParams.restorestate(direction, restartsize);
    algorithmTypes atype = dmrginp.algorithm_method();
    dmrginp.set_algorithm_method() = ONEDOT;
    if (mpigetrank()==0 && !RESTART && !FULLRESTART && dmrginp.get_sweep_type() != FULL) {
      for (int l=0; l<dmrginp.projectorStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
      }
      for (int l=0; l<dmrginp.baseStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
      }
    }
    dmrginp.set_algorithm_method() = atype;

    
    pout << "DONE COMPRESSING THE CORRECTION VECTOR"<<endl;
    pout << "NOW WE WILL OPTIMIZE THE RESPONSE WAVEFUNCTION"<<endl;
    //finally now calculate the response state
    if (dmrginp.get_sweep_type() == FULL)
      responseSweep(sweep_tol, dmrginp.targetState(), dmrginp.projectorStates(), dmrginp.baseStates());
    else
      responsepartialSweep(sweep_tol, dmrginp.targetState(), dmrginp.projectorStates(), dmrginp.baseStates());

  }
  else if (dmrginp.calc_type() == RESPONSE || dmrginp.calc_type() == RESPONSELCC || dmrginp.calc_type() == RESPONSEAAAV || dmrginp.calc_type() == RESPONSEAAAC)
  {
    //compressing the V|\Psi_0>, here \Psi_0 is the basestate and 
    //its product with V will have a larger bond dimension and is being compressed
    //it is called the target state
    //dmrginp.setimplicitTranspose() = false;


    sweepParams.restorestate(direction, restartsize);
    algorithmTypes atype = dmrginp.algorithm_method();
    dmrginp.set_algorithm_method() = ONEDOT;
    if (mpigetrank()==0 && !RESTART && !FULLRESTART && dmrginp.get_sweep_type() == FULL) {
      for (int l=0; l<dmrginp.projectorStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.projectorStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.projectorStates()[l]);
      }
      for (int l=0; l<dmrginp.baseStates().size(); l++) {
	Sweep::InitializeStateInfo(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::InitializeStateInfo(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, dmrginp.baseStates()[l]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, dmrginp.baseStates()[l]);
      }
    }
    dmrginp.set_algorithm_method() = atype;

    
    pout << "DONE COMPRESSING THE CORRECTION VECTOR"<<endl;
    pout << "NOW WE WILL OPTIMIZE THE RESPONSE WAVEFUNCTION"<<endl;
    //finally now calculate the response state
    if (dmrginp.get_sweep_type() == FULL)
      responseSweep(sweep_tol, dmrginp.targetState(), dmrginp.projectorStates(), dmrginp.baseStates());
    else
      responsepartialSweep(sweep_tol, dmrginp.targetState(), dmrginp.projectorStates(), dmrginp.baseStates());


  }
  else if (dmrginp.calc_type() == CALCOVERLAP)
  {
    pout.precision(12);
    if (mpigetrank() == 0) {
      for (int istate = 0; istate<dmrginp.nroots(); istate++) {
	bool direction;
	int restartsize;
	sweepParams.restorestate(direction, restartsize);
	Sweep::InitializeStateInfo(sweepParams, !direction, istate);
	Sweep::InitializeStateInfo(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, istate);
      }
      for (int istate = 0; istate<dmrginp.nroots(); istate++) 
	for (int j=istate; j<dmrginp.nroots() ; j++) {
	  int integralIndex = 0;
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, j, istate, integralIndex);
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, j, istate, integralIndex);
	}
      //Sweep::calculateAllOverlap(O);
    }
  }
  else if (dmrginp.calc_type() == CALCHAMILTONIAN)
  {
    pout.precision(12);

    for (int istate = 0; istate<dmrginp.nroots(); istate++) {
      bool direction;
      int restartsize;
      sweepParams.restorestate(direction, restartsize);
      
      if (mpigetrank() == 0) {
	Sweep::InitializeStateInfo(sweepParams, !direction, istate);
	Sweep::InitializeStateInfo(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, istate);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, istate);
      }
    }
    
    //Sweep::calculateHMatrixElements(H);
    pout << "overlap "<<endl<<O<<endl;
    pout << "hamiltonian "<<endl<<H<<endl;
  }
  else if (dmrginp.calc_type() == DMRG ||
	   dmrginp.calc_type() == ONEPDM ||
	   dmrginp.calc_type() == TWOPDM ||
	   dmrginp.calc_type() == THREEPDM ||
	   dmrginp.calc_type() == TRANSITION_ONEPDM ||
	   dmrginp.calc_type() == TRANSITION_TWOPDM ||
	   dmrginp.calc_type() == TRANSITION_THREEPDM)
  {
    if (dmrginp.get_sweep_type() != FULL)
	partialsweepDMRG(sweep_tol);
    else {
      if (RESTART && !FULLRESTART)
	restart(sweep_tol, reset_iter);
      else if (FULLRESTART) {
	fullrestartGenblock();
	reset_iter = true;
	sweepParams.restorestate(direction, restartsize);

	if (!direction) {
	  double last_fe = Sweep::do_one(sweepParams, false, direction, true, restartsize);
	}
	sweepParams.calc_niter();
	sweepParams.savestate(direction, restartsize);
	restart(sweep_tol, reset_iter);
      }
      else if (BACKWARD) {
	fullrestartGenblock();
	reset_iter = true;
	sweepParams.restorestate(direction, restartsize);
	sweepParams.calc_niter();
	sweepParams.savestate(direction, restartsize);
	restart(sweep_tol, reset_iter);
      }
      else 
	dmrg(sweep_tol);
    }
    if (dmrginp.calc_type() == ONEPDM) 
      Npdm::npdm(NPDM_ONEPDM);
    //Npdm(1, false, false);
    if (dmrginp.calc_type() == TRANSITION_ONEPDM) 
      Npdm::npdm(NPDM_ONEPDM, true);
    if (dmrginp.calc_type() == TWOPDM) 
      Npdm::npdm(NPDM_TWOPDM);
    if (dmrginp.calc_type() == TRANSITION_TWOPDM) 
      Npdm::npdm(NPDM_TWOPDM, true);
    if (dmrginp.calc_type() == THREEPDM) 
      Npdm::npdm(NPDM_THREEPDM);
    if (dmrginp.calc_type() == TRANSITION_THREEPDM) 
      Npdm::npdm(NPDM_THREEPDM,true);
  }
  else if (dmrginp.calc_type() ==FCI) {
    Sweep::fullci(sweep_tol);
  }    
  else if (dmrginp.calc_type() == TINYCALC) {
    Sweep::tiny(sweep_tol);
  }
  else if (dmrginp.calc_type() == RESTART_ONEPDM) {
    Npdm::npdm(NPDM_ONEPDM);
  }
  else if (dmrginp.calc_type() == RESTART_T_ONEPDM) {
    Npdm::npdm(NPDM_ONEPDM,true);
  }
  else if (dmrginp.calc_type() == RESTART_TWOPDM) {
    Npdm::npdm(NPDM_TWOPDM);
  }
  else if (dmrginp.calc_type() == RESTART_T_TWOPDM) {
    Npdm::npdm(NPDM_TWOPDM,true);
  }
  else if (dmrginp.calc_type() == RESTART_THREEPDM) {
    Npdm::npdm(NPDM_THREEPDM);
  }
  else if (dmrginp.calc_type() == RESTART_T_THREEPDM) {
    Npdm::npdm(NPDM_THREEPDM,true);
  }
  else {
    pout << "Invalid calculation types" << endl; abort();
  }
    /*
  case (TWOPDM):
    Npdm::npdm(NPDM_TWOPDM);
    break;

  case (THREEPDM):
    Npdm::npdm(NPDM_THREEPDM);
    break;

  case (FOURPDM):
    Npdm::npdm(NPDM_FOURPDM);
    break;

  case (NEVPT2PDM):
    Npdm::npdm(NPDM_NEVPT2);
    break;

  case(NEVPT2):
    nevpt2::nevpt2();
    break;
    
  case (RESTART_ONEPDM):
    Npdm::npdm(NPDM_ONEPDM,true);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true);
    }
    break;

  case (RESTART_TWOPDM):
    Npdm::npdm(NPDM_TWOPDM,true);
    break;
  case (RESTART_THREEPDM):
    Npdm::npdm(NPDM_THREEPDM,true);
    break;
  case (RESTART_FOURPDM):
    Npdm::npdm(NPDM_FOURPDM,true);
    break;
  case (RESTART_NEVPT2PDM):
    Npdm::npdm(NPDM_NEVPT2,true);
    break;
  case (TRANSITION_ONEPDM):
    Npdm::npdm(NPDM_ONEPDM,false,true);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true,true);      
    }
    break;
  case (TRANSITION_TWOPDM):
    Npdm::npdm(NPDM_TWOPDM,false,true);
    break;
  case (RESTART_T_ONEPDM):
    Npdm::npdm(NPDM_ONEPDM,true,true);
    if (dmrginp.hamiltonian() == BCS) {
      Npdm::npdm(NPDM_PAIRMATRIX,true,true);      
    }
    break;
  case (RESTART_T_TWOPDM):
    Npdm::npdm(NPDM_TWOPDM,true,true);
    break;
  case(RESTART_NEVPT2):
    nevpt2::nevpt2_restart();
    break;
    */

  cout.rdbuf(backup);
  double cputime = globaltimer.totalcputime();
  double walltime = globaltimer.totalwalltime();
  pout << setprecision(3) <<"\n\n\t\t\t BLOCK CPU  Time (seconds): " << cputime << endl;
  pout << setprecision(3) <<"\t\t\t BLOCK Wall Time (seconds): " << walltime << endl;

  delete [] stackmemory;
#ifndef SERIAL
  }

  //world.barrier();
  sleepBarrier(world, 0, 10);
  MPI_Comm_free(&Calc);
  sched_setaffinity(0, sizeof(oldmask), oldmask);
#endif

  return 0;
}


void calldmrg_(char* input, char* output) {
   int a;
   //a=calldmrg("dmrg.inp",0);//, output);
   a=calldmrg(input,0);//, output);
}


void fullrestartGenblock() {
  SweepParams sweepParams, sweepParamsTmp;
  bool direction; int restartsize;
//Temporary fix to restore sweep direction
//FIXME: NN wrote: please let me know if this makes some erroneous behaviors
  sweepParamsTmp.restorestate(direction, restartsize);

  sweepParams.set_sweep_iter() = 0;
  sweepParams.current_root() = -1;
//direction = true;
  restartsize = 0;

  SweepGenblock::do_one(sweepParams, false, !direction, RESTART, restartsize, -1, -1);
  
  sweepParams.restorestate(direction, restartsize);
  sweepParams.set_sweep_iter()=0;
  sweepParams.set_block_iter() = 0;
  
  sweepParams.savestate(direction, restartsize);
}  


void restart(double sweep_tol, bool reset_iter)
{
  double last_fe = 100.;
  double last_be = 100.;
  double old_fe = 0.;
  double old_be = 0.;
  bool direction;
  int restartsize;
  SweepParams sweepParams;
  bool dodiis = false;

  int domoreIter = 2;

  sweepParams.restorestate(direction, restartsize);

  if (!dmrginp.setStateSpecific()) {
    if(reset_iter) { //this is when you restart from the start of the sweep
      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
    }
    
    if (restartwarm)
      last_fe = Sweep::do_one(sweepParams, true, direction, true, restartsize);
    else
      last_fe = Sweep::do_one(sweepParams, false, direction, true, restartsize);
    
    
    while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	   (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
      {
	
	old_fe = last_fe;
	old_be = last_be;
	if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	  break;
	last_be = Sweep::do_one(sweepParams, false, !direction, false, 0);
	
	
	if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	  break;
	last_fe = Sweep::do_one(sweepParams, false, direction, false, 0);	
      }
  }
  else { //this is state specific calculation  
    const int nroots = dmrginp.nroots();

    bool direction;
    int restartsize;
    sweepParams.restorestate(direction, restartsize);

    //initialize state and canonicalize all wavefunctions
    int currentRoot = sweepParams.current_root();
    for (int i=0; i<nroots; i++) {
      sweepParams.current_root() = i;
      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
      }
    }

    //now generate overlaps with all the previous wavefunctions
    for (int i=0; i<currentRoot; i++) {
      sweepParams.current_root() = i;
      if (mpigetrank()==0) {
	for (int j=0; j<i; j++) {
	  int integralIndex = 0;
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	}
      }
    }
    sweepParams.current_root() = currentRoot;

    if (sweepParams.current_root() <0) {
      p1out << "This is most likely not a restart calculation and should be done without the restart command!!"<<endl;
      p1out << "Aborting!!"<<endl;
      exit(0);
    }
    pout << "RESTARTING STATE SPECIFIC CALCULATION OF STATE "<<sweepParams.current_root()<<" AT SWEEP ITERATION  "<<sweepParams.get_sweep_iter()<<endl;

    //this is so that the iteration is not one ahead after generate block for restart
    --sweepParams.set_sweep_iter(); sweepParams.savestate(direction, restartsize);
    for (int i=sweepParams.current_root(); i<nroots; i++) {
      sweepParams.current_root() = i;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	for (int j=0; j<i ; j++) {
	  int integralIndex = 0;
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	  Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
	}
      }
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, i, i);
      
      
      p1out << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<endl;
      dmrg_stateSpecific(sweep_tol, i);
      p1out << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<" FINSIHED"<<endl;

      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
      sweepParams.savestate(!direction, restartsize);
    }

    p1out << "ALL STATE SPECIFIC CALCUALTIONS FINISHED"<<endl;
  }


  if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()){
    pout << "\n\t\t\t Maximum sweep iterations achieved " << std::endl;
  }

}

void dmrg(double sweep_tol)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;

  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  // warm up sweep ...
  bool dodiis = false;

  int domoreIter = 0;
  bool direction;

  //this is regular dmrg calculation
  if(!dmrginp.setStateSpecific()) {
    sweepParams.current_root() = -1;
    last_fe = Sweep::do_one(sweepParams, true, true, false, 0);
    direction = false;
    while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	   (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = Sweep::do_one(sweepParams, false, false, false, 0);
      direction = true;
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      //For obtaining the extrapolated energy
      old_states=sweepParams.get_keep_states();
      new_states=sweepParams.get_keep_states_ls();
      
      last_fe = Sweep::do_one(sweepParams, false, true, false, 0);
      direction = false;
      
      new_states=sweepParams.get_keep_states();
      
      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      if (domoreIter == 2) {
	dodiis = true;
	break;
      }
      
    }
  }
  else { //this is state specific calculation  
    const int nroots = dmrginp.nroots();

    bool direction=true;
    int restartsize;
    //sweepParams.restorestate(direction, restartsize);
    //sweepParams.set_sweep_iter() = 0;
    //sweepParams.set_restart_iter() = 0;

    algorithmTypes atype;
    pout << "STARTING STATE SPECIFIC CALCULATION "<<endl;
    for (int i=0; i<nroots; i++) {
      atype = dmrginp.algorithm_method();
      dmrginp.set_algorithm_method() = ONEDOT;
      sweepParams.current_root() = i;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);

      }

      for (int j=0; j<i ; j++) {
	int integralIndex = 0;
	Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
      }
      dmrginp.set_algorithm_method() = atype;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, i, i);
      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
      sweepParams.savestate(!direction, restartsize);

      
      pout << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<endl;
      dmrg_stateSpecific(sweep_tol, i);
      pout << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<" FINSIHED"<<endl;
    }

    pout << "ALL STATE SPECIFIC CALCUALTIONS FINISHED"<<endl;
  }
}


void partialsweepDMRG(double sweep_tol)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;

  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  // warm up sweep ...
  bool dodiis = false;

  int domoreIter = 0;
  bool direction=true;

  int restartsize = 0;
  if (RESTART)
    sweepParams.restorestate(direction, restartsize);
  else if (FULLRESTART)
    fullrestartGenblock();

  //this is regular dmrg calculation
  if(!dmrginp.setStateSpecific()) {
    if (!( (RESTART) && sweepParams.get_sweep_iter() >= 2)) {
      sweepParams.current_root() = -1;
      dmrginp.get_sweep_type() = FULL;
      if (FULLRESTART)
	last_fe = Sweep::do_one(sweepParams, false, direction, false, 0);
      else
	last_fe = Sweep::do_one(sweepParams, true, direction, false, 0);
      last_be = Sweep::do_one(sweepParams, false, !direction, false, 0);
      dmrginp.get_sweep_type() = PARTIAL;
    }
    while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol) || 
	   (dmrginp.algorithm_method() == TWODOT_TO_ONEDOT && dmrginp.twodot_to_onedot_iter()+1 >= sweepParams.get_sweep_iter()) )
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = Sweep::do_one_partial(sweepParams, false, direction, false, 0);
      direction = !direction;
      //last_be = Sweep::do_one_partial(sweepParams, false, true, false, 0);
      //direction = true;
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      //For obtaining the extrapolated energy
      old_states=sweepParams.get_keep_states();
      new_states=sweepParams.get_keep_states_ls();
      
      last_fe = Sweep::do_one_partial(sweepParams, false, direction, false, 0);
      direction = !direction;
      //last_fe = Sweep::do_one_partial(sweepParams, false, false, false, 0);
      //direction = false;
      
      new_states=sweepParams.get_keep_states();
      
      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      if (domoreIter == 2) {
	dodiis = true;
	break;
      }
      
    }

    if (!direction)
      Sweep::do_one_partial(sweepParams, false, false, false, 0);
    //do two final full sweeps to store the rotation matrices
    dmrginp.get_sweep_type() = FULL;
    sweepParams.set_backward_starting_size() = 1;
    last_fe = Sweep::do_one(sweepParams, false, true, false, 0);
    last_be = Sweep::do_one(sweepParams, false, false, false, 0);
    dmrginp.get_sweep_type() = PARTIAL;

  }
  else { //this is state specific calculation  
    const int nroots = dmrginp.nroots();

    bool direction=true;
    int restartsize;
    //sweepParams.restorestate(direction, restartsize);
    //sweepParams.set_sweep_iter() = 0;
    //sweepParams.set_restart_iter() = 0;

    algorithmTypes atype;
    pout << "STARTING STATE SPECIFIC CALCULATION "<<endl;
    for (int i=0; i<nroots; i++) {
      atype = dmrginp.algorithm_method();
      dmrginp.set_algorithm_method() = ONEDOT;
      sweepParams.current_root() = i;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      if (mpigetrank()==0) {
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, i);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, direction, i);
	Sweep::InitializeStateInfo(sweepParams, !direction, i);

      }

      for (int j=0; j<i ; j++) {
	int integralIndex = 0;
	Sweep::InitializeOverlapSpinBlocks(sweepParams, direction, i, j, integralIndex);
	Sweep::InitializeOverlapSpinBlocks(sweepParams, !direction, i, j, integralIndex);
      }
      dmrginp.set_algorithm_method() = atype;

      p1out << "RUNNING GENERATE BLOCKS FOR STATE "<<i<<endl;

      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, i, i);
      sweepParams.set_sweep_iter() = 0;
      sweepParams.set_restart_iter() = 0;
      sweepParams.savestate(!direction, restartsize);

      
      pout << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<endl;
      dmrg_stateSpecific(sweep_tol, i);
      pout << "STATE SPECIFIC CALCULATION FOR STATE: "<<i<<" FINSIHED"<<endl;
    }

    pout << "ALL STATE SPECIFIC CALCUALTIONS FINISHED"<<endl;
  }
}


void responseSweep(double sweep_tol, int targetState, vector<int>& projectors, vector<int>& baseStates)
{
  double last_fe = 1.e6;
  double last_be = 1.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;

  bool direction, warmUp, restart;
  int restartSize=0;
  direction = true; //forward
  warmUp = true; //startup sweep
  restart = false; //not a restart

  sweepParams.current_root() = -1;

  algorithmTypes atype = dmrginp.algorithm_method();
  dmrginp.set_algorithm_method() = ONEDOT;

  //the baseState is the initial guess for the targetState
  if (FULLRESTART) {
    sweepParams.restorestate(direction, restartSize);
    direction = !direction;
    dmrginp.setGuessState() = targetState;
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
    bool tempdirection;
    sweepParams.restorestate(tempdirection, restartSize);
    sweepParams.calc_niter();
    sweepParams.set_sweep_iter() = 0;
    sweepParams.set_restart_iter() = 0;
    sweepParams.savestate(tempdirection, restartSize);
  }
  else if (RESTART) {
    dmrginp.set_algorithm_method() = atype;
    warmUp = false;
    restart = true;
    sweepParams.restorestate(direction, restartSize);
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
  }
  else 
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);

  dmrginp.set_algorithm_method() = atype;
  restart = false;
  restartSize = 0;
  warmUp = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;

      last_be = SweepResponse::do_one(sweepParams, warmUp, !direction, restart, restartSize, targetState, projectors, baseStates);
      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);

      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }
  
}

void responsepartialSweep(double sweep_tol, int targetState, vector<int>& projectors, vector<int>& baseStates)
{
  double last_fe = 1.e6;
  double last_be = 1.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;

  bool direction, warmUp, restart;
  int restartSize=0;
  direction = true; //forward
  warmUp = true; //startup sweep
  restart = false; //not a restart

  sweepParams.current_root() = -1;

  algorithmTypes atype = dmrginp.algorithm_method();
  dmrginp.set_algorithm_method() = ONEDOT;

  //the baseState is the initial guess for the targetState
  if (FULLRESTART) {
    sweepParams.restorestate(direction, restartSize);
    direction = !direction;
    dmrginp.setGuessState() = targetState;
    
    last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
    bool tempdirection;
    sweepParams.restorestate(tempdirection, restartSize);
    sweepParams.calc_niter();
    sweepParams.set_sweep_iter() = 0;
    sweepParams.set_restart_iter() = 0;
    sweepParams.savestate(tempdirection, restartSize);
  }
  else if (RESTART) {
    dmrginp.set_algorithm_method() = atype;
    warmUp = false;
    restart = true;
    sweepParams.restorestate(direction, restartSize);
    if (sweepParams.get_sweep_iter() < 1)
      last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
    else
      direction = !direction;
  }
  else {

    SpinQuantum oldq;
    //take the basestate wavefunction and expand it to accomodate large end block
    if (mpigetrank() == 0) {

      StackWavefunction w; StateInfo state;

      if (StackWavefunction::exists(-1)) {
	cout << "copying state -1 to "<<baseStates[0]<<endl;
	StackWavefunction::CopyState(-1, baseStates[0]); 
      }
      else
	StackWavefunction::CopyState(baseStates[0], -1); 
    }

    if (dmrginp.calc_type() == RESPONSEAAAV) {
      StackWavefunction w; StateInfo state;
      int start = dmrginp.getPartialSweep();
      int end   = dmrginp.last_site();
      StackSpinBlock site(start-1, start-1, 1, false);      	
      StackSpinBlock system = StackSpinBlock::buildBigEdgeBlock(start, end, false, true, 1, false);
      
      system.addAdditionalOps();
      StackSpinBlock newSystem;
      newSystem.default_op_components(false, false, true, false);
      newSystem.set_integralIndex() = 1;
      newSystem.setstoragetype(DISTRIBUTED_STORAGE);
      newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, site);
      
      
      vector<int> sites(2); sites[0] = 0; sites[1] = start-2;
      if (mpigetrank() == 0) {	
	w.LoadWavefunctionInfo(state, sites, baseStates[0], true);
	
	StackWavefunction w2; 
	w2.initialise(dmrginp.effective_molecule_quantum_vec(), *state.leftStateInfo, newSystem.get_stateInfo(), true);
	w2.UnCollectQuantaAlongColumns(*state.leftStateInfo, newSystem.get_stateInfo());
	for (int i=0; i<w.nrows(); i++)
	  for (int j=0; j<w.ncols(); j++) 
	    if (w.allowed(i,j)) {
	      int J = newSystem.get_stateInfo().unCollectedStateInfo->quantaMap(0,j)[0];
	      copy (w(i,j), w2(i,J));
	    }
	w2.CollectQuantaAlongColumns(*state.leftStateInfo, const_cast<StateInfo&>(*newSystem.get_stateInfo().unCollectedStateInfo));
	StateInfo bigstate;
	TensorProduct(*state.leftStateInfo, const_cast<StateInfo&>(newSystem.get_stateInfo()), bigstate,  PARTICLE_SPIN_NUMBER_CONSTRAINT);      
	w2.SaveWavefunctionInfo(bigstate, sites, baseStates[0]);
	w2.deallocate();
	w.deallocate();
      }

      //make edge blocks for all the combinations
      sites[0] = start-1; sites[1] = end-1;
      StackSpinBlock::store(false, sites, newSystem, targetState, baseStates[0]);  
      //system.deallocate(); system.clear();
	
      newSystem.set_integralIndex() = 0;
      //system = StackSpinBlock::buildBigEdgeBlock(start, end, false, true, 0, false);
      StackSpinBlock::store(false, sites, newSystem, targetState, baseStates[0]);  
      newSystem.deallocate(); newSystem.clear();
      system.removeAdditionalOps();
      system.deallocate();
      site.deallocate(); 
      
      {	
	StackSpinBlock site(start-1, start-1, 0, false);      	
	system = StackSpinBlock::buildBigEdgeBlock(start, end, false, true, 0, true);
	system.addAdditionalOps();
	StackSpinBlock newSystem;
	newSystem.default_op_components(false, false, true, true);
	newSystem.set_integralIndex() = 0;
	newSystem.setstoragetype(DISTRIBUTED_STORAGE);
	newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, site);
	StackSpinBlock::store(false, sites, newSystem, targetState, targetState);  
	StackSpinBlock::store(false, sites, newSystem, baseStates[0], baseStates[0]);  
	newSystem.deallocate();
	newSystem.clear();
	system.removeAdditionalOps();
	system.deallocate();
	site.deallocate(); 
      }

      //now take the expanded base wavefunction and canonicalize it 
    
      if (mpigetrank() == 0) {
	Sweep::InitializeStateInfoPartialSweep(sweepParams, true, baseStates[0]);
	Sweep::CanonicalizeWavefunctionPartialSweep(sweepParams, false, baseStates[0]);
	Sweep::CanonicalizeWavefunctionPartialSweep(sweepParams, true, baseStates[0]);
	Sweep::CanonicalizeWavefunctionPartialSweep(sweepParams, false, baseStates[0]);
      }
      
      //store the first block for all combinations
      StackSpinBlock forwardStart;
      InitBlocks::InitStartingBlock (forwardStart,true, targetState, targetState,
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     0, false, false, 0);
      sites[0] = 0; sites[1] = 0;
      StackSpinBlock::store(true, sites, forwardStart, targetState, targetState);  
      
      forwardStart.deallocate(); forwardStart.clear();
      InitBlocks::InitStartingBlock (forwardStart,true, baseStates[0], baseStates[0],
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     0, false, false, 0);
      sites[0] = 0; sites[1] = 0;
      StackSpinBlock::store(true, sites, forwardStart, baseStates[0], baseStates[0]);  
      
      forwardStart.deallocate(); forwardStart.clear();
      InitBlocks::InitStartingBlock (forwardStart,true, targetState, baseStates[0],
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     0, false, false, 1);
      sites[0] = 0; sites[1] = 0;
      StackSpinBlock::store(true, sites, forwardStart, targetState, baseStates[0]);  
      forwardStart.set_integralIndex() = 0;
      StackSpinBlock::store(true, sites, forwardStart, targetState, baseStates[0]);  
      forwardStart.deallocate(); forwardStart.clear();
      direction = false;
      
    }
    else if (dmrginp.calc_type() == RESPONSEAAAC) {
      if (mpigetrank() == 0)
	StackWavefunction::ChangeLastSite(dmrginp.last_site()-1, dmrginp.getPartialSweep()-1, baseStates[0]); 
      
      int start = dmrginp.getPartialSweep();
      int end   = dmrginp.last_site();
      StackSpinBlock site(start-1, start-1, 1, false);      
      StackSpinBlock system = StackSpinBlock::buildBigEdgeBlock(start, end, true, true, 1, false);
      
      system.addAdditionalOps(); system.set_loopblock(false);
      StackSpinBlock newSystem;
      newSystem.default_op_components(false, true, true, false);
      newSystem.set_integralIndex() = 1;
      newSystem.setstoragetype(DISTRIBUTED_STORAGE);
      newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, site);
      
      
      vector<int> sites(2); sites[0] = start-1; sites[1] = start-1;
      if (mpigetrank() == 0 ) {
	std::vector<SpinQuantum> quanta = newSystem.get_stateInfo().quanta;
	std::vector<Matrix> rotation(quanta.size());
	
	StateInfo lsi = *newSystem.get_stateInfo().leftStateInfo;
	StateInfo rsi = *newSystem.get_stateInfo().rightStateInfo;
	StateInfo si = newSystem.get_stateInfo();
	
	for (int i=0; i<si.quanta.size(); i++) {
	  const vector<int>& oldToNewStateI = si.oldToNewState[i];
	  int rowindex = 0;
	  for (int iSub =0; iSub < oldToNewStateI.size(); iSub++) {
	    int unCollectedI = oldToNewStateI[iSub];
	    int lindex = si.unCollectedStateInfo->leftUnMapQuanta[unCollectedI];
	    int rindex = si.unCollectedStateInfo->rightUnMapQuanta[unCollectedI];
	    
	    if (lindex == lsi.quanta.size()-1) {//this is the last index so all orbs are doubly occupied
	      if (rotation[i].Ncols() != 0) {
		cout << "something wrong in writing the matrix"<<endl;
		abort();
	      }
	      else {
		rotation[i].ReSize(si.quantaStates[i],1);rotation[i] = 0.0;
		rotation[i](rowindex+1,1) = 1.0;
	      }
	      break;
	    }
	    rowindex += si.unCollectedStateInfo->quantaStates[unCollectedI];	    
	  }
	}
	sites[0] = start-1; sites[1] = dmrginp.last_site()-1;
	SaveRotationMatrix(sites, rotation, baseStates[0]); 
      }

      //make edge blocks for all the combinations
      sites[0] = start; sites[1] = end-1;
      StackSpinBlock::store(false, sites, system, targetState, baseStates[0]);  
      //newSystem.deallocate();
      
      //system = StackSpinBlock::buildBigEdgeBlock(start, end, false, true, 0, false);
      system.set_integralIndex() = 0;
      StackSpinBlock::store(false, sites, system, targetState, baseStates[0]);  
      newSystem.deallocate();
      system.removeAdditionalOps();
      system.deallocate(); 
      system.clear();
      site.deallocate();
      
      system = StackSpinBlock::buildBigEdgeBlock(start, end, true, true, 0, true);
      StackSpinBlock::store(false, sites, system, targetState, targetState);  
      StackSpinBlock::store(false, sites, system, baseStates[0], baseStates[0]);  
      system.removeAdditionalOps();
      system.deallocate();
      system.clear();
      
      
      
      //store the first block for all combinations
      StackSpinBlock forwardStart;
      InitBlocks::InitStartingBlock (forwardStart,true, targetState, targetState,
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     0, false, false, 0);
      sites[0] = 0; sites[1] = 0;
      StackSpinBlock::store(true, sites, forwardStart, targetState, targetState);  
      
      forwardStart.deallocate(); forwardStart.clear();
      InitBlocks::InitStartingBlock (forwardStart,true, baseStates[0], baseStates[0],
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     0, false, false, 0);
      sites[0] = 0; sites[1] = 0;
      StackSpinBlock::store(true, sites, forwardStart, baseStates[0], baseStates[0]);  
      
      forwardStart.deallocate(); forwardStart.clear();
      InitBlocks::InitStartingBlock (forwardStart,true, targetState, baseStates[0],
				     sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), 
				     0, false, false, 1);
      sites[0] = 0; sites[1] = 0;
      StackSpinBlock::store(true, sites, forwardStart, targetState, baseStates[0]);  
      forwardStart.set_integralIndex() = 0;
      StackSpinBlock::store(true, sites, forwardStart, targetState, baseStates[0]);  
      forwardStart.deallocate(); forwardStart.clear();
      
      //now take the expanded base wavefunction and canonicalize it 	
      dmrginp.setPartialSweep() = dmrginp.setPartialSweep()+1;
      sweepParams.current_root() = baseStates[0];
      if (mpigetrank() == 0) {
	Sweep::InitializeStateInfoPartialSweep(sweepParams, false, baseStates[0]);
	Sweep::CanonicalizeWavefunctionPartialSweep(sweepParams, true, baseStates[0]);
	Sweep::CanonicalizeWavefunctionPartialSweep(sweepParams, false, baseStates[0]);
	Sweep::CanonicalizeWavefunctionPartialSweep(sweepParams, true, baseStates[0]);
      }
      direction = false;
    }


    dmrginp.setGuessState() = baseStates[0];

    //now run a warmup sweep followed by actual sweeps
    last_be = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, projectors, baseStates);
    dmrginp.set_algorithm_method() = atype;
    last_fe = SweepResponse::do_one(sweepParams, !warmUp, !direction, restart, restartSize, targetState, projectors, baseStates);
    direction = !direction;
  }

  dmrginp.set_algorithm_method() = atype;
  restart = false;
  restartSize = 0;
  warmUp = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;

      last_be = SweepResponse::do_one(sweepParams, warmUp, !direction, restart, restartSize, targetState, projectors, baseStates);
      direction = !direction;

      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      last_fe = SweepResponse::do_one(sweepParams, warmUp, !direction, restart, restartSize, targetState, projectors, baseStates);
      direction = !direction;

      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }

  
}


void compress(double sweep_tol, int targetState, int baseState)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;
  bool direction;

  sweepParams.current_root() = -1;
  //this is the warmup sweep, the baseState is used as the initial guess for the targetState
  last_fe = SweepCompress::do_one(sweepParams, true, true, false, 0, targetState, baseState);

  direction = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = SweepCompress::do_one(sweepParams, false, false, false, 0, targetState, baseState);
      direction = true;
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      last_fe = SweepCompress::do_one(sweepParams, false, true, false, 0, targetState, baseState);
      direction = false;
      
      
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }

  //we finally canonicalize the targetState
  //one has to canonicalize the wavefunction with atleast 3 sweeps, this is a quirk of the way 
  //we transform wavefunction
  algorithmTypes atype;
  atype = dmrginp.algorithm_method();
  dmrginp.set_algorithm_method() = ONEDOT;
  if (mpigetrank()==0) {
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);  
  }
  dmrginp.set_algorithm_method() = atype;
  
}

void dmrg_stateSpecific(double sweep_tol, int targetState)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  int ls_count=0;
  SweepParams sweepParams;
  int old_states=sweepParams.get_keep_states();
  int new_states;
  double old_error=0.0;
  double old_energy=0.0;
  // warm up sweep ...

  bool direction;
  int restartsize;
  sweepParams.restorestate(direction, restartsize);

  //initialize array of size m_maxiter or dmrginp.max_iter() for dw and energy
  sweepParams.current_root() = targetState;

  last_fe = Sweep::do_one(sweepParams, false, direction, true, restartsize);

  while ((fabs(last_fe - old_fe) > sweep_tol) || (fabs(last_be - old_be) > sweep_tol)  )
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()) 
	break;

      last_be = Sweep::do_one(sweepParams, false, !direction, false, 0);
      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;

      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;


      last_fe = Sweep::do_one(sweepParams, false, direction, false, 0);

      new_states=sweepParams.get_keep_states();


      pout << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;

    }
  pout << "Converged Energy  " << sweepParams.get_lowest_energy()[0]<< std::endl;
  if(dmrginp.max_iter() <= sweepParams.get_sweep_iter()) {
    
    pout << "Maximum sweep iterations achieved " << std::endl;
  }

  //one has to canonicalize the wavefunction with atleast 3 sweeps, this is a quirk of the way 
  //we transform wavefunction
  if (mpigetrank()==0) {
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, targetState);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, !direction, targetState);
    Sweep::InitializeStateInfo(sweepParams, direction, targetState);
    
  }

}


void Npdm(int pdm, bool restartpdm, bool transitionpdm)
{
  double sweep_tol = 1e-7;
  sweep_tol = dmrginp.get_sweep_tol();
  bool direction;
  int restartsize;
  bool direction_copy; int restartsize_copy;
  SweepParams sweepParams;
  SweepParams sweep_copy;
  
  if(sym == "dinfh") {
    pout << "Npdm not implemented with dinfh symmetry"<<endl;
    abort();
  }
  
  if (dmrginp.algorithm_method() == TWODOT ) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
    abort();
  }
  
  dmrginp.do_pdm() = true;
  
  // Screening can break things for NPDM (e.g. smaller operators won't be available from which to build larger ones etc...?)
  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops
  dmrginp.Sz() = dmrginp.total_spin_number().getirrep();
  sweep_copy.restorestate(direction_copy, restartsize_copy);


  {
    Timer timer;
    sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
    
    for (int state=0; state<dmrginp.nroots(); state++) {
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, state, state); //this will generate the cd operators                               
      if (pdm == 1) SweepOnepdm::do_one(sweepParams, false, direction, false, 0, state);     
      //if (pdm == 2) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state);     
      //else if (npdm_order == NPDM_TWOPDM) SweepTwopdm::do_one(sweepParams, false, direction, false, 0, state, state);
      else abort();
    }
    
  }
  
  sweep_copy.savestate(direction_copy, restartsize_copy);
  
}


/*
void restartResponseSweep(double sweep_tol, int targetState, int correctionVector, int baseState)
{
  double last_fe = 10.e6;
  double last_be = 10.e6;
  double old_fe = 0.;
  double old_be = 0.;
  SweepParams sweepParams;
  bool direction, warmUp=false, restart=true;
  int restartSize=0;

  sweepParams.restorestate(direction, restartSize);

  sweepParams.current_root() = -1;

  last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, correctionVector, baseState);

  warmUp = false;
  restart = false;
  while ( true)
    {
      old_fe = last_fe;
      old_be = last_be;
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      last_be = SweepResponse::do_one(sweepParams, warmUp, !direction, restart, restartSize, targetState, correctionVector, baseState);
      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	break;
      
      direction = true;
      last_fe = SweepResponse::do_one(sweepParams, warmUp, direction, restart, restartSize, targetState, correctionVector, baseState);

      
      p1out << "\t\t\t Finished Sweep Iteration "<<sweepParams.get_sweep_iter()<<endl;
      
    }
  
}
*/

