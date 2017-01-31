#include "global.h"
#include "IntegralMatrix.h"
#include "fciqmchelper.h"
#include "input.h"
#include "Stackspinblock.h"
#include "wrapper.h"
#include "rotationmat.h"
#include <sstream>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "sweepgenblock.h"
#include "sweep.h"
#include "fciqmchelper.h"

void ReadInput(char* conf);
namespace SpinAdapted{
MPS globalMPS;
}
using namespace SpinAdapted;

void initBoostMPI(int argc, char* argv[]) {
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
#endif
}
/*
void AddMPSs(int* states, double* scale, int nstates, int outstate) {
  AddMPS(states, scale, nstates, outstate);
}

void seedRandom(int seed) {
  srand(seed);
}

void ApplyCD(int i, int j) {
  MPS m(0);
  m.ApplyCD(i,j);
}

void CollapseToDeterminant(char* s, int stateIndex) {
  MPS m(stateIndex);
  m.CollapseToDeterminant(s);
}
*/

void ReadInputFromC(char* conf, int outputlevel) {
  ReadInput(conf);
  dmrginp.setOutputlevel() = outputlevel;
  dmrginp.initCumulTimer();
  MAX_THRD = dmrginp.thrds_per_node()[mpigetrank()];
  int mkl_thrd = dmrginp.mkl_thrds();

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

  dmrginp.matmultFlops.resize(numthrds, 0.);
  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();

}

void initializeGlobalMPS(int mpsindex) {
  SpinAdapted::globalMPS = MPS(mpsindex);
}
/*
void writeToDisk(unsigned long &occ, int length, int stateIndex)
{
  MPS m(&occ, length, stateIndex);
  //m.writeToDiskForDMRG(stateIndex);
}
*/
void evaluateOverlapAndHamiltonian(int state1, int state2, double* o, double* h)
{
  calcHamiltonianAndOverlap(state1, state2, *h, *o);
}

void readMPSFromDiskAndInitializeStaticVariables(bool initializeDotBlocks) {

  if (dmrginp.spinAdapted() && dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 ) {
    int i=0;
    StackSpinBlock s(i, i, 0, false);
    SpinQuantum sq = dmrginp.molecule_quantum();
    sq = SpinQuantum(sq.get_s().getirrep(), sq.get_s(), IrrepSpace(0));
    int qs = 1, ns = 1;
    StateInfo addstate(ns, &sq, &qs); 
    StackSpinBlock dummyblock(addstate, 0);
    StackSpinBlock newstartingBlock;
    newstartingBlock.set_integralIndex() = 0;
    newstartingBlock.default_op_components(false, true, true, false);
    newstartingBlock.setstoragetype(LOCAL_STORAGE);
    newstartingBlock.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, s, dummyblock);
    if (mpigetrank() == 0)
      MPS::siteBlocks.push_back(newstartingBlock); //alway make transpose operators as well
  }
  else {
    int i = 0;
    if (mpigetrank() == 0)
      MPS::siteBlocks.push_back(StackSpinBlock(i, i, 0, false)); //alway make transpose operators as well
  }

  if (mpigetrank() == 0) {
    if(!dmrginp.spinAdapted())
      MPS::sweepIters = dmrginp.last_site()/2-2;
    else
      MPS::sweepIters = dmrginp.last_site()-2;
    MPS::spinAdapted = false;
    if (initializeDotBlocks) {
      for (int i=1; i<MPS::sweepIters+2; i++) {
	  MPS::siteBlocks.push_back(StackSpinBlock(i, i, 0, false)); //alway make transpose operators as well
      }
    }
  }
#ifndef SERIAL
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, MPS::sweepIters, 0);
  boost::mpi::broadcast(world, MPS::spinAdapted, 0);
#endif

}


/*
void RDM(char* infile)
{
  setbuf(stdout, NULL);
  int msgsize=1000;
  char msgctr[msgsize];

  int nstates;
  std::vector<int> states;
  if (mpigetrank() == 0) {
    ifstream file(infile);
    int stateindex ;
    while(file >> stateindex) {
      states.push_back(stateindex);
      if (mpigetrank() == 0)
	printf("reading state %i\n", stateindex);
    }
    file.close();
  }
#ifndef SERIAL
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, states, 0);
#endif
  nstates = states.size();

  //read the sweepparam
  SweepParams sweepParams, sweep_copy;
  bool direction, direction_copy; int restartsize, restartsize_copy;
  sweep_copy.restorestate(direction_copy, restartsize_copy);

  if (dmrginp.algorithm_method() == TWODOT) {
    pout << "Npdm not allowed with twodot algorithm" << endl;
    abort();
  }

  boost::shared_ptr<Npdm::Npdm_driver_base> npdm_driver = boost::shared_ptr<Npdm::Npdm_driver_base>( new Npdm::Twopdm_driver( dmrginp.last_site() ) );

  dmrginp.setimplicitTranspose() = false;
  dmrginp.setStateSpecific() = true;
  dmrginp.do_pdm() = true;
  dmrginp.oneindex_screen_tol() = 0.0; //need to turn screening off for one index ops
  dmrginp.twoindex_screen_tol() = 0.0; //need to turn screening off for two index ops

  //calculate regular twordm
  for (int i=0; i<nstates; i++) {
    for (int j=i; j<i+1; j++) {

      dmrginp.set_fullrestart() = true;

      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      if (mpigetrank() == 0) {
	Sweep::InitializeStateInfo(sweepParams, direction, states[i]);
	Sweep::InitializeStateInfo(sweepParams, !direction, states[i]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, states[i]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, states[i]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, states[i]);
      }
      
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      if (mpigetrank() == 0) {
	Sweep::InitializeStateInfo(sweepParams, direction, states[j]);
	Sweep::InitializeStateInfo(sweepParams, !direction, states[j]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, states[j]);
	Sweep::CanonicalizeWavefunction(sweepParams, !direction, states[j]);
	Sweep::CanonicalizeWavefunction(sweepParams, direction, states[j]);
      }

      dmrginp.do_npdm_ops() = true;
      sweepParams = sweep_copy; direction = direction_copy; restartsize = restartsize_copy;
      //args (sweeparams, warmup, forward, restart, restartsize, statea, stateb)
      SweepGenblock::do_one(sweepParams, false, !direction, false, 0, states[i], states[j]); //this will generate the cd operators
      dmrginp.set_fullrestart() = false;


      //npdm_do_one_sweep(*npdm_driver, sweepParams, false, direction, false, 0, states[i], states[j]);
      SweepTwopdm::do_one(sweepParams, false, direction, false, 0, states[i], states[j]);
    }
  }


}
*/

void test(char* infile)
{
  setbuf(stdout, NULL);
  pout.precision(12);
  int msgsize=1000;
  char msgctr[msgsize];

  int nstates;
  std::vector<int> states;
  if (mpigetrank() == 0) {
    ifstream file(infile);
    int stateindex ;
    while(file >> stateindex) {
      states.push_back(stateindex);
      if (mpigetrank() == 0)
	printf("reading state %i\n", stateindex);
    }
    file.close();
  }
#ifndef SERIAL
  boost::mpi::communicator world;
  boost::mpi::broadcast(world, states, 0);
#endif
  nstates = states.size();


  std::vector< std::vector<double> > ham(nstates, std::vector<double>(nstates, 0.0));
  std::vector< std::vector<double> > Overlap(nstates, std::vector<double>(nstates, 0.0));


  for (int i=0; i<nstates; i++) {
    if(mpigetrank() == 0)
      printf("starting row : %i\n", i);
    for (int j=0; j<=i; j++) {
      double h=0,o=0;
      calcHamiltonianAndOverlap(states[i], states[j], h, o);
      ham[i][j] = h; ham[j][i] = h;
      Overlap[i][j] = o; Overlap[j][i] = o;
      if (mpigetrank() == 0) 
	printf("%i %i  %18.9e  %18.9e\n", i, j, h, o); 
    }
  }
  
  if(mpigetrank() == 0) {
    printf("printing hamiltonian\n");
    for (int i=0; i<nstates; i++) {
      for (int j=0; j<nstates; j++) 
	printf("%18.9e ", ham[i][j]);
      printf("\n");
    }

    /*
    printf("\n");
    printf("printing hamiltonian\n");
    for (int i=0; i<nstates; i++) {
      for (int j=0; j<nstates; j++) 
	printf("%18.9e ", ham[i][j]/sqrt(Overlap[i][i]*Overlap[j][j]));
      printf("\n");
    }
    */

    printf("\n");
    printf("printing overlap\n");
    for (int i=0; i<nstates; i++) {
      for (int j=0; j<nstates; j++) 
	printf("%18.9e ", Overlap[i][j]);
      printf("\n");
    }
  }
}

/*
void evaluateOverlapAndHamiltonian1(unsigned long *occ, int length, double* o, double* h) {
  MPS dmrgc(occ, length);
  dmrgc.writeToDiskForDMRG(2, false);
  calcHamiltonianAndOverlap(0, 2, *h, *o);
}
*/

void intFromString(unsigned long &occ, char* s) {
  occ = 0;
  long temp = 1;
  string ss(s);
  stringstream stream(ss);
  int n, i=0;
  while (stream >>n) {
    if (n==1)
      occ = occ | temp <<(63-i);
    i++;
  }
  return;
}
