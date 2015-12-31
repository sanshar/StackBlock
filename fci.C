/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include "Stackspinblock.h"
#include "Stackwavefunction.h"
#include "sweep.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace std;


void SpinAdapted::Sweep::fullci(double sweep_tol)
{
  int integralIndex = 0;
  SweepParams sweepParams;
  sweepParams.set_sweep_parameters();


  StackSpinBlock system, sysdot;
  InitBlocks::InitStartingBlock(system, true, 0, 0, sweepParams.get_forward_starting_size(),  sweepParams.get_backward_starting_size(), 0, false, true, integralIndex);
  int numsites = dmrginp.spinAdapted() ? dmrginp.last_site() : dmrginp.last_site()/2;
  int forwardsites = numsites/2+numsites%2;
  int backwardsites = numsites - forwardsites;
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));


  StackSpinBlock newSystem;
  for (int i=0; i<forwardsites-1; i++) {
    sysdot = StackSpinBlock(i+1, i+1, integralIndex, true);
    system.addAdditionalOps();
    newSystem.set_integralIndex() = integralIndex;
    if (i == forwardsites-2)
      newSystem.default_op_components(true, true, false, true);
    else
      newSystem.default_op_components(false, true, false, true);

    newSystem.setstoragetype(DISTRIBUTED_STORAGE);
    newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, sysdot);

    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysMemory = newSystem.memoryUsed();
    if (i != forwardsites-2) {
      if (i != 0) {
	newSystem.moveToNewMemory(system.getdata());
	Stackmem[0].deallocate(newSystem.getdata()+newSystem.memoryUsed(), memoryToFree);
      }
      system.clear();
      system = newSystem;
    }
  }

  StackSpinBlock environment, newEnvironment, envdot;
  InitBlocks::InitStartingBlock(environment, false, 0, 0, sweepParams.get_forward_starting_size(),  sweepParams.get_backward_starting_size(), 0, false, true, integralIndex);
  for (int i=0;i <backwardsites-1; i++) {
    envdot = StackSpinBlock(numsites-2-i, numsites-2-i, integralIndex, true);
    environment.addAdditionalOps();
    newEnvironment.set_integralIndex() = integralIndex;
    if (i == backwardsites-2)
      newEnvironment.default_op_components(true, false, true, true);
    else
      newEnvironment.default_op_components(false, false, true, true);
    newEnvironment.setstoragetype(DISTRIBUTED_STORAGE);
    newEnvironment.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, environment, envdot);

    if (i!=backwardsites-2) {
      if (i != 0) {
	long memoryToFree = newEnvironment.getdata() - environment.getdata();
	long newenvMemory = newEnvironment.memoryUsed();
	newEnvironment.moveToNewMemory(environment.getdata());
	Stackmem[0].deallocate(newEnvironment.getdata()+newEnvironment.memoryUsed(), memoryToFree);
      }
      environment.clear();
      environment = newEnvironment;
    }
  }

  pout <<"\t\t\t System Block :: "<< newSystem;
  pout <<"\t\t\t Environment Block :: "<< newEnvironment;
  newSystem.set_loopblock(true); newEnvironment.set_loopblock(false);
  StackSpinBlock big;
  InitBlocks::InitBigBlock(newSystem, newEnvironment, big); 


  int nroots = dmrginp.nroots(0);
  std::vector<StackWavefunction> solution(nroots);

  solution[0].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), false);
  solution[0].Clear();
  if (mpigetrank() == 0) {
    for (int i=1; i<nroots; i++) {
      solution[i].initialise(dmrginp.effective_molecule_quantum_vec(), big.get_leftBlock()->get_stateInfo(), big.get_rightBlock()->get_stateInfo(), false);
      solution[i].Clear();
    }
  }


  std::vector<double> energies(nroots);
  double tol = sweepParams.get_davidson_tol();

  pout << "\t\t\t Solving the Wavefunction "<<endl;
  int currentState = 0;
  std::vector<StackWavefunction> lowerStates;
  Solver::solve_wavefunction(solution, energies, big, tol, BASIC, false, true, false, false, sweepParams.get_additional_noise(), currentState, lowerStates);

  pout << "tensormultiply "<<*dmrginp.tensormultiply<<endl;

  for (int i=0; i<nroots; i++) {
    pout << "fullci energy "<< energies[i]<<endl;
  }
  if (!mpigetrank())
  {
#ifndef MOLPRO
    FILE* f = fopen("dmrg.e", "wb");
#else
    std::string efile;
    efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
    FILE* f = fopen(efile.c_str(), "wb");
#endif
    
    for(int j=0;j<nroots;++j) {
      double e = energies[j]; 
      fwrite( &e, 1, sizeof(double), f);
    }
    fclose(f);
  }


  if (mpigetrank() == 0) {
    for (int i=nroots-1; i>0; i--)
      solution[i].deallocate();
  }
  solution[0].deallocate();


}

void SpinAdapted::Sweep::tiny(double sweep_tol)
{
#ifndef SERIAL
  if(mpigetrank() == 0) {
#endif
    pout.precision(12);

  int nroots = dmrginp.nroots(0);
  SweepParams sweepParams;
  sweepParams.set_sweep_parameters();
  StackSpinBlock system(0,dmrginp.last_site()-1, 0, true);
  const StateInfo& sinfo = system.get_stateInfo();
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  for (int i=0; i<sinfo.totalStates; i++) {
    if (sinfo.quanta[i] == dmrginp.molecule_quantum()) {
      StackMatrix& h = system.get_op_array(HAM).get_element(0).at(0)->operator_element(i,i);
      DiagonalMatrix energies(h.Nrows()); energies = 0.0;
      diagonalise(h, energies);
      
      for (int x=0; x<nroots; x++) 
	pout << "fullci energy  "<< energies(x+1)<<endl;

      if (mpigetrank() == 0)
      {
#ifndef MOLPRO
	FILE* f = fopen("dmrg.e", "wb");
#else
	std::string efile;
	efile = str(boost::format("%s%s") % dmrginp.load_prefix() % "/dmrg.e" );
	FILE* f = fopen(efile.c_str(), "wb");
#endif
	
	for(int j=0;j<nroots;++j) {
	  double e = energies(j+1); 
	  fwrite( &e, 1, sizeof(double), f);
	}
	fclose(f);
      }

      return;
    }
  }

  pout << "The wavefunction symmetry is not possible with the orbitals supplied."<<endl;
  abort();
#ifndef SERIAL
  }
#endif
}

