#include "stackguess_wavefunction.h"
#include "global.h"
#include "solver.h"
#include "initblocks.h"
#include "MatrixBLAS.h"
#include <boost/format.hpp>
#ifndef SERIAL
#include "mpi.h"
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif
#include "rotationmat.h"
#include "Stackspinblock.h"
#include "Stackdensity.h"
#include "pario.h"
#include "davidson.h"
#include "sweep.h"
//#include "mps_mpi.h"
#include "type1.h"
#include "operatorfunctions.h"
#include "StateInfo.h"
#include <stdio.h>
#include <boost/filesystem.hpp>
#include "IntegralMatrix.h"

using namespace boost;
using namespace std;

double    block_time      ;
double    prepare_time    ;
double    multiply_time   ;
double    rotate_time     ;
double    dellocate_time  ;

namespace SpinAdapted{
std::vector<StackSpinBlock> siteBlocks_noDES;
}

double SpinAdapted::mps_nevpt::type1::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, perturber& pb, int baseState)
{
  int integralIndex = 0;
  StackSpinBlock system;
  system.nonactive_orb() = pb.orb();
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  std::vector<double> finalEnergy(nroots,-1.0e10);
  std::vector<double> finalEnergy_spins(nroots,0.);
  double finalError = 0.;

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  if (forward)
    p3out << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl;
  else
    p3out << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl;
  p3out << "\t\t\t ============================================================================ " << endl;

  InitBlocks::InitStartingBlock (system,forward, baseState, pb.wavenumber(), sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex, pb.braquanta, pb.ketquanta);
  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  p3out << "\t\t\t Starting block is :: " << endl << system << endl;

  StackSpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (restart)
  {
    if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
      dot_with_sys = false;
    if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
      dot_with_sys = false;
  }
  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      p3out << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      p3out << "\t\t\t ----------------------------" << endl;
	    if (forward)
      {
	      p3out << "\t\t\t Current direction is :: Forwards " << endl;
      }
	    else
      {
	      p3out << "\t\t\t Current direction is :: Backwards " << endl;
      }

      if (sweepParams.get_block_iter() != 0) 
	sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
       p3out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem; // new system after blocking and decimating
      newSystem.nonactive_orb() = pb.orb();

      //Need to substitute by:
     // if (warmUp )
     //   Startup(sweepParams, system, newSystem, dot_with_sys, pb.wavenumber(), baseState);
     // else {
     //   BlockDecimateAndCompress (sweepParams, system, newSystem, false, dot_with_sys, pb.wavenumber(), baseState);
     // }
      
      Timer timer;
    timer.start();
        BlockDecimateAndCompress (sweepParams, system, newSystem, warmUp, dot_with_sys,pb, baseState);
      //Need to substitute by?
    block_time += timer.elapsedcputime();


      system = newSystem;
      if (dmrginp.outputlevel() > 2){
	    p3out << system<<endl;
	    p3out << system.get_braStateInfo()<<endl;
	    system.printOperatorSummary();
      }
      
      //system size is going to be less than environment size
      if (forward && system.get_complementary_sites()[0] >= dmrginp.last_site()/2)
	    dot_with_sys = false;
      if (!forward && system.get_sites()[0]-1 < dmrginp.last_site()/2)
	    dot_with_sys = false;

      StackSpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState);	 	
      syssites = system.get_sites();
	    p3out << "\t\t\t saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      
#ifndef SERIAL
      mpi::communicator world;
      world.barrier();
#endif
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }

  //FIXME
  //It does not seem necessary.

  //when we are doing twodot, we still need to do the last sweep to make sure that the
  //correctionVector and base wavefunction are propogated correctly across sweeps
//  //especially when we switch from twodot to onedot algorithm
//  if (!sweepParams.get_onedot() && !warmUp) {
//      pout << "\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
//      pout << "\t\t\t ----------------------------" << endl;
//      if (dmrginp.outputlevel() > 0) {
//	    if (forward)
//	      pout << "\t\t\t Current direction is :: Forwards " << endl;
//	    else
//	      pout << "\t\t\t Current direction is :: Backwards " << endl;
//      }
//    sweepParams.set_onedot() = true;
//    sweepParams.set_env_add() = 0;
//    bool dot_with_sys = true;
//    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, baseState);
//    sweepParams.set_onedot() = false;
//    sweepParams.set_env_add() = 1;
//  }
//

  p3out << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  p3out << "\t\t\t Largest overlap for Sweep with " << sweepParams.get_keep_states() << " states is " << finalEnergy[0] << endl;
  sweepParams.set_largest_dw() = finalError;
  

  p3out << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalError;
}

void SpinAdapted::mps_nevpt::type1::BlockDecimateAndCompress (SweepParams &sweepParams, StackSpinBlock& system, StackSpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys, perturber& pb, int baseState)
{
  Timer timer;
    timer.start();
  int sweepiter = sweepParams.get_sweep_iter();
  if (dmrginp.outputlevel() > 0) {
    mcheck("at the start of block and decimate");
    p3out << "\t\t\t dot with system "<<dot_with_sys<<endl;
    p3out <<endl<< "\t\t\t Performing Blocking"<<endl;
  }
  // figure out if we are going forward or backwards
  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot;
  StackSpinBlock environment, environmentDot, newEnvironment;
  StackSpinBlock big;
  environment.nonactive_orb() = pb.orb();
  environmentDot.nonactive_orb() = pb.orb();
  newEnvironment.nonactive_orb() = pb.orb();
  int systemDotStart, systemDotEnd;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() -1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, 0, pb.orb(), false);
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, 0, pb.orb(), false);

  Sweep::makeSystemEnvironmentBigBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), pb.wavenumber(), baseState,pb.braquanta,pb.ketquanta);


  //analyse_operator_distribution(big);
  dmrginp.guessgenT -> stop();
  dmrginp.multiplierT -> start();
  std::vector<Matrix> rotatematrix;

  if (dmrginp.outputlevel() > 2)
    mcheck(""); 
  if (dmrginp.outputlevel() > 2) {
    if (!dot_with_sys && sweepParams.get_onedot())  { p3out << "\t\t\t System  Block"<<system;    }
    else p3out << "\t\t\t System  Block"<<newSystem;
    p3out << "\t\t\t Environment Block"<<newEnvironment<<endl;
    p3out << "\t\t\t Solving wavefunction "<<endl;
  }
  p3out <<"New System"<<endl;
  newSystem.printOperatorSummary();

  std::vector<StackWavefunction> solution; solution.resize(1);
  std::vector<StackWavefunction> outputState; outputState.resize(1);

  DiagonalMatrix e;


  //read the 0th wavefunction which we keep on the ket side because by default the ket stateinfo is used to initialize wavefunction
  //also when you use spinblock operators to multiply a state, it does so from the ket side i.e.  H|ket>
  //GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 0.0, baseState); 
	solution[0].initialise(pb.ketquanta, big.get_leftBlock()->get_ketStateInfo(), big.get_rightBlock()->get_ketStateInfo(), sweepParams.get_onedot());
	solution[0].Clear();
  GuessWave::guess_wavefunctions(solution[0], e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), baseState, dot_with_sys, 0.0); 

#ifndef SERIAL
  mpi::communicator world;
  broadcast(world, solution, 0);
#endif
  
  outputState[0].initialise(pb.braquanta,big.get_leftBlock()->get_braStateInfo(), big.get_rightBlock()->get_braStateInfo(),sweepParams.get_onedot());
  outputState[0].Clear();

    prepare_time +=timer.elapsedcputime();
    timer.start();
  if (pb.type() == TwoPerturbType::Va)
    big.multiplyCDD_sum(solution[0],&(outputState[0]),numthrds);
  if (pb.type() == TwoPerturbType::Vi)
    big.multiplyCCD_sum(solution[0],&(outputState[0]),numthrds);

    multiply_time +=timer.elapsedcputime();
  //davidson_f(solution[0], outputState[0]);
    timer.start();
  StackSpinBlock newbig;

  if (sweepParams.get_onedot() && !dot_with_sys)
  {
    InitBlocks::InitNewSystemBlock(system, systemDot, newSystem, baseState, pb.wavenumber(), systemDot.size(), dmrginp.direct(), system.get_integralIndex(), DISTRIBUTED_STORAGE, false, true,NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,pb.braquanta,pb.ketquanta);
    InitBlocks::InitBigBlock(newSystem, environment, newbig,pb.braquanta,pb.ketquanta); 

    StackWavefunction tempwave = outputState[0];
    GuessWave::onedot_shufflesysdot(big.get_braStateInfo(), newbig.get_braStateInfo(), outputState[0], tempwave);  
    outputState[0] = tempwave;

    tempwave = solution[0];
    GuessWave::onedot_shufflesysdot(big.get_ketStateInfo(), newbig.get_ketStateInfo(), solution[0], tempwave);  
    solution[0] = tempwave;

    big.get_rightBlock()->clear();
    big.clear();
  }
  else
    newbig = big;

//**********************
  StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo());
  bratracedMatrix.allocate(newSystem.get_braStateInfo());
  //FIXME
  //bratracedMatrix.makedensitymatrix(outputState, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
  bratracedMatrix.makedensitymatrix(outputState[0], newbig, 1.0);
  if (sweepParams.get_noise() > NUMERICAL_ZERO) {
    p3out << "adding noise  "<<trace(bratracedMatrix)<<"  "<<sweepiter<<"  "<<dmrginp.weights(sweepiter)[0]<<endl;

    double* backupData;
    if (mpigetrank() == 0) {
      double norm = trace(bratracedMatrix);
      backupData = Stackmem[omprank].allocate(bratracedMatrix.memoryUsed());
	    memset(backupData, 0, bratracedMatrix.memoryUsed()* sizeof(double));
      if (fabs(norm) > NUMERICAL_ZERO) {
	    DAXPY(bratracedMatrix.memoryUsed(), 1.0, bratracedMatrix.get_data(), 1, backupData, 1);
	    //DAXPY(bratracedMatrix.memoryUsed(), 1.0/norm, bratracedMatrix.get_data(), 1, backupData, 1);
      }
    }
    bratracedMatrix.Clear();
    //************************
    bratracedMatrix.add_onedot_noise(solution[0], newbig);

    if (mpigetrank() == 0) {
      double norm = trace(bratracedMatrix);
      if (fabs(norm) > 1e-10) {
	      DAXPY(bratracedMatrix.memoryUsed(), sweepParams.get_noise()/norm, bratracedMatrix.get_data(), 1, backupData, 1);
      }
      DCOPY(bratracedMatrix.memoryUsed(), &backupData[0], 1, bratracedMatrix.get_data(), 1);
      
      if (trace(bratracedMatrix) <1e-14) 
	      bratracedMatrix.SymmetricRandomise();
      
      p3out << "after noise  "<<trace(bratracedMatrix)<<"  "<<sweepParams.get_noise()<<endl;
      Stackmem[omprank].deallocate(backupData, bratracedMatrix.memoryUsed());
    }
  }

  //****************************
//  
//  StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo());
//  //long requiredData = getRequiredMemory(newSystem.get_braStateInfo(), bratracedMatrix.get_deltaQuantum());
//  //std::vector<double> data(requiredData, 0.0);
//  //bratracedMatrix.allocate(newSystem.get_braStateInfo(), &data[0]);
//  bratracedMatrix.allocate(newSystem.get_braStateInfo());
//
//  ////bratracedMatrix.makedensitymatrix(outputState, newbig, dmrginp.weights(sweepiter), 0.0, 0.0, true);
//  //bratracedMatrix.makedensitymatrix(outputState, newbig, std::vector<double>(1,1.0), 0.0, 0.0, true);
//  //if (sweepParams.get_noise() > NUMERICAL_ZERO) {
//  //  pout << "adding noise  "<<trace(bratracedMatrix)<<"  "<<sweepiter<<"  "<<dmrginp.weights(sweepiter)[0]<<endl;
//    bratracedMatrix.add_onedot_noise(solution[0], newbig);
//  //  if (trace(bratracedMatrix) <1e-14) 
//  //    bratracedMatrix.SymmetricRandomise();
//  //    
//  //  pout << "after noise  "<<trace(bratracedMatrix)<<"  "<<sweepParams.get_noise()<<endl;
//  //}
//
//
  std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
  LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseState);

  double braerror;
  if (!mpigetrank()) {
    braerror = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
  }

#ifndef SERIAL
  broadcast(world, ketrotateMatrix, 0);
  broadcast(world, brarotateMatrix, 0);
#endif

  p3out << "\t\t\t Total bra discarded weight "<<braerror<<endl<<endl;

  sweepParams.set_lowest_error() = braerror;

  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, pb.wavenumber());
  //FIXME
  //It is neccessary for twodot algorithm to save baseState wavefuntion.
  //I do not know why. 
  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
  outputState[0].SaveWavefunctionInfo (newbig.get_braStateInfo(), newbig.get_leftBlock()->get_sites(), pb.wavenumber());

  bratracedMatrix.deallocate();
  outputState[0].deallocate();
  solution[0].deallocate();


    rotate_time +=timer.elapsedcputime();
    timer.start();
  newEnvironment.removeAdditionalOps();
  newEnvironment.deallocate();
  environment.removeAdditionalOps();
  environment.deallocate();

  //TODO 
  //Why do I need this?
  //They should have been consistent.
//  solution[0].SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), baseState);
//  SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), ketrotateMatrix, baseState);

  p1out <<"\t\t\t Performing Renormalization "<<endl;
  newSystem.transform_operators(brarotateMatrix, ketrotateMatrix);

  {
    long memoryToFree = newSystem.getdata() - system.getdata();
    long newsysmem = newSystem.memoryUsed();
    newSystem.moveToNewMemory(system.getdata());
    Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
    //system.clear();
  }

    dellocate_time +=timer.elapsedcputime();
  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;
  p2out << "addnoise  S_0_opxop  S_1_opxop"<<endl;
  p2out << *dmrginp.addnoise<<" "<<*dmrginp.s0time<<" "<<*dmrginp.s1time<<endl;

  p2out << str(boost::format("%-40s - %-10.4f\n") % "Total memory" % (Stackmem[0].size*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Memory used" % (Stackmem[0].memused*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->system" % ((system.memoryUsed()+system.additionalMemoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->newSystem" % ((newSystem.memoryUsed()+newSystem.additionalMemoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->Envrionment" % ((environment.memoryUsed()+environment.additionalMemoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->newEnvrionment" % ((newEnvironment.memoryUsed()+newEnvironment.additionalMemoryUsed())*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "  |-->Memory left" % ((Stackmem[0].size-Stackmem[0].memused)*8/1.e9));
  p2out << str(boost::format("%-40s - %-10.4f\n") % "      |-->wavefunction"  %((solution[0].memoryUsed())*8/1.e9) );

}

void SpinAdapted::mps_nevpt::type1::cleanup(int left, int right, int cleanlevel)
{
  /*
  for (int site=0; site < dmrginp.last_site(); site++)
  {
    file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% 0 % "." % site % "-states" % left % "." % right % "-integral" %0 % "rank" % mpigetrank() % 0 % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());
    file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% 0 % "." % site % "-states" % left % "." % right % "-integral" %0 % "rank" % mpigetrank() % 0 % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());
    file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% (site+1) % "-" % (dmrginp.last_site()-1)% "-states" % left % "." % right % "-integral" %0 % "rank" % mpigetrank() % 0 % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());
    file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% (site+1) % "-" % (dmrginp.last_site()-1)% "-states" % left % "." % right % "-integral" %0 % "rank" % mpigetrank() % 0 % ".tmp" );
    if (boost::filesystem::exists(file)) remove(file.c_str());

  }
  */
  boost::filesystem::path p(dmrginp.save_prefix());
  boost::filesystem::directory_iterator iter(p);
  std::vector<boost::filesystem::path> deletefiles;
  while(iter!=boost::filesystem::directory_iterator())
  {
    if(iter->path().filename().string().find("Block-")==0)
      deletefiles.push_back(*iter);
    iter++;
  }
  for(int i=0;i<deletefiles.size();i++)
    boost::filesystem::remove(deletefiles[i]);
}


void SpinAdapted::mps_nevpt::type1::subspace_Vi(int baseState)
{
  
  double energy=0;
  double overlap=0;
  ViPerturber pb;
  MPS::siteBlocks.clear();
  int coresize = dmrginp.spinAdapted()? dmrginp.core_size():dmrginp.core_size()*2;
  int coreshift = dmrginp.spinAdapted()? dmrginp.act_size(): dmrginp.act_size()*2;
  //int virtsize = dmrginp.virt_size();
  //int virtshift = dmrginp.core_size()+dmrginp.act_size();
  for(int i=0; i< coresize; i++){
    double perturberEnergy=0;
    dmrginp.set_calc_type() = MPS_NEVPT;
    pb.init(i+coreshift);
    pout << "Begin Vi subspace with i = " << pb.orb(0)<<endl;
    SweepParams sweepParams;
    sweepParams.set_sweep_parameters();
    sweepParams.current_root() = baseState;
    //sweepParams.current_root() = -1;
    //double last_fe = Startup(sweepParams, true, true, false, 0, pb, baseState);
    Timer timer;
    Startup(sweepParams, true, pb, baseState);
    pout <<"Start up time :" << timer.elapsedwalltime();
    //sweepParams.current_root() = baseState;
    timer.start();
    while(true)
    {
      do_one(sweepParams, false, false, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	      break;
      do_one(sweepParams, false, true, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	      break;
    }
    cleanup(pb.wavenumber(), baseState);
    p1out <<"Sweep time :" << timer.elapsedwalltime()<<endl;

p1out << block_time     <<endl ;
p1out << prepare_time   <<endl ;
p1out << multiply_time  <<endl ;
p1out << rotate_time    <<endl ;
p1out << dellocate_time <<endl ;
 block_time     =0.0 ;
 prepare_time   =0.0 ;
 multiply_time  =0.0 ;
 rotate_time    =0.0 ;
 dellocate_time =0.0 ;

    MPS pbmps(pb.wavenumber());
    double o, h;
    dmrginp.set_calc_type() = DMRG;

    timer.start();
    calcHamiltonianAndOverlap(pb.wavenumber(), pb.wavenumber() ,h, o, pb.braquanta);

    p1out <<"Calculate Expectation time :" << timer.elapsedwalltime()<<endl;;



    if(!dmrginp.spinAdapted())
    {
      //In nonspinAdapted, alpha and beta have the results. Only one is neccessary. 
      o*=2;
      h*=2;
      i++;
    }
    if(o> NUMERICAL_ZERO){
      double fock =dmrginp.spinAdapted()? v_1[0](2*(i+coreshift),2*(i+coreshift)): v_1[0](i+coreshift,i+coreshift);
      //perturberEnergy = h/o+fock+perturber::CoreEnergy[0];
      perturberEnergy = h/o - fock;
      energy += o/(perturber::ZeroEnergy[baseState]- perturberEnergy) ;
      //overlap +=o;
      overlap += o/((perturber::ZeroEnergy[baseState]- perturberEnergy)*(perturber::ZeroEnergy[baseState]- perturberEnergy));
      p1out << "Zero Energy: " << perturber::ZeroEnergy[baseState]<<endl;
      p1out << "Norm of contracted MPO on MPS : " << o <<endl;
      p1out << "Amplitude : " << o/((perturber::ZeroEnergy[baseState]- perturberEnergy)*(perturber::ZeroEnergy[baseState]- perturberEnergy)) <<endl;
      p1out << "Ener(only CAS part) : " << h/o<<endl;
      p1out << "Energy : " << perturberEnergy<<endl;
      p1out << "Correction Energy: "<< o/(perturber::ZeroEnergy[baseState]- perturberEnergy)<<endl; 
    }
    else{
      p1out << "Amplitude : " << 0.0 <<endl;
      p1out << "Energy : " << 0.0<<endl;
    }



  }
  pout << "Nevpt2 correction to the energy for state 0 in subspace Vi is " << energy<<endl;;
  pout << "Nevpt2 Vi subspace perturber Amplitude : " << overlap<<endl;;
  //pout << "Core Energy of nevpt2 " <<perturber::CoreEnergy[0]<<endl;
  std::string file = str(boost::format("%s%s%d") % dmrginp.load_prefix() % "/Vi_" % baseState);
  std::fstream f(file,std::fstream::out);
  f << energy <<endl;
  f << overlap <<endl;
  f.close();
}

void SpinAdapted::mps_nevpt::type1::subspace_Va(int baseState)
{
  
  double energy=0;
  double overlap=0;
  VaPerturber pb;
  MPS::siteBlocks.clear();
  int virtsize = dmrginp.spinAdapted()? dmrginp.virt_size():dmrginp.virt_size()*2;
  int virtshift = dmrginp.spinAdapted()? dmrginp.core_size()+dmrginp.act_size(): (dmrginp.core_size()+dmrginp.act_size())*2;
  //int virtsize = dmrginp.virt_size();
  //int virtshift = dmrginp.core_size()+dmrginp.act_size();
  for(int i=0; i< virtsize; i++){
    double perturberEnergy=0;
    dmrginp.set_calc_type() = MPS_NEVPT;
    pb.init(i+virtshift);
    pout << "Begin Va subspace with a = " << pb.orb(0)<<endl;
    SweepParams sweepParams;
    sweepParams.set_sweep_parameters();
    sweepParams.current_root() = baseState;
    //sweepParams.current_root() = -1;
    //double last_fe = Startup(sweepParams, true, true, false, 0, pb, baseState);
    Timer timer;
    Startup(sweepParams, true, pb, baseState);
    pout <<"Start up time :" << timer.elapsedwalltime();
    //sweepParams.current_root() = baseState;
    timer.start();
    while(true)
    {
      do_one(sweepParams, false, false, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	      break;
      do_one(sweepParams, false, true, false, 0, pb, baseState);
      if(dmrginp.max_iter() <= sweepParams.get_sweep_iter())
	      break;
    }
    cleanup(pb.wavenumber(), baseState);
    p1out <<"Sweep time :" << timer.elapsedwalltime()<<endl;
   
    MPS pbmps(pb.wavenumber());
    double o, h;
    dmrginp.set_calc_type() = DMRG;

    timer.start();

    calcHamiltonianAndOverlap(pb.wavenumber(), pb.wavenumber() ,h, o, pb.braquanta);

    p1out <<"Calculate Expectation time :" << timer.elapsedwalltime()<<endl;

    if(!dmrginp.spinAdapted())
    {
      //In nonspinAdapted, alpha and beta have the results. Only one is neccessary. 
      o*=2;
      h*=2;
      i++;
    }
    if(o> NUMERICAL_ZERO){
      double fock =dmrginp.spinAdapted()? v_1[0](2*(i+virtshift),2*(i+virtshift)): v_1[0](i+virtshift,i+virtshift);
      //perturberEnergy = h/o+fock+perturber::CoreEnergy[0];
      perturberEnergy = h/o+fock;
      energy += o/(perturber::ZeroEnergy[baseState]- perturberEnergy) ;
      //overlap +=o;
      overlap += o/((perturber::ZeroEnergy[baseState]- perturberEnergy)*(perturber::ZeroEnergy[baseState]- perturberEnergy));
      p1out << "Norm of contracted MPO on MPS : " << o <<endl;
      p1out << "Amplitude : " << o/((perturber::ZeroEnergy[baseState]- perturberEnergy)*(perturber::ZeroEnergy[baseState]- perturberEnergy)) <<endl;
      p1out << "Ener(only CAS part) : " << h/o<<endl;
      p1out << "Energy : " << perturberEnergy<<endl;
      p1out << "Correction Energy: "<< o/(perturber::ZeroEnergy[baseState]- perturberEnergy)<<endl; 
    }
    else{
      p1out << "Amplitude : " << 0.0 <<endl;
      p1out << "Energy : " << 0.0<<endl;
    }



  }
  pout << "Nevpt2 correction to the energy for state 0 in subspace Va is " << energy<<endl;;
  pout << "Nevpt2 Va subspace perturber Amplitude : " << overlap<<endl;;
  //pout << "Core Energy of nevpt2 " <<perturber::CoreEnergy[0]<<endl;
  std::string file = str(boost::format("%s%s%d") % dmrginp.load_prefix() % "/Va_" % baseState);
  std::fstream f(file,std::fstream::out);
  f << energy <<endl;
  f << overlap <<endl;
  f.close();
}


  void SpinAdapted::mps_nevpt::type1::calcHamiltonianAndOverlap(int statea, int stateb, double& h, double& o, const vector<SpinQuantum>& contracted_quanta, int integralIndex) {

    StackSpinBlock system, siteblock;
    bool forward = true, restart=false, warmUp = false;
    int forward_starting_size=1, backward_starting_size=0, restartSize =0;
    InitBlocks::InitStartingBlock(system, forward, statea, stateb, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, integralIndex, contracted_quanta, contracted_quanta); 

    p2out << system<<endl;

    std::vector<Matrix> Rotationa, Rotationb;
    /*
    int sysquanta = system.get_stateInfo().quanta.size();
    if (mpigetrank() == 0) {
      Rotationa = statea.getSiteTensors(0);
      Rotationb = stateb.getSiteTensors(0);
      Rotationa.resize(sysquanta);
      Rotationb.resize(sysquanta);
    }

#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world, Rotationa, 0);
    mpi::broadcast(world, Rotationb, 0);
#endif

    system.transform_operators(const_cast<std::vector<Matrix>&>(Rotationa), 
			       const_cast<std::vector<Matrix>&>(Rotationb));
    */
    int sys_add = true; bool direct = true; 


    std::vector<int> rotSites(2,0);
    int sweepIters = dmrginp.spinAdapted() ? dmrginp.last_site() -2 : dmrginp.last_site()/2-2;
    int normToComp = sweepIters/2;

    for (int i=0; i<sweepIters-1; i++) {
      pout << i<<" out of "<<sweepIters-1<<"  ";
      StackSpinBlock newSystem;
      if (i>=normToComp && !system.has(CRE_DESCOMP))
	system.addAllCompOps();
      system.addAdditionalOps();

      StackSpinBlock dotsite(i+1, i+1, integralIndex, true);
      if (mpigetrank() == 0) {
	rotSites[1] = i+1;
	LoadRotationMatrix(rotSites, Rotationa, statea);
	LoadRotationMatrix(rotSites, Rotationb, stateb);
	//Rotationa = statea.getSiteTensors(i+1);
	//Rotationb = stateb.getSiteTensors(i+1);
      }

#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(calc, Rotationa, 0);
      mpi::broadcast(calc, Rotationb, 0);
#endif
      if (i <normToComp)  {
	pout << "norm ops "<<endl;
	InitBlocks::InitNewSystemBlock(system, dotsite, newSystem, statea, stateb, sys_add, direct, integralIndex, DISTRIBUTED_STORAGE, true, false, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, contracted_quanta, contracted_quanta);
      }
      else {
	pout << "comp ops "<<endl;
	InitBlocks::InitNewSystemBlock(system, dotsite, newSystem, statea, stateb, sys_add, direct, integralIndex, DISTRIBUTED_STORAGE, false, true, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, contracted_quanta, contracted_quanta);
      }
      newSystem.transform_operators(const_cast<std::vector<Matrix>&>(Rotationa), 
				    const_cast<std::vector<Matrix>&>(Rotationb));
      {
	long memoryToFree = newSystem.getdata() - system.getdata();
	long newsysmem = newSystem.memoryUsed();
	newSystem.moveToNewMemory(system.getdata());
	Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
      }
      system = newSystem;
    }

    StackSpinBlock newSystem, big;
    StackSpinBlock dotsite1(sweepIters, sweepIters, integralIndex, true);
    StackSpinBlock dotsite2(sweepIters+1, sweepIters+1, integralIndex, true);

    //For molecule has at most 4 orbitals, there is at most one iteration.
    //System does not have CompOps.
	  system.addAllCompOps();

    system.addAdditionalOps();
    InitBlocks::InitNewSystemBlock(system, dotsite1, newSystem, statea, stateb, sys_add, direct, integralIndex, DISTRIBUTED_STORAGE, false, true, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, contracted_quanta, contracted_quanta);
    
    
    newSystem.set_loopblock(false); system.set_loopblock(false);
    InitBlocks::InitBigBlock(newSystem, dotsite2, big, contracted_quanta, contracted_quanta); 
    
    StackWavefunction stateaw, statebw;
    StateInfo s;

    rotSites[1] += 1;
    stateaw.LoadWavefunctionInfo(s, rotSites, statea, true);
    statebw.LoadWavefunctionInfo(s, rotSites, stateb, true);

#ifndef SERIAL
    mpi::broadcast(calc, stateaw, 0);
    mpi::broadcast(calc, statebw, 0);
#endif
    if (mpigetrank() != 0) {
      double* dataa = Stackmem[omprank].allocate(stateaw.memoryUsed());
      stateaw.set_data(dataa);
      stateaw.allocateOperatorMatrix();
      double* datab = Stackmem[omprank].allocate(statebw.memoryUsed());
      statebw.set_data(datab);
      statebw.allocateOperatorMatrix();
    }
#ifndef SERIAL
    calc.barrier();
    MPI_Bcast(stateaw.get_data(), stateaw.memoryUsed(), MPI_DOUBLE, 0, Calc);
    MPI_Bcast(statebw.get_data(), statebw.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

    StackWavefunction temp; temp.initialise(stateaw);
    temp.Clear();
    

    big.multiplyH_2index(statebw, &temp, 1);

    if (mpigetrank() == 0)
      h = DotProduct(stateaw, temp);

    temp.Clear();

    big.multiplyOverlap(statebw, &temp, 1);
    if (mpigetrank() == 0)
      o = DotProduct(stateaw, temp);

#ifndef SERIAL
      mpi::communicator world;
    mpi::broadcast(calc, h, 0);
    mpi::broadcast(calc, o, 0);
#endif
    temp.deallocate();
    statebw.deallocate();
    stateaw.deallocate();

    return;
  }

/*
void SpinAdapted::mps_nevpt::type1::calcHamiltonianAndOverlap(const MPS& statea, double& h, double& o, perturber& pb) {
#ifndef SERIAL
  mpi::communicator world;
#endif


  StackSpinBlock system, siteblock;
  //system.nonactive_orb() =pb.orb();
  //siteblock.nonactive_orb() =pb.orb();
  bool forward = true, restart=false, warmUp = false;
  int leftState=0, rightState=0, forward_starting_size=1, backward_starting_size=1, restartSize =0;
  InitBlocks::InitStartingBlock(system, forward, leftState, rightState, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0,statea.getw().get_deltaQuantum(), statea.getw().get_deltaQuantum()); 

  if (dmrginp.outputlevel() > 2)
  {
    p3out << system<<endl;
    system.printOperatorSummary();
  }
  //system.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(0)), 
  //    		       const_cast<std::vector<Matrix>&>(statea.getSiteTensors(0)), false, false );

  int sys_add = true; bool direct = true; 
  //int sys_add = true; bool direct = false; 

  int normToComp = MPS::sweepIters/2;
  for (int i=0; i<MPS::sweepIters-1; i++) {
    StackSpinBlock newSystem;
  if (i>=normToComp && !system.has(CRE_DESCOMP))
	    system.addAllCompOps();
     system.addAdditionalOps();
    //TODO
    //After using CreDes operator at fist, results were wrong.
    if (i <normToComp && false)  {
	    InitBlocks::InitNewSystemBlock(system, MPS::siteBlocks_noDES[i+1], newSystem, leftState, rightState, sys_add, direct, 0, DISTRIBUTED_STORAGE, true, false, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum());
    }
    else {
	    InitBlocks::InitNewSystemBlock(system, MPS::siteBlocks_noDES[i+1], newSystem, leftState, rightState, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum());
    }
    p3out << newSystem<<endl;
    newSystem.printOperatorSummary();

    //InitBlocks::InitNewSystemBlock(system, MPS::siteBlocks_noDES[i+1], newSystem, leftState, rightState, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true,NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum());

    //newSystem.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)));
    newSystem.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), 
      			    const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), false);
    //newSystem.transform_operators(const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), 
    //  			    const_cast<std::vector<Matrix>&>(statea.getSiteTensors(i+1)), false );

      {
	      long memoryToFree = newSystem.getdata() - system.getdata();
	      long newsysmem = newSystem.memoryUsed();
	      newSystem.moveToNewMemory(system.getdata());
	      Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
      }
    system = newSystem;
  }

  StackSpinBlock newSystem, big;
  system.addAdditionalOps();
  //system.addAdditionalCompOps();
  //system.printOperatorSummary();
  //To set implicit_transpose to true.
  //The last site spinblock should have implicit_transpose true.
  InitBlocks::InitNewSystemBlock(system, MPS::siteBlocks_noDES[MPS::sweepIters], newSystem,  leftState, rightState, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true,NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum());
  
  newSystem.set_loopblock(false); system.set_loopblock(false);
  //newSystem.addAdditionalCompOps();
  //MPS::siteBlocks_noDES[MPS::sweepIters+1].set_loopblock(false);
  InitBlocks::InitBigBlock(newSystem, MPS::siteBlocks_noDES[MPS::sweepIters+1], big,statea.getw().get_deltaQuantum(),statea.getw().get_deltaQuantum()); 
  

  //FIXME
  //Assume statea.getw() has only one deltaquantum.
  //Spin Embeding for zero order wavefunction is needed.

  StackWavefunction temp;
  temp.initialise(statea.getw());
  temp.Clear();
  big.multiplyH_2index(const_cast<StackWavefunction&>(statea.getw()), &temp, 1);
  //big.multiplyH(const_cast<StackWavefunction&>(statea.getw()), &temp, 1);

  if (mpigetrank() == 0)
    h = DotProduct(statea.getw(), temp);

  temp.Clear();
  big.multiplyOverlap(const_cast<StackWavefunction&>(statea.getw()), &temp, 1);
  if (mpigetrank() == 0)
    o = DotProduct(statea.getw(), temp);
  if(dmrginp.spinAdapted())
  {

    double cg= 0.0;
    //TODO
    //Assume the zero order wavefunction has a spin zero. 
    //Spin Embeding must be used.
    SpinQuantum wQ= statea.getw().get_deltaQuantum(0);
    //cg*= dmrginp.get_ninej()(wQ.get_s().getirrep(), 1, dmrginp.effective_molecule_quantum().get_s().getirrep(), 
    //      					   pb.delta.get_s().getirrep(), pb.delta.get_s().getirrep(), 0,
    //                          dmrginp.effective_molecule_quantum().get_s().getirrep(), 0, dmrginp.effective_molecule_quantum().get_s().getirrep());
    //cg*= Symmetry::spatial_ninej(wQ.get_symm().getirrep(), -pb.delta.get_symm().getirrep(), dmrginp.effective_molecule_quantum().get_symm().getirrep(), 
    //      					   pb.delta.get_symm().getirrep(), -pb.delta.get_symm().getirrep(), 0,
    //                          dmrginp.effective_molecule_quantum().get_symm().getirrep(), 0, dmrginp.effective_molecule_quantum().get_symm().getirrep());
    cg += pow(clebsch(wQ.get_s().getirrep(),-1,1,1,0,0),2);
    cg += pow(clebsch(wQ.get_s().getirrep(),1,1,-1,0,0),2);
    //cout << "cg coefficient: " <<cg<<endl;
    h*= cg*cg;
    o*= cg*cg;
  }


#ifndef SERIAL
  mpi::broadcast(world, h, 0);
  mpi::broadcast(world, o, 0);
#endif

  cleanup(pb.wavenumber(), pb.wavenumber());
  return;
}

*/

void SpinAdapted::mps_nevpt::type1::Startup(const SweepParams &sweepParams, const bool &forward, perturber& pb, int baseState) {

#ifndef SERIAL
  mpi::communicator world;
#endif
  assert(forward);
  StackSpinBlock system;
  system.nonactive_orb() =pb.orb();
  bool restart=false, warmUp = false;
  int forward_starting_size=1, backward_starting_size=0, restartSize =0;
  InitBlocks::InitStartingBlock(system, forward, pb.wavenumber(), baseState, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0,pb.braquanta, pb.ketquanta); 

  StackSpinBlock::store (forward, system.get_sites(), system, pb.wavenumber(), baseState); // if restart, just restoring an existing block --

  for (int i=0; i<MPS::sweepIters; i++) {
    StackSpinBlock newSystem;
    StackSpinBlock dotSystem(i+1,i+1, 0, pb.orb(), false);


    //system.addAdditionalCompOps();
    //newSystem.default_op_components(true, system, dotSystem, true, true, false);
    newSystem.nevpt_op_components(dmrginp.direct(), system, dotSystem, pb);
    newSystem.setstoragetype(DISTRIBUTED_STORAGE);
    newSystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, dotSystem, pb.braquanta, pb.ketquanta);
    newSystem.printOperatorSummary();

    //TODO
    //Do not restore environment, becasue it contains so many uncessary operators. 
    //Waste of memory. 
    //StackSpinBlock Environment, big;
    //StackSpinBlock::restore (!forward, newSystem.get_complementary_sites() , Environment, baseState, baseState);
    //StackSpinBlock::restore (!forward, newSystem.get_complementary_sites() , Environment,sweepParams.current_root(),sweepParams.current_root());

    //big.BuildSumBlock(PARTICLE_SPIN_NUMBER_CONSTRAINT, newSystem, Environment, pb.braquanta, pb.ketquanta);

    //StateInfo envStateInfo;
    StateInfo ketStateInfo;
    StateInfo braStateInfo;
    StateInfo halfbraStateInfo;// It has the same left and right StateInfo as braStateInfo. However, its total quanta is pb.ketquanta.
    // It is used to project solution into to braStateInfo.

    StackWavefunction solution; 
    StackWavefunction outputState; 
    //StackWavefunction *outputState_array;

    //initiateMultiThread(outputState, outputState_array, numthrds);
    StackWavefunction solutionprojector; 
    solution.LoadWavefunctionInfo(ketStateInfo, newSystem.get_sites(), baseState, true);
    //StateInfo::restore(!forward, newSystem.get_complementary_sites(), ketStateInfo, baseState);
    #ifndef SERIAL
      broadcast(world, ketStateInfo, 0);
      broadcast(world, solution, 0);
    #endif
    outputState.initialise(pb.braquanta,newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), solution.get_onedot());
    outputState.Clear();
    solutionprojector.initialise(pb.ketquanta,newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), solution.get_onedot());
    solutionprojector.Clear();
    //TensorProduct (newSystem.get_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.braquanta[0], EqualQ, braStateInfo);
    //TODO
    //TensorProduct do not support const StateInfo&
    TensorProduct (newSystem.set_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.braquanta[0], EqualQ, braStateInfo);
    TensorProduct (newSystem.set_braStateInfo(), *(ketStateInfo.rightStateInfo), pb.ketquanta[0], EqualQ, halfbraStateInfo);

    //StateInfo::restore(forward, environmentsites, envStateInfo, baseState);

    //DiagonalMatrix e;
    //if(i == 0)
    //  GuessWave::guess_wavefunctions(solution, e, big, TRANSPOSE, true, true, 0.0, baseState); 
    //else
    //  GuessWave::guess_wavefunctions(solution, e, big, TRANSFORM, true, true, 0.0, baseState); 


    //SpinAdapted::operatorfunctions::Product(&newSystem, ccd, solution[0], &ketStateInfo, stateb.getw(), temp, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)), true, 1.0);

    

    boost::shared_ptr<StackSparseMatrix> O;
    if (pb.type() == TwoPerturbType::Va)
      O = newSystem.get_op_array(CDD_SUM).get_element(0).at(0);
    if (pb.type() == TwoPerturbType::Vi)
      O = newSystem.get_op_array(CCD_SUM).get_element(0).at(0);
    boost::shared_ptr<StackSparseMatrix> overlap = newSystem.get_op_array(OVERLAP).get_element(0).at(0);
    bool deallocate1 = O->memoryUsed() == 0 ? true : false; 
    if (deallocate1) {
      O->allocate(newSystem.get_braStateInfo(), newSystem.get_ketStateInfo());
      O->build(newSystem);
    }
    bool deallocate2 = overlap->memoryUsed() == 0 ? true : false; 
    if (deallocate2) {
      overlap->allocate(newSystem.get_braStateInfo(), newSystem.get_ketStateInfo());
      overlap->build(newSystem);
    }
    //cout <<"leftbra: "<< newSystem.get_braStateInfo()<<endl;
    //cout <<"right: "<< *(ketStateInfo.rightStateInfo)<<endl;
    //cout <<"leftket: "<< newSystem.get_ketStateInfo()<<endl;
    //cout <<"brastateinfo: " <<braStateInfo<<endl;
    //cout <<"halfbrastateinfo: " <<halfbraStateInfo<<endl;
    //cout <<"ketstateinfo: " <<ketStateInfo<<endl;
    SpinAdapted::operatorfunctions::TensorMultiply(*O, &braStateInfo, &ketStateInfo , solution, outputState, pb.delta, true, 1.0);
    SpinAdapted::operatorfunctions::TensorMultiply(*overlap, &halfbraStateInfo, &ketStateInfo , solution, solutionprojector, overlap->get_deltaQuantum(0), true, 1.0);
    if (deallocate2) overlap->deallocate();
    if (deallocate1) O->deallocate();
    //cout <<"outputState\n"<<outputState<<endl;
    //cout <<"Projector\n"<<solutionprojector<<endl;

    StackDensityMatrix bratracedMatrix(newSystem.get_braStateInfo());
    //long requiredData = getRequiredMemory(newSystem.get_braStateInfo(), bratracedMatrix.get_deltaQuantum());
    //std::vector<double> data(requiredData, 0.0);
    //bratracedMatrix.allocate(newSystem.get_braStateInfo(), &data[0]);
    bratracedMatrix.allocate(newSystem.get_braStateInfo());
    double norm = DotProduct(outputState, outputState);
    if(norm > NUMERICAL_ZERO)
      SpinAdapted::operatorfunctions::MultiplyWithOwnTranspose(outputState, bratracedMatrix, 0.5/norm);
    SpinAdapted::operatorfunctions::MultiplyWithOwnTranspose(solutionprojector, bratracedMatrix, 0.5);
    std::vector<Matrix> brarotateMatrix, ketrotateMatrix;
    //cout <<"bratracedMatrix"<<endl;
    //cout <<bratracedMatrix<<endl;
    //cout <<"keep_states"<<sweepParams.get_keep_states()<<"keep_qstates"<<sweepParams.get_keep_qstates()<<endl;
    LoadRotationMatrix (newSystem.get_sites(), ketrotateMatrix, baseState);
    double error;
    if (!mpigetrank())
      error = makeRotateMatrix(bratracedMatrix, brarotateMatrix, sweepParams.get_keep_states(), sweepParams.get_keep_qstates());
      //error = makeRotateMatrix(bratracedMatrix, brarotateMatrix, max(100,sweepParams.get_keep_states()/4), sweepParams.get_keep_qstates(),0.0);
    #ifndef SERIAL
    broadcast(world, ketrotateMatrix, 0);
    broadcast(world, brarotateMatrix, 0);
    #endif

    SaveRotationMatrix (newSystem.get_sites(), brarotateMatrix, pb.wavenumber());
    bratracedMatrix.deallocate();



    //cout <<"before renormalize: "<<endl<<newSystem.get_braStateInfo()<<endl;
    newSystem.transform_operators(brarotateMatrix,ketrotateMatrix);
    //cout <<"after renormalize: "<<endl<<newSystem.get_braStateInfo()<<endl;

    {
      long memoryToFree = newSystem.getdata() - system.getdata();
      long newsysmem = newSystem.memoryUsed();
      newSystem.moveToNewMemory(system.getdata());
      Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
      //system.clear();
    }

    StackSpinBlock::store (forward, newSystem.get_sites(), newSystem, pb.wavenumber(), baseState); // if restart, just restoring an existing block --
	  newSystem.printOperatorSummary();
    system=newSystem;
  }
  //TODO
  //It seems that there is no need to do Last Step of Sweep.
}

