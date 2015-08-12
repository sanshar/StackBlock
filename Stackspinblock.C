/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "Stackspinblock.h"
#include "IntegralMatrix.h"
#include "stackopxop.h"
#include "operatorfunctions.h"
#include "Stackwavefunction.h"
#include <boost/format.hpp>
#include "distribute.h"
#include "StackOperators.h"
#include "Stackdensity.h"
#include "csf.h"
#include "StateInfo.h"
#include <boost/serialization/array.hpp>
#include <boost/function.hpp>
#include <boost/functional.hpp>
#include <boost/bind.hpp>
#include <boostutils.h>


#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted{
using namespace operatorfunctions;

void StackSpinBlock::deallocate() 
{
  Stackmem[omprank].deallocate(data, totalMemory);
}

//the relevant data is in the brackets
//------- [---------------]       this is old data
//[---------]-------------        this is new data
//the user is incharge of making sure that there is enough space in the new data
//also some part of the new data might overlap with the old data (see the figure above)
// i.e. end of new data might point to the same memory location as the start of old data.
//so copy from begining and move to the end
void StackSpinBlock::moveToNewMemory(double* pData)
{
  double* oldData = data;
  data = pData;

  DCOPY(totalMemory, oldData, 1, data, 1);

  double* localdata = data; 
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
  {
    if(it->second->is_core()) {
      
      for (int i=0; i<it->second->get_size(); i++) {
	int vecsize = it->second->get_local_element(i).size();
	for (int j=0; j<vecsize; j++) {
	  it->second->get_local_element(i)[j]->set_data(localdata);
	  it->second->get_local_element(i)[j]->allocateOperatorMatrix();
	  localdata = localdata + it->second->get_local_element(i)[j]->memoryUsed();
	}
      }
    }
  }

}

void StackSpinBlock::printOperatorSummary()
{
#ifndef SERIAL
  mpi::communicator world;

  if (mpigetrank() != 0) {
    for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
      sendobject(it->second->get_size(), 0);
  }
  else {
    for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core())
      {
         p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  ";
      }
      else
      {
         p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  ";      
      }
      
      vector<int> numops(world.size(), 0);
      for (int proc = 0; proc <world.size(); proc++) {
         if (proc != 0) 
            receiveobject(numops[proc],proc);
         else 
            numops[proc] = it->second->get_size();
         p2out << " " << numops[proc]<<"  ";
      }
      p2out << endl;
    }
  }
#else
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::const_iterator it = ops.begin(); it != ops.end(); ++it)
  {
    if(it->second->is_core()) 
      p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  ";      
    else
      p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  ";      
    p2out << endl;
  }
#endif
  
}
ostream& operator<< (ostream& os, const StackSpinBlock& b)
{
  os << "\t\t\t Sites ::  ";
  for (int i = 0; i < b.sites.size(); ++i) { os << b.sites[i] << " "; } 
  
  if (dmrginp.outputlevel() > 1) {
    os << endl;
    os << b.braStateInfo;
    os << b.ketStateInfo;
  }
  else {
    os <<"    # states: "<<b.braStateInfo.totalStates;
    os <<"    # states: "<<b.ketStateInfo.totalStates<<endl;
  }
  return os;
}


StackSpinBlock::StackSpinBlock () : 
  additionalMemory(0),
  additionaldata(0),
  totalMemory(0),
  data(0),
  localstorage(false),
  name (rand()), 
  integralIndex(0),
  loopblock(false),
  direct(false), complementary(false), normal(true), leftBlock(0), rightBlock(0) { }

StackSpinBlock::StackSpinBlock(int start, int finish, int p_integralIndex, bool implicitTranspose, bool is_complement) :  
  name (rand()), 
  integralIndex(p_integralIndex),
  direct(false), leftBlock(0), rightBlock(0),  additionalMemory(0),
  additionaldata(0)
{
  complementary = is_complement;
  normal = !is_complement;

  //this is used to make dot block and we make the 
  //additional operators by default because they are cheap
  default_op_components(is_complement, implicitTranspose);

  std::vector<int> sites; 
  if (dmrginp.use_partial_two_integrals()) {
    if (start != finish) {
      pout << "Cannot use partial two electron integrals, when making spin block with more than two orbitals"<<endl;
      abort();
    }
    std::vector<int> o;
    for (int i=dmrginp.spatial_to_spin()[start]; i<dmrginp.spatial_to_spin()[start+1]; i+=2)
      o.push_back(i/2);
    twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
    twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
    mpi::communicator world;
    PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
    mpi::broadcast(world, ar, 0);
#endif
  }
  else
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  int lower = min(start, finish);
  int higher = max(start, finish);
  sites.resize(higher - lower + 1);
  for (int i=0; i < sites.size(); i++)
      sites[i] = lower + i;

  BuildTensorProductBlock(sites);
}

StackSpinBlock::StackSpinBlock (const StackSpinBlock& b) { *this = b; }

StackSpinBlock::StackSpinBlock(const StateInfo& s, int pintegralIndex)
{
  additionalMemory=0;
  additionaldata=0;
  data = 0;
  totalMemory = 0;
  braStateInfo = s;
  ketStateInfo = s;
  sites.resize(0);
  integralIndex = pintegralIndex;
}

void StackSpinBlock::BuildTensorProductBlock(std::vector<int>& new_sites)
{

  if (twoInt.get() == 0 && dmrginp.use_partial_two_integrals()) { //this is when dummy block is being added for non zero spin
    std::vector<int> o;
    for (int i=dmrginp.spatial_to_spin()[new_sites[0]]; i<dmrginp.spatial_to_spin()[new_sites[new_sites.size()-1]+1]; i+=2)
      o.push_back(i/2);
    twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
    twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
    mpi::communicator world;
    PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
    mpi::broadcast(world, ar, 0);
    //world.broadcast(twoInt);
#endif
  }
  else if (twoInt.get() == 0)
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  name = get_name();


  std::vector< std::vector<Csf> > ladders;
  std::vector< Csf > dets; 

  if (dmrginp.spinAdapted()) {
    sites = new_sites;
    dets = CSFUTIL::spinfockstrings(new_sites, ladders);
  }
  else {
    for (int i=0; i<new_sites.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[new_sites[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[new_sites[i]]+1 );
    }
    // dets is   {0, alpha, beta, alphabeta;}
    // in StateInfo.C, they are sorted, they become {0,beta,alpha,alphabeta}
    dets = CSFUTIL::spinfockstrings(new_sites);
    for (int j=0; j<dets.size(); j++)
      ladders.push_back(std::vector<Csf>(1,dets[j]));
  }

  braStateInfo = StateInfo(dets);
  ketStateInfo = StateInfo(dets);

  setstoragetype(LOCAL_STORAGE);
  complementary_sites = make_complement(sites);

  //this is where we are building blocks from sites
  //currently only used for building single site blocks
  //temporarily disable screening for single site blocks
  double twoindex_ScreenTol = dmrginp.twoindex_screen_tol();
  double oneindex_ScreenTol = dmrginp.oneindex_screen_tol();
  if (new_sites.size() == 1 ){
    dmrginp.twoindex_screen_tol() = 0.0;
    dmrginp.oneindex_screen_tol() = 0.0;
  }

  totalMemory = build_iterators();

  if (totalMemory != 0)
    data = Stackmem[omprank].allocate(totalMemory);
  if (new_sites.size() == 1 ) {
    dmrginp.twoindex_screen_tol() = twoindex_ScreenTol;
    dmrginp.oneindex_screen_tol() = oneindex_ScreenTol;
  }

  build_operators(dets, ladders);

}

std::vector<int> StackSpinBlock::make_complement(const std::vector<int>& sites)
{
  std::vector<int> complementary_sites;
  for (int i=0; i<dmrginp.last_site(); ++i)
    if (find(sites.begin(), sites.end(), i) == sites.end())
      complementary_sites.push_back(i);
  
  return complementary_sites;
  
}

long StackSpinBlock::build_iterators()
{
  dmrginp.builditeratorsT->start();
  long memoryRequired = 0;
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
  {
    if(it->second->is_core()) 
      memoryRequired += it->second->build_iterators(*this, true);
    else
      it->second->build_iterators(*this, false);
  }
  dmrginp.builditeratorsT->stop();
  return memoryRequired;
}


void StackSpinBlock::build_operators(std::vector< Csf >& dets, std::vector< std::vector<Csf> >& ladders)
{
  dmrginp.buildcsfops->start();
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  
  double* localdata = data; 
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(it->second->is_core()) {
      localdata = it->second->allocateOperators(braStateInfo, ketStateInfo, localdata);
      for (int i=0; i<it->second->get_size(); i++)
	for (int j=0; j<it->second->get_local_element(i).size(); j++) 
	  allops.push_back(it->second->get_local_element(i)[j]);
    }
  }
  
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<allops.size(); i++)
    allops[i]->buildUsingCsf(*this, ladders, dets);
  
  dmrginp.buildcsfops->stop();
}
  


void StackSpinBlock::build_operators()
{
  double* localdata = data; 
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    {
      if(it->second->is_core()) {
	localdata = it->second->allocateOperators(braStateInfo, ketStateInfo, localdata);
	it->second->build_operators(*this);
      }
    }
}



void StackSpinBlock::build_and_renormalise_operators(const std::vector<Matrix>& rotateMatrix, const StateInfo *newStateInfo)
{
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core()) {
      for (int i=0; i<it->second->get_size(); i++)
	for (int j=0; j<it->second->get_local_element(i).size(); j++)
	  allops.push_back(it->second->get_local_element(i)[j]);
    }
  }

  SplitStackmem();
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<allops.size(); i++)
    allops[i]->build_and_renormalise_transform(this, rotateMatrix, newStateInfo);
  MergeStackmem();
}


void StackSpinBlock::build_and_renormalise_operators(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core()) {
      it->second->build_and_renormalise_operators(*this, ot, leftMat, bra, rightMat, ket);
    }
  }
}

void StackSpinBlock::CleanUpOperators()
{
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
  {
      
    for (int i=0; i<it->second->get_size(); i++) {
      int vecsize = it->second->get_local_element(i).size();
      for (int j=0; j<vecsize; j++) {
	it->second->get_local_element(i)[j]->CleanUp();
	it->second->get_local_element(i)[j]->set_totalMemory() = 0;
      }
    }
  }

}


void StackSpinBlock::transform_operators(std::vector<Matrix>& rotateMatrix) 
{
  p1out << "\t\t\t Transforming to new basis " << endl;
  Timer transformtimer;

  StateInfo oldStateInfo = braStateInfo;
  std::vector<SpinQuantum> newQuanta;
  std::vector<int> newQuantaStates;
  std::vector<int> newQuantaMap;
  for (int Q = 0; Q < rotateMatrix.size (); ++Q)
  {
    if (rotateMatrix [Q].Ncols () != 0)
      {
	newQuanta.push_back (braStateInfo.quanta [Q]);
	newQuantaStates.push_back (rotateMatrix [Q].Ncols ());
	newQuantaMap.push_back (Q);
      }
  }
  StateInfo newStateInfo = StateInfo (newQuanta, newQuantaStates, newQuantaMap);

  //first find out how much memory is required and allocate all the operators
  long requiredMemory = 0;
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    requiredMemory += it->second->getRequiredMemory(newStateInfo, newStateInfo);
  totalMemory = requiredMemory;
  data = Stackmem[omprank].allocate(requiredMemory);
  double* localdata = data;
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    localdata = it->second->allocateOperators(newStateInfo, newStateInfo, localdata);

  build_and_renormalise_operators( rotateMatrix, &newStateInfo );

  braStateInfo = newStateInfo;
  braStateInfo.AllocatePreviousStateInfo ();
  *braStateInfo.previousStateInfo = oldStateInfo;
  ketStateInfo = braStateInfo;


  p3out << "\t\t\t total elapsed time " << globaltimer.totalwalltime() << " ... " 
       << globaltimer.elapsedwalltime() << endl;

  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (! it->second->is_core())
      ops[it->first]->set_core(true);

  this->direct = false;
  p3out << "\t\t\t transform time " << transformtimer.elapsedwalltime() << endl;

  leftBlock = 0;
  rightBlock = 0;
}


void StackSpinBlock::transform_operators(std::vector<Matrix>& leftrotateMatrix, std::vector<Matrix>& rightrotateMatrix, bool clearRightBlock, bool clearLeftBlock) 
{
  p1out << "\t\t\t Transforming to new basis " << endl;
  Timer transformtimer;

  StateInfo oldbraStateInfo=braStateInfo, oldketStateInfo=ketStateInfo;
  StateInfo newbraStateInfo, newketStateInfo;
  StateInfo::transform_state(leftrotateMatrix, braStateInfo, newbraStateInfo);
  StateInfo::transform_state(rightrotateMatrix, ketStateInfo, newketStateInfo);

  //first find out how much memory is required and allocate all the operators
  long requiredMemory = 0;
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    requiredMemory += it->second->getRequiredMemory(newbraStateInfo, newketStateInfo);
  totalMemory = requiredMemory;
  data = Stackmem[omprank].allocate(requiredMemory);
  double* localdata = data;
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    localdata = it->second->allocateOperators(newbraStateInfo, newketStateInfo, localdata);

  build_and_renormalise_operators( leftrotateMatrix, &newbraStateInfo, rightrotateMatrix, &newketStateInfo );

  braStateInfo = newbraStateInfo;
  braStateInfo.AllocatePreviousStateInfo ();
  *braStateInfo.previousStateInfo = oldbraStateInfo;

  ketStateInfo = newketStateInfo;
  ketStateInfo.AllocatePreviousStateInfo ();
  *ketStateInfo.previousStateInfo = oldketStateInfo;


  p3out << "\t\t\t total elapsed time " << globaltimer.totalwalltime() << " ... " 
       << globaltimer.elapsedwalltime() << endl;


  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (! it->second->is_core())
      ops[it->first]->set_core(true);

  this->direct = false;
  p3out << "\t\t\t transform time " << transformtimer.elapsedwalltime() << endl;

  if (leftBlock && clearLeftBlock)
    leftBlock->clear();
  if (rightBlock && clearRightBlock)
    rightBlock->clear();
  leftBlock = 0;
  rightBlock = 0;
}


void StackSpinBlock::clear()
{
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    it->second->clear();
}


void StackSpinBlock::recreateStateInfo(int condition)
{
  if (rightBlock == 0) return;
  ketStateInfo = StateInfo();
  braStateInfo = StateInfo();
  if(dmrginp.transition_diff_irrep()){
    if( condition== PARTICLE_SPIN_NUMBER_CONSTRAINT)
    // When bra and ket wavefuntion have different spatial or spin irrep,
    // Bra Stateinfo for the big block should not be used with quantum number of effective_molecule_quantum
      TensorProduct (leftBlock->braStateInfo, rightBlock->braStateInfo, dmrginp.bra_quantum(), EqualQ, braStateInfo);
    else if (condition== NO_PARTICLE_SPIN_NUMBER_CONSTRAINT) 
      TensorProduct (leftBlock->braStateInfo, rightBlock->braStateInfo, dmrginp.bra_quantum(), LessThanQ, braStateInfo);
    // When bra and ket wavefuntion have different spatial or spin irrep,
  }
 else {
   TensorProduct (leftBlock->braStateInfo, rightBlock->braStateInfo, braStateInfo, condition);
 }

  TensorProduct (leftBlock->ketStateInfo, rightBlock->ketStateInfo, ketStateInfo, condition);
}

void StackSpinBlock::collectQuanta()
{
  braStateInfo.CollectQuanta();
  ketStateInfo.CollectQuanta();
}

void StackSpinBlock::BuildSumBlockSkeleton(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, bool collectQuanta, StateInfo* compState)
{

  name = get_name();
  p1out << "\t\t\t Building Sum Block " << name << endl;
  leftBlock = &lBlock;
  rightBlock = &rBlock;

  sites.reserve (lBlock.sites.size () + rBlock.sites.size ());

  dmrginp.blockintegrals -> start();
  
  if (dmrginp.use_partial_two_integrals()) {
    if (rBlock.sites.size() == 1) {
      std::vector<int> o;
      for (int i=dmrginp.spatial_to_spin().at(rBlock.sites[0]); i<dmrginp.spatial_to_spin().at(rBlock.sites[0]+1); i+=2)
	o.push_back(i/2);
      twoInt = boost::shared_ptr<PartialTwoElectronArray> (new PartialTwoElectronArray(o));
      twoInt->Load(dmrginp.load_prefix(), integralIndex);
#ifndef SERIAL
      mpi::communicator world;
      PartialTwoElectronArray& ar = dynamic_cast<PartialTwoElectronArray&>(*twoInt.get());
      mpi::broadcast(world, ar, 0);
#endif

    }
    //pout << "Cannot use partial two electron integrals, when the dot block has more than one orbital"<<endl;
    //abort();

  }
  else 
    twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex],  boostutils::null_deleter());

  dmrginp.blockintegrals -> stop();

  dmrginp.blocksites -> start();

  sites = lBlock.sites;
  copy (rBlock.sites.begin(), rBlock.sites.end (), back_inserter (sites));
  sort(sites.begin(), sites.end());
  complementary_sites = make_complement(sites);
  p2out << "\t\t\t ";
  for (int i = 0; i < sites.size(); ++i) p2out << sites[i] << " ";
  p2out << endl;
  dmrginp.blocksites -> stop();

  dmrginp.statetensorproduct -> start();
  if(dmrginp.transition_diff_irrep()){
    if( condition== PARTICLE_SPIN_NUMBER_CONSTRAINT)
    // When bra and ket wavefuntion have different spatial or spin irrep,
    // Bra Stateinfo for the big block should not be used with quantum number of effective_molecule_quantum
      TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, dmrginp.bra_quantum(), EqualQ, braStateInfo);
    else if (condition== NO_PARTICLE_SPIN_NUMBER_CONSTRAINT) 
      TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, dmrginp.bra_quantum(), LessThanQ, braStateInfo,compState);
    // When bra and ket wavefuntion have different spatial or spin irrep,
  }
 else {
   TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, braStateInfo, condition, compState);
 }

  TensorProduct (lBlock.ketStateInfo, rBlock.ketStateInfo, ketStateInfo, condition, compState);
  dmrginp.statetensorproduct -> stop();

  dmrginp.statecollectquanta -> start();
  if (collectQuanta) {
  if (!( (dmrginp.hamiltonian() == BCS && condition == SPIN_NUMBER_CONSTRAINT)  ||
	 (dmrginp.hamiltonian() != BCS && condition == PARTICLE_SPIN_NUMBER_CONSTRAINT))) {
    braStateInfo.CollectQuanta();
    ketStateInfo.CollectQuanta();
  }
  }
  else {
    braStateInfo.hasCollectedQuanta = false;
    ketStateInfo.hasCollectedQuanta = false;
  }
  dmrginp.statecollectquanta -> stop();

}

void StackSpinBlock::BuildSumBlock(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, bool collectQuanta, StateInfo* compState)
{
  if (!(lBlock.integralIndex == rBlock.integralIndex && lBlock.integralIndex == integralIndex))  {
    pout << "The left, right and dot block should use the same integral indices"<<endl;
    pout << "ABORTING!!"<<endl;
    exit(0);
  }
  dmrginp.buildsumblock -> start();
  BuildSumBlockSkeleton(condition, lBlock, rBlock, collectQuanta, compState);

  totalMemory = build_iterators();
  if (totalMemory != 0)
    data = Stackmem[omprank].allocate(totalMemory);

  dmrginp.buildblockops -> start();
  build_operators();
  dmrginp.buildblockops -> stop();
  dmrginp.buildsumblock -> stop();
}


void StackSpinBlock::operator= (const StackSpinBlock& b)
{
  localstorage = b.localstorage;
  name = b.name;
  complementary = b.is_complementary();
  normal = b.is_normal();
  loopblock = b.is_loopblock();

  sites = b.sites;
  complementary_sites = b.complementary_sites;
  integralIndex = b.integralIndex;

  direct = b.is_direct();

  braStateInfo = b.braStateInfo;
  ketStateInfo = b.ketStateInfo;
  leftBlock = b.leftBlock;
  rightBlock = b.rightBlock;
  twoInt = b.twoInt;
  ops = b.ops;
  totalMemory = b.totalMemory;
  data = b.data;
  additionalMemory = b.additionalMemory;
  additionaldata = b.additionaldata;
}

void StackSpinBlock::initialise_op_array(opTypes optype, bool is_core)
{
  ops[optype] = make_new_stackop(optype, is_core);
  return;
}

void StackSpinBlock::multiplyOverlap(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{
  boost::shared_ptr<StackSparseMatrix> op = leftBlock->get_op_array(OVERLAP).get_local_element(0)[0];
  bool deallocate1 = op->memoryUsed() == 0 ? true : false; 
  op->allocate(leftBlock->get_braStateInfo(), leftBlock->get_ketStateInfo());
  op->build(*leftBlock);
  
  boost::shared_ptr<StackSparseMatrix> overlap = rightBlock->get_op_array(OVERLAP).get_local_element(0)[0];
  bool deallocate2 = overlap->memoryUsed() == 0 ? true : false; 
  overlap->allocate(rightBlock->get_braStateInfo(), rightBlock->get_ketStateInfo());
  overlap->build(*rightBlock);
  
  TensorMultiply(leftBlock, *op, *overlap, this, c, v, op->get_deltaQuantum(0) ,1.0);  // dmrginp.ef
  
  if (deallocate2) overlap->deallocate();
  if (deallocate1) op->deallocate();
  
}

int procWithMinOps(std::vector<boost::shared_ptr<StackSparseMatrix> >& allops)
{
  int size = 1;
#ifndef SERIAL
  boost::mpi::communicator world;
  size = world.size();
#endif
  std::vector<int> numOps(size, 0);

  numOps[mpigetrank()] = allops.size();
  int minproc = 0;
#ifndef SERIAL
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, &numOps[0], size, MPI_INT, MPI_SUM);
#endif

  for (int i=0; i< size; i++) {
    if (numOps[i] < numOps[minproc])
      minproc = i;
    pout << numOps[i]<<"  ";
  }
  pout << endl<<"  minproc "<< minproc<<endl;
  return minproc;
}

void StackSpinBlock::multiplyH(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{

  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));

  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<FUNCTOR2> allfuncs;

  //accumulate ham
  StackWavefunction* v_array; 
  initiateMultiThread(v, v_array, numthrds);

  FUNCTOR2 f4 = boost::bind(&stackopxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() ); 
  FUNCTOR2 f5 = boost::bind(&stackopxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() ); 
  FUNCTOR2 f6 = boost::bind(&stackopxop::cdxcdcomp, otherBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() );
  FUNCTOR2 f7 = boost::bind(&stackopxop::ddxcccomp, otherBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() );


  for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
    for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
      allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
      allfuncs.push_back(f4);
    }

  for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
    for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
      allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
      allfuncs.push_back(f5);
    }



  if (dmrginp.hamiltonian() != HUBBARD) {
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=0; j<loopBlock->get_op_array(CRE_DES).get_local_element(i).size(); j++) {
	allops.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	allfuncs.push_back(f6);
      }

    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++)
      for (int j=0; j<loopBlock->get_op_array(CRE_CRE).get_local_element(i).size(); j++) {
	allops.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	allfuncs.push_back(f7);
      }
  }

  int proc = procWithMinOps(allops);
  FUNCTOR2 f1 = boost::bind(&stackopxop::hamandoverlap, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum(), coreEnergy[integralIndex], proc);
  if (proc == mpigetrank()) {
    allops.push_back(rightBlock->get_op_rep(OVERLAP, hq)); allfuncs.push_back(f1);//this is just a placeholder function  
  }

  SplitStackmem();
  dmrginp.tensormultiply->start();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i<allops.size(); i++)  {
    allfuncs[i](allops[i]);
  }
  dmrginp.tensormultiply->stop();  

  MergeStackmem();
  accumulateMultiThread(v, v_array, numthrds);
  distributedaccumulate(*v);
}





void StackSpinBlock::diagonalH(DiagonalMatrix& e) const
{
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;
  
  DiagonalMatrix* e_array = new DiagonalMatrix[numthrds];
  for (int i=0; i<numthrds; i++)
    e_array[i] = e;

  if (mpigetrank() == 0) {
    for (int i=0; i<e.Nrows(); i++)
      e(i+1) += coreEnergy[integralIndex];
  }

  FUNCTOR2 f3 = boost::bind(&stackopxop::cdxcdcomp_d, otherBlock, _1, this, e_array);

  std::vector<boost::shared_ptr<StackSparseMatrix> > allops;
  std::vector<FUNCTOR2> allfuncs;


  if (dmrginp.hamiltonian() != HUBBARD) {
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=0; j<loopBlock->get_op_array(CRE_DES).get_local_element(i).size(); j++) {
	allops.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	allfuncs.push_back(f3);
      }
  }

  int proc = procWithMinOps(allops);

  FUNCTOR2 f1 = boost::bind(&stackopxop::ham_d, loopBlock, _1, this, e_array, proc);
  FUNCTOR2 f2 = boost::bind(&stackopxop::ham_d, otherBlock, _1, this, e_array, proc);

  if (proc == mpigetrank()) {
    allops.push_back(loopBlock->get_op_rep(HAM, hq)); allfuncs.push_back(f1); //the function is placeholder
    allops.push_back(otherBlock->get_op_rep(HAM, hq)); allfuncs.push_back(f2);//the function is placeholder
  }


  SplitStackmem();
  dmrginp.tensormultiply->start();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i<allops.size(); i++)  {
    allfuncs[i](allops[i]);
  }
  dmrginp.tensormultiply->stop();  
  MergeStackmem();

  for (int i=0; i<numthrds; i++)
    e += e_array[i];
  delete [] e_array;
  distributedaccumulate(e);

}

void StackSpinBlock::BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, bool random, const bool haveNormops)
{
  name = get_name();

  if (dmrginp.spinAdapted()) {
    sites = sts;
  }
  else {
    for (int i=0; i<sts.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]+1 );
    }
  }

  complementary_sites = make_complement(sites);

  assert (sites.size () > 0);
  sort (sites.begin (), sites.end ());

  //always have implicit transpose in this case
  default_op_components(!haveNormops, true);
  
  setstoragetype(DISTRIBUTED_STORAGE);


  std::vector< Csf > dets;
  std::vector< Csf > det_ex;
  Timer slatertimer;

  for (int i = 0; i < qnumbers.size (); ++i)
    {
      if (distribution [i] == 0 || (qnumbers [i].get_n() > 2*sites.size ())) continue;

      if(dmrginp.spinAdapted()) 
	det_ex = Csf::distribute (qnumbers [i].get_n(), qnumbers [i].get_s().getirrep(), IrrepVector(qnumbers [i].get_symm().getirrep(), 0) , sites [0],
				  sites [0] + sites.size (), dmrginp.last_site(), integralIndex);
      else
	det_ex = Csf::distributeNonSpinAdapted (qnumbers [i].get_n(), qnumbers [i].get_s().getirrep(), IrrepVector(qnumbers [i].get_symm().getirrep(), 0) , sites [0],
						sites [0] + sites.size (), dmrginp.last_site(), integralIndex);

      multimap <double, Csf > slater_emap;

      for (int j = 0; j < det_ex.size(); ++j) {
	slater_emap.insert (pair <double, Csf > (csf_energy (det_ex[j], integralIndex), det_ex[j]));
      }

      multimap <double, Csf >::iterator m = slater_emap.begin();
      int sz = det_ex.size();
      det_ex.resize (min (distribution [i], sz));
      for (int j = 0; j < det_ex.size(); ++j)
        {
          det_ex[j] = m->second;
          ++m;
        }

      copy (det_ex.begin(), det_ex.end(), back_inserter (dets));
    }

  std::multimap<SpinQuantum, Csf > tmp;
  for (int i = 0; i < dets.size (); ++i) {
    tmp.insert(pair<SpinQuantum, Csf > (SpinQuantum ( (dets [i]).n, (dets [i]).S, (dets[i]).sym_is()), dets[i]));
  }
  std::multimap<SpinQuantum, Csf >::iterator tmpiter = tmp.begin();
  dets.clear();
  for (; tmpiter != tmp.end();)
  {
    dets.push_back(tmpiter->second);
    tmpiter++;
  }

  braStateInfo = StateInfo (dets);
  ketStateInfo = StateInfo (dets);

  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  totalMemory = build_iterators();
  if (totalMemory != 0)
    data = Stackmem[omprank].allocate(totalMemory);
  pout << "Allocating "<<totalMemory<<" for the block "<<endl;

  p3out << "\t\t\t time in slater distribution " << slatertimer.elapsedwalltime() << endl;

  std::vector< std::vector<Csf> > ladders; ladders.resize(dets.size());

#pragma omp parallel for schedule(dynamic)
  for (int i=0; i< dets.size(); i++)
    ladders[i] = dets[i].spinLadder(min(2, dets[i].S.getirrep()));


  build_operators(dets, ladders);
  p3out << "\t\t\t time in slater operator build " << slatertimer.elapsedwalltime() << endl;


}

void StackSpinBlock::BuildSingleSlaterBlock(std::vector<int> sts) {
  name = get_name();
  if (dmrginp.spinAdapted()) {
    sites = sts;
  }
  else {
    for (int i=0; i<sts.size(); i++) {
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]   );
      sites.push_back( dmrginp.spatial_to_spin()[sts[i]]+1 );
    }
  }

  complementary_sites = make_complement(sites);
  assert (sites.size () > 0);
  sort (sites.begin (), sites.end ());

  int left = sites[0], right = sites[0] + sites.size(), edge = dmrginp.last_site();
  int n = 0, sp = 0;
  std::vector<bool> tmp(0);
  IrrepSpace irrep(0);

  if (dmrginp.spinAdapted()) {
    for (int orbI = left; orbI < right; ++orbI) {
      n += dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]] + dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1];
      sp += dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]] - dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1];

      // FIXME: NN wrote, follows don't work correctly for non-abelian symmetry
      if (dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]] == 1) {
        irrep = IrrepSpace(Symmetry::add(irrep.getirrep(),SymmetryOfSpatialOrb(orbI).getirrep())[0]);
      }
      if (dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1] == 1) {
        irrep = IrrepSpace(Symmetry::add(irrep.getirrep(),SymmetryOfSpatialOrb(orbI).getirrep())[0]);
      }
    }

    for (int i = 0; i < dmrginp.spatial_to_spin()[left]; ++i) {
      tmp.push_back(0);
    }
    for (int orbI = left; orbI < right; ++orbI) {
      tmp.push_back(dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]]);
      tmp.push_back(dmrginp.hf_occupancy()[dmrginp.spatial_to_spin()[orbI]+1]);
    }
    for (int i = 0; i < dmrginp.spatial_to_spin()[edge]-dmrginp.spatial_to_spin()[right]; ++i) {
      tmp.push_back(0);
    }
  } else {
    for (int i = 0; i < left; ++i) {
      tmp.push_back(0);
    }
    for (int orbI = left; orbI < right; ++orbI) {
      if (dmrginp.hf_occupancy()[orbI] == 1) {
        n += 1;
        sp += SpinOf(orbI);
      }
      tmp.push_back(dmrginp.hf_occupancy()[orbI]);
    }
    for (int i = 0; i < edge-right; ++i) {
      tmp.push_back(0);
    }
  }

  Slater new_det = Slater(Orbstring(tmp));
  map<Slater, double> m;
  m[new_det] = 1.0;

  Csf origin(m, n, SpinSpace(sp), sp, IrrepVector(irrep.getirrep(), 0));
  std::vector<Csf> dets(1, origin);
  braStateInfo = StateInfo (dets);
  ketStateInfo = StateInfo (dets);
  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[integralIndex], boostutils::null_deleter());

  totalMemory = build_iterators();
  if (totalMemory != 0)
    data = Stackmem[omprank].allocate(totalMemory);
  pout << "Allocating "<<totalMemory<<" for the block "<<endl;

  std::vector< std::vector<Csf> > ladders; ladders.resize(dets.size());
  for (int i=0; i< dets.size(); i++)
    ladders[i] = dets[i].spinLadder(min(2, dets[i].S.getirrep()));

  build_operators(dets, ladders);

}


void initialiseSingleSiteBlocks(std::vector<StackSpinBlock>& singleSiteBlocks, int integralIndex) {
  
  int niter; 
  if (dmrginp.spinAdapted())
    niter = dmrginp.last_site();
  else
    niter = dmrginp.last_site()/2; 

  singleSiteBlocks.resize(niter);
  int dotStart, dotEnd;
  for (int i=0; i<niter; i++) {
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));
    singleSiteBlocks[i] = StackSpinBlock(i, i, integralIndex, true);
  }
}


void StackSpinBlock::sendcompOps(StackOp_component_base& opcomp, int I, int J, int optype, int compsite)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    world.send(processorindex(compsite), optype+i*10+1000*J+100000*I, *oparray[i]);

    //now broadcast the data
    //MPI::COMM_WORLD.Bcast(oparray[i]->get_data(), oparray[i]->memoryUsed(), MPI_DOUBLE, trimap_2d(I, J, dmrginp.last_site()));
    MPI::COMM_WORLD.Send(oparray[i]->get_data(), oparray[i]->memoryUsed(), MPI_DOUBLE, processorindex(compsite), optype+i*10+1000*J+100000*I);
  }  
#endif
}

void StackSpinBlock::recvcompOps(StackOp_component_base& opcomp, int I, int J, int optype)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    world.recv(processorindex(trimap_2d(I, J, dmrginp.last_site())), optype+i*10+1000*J+100000*I, *oparray[i]);

    double* data = Stackmem[omprank].allocate(oparray[i]->memoryUsed());
    if (additionalMemory == 0) 
      additionaldata=data; 
    additionalMemory+=oparray[i]->memoryUsed();
    
    oparray[i]->set_data(data);
    oparray[i]->allocateOperatorMatrix();

    //now broadcast the data
    MPI::COMM_WORLD.Recv(oparray[i]->get_data(), oparray[i]->memoryUsed(), MPI_DOUBLE, processorindex(trimap_2d(I, J, dmrginp.last_site())), optype+i*10+1000*J+100000*I);
  }
#endif
}


void StackSpinBlock::removeAdditionalOps() 
{
  Stackmem[omprank].deallocate(additionaldata, additionalMemory);
}


void StackSpinBlock::addAdditionalOps()
{
#ifndef SERIAL
  dmrginp.datatransfer->start();
  boost::mpi::communicator world;
  if (world.size() == 1)
    return; //there is no need to have additional compops

  int length = dmrginp.last_site();

  //distribute cre to all processors
  if (!ops[CRE]->is_local()) {
    for(int i=0; i<get_sites().size(); i++) {
      if (ops[CRE]->has(sites[i])) {
        if (processorindex(sites[i]) != mpigetrank()) {
	  ops[CRE]->add_local_indices(sites[i]);
	}

        ops[CRE]->set_local() = true;
	boost::shared_ptr<StackSparseMatrix> op = ops[CRE]->get_element(sites[i])[0];

	//this only broadcasts the frame but no data
        mpi::broadcast(world, *op, processorindex(sites[i]));

	//now allocate the data
	if (processorindex(sites[i]) != mpigetrank()) {

	  double *data = Stackmem[omprank].allocate(op->memoryUsed());
	  op->set_data(data);
	  if (additionalMemory == 0) 
	    additionaldata = data;
	  additionalMemory+=op->memoryUsed();
	  op->allocateOperatorMatrix();
	}

	//now broadcast the data
	MPI::COMM_WORLD.Bcast(op->get_data(), op->memoryUsed(), MPI_DOUBLE, processorindex(sites[i]));
      }
    }
  }

  //distribute DES to all processors
  if (has(DES) && !ops[DES]->is_local()) {
    for(int i=0; i<get_sites().size(); i++) {
      if (ops[DES]->has(sites[i])) {
	if (processorindex(sites[i]) != mpigetrank()) ops[DES]->add_local_indices(sites[i]);

	ops[DES]->set_local() = true;
	boost::shared_ptr<StackSparseMatrix> op = ops[DES]->get_element(sites[i])[0];

	//this only broadcasts the frame but no data
        mpi::broadcast(world, *op, processorindex(sites[i]));

	//now allocate the data when it is not already there
	if (processorindex(sites[i]) != mpigetrank()) {
	  double *data = Stackmem[omprank].allocate(op->memoryUsed());
	  op->set_data(data);
	  if (additionalMemory == 0) 
	    additionaldata=data; 
	  additionalMemory+=op->memoryUsed();
	  
	  op->allocateOperatorMatrix();
	}

	//now broadcast the data
	MPI::COMM_WORLD.Bcast(op->get_data(), op->memoryUsed(), MPI_DOUBLE, processorindex(sites[i]));
      }
    }
  }
  

  vector<int> dotindice;
  dotindice.push_back((sites[0] == 0) ? complementary_sites[0] : complementary_sites[complementary_sites.size()-1]);
  if (!dmrginp.spinAdapted()) { // when non-spinadapted, sites are spin orbitals
    dotindice.push_back((sites[0] == 0) ? complementary_sites[1] : complementary_sites[complementary_sites.size()-2]);    
  }
  for (int idx = 0; idx < dotindice.size(); ++idx) {
    int dotopindex = dotindice[idx];
    int I = dotopindex;

    //CCDcompI should be broadcast to all procs
    if (has(CRE_CRE_DESCOMP)) {
      if( !ops[CRE_CRE_DESCOMP]->is_local() ) {
      int fromproc = processorindex(I);
      if( ops[CRE_CRE_DESCOMP]->has(I) ) {
	if (fromproc != mpigetrank()) ops[CRE_CRE_DESCOMP]->add_local_indices(I);
	
	std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[CRE_CRE_DESCOMP]->get_element(I);
	for (int iproc =0; iproc<oparray.size(); iproc++) {
	  //this only broadcasts the frame but no data
	  mpi::broadcast(world, *oparray[iproc], fromproc);
	  
	  //now allocate the data when it is not already there
	  if (fromproc != mpigetrank()) {
	    double *data = Stackmem[omprank].allocate(oparray[iproc]->memoryUsed());
	    oparray[iproc]->set_data(data);
	    if (additionalMemory == 0) 
	      additionaldata=data; 
	    additionalMemory+=oparray[iproc]->memoryUsed();
	    
	    oparray[iproc]->allocateOperatorMatrix();
	  }
	  
	  //now broadcast the data
	  MPI::COMM_WORLD.Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc);
	}
      }
      }
    }
    //CCDcompI should be broadcast to all procs
    if (has(CRE_DES_DESCOMP)) { 
      if (!ops[CRE_DES_DESCOMP]->is_local() ) {
      int fromproc = processorindex(I);
      if( ops[CRE_DES_DESCOMP]->has(I) ) {
	if (fromproc != mpigetrank()) ops[CRE_DES_DESCOMP]->add_local_indices(I);
	
	std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[CRE_DES_DESCOMP]->get_element(I);
	for (int iproc =0; iproc<oparray.size(); iproc++) {
	  //this only broadcasts the frame but no data
	  mpi::broadcast(world, *oparray[iproc], fromproc);
	  
	  //now allocate the data when it is not already there
	  if (fromproc != mpigetrank()) {
	    double *data = Stackmem[omprank].allocate(oparray[iproc]->memoryUsed());
	    oparray[iproc]->set_data(data);
	    if (additionalMemory == 0) 
	      additionaldata=data; 
	    additionalMemory+=oparray[iproc]->memoryUsed();
	    
	    oparray[iproc]->allocateOperatorMatrix();
	  }
	  
	  //now broadcast the data
	  MPI::COMM_WORLD.Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc);
	}
      }
      }
    }
    //CDII should be broadcast to all procs
    if (has(CRE_DESCOMP)) {
    if (!ops[CRE_DESCOMP]->is_local() ) {
      int fromproc = processorindex(trimap_2d(I, I, length));
      if( ops[CRE_DESCOMP]->has(I,I) ) {
	if (fromproc != mpigetrank()) ops[CRE_DESCOMP]->add_local_indices(I,I);
	
	std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[CRE_DESCOMP]->get_element(I,I);
	for (int iproc =0; iproc<oparray.size(); iproc++) {
	  //this only broadcasts the frame but no data
	  mpi::broadcast(world, *oparray[iproc], fromproc);
	  
	  //now allocate the data when it is not already there
	  if (fromproc != mpigetrank()) {
	    double *data = Stackmem[omprank].allocate(oparray[iproc]->memoryUsed());
	    oparray[iproc]->set_data(data);
	    if (additionalMemory == 0) 
	      additionaldata=data; 
	    additionalMemory+=oparray[iproc]->memoryUsed();
	    
	    oparray[iproc]->allocateOperatorMatrix();
	  }
	  
	  //now broadcast the data
	  MPI::COMM_WORLD.Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc);
	}
      }

      if( ops[DES_DESCOMP]->has(I,I) ) {
	if (fromproc != mpigetrank()) ops[DES_DESCOMP]->add_local_indices(I,I);
	
	std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[DES_DESCOMP]->get_element(I,I);
	for (int iproc =0; iproc<oparray.size(); iproc++) {
	  //this only broadcasts the frame but no data
	  mpi::broadcast(world, *oparray[iproc], fromproc);
	  
	  //now allocate the data when it is not already there
	  if (fromproc != mpigetrank()) {
	    double *data = Stackmem[omprank].allocate(oparray[iproc]->memoryUsed());
	    oparray[iproc]->set_data(data);
	    if (additionalMemory == 0) 
	      additionaldata=data; 
	    additionalMemory+=oparray[iproc]->memoryUsed();
	    
	    oparray[iproc]->allocateOperatorMatrix();
	  }
	  
	  //now broadcast the data
	  MPI::COMM_WORLD.Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc);
	}
      }
      if (has(DES)) {
	if( ops[CRE_CRECOMP]->has(I,I) ) {
	  if (fromproc != mpigetrank()) ops[CRE_CRECOMP]->add_local_indices(I,I);
	  
	  std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[CRE_CRECOMP]->get_element(I,I);
	  for (int iproc =0; iproc<oparray.size(); iproc++) {
	    //this only broadcasts the frame but no data
	    mpi::broadcast(world, *oparray[iproc], fromproc);
	    
	    //now allocate the data when it is not already there
	    if (fromproc != mpigetrank()) {
	      double *data = Stackmem[omprank].allocate(oparray[iproc]->memoryUsed());
	      oparray[iproc]->set_data(data);
	      if (additionalMemory == 0) 
		additionaldata=data; 
	      additionalMemory+=oparray[iproc]->memoryUsed();
	      
	      oparray[iproc]->allocateOperatorMatrix();
	    }
	    
	    //now broadcast the data
	    MPI::COMM_WORLD.Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc);
	  }
	}
	if( ops[DES_CRECOMP]->has(I,I) ) {
	  if (fromproc != mpigetrank()) ops[DES_CRECOMP]->add_local_indices(I,I);
	  
	  std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[DES_CRECOMP]->get_element(I,I);
	  for (int iproc =0; iproc<oparray.size(); iproc++) {
	    //this only broadcasts the frame but no data
	    mpi::broadcast(world, *oparray[iproc], fromproc);
	    
	    //now allocate the data when it is not already there
	    if (fromproc != mpigetrank()) {
	      double *data = Stackmem[omprank].allocate(oparray[iproc]->memoryUsed());
	      oparray[iproc]->set_data(data);
	      if (additionalMemory == 0) 
		additionaldata=data; 
	      additionalMemory+=oparray[iproc]->memoryUsed();
	      
	      oparray[iproc]->allocateOperatorMatrix();
	    }
	    
	    //now broadcast the data
	    MPI::COMM_WORLD.Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc);
	  }
	}
      }
    }
    }

    for (int i=0; i<complementary_sites.size(); i++) {
      int compsite = complementary_sites[i];
      if (std::find(dotindice.begin(), dotindice.end(), compsite) != dotindice.end())
        continue;
      int I = (compsite > dotopindex) ? compsite : dotopindex;
      int J = (compsite > dotopindex) ? dotopindex : compsite;

      if (processorindex(compsite) == processorindex(trimap_2d(I, J, length)) || ops[CRE_DESCOMP]->is_local())
        continue;

      if (processorindex(compsite) == mpigetrank()) {
        //this will potentially receive some ops        
        bool other_proc_has_ops = true;
        world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
        if (other_proc_has_ops) {
	  ops[CRE_DESCOMP]->add_local_indices(I, J);
	  recvcompOps(*ops[CRE_DESCOMP], I, J, CRE_DESCOMP);
        }

        other_proc_has_ops = true;
        world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
        if (other_proc_has_ops) {
	  ops[DES_DESCOMP]->add_local_indices(I, J);
	  recvcompOps(*ops[DES_DESCOMP], I, J, DES_DESCOMP);
        }

	if (has(DES)) {
	  other_proc_has_ops = true;
	  world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	  if (other_proc_has_ops) {
	    ops[CRE_CRECOMP]->add_local_indices(I, J);
	    recvcompOps(*ops[CRE_CRECOMP], I, J, CRE_CRECOMP);
	  }
	  other_proc_has_ops = true;
	  world.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	  if (other_proc_has_ops) {
	    ops[DES_CRECOMP]->add_local_indices(I, J);
	    recvcompOps(*ops[DES_CRECOMP], I, J, DES_CRECOMP);
	  }
	}
      } 
      else {
        //this will potentially send some ops
        if (processorindex(trimap_2d(I, J, length)) == mpigetrank()) {
	  bool this_proc_has_ops = ops[CRE_DESCOMP]->has_local_index(I, J);
	  world.send(processorindex(compsite), 0, this_proc_has_ops);
	  if (this_proc_has_ops) {
	    sendcompOps(*ops[CRE_DESCOMP], I, J, CRE_DESCOMP, compsite);
	  }
          this_proc_has_ops = ops[DES_DESCOMP]->has_local_index(I, J);
	  world.send(processorindex(compsite), 0, this_proc_has_ops);
	  if (this_proc_has_ops) {
	    sendcompOps(*ops[DES_DESCOMP], I, J, DES_DESCOMP, compsite);     
	  }
	  if (has(DES)) {
	    this_proc_has_ops = ops[CRE_CRECOMP]->has_local_index(I, J);
	    world.send(processorindex(compsite), 0, this_proc_has_ops);
	    if (this_proc_has_ops) {
	      sendcompOps(*ops[CRE_CRECOMP], I, J, CRE_CRECOMP, compsite);     
	    }
	    this_proc_has_ops = ops[DES_CRECOMP]->has_local_index(I, J);
	    world.send(processorindex(compsite), 0, this_proc_has_ops);
	    if (this_proc_has_ops) {
	      sendcompOps(*ops[DES_CRECOMP], I, J, DES_CRECOMP, compsite);     
	    }
	  }
	  
        } 
	else 
	  continue;
      }
    }
  }
  dmrginp.datatransfer->stop();
#endif
}



}
