/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include <sys/time.h>
#include "time.h"
#include "screen.h"
#include "Stackspinblock.h"
#include "MatrixBLAS.h"
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
  for (long i=0; i<totalMemory; i++)
    data[i] = oldData[i];
  //DCOPY(totalMemory, oldData, 1, data, 1);
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

      vector<int> numops(calc.size(), 0);
      for (int proc = 0; proc <calc.size(); proc++) {
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
    if(it->second->is_core()) {
      p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Core Operators  ";      
    }
    else {
      p2out << "\t\t\t " << it->second->size()<<" :  "<<it->second->get_op_string()<<"  Virtual Operators  ";      
    }
    p2out << endl;
  }
#endif
  
}
ostream& operator<< (ostream& os, const StackSpinBlock& b)
{
  os << "\t\t\t Sites ::  "<<b.sites[0]<<"--"<<b.sites[b.sites.size()-1]<<endl;
  
  if (dmrginp.outputlevel() > 5) {
    os << endl;
    os << b.braStateInfo;
    os << b.ketStateInfo;
  }
  else {
    os <<"\t\t\t # states: "<<b.braStateInfo.totalStates;
    os <<"\t\t\t # states: "<<b.ketStateInfo.totalStates<<endl;
  }
  return os;
}

  //thus is used to build the edge block for responseaaav and responseaaac
StackSpinBlock StackSpinBlock::buildBigEdgeBlock(int start, int end, bool haveNorm, bool haveComp, int p_integralIndex, bool implicitTranspose)
{
    if (dmrginp.calc_type() == RESPONSEAAAV) {
      StackSpinBlock system(end-1,end-1, p_integralIndex, implicitTranspose);
      SpinQuantum moleculeQ = dmrginp.molecule_quantum();
      dmrginp.set_molecule_quantum() = SpinQuantum(1, SpinSpace(0), IrrepSpace(0)); 

      for (int i=end-2; i >= start; i--) {
	StackSpinBlock newSystem;
	pout << i <<"  ";
	StackSpinBlock site(i, i, p_integralIndex, implicitTranspose);
	system.addAdditionalOps();
	newSystem.default_op_components(false, haveNorm, haveComp, implicitTranspose);
	newSystem.set_integralIndex() = p_integralIndex;
	newSystem.setstoragetype(DISTRIBUTED_STORAGE);
	newSystem.BuildSumBlock (PARTICLE_NUMBER_CONSTRAINT, system, site);
	{
	  long memoryToFree = newSystem.getdata() - system.getdata();
	  long newsysmem = newSystem.memoryUsed();
	  newSystem.moveToNewMemory(system.getdata());
	  Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
	}
	system.clear();
	system = newSystem;
      }
      pout << endl;
      dmrginp.set_molecule_quantum() = moleculeQ; 

      
      return system;
      
    }
    else {
      StackSpinBlock system(end-1,end-1, p_integralIndex, implicitTranspose);
      SpinQuantum moleculeQ = dmrginp.molecule_quantum();
      dmrginp.set_molecule_quantum() = SpinQuantum(3, SpinSpace(0), IrrepSpace(0)); 

      for (int i=end-2; i >= start; i--) {
	StackSpinBlock newSystem;
	pout << i <<"  ";
	StackSpinBlock site(i, i, p_integralIndex, implicitTranspose);
	system.addAdditionalOps();
	newSystem.default_op_components(false, haveNorm, haveComp, implicitTranspose);
	newSystem.set_integralIndex() = p_integralIndex;
	newSystem.setstoragetype(DISTRIBUTED_STORAGE);
	newSystem.BuildSumBlock (HOLE_NUMBER_CONSTRAINT, system, site);
	{
	  long memoryToFree = newSystem.getdata() - system.getdata();
	  long newsysmem = newSystem.memoryUsed();
	  newSystem.moveToNewMemory(system.getdata());
	  Stackmem[omprank].deallocate(newSystem.getdata()+newsysmem, memoryToFree);
	}
	system.clear();
	system = newSystem;
	dmrginp.set_molecule_quantum().particleNumber+=2;
      }
      pout << endl;
      dmrginp.set_molecule_quantum() = moleculeQ; 

      return system;      
    }

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
    mpi::broadcast(calc, ar, 0);
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

StackSpinBlock::StackSpinBlock(int start, int finish, int p_integralIndex, std::vector<int> non_active_orbs, bool implicitTranspose, bool is_complement) :  
  name (rand()), 
  integralIndex(p_integralIndex),
  nonactive_orbs(non_active_orbs),
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
    mpi::broadcast(calc, ar, 0);
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

void StackSpinBlock::set_twoInt(int pintegralIndex) {
  twoInt = boost::shared_ptr<TwoElectronArray>( &v_2[pintegralIndex], boostutils::null_deleter());
  integralIndex = pintegralIndex;
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
    mpi::broadcast(calc, ar, 0);
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
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allopsUsingCre;
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  
  double* localdata = data; 
  //first make cre

  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(it->second->is_core()) {
      localdata = it->second->allocateOperators(braStateInfo, ketStateInfo, localdata);
      for (int i=0; i<it->second->get_size(); i++)
	for (int j=0; j<it->second->get_local_element(i).size(); j++) {
	  allops.push_back(it->second->get_local_element(i)[j]);
	}
    }
  }
  std::vector<vector<vector<Csf> > > ompladders(numthrds, ladders);
  vector<vector<Csf> > ompdets(numthrds, dets);
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<allops.size(); i++) {
    allops[i]->buildUsingCsf(*this, ompladders[omprank], ompdets[omprank]);
  }
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

  //here enumerate is ordered so that ther cheapest renormalizations will be at the end
  //opTypes OPTYPES[14] {HAM, DES_DESCOMP, CRE_DESCOMP, CRE_CRE_DESCOMP, 
  //DES, DES_DES, CRE_CRECOMP, DES_CRE, DES_CRECOMP, CRE_DES_DESCOMP,  CRE, OVERLAP, CRE_CRE, CRE_DES}; 
  
  //for (int opindex = 0; opindex<14; opindex++) {
  //opTypes ot = OPTYPES[opindex];
  //std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.find(ot);
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core()) {
      for (int i=0; i<it->second->get_size(); i++)
	for (int j=0; j<it->second->get_local_element(i).size(); j++) 
	  allops.push_back(it->second->get_local_element(i)[j]);
    }
  }
  const int quantaSz = newStateInfo->quanta.size ();
  std::multimap<long, int> reorder;
  for (int i=0; i<quantaSz; i++)
    reorder.insert( std::pair<long, int >(newStateInfo->getquantastates(i), i));
  std::vector<int> reorderedVector(quantaSz, 0);
  int index = quantaSz-1;
  for (std::multimap<long, int >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
    reorderedVector[index] = it->second;
    index--;
  }  

  std::vector<long> backupMem(allops.size(), 0);
  for (int i=0; i<allops.size(); i++) {
    allops[i]->Clear();
    backupMem[i] = allops[i]->set_totalMemory();
    allops[i]->set_totalMemory() = 0;
  }

  const std::vector<int>& newQuantaMap = newStateInfo->newQuantaMap;
  
  dmrginp.parallelrenorm->start();
  SplitStackmem();
  //mpi::communicator world;

  //world.barrier();
  //for (int r = 0; r < mpigetsize(); ++r) {
  //if (r == mpigetrank()) {
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<allops.size()*reorderedVector.size(); i++) {
    int opindex = (i)%allops.size(), quantaindex = (i)/allops.size();

    std::vector<int>& colinds = allops[opindex]->getActiveCols(quantaindex);
    for(int qprimeindex = 0; qprimeindex <colinds.size(); qprimeindex++) {
      int qprime = colinds[qprimeindex];
      int Q = newQuantaMap[quantaindex], QPrime = newQuantaMap[qprime];

      double* data = Stackmem[omprank].allocate(get_ketStateInfo().getquantastates(Q)*get_ketStateInfo().getquantastates(QPrime));      
      StackMatrix m(data, get_braStateInfo().getquantastates(Q), get_ketStateInfo().getquantastates(QPrime));
      ::Clear(m);

      allops[opindex]->build(m, Q, QPrime, *this);

      MatrixRotate(rotateMatrix[Q], m, rotateMatrix[QPrime], 
		  allops[opindex]->operator_element(quantaindex, qprime));
      Stackmem[omprank].deallocate(data, m.Storage());
    }
  }
  //}
  //world.barrier();
  //}

  MergeStackmem();
  dmrginp.parallelrenorm->stop();

  for (int i=0; i<allops.size(); i++) {
    allops[i]->set_totalMemory() = backupMem[i];
  }

}


void StackSpinBlock::build_and_renormalise_operators(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
{
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<boost::shared_ptr<StackSparseMatrix> >  allopsOnDisk;
  std::vector<std::string >  fileNames;
  std::vector<opTypes> opTypesOnDisk;

  //these are for the three index operators which are build on the disk
  //and have special build_and_renormalise_operators
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core() && ot >= CRE_CRE_CRE && !(ot == RI_3INDEX || ot == RI_4INDEX) ) {
      for (int i=0; i<it->second->get_size(); i++) {
	for (int j=0; j<it->second->get_local_element(i).size(); j++) {
	  opTypesOnDisk.push_back(ot);
	  allopsOnDisk.push_back(it->second->get_local_element(i)[j]);
	  fileNames.push_back( it->second->get_filename());
	}
      }
    }
  }
  
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it) {
    opTypes ot = it->first;
    if(! it->second->is_core() && ot < CRE_CRE_CRE && !(ot == RI_3INDEX || ot == RI_4INDEX) ) {
      for (int i=0; i<it->second->get_size(); i++)
	for (int j=0; j<it->second->get_local_element(i).size(); j++) {
	  allops.push_back(it->second->get_local_element(i)[j]);
	}
    }
  }

  dmrginp.parallelrenorm->start();
  SplitStackmem();
  if (!dmrginp.do_npdm_in_core()){

    vector<StackSparseMatrix> tmp(numthrds);
#pragma omp parallel 
    {
      for (int I=0; I<allopsOnDisk.size()/numthrds+1; I++) {

	size_t mem = Stackmem[omprank].memused;
	double *ptr = Stackmem[omprank].data+mem;

	if (omprank == 0) {
	  //if the leftblock is not a dot block then check if the relevant operator is on leftblock
	  //if it is then preread it using only thread 0. Reading it using multiple threads is for some reason
	  //causing crashes
	  if (leftBlock->size() > 1) {
	    for (int i=0;i<numthrds; i++) {
	      if (I*numthrds+i < allopsOnDisk.size()) {
		allopsOnDisk[I*numthrds+i]->CleanUp();
		
		int o1 = allopsOnDisk[I*numthrds+i]->get_orbs(0),
		  o2 = allopsOnDisk[I*numthrds+i]->get_orbs(1),
		  o3 = allopsOnDisk[I*numthrds+i]->get_orbs(2);
		
		if (leftBlock->get_op_array(opTypesOnDisk[I*numthrds+i]).has(o1,o2,o3)) {
		  std::vector<boost::shared_ptr<StackSparseMatrix> > sysops = leftBlock->get_op_array(opTypesOnDisk[I*numthrds+i]).get_element(o1,o2,o3);
		  //loop over the sysdotops and pick the one corresponding to op
		  for (int jdx=0; jdx < sysops.size(); jdx++) {
		    boost::shared_ptr<StackSparseMatrix> sysop = sysops[jdx];
		    int len1 = sysop->get_filename().length(),
		      len2 = allopsOnDisk[I*numthrds+i]->get_filename().length();
		    char char1 = sysop->get_filename()[len1-1], char2 = allopsOnDisk[I*numthrds+i]->get_filename()[len2-1];
		    if ( char1 == char2) {
		      bool allocate = sysop->memoryUsed() == 0;
		      //if this is not in memory already, read it from the disk
		      if (allocate) {
			sysop->LoadThreadSafe(true);
			sysop->allocateOperatorMatrix();
			continue;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
#pragma omp barrier

	if (I*numthrds+omprank < allopsOnDisk.size()) {
	  int i = I*numthrds+omprank;
	  allopsOnDisk[i]->allocate(get_braStateInfo(), get_ketStateInfo());
	  allopsOnDisk[i]->build(*this);
	  
	  tmp[omprank] = *allopsOnDisk[i];
	  tmp[omprank].set_totalMemory() = 0; tmp[omprank].set_data(0); tmp[omprank].CleanUp();
	  tmp[omprank].allocate(*bra, *ket);
	  const std::vector<int>& lnewQuantaMap = bra->newQuantaMap;
	  const std::vector<int>& rnewQuantaMap = ket->newQuantaMap;
	  for (int newQ = 0; newQ < lnewQuantaMap.size(); newQ++)
	    for (int newQPrime = 0; newQPrime < rnewQuantaMap.size(); newQPrime++) {
	      if (tmp[omprank].allowed(newQ, newQPrime)) {
		int Q = lnewQuantaMap[newQ], QPrime = rnewQuantaMap[newQPrime];
		MatrixRotate(leftMat[Q], allopsOnDisk[i]->operator()(Q, QPrime), rightMat[QPrime], tmp[omprank](newQ, newQPrime));
	      }
	    }
	}

#pragma omp barrier
	if (omprank == 0) {
	  for (int i=0;i<numthrds; i++) {
	    if (I*numthrds+i < allopsOnDisk.size()) {
	      tmp[i].SaveThreadSafe();
	      tmp[i].CleanUp();
	      allopsOnDisk[I*numthrds+i]->CleanUp();
	    }
	  }
	}

#pragma omp barrier
	Stackmem[omprank].deallocate(ptr, Stackmem[omprank].memused-mem);

      }
    }
  }
  else {
#pragma omp parallel for schedule(dynamic) 
    for (int i=0; i<allopsOnDisk.size(); i++) {
	allopsOnDisk[i]->build_and_renormalise_transform(this, leftMat, bra, rightMat, ket);
    }
  }

#pragma omp parallel for schedule(dynamic) 
  for (int i=0; i<allops.size(); i++) {
    allops[i]->build_and_renormalise_transform(this, leftMat, bra, rightMat, ket);
  }

  MergeStackmem();

  dmrginp.parallelrenorm->stop();

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

  if (localdata != data+totalMemory) {
    pout << "Memory problem in transform_operators"<<endl;
    exit(0);
  }
  pout << "**** STACK MEMORY REMAINING ***** "<<1.0*(Stackmem[omprank].size-Stackmem[omprank].memused)*sizeof(double)/1.e9<<" GB"<<endl;
  
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
    if (dmrginp.do_npdm_in_core() || it->first < CRE_CRE_CRE)
      requiredMemory += it->second->getRequiredMemory(newbraStateInfo, newketStateInfo);
  totalMemory = requiredMemory;
  data = Stackmem[omprank].allocate(requiredMemory);
  double* localdata = data;
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = ops.begin(); it != ops.end(); ++it)
    if (dmrginp.do_npdm_in_core() || it->first < CRE_CRE_CRE)
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
  //p1out << "\t\t\t Building Sum Block " << name << endl;
  leftBlock = &lBlock;
  rightBlock = &rBlock;

  if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= lBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= rBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()!=0 ){
    if(lBlock.nonactive_orbs.size() != rBlock.nonactive_orbs.size()){
      pout << "Nonactive_orbs in left block and right block are different.";
      abort();
    }
    else{
      nonactive_orbs= rBlock.nonactive_orbs;
      for(int i=0; i<nonactive_orbs.size(); i++)
        if(lBlock.nonactive_orbs[i]!= rBlock.nonactive_orbs[i]){
          pout << "Nonactive_orbs in left block and right block are different.";
          abort();
        }
    }
  }

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
      mpi::broadcast(calc, ar, 0);
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
  //p2out << "\t\t\t ";
  //for (int i = 0; i < sites.size(); ++i) p2out << sites[i] << " ";
  //p2out << endl;
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

//Build Sum Block with different quanta num in bra and ket.
void StackSpinBlock::BuildSumBlockSkeleton(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum, bool collectQuanta)
{

  name = get_name();
  if (dmrginp.outputlevel() > 0) 
    pout << "\t\t\t Building Sum Block " << name << endl;
  leftBlock = &lBlock;
  rightBlock = &rBlock;

  if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= lBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()==0 )
    nonactive_orbs= rBlock.nonactive_orbs;
  else if(lBlock.nonactive_orbs.size()!=0 && rBlock.nonactive_orbs.size()!=0 ){
    if(lBlock.nonactive_orbs.size() != rBlock.nonactive_orbs.size()){
      pout << "Nonactive_orbs in left block and right block are different.";
      abort();
    }
    else{
      nonactive_orbs= rBlock.nonactive_orbs;
      for(int i=0; i<nonactive_orbs.size(); i++)
        if(lBlock.nonactive_orbs[i]!= rBlock.nonactive_orbs[i]){
          pout << "Nonactive_orbs in left block and right block are different.";
          abort();
        }
    }
  }


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

  sites.reserve (lBlock.sites.size () + rBlock.sites.size ());
  sites = lBlock.sites;
  copy (rBlock.sites.begin(), rBlock.sites.end (), back_inserter (sites));
  sort(sites.begin(), sites.end());
  complementary_sites = make_complement(sites);
  if (dmrginp.outputlevel() > 0) {
    pout << "\t\t\t ";
    for (int i = 0; i < sites.size(); ++i) pout << sites[i] << " ";
    pout << endl;
  }
  dmrginp.blocksites -> stop();

  dmrginp.statetensorproduct -> start();
  //TODO
  //It is very possible that several quantum numbers for bra.
  //For example, Perturber can have different Spatial Symmetry.
  //TODO
  //What does compState do? 
  //It is not used at all in the whole block code.
  if( condition== PARTICLE_SPIN_NUMBER_CONSTRAINT)
  {
    TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, braquantum[0], EqualQ, braStateInfo);
    TensorProduct (lBlock.ketStateInfo, rBlock.ketStateInfo, ketquantum[0], EqualQ, ketStateInfo);
  }
  else if (condition== NO_PARTICLE_SPIN_NUMBER_CONSTRAINT) 
  {
    TensorProduct (lBlock.braStateInfo, rBlock.braStateInfo, braquantum[0], LessThanQ, braStateInfo);
    TensorProduct (lBlock.ketStateInfo, rBlock.ketStateInfo, ketquantum[0], LessThanQ, ketStateInfo);
  }
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
  //To reduced the number of operators in onepdm calculation, two index operators are removed. However, they are still needed for single orbital sites.
  //For the sum block built from dummy sites (singlet embedding), there are no two index operators.
  if(dmrginp.do_npdm_ops() && (dmrginp.calc_type() == RESTART_ONEPDM || dmrginp.calc_type() == ONEPDM))
    if(sites.size()==1)
    {
      if (!is_direct()) {
        ops[CRE_DES] = make_new_stackop(CRE_DES, true);
        ops[CRE_CRE] = make_new_stackop(CRE_CRE, true);
        ops[DES_CRE] = make_new_stackop(DES_CRE, true);
        ops[DES_DES] = make_new_stackop(DES_DES, true);
      }
      else
      {
        ops[CRE_DES] = make_new_stackop(CRE_DES, false);
        ops[CRE_CRE] = make_new_stackop(CRE_CRE, false);
        ops[DES_CRE] = make_new_stackop(DES_CRE, false);
        ops[DES_DES] = make_new_stackop(DES_DES, false);

      }
    }

  totalMemory = build_iterators();
  if (totalMemory != 0)
    data = Stackmem[omprank].allocate(totalMemory);

  dmrginp.buildblockops -> start();
  build_operators();
  dmrginp.buildblockops -> stop();
  dmrginp.buildsumblock -> stop();
}

//Build Sum Block with different quanta num in bra and ket.
void StackSpinBlock::BuildSumBlock(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum, bool collectQuanta)
{
  if (!(lBlock.integralIndex == rBlock.integralIndex && lBlock.integralIndex == integralIndex))  {
    pout << "The left, right and dot block should use the same integral indices"<<endl;
    pout << "ABORTING!!"<<endl;
    exit(0);
  }
  dmrginp.buildsumblock -> start();
  BuildSumBlockSkeleton(condition, lBlock, rBlock, braquantum,ketquantum, collectQuanta);

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

void StackSpinBlock::multiplyCDD_sum(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{

  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<FUNCTOR2> allfuncs;
  StackWavefunction *v_array;

  initiateMultiThread(v, v_array, numthrds);
  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  SpinQuantum q= -getSpinQuantum( nonactive_orb(0));  
  FUNCTOR2 f = boost::bind(&stackopxop::CDDandoverlap, leftBlock, _1, this, boost::ref(c), v_array, q);
  if (mpigetsize()-1 == mpigetrank()) {
    allops.push_back(rightBlock->get_op_array(OVERLAP).get_element(0).at(0)); allfuncs.push_back(f);//this is just a placeholder function  
  }

  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  FUNCTOR2 f2 = boost::bind(&stackopxop::cdd_dxcdcomp, leftBlock, _1, this, boost::ref(c), v_array,q); 
  for (int i=0; i<rightBlock->get_op_array(DES).get_size(); i++)
    for (int j=0; j<rightBlock->get_op_array(DES).get_local_element(i).size(); j++) {
      allops.push_back(rightBlock->get_op_array(DES).get_local_element(i)[j]);
	    allfuncs.push_back(f2);
    }

  FUNCTOR2 f3 = boost::bind(&stackopxop::cdd_dxcdcomp, rightBlock, _1, this, boost::ref(c), v_array, q ); 
  for (int i=0; i<leftBlock->get_op_array(DES).get_size(); i++)
    for (int j=0; j<leftBlock->get_op_array(DES).get_local_element(i).size(); j++) {
      allops.push_back(leftBlock->get_op_array(DES).get_local_element(i)[j]);
	    allfuncs.push_back(f3);
    }
    
  FUNCTOR2 f4 = boost::bind(&stackopxop::cdd_cxddcomp, leftBlock, _1, this, boost::ref(c), v_array, q);
  for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
    for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
      allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	    allfuncs.push_back(f4);
    }
  
  FUNCTOR2 f5 = boost::bind(&stackopxop::cdd_cxddcomp, rightBlock, _1, this, boost::ref(c), v_array, q);
  for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
    for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
      allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	    allfuncs.push_back(f5);
    }
  
  SplitStackmem();
#pragma omp parallel for  schedule(dynamic) 
  for (int i = 0; i<allops.size(); i++)  {
      allfuncs[i](allops[i]);
  }

  MergeStackmem();
  accumulateMultiThread(v, v_array, numthrds);
  distributedaccumulate(*v);
  dmrginp.oneelecT -> stop();

}

void StackSpinBlock::multiplyCCD_sum(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{

  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<FUNCTOR2> allfuncs;
  StackWavefunction *v_array;

  initiateMultiThread(v, v_array, numthrds);

  dmrginp.oneelecT -> start();
  dmrginp.s0time -> start();

  SpinQuantum q= getSpinQuantum( nonactive_orb(0));  
  FUNCTOR2 f = boost::bind(&stackopxop::CCDandoverlap, leftBlock, _1, this, boost::ref(c), v_array, q);
  if (mpigetsize()-1 == mpigetrank()) {
    allops.push_back(rightBlock->get_op_array(OVERLAP).get_element(0).at(0)); allfuncs.push_back(f);//this is just a placeholder function  
  }


  dmrginp.s0time -> stop();
#ifndef SERIAL
  boost::mpi::communicator world;
  int size = world.size();
#endif

  dmrginp.s1time -> start();
  FUNCTOR2 f2 = boost::bind(&stackopxop::ccd_cxcdcomp, leftBlock, _1, this, boost::ref(c), v_array, q); 
  for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
    for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
      allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	    allfuncs.push_back(f2);
    }

  FUNCTOR2 f3 = boost::bind(&stackopxop::ccd_cxcdcomp, rightBlock, _1, this, boost::ref(c), v_array, q ); 
  for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
    for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
      allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	    allfuncs.push_back(f3);
    }
    
  FUNCTOR2 f4 = boost::bind(&stackopxop::ccd_dxcccomp, leftBlock, _1, this, boost::ref(c), v_array, q);
  for (int i=0; i<rightBlock->get_op_array(DES).get_size(); i++)
    for (int j=0; j<rightBlock->get_op_array(DES).get_local_element(i).size(); j++) {
      allops.push_back(rightBlock->get_op_array(DES).get_local_element(i)[j]);
	    allfuncs.push_back(f4);
    }
  
  FUNCTOR2 f5 = boost::bind(&stackopxop::ccd_dxcccomp, rightBlock, _1, this, boost::ref(c), v_array, q);
  for (int i=0; i<leftBlock->get_op_array(DES).get_size(); i++)
    for (int j=0; j<leftBlock->get_op_array(DES).get_local_element(i).size(); j++) {
      allops.push_back(leftBlock->get_op_array(DES).get_local_element(i)[j]);
	    allfuncs.push_back(f5);
    }
  
  SplitStackmem();
#pragma omp parallel for  schedule(dynamic) 
  for (int i = 0; i<allops.size(); i++)  {
      allfuncs[i](allops[i]);
  }
  MergeStackmem();

  accumulateMultiThread(v, v_array, numthrds);
  distributedaccumulate(*v);

  dmrginp.oneelecT -> stop();

}

int procWithMinOps(std::vector<boost::shared_ptr<StackSparseMatrix> >& allops)
{
  int size = 1;
#ifndef SERIAL
  boost::mpi::communicator world;
  size = calc.size();
#endif
  std::vector<int> numOps(size, 0);

  numOps[mpigetrank()] = allops.size();
  int minproc = 0;
#ifndef SERIAL
  MPI_Allreduce(MPI_IN_PLACE, &numOps[0], size, MPI_INT, MPI_SUM, Calc);
#endif


  for (int i=0; i< size; i++) {
    if (numOps[i] < numOps[minproc])
      minproc = i;
  }
  return minproc;
}

void StackSpinBlock::multiplyH(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{

  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));

  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  StackWavefunction C1, C2;
  C1.initialise(c);
  DCOPY(c.memoryUsed(), c.get_data(), 1, C1.get_data(), 1); 
  C2.initialise(c);
  DCOPY(c.memoryUsed(), c.get_data(), 1, C2.get_data(), 1); 

  //uncollect C[1]
  if (loopBlock == get_leftBlock()) C1.UnCollectQuantaAlongRows(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
  else C1.UnCollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
  if (otherBlock->get_rightBlock() != 0) {
    if (loopBlock == get_leftBlock()) C2.UnCollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
    else C2.UnCollectQuantaAlongRows(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
  }

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<FUNCTOR2> allfuncs;
  std::vector<boost::shared_ptr<StackSparseMatrix> > allops2;
  std::vector<FUNCTOR3> allfuncs2;
  std::vector<boost::shared_ptr<StackSparseMatrix> > allops3;
  std::vector<FUNCTOR3> allfuncs3;

  //accumulate ham
  StackWavefunction* v_array; 
  initiateMultiThread(v, v_array, numthrds);

  FUNCTOR2 f4, f5;
  FUNCTOR3 f4b, f5b;

  int unCollectIndex = 0, collectedIndex =0;
  //finally we will collect the V again and use the twoindex functions
  FUNCTOR2 f1 = boost::bind(&stackopxop::hamandoverlap, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum(), coreEnergy[integralIndex], mpigetsize()-1);
  if (mpigetsize()-1 == mpigetrank()) {
    allops.push_back(rightBlock->get_op_array(OVERLAP).get_element(0).at(0)); allfuncs.push_back(f1);//this is just a placeholder function  
  }
  collectedIndex = allfuncs.size();

  //first line up the functions that use uncollected wavefunctions
  if (loopBlock == rightBlock) {
    if (otherBlock->get_rightBlock() != 0)
      f5 = boost::bind(&stackopxop::cxcddcomp_3index, rightBlock, _1, this, boost::ref(C2), v_array, dmrginp.effective_molecule_quantum());
    else
      f5 = boost::bind(&stackopxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum());

    f4b = boost::bind(&stackopxop::cxcddcomp_3indexElement, leftBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum() );
    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	allfuncs.push_back(f5);
      }
  }
  else {
    if (otherBlock->get_rightBlock() != 0)
      f4 = boost::bind(&stackopxop::cxcddcomp_3index, leftBlock, _1, this, boost::ref(C2), v_array, dmrginp.effective_molecule_quantum() ); 
    else
      f4 = boost::bind(&stackopxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() ); 
    f5b = boost::bind(&stackopxop::cxcddcomp_3indexElement, rightBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum()); 

    for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	allfuncs.push_back(f4);
      }    
  }

  if (otherBlock->get_rightBlock() != 0)
    unCollectIndex = allfuncs.size();
  else {
    collectedIndex = allfuncs.size();
    unCollectIndex = allfuncs.size();
  }

  if (loopBlock == rightBlock) {
    for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	allops2.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	allfuncs2.push_back(f4b);
      }
  }
  else {
    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	allops2.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	allfuncs2.push_back(f5b);
      }
  }

  
  FUNCTOR3 f6b = boost::bind(&stackopxop::cdxcdcomp_3indexElement, otherBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum());
  FUNCTOR3 f7b = boost::bind(&stackopxop::ddxcccomp_3indexElement, otherBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum() );


  //all these will use the threeindex functions
  if (dmrginp.hamiltonian() != HUBBARD) {
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=1; j<loopBlock->get_op_array(CRE_DES).get_local_element(i).size(); j++) {
	allops2.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	allfuncs2.push_back(f6b);
      }

    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++)
      for (int j=1; j<loopBlock->get_op_array(CRE_CRE).get_local_element(i).size(); j++) {
	allops2.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	allfuncs2.push_back(f7b);
      }
  }
  if (dmrginp.hamiltonian() != HUBBARD) {
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=0; j<1; j++) {
	allops2.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	allfuncs2.push_back(f6b);
      }
    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++)
      for (int j=0; j<1; j++) {
	allops2.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	allfuncs2.push_back(f7b);
      }
  }

  std::vector<int> reorderedVector;

  if (loopBlock == leftBlock) {
    const int rightKetOpSz = get_rightBlock()->get_ketStateInfo().quanta.size ();
    std::multimap<long, int> reorder;
    for (int i=0; i<rightKetOpSz; i++)
      reorder.insert( std::pair<long, int >(get_rightBlock()->get_ketStateInfo().getquantastates(i), i));
    reorderedVector.resize(rightKetOpSz, 0);
    int index = rightKetOpSz-1;
    for (std::multimap<long, int >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
      reorderedVector[index] = it->second;
      index--;
    }  
  }
  else {
    const int leftKetOpSz = get_leftBlock()->get_ketStateInfo().quanta.size ();
    std::multimap<long, int> reorder;
    for (int i=0; i<leftKetOpSz; i++)
      reorder.insert( std::pair<long, int >(get_leftBlock()->get_ketStateInfo().getquantastates(i), i));
    reorderedVector.resize(leftKetOpSz, 0);
    int index = leftKetOpSz-1;
    for (std::multimap<long, int >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
      reorderedVector[index] = it->second;
      index--;
    }  
  }

  dmrginp.matmultNum = 0;

  struct timeval start, end;
  gettimeofday(&start, NULL);

  SplitStackmem();
  dmrginp.tensormultiply->start();
  std::vector<int> collected(numthrds, 0);
  std::vector<int> numops(numthrds, 0);
#pragma omp parallel for  schedule(dynamic)
  for (int i = 0; i<allops.size()+(allops2.size()+allops3.size())*reorderedVector.size(); i++)  {

    if (i>=collectedIndex && i<unCollectIndex && collected[omprank] == 0) {
      if (loopBlock == get_leftBlock()) v_array[omprank].UnCollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
      else  v_array[omprank].UnCollectQuantaAlongRows(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
      collected[omprank] = 1;
    }
    else if (i>=unCollectIndex &&  collected[omprank]==0) {
      if (loopBlock == get_leftBlock()) v_array[omprank].UnCollectQuantaAlongRows(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
      else  v_array[omprank].UnCollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
      collected[omprank] = 2;
    }
    else if (i>=unCollectIndex &&  collected[omprank]==1) {
      if (loopBlock == get_leftBlock()) {
	v_array[omprank].CollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), *otherBlock->get_ketStateInfo().unCollectedStateInfo);
	v_array[omprank].UnCollectQuantaAlongRows(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
      }
      else {
	v_array[omprank].CollectQuantaAlongRows(*otherBlock->get_ketStateInfo().unCollectedStateInfo, loopBlock->get_ketStateInfo());
	v_array[omprank].UnCollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
      }
      collected[omprank] = 2;
    }

    if (i<allops.size()) {
      allfuncs[i](allops[i]);
    }
    else if (i < allops.size()+allops2.size()*reorderedVector.size()) {
      int opindex = (i-allops.size())%allops2.size(), quantaindex = (i-allops.size())/allops2.size();
      allfuncs2[opindex](allops2[opindex], reorderedVector[quantaindex]);
    }
    else {
      int opindex = (i-allops.size()-allops2.size()*reorderedVector.size())%allops3.size(), 
	quantaindex = (i-allops.size()-allops2.size()*reorderedVector.size())/allops3.size();
      allfuncs3[opindex](allops3[opindex], reorderedVector[quantaindex]);
    }
  }
  
  for (int i = 0; i<numthrds; i++)  {
    if (collected[i]==2) {
      if (loopBlock == get_leftBlock()) v_array[i].CollectQuantaAlongRows(*loopBlock->get_ketStateInfo().unCollectedStateInfo, otherBlock->get_ketStateInfo());
      else v_array[i].CollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), *loopBlock->get_ketStateInfo().unCollectedStateInfo);
      collected[i] = 0;
    }
    if (collected[i]==1) {
      if (loopBlock == get_leftBlock()) v_array[i].CollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), *otherBlock->get_ketStateInfo().unCollectedStateInfo);
      else v_array[i].CollectQuantaAlongRows(*otherBlock->get_ketStateInfo().unCollectedStateInfo, loopBlock->get_ketStateInfo());
      collected[i] = 0;
    }
  }

  dmrginp.tensormultiply->stop();  

  gettimeofday(&end, NULL);
  //pout <<"cd/cc "<< *dmrginp.cdtime<<"  "<<*dmrginp.cctime<<"  "<<*dmrginp.guesswf;
  //pout <<"  "<< end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec - start.tv_usec)<<endl;

  dmrginp.cdtime->reset();
  dmrginp.cctime->reset();

  MergeStackmem();
  accumulateMultiThread(v, v_array, numthrds);
  distributedaccumulate(*v);

  C2.deallocate();
  C1.deallocate();
}
/*
void StackSpinBlock::multiplyH(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{

  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));

  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;

  StackWavefunction C1, C2;
  C1.initialise(c);
  DCOPY(c.memoryUsed(), c.get_data(), 1, C1.get_data(), 1); 
  C2.initialise(c);
  DCOPY(c.memoryUsed(), c.get_data(), 1, C2.get_data(), 1); 

  //uncollect C[1]
  if (loopBlock == get_leftBlock()) C1.UnCollectQuantaAlongRows(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
  else C1.UnCollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
  if (otherBlock->get_rightBlock() != 0) {
    if (loopBlock == get_leftBlock()) C2.UnCollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
    else C2.UnCollectQuantaAlongRows(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
  }

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<FUNCTOR2> allfuncs;
  std::vector<boost::shared_ptr<StackSparseMatrix> > allops2;
  std::vector<FUNCTOR3> allfuncs2;

  //accumulate ham
  StackWavefunction* v_array; 
  initiateMultiThread(v, v_array, numthrds);

  FUNCTOR2 f4, f5;
  FUNCTOR3 f4b, f5b;

  int unCollectIndex = 0, collectedIndex =0;
  //finally we will collect the V again and use the twoindex functions
  FUNCTOR2 f1 = boost::bind(&stackopxop::hamandoverlap, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum(), coreEnergy[integralIndex], mpigetsize()-1);
  if (mpigetsize()-1 == mpigetrank()) {
    allops.push_back(rightBlock->get_op_array(OVERLAP).get_element(0).at(0)); allfuncs.push_back(f1);//this is just a placeholder function  
  }
  collectedIndex = allfuncs.size();

  //first line up the functions that use uncollected wavefunctions
  if (loopBlock == rightBlock) {
    if (otherBlock->get_rightBlock() != 0)
      f5 = boost::bind(&stackopxop::cxcddcomp_3index, rightBlock, _1, this, boost::ref(C2), v_array, dmrginp.effective_molecule_quantum());
    else
      f5 = boost::bind(&stackopxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum());

    f4b = boost::bind(&stackopxop::cxcddcomp_3indexElement, leftBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum() );

    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f5);
      }
  } else {
    // loopBlock is leftBlock
    if (otherBlock->get_rightBlock() != 0)
      f4 = boost::bind(&stackopxop::cxcddcomp_3index, leftBlock, _1, this, boost::ref(C2), v_array, dmrginp.effective_molecule_quantum() ); 
    else
      f4 = boost::bind(&stackopxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() );

    f5b = boost::bind(&stackopxop::cxcddcomp_3indexElement, rightBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum()); 

    for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f4);
      }    
  }

  if (otherBlock->get_rightBlock() != 0)
    unCollectIndex = allfuncs.size();
  else {
    collectedIndex = allfuncs.size();
    unCollectIndex = allfuncs.size();
  }

  if (loopBlock == rightBlock) {
    for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops2.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs2.push_back(f4b);
      }
  }
  else {
    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops2.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs2.push_back(f5b);
      }
  }

  FUNCTOR3 f6b = boost::bind(&stackopxop::cdxcdcomp_3indexElement, otherBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum());
  FUNCTOR3 f7b = boost::bind(&stackopxop::ddxcccomp_3indexElement, otherBlock, _1, this, boost::ref(C1), v_array, _2, dmrginp.effective_molecule_quantum() );


  //all these will use the threeindex functions
  if (dmrginp.hamiltonian() != HUBBARD) {
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=1; j<loopBlock->get_op_array(CRE_DES).get_local_element(i).size(); j++) {
	      allops2.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	      allfuncs2.push_back(f6b);
      }

    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++)
      for (int j=1; j<loopBlock->get_op_array(CRE_CRE).get_local_element(i).size(); j++) {
	      allops2.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	      allfuncs2.push_back(f7b);
      }
    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++)
      for (int j=0; j<1; j++) {
	      allops2.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	      allfuncs2.push_back(f7b);
      }
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=0; j<1; j++) {
	      allops2.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	      allfuncs2.push_back(f6b);
      }
  }

  std::vector<int> reorderedVector;

  if (loopBlock == leftBlock) {
    const int rightKetOpSz = get_rightBlock()->get_ketStateInfo().quanta.size ();
    std::multimap<long, int> reorder;
    for (int i=0; i<rightKetOpSz; i++)
      reorder.insert( std::pair<long, int >(get_rightBlock()->get_ketStateInfo().getquantastates(i), i));
    reorderedVector.resize(rightKetOpSz, 0);
    int index = rightKetOpSz-1;
    for (std::multimap<long, int >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
      reorderedVector[index] = it->second;
      index--;
    }  
  }
  else {
    const int leftKetOpSz = get_leftBlock()->get_ketStateInfo().quanta.size ();
    std::multimap<long, int> reorder;
    for (int i=0; i<leftKetOpSz; i++)
      reorder.insert( std::pair<long, int >(get_leftBlock()->get_ketStateInfo().getquantastates(i), i));
    reorderedVector.resize(leftKetOpSz, 0);
    int index = leftKetOpSz-1;
    for (std::multimap<long, int >::iterator it = reorder.begin(); it!=reorder.end(); it++) {
      reorderedVector[index] = it->second;
      index--;
    }  
  }

  dmrginp.matmultNum = 0;

  struct timeval start, end;
  gettimeofday(&start, NULL);

  SplitStackmem();
  dmrginp.tensormultiply->start();
  std::vector<int> collected(numthrds, 0);
  std::vector<int> numops(numthrds, 0);
  //pout << "allops.size() " << allops.size() << endl;
  //pout << "allops2.size() " << allops2.size() << endl;
#pragma omp parallel for  schedule(dynamic) 
  for (int i = 0; i<allops.size()+allops2.size()*reorderedVector.size(); i++)  {

    if (i>=collectedIndex && i<unCollectIndex && collected[omprank] == 0) {
      if (loopBlock == get_leftBlock()) v_array[omprank].UnCollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
      else  v_array[omprank].UnCollectQuantaAlongRows(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
      collected[omprank] = 1;
    }
    else if (i>=unCollectIndex &&  collected[omprank]==0) {
      if (loopBlock == get_leftBlock()) v_array[omprank].UnCollectQuantaAlongRows(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
      else  v_array[omprank].UnCollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
      collected[omprank] = 2;
    }
    else if (i>=unCollectIndex &&  collected[omprank]==1) {
      if (loopBlock == get_leftBlock()) {
	v_array[omprank].CollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), *otherBlock->get_ketStateInfo().unCollectedStateInfo);
	v_array[omprank].UnCollectQuantaAlongRows(loopBlock->get_ketStateInfo(), otherBlock->get_ketStateInfo());
      }
      else {
	v_array[omprank].CollectQuantaAlongRows(*otherBlock->get_ketStateInfo().unCollectedStateInfo, loopBlock->get_ketStateInfo());
	v_array[omprank].UnCollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), loopBlock->get_ketStateInfo());
      }
      collected[omprank] = 2;
    }

    if (i<allops.size()) {
      allfuncs[i](allops[i]);
    }
    else {
      int opindex = (i-allops.size())%allops2.size(), quantaindex = (i-allops.size())/allops2.size();
      allfuncs2[opindex](allops2[opindex], reorderedVector[quantaindex]);
    }
  }

  
  for (int i = 0; i<numthrds; i++)  {
    if (collected[i]==2) {
      if (loopBlock == get_leftBlock()) v_array[i].CollectQuantaAlongRows(*loopBlock->get_ketStateInfo().unCollectedStateInfo, otherBlock->get_ketStateInfo());
      else v_array[i].CollectQuantaAlongColumns(otherBlock->get_ketStateInfo(), *loopBlock->get_ketStateInfo().unCollectedStateInfo);
      collected[i] = 0;
    }
    if (collected[i]==1) {
      if (loopBlock == get_leftBlock()) v_array[i].CollectQuantaAlongColumns(loopBlock->get_ketStateInfo(), *otherBlock->get_ketStateInfo().unCollectedStateInfo);
      else v_array[i].CollectQuantaAlongRows(*otherBlock->get_ketStateInfo().unCollectedStateInfo, loopBlock->get_ketStateInfo());
      collected[i] = 0;
    }
  }

  dmrginp.tensormultiply->stop();  

  gettimeofday(&end, NULL);
  //pout <<"cd/cc "<< *dmrginp.cdtime<<"  "<<*dmrginp.cctime<<"  "<<*dmrginp.guesswf;
  //pout <<"  "<< end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec - start.tv_usec)<<endl;

  dmrginp.cdtime->reset();
  dmrginp.cctime->reset();

  MergeStackmem();
  accumulateMultiThread(v, v_array, numthrds);
  distributedaccumulate(*v);

  C2.deallocate();
  C1.deallocate();
}
*/

void StackSpinBlock::multiplyH_2index(StackWavefunction& c, StackWavefunction* v, int num_threads) const
{
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));
  StackSpinBlock* loopBlock=(leftBlock->is_loopblock()) ? leftBlock : rightBlock;
  StackSpinBlock* otherBlock = loopBlock == leftBlock ? rightBlock : leftBlock;


  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;
  std::vector<FUNCTOR2> allfuncs;

  //accumulate ham
  StackWavefunction* v_array; 
  initiateMultiThread(v, v_array, numthrds);

  FUNCTOR2 f4, f5;

  //finally we will collect the V again and use the twoindex functions
  FUNCTOR2 f1 = boost::bind(&stackopxop::hamandoverlap, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum(), coreEnergy[integralIndex], mpigetsize()-1);
  if (mpigetsize()-1 == mpigetrank()) {
    allops.push_back(rightBlock->get_op_array(OVERLAP).get_element(0).at(0)); allfuncs.push_back(f1);//this is just a placeholder function  
  }

  if (loopBlock == rightBlock) {
    f5 = boost::bind(&stackopxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum());
    f4 = boost::bind(&stackopxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() );
    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f5);
      }
  }
  else {
    f4 = boost::bind(&stackopxop::cxcddcomp, leftBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() ); 
    f5 = boost::bind(&stackopxop::cxcddcomp, rightBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum()); 

    for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
        allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f4);
      }    
  }

  if (loopBlock == rightBlock) {
    for (int i=0; i<rightBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<rightBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops.push_back(rightBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f4);
      }
  }
  else {
    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	      allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f5);
      }
  }

  FUNCTOR2 f6 = boost::bind(&stackopxop::cdxcdcomp, otherBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum());
  FUNCTOR2 f7 = boost::bind(&stackopxop::ddxcccomp, otherBlock, _1, this, boost::ref(c), v_array, dmrginp.effective_molecule_quantum() );

  //all these will use the threeindex functions
  if (dmrginp.hamiltonian() != HUBBARD) {
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++)
      for (int j=1; j<loopBlock->get_op_array(CRE_DES).get_local_element(i).size(); j++) {
	      allops.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	      allfuncs.push_back(f6);
      }

    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++)
      for (int j=1; j<loopBlock->get_op_array(CRE_CRE).get_local_element(i).size(); j++) {
	      allops.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	      allfuncs.push_back(f7);
      }
    for (int i=0; i<loopBlock->get_op_array(CRE_CRE).get_size(); i++) {
	    allops.push_back(loopBlock->get_op_array(CRE_CRE).get_local_element(i)[0]);
	    allfuncs.push_back(f7);
    }
    for (int i=0; i<loopBlock->get_op_array(CRE_DES).get_size(); i++) {
	    allops.push_back(loopBlock->get_op_array(CRE_DES).get_local_element(i)[0]);
      allfuncs.push_back(f6);
    }
  }
  dmrginp.matmultNum = 0;

  struct timeval start, end;
  gettimeofday(&start, NULL);

  SplitStackmem();
  dmrginp.tensormultiply->start();
  std::vector<int> collected(numthrds, 0);
  std::vector<int> numops(numthrds, 0);

  //pout << "allops.size() " << allops.size() << endl;
#pragma omp parallel for  schedule(dynamic) 
  for (int i = 0; i<allops.size(); i++)  {
      allfuncs[i](allops[i]);
  }

  dmrginp.tensormultiply->stop();  

  gettimeofday(&end, NULL);
  //pout <<"cd/cc "<< *dmrginp.cdtime<<"  "<<*dmrginp.cctime<<"  "<<*dmrginp.guesswf;
  //pout <<"  "<< end.tv_sec-start.tv_sec + 1e-6*(end.tv_usec - start.tv_usec)<<endl;

  dmrginp.cdtime->reset();
  dmrginp.cctime->reset();

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
    allops.push_back(loopBlock->get_op_array(HAM).get_element(0).at(0)); allfuncs.push_back(f1); //the function is placeholder
    allops.push_back(otherBlock->get_op_array(HAM).get_element(0).at(0)); allfuncs.push_back(f2);//the function is placeholder
  }


  SplitStackmem();
  //dmrginp.tensormultiply->start();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i<allops.size(); i++)  {
    allfuncs[i](allops[i]);
  }
  //dmrginp.tensormultiply->stop();  
  MergeStackmem();

  for (int i=0; i<numthrds; i++)
    e += e_array[i];
  delete [] e_array;
  distributedaccumulate(e);

}

void StackSpinBlock::BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, bool haveCompops, const bool haveNormops)
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
  default_op_components(false, haveNormops, haveCompops, true);
  
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

  //p3out << "\t\t\t time in slater distribution " << slatertimer.elapsedwalltime() << endl;

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
    if (dmrginp.calc_type() == COMPRESS || dmrginp.calc_type() == RESPONSE)
      singleSiteBlocks[i] = StackSpinBlock(i, i, integralIndex, false);
    else
      singleSiteBlocks[i] = StackSpinBlock(i, i, integralIndex, true);
  }
}

}
