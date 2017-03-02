/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include "Stackspinblock.h"
#include "IntegralMatrix.h"
#include "screen.h"
#include "stackopxop.h"
#include "operatorfunctions.h"
#include "Stackwavefunction.h"
#include <boost/format.hpp>
#include <sstream>
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
#include <stdio.h>
#include <tuple>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted{
using namespace operatorfunctions;

void StackSpinBlock::make_iterator(StackSpinBlock& b, opTypes op, int* data, int index, int numIndices) {

  std::vector<int> oneindex (numIndices, 0.0);
  std::vector<std::pair<int, int> > twoindex(numIndices, std::pair<int, int>(0,0));
  std::map<std::tuple<int, int, int>, int > threeindex;
  
  if (index == 1) {
    for (int i=0; i<numIndices; i++) {
      oneindex[i] = data[i];
    }
  }
  else if (index == 2){
    for (int i=0; i<numIndices; i++) {
      twoindex[i].first = data[2*i];
      twoindex[i].second = data[2*i+1];
    }
  }
  else if (index == 3){
    for (int i=0; i<numIndices; i++) 
      threeindex[std::make_tuple(data[4*i], data[4*i+1], data[4*i+2])] = data[4*i+3];
  }

  b.ops[op]->build_iterators(b, oneindex, twoindex, threeindex);
  
} 

void StackSpinBlock::restore (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char* name)
{
  dmrginp.diski->start();
  Timer disktimer;
  std::string file[numthrds];

  for (int i=0; i<numthrds; i++) {
    if (forward)
      file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % i % ".tmp" );
    else
      file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % i % ".tmp" );
  }

  p1out << "\t\t\t Restoring block file :: " << file[0] << endl;


  int lstate =  left;
  int rstate =  right;

  if (mpigetrank() == 0) {
    StateInfo::restore(forward, sites, b.braStateInfo, lstate);
    StateInfo::restore(forward, sites, b.ketStateInfo, rstate);
  }
  
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(calc, b.braStateInfo, 0);
  mpi::broadcast(calc, b.ketStateInfo, 0);
#endif


  dmrginp.rawdatai->start();
  
  int* initialData = new int[31];
  int allindexsize;
  
  FILE *fp[numthrds];
  fp[0] = fopen(file[0].c_str(), "rb");

  fread(initialData, sizeof(int), 31, fp[0]);
  fread(&allindexsize, sizeof(int), 1, fp[0]);
  int* allindices = new int[allindexsize];//large data
  fread(allindices, sizeof(int), allindexsize, fp[0]);

  fread(&(b.totalMemory), sizeof(long), 1, fp[0]);
  b.data = Stackmem[omprank].allocate(b.totalMemory);

  double walltime = globaltimer.totalwalltime();
  fread(b.data, sizeof(double), b.totalMemory, fp[0]);
  fclose(fp[0]);

  pout << str( boost::format("Read  %-10.4fG of data in  %-10.4f s\n") % (b.totalMemory*sizeof(double)/1.e9) % (globaltimer.totalwalltime()-walltime));

  dmrginp.rawdatai->stop();

  int dataperthrd = b.totalMemory/numthrds;
  int dataonlastthrd = b.totalMemory/numthrds + b.totalMemory%numthrds;

    
  //nowunpack the first 25 integers
  b.localstorage         = initialData[0] == 1 ? true : false;
  b.name                 = initialData[1];
  b.complementary        = initialData[2] == 1 ? true : false;
  b.normal               = initialData[3] == 1 ? true : false;
  b.direct               = initialData[4] == 1 ? true : false;
  b.loopblock            = initialData[5] == 1 ? true : false;
  int numsites           = initialData[6];
  int compsites          = initialData[7];
  b.integralIndex        = initialData[8];
  int numcre             = initialData[9];
  int numdes             = initialData[10];
  int numcredes          = initialData[11];
  int numdescre          = initialData[12];
  int numcrecre          = initialData[13];
  int numdesdes          = initialData[14];
  int numcredescomp      = initialData[15];
  int numdescrecomp      = initialData[16];
  int numdesdescomp      = initialData[17];
  int numcrecrecomp      = initialData[18];
  int numcrecredescomp   = initialData[19];
  int numcredesdescomp   = initialData[20];
  int numham             = initialData[21];
  int numoverlap         = initialData[22];
  int numcdd_sum         = initialData[23];
  int numcdd_cd          = initialData[24];
  int numcdd_dd          = initialData[25];
  int numccd_sum         = initialData[26];
  int numccd_cd          = initialData[27];
  int numccd_cc          = initialData[28];
  int numri3index        = initialData[29];
  int numri4index        = initialData[30];


  dmrginp.readmakeiter->start();

  b.sites.resize(numsites, 0);
  b.complementary_sites.resize(compsites, 0);

  int index = 0;
  for (int i=0; i<numsites; i++)
    b.sites[i] = allindices[index+i];
  index += numsites;

  for (int i=0; i<compsites; i++)
    b.complementary_sites[i] = allindices[index+i];
  index += compsites;


  if (numham           != -1) {b.ops[HAM]           =   make_new_stackop(HAM, true);}
  if (numcre           != -1) {b.ops[CRE]           =   make_new_stackop(CRE, true);}
  if (numcrecre        != -1) {b.ops[CRE_CRE]           =   make_new_stackop(CRE_CRE, true);}
  if (numdesdescomp    != -1) {b.ops[DES_DESCOMP]           =   make_new_stackop(DES_DESCOMP, true);}
  if (numcredes        != -1) {b.ops[CRE_DES]           =   make_new_stackop(CRE_DES, true);}
  if (numcredescomp    != -1) {b.ops[CRE_DESCOMP]           =   make_new_stackop(CRE_DESCOMP, true);}
  if (numcrecredescomp != -1) {b.ops[CRE_CRE_DESCOMP]           =   make_new_stackop(CRE_CRE_DESCOMP, true);}
  if (numdes           != -1) {b.ops[DES]           =   make_new_stackop(DES, true);}
  if (numdesdes        != -1) {b.ops[DES_DES]           =   make_new_stackop(DES_DES, true);}
  if (numcrecrecomp    != -1) {b.ops[CRE_CRECOMP]           =   make_new_stackop(CRE_CRECOMP, true);}
  if (numdescre        != -1) {b.ops[DES_CRE]           =   make_new_stackop(DES_CRE, true);}
  if (numdescrecomp    != -1) {b.ops[DES_CRECOMP]           =   make_new_stackop(DES_CRECOMP, true);}
  if (numcredesdescomp != -1) {b.ops[CRE_DES_DESCOMP]           =   make_new_stackop(CRE_DES_DESCOMP, true);}
  if (numoverlap       != -1) {b.ops[OVERLAP]           =   make_new_stackop(OVERLAP, true);}
  if (numcdd_sum       != -1) {b.ops[CDD_SUM]           =   make_new_stackop(CDD_SUM, true);}
  if (numcdd_cd        != -1) {b.ops[CDD_CRE_DESCOMP]   =   make_new_stackop(CDD_CRE_DESCOMP, true);}
  if (numcdd_dd        != -1) {b.ops[CDD_DES_DESCOMP]   =   make_new_stackop(CDD_DES_DESCOMP, true);}
  if (numccd_sum       != -1) {b.ops[CCD_SUM]           =   make_new_stackop(CCD_SUM, true);}
  if (numccd_cd        != -1) {b.ops[CCD_CRE_DESCOMP]   =   make_new_stackop(CCD_CRE_DESCOMP, true);}
  if (numccd_cc        != -1) {b.ops[CCD_CRE_CRECOMP]   =   make_new_stackop(CCD_CRE_CRECOMP, true);}
  if (numri3index      != -1) {b.ops[RI_3INDEX]         =   make_new_stackop(RI_3INDEX, true);}
  if (numri4index      != -1) {b.ops[RI_4INDEX]         =   make_new_stackop(RI_4INDEX, true);}

  if (b.localstorage) b.setstoragetype(LOCAL_STORAGE);
  else b.setstoragetype(DISTRIBUTED_STORAGE);

  //this should be in the same order as opTypes are written in enumerator.h file
  //
  if (numham           != -1) { make_iterator(b, HAM,             &allindices[index], 1,  numham);           index+=numham;} 
  if (numcre           != -1) { make_iterator(b, CRE,             &allindices[index], 1,  numcre);           index+=numcre;}
  if (numcrecre        != -1) { make_iterator(b, CRE_CRE,         &allindices[index], 2,  numcrecre);        index+=2*numcrecre;}
  if (numdesdescomp    != -1) { make_iterator(b, DES_DESCOMP,     &allindices[index], 2,  numdesdescomp);    index+=2*numdesdescomp;}
  if (numcredes        != -1) { make_iterator(b, CRE_DES,         &allindices[index], 2,  numcredes);        index+=2*numcredes;}
  if (numcredescomp    != -1) { make_iterator(b, CRE_DESCOMP,     &allindices[index], 2,  numcredescomp);    index+=2*numcredescomp;}
  if (numcrecredescomp != -1) { make_iterator(b, CRE_CRE_DESCOMP, &allindices[index], 1,  numcrecredescomp); index+=numcrecredescomp;}
  if (numdes           != -1) { make_iterator(b, DES,             &allindices[index], 1,  numdes);           index+=numdes;}
  if (numdesdes        != -1) { make_iterator(b, DES_DES,         &allindices[index], 2,  numdesdes);        index+=2*numdesdes;}
  if (numcrecrecomp    != -1) { make_iterator(b, CRE_CRECOMP,     &allindices[index], 2,  numcrecrecomp);    index+=2*numcrecrecomp;}
  if (numdescre        != -1) { make_iterator(b, DES_CRE,         &allindices[index], 2,  numdescre);        index+=2*numdescre;}
  if (numdescrecomp    != -1) { make_iterator(b, DES_CRECOMP,     &allindices[index], 2,  numdescrecomp);    index+=2*numdescrecomp;}
  if (numcredesdescomp != -1) { make_iterator(b, CRE_DES_DESCOMP, &allindices[index], 1,  numcredesdescomp); index+=numcredesdescomp;}
  if (numoverlap       != -1) { make_iterator(b, OVERLAP,         &allindices[index], 1,  numoverlap);       index+=numoverlap;}
  if (numcdd_sum       != -1) { make_iterator(b, CDD_SUM,         &allindices[index], 1,  numcdd_sum);       index+=numcdd_sum;}
  if (numcdd_cd        != -1) { make_iterator(b, CDD_CRE_DESCOMP, &allindices[index], 1,  numcdd_cd );       index+=numcdd_cd ;}
  if (numcdd_dd        != -1) { make_iterator(b, CDD_DES_DESCOMP, &allindices[index], 1,  numcdd_dd );       index+=numcdd_dd ;}
  if (numccd_sum       != -1) { make_iterator(b, CCD_SUM,         &allindices[index], 1,  numccd_sum);       index+=numccd_sum;}
  if (numccd_cd        != -1) { make_iterator(b, CCD_CRE_DESCOMP, &allindices[index], 1,  numccd_cd );       index+=numccd_cd ;}
  if (numccd_cc        != -1) { make_iterator(b, CCD_CRE_CRECOMP, &allindices[index], 1,  numccd_cc );       index+=numccd_cc ;}
  if (numri3index      != -1) { b.ops[RI_3INDEX]->build_iterators(b, false);}
  if (numri4index      != -1) { b.ops[RI_4INDEX]->build_iterators(b, false);}
  dmrginp.readmakeiter->stop();


  dmrginp.rawdatai->start();
  dmrginp.rawdatai->stop();

  dmrginp.readallocatemem->start();
  double* localdata = b.data; 
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = b.ops.begin(); it != b.ops.end(); ++it)
  {
    if(it->second->is_core() && it->first != RI_3INDEX && it->first != RI_4INDEX) {

      for (int i=0; i<it->second->get_size(); i++) {
	int vecsize = it->second->get_local_element(i).size();
	for (int j=0; j<vecsize; j++) {
	  it->second->get_local_element(i)[j]->set_data(localdata);
	  localdata = it->second->get_local_element(i)[j]->allocate(b.braStateInfo, b.ketStateInfo, localdata);
	}
      }
    }
  }
  dmrginp.readallocatemem->stop();
  dmrginp.diski->stop();

  delete [] allindices;
  delete [] initialData;

}
  
void StackSpinBlock::store (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char *name)
{
  dmrginp.disko->start();
  Timer disktimer;
  std::string file[numthrds];

  for (int i=0; i<numthrds; i++) {
    if (dmrginp.spinAdapted()) {
      if (forward)
	file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % i % ".tmp" );
      else
	file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % i % ".tmp" );
    }
    else {
      if (forward)
	file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% (sites[0]/2) % "." % (sites[sites.size()-1]/2) % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % i % ".tmp" );
      else
	file[i] = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% (sites[0]/2) % "." % (sites[sites.size()-1]/2) % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % i % ".tmp" );
    }
  }
  
  p1out << "\t\t\t Saving block file :: " << file[0] << endl;

  int lstate =  left;
  int rstate =  right;
  
  if (mpigetrank()==0) {
    StateInfo::store(forward, sites, b.braStateInfo, lstate);
    StateInfo::store(forward, sites, b.ketStateInfo, rstate);
  }


  int* initialData = new int[31];

  //nowunpack the first 23 integers
  initialData[0]    = b.localstorage      == 1 ? true : false;
  initialData[1]    = b.name ;
  initialData[2]    = b.complementary     == 1 ? true : false;
  initialData[3]    = b.normal            == 1 ? true : false;
  initialData[4]    = b.direct            == 1 ? true : false;
  initialData[5]    = b.loopblock         == 1 ? true : false;
  initialData[6]    = b.sites.size();
  initialData[7]    = b.complementary_sites.size();
  initialData[8]    = b.integralIndex ;
  initialData[9]    = b.has(CRE)                 ?  b.ops[CRE]->size()             : -1 ;
  initialData[10]   = b.has(DES)                 ?  b.ops[DES]->size()             : -1;
  initialData[11]   = b.has(CRE_DES)             ?  b.ops[CRE_DES]->size()         : -1;
  initialData[12]   = b.has(DES_CRE)             ?  b.ops[DES_CRE]->size()         : -1;
  initialData[13]   = b.has(CRE_CRE)             ?  b.ops[CRE_CRE]->size()         : -1;
  initialData[14]   = b.has(DES_DES)             ?  b.ops[DES_DES]->size()         : -1;
  initialData[15]   = b.has(CRE_DESCOMP)         ?  b.ops[CRE_DESCOMP]->size()     : -1;
  initialData[16]   = b.has(DES_CRECOMP)         ?  b.ops[DES_CRECOMP]->size()     : -1;
  initialData[17]   = b.has(DES_DESCOMP)         ?  b.ops[DES_DESCOMP]->size()     : -1;
  initialData[18]   = b.has(CRE_CRECOMP)         ?  b.ops[CRE_CRECOMP]->size()     : -1;
  initialData[19]   = b.has(CRE_CRE_DESCOMP)     ?  b.ops[CRE_CRE_DESCOMP]->size() : -1;
  initialData[20]   = b.has(CRE_DES_DESCOMP)     ?  b.ops[CRE_DES_DESCOMP]->size() : -1;
  initialData[21]   = b.has(HAM)                 ?  b.ops[HAM]->size()             : -1;
  initialData[22]   = b.has(OVERLAP)             ?  b.ops[OVERLAP]->size()         : -1;
  initialData[23]   = b.has(CDD_SUM)             ?  b.ops[CDD_SUM]->size()         : -1;
  initialData[24]   = b.has(CDD_CRE_DESCOMP)     ?  b.ops[CDD_CRE_DESCOMP]->size() : -1;
  initialData[25]   = b.has(CDD_DES_DESCOMP)     ?  b.ops[CDD_DES_DESCOMP]->size() : -1;
  initialData[26]   = b.has(CCD_SUM)             ?  b.ops[CCD_SUM]->size()         : -1;
  initialData[27]   = b.has(CCD_CRE_DESCOMP)     ?  b.ops[CCD_CRE_DESCOMP]->size() : -1;
  initialData[28]   = b.has(CCD_CRE_CRECOMP)     ?  b.ops[CCD_CRE_CRECOMP]->size() : -1;
  initialData[29]   = b.has(RI_3INDEX)           ?  b.ops[RI_3INDEX]->size()       : -1;
  initialData[30]   = b.has(RI_4INDEX)           ?  b.ops[RI_4INDEX]->size()       : -1;

  std::vector<int> allindices;
  allindices.insert(allindices.end(), b.sites.begin(), b.sites.end());
  allindices.insert(allindices.end(), b.complementary_sites.begin(), b.complementary_sites.end());

  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = b.ops.begin(); it != b.ops.end(); ++it)
  {
    if (it->first != RI_3INDEX && it->first != RI_4INDEX) {
      std::vector<int> indices = it->second->get_global_array();
      allindices.insert(allindices.end(), indices.begin(), indices.end());
    }
  }

  dmrginp.rawdatao->start();
  FILE *fp[numthrds];
  fp[0] = fopen(file[0].c_str(), "wb");
  int size = allindices.size();
  fwrite(initialData, sizeof(int), 31, fp[0]);
  fwrite(&size, sizeof(int), 1, fp[0]);
  fwrite(&allindices[0], sizeof(int), allindices.size(), fp[0]);

  

  fwrite(&b.totalMemory, sizeof(long), 1, fp[0]);

  int dataperthrd = b.totalMemory/numthrds;
  int dataonlastthrd = b.totalMemory/numthrds + b.totalMemory%numthrds;

  double walltime = globaltimer.totalwalltime();
  fwrite(b.data, sizeof(double), b.totalMemory, fp[0]);
  pout <<str( boost::format("Wrote  %-10.4fG of data in  %-10.4f s\n") % (b.totalMemory*sizeof(double)/1.e9) % (globaltimer.totalwalltime()-walltime) );

  fclose(fp[0]);


  dmrginp.rawdatao->stop();

  delete [] initialData;
  dmrginp.disko->stop();

}


void StackSpinBlock::sendcompOps(StackOp_component_base& opcomp, int I, int J, int optype, int compsite)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    calc.send(processorindex(compsite), optype+i*10+1000*J+100000*I, *oparray[i]);

    //now broadcast the data
    //MPI::COMM_WORLD.Bcast(oparray[i]->get_data(), oparray[i]->memoryUsed(), MPI_DOUBLE, trimap_2d(I, J, dmrginp.last_site()));
    MPI_Send(oparray[i]->get_data(), oparray[i]->memoryUsed(), MPI_DOUBLE, processorindex(compsite), optype+i*10+1000*J+100000*I, Calc);
  }  
#endif
}

void StackSpinBlock::recvcompOps(StackOp_component_base& opcomp, int I, int J, int optype)
{
#ifndef SERIAL
  boost::mpi::communicator world;
  std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = opcomp.get_element(I,J);
  for(int i=0; i<oparray.size(); i++) {
    calc.recv(processorindex(trimap_2d(I, J, dmrginp.last_site())), optype+i*10+1000*J+100000*I, *oparray[i]);

    double* data = Stackmem[omprank].allocate(oparray[i]->memoryUsed());
    if (additionalMemory == 0) 
      additionaldata=data; 
    additionalMemory+=oparray[i]->memoryUsed();
    
    oparray[i]->set_data(data);
    oparray[i]->allocateOperatorMatrix();

    //now broadcast the data
    MPI_Recv(oparray[i]->get_data(), oparray[i]->memoryUsed(), MPI_DOUBLE, processorindex(trimap_2d(I, J, dmrginp.last_site())), optype+i*10+1000*J+100000*I, Calc,MPI_STATUS_IGNORE);
  }
#endif
}


void StackSpinBlock::removeAdditionalOps() 
{
  Stackmem[omprank].deallocate(additionaldata, additionalMemory);
}

void StackSpinBlock::addOneIndexNormOps() {
  // used in NPDM code
#ifndef SERIAL
  boost::mpi::communicator world;
  if (calc.size() == 1)
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
	mpi::broadcast(calc, *op, processorindex(sites[i]));
	
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
	MPI_Bcast(op->get_data(), op->memoryUsed(), MPI_DOUBLE, processorindex(sites[i]), Calc);
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
	mpi::broadcast(calc, *op, processorindex(sites[i]));
	
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
	MPI_Bcast(op->get_data(), op->memoryUsed(), MPI_DOUBLE, processorindex(sites[i]), Calc);
      }
    }
  }
#endif
}

void StackSpinBlock::addOneIndexOps() 
{
#ifndef SERIAL
  boost::mpi::communicator world;
  if (calc.size() == 1)
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
	mpi::broadcast(calc, *op, processorindex(sites[i]));
	
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
	MPI_Bcast(op->get_data(), op->memoryUsed(), MPI_DOUBLE, processorindex(sites[i]), Calc);
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
	mpi::broadcast(calc, *op, processorindex(sites[i]));
	
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
	MPI_Bcast(op->get_data(), op->memoryUsed(), MPI_DOUBLE, processorindex(sites[i]), Calc);
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
	    mpi::broadcast(calc, *oparray[iproc], fromproc);
	    
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
	    MPI_Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc, Calc);
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
	    mpi::broadcast(calc, *oparray[iproc], fromproc);
	    
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
	    MPI_Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc, Calc);
	  }
	}
      }
    }
  }
#endif
}

void StackSpinBlock::messagePassTwoIndexOps()
{
#ifndef SERIAL    
  boost::mpi::communicator world;
  int length = dmrginp.last_site();
  vector<int> dotindice;
  dotindice.push_back((sites[0] == 0) ? complementary_sites[0] : complementary_sites[complementary_sites.size()-1]);
  if (!dmrginp.spinAdapted()) { // when non-spinadapted, sites are spin orbitals
    dotindice.push_back((sites[0] == 0) ? complementary_sites[1] : complementary_sites[complementary_sites.size()-2]);    
  }
  for (int idx = 0; idx < dotindice.size(); ++idx) {
    int dotopindex = dotindice[idx];
    int I = dotopindex;
    
    for (int idx1 = idx; idx1 < dotindice.size(); ++idx1) {
      //CDII should be broadcast to all procs
      int J = dotindice[idx1];
      if (J > I) {
        int K = J; J = I; I = K;
      }
      if (has(CRE_DESCOMP)) {
        if (!ops[CRE_DESCOMP]->is_local() ) {
	        int fromproc = processorindex(trimap_2d(I, J, length));

	        //take the CDComp(I,I) where I=dot index and broadcast it
	        if( ops[CRE_DESCOMP]->has(I,J) ) {
	          if (fromproc != mpigetrank()) ops[CRE_DESCOMP]->add_local_indices(I,J);
	          
	          std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[CRE_DESCOMP]->get_element(I,J);
	          for (int iproc =0; iproc<oparray.size(); iproc++) {
	            //this only broadcasts the frame but no data
	            mpi::broadcast(calc, *oparray[iproc], fromproc);
	            
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
	            MPI_Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc, Calc);
	          }
	        }
	
	        //take the DDComp(I,I) where I=dot index and broadcast it
	        if( ops[DES_DESCOMP]->has(I,J) ) {
	          if (fromproc != mpigetrank()) ops[DES_DESCOMP]->add_local_indices(I,J);
	          
	          std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[DES_DESCOMP]->get_element(I,J);
	          for (int iproc =0; iproc<oparray.size(); iproc++) {
	            //this only broadcasts the frame but no data
	            mpi::broadcast(calc, *oparray[iproc], fromproc);
	            
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
	            MPI_Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc, Calc);
	          }
	        }
	
	        //take the DCComp(I,I) and CCComp(I,I) where I=dot index and broadcast it
	        if (has(DES) ) {
	          if( ops[CRE_CRECOMP]->has(I,J) ) {
	            if (fromproc != mpigetrank()) ops[CRE_CRECOMP]->add_local_indices(I,J);
	            
	            std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[CRE_CRECOMP]->get_element(I,J);
	            for (int iproc =0; iproc<oparray.size(); iproc++) {
	              //this only broadcasts the frame but no data
	              mpi::broadcast(calc, *oparray[iproc], fromproc);
	              
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
	              MPI_Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc, Calc);
	            }
	          }
	          if( ops[DES_CRECOMP]->has(I,J) ) {
	            if (fromproc != mpigetrank()) ops[DES_CRECOMP]->add_local_indices(I,J);
	            
	            std::vector<boost::shared_ptr<StackSparseMatrix> > oparray = ops[DES_CRECOMP]->get_element(I,J);
	            for (int iproc =0; iproc<oparray.size(); iproc++) {
	              //this only broadcasts the frame but no data
	              mpi::broadcast(calc, *oparray[iproc], fromproc);
	              
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
	              MPI_Bcast(oparray[iproc]->get_data(), oparray[iproc]->memoryUsed(), MPI_DOUBLE, fromproc, Calc);
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
	          calc.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	          if (other_proc_has_ops) {
	            ops[CRE_DESCOMP]->add_local_indices(I, J);
	            recvcompOps(*ops[CRE_DESCOMP], I, J, CRE_DESCOMP);
	          }
	          
	          other_proc_has_ops = true;
	          calc.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	          if (other_proc_has_ops) {
	            ops[DES_DESCOMP]->add_local_indices(I, J);
	            recvcompOps(*ops[DES_DESCOMP], I, J, DES_DESCOMP);
	          }
	          
	          if (has(DES)) {
	            other_proc_has_ops = true;
	            calc.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	            if (other_proc_has_ops) {
	      	      ops[CRE_CRECOMP]->add_local_indices(I, J);
	      	      recvcompOps(*ops[CRE_CRECOMP], I, J, CRE_CRECOMP);
	            }
	            other_proc_has_ops = true;
	            calc.recv(processorindex(trimap_2d(I, J, length)), 0, other_proc_has_ops);
	            if (other_proc_has_ops) {
	      	      ops[DES_CRECOMP]->add_local_indices(I, J);
	      	      recvcompOps(*ops[DES_CRECOMP], I, J, DES_CRECOMP);
	            }
	          }
	        } else {
	          //this will potentially send some ops
	          if (processorindex(trimap_2d(I, J, length)) == mpigetrank()) {
	            bool this_proc_has_ops = ops[CRE_DESCOMP]->has_local_index(I, J);
	            calc.send(processorindex(compsite), 0, this_proc_has_ops);
	            if (this_proc_has_ops) {
	      	      sendcompOps(*ops[CRE_DESCOMP], I, J, CRE_DESCOMP, compsite);
	            }
	            this_proc_has_ops = ops[DES_DESCOMP]->has_local_index(I, J);
	            calc.send(processorindex(compsite), 0, this_proc_has_ops);
	            if (this_proc_has_ops) {
	      	      sendcompOps(*ops[DES_DESCOMP], I, J, DES_DESCOMP, compsite);     
	            }
	            if (has(DES)) {
	      	      this_proc_has_ops = ops[CRE_CRECOMP]->has_local_index(I, J);
	      	      calc.send(processorindex(compsite), 0, this_proc_has_ops);
	      	      if (this_proc_has_ops) {
	      	        sendcompOps(*ops[CRE_CRECOMP], I, J, CRE_CRECOMP, compsite);     
	      	      }
	      	      this_proc_has_ops = ops[DES_CRECOMP]->has_local_index(I, J);
	      	      calc.send(processorindex(compsite), 0, this_proc_has_ops);
	      	      if (this_proc_has_ops) {
	      	        sendcompOps(*ops[DES_CRECOMP], I, J, DES_CRECOMP, compsite);     
	      	      }
	            }
	          } else continue;
	        }
	      }
      }
    }
  }
#endif
}


  //In loop block you never have any comp ops, but we need a few
  //to make CCDcomp, these operators are created in parallel and passed around
void StackSpinBlock::formTwoIndexOps() {

#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  //if dont have 2 index normal ops then return
  if (!has(CRE_DES)) return;

  int length = dmrginp.last_site();

  //add CRE_DESCOMP and DES_DESCOMP ops
  ops[CRE_DESCOMP] = make_new_stackop(CRE_DESCOMP, true);
  ops[DES_DESCOMP] = make_new_stackop(DES_DESCOMP, true);
  ops[CRE_DESCOMP]->set_local()=false;
  ops[DES_DESCOMP]->set_local()=false;

  if (has(DES)) {
    ops[DES_CRECOMP] = make_new_stackop(DES_CRECOMP, true);
    ops[CRE_CRECOMP] = make_new_stackop(CRE_CRECOMP, true);
    ops[DES_CRECOMP]->set_local()=false;
    ops[CRE_CRECOMP]->set_local()=false;
  }

  vector<int> dotindice;
  dotindice.push_back((sites[0] == 0) ? complementary_sites[0] : complementary_sites[complementary_sites.size()-1]);
  if (!dmrginp.spinAdapted()) { // when non-spinadapted, sites are spin orbitals
    dotindice.push_back((sites[0] == 0) ? complementary_sites[1] : complementary_sites[complementary_sites.size()-2]);    
  }


  const double screen_tol = dmrginp.twoindex_screen_tol();
  //build iterators
  std::vector<int> cinds;
  std::map<std::tuple<int, int, int>, int> tuple;
  std::vector<std::pair<int, int> > cdpair, ddpair;
  for (int idx = 0; idx < dotindice.size(); ++idx) 
    for (int i=0; i<complementary_sites.size(); i++) {
      int dotopindex = dotindice[idx];
      int compsite = complementary_sites[i];
      if (find(dotindice.begin(), dotindice.end(), compsite) != dotindice.end() && dotopindex > compsite) {
	continue;
      }
      int I = (compsite > dotopindex) ? compsite : dotopindex;
      int J = (compsite > dotopindex) ? dotopindex : compsite;
      if (dmrginp.use_partial_two_integrals() || screen_cd_interaction(I, J, sites, *get_twoInt(), screen_tol))
	cdpair.push_back(std::pair<int, int>(I, J));

      if (dmrginp.use_partial_two_integrals() || screen_dd_interaction(I, J, sites, *get_twoInt(), screen_tol))
	ddpair.push_back(std::pair<int, int>(I, J));
    }
  ops[CRE_DESCOMP]->build_iterators(*this, cinds, cdpair, tuple);
  ops[DES_DESCOMP]->build_iterators(*this, cinds, ddpair, tuple);
  if (has(DES)) {
    ops[DES_CRECOMP]->build_iterators(*this, cinds, cdpair, tuple);
    ops[CRE_CRECOMP]->build_iterators(*this, cinds, ddpair, tuple);
  }
  
  for (int idx = 0; idx < dotindice.size(); ++idx) 
  for (int i=0; i<complementary_sites.size(); i++) {
    int dotopindex = dotindice[idx];
    int compsite = complementary_sites[i];
    if (find(dotindice.begin(), dotindice.end(), compsite) != dotindice.end() && dotopindex > compsite) {
      continue;
    }
    int I = (compsite > dotopindex) ? compsite : dotopindex;
    int J = (compsite > dotopindex) ? dotopindex : compsite;
    if (dmrginp.use_partial_two_integrals() || screen_cd_interaction(I, J, sites, *get_twoInt(), screen_tol)) {
      
      //I and J are both dot indices
      if (find (dotindice.begin(), dotindice.end(), I) != dotindice.end() &&
	  find (dotindice.begin(), dotindice.end(), J) != dotindice.end()) {
	if (!ops[CRE_DESCOMP]->has_local_index(I,J)) 
	  ops[CRE_DESCOMP]->add_local_indices(I,J);
	if (has(DES)) {
	  if (!ops[DES_CRECOMP]->has_local_index(I,J))
	    ops[DES_CRECOMP]->add_local_indices(I,J);
	}
      }
      else {
	//the CCD_comp(compsite) needs CD_comp(compsite, dotsite)
	if (processorindex(compsite) == mpigetrank() && !ops[CRE_DESCOMP]->has_local_index(I,J)) 
	  ops[CRE_DESCOMP]->add_local_indices(I,J);
	if (processorindex(compsite) != mpigetrank() && ops[CRE_DESCOMP]->has_local_index(I,J)) 
	  ops[CRE_DESCOMP]->remove_local_indices(I,J);

	if (has(DES)) {
	  if (processorindex(compsite) == mpigetrank() && !ops[DES_CRECOMP]->has_local_index(I,J)) 
	    ops[DES_CRECOMP]->add_local_indices(I,J);
	  if (processorindex(compsite) != mpigetrank() && ops[DES_CRECOMP]->has_local_index(I,J)) 
	    ops[DES_CRECOMP]->remove_local_indices(I,J);
	}

      }
    }
    
    if (dmrginp.use_partial_two_integrals() || screen_dd_interaction(I, J, sites, *get_twoInt(), screen_tol)) {
      
      //I and J are both dot indices
      if (find (dotindice.begin(), dotindice.end(), I) != dotindice.end() &&
	  find (dotindice.begin(), dotindice.end(), J) != dotindice.end()) {
	if (!ops[DES_DESCOMP]->has_local_index(I,J)) 
	  ops[DES_DESCOMP]->add_local_indices(I,J);
	if (has(DES)) { 
	  if (!ops[CRE_CRECOMP]->has_local_index(I,J))
	    ops[CRE_CRECOMP]->add_local_indices(I,J);
	}
      }
      else {
	if (processorindex(compsite) == mpigetrank() && !ops[DES_DESCOMP]->has_local_index(I,J)) 
	  ops[DES_DESCOMP]->add_local_indices(I,J);
	if (processorindex(compsite) != mpigetrank() && ops[DES_DESCOMP]->has_local_index(I,J)) 
	  ops[DES_DESCOMP]->remove_local_indices(I,J);

	if (has(DES)) {
	  if (processorindex(compsite) == mpigetrank() && !ops[CRE_CRECOMP]->has_local_index(I,J)) 
	    ops[CRE_CRECOMP]->add_local_indices(I,J);
	  if (processorindex(compsite) != mpigetrank() && ops[CRE_CRECOMP]->has_local_index(I,J)) 
	    ops[CRE_CRECOMP]->remove_local_indices(I,J);
	}
      }
    }

  }

  vector<opTypes> opTypevec;
  opTypevec.push_back(CRE_DESCOMP);
  opTypevec.push_back(DES_DESCOMP);
  if (has(DES)) {
    opTypevec.push_back(DES_CRECOMP);
    opTypevec.push_back(CRE_CRECOMP);
  }
  for (int optypeindex=0; optypeindex<opTypevec.size(); optypeindex++) {
    opTypes optype = opTypevec[optypeindex];

#ifndef SERIAL
    for (int proc=0; proc<calc.size(); proc++) {
#else
    {
      int proc = 0;
#endif

      //for each processor loop over all its CDcomp operators and place them in ops_On_proc
      int arraysize = 0;
      std::vector<boost::shared_ptr<StackSparseMatrix> > ops_On_proc;
      if (mpigetrank() == proc) {
	for (int i=0; i<ops[optype]->get_size(); i++) {
	  int vecsize = ops[optype]->get_local_element(i).size();
	  for (int j=0; j< vecsize; j++) 
	    ops_On_proc.push_back(ops[optype]->get_local_element(i)[j]);
	}
	arraysize = ops_On_proc.size();
      }

#ifndef SERIAL
      //broadcast this ops_On_proc, so each processor has the shell of the operators
      mpi::broadcast(calc, arraysize, proc);
      if (mpigetrank() != proc) {
	for (int i=0; i<arraysize; i++) {
	  if (optype == CRE_DESCOMP)
	    ops_On_proc.push_back( boost::shared_ptr<StackCreDesComp>(new StackCreDesComp));
	  else if (optype == DES_DESCOMP)
	    ops_On_proc.push_back( boost::shared_ptr<StackDesDesComp>(new StackDesDesComp));
	  else if (optype == DES_CRECOMP)
	    ops_On_proc.push_back( boost::shared_ptr<StackDesCreComp>(new StackDesCreComp));
	  else if (optype == CRE_CRECOMP)
	    ops_On_proc.push_back( boost::shared_ptr<StackCreCreComp>(new StackCreCreComp));
	}
      }
#endif
      
      for (int i=0; i<ops_On_proc.size(); i++) {
#ifndef SERIAL
	mpi::broadcast(calc, *ops_On_proc[i], proc);
#endif
	if (ops_On_proc[i]->memoryUsed() == 0)
	  ops_On_proc[i]->allocate(braStateInfo, ketStateInfo);
	else
	  ops_On_proc[i]->Clear();
      }

      //make these comp operators from normal operators using multithreading
      //SplitStackmem();
#pragma omp parallel for schedule(dynamic)
      for (int i=0; i<ops_On_proc.size(); i++) {
	if (optype == CRE_DESCOMP)
	  ((StackCreDesComp&)(*ops_On_proc[i])).buildfromCreDes(*this);
	else if (optype == DES_DESCOMP)
	  ((StackDesDesComp&)(*ops_On_proc[i])).buildfromDesDes(*this);
	else if (optype == DES_CRECOMP)
	  ((StackDesCreComp&)(*ops_On_proc[i])).buildfromDesCre(*this);
	else if (optype == CRE_CRECOMP)
	  ((StackCreCreComp&)(*ops_On_proc[i])).buildfromCreCre(*this);
      }

      for (int i=0; i<ops_On_proc.size(); i++) {

#ifndef SERIAL
	MPI_Allreduce(MPI_IN_PLACE, ops_On_proc[i]->get_data(), ops_On_proc[i]->memoryUsed(), MPI_DOUBLE, MPI_SUM, Calc);
#endif
	
	if (mpigetrank() == proc) {
	  if (additionalMemory == 0)
	    additionaldata = ops_On_proc[i]->get_data();
	  additionalMemory += ops_On_proc[i]->memoryUsed();
	}
      }

      for (int i=ops_On_proc.size()-1; i>=0; i--) {
	if (mpigetrank() != proc)
	  ops_On_proc[i]->deallocate();
      }

      //MergeStackmem();
      //calc.barrier();
    }
  }
}


  //If Cre is not local then all procs get all copies of cre
  //and we also perform messagepasstwoindex ops or form two index ops
  //depending on whether two comp ops are already present
void StackSpinBlock::addAdditionalOps()
{
  dmrginp.datatransfer->start();

  addOneIndexOps();

  //we Already have Comp operators, now we just have to spread it around
  if ( has(CRE_DESCOMP))
    messagePassTwoIndexOps();
  else 
    formTwoIndexOps();
    
  dmrginp.datatransfer->stop();
}

  //mid way during the sweep one has to go from having normal ops to
  //comp ops and they are generated in thsi function
void StackSpinBlock::addAllCompOps() {

#ifndef SERIAL
  boost::mpi::communicator world;
#endif
  //if dont have 2 index normal ops then return
  if (!has(CRE_DES)) return;

  int length = dmrginp.last_site();

  //add CRE_DESCOMP and DES_DESCOMP ops
  ops[CRE_DESCOMP] = make_new_stackop(CRE_DESCOMP, true);
  ops[DES_DESCOMP] = make_new_stackop(DES_DESCOMP, true);
  ops[CRE_DESCOMP]->set_local()=false;
  ops[DES_DESCOMP]->set_local() = false;

  ops[CRE_DESCOMP]->build_iterators(*this, false);
  ops[DES_DESCOMP]->build_iterators(*this, false);

  if( has(DES)) {
    ops[DES_CRECOMP] = make_new_stackop(DES_CRECOMP, true);
    ops[CRE_CRECOMP] = make_new_stackop(CRE_CRECOMP, true);
    ops[DES_CRECOMP]->set_local()=false;
    ops[CRE_CRECOMP]->set_local() = false;
    
    ops[DES_CRECOMP]->build_iterators(*this, false);
    ops[CRE_CRECOMP]->build_iterators(*this, false);
  }

  vector<opTypes> opTypevec;
  opTypevec.push_back(CRE_DESCOMP);
  opTypevec.push_back(DES_DESCOMP);
  if (has(DES)) {
    opTypevec.push_back(CRE_CRECOMP);
    opTypevec.push_back(DES_CRECOMP);
  }
  // CDcomp
  //loop over all processors
  for (int optypeindex=0; optypeindex<opTypevec.size(); optypeindex++) {
    opTypes optype = opTypevec[optypeindex];

#ifndef SERIAL
    for (int proc=0; proc<calc.size(); proc++) {
#else
    {
      int proc = 0;
#endif

      //for each processor loop over all its CDcomp operators and place them in ops_On_proc
      int arraysize = 0;
      std::vector<boost::shared_ptr<StackSparseMatrix> > ops_On_proc;
      if (mpigetrank() == proc) {
	for (int i=0; i<ops[optype]->get_size(); i++) {
	  int vecsize = ops[optype]->get_local_element(i).size();
	  for (int j=0; j< vecsize; j++) 
	    ops_On_proc.push_back(ops[optype]->get_local_element(i)[j]);
	}
	arraysize = ops_On_proc.size();
      }

#ifndef SERIAL
      //broadcast this ops_On_proc, so each processor has the shell of the operators
      mpi::broadcast(calc, arraysize, proc);
      if (mpigetrank() != proc) 
	for (int i=0; i<arraysize; i++)
	  ops_On_proc.push_back( boost::shared_ptr<StackSparseMatrix>(new StackSparseMatrix));
      
      for (int i=0; i<ops_On_proc.size(); i++) {
	mpi::broadcast(calc, *ops_On_proc[i], proc);
	if (ops_On_proc[i]->memoryUsed() == 0)
	  ops_On_proc[i]->allocate(braStateInfo, ketStateInfo);
      }
#else
      for (int i=0; i<ops_On_proc.size(); i++) {
	if (ops_On_proc[i]->memoryUsed() == 0)
	  ops_On_proc[i]->allocate(braStateInfo, ketStateInfo);
      }
#endif

      //make these comp operators from normal operators using multithreading
      //SplitStackmem();
#pragma omp parallel for schedule(dynamic)
      for (int i=0; i<ops_On_proc.size(); i++) {
	if (optype == CRE_DESCOMP)
	  ((StackCreDesComp&)(*ops_On_proc[i])).buildfromCreDes(*this);
	else if (optype == DES_DESCOMP)
	  ((StackDesDesComp&)(*ops_On_proc[i])).buildfromDesDes(*this);
	else if (optype == DES_CRECOMP)
	  ((StackDesCreComp&)(*ops_On_proc[i])).buildfromDesCre(*this);
	else if (optype == CRE_CRECOMP)
	  ((StackCreCreComp&)(*ops_On_proc[i])).buildfromCreCre(*this);
      }

      for (int i=0; i<ops_On_proc.size(); i++) {
#ifndef SERIAL
	MPI_Allreduce(MPI_IN_PLACE, ops_On_proc[i]->get_data(), ops_On_proc[i]->memoryUsed(), MPI_DOUBLE, MPI_SUM, Calc);
#endif
	if (mpigetrank() == proc) {
	  if (additionalMemory == 0)
	    additionaldata = ops_On_proc[i]->get_data();
	  additionalMemory += ops_On_proc[i]->memoryUsed();
	}
      }

      for (int i=ops_On_proc.size()-1; i>=0; i--) {
	if (mpigetrank() != proc)
	  ops_On_proc[i]->deallocate();
      }
      //MergeStackmem();
    }
  }

  
}


}
