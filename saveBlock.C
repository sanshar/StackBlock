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

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted{
using namespace operatorfunctions;

void StackSpinBlock::make_iterator(StackSpinBlock& b, opTypes op, int* data, bool oneIndex, int numIndices) {

  b.ops[op] = make_new_stackop(op, true);
  std::vector<int> oneindex (numIndices, 0.0);
  std::vector<std::pair<int, int> > twoindex(numIndices, std::pair<int, int>(0,0));
  
  if (oneIndex) {
    for (int i=0; i<numIndices; i++) {
      oneindex[i] = data[i];
    }
  }
  else {
    for (int i=0; i<numIndices; i++) {
      twoindex[i].first = data[2*i];
      twoindex[i].second = data[2*i+1];
    }
  }

  b.ops[op]->build_iterators(b, oneindex, twoindex);
  
} 

void StackSpinBlock::restore (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char* name)
{
  dmrginp.diski->start();
  Timer disktimer;
  std::string file;

  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  
  p1out << "\t\t\t Restoring block file :: " << file << endl;


  int lstate =  left;
  int rstate =  right;

  if (mpigetrank() == 0) {
    StateInfo::restore(forward, sites, b.braStateInfo, lstate);
    StateInfo::restore(forward, sites, b.ketStateInfo, rstate);
  }
  
#ifndef SERIAL
  mpi::communicator world;
  mpi::broadcast(world, b.braStateInfo, 0);
  mpi::broadcast(world, b.ketStateInfo, 0);
#endif


  int* initialData = new int[23];
  int allindexsize;

  dmrginp.rawdatai->start();
  FILE *fp = fopen(file.c_str(), "rb");
  fread(initialData, sizeof(int), 23, fp);
  fread(&allindexsize, sizeof(int), 1, fp);
  int* allindices = new int[allindexsize];//large data
  fread(allindices, sizeof(int), allindexsize, fp);
  fread(&(b.totalMemory), sizeof(long), 1, fp);
  b.data = Stackmem[omprank].allocate(b.totalMemory);
  fread(b.data, sizeof(double), b.totalMemory, fp);
  fclose(fp);
  dmrginp.rawdatai->stop();

  //nowunpack the first 23 integers
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

  //this should be in the same order as opTypes are written in enumerator.h file
  //
  if (numham           != -1) { make_iterator(b, HAM,             &allindices[index], true,  numham);           index+=numham;} 
  if (numcre           != -1) { make_iterator(b, CRE,             &allindices[index], true,  numcre);           index+=numcre;}
  if (numcrecre        != -1) { make_iterator(b, CRE_CRE,         &allindices[index], false, numcrecre);        index+=2*numcrecre;}
  if (numdesdescomp    != -1) { make_iterator(b, DES_DESCOMP,     &allindices[index], false, numdesdescomp);    index+=2*numdesdescomp;}
  if (numcredes        != -1) { make_iterator(b, CRE_DES,         &allindices[index], false, numcredes);        index+=2*numcredes;}
  if (numcredescomp    != -1) { make_iterator(b, CRE_DESCOMP,     &allindices[index], false, numcredescomp);    index+=2*numcredescomp;}
  if (numcrecredescomp != -1) { make_iterator(b, CRE_CRE_DESCOMP, &allindices[index], true,  numcrecredescomp); index+=numcrecredescomp;}
  if (numdes           != -1) { make_iterator(b, DES,             &allindices[index], true,  numdes);           index+=numdes;}
  if (numdesdes        != -1) { make_iterator(b, DES_DES,         &allindices[index], false, numdesdes);        index+=2*numdesdes;}
  if (numcrecrecomp    != -1) { make_iterator(b, CRE_CRECOMP,     &allindices[index], false, numcrecrecomp);    index+=2*numcrecrecomp;}
  if (numdescre        != -1) { make_iterator(b, DES_CRE,         &allindices[index], false, numdescre);        index+=2*numdescre;}
  if (numdescrecomp    != -1) { make_iterator(b, DES_CRECOMP,     &allindices[index], false, numdescrecomp);    index+=2*numdescrecomp;}
  if (numcredesdescomp != -1) { make_iterator(b, CRE_DES_DESCOMP, &allindices[index], true,  numcredesdescomp); index+=numcredesdescomp;}
  if (numoverlap       != -1) { make_iterator(b, OVERLAP,         &allindices[index], true,  numoverlap);       index+=numoverlap;}
  dmrginp.readmakeiter->stop();


  dmrginp.rawdatai->start();
  dmrginp.rawdatai->stop();

  dmrginp.readallocatemem->start();
  double* localdata = b.data; 
  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = b.ops.begin(); it != b.ops.end(); ++it)
  {
    if(it->second->is_core()) {
      
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
  std::string file;

  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  
  p1out << "\t\t\t Saving block file :: " << file << endl;

  int lstate =  left;
  int rstate =  right;
  
  if (mpigetrank()==0) {
    StateInfo::store(forward, sites, b.braStateInfo, lstate);
    StateInfo::store(forward, sites, b.ketStateInfo, rstate);
  }


  int* initialData = new int[23];

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

  std::vector<int> allindices;
  allindices.insert(allindices.end(), b.sites.begin(), b.sites.end());
  allindices.insert(allindices.end(), b.complementary_sites.begin(), b.complementary_sites.end());

  for (std::map<opTypes, boost::shared_ptr< StackOp_component_base> >::iterator it = b.ops.begin(); it != b.ops.end(); ++it)
  {
    std::vector<int> indices = it->second->get_global_array();
    allindices.insert(allindices.end(), indices.begin(), indices.end());
  }

  dmrginp.rawdatao->start();
  FILE *fp = fopen(file.c_str(), "wb");
  int size = allindices.size();
  fwrite(initialData, sizeof(int), 23, fp);
  fwrite(&size, sizeof(int), 1, fp);
  fwrite(&allindices[0], sizeof(int), allindices.size(), fp);
  fwrite(&b.totalMemory, sizeof(long), 1, fp);
  fwrite(b.data, sizeof(double), b.totalMemory, fp);
  fclose(fp);
  dmrginp.rawdatao->stop();

  delete [] initialData;
  dmrginp.disko->stop();

}


}
