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

std::string StackSpinBlock::restore (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char* name)
{
  dmrginp.diskio->start();
  Timer disktimer;
  std::string file;

  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  
  p1out << "\t\t\t Restoring block file :: " << file << endl;

  std::ifstream ifs(file.c_str(), std::ios::binary);

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

  b.Load (ifs);
  ifs.close();

  

  dmrginp.diskio->stop();

  return file;
}
  
void StackSpinBlock::store (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char *name)
{
  dmrginp.diskio->start();
  Timer disktimer;
  std::string file;

  if (forward)
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-f-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  else
    file = str(boost::format("%s%s%d%s%d%s%d%s%d%s%d%s%d%s") % dmrginp.save_prefix() % "/Block-b-sites-"% sites[0] % "." % sites[sites.size()-1] % "-states" % left % "." % right % "-integral" %b.integralIndex % "rank" % mpigetrank() % ".tmp" );
  
  p1out << "\t\t\t Saving block file :: " << file << endl;
  
  
  std::ofstream ofs(file.c_str(), std::ios::binary);
  
  int lstate =  left;
  int rstate =  right;
  
  if (mpigetrank()==0) {
    StateInfo::store(forward, sites, b.braStateInfo, lstate);
    StateInfo::store(forward, sites, b.ketStateInfo, rstate);
  }


  b.Save (ofs);
  ofs.close(); 

  dmrginp.diskio->stop();

  //p1out << "\t\t\t block save disk time " << disktimer.elapsedwalltime() << " " << disktimer.elapsedcputime() << endl;
}

void StackSpinBlock::Save (std::ofstream &ofs)
{
  boost::archive::binary_oarchive save_block(ofs);
  save_block << *this;

  save_block << totalMemory;
  save_block << boost::serialization::make_array<double>(data, totalMemory);
  //for (int i=0; i<totalMemory; i++) 
  //save_block << data[i];

}

//helper function
void StackSpinBlock::Load (std::ifstream & ifs)
{

  boost::archive::binary_iarchive load_block(ifs);
  load_block >> *this;

  load_block >> totalMemory;
  data = Stackmem[omprank].allocate(totalMemory);
  load_block >> boost::serialization::make_array<double>(data, totalMemory);
  //for (int i=0; i<totalMemory; i++) 
  //load_block >> data[i];

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
}
