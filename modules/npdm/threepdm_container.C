/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "global.h"
#include <boost/format.hpp>
#include "threepdm_container.h"
#include "npdm_permutations.h"
#include <boost/filesystem.hpp>
#include <boost/range/algorithm.hpp>
#include <math.h>  
#include <boost/shared_ptr.hpp>
#include "IntegralMatrix.h"
#include "boostutils.h"

namespace SpinAdapted{
namespace Npdm{

//===========================================================================================================================================================

Threepdm_container::Threepdm_container( int sites )
{
  if ( dmrginp.store_spinpdm() ) {
    if(dmrginp.spinAdapted()) {
      double *data = Stackmem[omprank].allocate(pow(2*sites,6));
      threepdm.resize(2*sites,2*sites,2*sites,2*sites,2*sites,2*sites, data);
    }
    else {
      double* data = Stackmem[omprank].allocate(pow(sites,6));
      threepdm.resize(sites,sites,sites,sites,sites,sites, data);
    }
    threepdm.Clear();
  } 
  if ( !dmrginp.spatpdm_disk_dump() ) {
    if(dmrginp.spinAdapted()) {
      size_t s = sites;
      size_t len = s*s;
      double* data = Stackmem[omprank].allocate( (len*len*len + 3*len*len + 2*len)/6);
      pout << "allocating "<<len<<" doubles "<<endl;
      spatial_threepdm.resize(s, data);
    }
    else {
      size_t s = sites/2;
      size_t len = s*s;
      double* data = Stackmem[omprank].allocate((len*len*len + 3*len*len + 2*len)/6);
      spatial_threepdm.resize(s, data);
    }
    spatial_threepdm.Clear();
  }

  else{
    elements_stride_.resize(6);
    for(int i=0; i<6; i++ )
      elements_stride_[i]=pow(sites,5-i);

    char file[5000];
    sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",mpigetrank(),".tmp");
    //std::ofstream spatpdm_disk(file, std::ios::binary);
    spatpdm_disk=fopen(file,"wb");
    //32M buffer;
    setvbuf(spatpdm_disk,NULL,_IOFBF,1024*1024*32);
    //spatpdm_disk.open(file, std::ios::binary);
    //spatpdm_disk.close();
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdms(const int& i, const int& j, int integralIndex)
{
#ifndef SERIAL
  calc.barrier();
#endif
  Timer timer;
  if ( dmrginp.store_spinpdm() ) {
    accumulate_npdm();
    save_npdm_binary(i, j);
    save_npdm_text(i, j);
  }
  if ( !dmrginp.spatpdm_disk_dump() ) {
    accumulate_spatial_npdm();
    save_spatial_npdm_text(i, j, integralIndex);
  }
  save_spatial_npdm_binary(i, j);

#ifndef SERIAL
  calc.barrier();
#endif
  double cputime =timer.elapsedcputime(); 
  p3out << "3PDM save full array time " << timer.elapsedwalltime() << " " << cputime << endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_text(const int &i, const int &j)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << threepdm.dim1() << endl;

    double trace = 0.0;
    for(int i=0; i<threepdm.dim1(); ++i)
      for(int j=0; j<threepdm.dim2(); ++j)
        for(int k=0; k<threepdm.dim3(); ++k)
          for(int l=0; l<threepdm.dim4(); ++l)
            for(int m=0; m<threepdm.dim5(); ++m)
              for(int n=0; n<threepdm.dim6(); ++n) {
                if ( abs(threepdm(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += threepdm(i,j,k,l,m,n);
                }
              }
    ofs.close();
    pout << "Spin-orbital 3PDM trace = " << trace << "\n";
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_text(const int &i, const int &j, int integralIndex)
{
  if( mpigetrank() == 0)
  {
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.", i, j,".txt");
    ofstream ofs(file);
    ofs << spatial_threepdm.dim1() << endl;

    double energy1 = 0.0, energy2=0.0;
    double trace = 0.0;
    double nelec = dmrginp.total_particle_number();
    size_t dim1 = spatial_threepdm.dim1();
    size_t dim2 = dim1*dim1, dim3=dim2*dim1;
    double* twordm = Stackmem[omprank].allocate(dim1*dim1*dim1*dim1);
    double* onerdm = Stackmem[omprank].allocate(dim1*dim1);
    memset(twordm, 0, sizeof(double)*dim3*dim1);
    memset(onerdm, 0, sizeof(double)*dim2);

    for(int i=0; i<spatial_threepdm.dim1(); ++i)
      for(int j=0; j<spatial_threepdm.dim2(); ++j)
        for(int k=0; k<spatial_threepdm.dim3(); ++k) 
          for(int l=0; l<spatial_threepdm.dim4(); ++l)
            for(int m=0; m<spatial_threepdm.dim5(); ++m)
              for(int n=0; n<spatial_threepdm.dim6(); ++n) {
                if ( abs(spatial_threepdm(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
                  ofs << boost::format("%d %d %d %d %d %d %20.14e\n") % i % j % k % l % m % n % spatial_threepdm(i,j,k,l,m,n);
                  if ( (i==n) && (j==m) && (k==l) ) trace += spatial_threepdm(i,j,k,l,m,n);
                }
              }

    ofs.close();

    pout << "Spatial      3PDM trace  = " << trace << "\n";
    
    Stackmem[omprank].deallocate(onerdm, dim2);
    Stackmem[omprank].deallocate(twordm, dim3*dim1);
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_npdm_binary(const int &I, const int &J)
{
  if( mpigetrank() == 0)
  {
    
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm.", I, J,".bin");
    std::ofstream ofs(file, std::ios::binary);
    ofs.write((char*)(threepdm.data), sizeof(double)*threepdm.get_size());
    ofs.close();
  }
  
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::external_sort_index(const int &i, const int &j)
{

  boost::sort(nonspin_batch);
#ifndef SERIAL
  boost::mpi::communicator world;
  if(mpigetrank() != 0){
      calc.send(0,0, nonspin_batch);
  }
  else{
    //Store index from different processors on disk of root node.
    for(int p=0; p< calc.size();p++){
      if(p!=0){
        calc.recv(p,0, nonspin_batch);
      }
      char file[5000];
      //batch_index tmpbuffer[1000000];
      //std::copy(nonspin_batch.begin(),nonspin_batch.end(),tmpbuffer);
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.", i, j,p,".bin");
      FILE* inputfile=fopen(file,"wb");
      fwrite(&nonspin_batch[0],sizeof(Sortpdm::batch_index),nonspin_batch.size(),inputfile);
      fclose(inputfile);
      nonspin_batch.clear();
      if(calc.size() == 0) return;
    }
    //external sort nonspin_batch
    //TODO
    //It is not parallel.
    char outfilename[5000];
    sprintf (outfilename, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.",i,j,".bin");
    FILE* outputfile = fopen(outfilename,"wb");
    long sorting_buff= 1024*1024*(32/calc.size());
    //For batch_index, the sorting buff is about 96M/calc.size();
    std::vector<Sortpdm::cache<Sortpdm::batch_index>> filecache;
    for(int p=0; p< calc.size();p++){
      char file[5000];
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.", i, j,p,".bin");
      Sortpdm::cache<Sortpdm::batch_index> tmpcache( file, sorting_buff);
      filecache.push_back(tmpcache);
    }
    long outputbuffsize=sorting_buff*4;
    long outputbuff_position = 0;
    Sortpdm::batch_index* outputbuff = new Sortpdm::batch_index[outputbuffsize];
    for(;;){
      // select the smallest one in the current positions of different caches.
      int smallest = 0;
      for(int i=1 ; i< filecache.size(); i++){
        if(filecache[i].value() < filecache[smallest].value())
          smallest =i;
      }
      outputbuff[outputbuff_position++]=filecache[smallest].value();

      if(!filecache[smallest].forward()){
        filecache[smallest].clear();
        filecache.erase(filecache.begin()+smallest);
        if (filecache.size()==0) {
        fwrite(outputbuff,sizeof(Sortpdm::batch_index),outputbuff_position,outputfile);
        break;
        }
      }

      if(outputbuff_position == outputbuffsize){
        fwrite(outputbuff,sizeof(Sortpdm::batch_index),outputbuffsize,outputfile);
        outputbuff_position=0;
      }

      }
      fclose(outputfile);
      //Finish external sort of index.

      //Clean up.
      for(int p=0; p< calc.size();p++){
        char file[5000];
        sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm_index.", i, j,p,".bin");
        boost::filesystem::remove(file);
      }
    }
#else
    char file[5000];
    sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/threepdm_index.", i, j,".bin");
    FILE* outfile=fopen(file,"wb");
    fwrite(&nonspin_batch[0],sizeof(Sortpdm::batch_index),nonspin_batch.size(),outfile);
    nonspin_batch.clear();
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::save_spatial_npdm_binary(const int &I, const int &J)
{ 
  if(!dmrginp.spatpdm_disk_dump())
  {
    if( mpigetrank() == 0)
    {
      int integralIndex = 0;
      double energy1 = 0.0, energy2=0.0;
      double trace = 0.0;
      double nelec = dmrginp.total_particle_number();
      size_t dim1 = spatial_threepdm.dim1();
      size_t dim2 = dim1*dim1, dim3=dim2*dim1;
      double* twordm = Stackmem[omprank].allocate(dim1*dim1*dim1*dim1);
      double* onerdm = Stackmem[omprank].allocate(dim1*dim1);
      memset(twordm, 0, sizeof(double)*dim3*dim1);
      memset(onerdm, 0, sizeof(double)*dim2);
      
      for(int i=0; i<spatial_threepdm.dim1(); ++i)
	for(int j=0; j<spatial_threepdm.dim2(); ++j)
	  for(int k=0; k<spatial_threepdm.dim3(); ++k) 
	    for(int l=0; l<spatial_threepdm.dim4(); ++l)
	      for(int m=0; m<spatial_threepdm.dim5(); ++m)
		for(int n=0; n<spatial_threepdm.dim6(); ++n) {
		  if ( abs(spatial_threepdm(i,j,k,l,m,n)) > NUMERICAL_ZERO ) {
		    if (i==n && j==m && k==l) trace  += spatial_threepdm(i,j,k, l,m,n);
		    spatial_threepdm(k,j,i, n,m,l) = spatial_threepdm(i,j,k, l,m,n);
		    spatial_threepdm(i,k,j, m,l,n) = spatial_threepdm(i,j,k, l,m,n);
		    spatial_threepdm(j,i,k, l,n,m) = spatial_threepdm(i,j,k, l,m,n);
		    
		    spatial_threepdm(k,i,j, m,n,l) = spatial_threepdm(i,j,k, l,m,n);
		    spatial_threepdm(j,k,i, n,l,m) = spatial_threepdm(i,j,k, l,m,n);
		    spatial_threepdm(i,j,k, l,m,n) = spatial_threepdm(i,j,k, l,m,n);
		  }
		}
      
      
      for(int i=0; i<spatial_threepdm.dim1(); ++i)
	for(int j=0; j<spatial_threepdm.dim2(); ++j)
	  for(int k=0; k<spatial_threepdm.dim3(); ++k) 
	    for(int l=0; l<spatial_threepdm.dim4(); ++l)
	      for(int m=0; m<spatial_threepdm.dim5(); ++m) {
		twordm[i*dim3+j*dim2+k*dim1+l] += spatial_threepdm(i,j,m,m,k,l)/(nelec-2.);
	      }
      
      const std::vector<int>& ro = dmrginp.reorder_vector();
      for (int i=0; i<dim1; i++)
	for (int j=0; j<dim1; j++)
	  for (int k=0; k<dim1; k++) {
	    onerdm[i*dim1+j] += twordm[i*dim3+k*dim2+k*dim1+j]/(nelec-1.);
	    for (int l=0; l<dim1; l++) {
	      energy2 += v_2[integralIndex](2*i,2*j,2*l,2*k)*twordm[ro.at(i)*dim3+ro.at(j)*dim2+ro.at(k)*dim1+ro.at(l)]*0.5;
	    }
	  }
      for (int i=0; i<dim1; i++)
	for (int j=0; j<dim1; j++)
	  energy1 += v_1[integralIndex](2*i,2*j)*onerdm[ro.at(i)*dim1+ro.at(j)];
      
      pout << "Spatial      3PDM trace  = " << trace << "\n";
      pout << "Spatial      3PDM Energy = " << energy1+energy2+coreEnergy[integralIndex] << "\n";
      
      Stackmem[omprank].deallocate(onerdm, dim2);
      Stackmem[omprank].deallocate(twordm, dim3*dim1);
      
      char file[5000];
      sprintf (file, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",I,J,".bin");
#ifndef MOLCAS
      std::ofstream ofs(file, std::ios::binary);
      //boost::archive::binary_oarchive save(ofs);
      ofs.write( (char*)(spatial_threepdm.data), sizeof(double)*spatial_threepdm.get_size());
      ofs.close();
#else
      FILE* f = fopen(file,"wb");
      fwrite( (char*)(spatial_threepdm.data),sizeof(double),spatial_threepdm.get_size(),f);
      fclose(f);
#endif
    }
  }
  else{
    fclose(spatpdm_disk);
    //When spatpdm_disk is opened, the state numbers, i and j, are not known. 
    char oldfile[5000];
    sprintf (oldfile, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",mpigetrank(),".tmp");
    char newfile[5000];
    sprintf (newfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",I,J,mpigetrank(),".bin");
    boost::filesystem::rename(oldfile,newfile);
    if(dmrginp.pdm_unsorted())
      external_sort_index(I,J);
    
    else{
      char file[5000];
      sprintf (file, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",I,J,mpigetrank(),".bin");
      char finalfile[5000];
      sprintf (finalfile, "%s%s%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",I,J,".bin");
#ifndef SERIAL
      boost::mpi::communicator world;
      calc.barrier();
      Timer timer1;
      char tmpfile[5000];
      sprintf (tmpfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",I,J,mpigetrank(),".tmp");
      char sortedfile[5000];
      sprintf (sortedfile, "%s%s%d.%d.%d%s", dmrginp.save_prefix().c_str(),"/spatial_threepdm.",I,J,mpigetrank(),".bin");
      Sortpdm::partition_data<Sortpdm::index_element>((long)pow(dmrginp.last_site(),6),file,tmpfile);
      //TODO
      //tmpfile and sortedfile can be the same file. Because tmpfile is divided into many small files. It can be overwritten by sortedfile.
      //However, when they are different, it is a little faster. Maybe compiler can do some optimizations. 
      Sortpdm::externalsort<Sortpdm::index_element>(tmpfile,sortedfile,(long)pow(dmrginp.last_site(),6));
      calc.barrier();
      double cputime = timer1.elapsedcputime();
      p3out << "3PDM parallel external sort time " << timer1.elapsedwalltime() << " " << cputime << endl;
      Timer timer;
      Sortpdm::mergefile(sortedfile);
      calc.barrier();
      if(mpigetrank()==0) boost::filesystem::rename(sortedfile,finalfile);
      boost::filesystem::remove(tmpfile);
      cputime = timer.elapsedcputime();
      p3out << "3PDM merge sorted file time " << timer.elapsedwalltime() << " " << cputime << endl;
#else
      Timer timer2;
      Sortpdm::externalsort<Sortpdm::index_element>(file,finalfile,(long)pow(dmrginp.last_site(),6));
      boost::filesystem::remove(file);
      double cputime = timer2.elapsedcputime();
      p3out << "3PDM external sort time " << timer2.elapsedwalltime() << " " << cputime << endl;
#endif
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::load_npdm_binary(const int &i, const int &j) { abort(); }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::accumulate_npdm()
{
#ifndef SERIAL
  array_6d<double> tmp_recv;
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<calc.size(); ++p) {
      double *pProcData = Stackmem[omprank].allocate(threepdm.get_size());
      MPI_Recv(pProcData, threepdm.get_size(), MPI_DOUBLE, p, p,  Calc, MPI_STATUS_IGNORE);

      for (int i=0; i<threepdm.get_size(); i++) 
	if (abs(pProcData[i]) > NUMERICAL_ZERO) 
	  threepdm.data[i] = pProcData[i];
      Stackmem[omprank].deallocate(pProcData, threepdm.get_size());
    }
  }
  else 
  {
    MPI_Send(spatial_threepdm.data, spatial_threepdm.get_size(), MPI_DOUBLE, 0, mpigetrank(),  Calc);
    //calc.send(0, mpigetrank(), spatial_threepdm);
  }
#endif
} 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::accumulate_spatial_npdm()
{
#ifndef SERIAL
  mpi::communicator world;
  if( mpigetrank() == 0)
  {
    for(int p=1; p<calc.size(); ++p) {
      double *pProcData = Stackmem[omprank].allocate(spatial_threepdm.get_size());
      size_t datadim = spatial_threepdm.get_size();
      size_t maxint = 26843540;//mpi cannot transfer more than these number of doubles
      size_t maxiter = datadim/maxint;
      for (int i=0; i<maxiter; i++)
	MPI_Recv(pProcData+i*maxint, maxint, MPI_DOUBLE, p, p*100000+i,  Calc, MPI_STATUS_IGNORE);
      MPI_Recv(pProcData+maxiter*maxint, datadim-maxiter*maxint, MPI_DOUBLE, p, p*100000+maxiter,  Calc, MPI_STATUS_IGNORE);

      for (int i=0; i<spatial_threepdm.get_size(); i++) 
	if (abs(pProcData[i]) > NUMERICAL_ZERO) 
	  spatial_threepdm.data[i] = pProcData[i];
      Stackmem[omprank].deallocate(pProcData, spatial_threepdm.get_size());
    }
  }
  else 
  {
    size_t datadim = spatial_threepdm.get_size();
    size_t maxint = 26843540;//mpi cannot transfer more than these number of doubles
    size_t maxiter = datadim/maxint;
    double* pProcData = spatial_threepdm.data;
    for (int i=0; i<maxiter; i++)
      MPI_Send(pProcData+i*maxint, maxint, MPI_DOUBLE, 0, mpigetrank()*100000+i,  Calc);
    MPI_Send(pProcData+maxiter*maxint, datadim-maxiter*maxint, MPI_DOUBLE, 0, mpigetrank()*100000+maxiter,  Calc);
      //MPI_Send(spatial_threepdm.data, spatial_threepdm.get_size(), MPI_DOUBLE, 0, mpigetrank(),  Calc);
    //calc.send(0, mpigetrank(), spatial_threepdm);
  }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void Threepdm_container::update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{

  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    double val = it->second;
    if ( abs(val) < NUMERICAL_ZERO ) continue;

    assert( (it->first).size() == 6 );
    int i0 = (it->first)[0];
    int j0 = (it->first)[1];
    int k0 = (it->first)[2];
    int l0 = (it->first)[3];
    int m0 = (it->first)[4];
    int n0 = (it->first)[5];
    int i = ro.at(i0/2)*2 + i0%2;
    int j = ro.at(j0/2)*2 + j0%2;
    int k = ro.at(k0/2)*2 + k0%2;
    int l = ro.at(l0/2)*2 + l0%2;
    int m = ro.at(m0/2)*2 + m0%2;
    int n = ro.at(n0/2)*2 + n0%2;

    //if ( abs(val) > 1e-8 ) {
    //  pout << "so-threepdm val: i,j,k,l,m,n = " 
    //       << i << "," << j << "," << k << "," << l << "," << m << "," << n
    //       << "\t\t" << val << endl;
    //}

    // Test for duplicates
    threepdm(i,j,k,l,m,n) = val;
  }

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// This routine assumes that no spin-orbital indices are generated more than once

void Threepdm_container::update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch )
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      // Spin indices
      int i = (it->first)[0];
      int j = (it->first)[1];
      int k = (it->first)[2];
      int l = (it->first)[3];
      int m = (it->first)[4];
      int n = (it->first)[5];

//      if ( i%2 != n%2 ) continue;
//      if ( j%2 != m%2 ) continue;
//      if ( k%2 != l%2 ) continue;

//      spatial_threepdm( ro.at(i/2), ro.at(j/2), ro.at(k/2), ro.at(l/2), ro.at(m/2), ro.at(n/2) ) += it->second;
      spatial_threepdm( ro.at(i), ro.at(j), ro.at(k), ro.at(l), ro.at(m), ro.at(n) ) = it->second;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

long Threepdm_container::oneindex_spin(const std::vector<int> & orbital_element_index)
{
  // Take into account orbital reordering
  const std::vector<int>& ro = dmrginp.reorder_vector();

  assert( orbital_element_index.size() == 6);
  long linearindex=0;
  for(int i=0; i< 6; i++)
    linearindex+=(ro.at(orbital_element_index[i]))*elements_stride_[i];
  return linearindex;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------


void Threepdm_container::dump_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch)
{
  long spatpdm_disk_position= ftell(spatpdm_disk);
  std::map < long, double>  index_and_elements;
  for (auto it = spin_batch.begin(); it != spin_batch.end(); ++it) {
    assert( (it->first).size() == 6 );

    // Store significant elements only
    if ( abs(it->second) > NUMERICAL_ZERO ) {
      long linearindex = oneindex_spin(it->first);
      index_and_elements.insert(std::pair<long,double>(linearindex,it->second));
    }
  }
  if(index_and_elements.size()==0) return;
  Sortpdm::index_element index_elements[72];
  assert(72>= index_and_elements.size());
  int i=0;
  for(auto it = index_and_elements.begin(); it!=index_and_elements.end();it++){
    index_elements[i].index=it->first;
    index_elements[i].element=it->second;
    i++;
  }
  fwrite(index_elements,sizeof(Sortpdm::index_element),index_and_elements.size(),spatpdm_disk);
  if(dmrginp.pdm_unsorted()){
  Sortpdm::batch_index onerecord(index_and_elements.begin()->first,spatpdm_disk_position/sizeof(Sortpdm::index_element),index_and_elements.size(),mpigetrank());
  nonspin_batch.push_back(onerecord);
  }


  //char file[5000];
  //sprintf (file, "%s%s%d%s", dmrginp.save_prefix().c_str(),"/testspatial_threepdm.",mpigetrank(),".bin");
  //std::ofstream spatpdm_disk(file, std::ios::ate| std::ios::binary);
  //spatpdm_disk.open(file, std::ios::binary);
  //boost::archive::binary_oarchive save(spatpdm_disk);
  //for(auto it= index_elements.begin(); it!=index_elements.end(); it++)
  //{
  //  //pout <<"element:  "<<it-> index<< "\t\t"<<it->element<<endl;
  //  save << *it;
  //  //spatpdm_disk << it->index;
  //  //spatpdm_disk << it->element;

  //}
  //spatpdm_disk.close();
  
#if 0
#ifndef SERIAL
  //mpi 
  boost::mpi::communicator world;
  std::vector<std::map < long, double>> indexelement_array;

  // dump std::map of index and elements into disk;
  if(mpigetrank()==0){
    boost::mpi::gather(calc, index_and_elements, indexelement_array, 0);
    std::map < long, double> elements;
    for(auto it = indexelement_array.begin(); it!=indexelement_array.end();it++)
      elements.insert(it->begin(),it->end());

    for(auto it = elements.begin(); it!=elements.end();it++)
    {
      spatialpdm_disk.seekp(it->first*sizeof(double),ios::beg);
        spatialpdm_disk << it->second;
    }
  }
  else
    boost::mpi::gather(calc, index_and_elements, 0);

#else
    for(auto it = index_and_elements.begin(); it!=index_and_elements.end();it++)
    {
     // spatialpdm_disk.seekp(it->first*sizeof(double),ios_base::beg);
      spatialpdm_disk.seekp((it->first)*sizeof(double));
      spatialpdm_disk << it->second;
    }
#endif
#endif
}


void Threepdm_container::store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements)
{
  Threepdm_permutations perm;
  if(dmrginp.spinAdapted())
  {

    std::vector< std::pair< std::vector<int>, double > > spatial_batch;
    perm.get_spatial_batch(new_spin_orbital_elements,spatial_batch);
    if(dmrginp.spatpdm_disk_dump() ){
      dump_to_disk(spatial_batch);
    }
    else update_full_spatial_array(spatial_batch);
    if( dmrginp.store_spinpdm())
    {
      std::vector< std::pair< std::vector<int>, double > > spin_batch;
      perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
      update_full_spin_array( spin_batch );
    }
  }
  else{
    std::vector< std::pair< std::vector<int>, double > > spin_batch;
    perm.process_new_elements( new_spin_orbital_elements, nonredundant_elements, spin_batch );
    update_full_spin_array( spin_batch );
  }
}

//===========================================================================================================================================================

}
}


