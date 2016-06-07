/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef THREEPDM_CONTAINER_H
#define THREEPDM_CONTAINER_H

#include "multiarray.h"
#include "npdm_container.h"
#include "externalsort.h"


//Internally the 3RDM is stored as (i,j,k,l,m,n) = E^{i,j,k}_{n,m,l}
namespace SpinAdapted{
namespace Npdm{

  using namespace Sortpdm;
//===========================================================================================================================================================

class Threepdm_container : public Npdm_container {

  public:
    Threepdm_container( int sites );
//FIXME destructor?


    void save_npdms(const int &i, const int &j, int integralIndex=0);
    void store_npdm_elements( const std::vector< std::pair< std::vector<int>, double > > & new_spin_orbital_elements );
    void clear() { threepdm.Clear(); spatial_threepdm.Clear(); nonredundant_elements.clear(); }

    array_6d_3rdm<double>& get_spatial_threepdm() { assert(dmrginp.store_spinpdm()); return spatial_threepdm; }
    virtual ~Threepdm_container() 
    {
      if(spatial_threepdm.get_size() != 0) 
	Stackmem[omprank].deallocate(spatial_threepdm.data, spatial_threepdm.get_size());

      if (threepdm.get_size() != 0)
	Stackmem[omprank].deallocate(threepdm.data, threepdm.get_size());
    }

  private:
    // Vector to store nonredundant spin-orbital elements only
    std::vector< std::pair< std::vector<int>, double > > nonredundant_elements;
    // Optional arrays to store the full spin and/or spatial PDMs in core if memory allows.
    array_6d<double> threepdm;
    array_6d_3rdm<double> spatial_threepdm;
    std::vector<int> elements_stride_;
    std::vector<batch_index> nonspin_batch;
    FILE* spatpdm_disk;
    FILE* batch_index_file;
    //long spatpdm_disk_position;

    void save_npdm_text(const int &i, const int &j);
    void save_npdm_binary(const int &i, const int &j);
    void save_spatial_npdm_text(const int &i, const int &j, int integralIndex=0);
    void save_spatial_npdm_binary(const int &i, const int &j);
    void load_npdm_binary(const int &i, const int &j);
    void accumulate_npdm();
    void accumulate_spatial_npdm();
    void external_sort_index(const int &i, const int &j);
  
    void update_full_spin_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    void update_full_spatial_array( std::vector< std::pair< std::vector<int>, double > >& spin_batch );
    long oneindex_spin(const std::vector<int> & orbital_element_index);
    void dump_to_disk(std::vector< std::pair< std::vector<int>, double > > & spin_batch);


    template<class U> friend class Npdm::Sortpdm::cache;
};

//===========================================================================================================================================================

}
}

#endif
