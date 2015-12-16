/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef INTEGRAL_ARRAY_HEADER
#define INTEGRAL_ARRAY_HEADER
#define WANT_STREAM
#include <iostream>
#include <fstream>
#include <cmath>
//#include "parser.h"
#include <newmat.h>
#include <newmatio.h>
#include <boost/filesystem.hpp>
#include <multiarray.h>
#include <newmatutils.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "blas_calls.h"
#include <boost/algorithm/string.hpp>
#include <string>
#include <algorithm>
#include "input.h"
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "global.h"

using namespace boost;

namespace SpinAdapted{
  
class OneElectronArray
  {
  private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & dim & rhf & dummyZero & bin & Occnum & isCalcLCC & isH0;
    }
    //shared_mem_vector rep;
    double* rep;
    double dummyZero; /* this variable holds the value 0, for use with rhf integrals */
  public:
    bool isCalcLCC; //if the calculation is of LCC type then we dont need separate integrals for H0 and V
    bool isH0; //if this is H0 then make sure that all integrals are zero if dn = 0
    long dim;
    bool rhf;
    bool bin; //if the file is binary
    vector<double> Occnum;
    
  OneElectronArray() : isCalcLCC(false), dim(0), dummyZero(0.0), rhf(false), bin(false), rep(0)
    {}
  OneElectronArray(bool pLCC, bool pisH0) : isCalcLCC(pLCC), isH0(pisH0), dim(0), dummyZero(0.0), rhf(false), bin(false), rep(0)
    {}
    /*
  OneElectronArray(int n, bool rhf_=false, bool bin_=false):dummyZero(0.0), rhf(rhf_), bin(bin_), rep(0)
    {
      ReSize(n);
    }
    */
    bool& set_isH0() {return isH0;}
    double*& set_data() {return rep;}
    double& operator()(int i, int j);
    double operator()(int i, int j) const;
    void ReSize(int n);
    long NOrbs() const {
      return dim;
    }


    void ReadOccNum(ifstream& inFile);
    
    void ReadFromDumpFile(ifstream& dumpFile, int norbs);

    void DumpToFile(ofstream& dumpFile) const;

  };


// ***************************************************************************
// class TwoElectronArray
// ***************************************************************************
/*!
@brief Class for real Hermitian two-electron quantities such as two-e integrals and two-particle r.d.m.

- Uses (12|1'2') format.		

- Hermitian symmetry (ij|kl) = (kl|ij) is AUTOMATICALLY enforced in the underlying storage.

- Additional permutational symmetry (appropriate for 2-e integrals, is
used by default e.g. (i(1)j(2)|k(1)l(2)) = (k(1)j(2)|i(1)l(2)) but
can be turned off through the flag permSymm (default is true).

- The convention for indexing and spins is: odd-index = alpha, even-index = beta.
e.g. twoe(0, 2, 0, 2) would be a all-beta type quantity.

- For spin-restricted quantities, the underlying storage (in rep) can be reduced
by setting the rhf flag. However, to the external user the indexing still behaves
as unrestricted. Default is false.
*/
class TwoElectronArray // 2e integral, notation (12|12), symmetric matrix
  {
  private:

    double* rep; //!< underlying storage enforcing real hermitian symmetry.

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & dim;
      ar & indexMap;
      ar & rhf;
      ar & bin;
      ar & permSymm;
      ar & matDim;
      ar & isCalcLCC;
      ar & isH0;
    }


  public:
    long dim;  //!< #spinorbs
    array_2d<int> indexMap;
    double dummyZero; // dummy variable to hold zero for RHF integrals
    bool isCalcLCC; //if the calculation is of LCC type then we dont need separate integrals for H0 and V
    bool isH0; //if this is H0 then make sure that all integrals are zero if dn = 0

    // default is unrestrictedPermSymm
    enum TwoEType { unrestrictedPermSymm, restrictedPermSymm,
                    unrestrictedNonPermSymm, restrictedNonPermSymm };
    bool rhf; //!< reduced storage for restricted quantities.
    bool permSymm; //!< enforce additional 2e- permutational symmetry.
    bool bin; //if binary file is to be read
    long matDim; 
  public:

  TwoElectronArray() : isCalcLCC(false), dim(0), rhf(false), permSymm(true), bin(false), dummyZero(0.0)
    {}

    TwoElectronArray(bool pLCC, bool pisH0) : isCalcLCC(pLCC), isH0(pisH0), dim(0), rhf(false), permSymm(true), bin(false), dummyZero(0.0)
    {}

    TwoElectronArray(TwoEType twoetype);
    
    explicit TwoElectronArray(long n, TwoEType twoetype)
    {
      *this = TwoElectronArray(twoetype);
      ReSize(n);
    }

    int NOrbs() const
      {
        return dim;
      }
    array_2d<int>& GetMap()
    {
      return indexMap;
    }

    double*& set_data() {return rep;}
    void ReSize(int n);

    int MapIndices(int n);

    virtual void Load(std::string prefix, int index) {}
    virtual void Save(std::string prefix, int index) {}
    virtual double& operator()(int i, int j, int k, int l);
    virtual double operator () (int i, int j, int k, int l) const;

    void ReadFromDumpFile(ifstream& dumpFile, int norbs);
    
    void DumpToFile(ofstream& dumpFile) const;

  };



class PartialTwoElectronArray : public TwoElectronArray // 2e integral, notation (OrbIndex2|12), where OrbIndex is given
{
  private:

    std::vector<std::vector<double> > rep; //!< underlying storage 
    std::vector<int> OrbIndex; //this is the orb i
    int dim;  //!< #spinorbs
    double dummyZero;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version)
    {
      ar & rep & dim & OrbIndex;
    }


  public:

    PartialTwoElectronArray(std::vector<int>& pOrbIndex) : OrbIndex(pOrbIndex)
    {}

    PartialTwoElectronArray(std::vector<int>& pOrbIndex, int pdim) : OrbIndex(pOrbIndex), dim(pdim)
    {}

    void populate(TwoElectronArray& v2);

    void setOrbIndex(std::vector<int>& pOrbIndex) {OrbIndex = pOrbIndex;}

    void Load(std::string prefix, int index=0);

    void Save(std::string prefix, int index=0);

    double operator () (int i, int j, int k, int l) const;
 
    double& operator () (int i, int j, int k, int l);

  };

// ****************************************************
// class PairArray
// ****************************************************
/*!
@brief Class for electron-pairing quantities such as pairing potential and pairing matrix in particle number nonconserving calculations. Only singlet pairing is included

- pairing part in Hamiltonian is given by \Delta_{ij}a_{ia}^\dagger a_{jb}^\dagger+c.c.

- pairing potential \Delta can be symmetric or nonsymmetric, depends on spin-restricted or unrestricted calculation

- default is unrestricted

*/


class PairArray {
  private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive &ar, const unsigned int version) {
      ar & dim & rhf & dummyZero & bin;
    }
    double* rep;
    double dummyZero;
  public:
    bool rhf;
    bool bin;
    long dim;

  public:
    PairArray(): dim(0), rhf(false), bin(false), dummyZero(0.0), rep(0) {}
    //PairArray(int n, bool rhf_=false, bool bin_=false): rhf(rhf_), bin(bin_), dummyZero(0.0), rep(0) {
    //  ReSize(n);
    //}

    double*& set_data() {return rep;}    
    double operator() (int i, int j) const;    
    double& operator() (int i, int j);
    
    void ReSize(int n);

    int NOrbs() const {
      return dim;
    }

    void ReadFromDumpFile(ifstream& dumpFile, int norbs) {
      cerr << "PairArray::ReadFromDumpFile not implemented yet!";
      abort();
    }
    
    void DumpToFile(ofstream& dumpFile) const {
      cerr << "PairArray::DumpToFile not implemented yet!";
      abort();
    }

  };

// ****************************************************
// class CCCCArray
// ****************************************************
/*!
@brief Class for CCCC type quantities

- The Hamiltonian is given by
    0.25*w_{ijkl}C_{ia}C_{ja}C_{kb}C_{lb} + c.c.
  where C means creation operator

- Antisymmetry is enforced:
    w_{ijkl}=-w_{jikl}=-w_{ijlk}=w_{jilk}
  therefore, the underlying storage has i>j, l>k, composed indice ij=i*(i-1)/2+j, lk=l*(l-1)/2+k, and the storage is a matrix with (ij, lk) as indice

- Spin symmetry is controled by rhf option. if rhf = true, additionaly
    w_{ijkl}=w_{lkji}
  the underlying storage bocomes a symmetric matrix wrt ij and lk
*/


class CCCCArray {
  private:

    double* rep;
    double dummyZero;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, const unsigned int version)
    {
      ar & dim;
      ar & indexMap;
      ar & rhf;
      ar & bin;
      ar & matDim;
    }

  public:
    long dim;
    array_2d<int> indexMap;

    bool rhf;
    bool bin;
    long matDim; 

    // constructors
    CCCCArray(): dim(0), rhf(false), bin(false), dummyZero(0.) {}
    
    explicit CCCCArray(int n, bool _rhf) : rhf(_rhf), bin(false), dummyZero(0.) {
      ReSize(n);
    }

    int NOrbs() const { return dim;}
    array_2d<int>& GetMap() {  return indexMap;}
    double*& set_data() {return rep;}
    void ReSize(int n);

    int MapIndices(int n);
    
    virtual double operator() (int i, int j, int k, int l) const;
    virtual void set(int i, int j, int k, int l, double value);
    

    void ReadFromDumpFile(ifstream& dumpFile, int norbs) {
      cerr << "CCCCArray::ReadFromDumpFile not implemented yet!";
      abort();
    }
    
    void DumpToFile(ofstream& dumpFile) const {
      cerr << "CCCCArray::DumpToFile not implemented yet!";
      abort();
    }
};

// ****************************************************
// class CCCDArray
// ****************************************************
/*!
@brief Class for CCCD type quantities

- The Hamiltonian is given by
    1/6*w_{ijkl}C_iC_jC_kD_l
  where C means creation operator and D means destruction operator

- Antisymmetry is enforced:
    w_{ijkl}=-w_{ikjl}=-w_{jikl}=w_{jkil}=-w_{kjil}=w_{kijl}
  Considering spin, we only store w_{ijkl} where i>j,k,l and spin of i,j,l the same, while k the opposite. Those which cannot be deduced from this type is 0. composed indice ij=i*(i-1)/2+j, kl=k*n+l and the storage is two matrices with (ij, kl) as indice

- Spin symmetry is controled by rhf option. if rhf = true, additionaly
    w_{ia,ja,kb,la}=-w_{ib,jb,ka,lb}
  the underlying storage bocomes only one matrix
*/


class CCCDArray {
  private:
    //Matrix repA, repB;
    double *rep;
    double dummyZero;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, const unsigned int version)
    {
      ar & dim;
      ar & indexMap;
      ar & rhf;
      ar & bin;
      ar & matDim;
    }
  
  public:
    long dim;
    array_2d<int> indexMap;    
    
    bool rhf;
    bool bin;
    long matDim;

    // constructors
    CCCDArray(): dim(0), rhf(false), bin(false), dummyZero(0.) {}
    
    explicit CCCDArray(int n, bool _rhf) : rhf(_rhf), bin(false), dummyZero(0.) {
      ReSize(n);
    }

    int NOrbs() const { return dim;}
    array_2d<int>& GetMap() {  return indexMap;}
    double*& set_data() {return rep;}
    void ReSize(int n);
    
    int MapIndices(int n);

    virtual double operator() (int i, int j, int k, int l) const;    
    virtual void set(int i, int j, int k, int l, double value);

    void ReadFromDumpFile(ifstream& dumpFile, int norbs) {
      cerr << "CCCDArray::ReadFromDumpFile not implemented yet!";
      abort();
    }
    
    void DumpToFile(ofstream& dumpFile) const {
      cerr << "CCCDArray::DumpToFile not implemented yet!";
      abort();
    }
    
};

} // namespace

#endif
