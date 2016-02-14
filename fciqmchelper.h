#ifndef FCIQMC_HELPER_HEADER_H
#define FCIQMC_HELPER_HEADER_H

#include "global.h"
#include <vector>
#include "boost/shared_ptr.hpp"
#include "StateInfo.h"
#include "StackOperators.h"
#include "Stackwavefunction.h"
#include "ObjectMatrix.h"

namespace SpinAdapted{



//the MPS is stored in the left canonical form
//LLLLL..LC 

class MPS{

 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & stateInfos & stateIndex;
    }

  std::vector< StateInfo> stateInfos;
  int stateIndex;
  //std::vector< std::vector<Matrix> > SiteTensors; //these are the L matrices
  //StackWavefunction w; //the last wavefunction

  void Init(std::vector<bool>& occ, int stateIndex);
 public:
  static int sweepIters;
  static bool spinAdapted;
  static std::vector<StackSpinBlock> siteBlocks;

  MPS() {};
  MPS(int stateindex); 
  MPS(std::vector<bool>& occ, int stateIndex);
  MPS(ulong* occnum, int length, int stateIndex);

  void ApplyCD(int i, int j);
  //std::vector<Matrix>& getSiteTensors(int i) {return SiteTensors[i];}
  //const std::vector<Matrix>& getSiteTensors(int i) const {return SiteTensors[i];}
  //const StackWavefunction& getw() const {return w;}
  //void scale(double r) {Scale(r, w);}
  //void normalize() {int success; w.Normalise(&success);}
  //double get_coefficient(const vector<bool>& occ_strings);
  //void writeToDiskForDMRG(int state, bool writeStateAverage=false);
  void CollapseToDeterminant(char* s);
};


 //statea is multiplied with Operator O|Mpsa> and then we compress it to get stateb
 //void compressOperatorTimesMPS(const MPS& statea, MPS& stateb);

 //calculate overlap between a and b <Mpsa|Mpsb>
 //double calculateOverlap (const MPS& a, const MPS& b);

 //calculate hamiltonian matrix between a and b <Mpsa|H|Mpsb>
 void Canonicalize(const bool& forward, int stateIndex);
 void calcHamiltonianAndOverlap(int statea, int stateb, double& h, double& o, bool sameStates=false) ;
 void AddMPS(int* state, double* scale, int nstates, int stateout);
 //void AddMPS(int statea, int stateb, int stateout);
}

#endif
