/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_STACKDENSITY_HEADER
#define SPIN_STACKDENSITY_HEADER
#include "StackOperators.h"
#include "StateInfo.h"

namespace SpinAdapted{
  class StackSpinBlock;
  class StackWavefunction;
class StackDensityMatrix : public SpinAdapted::StackSparseMatrix
{
public:
  void add_onedot_noise(StackWavefunction& wave_solutions, StackSpinBlock& big, bool act2siteops = true);
  void add_twodot_noise(const StackSpinBlock &big, const double noise);
  StackDensityMatrix()
  {
    conj = 'n';
    orbs = std::vector<int>();
    initialised = true;
    fermion = false;
    deltaQuantum.assign(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
  }
  StackDensityMatrix(const StateInfo& s) {
    conj = 'n';
    orbs = std::vector<int>();
    initialised = true;
    fermion = false;
    deltaQuantum.assign(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));   
    if (dmrginp.hamiltonian() == BCS) {
      int n_min = 100000, n_max = 0;
      for (int i = 0; i < s.quanta.size(); ++i) {
        if (s.quanta[i].get_n() > n_max)
          n_max = s.quanta[i].get_n();
        if (s.quanta[i].get_n() < n_min)
          n_min = s.quanta[i].get_n();
      }
      for (int dn=2; dn<=n_max-n_min; dn+=2) {
        deltaQuantum.push_back(SpinQuantum(dn, SpinSpace(0), IrrepSpace(0)));
        deltaQuantum.push_back(SpinQuantum(-dn, SpinSpace(0), IrrepSpace(0)));
      }
    }
  }
  void makedensitymatrix(std::vector<StackWavefunction>& wave_solutions, StackSpinBlock &big, const std::vector<double> &wave_weights,
			 const double noise, const double additional_noise, bool warmup);
  void makedensitymatrix(StackWavefunction& wave_solutions, StackSpinBlock &big, const double &wave_weight);
  StackDensityMatrix& operator+=(const StackDensityMatrix& other);

  void build(const StackSpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b=0) {return 0.0;}
};
}

#endif
