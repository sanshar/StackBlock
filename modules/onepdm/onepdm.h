/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef ONEPDM_HEADER_H
#define ONEPDM_HEADER_H
#include <vector>
#include <newmat.h>

namespace SpinAdapted{
  class StackWavefunction;
  class StackSpinBlock;

void compute_onepdm(std::vector<StackWavefunction>& solutions, const StackSpinBlock& system, const StackSpinBlock& systemDot, const StackSpinBlock& newSystem, const StackSpinBlock& newEnvironment, const StackSpinBlock& big, const int numprocs);

void compute_one_pdm_0_2(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_one_pdm_2_0_0(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_one_pdm_0_2_0(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_one_pdm_1_1_0(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_one_pdm_1_1(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);

/*
void compute_pair_0_2(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_pair_2_0_0(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_pair_0_2_0(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_pair_1_1_0(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
void compute_pair_1_1(StackWavefunction& wave1, StackWavefunction& wave2, const StackSpinBlock& big, Matrix& onepdm);
*/
void save_onepdm_spatial_text(const Matrix& onepdm, const int &i, const int &j);
void save_onepdm_spatial_binary(const Matrix& onepdm, const int &i, const int &j);
void save_onepdm_binary(const Matrix& onepdm, const int &i, const int &j);
void load_onepdm_binary(Matrix& onepdm, const int &i, const int &j);
void save_onepdm_text(const Matrix& onepdm, const int &i, const int &j);
void accumulate_onepdm(Matrix& onepdm);

void save_pairmat_binary(const Matrix& pairmat, const int &i, const int &j);
void load_pairmat_binary(Matrix& pairmat, const int &i, const int &j);
void save_pairmat_text(const Matrix& onepdm, const int &i, const int &j);

std::vector<int> distribute_procs(const int numprocs, const int numjobs);
}
#endif
