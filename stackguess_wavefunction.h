#ifndef SPIN_STACKGUESS_WAVEFUNCTION_HEADER
#define SPIN_STACKGUESS_WAVEFUNCTION_HEADER
#include "ObjectMatrix.h"
#include "rotationmat.h"
#include <vector>

namespace SpinAdapted{
  class StateInfo;
  class StackSpinBlock;
  class StackWavefunction;

namespace GuessWave
{
  void TransformLeftBlock(StackWavefunction& oldwavefunction, const StateInfo& newstateinfo, const std::vector<Matrix>& RotationMatrix, StackWavefunction& tempoldWave);
  void TransformRightBlock(const StackWavefunction& tempnewWave, const StateInfo& tempoldStateInfo, const std::vector<Matrix>& RotationMatrix, StackWavefunction& trial);
  void guess_wavefunctions(std::vector<StackWavefunction>& solution, DiagonalMatrix& e, const StackSpinBlock &big, 
			   const guessWaveTypes &guesswavetype, const bool &onedot,  const bool& transpose_guess_wave, int nroots, double additional_noise=0.0, int currentState=0);
  void guess_wavefunctions(StackWavefunction& solution, DiagonalMatrix& e, const StackSpinBlock &big,
			   const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave, double additional_noise=0.0);

  void guess_wavefunctions(StackWavefunction& solution, DiagonalMatrix& e, const StackSpinBlock &big,
			   const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave, double additional_noise, bool ket);

  //onedot transpose wave guess
  void transpose_previous_wavefunction(StackWavefunction& trial, const StateInfo& stateInfo, const std::vector<int>& rightsites, const std::vector<int> &dotsites, const int state, const bool &onedot, const bool& transpose_guess_wave);
  void transpose_previous_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave, bool ket=true);
  void onedot_transpose_wavefunction(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
                                      const StackWavefunction& guesswf, StackWavefunction& transposewf);
  void onedot_threeindex_to_twoindex_wavefunction(const StateInfo& twostateinfo, const ObjectMatrix3D< std::vector<Matrix> >& 
						  threewavefunction, StackWavefunction& twowavefunction, const StateInfo& guessstateinfo);
  void onedot_twoindex_to_threeindex_wavefunction(const StateInfo& stateinfo, const StackWavefunction& twowavefunction, 
						  ObjectMatrix3D< std::vector<Matrix> > & threewavefunction);

  void onedot_shufflesysdot(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
			    const StackWavefunction& guesswf, StackWavefunction& transposewf);

  void onedot_twoindex_to_threeindex_shufflesysdot(const StateInfo& stateinfo, const StackWavefunction& twowavefunction, 
						   ObjectMatrix3D< vector<Matrix> >& threewavefunction);
  void transform_previous_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state, const bool &onedt, 
				       const bool& transpose_guess_wave);

  void transform_previous_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state, const bool &onedt, 
				       const bool& transpose_guess_wave, bool ket);
  void transform_previous_wavefunction(StackWavefunction& trial, const StateInfo& stateInfo, const std::vector<int> &leftsites, const std::vector<int> &rightsites, const int state, const bool &onedot, const bool& transpose_guess_wave);

  //ondeot transform wave guess                                                                                  
  void onedot_transform_wavefunction(const StateInfo& oldstateinfo, const StateInfo& newstateinfo, const StackWavefunction& oldwavefunction,
				     const std::vector<Matrix>& inverseLeftRotationMatrix,
                                     const std::vector<Matrix>& rightRotationMatrix, StackWavefunction& newwavefunction, 
				     const bool& transpose_guess_wave, bool ket = true);
  void basic_guess_wavefunction(DiagonalMatrix& e, StackWavefunction& trial, const StateInfo *stateinfo, const int state);
  void transform_previous_twodot_to_onedot_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state);

}
}

#endif
