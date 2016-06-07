/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "MatrixBLAS.h"
#include "operatorfunctions.h"
#include "Stackspinblock.h"
#include "StateInfo.h"
#include "Stackwavefunction.h"
#include "stackguess_wavefunction.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
#include "pario.h"

namespace SpinAdapted{
void GuessWave::TransformLeftBlock(StackWavefunction& oldwavefunction, const StateInfo& newstateinfo, const std::vector<Matrix>& RotationMatrix, StackWavefunction& tempoldWave)
{  
  for (int a=0; a<tempoldWave.nrows(); a++)
    for (int b=0; b<tempoldWave.ncols(); b++)
    {      
      int olda = newstateinfo.leftStateInfo->leftStateInfo->newQuantaMap[a];
      if (oldwavefunction.allowed(olda, b)) {
	  StackMatrix& lM = oldwavefunction(olda, b);
	Matrix tM = RotationMatrix[olda];
	StackMatrix& nM = tempoldWave.operator_element(a, b);
	MatrixMultiply(tM, 'n', lM, 'n', nM, 1.0);
      }
    }
}

void GuessWave::TransformRightBlock(const StackWavefunction& tempnewWave, const StateInfo& oldStateInfo, const std::vector<Matrix>& RotationMatrix, StackWavefunction& trial)
{

  for (int a=0; a<tempnewWave.nrows(); a++)
    for (int b=0; b<tempnewWave.ncols(); b++)
    {      
      int transB = oldStateInfo.rightStateInfo->leftStateInfo->newQuantaMap[b];
      if (tempnewWave.allowed(a,b)) {
	StackMatrix& nM = trial.operator_element(a, transB);
	const StackMatrix& oM = tempnewWave.operator_element(a, b);

	Matrix rM = RotationMatrix[transB];
	//rM = rM.t();
	MatrixMultiply(oM, 'n', rM, 't', nM, 1.0);
      }
    }

}




void GuessWave::transpose_previous_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave, bool ket)
{
  StateInfo oldStateInfo;
  StackWavefunction oldWave, temptrial;  ///ADD THE data


  p2out << "\t\t\t Transposing previous wavefunction" << endl;
  if(!onedot)
  {
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_rightBlock()->get_sites(), state, true);

    temptrial.initialise(trial); temptrial.Clear();
    temptrial.copyData(trial);
    //copy(trial.get_operatorMatrix(), temptrial.get_operatorMatrix());

    if (oldWave.get_onedot()) {
      temptrial.deallocate();
      temptrial.initialise(oldWave.get_deltaQuantum(), *oldStateInfo.rightStateInfo, *oldStateInfo.leftStateInfo, false);
      temptrial.Clear();
    }
    for (int i = 0; i < temptrial.nrows(); ++i)
      for (int j = 0; j < temptrial.ncols(); ++j)
        if (temptrial.allowed(i, j))
        {
          assert(oldWave.allowed(j, i));
          assert(temptrial(i, j).Nrows() == oldWave(j, i).Ncols());
          assert(temptrial(i, j).Ncols() == oldWave(j, i).Nrows());
          const StateInfo* s = &big.get_stateInfo();
	  StackMatrix& tmp = oldWave(j, i);

	  for (int r = 0; r<tmp.Nrows(); r++)
	    for (int c = 0; c<tmp.Ncols(); c++)
	      temptrial(i,j)(c+1, r+1) = tmp(r+1, c+1);
	  int parity = getCommuteParity(oldStateInfo.rightStateInfo->quanta[i],
					oldStateInfo.leftStateInfo->quanta[j],
					oldWave.get_deltaQuantum(0));
	  if (parity == -1)
	    MatrixScale(-1.0, temptrial(i,j));

        }
    if (oldWave.get_onedot()) {
      //the previous wavefunction was onedot and now make this previous onedot into two dot
      //DCOPY(temptrial.memoryUsed(), const_cast<double*>(temptrial.get_data()), 1, oldWave.get_data(), 1);

      std::vector<Matrix> rightRotationMatrix;
      LoadRotationMatrix(big.get_leftBlock()->get_sites(), rightRotationMatrix, state);
      trial.initialise(oldWave.get_deltaQuantum(), *big.get_stateInfo().leftStateInfo, *big.get_stateInfo().rightStateInfo, 
		       trial.get_onedot(), trial.get_data(), trial.memoryUsed());
      trial.Clear();
      for (int i=0; i<temptrial.nrows(); i++)
      for (int j=0; j<temptrial.ncols(); j++) {
      if (temptrial.allowed(i,j)) {
	MatrixMultiply(rightRotationMatrix[oldStateInfo.rightStateInfo->newQuantaMap[i]], 'n', 
		       temptrial(i,j), 'n', trial(oldStateInfo.rightStateInfo->newQuantaMap[i], j), 1.0);
      }
      }
    }
    else 
      trial.copyData(temptrial);
    //copy(temptrial.get_operatorMatrix(), trial.get_operatorMatrix());

    oldStateInfo.Free();
    temptrial.deallocate();
    oldWave.deallocate();
  }
  else
  {
    //this is because we are doing transpose before the middle of the sweep (usually it is done at the very end)
    //this can happen when one is doing partial sweep
    if (big.get_leftBlock()->get_leftBlock() == 0) {
      vector<int> wfsites = big.get_rightBlock()->get_sites();
      wfsites.insert(wfsites.end(), big.get_rightBlock()->get_rightBlock()->get_sites().begin(), big.get_rightBlock()->get_rightBlock()->get_sites().end());
      sort(wfsites.begin(), wfsites.end());
      oldWave.LoadWavefunctionInfo(oldStateInfo, wfsites, state, true); 

      //the oldwav is  [s.] [e]
      //we want to convert it to [e][s.]
      for (int i=0; i<oldWave.nrows(); i++)
	for (int j=0; j<oldWave.ncols(); j++)
	  if (oldWave.allowed(i,j)) {
	    assert(trial.allowed(j,i));
	    int parity = getCommuteParity(oldStateInfo.leftStateInfo->quanta[i], oldStateInfo.rightStateInfo->quanta[j], oldWave.get_deltaQuantum(0));
	    for (int x=0; x<oldWave(i,j).Nrows(); x++)
	      for (int y=0; y<oldWave(i,j).Ncols(); y++)
		trial(j,i)(y+1,x+1) = parity*oldWave(i,j)(x+1,y+1);
		  //MatrixScaleAdd(parity*1.0, oldWave(i,j), trial(j,i));

	  }
      oldStateInfo.Free();
      oldWave.deallocate();

    }
    else {
      vector<int> wfsites = big.get_rightBlock()->get_sites();
      wfsites.insert(wfsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
      sort(wfsites.begin(), wfsites.end());
      oldWave.LoadWavefunctionInfo(oldStateInfo, wfsites, state, true);
      if(ket)
	onedot_transpose_wavefunction(oldStateInfo, big.get_stateInfo(), oldWave, trial);
      else
	onedot_transpose_wavefunction(oldStateInfo, big.get_braStateInfo(), oldWave, trial);
      oldStateInfo.Free();
      oldWave.deallocate();
    }
  }

}

void GuessWave::transpose_previous_wavefunction(StackWavefunction& trial, const StateInfo& stateInfo, const std::vector<int>& rightsites, const std::vector<int> &dotsites, const int state, const bool &onedot, const bool& transpose_guess_wave)
{
  StateInfo oldStateInfo;
  StackWavefunction oldWave;

  vector<int> wfsites = rightsites;
  if (dmrginp.spinAdapted())
    wfsites.insert(wfsites.end(), dotsites.begin(), dotsites.end());
  else {
    wfsites.push_back(2*dotsites[0]);
    wfsites.push_back(2*dotsites[0]+1);
  }
  sort(wfsites.begin(), wfsites.end());
  oldWave.LoadWavefunctionInfo(oldStateInfo, wfsites, state, true);

  oldWave.set_deltaQuantum() = trial.set_deltaQuantum();

  StateInfo intermediate1, intermediate2;
  TensorProduct(*stateInfo.rightStateInfo, *stateInfo.leftStateInfo->rightStateInfo, intermediate1, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
  intermediate1.CollectQuanta();
  TensorProduct(intermediate1, *stateInfo.leftStateInfo->leftStateInfo, intermediate2, PARTICLE_SPIN_NUMBER_CONSTRAINT);
  onedot_transpose_wavefunction(intermediate2, stateInfo, oldWave, trial);
  oldStateInfo.Free();
  oldWave.deallocate();
}


/*!                                                                                                                                       
  @brief Transpose wavefunction from [s.][e] to [e.][s] form.                                                                             
*/

void GuessWave::onedot_transpose_wavefunction(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
					      const StackWavefunction& guesswf, StackWavefunction& transposewf)
{
  ObjectMatrix3D< vector<Matrix> > threewave;
  // first convert into three index wavefunction [s.][e] -> [s].[e]                                                                       
  onedot_twoindex_to_threeindex_wavefunction(guessstateinfo, guesswf, threewave);
  ObjectMatrix3D< vector<Matrix> > threewavetranspose(threewave.NDim2(), threewave.NDim1(), threewave.NDim0());

  // now perform transpose. [s].[e] -> [e].[s], so we need to move state |s>|.> past state |e>                                            
  // keeping track of the parity change                                                                                                   
  for (int a = 0; a < threewave.NDim0(); ++a)
    for (int b = 0; b < threewave.NDim1(); ++b)
      for (int c = 0; c < threewave.NDim2(); ++c)
	if (threewave(a, b, c).size () != 0)
	  {
	    threewavetranspose(c, b, a).resize( threewave(a, b, c).size());
	    for (int i=0; i< threewave(a,b,c).size(); i++) 
	    {
	      if (threewave(a,b,c)[i].Ncols() == 0)
		continue;
	      //int ab = guessstateinfo.leftStateInfo->unCollectedStateInfo->quantaMap(a , b)[i];
	      //int absize = guessstateinfo.leftStateInfo->unCollectedStateInfo->quantaMap(a , b).size();
	      copy(threewave(a, b, c)[i].t(), threewavetranspose(c, b, a)[i]);

	      // from |s.e> configuration, determine parity for |e.s> form                                                              

	      SpinQuantum Aq = guessstateinfo.leftStateInfo->leftStateInfo->quanta[a];
	      SpinQuantum Bq = guessstateinfo.leftStateInfo->rightStateInfo->quanta[b];
	      SpinQuantum Cq = guessstateinfo.rightStateInfo->quanta[c];
	      SpinQuantum ABq = (Aq+Bq)[i];
	      //SpinQuantum ABq = guessstateinfo.leftStateInfo->unCollectedStateInfo->quanta[ab];
	      
	      // first, from |s.e> -> |.se>	      
	      int parity1 = getCommuteParity(Aq, Bq, ABq);		

	      // next, |.se> -> |e.s>
	      int parity = parity1*getCommuteParity(ABq, Cq, transposewf.get_deltaQuantum(0));

	      if (parity == -1) 
		threewavetranspose(c, b, a)[i] *= -1.;
	    }
	  }

  onedot_threeindex_to_twoindex_wavefunction(transposestateinfo, threewavetranspose,
					     transposewf, *guessstateinfo.leftStateInfo->unCollectedStateInfo);
}



/*!                                                                                                                                         
@brief Convert wavefunction in three block form (with single dot) [s].[e] to two block form  [s.] [e].
  Note it always group the dot with [s], not with [e].
  Three block form means that the total wavefunction is written explicitly as
  |psi> = c_{s,dot,e} |s>|dot>|e>
  Two block form means that the total wavefunction is written as
  |psi> = c_{sdot,e} |sdot>|e>                                  
  @param twostateinfo: StateInfo for output form [s.][e]                                                                                   
  @param threewavefunction: input wavefunction [s].[e] notation
  @param twowavefunction: output wavefunction in [s.][e] notation                                                                         
*/
void GuessWave::onedot_threeindex_to_twoindex_wavefunction(const StateInfo& twostateinfo,
							   const ObjectMatrix3D< std::vector<Matrix> >& threewavefunction, 
							   StackWavefunction& twowavefunction, const StateInfo& prevUnCollectedSI)
{
  int aSz = threewavefunction.NDim0();
  int bSz = threewavefunction.NDim1();
  int cSz = threewavefunction.NDim2();

  const StateInfo& uncollectedstateinfo = *(twostateinfo.leftStateInfo->unCollectedStateInfo);

  twowavefunction.initialise(twowavefunction.get_deltaQuantum(), uncollectedstateinfo, *twostateinfo.rightStateInfo,
			     twowavefunction.get_onedot(), twowavefunction.get_data(), twowavefunction.memoryUsed()); // twowavefunction should already know its own quantum number
  twowavefunction.Clear();

  for (int a = 0; a < aSz; ++a)
    for (int b = 0; b < bSz; ++b)
      for (int c = 0; c < cSz; ++c)
	if (threewavefunction (a, b, c).size() != 0 && uncollectedstateinfo.allowedQuanta(a, b))
	  for (int i =0; i<uncollectedstateinfo.quantaMap(a,b).size(); i++)
	  {
	    int ab = uncollectedstateinfo.quantaMap (a, b)[i];

	    if (!twowavefunction.allowed(ab,c)) continue;
	    int A, B, AB, C, J, CB;
	    A = twostateinfo.leftStateInfo->leftStateInfo->quanta[a].get_s().getirrep(); 
	    B = twostateinfo.leftStateInfo->rightStateInfo->quanta[b].get_s().getirrep();
	    AB = uncollectedstateinfo.quanta[ab].get_s().getirrep();
	    C =  twostateinfo.rightStateInfo->quanta[c].get_s().getirrep(); 
	    J = twowavefunction.get_deltaQuantum(0).get_s().getirrep(); 

	    int Al, Bl, ABl, Cl, Jl, CBl;
	    Al = twostateinfo.leftStateInfo->leftStateInfo->quanta[a].get_symm().getirrep(); 
	    Bl = twostateinfo.leftStateInfo->rightStateInfo->quanta[b].get_symm().getirrep();
	    ABl = uncollectedstateinfo.quanta[ab].get_symm().getirrep();
	    Cl =  twostateinfo.rightStateInfo->quanta[c].get_symm().getirrep(); 
	    Jl = twowavefunction.get_deltaQuantum(0).get_symm().getirrep(); 

	    for (int j=0; j<prevUnCollectedSI.quantaMap(c, b).size(); j++) {
	      int cb = prevUnCollectedSI.quantaMap(c,b)[j];
	      CB = prevUnCollectedSI.quanta[cb].get_s().getirrep();
	      CBl = prevUnCollectedSI.quanta[cb].get_symm().getirrep();

	      int insertionNum = prevUnCollectedSI.quanta[cb].insertionNum(twostateinfo.leftStateInfo->rightStateInfo->quanta[b], twostateinfo.rightStateInfo->quanta[c]);
	      double scale = 1.0;
	      if(dmrginp.spinAdapted()) {
		scale = sixj(A, B, AB, C, J, CB);
		scale *= pow((1.0*AB+1.0)*(1.0*CB+1.0), 0.5)* pow(-1.0, static_cast<int>((A +B +J +C)/2)); 
	      }
	      scale *= Symmetry::spatial_sixj(Al, Bl, ABl, Cl, Jl, CBl);

	      if (threewavefunction(a, b, c)[insertionNum].Nrows() != 0)  {
		MatrixScaleAdd(scale, threewavefunction (a, b, c)[insertionNum], twowavefunction.operator_element (ab, c));
	      }
	    }
	  }

  // from tensor product form |a>|b> group together blocks with same quantum numbers                                                      
  twowavefunction.CollectQuantaAlongRows(uncollectedstateinfo, *twostateinfo.rightStateInfo);
}

void GuessWave::basic_guess_wavefunction(DiagonalMatrix& e, StackWavefunction& trial, const StateInfo *stateinfo, const int state)
{
  p2out << "\t\t\t No trial vector" << endl;
  multimap<double, int> e_sort;
  for (int i = 0; i < e.Nrows (); ++i) {
    e_sort.insert (pair<double, int> (e (i+1), i+1));
  }

  int states = stateinfo->totalStates;
  RowVector trialvector(states);
  trialvector = 0.;

  int index = 0;
  for (multimap<double, int>::iterator it = e_sort.begin(); it != e_sort.end(); it++) {
    if (index == state) {
      trialvector(it->second) = 1.;
      break;
    }
    index++;
  }
  trial.CollectFrom(trialvector);
}


void GuessWave::guess_wavefunctions(StackWavefunction& solution, DiagonalMatrix& e, const StackSpinBlock &big,
				    const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave,
				    double additional_noise)
{
#ifndef SERIAL
  mpi::communicator world;
#endif

  if (!mpigetrank())
  {
    switch(guesswavetype)
    {
    case TRANSFORM:
      transform_previous_wavefunction(solution, big, state, onedot, transpose_guess_wave);
      break;
    case BASIC:
      basic_guess_wavefunction(e, solution, &big.get_stateInfo(), state);
      break;
    case TRANSPOSE:
      transpose_previous_wavefunction(solution, big, state, onedot, transpose_guess_wave);
      break;
    }

    double norm = DotProduct(solution, solution);

    if (guesswavetype == BASIC) {
      StackWavefunction noiseMatrix; noiseMatrix.initialise( solution);
      noiseMatrix.Randomise();
      double norm = DotProduct(noiseMatrix, noiseMatrix);
      if (abs(norm) >= 1e-14) {
	ScaleAdd(1e-3/sqrt(norm), noiseMatrix, solution);
      }
      Normalise(solution);
      noiseMatrix.deallocate();
    }

    norm = DotProduct(solution, solution);
    p2out << "\t\t\t Norm of wavefunction :: "<<norm<<endl;
    
  }
}


void GuessWave::guess_wavefunctions(StackWavefunction& solution, DiagonalMatrix& e, const StackSpinBlock &big,
				    const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave,
				    double additional_noise, bool ket)
// ket determines use braSateinfo or ketStateinfo to guess_wavefunctions.
{

  if (!mpigetrank())
  {
    switch(guesswavetype)
    {
    case TRANSFORM:
      transform_previous_wavefunction(solution, big, state, onedot, transpose_guess_wave,ket);
      break;
    case BASIC:
      basic_guess_wavefunction(e, solution, &big.get_stateInfo(), state);
      break;
    case TRANSPOSE:
      transpose_previous_wavefunction(solution, big, state, onedot, transpose_guess_wave,ket);
      break;
    }

    double norm = DotProduct(solution, solution);

    if (guesswavetype == BASIC) {
      StackWavefunction noiseMatrix; noiseMatrix.initialise( solution);
      noiseMatrix.Randomise();
      double norm = DotProduct(noiseMatrix, noiseMatrix);
      if (abs(norm) >= 1e-14) {
	ScaleAdd(1e-6/sqrt(norm), noiseMatrix, solution);
      }
      noiseMatrix.deallocate();
    }

    Normalise(solution);
    norm = DotProduct(solution, solution);
    p2out << "\t\t\t Norm of wavefunction :: "<<norm<<endl;
    
  }
}

void GuessWave::guess_wavefunctions(std::vector<StackWavefunction>& solution, DiagonalMatrix& e, const StackSpinBlock &big, 
				    const guessWaveTypes &guesswavetype, const bool &onedot, const bool &transpose_guess_wave, int nroots,
				    double additional_noise, int currentState)
{

  for(int i=0;i<nroots;++i) {
    int state = (dmrginp.setStateSpecific() || dmrginp.calc_type() == COMPRESS) ? currentState : i;
    guess_wavefunctions(solution[i], e, big, guesswavetype, onedot, state, transpose_guess_wave, additional_noise);
  }
}

/*!  
  @brief Convert wavefunction in twoblock form (with single dot) [s.] [e] to three block form [s].[e] 

  Three block form means that the total wavefunction is written explicitly as
  |psi> = c_{s,dot,e} |s>|dot>|e>

  Two block form means that the total wavefunction is written as
  |psi> = c_{sdot,e} |sdot>|e>
  
  @param stateinfo: StateInfo for [s.][e] superblock
  @param twowavefunction: input wavefunction in [s.][e] notation
  @param threewavefunction: output wavefunction [s].[e] notation
 */
void GuessWave::onedot_twoindex_to_threeindex_wavefunction(const StateInfo& stateinfo, const StackWavefunction& twowavefunction, 
							   ObjectMatrix3D< vector<Matrix> >& threewavefunction)
{
  int NDimsys = stateinfo.leftStateInfo->leftStateInfo->quanta.size();
  int NDimdot = stateinfo.leftStateInfo->rightStateInfo->quanta.size();
  int NDimenv = stateinfo.rightStateInfo->quanta.size();

  threewavefunction.ReSize(NDimsys, NDimdot, NDimenv);

  StackWavefunction uncollectedwf; uncollectedwf.initialise(twowavefunction);
  uncollectedwf.copyData(twowavefunction);
  //copy(twowavefunction.get_operatorMatrix(), uncollectedwf.get_operatorMatrix());
  //DCOPY(uncollectedwf.memoryUsed(), const_cast<double*>(twowavefunction.get_data()), 1, uncollectedwf.get_data(), 1);
  uncollectedwf.UnCollectQuantaAlongRows(*stateinfo.leftStateInfo, *stateinfo.rightStateInfo);

  const StateInfo& uncollectedstateinfo = *(stateinfo.leftStateInfo->unCollectedStateInfo);
  for (int ab = 0; ab < uncollectedwf.nrows(); ++ab)
    for (int c = 0; c < uncollectedwf.ncols(); ++c)
      if (uncollectedwf.allowed(ab, c))
	{
	  int a = uncollectedstateinfo.leftUnMapQuanta[ab];
	  int b = uncollectedstateinfo.rightUnMapQuanta[ab];
	  int AB = uncollectedstateinfo.quanta[ab].get_s().getirrep();
	  int A = stateinfo.leftStateInfo->leftStateInfo->quanta[a].get_s().getirrep();
	  int B = stateinfo.leftStateInfo->rightStateInfo->quanta[b].get_s().getirrep();
	  int insertionNum = uncollectedstateinfo.quanta[ab].insertionNum(stateinfo.leftStateInfo->leftStateInfo->quanta[a],  stateinfo.leftStateInfo->rightStateInfo->quanta[b]);
	  vector<SpinQuantum> spq = stateinfo.leftStateInfo->leftStateInfo->quanta[a]+  stateinfo.leftStateInfo->rightStateInfo->quanta[b];
	  if (threewavefunction(a, b, c).size() != spq.size())
	    threewavefunction(a, b, c).resize(spq.size());

	  copy(uncollectedwf(ab, c), threewavefunction(a, b, c)[insertionNum]);
	}
  uncollectedwf.deallocate();
}


void GuessWave::onedot_shufflesysdot(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
				     const StackWavefunction& guesswf, StackWavefunction& transposewf)
{
  ObjectMatrix3D< vector<Matrix> > threewave;
  // first convert into three index wavefunction [s][e.] -> [s].[e]                                                         

  onedot_twoindex_to_threeindex_shufflesysdot(guessstateinfo, guesswf, threewave);
  onedot_threeindex_to_twoindex_wavefunction(transposestateinfo, threewave, transposewf, *guessstateinfo.rightStateInfo->unCollectedStateInfo);
}

void GuessWave::onedot_twoindex_to_threeindex_shufflesysdot(const StateInfo& stateinfo, const StackWavefunction& twowavefunction, 
							    ObjectMatrix3D< vector<Matrix> >& threewavefunction)
{
  int NDimsys = stateinfo.leftStateInfo->quanta.size();
  int NDimdot = stateinfo.rightStateInfo->rightStateInfo->quanta.size();
  int NDimenv = stateinfo.rightStateInfo->leftStateInfo->quanta.size();

  threewavefunction.ReSize(NDimsys, NDimdot, NDimenv);

  StackWavefunction uncollectedwf; uncollectedwf.initialise(twowavefunction);
  uncollectedwf.copyData(twowavefunction);
  //copy(twowavefunction.get_operatorMatrix(), uncollectedwf.get_operatorMatrix());
  //DCOPY(uncollectedwf.memoryUsed(), const_cast<double*>(twowavefunction.get_data()), 1, uncollectedwf.get_data(), 1);
  uncollectedwf.UnCollectQuantaAlongColumns(*stateinfo.leftStateInfo, *stateinfo.rightStateInfo);
  const StateInfo& uncollectedstateinfo = *(stateinfo.rightStateInfo->unCollectedStateInfo);

  for (int a = 0; a < uncollectedwf.nrows(); ++a)
    for (int bc = 0; bc < uncollectedwf.ncols(); ++bc)
      if (uncollectedwf.allowed(a, bc))
        {
          int b = uncollectedstateinfo.rightUnMapQuanta[bc];
          int c = uncollectedstateinfo.leftUnMapQuanta[bc];
          bool bodd = IsFermion(stateinfo.rightStateInfo->rightStateInfo->quanta[b]);
          bool codd = IsFermion(stateinfo.rightStateInfo->leftStateInfo->quanta[c]);
	  //int parity = (bodd&codd) ? -1 : 1;
	  int j1 = stateinfo.rightStateInfo->leftStateInfo->quanta[c].get_s().getirrep();
	  int j2 = stateinfo.rightStateInfo->rightStateInfo->quanta[b].get_s().getirrep();
	  int J = uncollectedstateinfo.quanta[bc].get_s().getirrep();
	  
	  int parity = getCommuteParity(stateinfo.rightStateInfo->leftStateInfo->quanta[c],
				    stateinfo.rightStateInfo->rightStateInfo->quanta[b],
				    uncollectedstateinfo.quanta[bc]);

	  //parity *= pow(-1.0, static_cast<int>( (3*j1 - j2 + J)/2));
	  int insertionNum = uncollectedstateinfo.quanta[bc].insertionNum(stateinfo.rightStateInfo->rightStateInfo->quanta[b], 
									  stateinfo.rightStateInfo->leftStateInfo->quanta[c]);
	  vector<SpinQuantum> spq = stateinfo.rightStateInfo->rightStateInfo->quanta[b]+ stateinfo.rightStateInfo->leftStateInfo->quanta[c];

	  if(threewavefunction(a, b, c).size() != spq.size())
	    threewavefunction(a, b, c).resize(spq.size());
          copy(uncollectedwf(a, bc), threewavefunction(a, b, c)[insertionNum]);
          if (parity == -1)
	    threewavefunction(a,b,c)[insertionNum] *= -1.0;

        }
  uncollectedwf.deallocate();
}


void GuessWave::transform_previous_wavefunction(StackWavefunction& trial, const StateInfo& stateInfo, const std::vector<int> &leftsites, const std::vector<int> &rightsites, const int state, const bool &onedot, const bool& transpose_guess_wave)
{
  p2out << "\t\t\t Transforming previous wavefunction " << endl;
  
  ObjectMatrix3D< vector<Matrix> > oldTrialWavefunction;
  ObjectMatrix3D< vector<Matrix> > newTrialWavefunction;
  StateInfo oldStateInfo;
  StackWavefunction oldWave;
  DiagonalMatrix D;
  Matrix U;
  Matrix V;
  std::vector<Matrix> leftRotationMatrix;

  oldWave.LoadWavefunctionInfo (oldStateInfo, leftsites, state, true);
  LoadRotationMatrix (leftsites, leftRotationMatrix, state);


  for (int q = 0; q < leftRotationMatrix.size (); ++q)
  {
    if (leftRotationMatrix [q].Nrows () > 0)
    {
      leftRotationMatrix[q] = leftRotationMatrix[q].t();
    }
  }

  std::vector<Matrix> rightRotationMatrix;
  LoadRotationMatrix(rightsites, rightRotationMatrix, state);
    
  onedot_transform_wavefunction(oldStateInfo, stateInfo, oldWave, leftRotationMatrix, rightRotationMatrix, trial, transpose_guess_wave);

  double norm = DotProduct(trial, trial);
  oldStateInfo.Free();
  oldWave.deallocate();
}



void GuessWave::transform_previous_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave)
{
  p2out << "\t\t\t Transforming previous wavefunction " << endl;
  
  ObjectMatrix3D< vector<Matrix> > oldTrialWavefunction;
  ObjectMatrix3D< vector<Matrix> > newTrialWavefunction;
  StateInfo oldStateInfo;
  StackWavefunction oldWave;
  DiagonalMatrix D;
  Matrix U;
  Matrix V;
  std::vector<Matrix> leftRotationMatrix;
  if (transpose_guess_wave || !onedot){
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_leftBlock()->get_sites(), state, true);
    LoadRotationMatrix (big.get_leftBlock()->get_leftBlock()->get_sites(), leftRotationMatrix, state);
  }
  else{
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_sites(), state, true);
    LoadRotationMatrix (big.get_leftBlock()->get_sites(), leftRotationMatrix, state);
  }

  for (int q = 0; q < leftRotationMatrix.size (); ++q)
  {
    if (leftRotationMatrix [q].Nrows () > 0)
    {

      leftRotationMatrix[q] = leftRotationMatrix[q].t();
    }
  }

  std::vector<Matrix> rightRotationMatrix;
  if(!onedot)
  {
    StackWavefunction tempoldWave;
    tempoldWave.initialise(oldWave.get_deltaQuantum(), *big.get_stateInfo().leftStateInfo->leftStateInfo, *oldStateInfo.rightStateInfo, false); 
    tempoldWave.Clear();

    TransformLeftBlock(oldWave, big.get_stateInfo(), leftRotationMatrix, tempoldWave);

    StateInfo tempoldStateInfo;
    if (dmrginp.hamiltonian() == BCS)
      TensorProduct (*(big.get_stateInfo().leftStateInfo->leftStateInfo), *oldStateInfo.rightStateInfo, tempoldStateInfo,
		   SPIN_NUMBER_CONSTRAINT);
    else
      TensorProduct (*(big.get_stateInfo().leftStateInfo->leftStateInfo), *oldStateInfo.rightStateInfo, tempoldStateInfo,
		   PARTICLE_SPIN_NUMBER_CONSTRAINT);

    tempoldStateInfo.CollectQuanta();

    StackWavefunction tempnewWave;
    tempnewWave.initialise(oldWave.get_deltaQuantum(), *big.get_stateInfo().leftStateInfo, *oldStateInfo.rightStateInfo->leftStateInfo, oldWave.get_onedot()); 
    tempnewWave.Clear();
    StateInfo tempnewStateInfo;
    if (dmrginp.hamiltonian() == BCS)
      TensorProduct (*(big.get_stateInfo().leftStateInfo), *oldStateInfo.rightStateInfo->leftStateInfo, tempnewStateInfo,
		   SPIN_NUMBER_CONSTRAINT);
    else
      TensorProduct (*(big.get_stateInfo().leftStateInfo), *oldStateInfo.rightStateInfo->leftStateInfo, tempnewStateInfo,
		   PARTICLE_SPIN_NUMBER_CONSTRAINT);

    
    tempnewStateInfo.CollectQuanta();
    onedot_shufflesysdot(tempoldStateInfo, tempnewStateInfo, tempoldWave, tempnewWave);

    LoadRotationMatrix (big.get_rightBlock()->get_sites(), rightRotationMatrix, state);

    TransformRightBlock(tempnewWave, oldStateInfo, rightRotationMatrix, trial);

    tempnewWave.deallocate();
    tempoldWave.deallocate();
  }
  else
  {
    vector<int> rotsites;
    if (transpose_guess_wave) {
      rotsites = big.get_rightBlock()->get_sites();
      rotsites.insert(rotsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
      sort(rotsites.begin(), rotsites.end());
      LoadRotationMatrix(rotsites, rightRotationMatrix, state);
    }
    else {
      rotsites = big.get_rightBlock()->get_sites();
      LoadRotationMatrix(rotsites, rightRotationMatrix, state);
    }
    onedot_transform_wavefunction(oldStateInfo, big.get_stateInfo(), oldWave, leftRotationMatrix, rightRotationMatrix, trial, transpose_guess_wave);
  }

  oldStateInfo.Free();
  oldWave.deallocate();
  double norm = DotProduct(trial, trial);
}

void GuessWave::transform_previous_twodot_to_onedot_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state)
{
  p2out << "\t\t\t Transforming previous wavefunction " << endl;
  
  ObjectMatrix3D< vector<Matrix> > oldTrialWavefunction;
  ObjectMatrix3D< vector<Matrix> > newTrialWavefunction;
  StateInfo oldStateInfo;
  StackWavefunction oldWave;
  DiagonalMatrix D;
  Matrix U;
  Matrix V;
  std::vector<Matrix> leftRotationMatrix;

  oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_leftBlock()->get_sites(), state, true);
  LoadRotationMatrix (big.get_leftBlock()->get_leftBlock()->get_sites(), leftRotationMatrix, state);

  for (int q = 0; q < leftRotationMatrix.size (); ++q)
    if (leftRotationMatrix [q].Nrows () > 0)
      leftRotationMatrix[q] = leftRotationMatrix[q].t();



  int aSz = big.get_stateInfo().leftStateInfo->leftStateInfo->quanta.size (); 
  int cSz = oldStateInfo.rightStateInfo->quanta.size();

  StackWavefunction tmpwavefunction;
  tmpwavefunction.initialise(dmrginp.effective_molecule_quantum_vec(), *big.get_stateInfo().leftStateInfo->leftStateInfo, *oldStateInfo.rightStateInfo, false);
  tmpwavefunction.Clear();

  //now contract the [s.] => [s] using the left rotation matrix
  for (int a = 0; a < aSz; ++a) //aSz // aSz <= oldASz
    for (int c = 0; c < cSz; ++c) //cSz
    {
      int oldA = big.get_stateInfo().leftStateInfo->leftStateInfo->newQuantaMap [a];
      
      if (oldWave.allowed(oldA, c)) {
	StackMatrix& tM = oldWave(oldA,c); //tmp (oldA, b, c);
	const Matrix& lM = leftRotationMatrix [oldA];
	StackMatrix& nM = tmpwavefunction.operator_element(a, c);
	//nM.ReSize (lM.Nrows (), tM.Ncols ()); //***********
	//SpinAdapted::Clear (nM);
	MatrixMultiply (lM, 'n', tM, 'n', nM, 1.); 
      }
    }

  StateInfo tempoldStateInfo;
  trial.set_deltaQuantum() = dmrginp.effective_molecule_quantum_vec();
  trial.set_initialised() = true;
  trial.set_fermion() = false;
  trial.set_onedot(true);

  TensorProduct( *big.get_stateInfo().leftStateInfo->leftStateInfo, *oldStateInfo.rightStateInfo,  tempoldStateInfo, PARTICLE_SPIN_NUMBER_CONSTRAINT);
  onedot_shufflesysdot(  tempoldStateInfo, big.get_stateInfo(), tmpwavefunction, trial);

  tmpwavefunction.deallocate();
  oldWave.deallocate();
  oldStateInfo.Free();

  double norm = DotProduct(trial, trial);
}


void GuessWave::transform_previous_wavefunction(StackWavefunction& trial, const StackSpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave,bool ket)
// ket determines use braSateinfo or ketStateinfo to guess_wavefunctions.
{
  p2out << "\t\t\t Transforming previous wavefunction " << endl;
  
  ObjectMatrix3D< vector<Matrix> > oldTrialWavefunction;
  ObjectMatrix3D< vector<Matrix> > newTrialWavefunction;
  StateInfo oldStateInfo;
  StackWavefunction oldWave;
  DiagonalMatrix D;
  Matrix U;
  Matrix V;
  std::vector<Matrix> leftRotationMatrix;
  if (transpose_guess_wave || !onedot){
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_leftBlock()->get_sites(), state, true);
    LoadRotationMatrix (big.get_leftBlock()->get_leftBlock()->get_sites(), leftRotationMatrix, state);
  }
  else{
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_sites(), state, true);
    LoadRotationMatrix (big.get_leftBlock()->get_sites(), leftRotationMatrix, state);
  }

  for (int q = 0; q < leftRotationMatrix.size (); ++q)
  {
    if (leftRotationMatrix [q].Nrows () > 0)
    {

      leftRotationMatrix[q] = leftRotationMatrix[q].t();
    }
  }

  std::vector<Matrix> rightRotationMatrix;
  if(!onedot)
  {
    StackWavefunction tempoldWave;
    tempoldWave.initialise(oldWave.get_deltaQuantum(), *big.get_stateInfo().leftStateInfo->leftStateInfo, *oldStateInfo.rightStateInfo, oldWave.get_onedot()); 
    tempoldWave.Clear();
    TransformLeftBlock(oldWave, big.get_stateInfo(), leftRotationMatrix, tempoldWave);

    StateInfo tempoldStateInfo;
    if (dmrginp.hamiltonian() == BCS)
      TensorProduct (*(big.get_stateInfo().leftStateInfo->leftStateInfo), *oldStateInfo.rightStateInfo, tempoldStateInfo,
		   SPIN_NUMBER_CONSTRAINT);
    else
      TensorProduct (*(big.get_stateInfo().leftStateInfo->leftStateInfo), *oldStateInfo.rightStateInfo, tempoldStateInfo,
		   PARTICLE_SPIN_NUMBER_CONSTRAINT);

    tempoldStateInfo.CollectQuanta();

    StackWavefunction tempnewWave;
    tempnewWave.initialise(oldWave.get_deltaQuantum(), *big.get_stateInfo().leftStateInfo, *oldStateInfo.rightStateInfo->leftStateInfo, oldWave.get_onedot()); 
    tempnewWave.Clear();
    StateInfo tempnewStateInfo;
    if (dmrginp.hamiltonian() == BCS)
      TensorProduct (*(big.get_stateInfo().leftStateInfo), *oldStateInfo.rightStateInfo->leftStateInfo, tempnewStateInfo,
		   SPIN_NUMBER_CONSTRAINT);
    else
      TensorProduct (*(big.get_stateInfo().leftStateInfo), *oldStateInfo.rightStateInfo->leftStateInfo, tempnewStateInfo,
		   PARTICLE_SPIN_NUMBER_CONSTRAINT);

    
    tempnewStateInfo.CollectQuanta();
    onedot_shufflesysdot(tempoldStateInfo, tempnewStateInfo, tempoldWave, tempnewWave);

    LoadRotationMatrix (big.get_rightBlock()->get_sites(), rightRotationMatrix, state);

    //trial.AllowQuantaFor(*big.get_stateInfo().leftStateInfo, *big.get_stateInfo().rightStateInfo, oldWave.get_deltaQuantum()); //**********
    TransformRightBlock(tempnewWave, oldStateInfo, rightRotationMatrix, trial);

    tempnewWave.deallocate();
    tempoldWave.deallocate();
  }
  else
  {
    vector<int> rotsites;
    if (transpose_guess_wave) {
      rotsites = big.get_rightBlock()->get_sites();
      rotsites.insert(rotsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
      sort(rotsites.begin(), rotsites.end());
      LoadRotationMatrix(rotsites, rightRotationMatrix, state);
    }
    else {
      rotsites = big.get_rightBlock()->get_sites();
      LoadRotationMatrix(rotsites, rightRotationMatrix, state);
    }
    if(ket)
      onedot_transform_wavefunction(oldStateInfo, big.get_stateInfo(), oldWave, leftRotationMatrix, rightRotationMatrix, trial, transpose_guess_wave,ket);
    else 
      onedot_transform_wavefunction(oldStateInfo, big.get_braStateInfo(), oldWave, leftRotationMatrix, rightRotationMatrix, trial, transpose_guess_wave,ket);
  }

  oldStateInfo.Free();
  oldWave.deallocate();

  double norm = DotProduct(trial, trial);
}

/*!  
  @brief Transform wavefunction from previous block configuration 
  [s.][e] -> [s'.][e'] (where s' = R[s.]; [R denotes renormalise] i.e. next block configuration)
  
  @param oldstateinfo: StateInfo for input configuration [s.][e] 
  @param newstateinfo: StateInfo for output configuration [s'.][e'] 
  @param oldWavefunction: input wavefunction
  @param inverseLeftRotationMatrix: inverse of transformation  from [s.]->[s']
  @param rightRotationMatrix: transformation from [.e'] -> [e]
  @param newwavefunction: output wavefunction. Should already have quantum number assigned.
 */


void GuessWave::onedot_transform_wavefunction(const StateInfo& oldstateinfo, const StateInfo& newstateinfo, 
					      const StackWavefunction& oldwavefunction, const vector<Matrix>& leftRotationMatrix, 
					      const vector<Matrix>& rightRotationMatrix, StackWavefunction& newwavefunction, 
					      const bool& transpose_guess_wave, bool ket)
{
  assert (oldwavefunction.get_deltaQuantum() == newwavefunction.get_deltaQuantum());
  Timer transform;
  int oldASz = oldstateinfo.leftStateInfo->quanta.size ();
  int oldCSz = oldstateinfo.rightStateInfo->quanta.size ();
  // Old wavefunction is in [s.][e] configuration
  // Then [s.] -> [s'], and [.e'] -> [e]
  // Make wavefunction in [s'][.e'] configuration 
  // Then later transfer dot, [s'][.e'] -> [s'.][e']
  StateInfo newenvstateinfo;
  ObjectMatrix<char> preparedQuanta;
  StackWavefunction tmpwavefunction;
  int aSz;
  int cSz;
  if (transpose_guess_wave){
    if(dmrginp.transition_diff_irrep() && (ket == false))
      TensorProduct (*(newstateinfo.rightStateInfo), *(newstateinfo.leftStateInfo->rightStateInfo), dmrginp.bra_quantum(), LessThanQ, newenvstateinfo);
    else
      TensorProduct (*(newstateinfo.rightStateInfo), *(newstateinfo.leftStateInfo->rightStateInfo), newenvstateinfo,
		   NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newenvstateinfo.CollectQuanta();
    tmpwavefunction.initialise (newwavefunction.get_deltaQuantum(), *(newstateinfo.leftStateInfo->leftStateInfo), newenvstateinfo, newwavefunction.get_onedot());
    tmpwavefunction.Clear();
    aSz = newstateinfo.leftStateInfo->leftStateInfo->quanta.size (); 
    cSz = newenvstateinfo.quanta.size();
  }
  else {
    tmpwavefunction.initialise (newwavefunction.get_deltaQuantum(), *(newstateinfo.leftStateInfo), *(newstateinfo.rightStateInfo), newwavefunction.get_onedot());
    tmpwavefunction.Clear();
    aSz = newstateinfo.leftStateInfo->quanta.size (); 
    cSz = newstateinfo.rightStateInfo->quanta.size();
  }

  ObjectMatrix<Matrix> tmp(oldASz, cSz); //cSz
  for (int a = 0; a < oldASz; ++a)
    for (int c = 0; c < oldCSz; ++c) // oldCSz <= cSz
    {
      if (oldwavefunction.allowed(a,c)) {
	const StackMatrix& oM = oldwavefunction.operator_element(a, c);

	int transC = oldstateinfo.rightStateInfo->newQuantaMap [c];

	Matrix& tM = tmp (a, transC);
	Matrix rM = rightRotationMatrix[transC];
	tM.ReSize (oM.Nrows (), rM.Nrows ());
	
	double parity = 1.;
	
	SpinAdapted::Clear (tM);
	// we have the rotation matrices C_r'r, C_l'l, and we are transforming the wavefunction
	// coefficient d_lr -> consequently, we need to multiply by the transpose, but NOT hermitian
	// conjugate, of rM
	rM = rM.t();
	MatrixMultiply (oM, 'n', rM, 'n', tM, parity); // phi -> c psi : phi d -> psi c c^t d : consequently compute c^t d 
      }
    }

  for (int c = 0; c < cSz; ++c) //cSz
    for (int a = 0; a < aSz; ++a) //aSz // aSz <= oldASz
    {
      int oldA = 0;
      if (transpose_guess_wave)
	oldA = newstateinfo.leftStateInfo->leftStateInfo->newQuantaMap [a];
      else
	oldA = newstateinfo.leftStateInfo->newQuantaMap [a];

      Matrix& tM = tmp(oldA,c); //tmp (oldA, b, c);
      if (tM.Ncols () != 0) // this quanta combination is not allowed
      {
	//assert (newstateinfo.leftStateInfo->leftStateInfo->quanta [a] == oldstateinfo.leftStateInfo->quanta [oldA]);
	const Matrix& lM = leftRotationMatrix [oldA];
	//assert (newstateinfo.leftStateInfo->leftStateInfo->quantaStates [a] == lM.Nrows ());
	StackMatrix& nM = tmpwavefunction.operator_element(a, c);//tmpwavefunction (a, b, c);
	MatrixMultiply (lM, 'n', tM, 'n', nM, 1.); 
      }
    }

  if (!transpose_guess_wave){
    newwavefunction.copyData(tmpwavefunction);
    //copy(tmpwavefunction.get_operatorMatrix(), newwavefunction.get_operatorMatrix());
    //newwavefunction = tmpwavefunction;
    newwavefunction.set_onedot(oldwavefunction.get_onedot());
    tmpwavefunction.deallocate();
    return;
  }

  // Now, wavefunction is in [s'][e'.] config, change to [s'.][e'] config.
  StateInfo tempoldStateInfo;
  if (dmrginp.hamiltonian() == BCS)
    TensorProduct (*(newstateinfo.leftStateInfo->leftStateInfo), newenvstateinfo, tempoldStateInfo,
		SPIN_NUMBER_CONSTRAINT);
  else if (dmrginp.transition_diff_irrep() && (ket == false))
    TensorProduct (*(newstateinfo.leftStateInfo->leftStateInfo), newenvstateinfo, dmrginp.bra_quantum(), EqualQ, tempoldStateInfo);
  else
    TensorProduct (*(newstateinfo.leftStateInfo->leftStateInfo), newenvstateinfo, tempoldStateInfo,
		 PARTICLE_SPIN_NUMBER_CONSTRAINT);
  tempoldStateInfo.CollectQuanta();


  onedot_shufflesysdot(  tempoldStateInfo, newstateinfo, tmpwavefunction, newwavefunction);
  tmpwavefunction.deallocate();
}

}
