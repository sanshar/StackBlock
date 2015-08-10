/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "StackOperators.h"
#include "Stackwavefunction.h"
#include "Stackspinblock.h"
#include "SpinQuantum.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>
#include "pario.h"


long SpinAdapted::getRequiredMemoryForWavefunction(const StateInfo& sr, const StateInfo& sc, const std::vector<SpinQuantum>& q) {

  long memory = 0;

  for (int i=0; i < sr.quanta.size(); i++)
    for (int j=0; j<sc.quanta.size(); j++) {
      bool allowedcoupling = false;
      for (int k=0; k<q.size(); k++) {
	if (q[k].allow(sr.quanta[i], sc.quanta[j])) {
	  allowedcoupling = true;
	  break;
	}
      }
      if (allowedcoupling)
	memory += sr.quantaStates[i]*sc.quantaStates[j]+ CACHEBUFFER;
    }

  return memory;
}


void SpinAdapted::StackWavefunction::initialise(const StackWavefunction& w)
{
  *this = w;
  data = Stackmem.allocate(totalMemory);

  long index = 0;
  for (int i = 0; i<nonZeroBlocks.size(); i++) {
    int lQ = nonZeroBlocks[i].first.first, rQ = nonZeroBlocks[i].first.second;
    operatorMatrix(lQ,rQ).allocate(&data[index], operatorMatrix(lQ, rQ).Nrows(), operatorMatrix(lQ, rQ).Ncols());
    index += operatorMatrix(lQ, rQ).Nrows()* operatorMatrix(lQ, rQ).Ncols()+ CACHEBUFFER;
    nonZeroBlocks[i].second= operatorMatrix(lQ,rQ);
  }
  if (index != totalMemory) exit(0);
}

void SpinAdapted::StackWavefunction::initialise(const vector<SpinQuantum>& dQ, const StateInfo& sl, const StateInfo& sr, const bool &onedot_)
{
  long totalMemory = getRequiredMemoryForWavefunction(sl, sr, dQ);
  double* data = Stackmem.allocate(totalMemory);
  //memset(data, 0, totalMemory*sizeof(double));
  initialise(dQ, sl, sr, onedot_, data, totalMemory);  
}


void SpinAdapted::StackWavefunction::initialise(const vector<SpinQuantum>& dQ, const StateInfo& sl, const StateInfo& sr, const bool &onedot_, double* pData, long ptotalMemory)
{
  // initialized a ket wavefunction
  operatorMatrix.resize(0,0);
  data = pData;
  totalMemory = ptotalMemory;
  initialised = true;
  fermion = false;
  built = true;
  deltaQuantum = dQ;
  onedot = onedot_; 


  nonZeroBlocks.resize(0);
  mapToNonZeroBlocks.clear();
  rowCompressedForm.clear();rowCompressedForm.resize(sl.quanta.size(), vector<int>());
  colCompressedForm.clear();colCompressedForm.resize(sr.quanta.size(), vector<int>());
  operatorMatrix.resize(sl.quanta.size (), sr.quanta.size ());
  allowedQuantaMatrix.resize(sl.quanta.size (), sr.quanta.size ());

  long index = 0;
  for (int lQ = 0; lQ < sl.quanta.size (); ++lQ)
    for (int rQ = 0; rQ < sr.quanta.size (); ++rQ) {
      allowedQuantaMatrix(lQ, rQ) = false;
      for (int i = 0; i < deltaQuantum.size(); ++i)
        if (deltaQuantum[i].allow(sl.quanta [lQ] , sr.quanta [rQ])) {
          allowedQuantaMatrix(lQ, rQ) = true;
	  rowCompressedForm[lQ].push_back(rQ);
	  colCompressedForm[rQ].push_back(lQ);
          break;
        }

      if (allowedQuantaMatrix(lQ, rQ))
      {
	operatorMatrix(lQ,rQ).allocate(&data[index], sl.quantaStates [lQ], sr.quantaStates [rQ]);
	index += sl.quantaStates [lQ]* sr.quantaStates [rQ]+ CACHEBUFFER;
	nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(lQ,rQ), operatorMatrix(lQ,rQ)));
	mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(lQ, rQ), nonZeroBlocks.size()-1));
      }

    }

}


/*
void SpinAdapted::StackWavefunction::FlattenInto (Matrix& C)
{
  C.ReSize(totalMemory,1);
  DCOPY(totalMemory, data, 1, C.Store(), 1);
}
*/

void SpinAdapted::StackWavefunction::CollectFrom (const RowVector& C)
{
  assert(totalMemory == C.Storage());
  long index = 0;
  for (int lQ = 0; lQ < nrows(); ++lQ)
  for (int rQ = 0; rQ < ncols(); ++rQ)
  if (allowed(lQ, rQ)) {
    int copysize = operatorMatrix(lQ, rQ).Nrows()* operatorMatrix(lQ, rQ).Ncols();
    DCOPY(copysize, C.Store()+(index), 1, operatorMatrix(lQ, rQ).Store(), 1);
    index += copysize;
  }
}


  
void SpinAdapted::StackWavefunction::SaveWavefunctionInfo (const StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num)
{
  char file [5000];
  int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", wave_num, ".tmp");
  p1out << "\t\t\t Saving Wavefunction " << file << endl;
  if (mpigetrank() == 0)
    {
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_wave(ofs);
      save_wave << onedot << waveInfo << *waveInfo.leftStateInfo << *(waveInfo.leftStateInfo->leftStateInfo);
      save_wave << *(waveInfo.leftStateInfo->rightStateInfo) << *waveInfo.rightStateInfo;
      if(!onedot)
	save_wave << *(waveInfo.rightStateInfo->leftStateInfo) << *(waveInfo.rightStateInfo->rightStateInfo);

      this->Save (ofs);
      ofs.close();
    }

}

void SpinAdapted::StackWavefunction::LoadWavefunctionInfo (StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num, bool allocateData)
{
  char file [5000];
  int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
  sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.load_prefix().c_str(), "/wave-", first, "-", last, ".", mpigetrank(), ".", wave_num, ".tmp");
  p1out << "\t\t\t Loading Wavefunction " << file << endl;
  waveInfo.Allocate ();
  if (mpigetrank() == 0)
    {
      std::ifstream ifs(file, std::ios::binary);
      boost::archive::binary_iarchive load_wave(ifs);
      load_wave >> onedot >> waveInfo >> *waveInfo.leftStateInfo >> *(waveInfo.leftStateInfo->leftStateInfo)
		>> *(waveInfo.leftStateInfo->rightStateInfo) >> *waveInfo.rightStateInfo;
      if(!onedot)
	load_wave >> *(waveInfo.rightStateInfo->leftStateInfo) >> *(waveInfo.rightStateInfo->rightStateInfo);

      this->Load (ifs, allocateData);
      ifs.close();
      allocateOperatorMatrix(); 
    }
}

void SpinAdapted::StackWavefunction::CollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol)
{
  //mdebugcheck("before collectquantaalongrows");
  try
    {
      Timer ctimer;

      StateInfo tmpState = sRow;
      tmpState.CollectQuanta ();

      StackWavefunction tmpOper;
      tmpOper.initialise(deltaQuantum, tmpState, sCol, onedot);

      rowCompressedForm.clear();rowCompressedForm.resize(tmpState.quanta.size(), vector<int>());
      colCompressedForm.clear();colCompressedForm.resize(sCol.quanta.size(), vector<int>());
      nonZeroBlocks.resize(0);
      mapToNonZeroBlocks.clear();


      ObjectMatrix<StackMatrix*> matRef;
      for (int i = 0; i < tmpState.quanta.size (); ++i)
	for (int j = 0; j < sCol.quanta.size (); ++j) {
	  std::vector<int> dum (1); dum [0] = j;
	  if (tmpOper.allowed(i, j)) {
	    OperatorMatrixReference (matRef, tmpState.oldToNewState [i], dum);
	    CatenateProduct (matRef, tmpOper.operator_element(i,j));

	    rowCompressedForm[i].push_back(j);
	    colCompressedForm[j].push_back(i);
	    nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(i,j), tmpOper.operatorMatrix(i,j)));
	    mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(i, j), nonZeroBlocks.size()-1));
	  }
	}

      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      operatorMatrix = tmpOper.operatorMatrix;

      allocateOperatorMatrix(); 
      copy(tmpOper.operatorMatrix, operatorMatrix);
      //DCOPY(totalMemory, tmpOper.data, 1, data, 1);

      tmpOper.deallocate();
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongrows");
}



void SpinAdapted::StackWavefunction::UnCollectQuantaAlongRows (const StateInfo& sRow, const StateInfo& sCol)
{
  try
    {
      rowCompressedForm.clear();rowCompressedForm.resize(sRow.unCollectedStateInfo->quanta.size(), vector<int>());
      colCompressedForm.clear();colCompressedForm.resize(sCol.quanta.size(), vector<int>());
      nonZeroBlocks.resize(0);
      mapToNonZeroBlocks.clear();

      StackWavefunction tmpOper;
      tmpOper.initialise(deltaQuantum, *sRow.unCollectedStateInfo, sCol, onedot);
      tmpOper.Clear();

      for (int i = 0; i < sRow.quanta.size (); ++i)
	{
	  const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
	  int firstRow = 0;
	  for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	    {
	      int unCollectedI = oldToNewStateI [iSub];
	      int lastRowSize = sRow.unCollectedStateInfo->quantaStates [unCollectedI];
	      for (int j = 0; j < sCol.quanta.size (); ++j)
		if (tmpOper.allowedQuantaMatrix (unCollectedI, j)) {
		  for (int row = 0; row<lastRowSize; row++) 
		    for (int col = 0; col <tmpOper.operator_element(unCollectedI, j).Ncols(); col++) {
		      tmpOper.operator_element(unCollectedI, j)(row+1, col+1) = operator_element(i, j)( row + firstRow + 1, col+1);
		    }


		  rowCompressedForm[unCollectedI].push_back(j);
		  colCompressedForm[j].push_back(unCollectedI);
		  nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(unCollectedI,j), tmpOper.operatorMatrix(unCollectedI, j)));
		  mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(unCollectedI, j), nonZeroBlocks.size()-1));
		}
	      firstRow += lastRowSize;
	    }
	}
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      operatorMatrix = tmpOper.operatorMatrix;

      allocateOperatorMatrix(); 
      copy(tmpOper.operatorMatrix, operatorMatrix);
      
      tmpOper.deallocate();

    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
}



void SpinAdapted::StackWavefunction::CollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol)
{
  //mdebugcheck("before collectquantaalongcolumns");
  try
    {
      StateInfo tmpState = sCol;
      tmpState.CollectQuanta ();

      StackWavefunction tmpOper;
      tmpOper.initialise (deltaQuantum, sRow, tmpState, onedot);
      tmpOper.Clear();

      rowCompressedForm.clear();rowCompressedForm.resize(sRow.quanta.size(), vector<int>());
      colCompressedForm.clear();colCompressedForm.resize(tmpState.quanta.size(), vector<int>());
      nonZeroBlocks.resize(0);
      mapToNonZeroBlocks.clear();

      ObjectMatrix<StackMatrix*> matRef;
      for (int i = 0; i < sRow.quanta.size (); ++i)
	for (int j = 0; j < tmpState.quanta.size (); ++j)
	  {
	    std::vector<int> dum (1); dum [0] = i;
	    if (tmpOper.allowed(i, j))
	      {
		OperatorMatrixReference (matRef, dum, tmpState.oldToNewState [j]);
		CatenateProduct (matRef, tmpOper.operator_element(i,j));

		rowCompressedForm[i].push_back(j);
		colCompressedForm[j].push_back(i);
		nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(i,j), tmpOper.operatorMatrix(i,j)));
		mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(i, j), nonZeroBlocks.size()-1));
	      }
	  }
      allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
      operatorMatrix = tmpOper.operatorMatrix;

      allocateOperatorMatrix(); 
      copy(tmpOper.operatorMatrix, operatorMatrix);
      //DCOPY(totalMemory, tmpOper.data, 1, data, 1);

      tmpOper.deallocate();
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }
  //mdebugcheck("after collectquantaalongcolumns");
}


void  SpinAdapted::StackWavefunction::UnCollectQuantaAlongColumns (const StateInfo& sRow, const StateInfo& sCol)
{
  //try                                                                                                                                   
  // {                                                                                                                                    
  rowCompressedForm.clear();rowCompressedForm.resize(sRow.quanta.size(), vector<int>());
  colCompressedForm.clear();colCompressedForm.resize(sCol.unCollectedStateInfo->quanta.size(), vector<int>());
  nonZeroBlocks.resize(0);
  mapToNonZeroBlocks.clear();

  StackWavefunction tmpOper;
  tmpOper.initialise(deltaQuantum, sRow, *sCol.unCollectedStateInfo, onedot);
  tmpOper.Clear();

  for (int j = 0; j < sRow.quanta.size (); ++j)
  for (int i = 0; i < sCol.quanta.size (); ++i)
  {
    const std::vector<int>& oldToNewStateI = sCol.oldToNewState [i];
    int firstCol = 0;
    for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
    {
      int unCollectedI = oldToNewStateI [iSub];
      int lastColSize = sCol.unCollectedStateInfo->quantaStates [unCollectedI];
      if (tmpOper.allowed(j, unCollectedI)){
	for (int row = 0; row<tmpOper.operator_element(j, unCollectedI).Nrows(); row++) 
	for (int col = 0; col <lastColSize; col++) 
	  tmpOper.operator_element(j, unCollectedI)(row+1, col+1) = operator_element(j, i)( row + 1, col+1+firstCol);
	      	      
	rowCompressedForm[j].push_back(unCollectedI);
	colCompressedForm[unCollectedI].push_back(j);
	nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(j, unCollectedI), tmpOper.operatorMatrix(j, unCollectedI)));
	mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(j, unCollectedI), nonZeroBlocks.size()-1));
	
      }
      firstCol += lastColSize;
    }
  }
  allowedQuantaMatrix = tmpOper.allowedQuantaMatrix;
  operatorMatrix = tmpOper.operatorMatrix;
  allocateOperatorMatrix(); 
  copy(tmpOper.operatorMatrix, operatorMatrix);

  
  tmpOper.deallocate();
}




SpinAdapted::StackWavefunction& SpinAdapted::StackWavefunction::operator+=(const StackWavefunction& other)
{
  for (int i = 0; i < nrows(); ++i)
    for (int j = 0; j < ncols(); ++j)
      if (allowed(i, j))
        {
          assert(other.allowed(i, j));
          MatrixScaleAdd(1., other.operator_element(i, j), operator_element(i, j));
        }
  return *this;
}
