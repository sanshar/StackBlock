/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "StackOperators.h"
#include "rotationmat.h"
#include "pario.h"
#include "MatrixBLAS.h"
#include <include/sortutils.h>
#include <boost/serialization/vector.hpp>
#include "pario.h"
#include "cmath"
using namespace boost;
using namespace std;

void SpinAdapted::UnCollectQuantaAlongRows(std::vector<Matrix>& rotation, const StateInfo&  sRow) {

   try
     {
       std::vector<Matrix> tmpOper(sRow.unCollectedStateInfo->quanta.size());
       for (int i=0; i<sRow.quanta.size(); i++) {
	 if (rotation[i].Ncols() == 0) continue;
	 const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
	 for (int iSub = 0; iSub <oldToNewStateI.size(); iSub++) {
	   int unCollectedI = oldToNewStateI[iSub];
	   tmpOper[unCollectedI] = Matrix(sRow.unCollectedStateInfo->quantaStates[unCollectedI], rotation[i].Ncols());
	   tmpOper[unCollectedI] = 0.0;
	 }
       }

       for (int i = 0; i < sRow.quanta.size (); ++i)
	 {
	   if (rotation[i].Ncols() == 0) continue;
	   const std::vector<int>& oldToNewStateI = sRow.oldToNewState [i];
	   int firstRow = 0;
	   for (int iSub = 0; iSub < oldToNewStateI.size (); ++iSub)
	     {
	       int unCollectedI = oldToNewStateI [iSub];
	       int lastRowSize = sRow.unCollectedStateInfo->quantaStates [unCollectedI];
	       //for (int j = 0; j < sCol.quanta.size (); ++j)
	       //if (sRow.unCollectedStateInfo->quanta[unCollectedI] == sCol.quanta[j]) {
	       for (int row = 0; row<lastRowSize; row++) 
		 for (int col = 0; col <tmpOper[unCollectedI].Ncols(); col++) {
		   tmpOper[unCollectedI](row+1, col+1) = rotation[i](row+firstRow+1, col+1);
		   //tmpOper.operator_element(unCollectedI, j)(row+1, col+1) = operatorMatrix(i, j)( row + firstRow + 1, col+1);
		 }
	       firstRow += lastRowSize;
	     }
	 }
       rotation = tmpOper;
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }


}


void SpinAdapted::CollectQuantaAlongRows(std::vector<Matrix>& rotation, const StateInfo&  sRow) {

   try
     {
       StateInfo tmpState = sRow;
       tmpState.CollectQuanta ();

       std::vector<Matrix> tmpOper(tmpState.quanta.size());
       for (int i=0; i<tmpState.quanta.size(); i++) {
	 const std::vector<int>& oldToNewStateI = tmpState.oldToNewState [i];
	 for (int iSub = 0; iSub <oldToNewStateI.size(); iSub++) {
	   int unCollectedI = oldToNewStateI[iSub];
	   if (rotation[unCollectedI].Ncols() != 0) {
	     tmpOper[i] = Matrix(tmpState.quantaStates[i], rotation[unCollectedI].Ncols());
	     tmpOper[i] = 0.0;
	     break;
	   }
	 }
       }

       ObjectMatrix<Matrix*> matRef;
       for (int i = 0; i < tmpState.quanta.size (); ++i)
       {
	 if (tmpOper[i].Ncols() == 0) continue;
	 std::vector<int> dum (1); 

	 int rows = tmpState.oldToNewState[i].size ();
	 int cols = dum.size ();
	 matRef.ReSize (rows, cols);
	 //cout << "***** "<<i<<"  "<<rows<<"  "<<tmpState.oldToNewState[0][0]<<endl;
	 std::vector<Matrix> dummyMats(rows);
	 for (int x = 0; x < rows; ++x) 
	   matRef (x,0) = &rotation[tmpState.oldToNewState [i][x]];

	 //cout << tmpOper[i]<<endl;
	 CatenateProduct (matRef, tmpOper[i]);
       }
       rotation = tmpOper;
    }
  catch (Exception)
    {
      Exception::what ();
      abort ();
    }


}

int SpinAdapted::PickSingleNumberAtRandom(std::vector<double>& d)
{ 
  int totalLength = d.size();

  vector<double> dCumulative(totalLength, 0.0);
  int index = 0;
  double oldCumulative = 0.0;
  for (int i=0; i<d.size(); i++) {
    dCumulative[i] = fabs(d[i]) + oldCumulative;
    oldCumulative = dCumulative[i];
  }

  double random = ((double) rand() / (RAND_MAX));

  index=-1;
  for (int i=0; i<totalLength; i++)
    if (dCumulative[i]/oldCumulative > random) {
      index = i;
      break;
    }

  if (index == -1)
    index = totalLength-1;
       
  int counter = 0;
  for (int i=0; i<d.size(); i++) {
    if (i == index)
      d[i] = 1.0*((d[i] > 0) ? 1 : ((d[i] < 0) ? -1 : 0));
    else
      d[i] = 0.0;
  }
  return index;
}

int SpinAdapted::PickSingleVectorAtRandom(std::vector<DiagonalMatrix>& d)
{ 
  int totalLength = 0;
  for (int i=0; i<d.size(); i++)
    totalLength += d[i].Ncols();

  vector<double> dCumulative(totalLength, 0.0);
  int index = 0;
  double oldCumulative = 0.0;
  for (int i=0; i<d.size(); i++) {
    for (int j=0; j<d[i].Ncols(); j++) {
      dCumulative[index] = fabs(d[i](j+1)) + oldCumulative;
      oldCumulative = dCumulative[index];
      index++;
    }
  }

  double random = ((double) rand() / (RAND_MAX));

  index=-1;
  for (int i=0; i<totalLength; i++)
    if (dCumulative[i]/oldCumulative > random) {
      index = i;
      break;
    }

  if (index == -1)
    index = totalLength-1;
       
  int counter = 0;
  for (int i=0; i<d.size(); i++) 
    for (int j=0; j<d[i].Ncols(); j++) {
      if (counter == index)
	d[i](j+1) = 1.0;
      else
	d[i](j+1) = 0.0;
      counter++;
    }
  return index;
}

void SpinAdapted::SaveRotationMatrix (const std::vector<int>& sites, const std::vector<Matrix>& m1, int state)
{
  dmrginp.diskwo->start();
  Timer disktimer;
  int rank = mpigetrank();
  if (rank == 0)
    {

      char file [5000];
      int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
      if (state == -1)
	sprintf (file, "%s%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state_average.tmp");
      else
	sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state",state, ".tmp");
      p1out << "\t\t\t Saving Rotation Matrix :: " << file << endl;
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_mat(ofs);
      save_mat << m1;
      ofs.close();
    }
  dmrginp.diskwo->stop();
}

void SpinAdapted::LoadRotationMatrix (const std::vector<int>& sites, std::vector<Matrix>& m1, int state)
{
  dmrginp.diskwi->start();
  Timer disktimer;
  int rank = mpigetrank();
  if (rank == 0)
  {
    char file [5000];
    int first = min(sites[0], *sites.rbegin()), last = max(sites[0], *sites.rbegin());
    if(state == -1)
      sprintf (file, "%s%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state_average.tmp");
    else
      sprintf (file, "%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/Rotation-", first, "-", last, ".", mpigetrank(),".state",state, ".tmp");
    p1out << "\t\t\t Loading Rotation Matrix :: " << file << endl;
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> m1;
    ifs.close();
  }
  dmrginp.diskwi->stop();
}

void SpinAdapted::allocate(const StateInfo& row, const StateInfo& col, std::vector<Matrix>& rotations)
{
  rotations.resize(row.quanta.size());
  for (int i=0; i<row.quanta.size(); i++) {
    int nrows = row.quantaStates[i];
    int ncols = 0;
    for (int j=0; j<col.quanta.size(); j++) 
      if (col.quanta[j] == row.quanta[i])
	ncols += col.quantaStates[j];

    rotations[i].ReSize(nrows, ncols);
  }
}


bool SpinAdapted::can_connect(int n, int spin, int right_block_size)
{
  for(int alpha=0;alpha<=min(right_block_size/2,dmrginp.total_particle_number()-n);++alpha)
  {
    int beta = dmrginp.total_particle_number() - n - alpha;
    if(dmrginp.total_spin_number().getirrep() - (alpha - beta) == spin && beta <= right_block_size/2)
      return true;
  }
  return false;
}

int Binom(int n, int k)
{
  vector<int> b(n+1);
  b[0] = 1;
  for(int i=1;i<=n;++i)
  {
    b[i] = 1;
    for(int j=i-1;j>0;--j)
      if(INT_MAX - b[j-1] > 0)
        b[j] += b[j-1];
      else
        b[j] = INT_MAX;
  }
  return b[k];
}

int SpinAdapted::get_total_states(const int &this_size, const int &other_size)
{
  int maxj = this_size/2+1;
  vector<int> statesWithJ(maxj, 0);
  vector<int> temp(maxj, 0);

  statesWithJ[0] = 2; statesWithJ[1] = 1;
  for (int i=1; i<this_size/2; i++)
  {
    temp = statesWithJ;
    for (int j=0; j<maxj; j++)
      statesWithJ[j] =0;
    for (int j=0; j<maxj; j++)
    {
      if (INT_MAX - 2*temp[j] > statesWithJ[j])
	statesWithJ[j] += 2*temp[j];
      else statesWithJ[j] = INT_MAX;

      if (j+1 <maxj)
      {
	if (INT_MAX - temp[j] > statesWithJ[j+1])
	  statesWithJ[j+1] += temp[j];
	else
	  statesWithJ[j+1] = INT_MAX;
      }
      if (j-1 >= 0)
      {
	if (INT_MAX - temp[j] > statesWithJ[j-1])
	  statesWithJ[j-1] += temp[j];
	else
	  statesWithJ[j-1] = INT_MAX;
      }
    }
  }
  double retval = 0;
  for (int i=0; i<maxj; i++)
  {
    if (INT_MAX - retval > statesWithJ[i])
	retval += statesWithJ[i];
    else
      retval = INT_MAX;
  }
  return retval;
}


double SpinAdapted::assign_matrix_by_dm(std::vector<Matrix>& rotatematrix, std::vector<DiagonalMatrix>& eigenmatrix, StackSparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta, int totalstatesbydm, int totalstatesbyquanta, int left_block_size, int right_block_size)
{
  const int min_states = totalstatesbydm;

  

  p2out << " \t\t\t assigning a total of " << min_states << " states using the dm alone " << endl;
  double totalnorm = 0.;
  rotatematrix.resize(eigenmatrix.size());

  for (int i = 0; i < totalstatesbydm; ++i)
    {
      int q = inorderwts[i].first;
      int qs = inorderwts[i].second;

      if( eigenmatrix[q].element(qs, qs) > 1.e-13)
      {
        if (rotatematrix[q].Ncols() == 0)
        {
	  rotatematrix[q].ReSize(transformmatrix(q,q).Nrows(), 1);
	  for (int i=0; i<transformmatrix(q,q).Nrows(); i++)
	    rotatematrix[q](i+1,1) = transformmatrix(q,q)(i+1, qs+1);
        }
        else
        {
	  Matrix bkp = rotatematrix[q];
	  rotatematrix[q].ReSize(bkp.Nrows(), bkp.Ncols()+1);

	  for (int i=0; i<bkp.Nrows(); i++)
	    for (int j=0; j<bkp.Ncols(); j++)
	      rotatematrix[q](i+1, j+1) = bkp(i+1, j+1);
	  
	  for (int i=0; i<transformmatrix(q,q).Nrows(); i++)
	    rotatematrix[q](i+1,bkp.Ncols()+1) = transformmatrix(q,q)(i+1, qs+1);
        }
        vector<int>::iterator findit = find(wtsbyquanta[q].begin(), wtsbyquanta[q].end(), qs);
        if (findit == wtsbyquanta[q].end()) { pout << " error in assign matrix " << endl; abort(); }
        wtsbyquanta[q].erase(findit);
        totalnorm += eigenmatrix[q].element(qs, qs);
      }
    }

  p2out << " \t\t\t assigning a total of " << totalstatesbyquanta << " states using quanta selection " << " for a norm of " << totalnorm << endl;

  int assignedbyq = 0;
  
  int nquanta = rotatematrix.size();
 
  int totalstatesleft = 0;
  for (int i = 0; i < nquanta; ++i)
    {
      totalstatesleft += wtsbyquanta[i].size();
    }

  p2out << " \t\t\t a total of " << totalstatesleft << " to be assigned " << endl;
  
  // now sort quanta in order of importance
  vector<double> totalquantaweights(nquanta);
  for (int q = 0; q < totalquantaweights.size(); ++q)
  {
    for (int qs = 0; qs < wtsbyquanta[q].size(); ++qs)
      totalquantaweights[q] += eigenmatrix[q].element(qs, qs);
  }
  vector<int> inorderquanta(nquanta);
  sort_data_to_indices(totalquantaweights, inorderquanta);
  reverse(inorderquanta.begin(), inorderquanta.end());  


  // reorder modified wtsbyquanta into a usable form
  vector<pair<int, int> > linearwtsbyquanta;
  int qspointer = 0;
  
  while(totalstatesleft)
  {
    for (int i = 0; i < nquanta; ++i)
	  {
	    int q = inorderquanta[i];
	    if (qspointer < wtsbyquanta[q].size())
	    {
	      linearwtsbyquanta.push_back(make_pair(q, wtsbyquanta[q][qspointer]));
	      --totalstatesleft;
	    }
	  }
    ++qspointer;
  }

  for (int i = 0; i < totalstatesbyquanta; ++i)
  {
    int q = linearwtsbyquanta[i].first;
    int qs = linearwtsbyquanta[i].second;
    if( eigenmatrix[q].element(qs, qs) > 1.e-13)
    {
      if (rotatematrix[q].Ncols() == 0)
      {
	rotatematrix[q].ReSize(transformmatrix(q,q).Nrows(), 1);
	for (int i=0; i<transformmatrix(q,q).Nrows(); i++)
	  rotatematrix[q](i+1,1) = transformmatrix(q,q)(i+1, qs+1);
	//rotatematrix[q] = transformmatrix(q, q).Column(qs + 1);
      }
      else
      {
	Matrix bkp = rotatematrix[q];
	rotatematrix[q].ReSize(rotatematrix[q].Nrows(), rotatematrix[q].Ncols()+1);
	for (int i=0; i<bkp.Nrows(); i++)
	  for (int j=0; j<bkp.Ncols(); j++)
	    rotatematrix[q](i+1, j+1) = bkp(i+1, j+1);
	
	for (int i=0; i<transformmatrix(q,q).Nrows(); i++)
	  rotatematrix[q](i+1,bkp.Ncols()+1) = transformmatrix(q,q)(i+1, qs+1);
	//rotatematrix[q] |= transformmatrix(q, q).Column(qs + 1);      
      }
      vector<int>::iterator findit = find(wtsbyquanta[q].begin(), wtsbyquanta[q].end(), qs);
      if (findit == wtsbyquanta[q].end()) { pout << " error in assign matrix " << endl; abort(); }
      wtsbyquanta[q].erase(findit);
      totalnorm += eigenmatrix[q].element(qs, qs);
    }
  }
  
  double norm = 0.;
  for(int i=0;i<eigenmatrix.size();++i)
    for(int j=0;j<eigenmatrix[i].Nrows();++j)
      norm += eigenmatrix[i].element(j, j);
  p2out << " \t\t\t total norm: " << norm <<"  norm after truncation: "<<totalnorm<< endl;

  return norm-totalnorm;

  //return (1. - totalnorm/norm);
}


void SpinAdapted::diagonalise_dm(StackSparseMatrix& tracedMatrix, std::vector<DiagonalMatrix>& eigenMatrix)
{
  int nquanta = tracedMatrix.nrows();
  eigenMatrix.resize(nquanta);
  vector<double> totalquantaweights(nquanta);
  for (int tQ = 0; tQ < nquanta; ++tQ)
    {
      int nStates = tracedMatrix.operator_element(tQ, tQ).Nrows ();
      DiagonalMatrix weights (nStates);
#ifdef USELAPACK
      diagonalise(tracedMatrix.operator_element(tQ, tQ), weights);
#else
      SymmetricMatrix dM (nStates);
      dM << tracedMatrix.operator_element(tQ,tQ);
      EigenValues (dM, weights, transformMatrix.operator_element(tQ,tQ));
#endif
      for(int i=0;i<weights.Nrows();++i)
        if(weights.element(i,i) < 1.e-14)
          weights.element(i,i) = 0.;
      eigenMatrix [tQ] = weights;
    }  
}


void SpinAdapted::svd_densitymat(StackSparseMatrix& wavefn, std::vector<Matrix>& U, std::vector<DiagonalMatrix>& eigenMatrix, 
				 std::vector<Matrix>& V) {

  int nquanta = wavefn.nrows();
  eigenMatrix.resize(nquanta);
  V.resize(nquanta);
  U.resize(nquanta);
  vector<double> totalquantaweights(nquanta);
  for (int tQ = 0; tQ < nquanta; ++tQ) 
  for (int rQ = 0; rQ < wavefn.ncols(); ++rQ) 
  if  (wavefn.allowed(tQ, rQ)) {
    int nStates = wavefn.operator_element(tQ, rQ).Nrows ();
    Matrix wftQ; copy(wavefn.operator_element(tQ, rQ), wftQ);
#ifdef USELAPACK
    svd(wftQ, eigenMatrix[tQ], U[tQ], V[tQ]);
#else
    cerr << "not working without lapack"<<endl;
    exit(0);
#endif
  }
}

void SpinAdapted::svd_densitymat(StackSparseMatrix& tracedMatrix, std::vector<DiagonalMatrix>& eigenMatrix) {
  // SVD of matrix M=(A,B,C)=USV^T
  // since MM^T=AA^T+BB^T+CC^T=USS^TU^T, we don't have to explicitly construct M
  int nquanta = tracedMatrix.nrows();
  eigenMatrix.resize(nquanta);
  vector<double> totalquantaweights(nquanta);
  for (int tQ = 0; tQ < nquanta; ++tQ) {
    int nStates = tracedMatrix.operator_element(tQ, tQ).Nrows ();
    DiagonalMatrix weights(nStates);
    std::vector<double> data(nStates*nStates, 0.0);
    StackMatrix M(&data[0], nStates, nStates);
    for (int sQ = 0; sQ < nquanta; ++sQ)
      if (tracedMatrix.allowed(tQ, sQ))
	MatrixMultiply (tracedMatrix.operator_element(tQ, tQ), 'n', tracedMatrix.operator_element(tQ, tQ), 't',
			M, 1.0);

#ifdef USELAPACK
    diagonalise(M, weights);
#else
    SymmetricMatrix dM (nStates);
    dM << M;
    EigenValues(dM, weights, transformMatrix.operator_element(tQ,tQ));
#endif

    copy(M, tracedMatrix.operator_element(tQ, tQ));

    for (int i=0; i<weights.Nrows(); ++i) {
      if(weights.element(i,i) < 1.e-28)
        weights.element(i,i) = 0.;
      else
        weights.element(i,i) = sqrt(weights.element(i,i));
    }
    eigenMatrix[tQ] = weights;    
  }
}

void SpinAdapted::sort_weights(std::vector<DiagonalMatrix>& eigenMatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& weightsbyquanta)
{
  // first sort weights with a multimap
  multimap<double, pair<int,int> > weightmap;
  int nquanta = eigenMatrix.size();
  vector<double> totalquantaweights(nquanta);

  for (int q = 0; q < nquanta; ++q)
  {
    //if (q == hfQuantaindex) continue; //the hartree fock quanta are always included first during warmup
    for (int qs = 0; qs < eigenMatrix[q].Nrows(); ++qs)
      {
	weightmap.insert (pair <double, pair<int,int> > (eigenMatrix[q].element(qs, qs), pair<int,int> (q, qs)));
	totalquantaweights[q] += eigenMatrix[q].element(qs, qs);

      }
  }

  multimap<double, pair<int,int> >::reverse_iterator w = weightmap.rbegin();


  // now put all the sorted indices in
  while (w != weightmap.rend())
    {
      inorderwts.push_back(make_pair(w->second.first, w->second.second));
      ++w;
    }

  // sort quantas by weight
  weightsbyquanta.resize(nquanta);
  for (int q = 0; q < nquanta; ++q)
    for (int qs = 0; qs < eigenMatrix[q].Nrows(); ++qs)
      weightsbyquanta[q].push_back(eigenMatrix[q].Nrows() - qs - 1);
}

