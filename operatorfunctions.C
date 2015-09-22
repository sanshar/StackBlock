/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "StackMatrix.h"
#include "timer.h"
#include "Stackspinblock.h"
#include "MatrixBLAS.h"
#include <math.h>
#include "global.h"
#include <omp.h>
#include <iostream>
#include <map>
#include <vector>
#include <newmat.h>
#include "StateInfo.h"
#include "StackBaseOperator.h"
#include "operatorfunctions.h"
#include "Stackspinblock.h"
#include "Stackwavefunction.h"
#include "couplingCoeffs.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pario.h"


  //TENSOR TRACE A x I  ->  C
void SpinAdapted::operatorfunctions::TensorTraceElement(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, StackMatrix& cel, int cq, int cqprime, double scale)
{
  if (fabs(scale) < TINY) return;
  assert(c.allowed(cq, cqprime));

  if (cstateinfo->hasCollectedQuanta) {
    int aq, aqprime, bq, bqprime, bstates;
    const char conjC = (ablock == cblock->get_leftBlock()) ? 'n' : 't';
    
    const std::vector<int> oldToNewI = cstateinfo->oldToNewState.at(cq);
    const std::vector<int> oldToNewJ = cstateinfo->oldToNewState.at(cqprime);
    
    const StateInfo* rS = cstateinfo->rightStateInfo, *lS = cstateinfo->leftStateInfo;
    int rowstride =0, colstride = 0;

    for (int oldi =0; oldi < oldToNewI.size(); oldi++) {
      colstride = 0;
      for (int oldj = 0; oldj < oldToNewJ.size(); oldj++)
	{
	  if (conjC == 'n')
	    {
	      aq = cstateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	      aqprime = cstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	      bq = cstateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	      bqprime = cstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	      bstates = cstateinfo->rightStateInfo->getquantastates(bq);
	    }
	  else 
	    {
	      aq = cstateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	      aqprime = cstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	      bq = cstateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	      bqprime = cstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	      bstates = cstateinfo->leftStateInfo->getquantastates(bq);
	    }
	  
	  if (a.allowed(aq, aqprime) && (bq == bqprime))
	    {
	      DiagonalMatrix unitMatrix(bstates);
	      unitMatrix = 1.;
	      
	      Matrix unity(bstates, bstates);
	      unity = unitMatrix;
	      
	      if (conjC == 'n')
		{
		  double scaleb = dmrginp.get_ninej()(lS->quanta[aqprime].get_s().getirrep() , rS->quanta[bqprime].get_s().getirrep(), cstateinfo->quanta[cqprime].get_s().getirrep(), 
						      a.get_spin().getirrep(), 0, c.get_spin().getirrep(),
						      lS->quanta[aq].get_s().getirrep() , rS->quanta[bq].get_s().getirrep(), cstateinfo->quanta[cq].get_s().getirrep());
		  
		  scaleb *= Symmetry::spatial_ninej(lS->quanta[aqprime].get_symm().getirrep() , rS->quanta[bqprime].get_symm().getirrep(), cstateinfo->quanta[cqprime].get_symm().getirrep(), 
						    a.get_symm().getirrep(), 0, c.get_symm().getirrep(),
						    lS->quanta[aq].get_symm().getirrep() , rS->quanta[bq].get_symm().getirrep(), cstateinfo->quanta[cq].get_symm().getirrep());
		  
		  MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scale, unity, 'n', scaleb, 
				       cel, rowstride, colstride);
		}
	      else {
		double scaleb = dmrginp.get_ninej()(lS->quanta[bqprime].get_s().getirrep(), rS->quanta[aqprime].get_s().getirrep() , cstateinfo->quanta[cqprime].get_s().getirrep(), 
						    0, a.get_spin().getirrep(), c.get_spin().getirrep(),
						    lS->quanta[bq].get_s().getirrep(), rS->quanta[aq].get_s().getirrep() , cstateinfo->quanta[cq].get_s().getirrep());
		scaleb *= Symmetry::spatial_ninej(lS->quanta[bqprime].get_symm().getirrep() , rS->quanta[aqprime].get_symm().getirrep(), cstateinfo->quanta[cqprime].get_symm().getirrep(), 
						  0, a.get_symm().getirrep(), c.get_symm().getirrep(),
						  lS->quanta[bq].get_symm().getirrep() , rS->quanta[aq].get_symm().getirrep(), cstateinfo->quanta[cq].get_symm().getirrep());
		if (a.get_fermion() && IsFermion (cstateinfo->leftStateInfo->quanta[bqprime]) ) scaleb *= -1.;

		MatrixTensorProduct (unity, 'n', scaleb, a.operator_element(aq, aqprime), a.conjugacy(), scale, 
				     cel, rowstride, colstride);
	      }
	    }
	  colstride += cstateinfo->unCollectedStateInfo->quantaStates[ oldToNewJ[oldj] ];
	  
	}
      rowstride += cstateinfo->unCollectedStateInfo->quantaStates[ oldToNewI[oldi] ];
      
    }
  }
  else {
    int aq, aqprime, bq, bqprime, bstates;
    const char conjC = (ablock == cblock->get_leftBlock()) ? 'n' : 't';
    
    //const std::vector<int> oldToNewI = cstateinfo->oldToNewState.at(cq);
    //const std::vector<int> oldToNewJ = cstateinfo->oldToNewState.at(cqprime);
    
    const StateInfo* rS = cstateinfo->rightStateInfo, *lS = cstateinfo->leftStateInfo;
    int rowstride =0, colstride = 0;
    
    if (conjC == 'n')
    {
      aq = cstateinfo->leftUnMapQuanta[ cq];
      aqprime = cstateinfo->leftUnMapQuanta[ cqprime ];
      bq = cstateinfo->rightUnMapQuanta[ cq ];
      bqprime = cstateinfo->rightUnMapQuanta[ cqprime ];
      bstates = cstateinfo->rightStateInfo->getquantastates(bq);
    }
    else 
    {
      aq = cstateinfo->rightUnMapQuanta[ cq ];
      aqprime = cstateinfo->rightUnMapQuanta[ cqprime ];
      bq = cstateinfo->leftUnMapQuanta[ cq ];
      bqprime = cstateinfo->leftUnMapQuanta[ cqprime ];
      bstates = cstateinfo->leftStateInfo->getquantastates(bq);
    }
    
    if (a.allowed(aq, aqprime) && (bq == bqprime))
    {
      DiagonalMatrix unitMatrix(bstates);
      unitMatrix = 1.;
      
      Matrix unity(bstates, bstates);
      unity = unitMatrix;
      
      if (conjC == 'n')
      {
	double scaleb = dmrginp.get_ninej()(lS->quanta[aqprime].get_s().getirrep() , rS->quanta[bqprime].get_s().getirrep(), cstateinfo->quanta[cqprime].get_s().getirrep(), 
					    a.get_spin().getirrep(), 0, c.get_spin().getirrep(),
					    lS->quanta[aq].get_s().getirrep() , rS->quanta[bq].get_s().getirrep(), cstateinfo->quanta[cq].get_s().getirrep());
	
	scaleb *= Symmetry::spatial_ninej(lS->quanta[aqprime].get_symm().getirrep() , rS->quanta[bqprime].get_symm().getirrep(), cstateinfo->quanta[cqprime].get_symm().getirrep(), 
					  a.get_symm().getirrep(), 0, c.get_symm().getirrep(),
					  lS->quanta[aq].get_symm().getirrep() , rS->quanta[bq].get_symm().getirrep(), cstateinfo->quanta[cq].get_symm().getirrep());
	
	MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scale, unity, 'n', scaleb, 
			     cel, rowstride, colstride);
      }
      else {
	double scaleb = dmrginp.get_ninej()(lS->quanta[bqprime].get_s().getirrep(), rS->quanta[aqprime].get_s().getirrep() , cstateinfo->quanta[cqprime].get_s().getirrep(), 
					    0, a.get_spin().getirrep(), c.get_spin().getirrep(),
					    lS->quanta[bq].get_s().getirrep(), rS->quanta[aq].get_s().getirrep() , cstateinfo->quanta[cq].get_s().getirrep());
	scaleb *= Symmetry::spatial_ninej(lS->quanta[bqprime].get_symm().getirrep() , rS->quanta[aqprime].get_symm().getirrep(), cstateinfo->quanta[cqprime].get_symm().getirrep(), 
					  0, a.get_symm().getirrep(), c.get_symm().getirrep(),
					  lS->quanta[bq].get_symm().getirrep() , rS->quanta[aq].get_symm().getirrep(), cstateinfo->quanta[cq].get_symm().getirrep());
	if (a.get_fermion() && IsFermion (cstateinfo->leftStateInfo->quanta[bqprime]) ) scaleb *= -1.;
	MatrixTensorProduct (unity, 'n', scaleb, a.operator_element(aq, aqprime), a.conjugacy(), scale, 
			     cel, rowstride, colstride);
      }
    }
  }
}
  
void SpinAdapted::operatorfunctions::TensorTrace (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix* cDiagonal, Real scale)
{
  if (fabs(scale) < TINY) return;
  assert (a.get_initialised());
  
  
  const char conjC = (ablock == cblock->get_leftBlock()) ? 'n' : 't';
  
  const int aSz = ablock->get_stateInfo().quanta.size ();
  const int bSz = (conjC == 'n') ? cblock->get_stateInfo().rightStateInfo->quanta.size () : cblock->get_stateInfo().leftStateInfo->quanta.size ();
  const StateInfo& s = cblock->get_stateInfo();
  
  const StateInfo* lS = s.leftStateInfo, *rS = s.rightStateInfo;

  int OMPRANK = omprank; 

  for (int aQ = 0; aQ < aSz; ++aQ)
    if (a.allowed(aQ, aQ))
      for (int bQ = 0; bQ < bSz; ++bQ)
	if (s.allowedQuanta (aQ, bQ, conjC))
	  {
	    int cQ = s.quantaMap (aQ, bQ, conjC)[0];
	    for (int cQState = 0; cQState < s.quantaStates [cQ]; ++cQState)
	      {
		Real scaleB = 1.0;
		int aQState;
		int bQState;
		
		if (conjC == 'n')
		  {
		    s.UnMapQuantumState (cQState, s.rightStateInfo->quantaStates [bQ], aQState, bQState);
		    scaleB *= dmrginp.get_ninej()(lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
						  a.get_spin().getirrep(), 0, 0,
						  lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
		    scaleB *= Symmetry::spatial_ninej(lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
						      a.get_symm().getirrep(), 0, 0,
						      lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		  }
		else
		  {
		    scaleB *= dmrginp.get_ninej()(lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
						  0, a.get_spin().getirrep(), 0,
						  lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
		    scaleB *= Symmetry::spatial_ninej(lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
						      0, a.get_symm().getirrep(), 0,
						      lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		    if (a.get_fermion()&& IsFermion(lS->quanta[bQ])) scaleB *= -1.0;
		    s.UnMapQuantumState (cQState, s.rightStateInfo->quantaStates [aQ], bQState, aQState);
		  }
		long dindex = s.unBlockedIndex[cQ] + cQState + 1;
		cDiagonal[OMPRANK](dindex) += scale * scaleB * a.operator_element(aQ, aQ)(aQState + 1, aQState + 1); 
		
	      }
	  }
  
}


void SpinAdapted::operatorfunctions::TensorProduct (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock* cblock, const StateInfo* cstateinfo, DiagonalMatrix* cDiagonal, double scale)
{
  if (fabs(scale) < TINY) return;
  const int aSz = a.nrows();
  const int bSz = b.nrows();
  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';
  const StackSpinBlock* bblock = (cblock->get_leftBlock() == ablock) ? cblock->get_rightBlock() : cblock->get_leftBlock();
  const StateInfo& s = cblock->get_stateInfo();
  const StateInfo* lS = s.leftStateInfo, *rS = s.rightStateInfo;

  int OMPRANK = omprank;

  for (int aQ = 0; aQ < aSz; ++aQ)
    if (a.allowed(aQ, aQ))
      for (int bQ = 0; bQ < bSz; ++bQ)
	if (b.allowed(bQ, bQ))
	  if (s.allowedQuanta (aQ, bQ, conjC))
	  {
	    int cQ = s.quantaMap (aQ, bQ, conjC)[0];
	    Real scaleA = scale;
	    Real scaleB = 1;
	    if (conjC == 'n')
	      {
		scaleB *= dmrginp.get_ninej()(lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
					      a.get_spin().getirrep(), b.get_spin().getirrep(), 0,
					      lS->quanta[aQ].get_s().getirrep() , rS->quanta[bQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
		scaleB *= Symmetry::spatial_ninej(lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
				     a.get_symm().getirrep(), b.get_symm().getirrep(), 0,
				     lS->quanta[aQ].get_symm().getirrep() , rS->quanta[bQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		
		if (b.get_fermion() && IsFermion (lS->quanta [aQ])) scaleB *= -1.0;
		for (int aQState = 0; aQState < lS->quantaStates[aQ] ; aQState++)
		  MatrixDiagonalScale(a.operator_element(aQ, aQ)(aQState+1, aQState+1)*scaleA*scaleB, b.operator_element(bQ, bQ), 
				      cDiagonal[OMPRANK].Store()+s.unBlockedIndex[cQ]+aQState*rS->quantaStates[bQ]);

	      }
	    else
	      {
		scaleB *= dmrginp.get_ninej()(lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep(), 
					      b.get_spin().getirrep(), a.get_spin().getirrep(), 0,
					      lS->quanta[bQ].get_s().getirrep() , rS->quanta[aQ].get_s().getirrep(), cstateinfo->quanta[cQ].get_s().getirrep());
		scaleB *= Symmetry::spatial_ninej(lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep(), 
				     b.get_symm().getirrep(), a.get_symm().getirrep(), 0,
				     lS->quanta[bQ].get_symm().getirrep() , rS->quanta[aQ].get_symm().getirrep(), cstateinfo->quanta[cQ].get_symm().getirrep());
		
		if (a.get_fermion()&& IsFermion(lS->quanta[bQ])) scaleB *= -1.0;
		for (int bQState = 0; bQState < lS->quantaStates[bQ] ; bQState++)
		  MatrixDiagonalScale(b.operator_element(bQ, bQ)(bQState+1, bQState+1)*scaleA*scaleB, a.operator_element(aQ, aQ), 
				      cDiagonal[OMPRANK].Store()+s.unBlockedIndex[cQ]+bQState*rS->quantaStates[aQ]);
	      }
	  }
}


void SpinAdapted::operatorfunctions::TensorTrace(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, double scale, int num_thrds) {
  
  if (fabs(scale) < TINY) return;
  assert (a.get_initialised() && c.get_initialised());
  
  std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();

  int quanta_thrds = dmrginp.quanta_thrds();
#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int cq = nonZeroBlocks[index].first.first, cqprime = nonZeroBlocks[index].first.second;
    TensorTraceElement(ablock, a, cblock, cstateinfo, c, nonZeroBlocks[index].second, cq, cqprime, scale);
  }

}



void SpinAdapted::operatorfunctions::TensorProductElement(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, StackMatrix& cel, int cq, int cqprime, double scale)
{
  //cstateinfo is not used
  //This function can be used for different bra and ket stateinfo.
  if (fabs(scale) < TINY) return;
  assert (a.get_initialised());
  assert (b.get_initialised());
  assert (c.get_initialised());
 

  if (cstateinfo->hasCollectedQuanta) {
    const StateInfo *ketstateinfo = &cblock->get_ketStateInfo(), 
      *brastateinfo = &cblock->get_braStateInfo();
    
    const StackSpinBlock* bblock = (cblock->get_leftBlock() == ablock) ? cblock->get_rightBlock() : cblock->get_leftBlock();
    
    const std::vector<int>& oldToNewI = brastateinfo->oldToNewState.at(cq);
    const std::vector<int>& oldToNewJ = ketstateinfo->oldToNewState.at(cqprime);
    
    const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';
    
    const StateInfo* lbraS = brastateinfo->leftStateInfo, *rbraS = brastateinfo->rightStateInfo;
    const StateInfo* lketS = ketstateinfo->leftStateInfo, *rketS = ketstateinfo->rightStateInfo;
    int rowstride = 0, colstride = 0;
    
    int aq, aqprime, bq, bqprime;
    
    //pout << "old to new size "<<oldToNewI.size()<<" "<<oldToNewJ.size()<<endl;
    for (int oldi =0; oldi < oldToNewI.size(); oldi++) {
      colstride = 0;
      for (int oldj = 0; oldj < oldToNewJ.size(); oldj++){
	if (conjC == 'n'){
	  aq = brastateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	  aqprime = ketstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	  bq = brastateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	  bqprime = ketstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	}
	else {
	  aq = brastateinfo->rightUnMapQuanta[ oldToNewI[oldi] ];
	  aqprime = ketstateinfo->rightUnMapQuanta[ oldToNewJ[oldj] ];
	  bq = brastateinfo->leftUnMapQuanta[ oldToNewI[oldi] ];
	  bqprime = ketstateinfo->leftUnMapQuanta[ oldToNewJ[oldj] ];
	}
  
	Real scaleA = scale;
	Real scaleB = 1.0;
	if (a.allowed(aq, aqprime) && b.allowed(bq, bqprime)){
	  if (conjC == 'n') {
	    scaleB = dmrginp.get_ninej()(lketS->quanta[aqprime].get_s().getirrep() , rketS->quanta[bqprime].get_s().getirrep(), ketstateinfo->quanta[cqprime].get_s().getirrep(), 
					 a.get_spin().getirrep(), b.get_spin().getirrep(), c.get_spin().getirrep(),
					 lbraS->quanta[aq].get_s().getirrep() , rbraS->quanta[bq].get_s().getirrep(), brastateinfo->quanta[cq].get_s().getirrep());
	    scaleB *= Symmetry::spatial_ninej(lketS->quanta[aqprime].get_symm().getirrep() , rketS->quanta[bqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
					      a.get_symm().getirrep(), b.get_symm().getirrep(), c.get_symm().getirrep(),
					      lbraS->quanta[aq].get_symm().getirrep() , rbraS->quanta[bq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	    scaleB *= b.get_scaling(rbraS->quanta[bq], rketS->quanta[bqprime]);
	    scaleA *= a.get_scaling(lbraS->quanta[aq], lketS->quanta[aqprime]);
	    if (b.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta [aqprime])) scaleB *= -1;
	    MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scaleA, 
				 b.operator_element(bq, bqprime), b.conjugacy(), scaleB, cel,rowstride, colstride);
	  }
	  else {
	    scaleB = dmrginp.get_ninej()(lketS->quanta[bqprime].get_s().getirrep(), rketS->quanta[aqprime].get_s().getirrep() , ketstateinfo->quanta[cqprime].get_s().getirrep(), 
					 b.get_spin().getirrep(), a.get_spin().getirrep(), c.get_spin().getirrep(),
					 lbraS->quanta[bq].get_s().getirrep(), rbraS->quanta[aq].get_s().getirrep() , brastateinfo->quanta[cq].get_s().getirrep());
	    scaleB *= Symmetry::spatial_ninej(lketS->quanta[bqprime].get_symm().getirrep() , rketS->quanta[aqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
					      b.get_symm().getirrep(), a.get_symm().getirrep(), c.get_symm().getirrep(),
					      lbraS->quanta[bq].get_symm().getirrep() , rbraS->quanta[aq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	    scaleB *= b.get_scaling(lbraS->quanta[bq], lketS->quanta[bqprime]);
	    scaleA *= a.get_scaling(rbraS->quanta[aq], rketS->quanta[aqprime]);
	    if (a.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta[bqprime]) ) scaleB *= -1.;
	    
	    MatrixTensorProduct (b.operator_element(bq, bqprime), b.conjugacy(), scaleB, 
				 a.operator_element(aq, aqprime), a.conjugacy(), scaleA, cel, rowstride, colstride);
	  }
	}
	colstride += ketstateinfo->unCollectedStateInfo->quantaStates[ oldToNewJ[oldj] ];
	
      }
      rowstride += brastateinfo->unCollectedStateInfo->quantaStates[ oldToNewI[oldi] ];
      
    }
    
    
  }
  else {
    const StateInfo *ketstateinfo = &cblock->get_ketStateInfo(), 
      *brastateinfo = &cblock->get_braStateInfo();
    
    const StackSpinBlock* bblock = (cblock->get_leftBlock() == ablock) ? cblock->get_rightBlock() : cblock->get_leftBlock();
    
    //const std::vector<int>& oldToNewI = brastateinfo->oldToNewState.at(cq);
    //const std::vector<int>& oldToNewJ = ketstateinfo->oldToNewState.at(cqprime);
    
    const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';
    
    const StateInfo* lbraS = brastateinfo->leftStateInfo, *rbraS = brastateinfo->rightStateInfo;
    const StateInfo* lketS = ketstateinfo->leftStateInfo, *rketS = ketstateinfo->rightStateInfo;
    int rowstride = 0, colstride = 0;
    
    int aq, aqprime, bq, bqprime;
    
    //pout << "old to new size "<<oldToNewI.size()<<" "<<oldToNewJ.size()<<endl;
    if (conjC == 'n'){
      aq = brastateinfo->leftUnMapQuanta[ cq ];
      aqprime = ketstateinfo->leftUnMapQuanta[ cqprime ];
      bq = brastateinfo->rightUnMapQuanta[ cq ];
      bqprime = ketstateinfo->rightUnMapQuanta[ cqprime ];
    }
    else {
      aq = brastateinfo->rightUnMapQuanta[ cq ];
      aqprime = ketstateinfo->rightUnMapQuanta[ cqprime ];
      bq = brastateinfo->leftUnMapQuanta[ cq ];
      bqprime = ketstateinfo->leftUnMapQuanta[ cqprime ];
    }
    
    Real scaleA = scale;
    Real scaleB = 1.0;
    if (a.allowed(aq, aqprime) && b.allowed(bq, bqprime)){
      if (conjC == 'n') {
	scaleB = dmrginp.get_ninej()(lketS->quanta[aqprime].get_s().getirrep() , rketS->quanta[bqprime].get_s().getirrep(), ketstateinfo->quanta[cqprime].get_s().getirrep(), 
				     a.get_spin().getirrep(), b.get_spin().getirrep(), c.get_spin().getirrep(),
				     lbraS->quanta[aq].get_s().getirrep() , rbraS->quanta[bq].get_s().getirrep(), brastateinfo->quanta[cq].get_s().getirrep());
	scaleB *= Symmetry::spatial_ninej(lketS->quanta[aqprime].get_symm().getirrep() , rketS->quanta[bqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
					  a.get_symm().getirrep(), b.get_symm().getirrep(), c.get_symm().getirrep(),
					  lbraS->quanta[aq].get_symm().getirrep() , rbraS->quanta[bq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	scaleB *= b.get_scaling(rbraS->quanta[bq], rketS->quanta[bqprime]);
	scaleA *= a.get_scaling(lbraS->quanta[aq], lketS->quanta[aqprime]);
	if (b.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta [aqprime])) scaleB *= -1;
	MatrixTensorProduct (a.operator_element(aq, aqprime), a.conjugacy(), scaleA, 
			     b.operator_element(bq, bqprime), b.conjugacy(), scaleB, cel,rowstride, colstride);
      }
      else {
	scaleB = dmrginp.get_ninej()(lketS->quanta[bqprime].get_s().getirrep(), rketS->quanta[aqprime].get_s().getirrep() , ketstateinfo->quanta[cqprime].get_s().getirrep(), 
				     b.get_spin().getirrep(), a.get_spin().getirrep(), c.get_spin().getirrep(),
				     lbraS->quanta[bq].get_s().getirrep(), rbraS->quanta[aq].get_s().getirrep() , brastateinfo->quanta[cq].get_s().getirrep());
	scaleB *= Symmetry::spatial_ninej(lketS->quanta[bqprime].get_symm().getirrep() , rketS->quanta[aqprime].get_symm().getirrep(), ketstateinfo->quanta[cqprime].get_symm().getirrep(), 
					  b.get_symm().getirrep(), a.get_symm().getirrep(), c.get_symm().getirrep(),
					  lbraS->quanta[bq].get_symm().getirrep() , rbraS->quanta[aq].get_symm().getirrep(), brastateinfo->quanta[cq].get_symm().getirrep());
	scaleB *= b.get_scaling(lbraS->quanta[bq], lketS->quanta[bqprime]);
	scaleA *= a.get_scaling(rbraS->quanta[aq], rketS->quanta[aqprime]);
	if (a.get_fermion() && IsFermion (ketstateinfo->leftStateInfo->quanta[bqprime]) ) scaleB *= -1.;
	
	MatrixTensorProduct (b.operator_element(bq, bqprime), b.conjugacy(), scaleB, 
			     a.operator_element(aq, aqprime), a.conjugacy(), scaleA, cel, rowstride, colstride);
      }
  
    }



  }
}




void SpinAdapted::operatorfunctions::TensorProduct (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, const StateInfo *cstateinfo, StackSparseMatrix& c, double scale, int num_thrds)
{
  if (fabs(scale) < TINY) return;
  int rows = c.nrows();
  int cols = c.ncols();

  int quanta_thrds = dmrginp.quanta_thrds();
#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int cq = 0; cq < rows; ++cq)
  for (int cqprime = 0; cqprime < cols; ++cqprime)
  if (c.allowed(cq, cqprime)) {
      TensorProductElement(ablock, a, b, cblock, cstateinfo, c, c.operator_element(cq, cqprime), cq, cqprime, scale);
  }


  /*
  std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();
#pragma omp parallel for schedule(dynamic)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int cq = nonZeroBlocks[index].first.first, cqprime = nonZeroBlocks[index].first.second;
    if (c.conjugacy() == 'n')
      TensorProductElement(ablock, a, b, cblock, cstateinfo, c, c.operator_element(cq, cqprime), cq, cqprime, scale);
    else
      TensorProductElement(ablock, a, b, cblock, cstateinfo, c, c.operator_element(cq, cqprime), cqprime, cq, scale);
  }
  */
}




void SpinAdapted::operatorfunctions::TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale, int num_thrds)
{
  // cannot be used for situation with different bra and ket
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;


  assert (cblock->get_leftBlock() == ablock || cblock->get_rightBlock() == ablock);
  if (cblock->get_leftBlock() == ablock)
  {
    const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = a.get_nonZeroBlocks();
    //#pragma omp parallel for schedule(dynamic)
    for (int index = 0; index<nonZeroBlocks.size(); index++) {
      int lQ = nonZeroBlocks[index].first.first, lQPrime = nonZeroBlocks[index].first.second;
      for (int rQ = 0; rQ < rightKetOpSz; ++rQ) 
	if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
	{
	  double fac=scale;
	  fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
				     a.get_spin().getirrep(), 0, a.get_spin().getirrep(),
				     lbraS->quanta[lQ].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
	  fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
					 a.get_symm().getirrep(), 0, a.get_symm().getirrep(),
					 lbraS->quanta[lQ].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
	  fac *= a.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	  MatrixMultiply (a.operator_element(lQ, lQPrime), a.conjugacy(), c.operator_element(lQPrime, rQ), c.conjugacy(),
			  v.operator_element(lQ, rQ), fac);
	}
      
    }
  }
  else
  {
    const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = a.get_nonZeroBlocks();
    //#pragma omp parallel for schedule(dynamic)
    for (int index = 0; index<nonZeroBlocks.size(); index++) {
      int rQ = nonZeroBlocks[index].first.first, rQPrime = nonZeroBlocks[index].first.second;
      for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime) 
	if (v.allowed(lQPrime, rQ) && c.allowed(lQPrime, rQPrime)) {
	  double fac = scale;
	  fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
				     0, a.get_spin().getirrep(), a.get_spin().getirrep(),
				     lketS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
	  fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
					 0, a.get_symm().getirrep(), a.get_symm().getirrep(),
					 lketS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
	  fac *= a.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	  double parity = a.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
	  
	  MatrixMultiply (c.operator_element(lQPrime, rQPrime), c.conjugacy(),
			  a.operator_element(rQ, rQPrime), TransposeOf(a.conjugacy()), v.operator_element(lQPrime, rQ), fac*parity);
	}
      
    }
  }
}


void SpinAdapted::operatorfunctions::TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale)
{
  long starttime = globaltimer.totalwalltime();

    // can be used for situation with different bra and ket
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *rbraS = cblock->get_braStateInfo().rightStateInfo;
  const StateInfo* lketS = cblock->get_ketStateInfo().leftStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';

  const StackSparseMatrix& leftOp = (conjC == 'n') ? a : b; // an ugly hack to support the release memory optimisation
  const StackSparseMatrix& rightOp = (conjC == 'n') ? b : a;
  const char leftConj = (conjC == 'n') ? a.conjugacy() : b.conjugacy();
  const char rightConj = (conjC == 'n') ? b.conjugacy() : a.conjugacy();

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[0].get_nonZeroBlocks();

  long maxlen = 0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
      if (maxlen < lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime))
	maxlen = lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime);

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int lQ = nonZeroBlocks[index].first.first, rQ = nonZeroBlocks[index].first.second;

    const std::vector<int>& colinds = rightOp.getActiveCols(rQ);
    for (int rrop=0; rrop <colinds.size(); rrop ++) {
      int rQPrime = colinds[rrop];
      
	const std::vector<int>& rowinds = c.getActiveRows(rQPrime);
	for (int l = 0; l < rowinds.size(); l++) {
	  int lQPrime = rowinds[l];
	  if (leftOp.allowed(lQ, lQPrime) ) {

	    StackMatrix m(dataArray[omprank], lketS->getquantastates(lQPrime), rbraS->getquantastates(rQ));
	    
	    double factor = scale*leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);	      
	    factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					  leftOp.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					  lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	    factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
					      leftOp.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
					      lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v[OMPRANK].get_symm().getirrep());
	    int parity = rightOp.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
	    factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);

	    MatrixMultiply (c.operator_element(lQPrime, rQPrime), 'n', rightOp.operator_element(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), 
			    m, 1.0, 0.);	      
	    MatrixMultiply (leftOp.operator()(lQ, lQPrime), leftConj, m, 'n',  v[OMPRANK].operator_element(lQ, rQ), factor*parity);

	  }
	}
    }
  }
  
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
   Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }

  

}


void SpinAdapted::operatorfunctions::TensorMultiplysplitLeft(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale)
{
  long starttime = globaltimer.totalwalltime();

    // can be used for situation with different bra and ket
  const int unCollectedleftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().unCollectedStateInfo->quanta.size ();
  const int unCollectedleftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().unCollectedStateInfo->quanta.size ();
  const int leftBraOpSz = cblock->get_leftBlock()->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int dotBraOpSz = cblock->get_leftBlock()->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int dotKetOpSz = cblock->get_leftBlock()->get_rightBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[0].get_nonZeroBlocks();

  long maxlen = 0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
      if (maxlen < lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime))
	maxlen = lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime);

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }

  //#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int luncollectedQ = nonZeroBlocks[index].first.first, rQ = nonZeroBlocks[index].first.second;
    int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];

    const std::vector<int>& colinds = rightOp.getActiveCols(rQ);
    for (int rrop=0; rrop <colinds.size(); rrop ++) {
      int rQPrime = colinds[rrop];
      
	const std::vector<int>& rowinds = c.getActiveRows(rQPrime);
	for (int l = 0; l < rowinds.size(); l++) {
	  int luncollectedQPrime = rowinds[l];
	  int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];
	  if (dotOp.allowed(dotQ, dotQPrime, LEFTOP.conjugacy()) && leftOp.allowed(lQ, lQPrime, LEFTOP.conjugacy())) {

	    StackMatrix m(dataArray[omprank], unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQ));
	    
	    double factor = scale*LEFTOP.get_scaling(unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]);	      
	    factor *= dmrginp.get_ninej()(unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					  LEFTOP.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					  unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	    factor *= Symmetry::spatial_ninej(unCollectedlketS->quanta[luncollectedQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
					      LEFTOP.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
					      unCollectedlbraS->quanta[luncollectedQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v[OMPRANK].get_symm().getirrep());

	    double scaleB = 1.0;
	    if (LEFTOP.conjugacy() == 'n') {
	      scaleB = dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), 
					   leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
					   lbraS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep());
	      scaleB *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , dotketS->quanta[dotQPrime].get_symm().getirrep(), unCollectedlketS->quanta[luncollectedQPrime].get_symm().getirrep(), 
					    leftOp.get_symm().getirrep(), dotOp.get_symm().getirrep(), LEFTOP.get_symm().getirrep(),
					    lbraS->quanta[lQ].get_symm().getirrep() , dotketS->quanta[dotQ].get_symm().getirrep(), unCollectedlbraS->quanta[luncollectedQ].get_symm().getirrep());
	      scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	      scaleB *= leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	      scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	      if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQPrime])) scaleB *= -1;
	    }
	    else {
	      scaleB = dmrginp.get_ninej()(lketS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQ].get_s().getirrep(), 
					   leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
					   lbraS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQPrime].get_s().getirrep());
	      scaleB *= Symmetry::spatial_ninej(lketS->quanta[lQ].get_symm().getirrep() , dotketS->quanta[dotQ].get_symm().getirrep(), unCollectedlketS->quanta[luncollectedQ].get_symm().getirrep(), 
					    leftOp.get_symm().getirrep(), dotOp.get_symm().getirrep(), LEFTOP.get_symm().getirrep(),
					    lbraS->quanta[lQPrime].get_symm().getirrep() , dotketS->quanta[dotQPrime].get_symm().getirrep(), unCollectedlbraS->quanta[luncollectedQPrime].get_symm().getirrep());
	      scaleB *= dotOp.operator_element(dotQPrime, dotQ)(1,1);
	      scaleB *= leftOp.get_scaling(lbraS->quanta[lQPrime], lketS->quanta[lQ]);
	      scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQPrime], dotketS->quanta[dotQ]);
	      if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQ])) scaleB *= -1;
	    }

	    int parity = rightOp.get_fermion() && IsFermion(unCollectedlketS->quanta[luncollectedQPrime]) ? -1 : 1;
	    factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);


	    MatrixMultiply (c.operator_element(luncollectedQPrime, rQPrime), 'n', rightOp.operator_element(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), 
			    m, 1.0, 0.);	      
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime, LEFTOP.conjugacy()), leftOp.conjugacy()=='n' ? LEFTOP.conjugacy() : TransposeOf(LEFTOP.conjugacy()), m, 'n',  v[OMPRANK].operator_element(luncollectedQ, rQ), factor*parity*scaleB);

	  }
	}
    }
  }
  
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
   Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }

  

}



void SpinAdapted::operatorfunctions::TensorMultiplysplitRight(const StackSparseMatrix& leftOp, const StackSparseMatrix& rightOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale)
{
  long starttime = globaltimer.totalwalltime();

    // can be used for situation with different bra and ket
  const int unCollectedrightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().unCollectedStateInfo->quanta.size ();
  const int unCollectedrightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().unCollectedStateInfo->quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int dotBraOpSz = cblock->get_rightBlock()->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int dotKetOpSz = cblock->get_rightBlock()->get_rightBlock()->get_ketStateInfo().quanta.size ();
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();

  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[0].get_nonZeroBlocks();

  long maxlen = 0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
      if (maxlen < lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime))
	maxlen = lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime);

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }

  //#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int lQ = nonZeroBlocks[index].first.first, runcollectedQ = nonZeroBlocks[index].first.second;
    int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];

    const std::vector<int>& colinds = leftOp.getActiveCols(lQ);
    for (int lrop=0; lrop <colinds.size(); lrop ++) {
      int lQPrime = colinds[lrop];
      
	const std::vector<int>& colinds2 = c.getActiveCols(lQPrime);
	for (int r = 0; r < colinds2.size(); r++) {
	  int runcollectedQPrime = colinds2[r];
	  int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

	  if (dotOp.allowed(dotQ, dotQPrime, RIGHTOP.conjugacy()) && rightOp.allowed(rQ, rQPrime, RIGHTOP.conjugacy())) {

	    StackMatrix m(dataArray[omprank], lketS->getquantastates(lQPrime), unCollectedrbraS->getquantastates(runcollectedQ));
	    
	    double factor = scale*leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);	      
	    factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					  leftOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), opQ.get_s().getirrep(),
					  lbraS->quanta[lQ].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	    factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_symm().getirrep() , c.get_deltaQuantum(0).get_symm().getirrep(), 
					  leftOp.get_symm().getirrep(), RIGHTOP.get_symm().getirrep(), opQ.get_symm().getirrep(),
					  lbraS->quanta[lQ].get_symm().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_symm().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_symm().getirrep());


	    double scaleB = 1.0;
	    if (RIGHTOP.conjugacy() == 'n') {
	      scaleB = dmrginp.get_ninej()(rketS->quanta[rQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep(), 
						   rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
						   rbraS->quanta[rQ].get_s().getirrep() , dotbraS->quanta[dotQ].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep());
	      scaleB *= Symmetry::spatial_ninej(rketS->quanta[rQPrime].get_symm().getirrep() , dotketS->quanta[dotQPrime].get_symm().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_symm().getirrep(), 
					     rightOp.get_symm().getirrep(), dotOp.get_symm().getirrep(), RIGHTOP.get_symm().getirrep(),
					     rbraS->quanta[rQ].get_symm().getirrep() , dotketS->quanta[dotQ].get_symm().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_symm().getirrep());
	      scaleB *= dotOp.operator_element(dotQ, dotQPrime, RIGHTOP.conjugacy())(1,1);
	      scaleB *= rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	      scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	      if (dotOp.get_fermion() && IsFermion(rketS->quanta[rQPrime])) scaleB *= -1;
	    }
	    else {
	      scaleB = dmrginp.get_ninej()(rketS->quanta[rQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQ].get_s().getirrep(), 
						   rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
						   rbraS->quanta[rQPrime].get_s().getirrep() , dotbraS->quanta[dotQPrime].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQPrime].get_s().getirrep());
	      scaleB *= Symmetry::spatial_ninej(rketS->quanta[rQ].get_symm().getirrep() , dotketS->quanta[dotQ].get_symm().getirrep(), unCollectedrketS->quanta[runcollectedQ].get_symm().getirrep(), 
					     rightOp.get_symm().getirrep(), dotOp.get_symm().getirrep(), RIGHTOP.get_symm().getirrep(),
					     rbraS->quanta[rQPrime].get_symm().getirrep() , dotketS->quanta[dotQPrime].get_symm().getirrep(), unCollectedrbraS->quanta[runcollectedQPrime].get_symm().getirrep());
	      scaleB *= dotOp.operator_element(dotQ, dotQPrime, RIGHTOP.conjugacy())(1,1);
	      scaleB *= rightOp.get_scaling(rbraS->quanta[rQPrime], rketS->quanta[rQ]);
	      scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQPrime], dotketS->quanta[dotQ]);
	      if (dotOp.get_fermion() && IsFermion(rketS->quanta[rQ])) scaleB *= -1;
	    }

	    int parity = RIGHTOP.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
	    factor *=  RIGHTOP.get_scaling(unCollectedrbraS->quanta[runcollectedQ], unCollectedrketS->quanta[runcollectedQPrime]);

	    MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', rightOp.operator_element(rQ, rQPrime, RIGHTOP.conjugacy()), 
			    rightOp.conjugacy() == 'n' ? TransposeOf(RIGHTOP.conjugacy()) : RIGHTOP.conjugacy(),  m, 1.0, 0.);	      
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftOp.conjugacy(), m, 'n',  v[OMPRANK].operator_element(lQ, runcollectedQ),
			    factor*parity*scaleB);

	  }
	}
    }
  }
  
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
   Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }

  

}


void SpinAdapted::operatorfunctions::MultiplyWithOwnTranspose(const StackSparseMatrix& a, StackSparseMatrix& c, Real scale)
{
  if (fabs(scale) < TINY) return;
  const int aSz = a.nrows();
  const int aSzPrime = a.ncols();
  
    assert (c.nrows() == a.nrows() &&
	    c.ncols() == a.nrows());
    
  int quanta_thrds = dmrginp.quanta_thrds();
#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
    for (int aQ = 0; aQ < aSz; ++aQ)
      for (int aQPrime = 0; aQPrime < aSzPrime; ++aQPrime)
	for (int bQPrime = 0; bQPrime < aSz; ++bQPrime)
	  {
	    if (a.allowed(aQ, aQPrime) && a.allowed(bQPrime, aQPrime) && c.allowed(aQ, bQPrime) ) {
	      
	      MatrixMultiply (a.operator_element(aQ, aQPrime), 'n', a.operator_element(bQPrime, aQPrime), 't',
			      c.operator_element(aQ, bQPrime), scale);
	    }
	  }
}

void SpinAdapted::operatorfunctions::Product (const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, StackSparseMatrix& c, double scale)
{
  const StateInfo* astate = &ablock->get_stateInfo(); 
  if (fabs(scale) < TINY) return;
  int rows = c.nrows();
  for (int cq = 0; cq < rows; ++cq)
    for (int cqprime = 0; cqprime < rows; ++cqprime)
      if (c.allowed(cq, cqprime))
	for (int aprime = 0; aprime < rows; aprime++)
	  if (a.allowed(cq, aprime) && b.allowed(aprime, cqprime))
	  {
	    int apj = astate->quanta[aprime].get_s().getirrep(), cqj = astate->quanta[cq].get_s().getirrep(), cqpj = astate->quanta[cqprime].get_s().getirrep();
	    double factor = a.get_scaling(astate->quanta[cq], astate->quanta[aprime]);
	    factor *= b.get_scaling(astate->quanta[aprime], astate->quanta[cqprime]);
	    if(dmrginp.spinAdapted()){
	      
	      factor *= racah(cqpj, b.get_spin().getirrep(), cqj, a.get_spin().getirrep(), apj, c.get_spin().getirrep()) * pow( (1.0*c.get_spin().getirrep()+1.0)*(1.0*apj+1.0), 0.5 )
		*pow(-1.0, static_cast<int>((b.get_spin().getirrep()+a.get_spin().getirrep()-c.get_spin().getirrep())/2.0));
	    }
	    MatrixMultiply(a.operator_element(cq, aprime), a.conjugacy(), b.operator_element(aprime, cqprime), b.conjugacy(),
			   c.operator_element(cq, cqprime), scale*factor, 1.0);
	    
	  }
}







void SpinAdapted::operatorfunctions::OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const StackSparseMatrix& op1, StackSparseMatrix& op2)
{
  const StateInfo& s = b.get_stateInfo();
  for (int lQ = 0; lQ< op2.nrows(); lQ++)
    for (int rQ = 0; rQ<op2.ncols(); rQ++)
      if (op2.allowed(lQ, rQ) && op1.allowed(lQ,rQ))
      {
	double factor = op1.get_scaling(s.quanta[lQ], s.quanta[rQ]);
	MatrixScaleAdd(scaleV*factor, op1.operator_element(lQ,rQ), op2.operator_element(lQ,rQ));
      }

}

/*
void SpinAdapted::operatorfunctions::MultiplyProduct(const StackSparseMatrix& a, const StackSparseMatrix& b, StackSparseMatrix& c, Real scale)
{
  if (fabs(scale) < TINY) return;
  const int aSz = a.nrows();
  const int aSzPrime = a.ncols();
  const int bSzPrime = b.ncols();

  assert (a.ncols() == b.nrows() && c.nrows() == a.nrows() &&
          c.ncols() == b.ncols());

  for (int aQ = 0; aQ < aSz; ++aQ)
    for (int aQPrime = 0; aQPrime < aSzPrime; ++aQPrime)
      for (int bQPrime = 0; bQPrime < bSzPrime; ++bQPrime)
        {
          if (a.allowed(aQ, aQPrime) && b.allowed(aQPrime, bQPrime) && c.allowed(aQ, bQPrime) ) {

	    MatrixMultiply (a.operator_element(aQ, aQPrime), a.conjugacy(), b.operator_element(aQPrime, bQPrime), b.conjugacy(),
			    c.operator_element(aQ, bQPrime), scale);
	  }
        }
}

*/

  

	      


