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


bool allowed(const std::vector<SpinQuantum>& dvec, const SpinQuantum& braQ, const SpinQuantum& ketQ) {
  for (int k =0; k<dvec.size(); k++) {
    if (braQ.allow(dvec[k], ketQ)) {
      return true;
    }
  }
  return false;
}

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
      for (int lQ = 0; lQ < leftBraOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQ, lQPrime))
	      {
		const StackMatrix& aop = a.operator_element(lQ, lQPrime);
		for (int rQ = 0; rQ < rightKetOpSz; ++rQ) 
		  if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
		    {
		      double fac=scale;
		      fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						 a.get_spin().getirrep(), 0, a.get_spin().getirrep(),
						 lbraS->quanta[lQ].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		      //fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
		      //a.get_symm().getirrep(), 0, a.get_symm().getirrep(),
		      //lbraS->quanta[lQ].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		      fac *= a.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
		      MatrixMultiply (aop, a.conjugacy(), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }
		
	      }
	  }
      }
    }
  else
    {
      for (int rQ = 0; rQ < rightBraOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightKetOpSz; ++rQPrime)
	  if (a.allowed(rQ, rQPrime))
	    {
	      const StackMatrix& aop = a.operator_element(rQ, rQPrime);
	      for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime) 
		if (v.allowed(lQPrime, rQ) && c.allowed(lQPrime, rQPrime)) {
		  double fac = scale;
		  fac *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					     0, a.get_spin().getirrep(), a.get_spin().getirrep(),
					     lketS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  //fac *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
		  //0, a.get_symm().getirrep(), a.get_symm().getirrep(),
		  //lketS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  fac *= a.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
		  double parity = a.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
		  
		  MatrixMultiply (c.operator_element(lQPrime, rQPrime), c.conjugacy(),
				  aop, TransposeOf(a.conjugacy()), v.operator_element(lQPrime, rQ), fac*parity);
		}
	      
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

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[omprank].get_nonZeroBlocks();

  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

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


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();


  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int luncollectedQPrime = nonZeroBlocks[index].first.first, rQPrime = nonZeroBlocks[index].first.second;
    int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

    const std::vector<int>& rowinds2 = rightOp.getActiveRows(rQPrime);
    for (int rrop=0; rrop <rowinds2.size(); rrop ++) {
      int rQ = rowinds2[rrop];

      StackMatrix m(dataArray[omprank], unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQ));
      ::Clear(m);
      
      const std::vector<int>& rowinds = v[OMPRANK].getActiveRows(rQ);
      for (int l = 0; l < rowinds.size(); l++) {
	int luncollectedQ = rowinds[l];
	int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
	if (dotOp.allowed(dotQ, dotQPrime, LEFTOP.conjugacy()) && leftOp.allowed(lQ, lQPrime, LEFTOP.conjugacy()) &&
	    allowed(LEFTOP.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime])) {


	  
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
	  if (fabs(factor*parity*scaleB) < TINY) continue; 

	  if (MatrixDotProduct(m,m) < TINY) {
	    if (rightOp.opName() == "OVERLAP")
	      copy(c.operator_element(luncollectedQPrime, rQPrime), m);
	    else
	      MatrixMultiply (c.operator_element(luncollectedQPrime, rQPrime), 'n', rightOp.operator_element(rQ, rQPrime), 
			      TransposeOf(rightOp.conjugacy()), m, 1.0, 0.);
	  }

	  if (leftOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(luncollectedQ, rQ));
	  else
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime, LEFTOP.conjugacy()), 
			    leftOp.conjugacy()=='n' ? LEFTOP.conjugacy() : TransposeOf(LEFTOP.conjugacy()), 
			    m, 'n',  v[OMPRANK].operator_element(luncollectedQ, rQ), factor*parity*scaleB);
	  
	}
      }
    }
  }
  
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
   Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }

  

}

void SpinAdapted::operatorfunctions::TensorMultiplyleftdot(const StackSparseMatrix& leftOp, StateInfo *cstate, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale)
{
  long starttime = globaltimer.totalwalltime();

    // can be used for situation with different bra and ket
  const int unCollectedleftBraOpSz = cstate->leftStateInfo->unCollectedStateInfo->quanta.size ();
  const int unCollectedleftKetOpSz = cstate->leftStateInfo->unCollectedStateInfo->quanta.size ();
  const int leftBraOpSz = cstate->leftStateInfo->leftStateInfo->quanta.size ();
  const int leftKetOpSz = cstate->leftStateInfo->leftStateInfo->quanta.size ();
  const int dotBraOpSz = cstate->leftStateInfo->rightStateInfo->quanta.size ();
  const int dotKetOpSz = cstate->leftStateInfo->rightStateInfo->quanta.size ();
  const int rightBraOpSz = cstate->rightStateInfo->quanta.size ();
  const int rightKetOpSz = cstate->rightStateInfo->quanta.size ();

  const boost::shared_ptr<StateInfo> unCollectedlbraS = cstate->leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cstate->leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cstate->leftStateInfo->leftStateInfo;
  const StateInfo* lketS = cstate->leftStateInfo->leftStateInfo;
  const StateInfo* dotbraS = cstate->leftStateInfo->rightStateInfo;
  const StateInfo* dotketS = cstate->leftStateInfo->rightStateInfo;
  const StateInfo* rbraS = cstate->rightStateInfo;
  const StateInfo* rketS = cstate->rightStateInfo;

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();

  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }
  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int luncollectedQPrime = nonZeroBlocks[index].first.first, rQPrime = nonZeroBlocks[index].first.second;
    int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

    int rQ = rQPrime;
      
    const std::vector<int>& rowinds = v[OMPRANK].getActiveRows(rQ);
    for (int l = 0; l < rowinds.size(); l++) {
      int luncollectedQ = rowinds[l];
      int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
      if (dotQ== dotQPrime && leftOp.allowed(lQ, lQPrime) &&
	    allowed(leftOp.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime])) {


	  
	double factor = scale*leftOp.get_scaling(unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]);	      
	factor *= dmrginp.get_ninej()(unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
				      leftOp.get_spin().getirrep(), 0, leftOp.get_spin().getirrep(),
				      unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  
	double scaleB = 1.0;

	scaleB = dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), 
				     leftOp.get_spin().getirrep(), 0, leftOp.get_spin().getirrep(),
				     lbraS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep());

	scaleB *= leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);

	MatrixScaleAdd(factor*scaleB*leftOp.operator_element(lQ, lQPrime)(1,1), c.operator_element(luncollectedQPrime, rQPrime), v[OMPRANK].operator_element(luncollectedQ, rQ));
	
	
      }
    }
  }
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
    Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }
  
  
  
}



void SpinAdapted::operatorfunctions::TensorMultiplydotop(const StackSparseMatrix& dotOp, StateInfo *cstate, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale)
{
  long starttime = globaltimer.totalwalltime();

    // can be used for situation with different bra and ket
  const int unCollectedleftBraOpSz = cstate->leftStateInfo->unCollectedStateInfo->quanta.size ();
  const int unCollectedleftKetOpSz = cstate->leftStateInfo->unCollectedStateInfo->quanta.size ();
  const int leftBraOpSz = cstate->leftStateInfo->leftStateInfo->quanta.size ();
  const int leftKetOpSz = cstate->leftStateInfo->leftStateInfo->quanta.size ();
  const int dotBraOpSz = cstate->leftStateInfo->rightStateInfo->quanta.size ();
  const int dotKetOpSz = cstate->leftStateInfo->rightStateInfo->quanta.size ();
  const int rightBraOpSz = cstate->rightStateInfo->quanta.size ();
  const int rightKetOpSz = cstate->rightStateInfo->quanta.size ();

  const boost::shared_ptr<StateInfo> unCollectedlbraS = cstate->leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cstate->leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cstate->leftStateInfo->leftStateInfo;
  const StateInfo* lketS = cstate->leftStateInfo->leftStateInfo;
  const StateInfo* dotbraS = cstate->leftStateInfo->rightStateInfo;
  const StateInfo* dotketS = cstate->leftStateInfo->rightStateInfo;
  const StateInfo* rbraS = cstate->rightStateInfo;
  const StateInfo* rketS = cstate->rightStateInfo;


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();


  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;
  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }

  SpinQuantum hq(0,SpinSpace(0),IrrepSpace(0));

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int luncollectedQPrime = nonZeroBlocks[index].first.first, rQPrime = nonZeroBlocks[index].first.second;
    int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

    int rQ = rQPrime;

    StackMatrix m(dataArray[omprank], unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQ));
    ::Clear(m);
      
    const std::vector<int>& rowinds = v[OMPRANK].getActiveRows(rQ);
    for (int l = 0; l < rowinds.size(); l++) {
      int luncollectedQ = rowinds[l];
      int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
      if (dotOp.allowed(dotQ, dotQPrime) && lQ == lQPrime &&
	  allowed(dotOp.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime])) {
	
	double factor = scale*dotOp.get_scaling(unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]);	      
	factor *= dmrginp.get_ninej()(unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
				      dotOp.get_spin().getirrep(), 0, dotOp.get_spin().getirrep(),
				      unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	
	double scaleB = 1.0;
	
	scaleB = dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), 
				     0, dotOp.get_spin().getirrep(), dotOp.get_spin().getirrep(),
				     lbraS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep());
	
	scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQPrime])) scaleB *= -1;
	
	MatrixScaleAdd(factor*scaleB, c.operator_element(luncollectedQPrime, rQPrime), v[OMPRANK].operator_element(luncollectedQ, rQ));
	
      }
    }
  }

  
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
    Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }
  
  

}


void SpinAdapted::operatorfunctions::TensorMultiplysplitLeftElement(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int rQPrime, double scale)
{
  long starttime = globaltimer.totalwalltime();

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();
  const std::vector<int>& rowindsc = c.getActiveRows(rQPrime);
  if (rowindsc.size() == 0) return;
  const std::vector<int>& rowindsr = rightOp.getActiveRows(rQPrime);
  if (rowindsr.size() == 0) return;

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


  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray;
  dataArray = Stackmem[OMPRANK].allocate(maxlen);


  for (int rrop = 0; rrop<rowindsr.size(); rrop++) {
    int rQ = rowindsr[rrop];
    bool deallocate = rightOp.memoryUsed() == 0 ? true : false;
    StackMatrix ropm;
    //make the rightop element
    double* ropdata = deallocate ? Stackmem[omprank].allocate(rbraS->quantaStates[rQ]* rketS->quantaStates[rQPrime]) :  0;
    if (rightOp.conjugacy() == 'n') ropm = deallocate ? StackMatrix(ropdata, rbraS->quantaStates[rQ], rketS->quantaStates[rQPrime]) : StackMatrix();
    else ropm = deallocate ? StackMatrix(ropdata, rketS->quantaStates[rQPrime], rbraS->quantaStates[rQ]) : StackMatrix();
    
    memset(ropdata, 0, ropm.Storage() * sizeof(double));
    const_cast<StackSparseMatrix&>(rightOp).build(ropm, rQ, rQPrime, *cblock->get_rightBlock());


    for (int rc=0; rc <rowindsc.size(); rc ++) {
      int luncollectedQPrime = rowindsc[rc];

      int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

      StackMatrix m(dataArray, unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQ));
      ::Clear(m);
      
      const std::vector<int>& rowinds = v[OMPRANK].getActiveRows(rQ);
      for (int l = 0; l < rowinds.size(); l++) {
      int luncollectedQ = rowinds[l];
      int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
      if (dotOp.allowed(dotQ, dotQPrime, LEFTOP.conjugacy()) && leftOp.allowed(lQ, lQPrime, LEFTOP.conjugacy()) &&
	  allowed(LEFTOP.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime])) {
	
	
	
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
	if (fabs(factor*parity*scaleB) < TINY) continue; 
	
	if (MatrixDotProduct(m,m) < TINY) {
	  if (rightOp.opName() == "OVERLAP")
	    copy(c.operator_element(luncollectedQPrime, rQPrime), m);
	  else {

	    
	    MatrixMultiply (c.operator_element(luncollectedQPrime, rQPrime), 'n', ropm, 
			    TransposeOf(rightOp.conjugacy()), m, 1.0, 0.);
	  }
	}
	
	if (leftOp.opName() == "OVERLAP")
	  MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(luncollectedQ, rQ));
	else
	  MatrixMultiply (leftOp.operator_element(lQ, lQPrime, LEFTOP.conjugacy()), 
			  leftOp.conjugacy()=='n' ? LEFTOP.conjugacy() : TransposeOf(LEFTOP.conjugacy()), 
			  m, 'n',  v[OMPRANK].operator_element(luncollectedQ, rQ), factor*parity*scaleB);
	
      }
      }
      
    }
  
    if (deallocate) Stackmem[omprank].deallocate(ropm.Store(), ropm.Storage());
  }
  Stackmem[OMPRANK].deallocate(dataArray, maxlen);
}

void SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitLeftElement(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int rQPrime, double scale)
{
  long starttime = globaltimer.totalwalltime();

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();
  const std::vector<int>& rowindsc = c.getActiveRows(rQPrime);
  if (rowindsc.size() == 0) return;
  //const std::vector<int>& rowindsr = rightOp.getActiveRows(rQPrime);
  //if (rowindsr.size() == 0) return;

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



  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;
  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray1, *dataArray2;
  dataArray1 = Stackmem[OMPRANK].allocate(maxlen);


  for (int rQ = 0; rQ< rightBraOpSz; rQ++)
  if (allowed(rightOp.get_deltaQuantum(), rbraS->quanta[rQ], rketS->quanta[rQPrime])) {

    bool deallocate = rightOp.memoryUsed() == 0 ? true : false;
    StackMatrix ropm;
    //make the rightop element
    double* ropdata = deallocate ? Stackmem[omprank].allocate(rbraS->quantaStates[rQ]* rketS->quantaStates[rQPrime]) :  0;
    ropm = deallocate ? StackMatrix(ropdata, rbraS->quantaStates[rQ], rketS->quantaStates[rQPrime]) : StackMatrix();    
    memset(ropdata, 0, ropm.Storage() * sizeof(double));
    const_cast<StackSparseMatrix&>(rightOp).build(ropm, rQ, rQPrime, *cblock->get_rightBlock());


    for (int rc=0; rc <rowindsc.size(); rc ++) {
      int luncollectedQPrime = rowindsc[rc];

      int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

      StackMatrix m(dataArray1, unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQ));
      ::Clear(m);
      
      const std::vector<int>& rowinds = v[OMPRANK].getActiveRows(rQ);
      for (int l = 0; l < rowinds.size(); l++) {
	int luncollectedQ = rowinds[l];
	int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
	
	//NONTRANSPOSE
	if (dotOp.allowed(dotQ, dotQPrime, LEFTOP.conjugacy()) && leftOp.allowed(lQ, lQPrime, LEFTOP.conjugacy()) &&
	    allowed(LEFTOP.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime])) {
	
	  
	  
	  double factor = scale*LEFTOP.get_scaling(unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]);	      
	  factor *= dmrginp.get_ninej()(unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					LEFTOP.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  double scaleB = 1.0;
	  
	  if (LEFTOP.conjugacy() == 'n') {
	    scaleB = dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), 
					 leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
					 lbraS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep());
	    scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	    scaleB *= leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	    scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	    if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQPrime])) scaleB *= -1;
	  }
	  else {
	    scaleB = dmrginp.get_ninej()(lketS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQ].get_s().getirrep(), 
					 leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
					 lbraS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQPrime].get_s().getirrep());
	    scaleB *= dotOp.operator_element(dotQPrime, dotQ)(1,1);
	    scaleB *= leftOp.get_scaling(lbraS->quanta[lQPrime], lketS->quanta[lQ]);
	    scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQPrime], dotketS->quanta[dotQ]);
	    if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQ])) scaleB *= -1;
	  }
	  
	  int parity = rightOp.get_fermion() && IsFermion(unCollectedlketS->quanta[luncollectedQPrime]) ? -1 : 1;
	  factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	  if (fabs(factor*parity*scaleB) < TINY) continue; 
	  
	  if (MatrixDotProduct(m,m) < TINY) {
	    if (rightOp.opName() == "OVERLAP")
	      copy(c.operator_element(luncollectedQPrime, rQPrime), m);
	    else {	    
	      MatrixMultiply (c.operator_element(luncollectedQPrime, rQPrime), 'n', ropm, 
			      TransposeOf(rightOp.conjugacy()), m, 1.0, 0.);
	    }
	  }
	  
	  if (leftOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(luncollectedQ, rQ));
	  else
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime, LEFTOP.conjugacy()), 
			    leftOp.conjugacy()=='n' ? LEFTOP.conjugacy() : TransposeOf(LEFTOP.conjugacy()), 
			    m, 'n',  v[OMPRANK].operator_element(luncollectedQ, rQ), factor*parity*scaleB);
	}
      }
    }	

    const std::vector<int>& rowindsc2 = c.getActiveRows(rQ);
    if (rowindsc2.size() == 0) {
      if (deallocate) Stackmem[omprank].deallocate(ropdata, rbraS->quantaStates[rQ]* rketS->quantaStates[rQPrime]);
      continue;
    }

    for (int rc=0; rc <rowindsc2.size(); rc ++) {
      int luncollectedQPrime = rowindsc2[rc];

      int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

      StackMatrix m(dataArray1, unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQPrime));
      ::Clear(m);
	
      const std::vector<int>& rowinds = v[OMPRANK].getActiveRows(rQPrime);
      for (int l = 0; l < rowinds.size(); l++) {
	int luncollectedQ = rowinds[l];
	int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
	//TRANSPOSE
	if (dotOp.allowed(dotQ, dotQPrime, TransposeOf(LEFTOP.conjugacy())) && leftOp.allowed(lQ, lQPrime, TransposeOf(LEFTOP.conjugacy())) &&
	    allowed(LEFTOP.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQPrime], unCollectedlketS->quanta[luncollectedQ])) {

	  double factor = scale;
	  if (LEFTOP.conjugacy() == 'n') factor*=getStandAlonescaling(-LEFTOP.get_deltaQuantum(0), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]);	      
	  factor *= dmrginp.get_ninej()(unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep(),c.get_deltaQuantum(0).get_s().getirrep(), 
					LEFTOP.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep(), rbraS->quanta[rQPrime].get_s().getirrep(),v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  
	  double scaleB = 1.0;
	  
	  if (LEFTOP.conjugacy() == 't') {
	    scaleB = dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), 
					 leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
					 lbraS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep());
	    scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	    scaleB *= leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	    scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	    if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQPrime])) scaleB *= -1;
	  }
	  else {
	    scaleB = dmrginp.get_ninej()(lketS->quanta[lQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedlketS->quanta[luncollectedQ].get_s().getirrep(), 
					 leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
					 lbraS->quanta[lQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedlbraS->quanta[luncollectedQPrime].get_s().getirrep());
	    scaleB *= dotOp.operator_element(dotQPrime, dotQ)(1,1);
	    scaleB *= leftOp.get_scaling(lbraS->quanta[lQPrime], lketS->quanta[lQ]);
	    scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQPrime], dotketS->quanta[dotQ]);
	    if (dotOp.get_fermion() && IsFermion(lketS->quanta[lQ])) scaleB *= -1;
	  }
	  
	  int parity = rightOp.get_fermion() && IsFermion(unCollectedlketS->quanta[luncollectedQPrime]) ? -1 : 1;
	  if (!dmrginp.spinAdapted() && LEFTOP.get_fermion()) parity *= -1;

	  factor *=  getStandAlonescaling(-rightOp.get_deltaQuantum(0), rbraS->quanta[rQPrime], rketS->quanta[rQ]);
	  if (fabs(factor*parity*scaleB) < TINY) continue; 
	  
	  if (rightOp.opName() == "OVERLAP")
	    copy(c.operator_element(luncollectedQPrime, rQ), m);
	  else {	    
	    MatrixMultiply (c.operator_element(luncollectedQPrime, rQ), 'n', ropm, 
			    'n', m, 1.0, 0.);
	  }
	  
	  
	  if (leftOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(luncollectedQ, rQPrime));
	  else
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime, TransposeOf(LEFTOP.conjugacy())), 
			    leftOp.conjugacy()=='n' ? TransposeOf(LEFTOP.conjugacy()) : LEFTOP.conjugacy(), 
			    m, 'n',  v[OMPRANK].operator_element(luncollectedQ, rQPrime), factor*parity*scaleB);
	  
	}
      }
    }  

    if (deallocate) Stackmem[omprank].deallocate(ropdata, rbraS->quantaStates[rQ]* rketS->quantaStates[rQPrime]);
    
  }
  Stackmem[OMPRANK].deallocate(dataArray1, maxlen);

}


void SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElement(const StackSparseMatrix& leftOp, const StackSparseMatrix& rightOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int lQPrime, double scale, bool doTranspose)
{
  long starttime = globaltimer.totalwalltime();

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();
  const std::vector<int>& colindsc = c.getActiveCols(lQPrime);
  if (colindsc.size() == 0) return;

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

  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray1;
  dataArray1 = Stackmem[OMPRANK].allocate(maxlen);

  for (int lQ = 0; lQ< leftBraOpSz; lQ++)
  if (allowed(leftOp.get_deltaQuantum(), lbraS->quanta[lQ], lketS->quanta[lQPrime])) {

    bool deallocate = leftOp.memoryUsed() == 0 ? true : false;
    StackMatrix lopm;
    //make the rightop element
    double* lopdata = deallocate ? Stackmem[omprank].allocate(lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]) :  0;
    if (leftOp.conjugacy() =='n') 
      lopm = deallocate ? StackMatrix(lopdata, lbraS->quantaStates[lQ], lketS->quantaStates[lQPrime]) : StackMatrix();    
    else
      lopm = deallocate ? StackMatrix(lopdata, lbraS->quantaStates[lQPrime], lketS->quantaStates[lQ]) : StackMatrix();    
    memset(lopdata, 0, lopm.Storage() * sizeof(double));
    const_cast<StackSparseMatrix&>(leftOp).build(lopm, lQ, lQPrime, *cblock->get_leftBlock());

    for (int rc=0; rc <colindsc.size(); rc ++) {
      int runcollectedQPrime = colindsc[rc];

      int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];
      StackMatrix m(dataArray1, lbraS->getquantastates(lQ), unCollectedrketS->getquantastates(runcollectedQPrime));
      ::Clear(m);
      
      const std::vector<int>& colinds2 = v[OMPRANK].getActiveCols(lQ);
      for (int r = 0; r < colinds2.size(); r++) {
	int runcollectedQ = colinds2[r];
	int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];

	//NONTRANSPOSE
	if (dotOp.allowed(dotQ, dotQPrime, RIGHTOP.conjugacy()) && rightOp.allowed(rQ, rQPrime, RIGHTOP.conjugacy()) &&
	    allowed(RIGHTOP.get_deltaQuantum(), unCollectedrbraS->quanta[runcollectedQ], unCollectedrketS->quanta[runcollectedQPrime])) {


	  double factor = scale*leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);	      
	  factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					leftOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), opQ.get_s().getirrep(),
					lbraS->quanta[lQ].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  
	  
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
	  
	  int parity = 1;
	  if (leftOp.conjugacy() == 'n') parity*=RIGHTOP.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
	  else parity*=RIGHTOP.get_fermion() && IsFermion(lketS->quanta[lQ]) ? -1 : 1;
	  factor *=  RIGHTOP.get_scaling(unCollectedrbraS->quanta[runcollectedQ], unCollectedrketS->quanta[runcollectedQPrime]);
	  if (fabs(factor*parity*scaleB) < TINY) continue; 

	  if (MatrixDotProduct(m,m) < TINY) {
	    if (leftOp.opName() == "OVERLAP")
	      copy(c.operator_element(lQPrime, runcollectedQPrime), m);
	    else
	      MatrixMultiply (lopm, leftOp.conjugacy(), c.operator_element(lQPrime, runcollectedQPrime),
			      'n', m, 1.0, 0.0);
	  }

	  if (rightOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(lQ, runcollectedQ));
	  else
	    MatrixMultiply (m, 'n', rightOp.operator_element(rQ, rQPrime, RIGHTOP.conjugacy()), 
			    rightOp.conjugacy() == 'n' ? TransposeOf(RIGHTOP.conjugacy()) : RIGHTOP.conjugacy(), 
			    v[OMPRANK].operator_element(lQ, runcollectedQ), factor*parity*scaleB);	      

	}
      }
    }

    if (!doTranspose) {
      if (deallocate) Stackmem[omprank].deallocate(lopdata, lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]);
      continue;
    }
    const std::vector<int>& colindsc2 = c.getActiveCols(lQ);

    if (colindsc2.size() == 0) {
      if (deallocate) Stackmem[omprank].deallocate(lopdata, lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]);
      continue;
    }

    for (int rc=0; rc <colindsc2.size(); rc ++) {
      int runcollectedQPrime = colindsc2[rc];

      int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];
      StackMatrix m(dataArray1, lbraS->getquantastates(lQPrime), unCollectedrketS->getquantastates(runcollectedQPrime));
      ::Clear(m);
      
      const std::vector<int>& colinds2 = v[OMPRANK].getActiveCols(lQPrime);
      for (int r = 0; r < colinds2.size(); r++) {
	int runcollectedQ = colinds2[r];
	int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];

	//TRANSPOSE
	if (dotOp.allowed(dotQ, dotQPrime, TransposeOf(RIGHTOP.conjugacy())) && rightOp.allowed(rQ, rQPrime, TransposeOf(RIGHTOP.conjugacy())) &&
	    allowed(RIGHTOP.get_deltaQuantum(), unCollectedrbraS->quanta[runcollectedQPrime], unCollectedrketS->quanta[runcollectedQ])) {


	  //const_cast<StackSparseMatrix&>(leftOp).set_conjugacy('t');
	  double factor = scale*getStandAlonescaling(-leftOp.get_deltaQuantum(0), lbraS->quanta[lQPrime], lketS->quanta[lQ]);
	  //const_cast<StackSparseMatrix&>(leftOp).set_conjugacy('n');
	  factor *= dmrginp.get_ninej()(lketS->quanta[lQ].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					leftOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), opQ.get_s().getirrep(),
					lbraS->quanta[lQPrime].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  
	  
	  double scaleB = 1.0;
	  
	  if (RIGHTOP.conjugacy() == 't') {
	    scaleB = dmrginp.get_ninej()(rketS->quanta[rQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep(), 
					 rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
					 rbraS->quanta[rQ].get_s().getirrep() , dotbraS->quanta[dotQ].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep());
	    
	    scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	    scaleB *= rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	    scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	    if (dotOp.get_fermion() && IsFermion(rketS->quanta[rQPrime])) scaleB *= -1;
	  }
	  else {
	    scaleB = dmrginp.get_ninej()(rketS->quanta[rQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQ].get_s().getirrep(), 
					 rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
					 rbraS->quanta[rQPrime].get_s().getirrep() , dotbraS->quanta[dotQPrime].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQPrime].get_s().getirrep());
	    scaleB *= dotOp.operator_element(dotQPrime, dotQ)(1,1);
	    scaleB *= rightOp.get_scaling(rbraS->quanta[rQPrime], rketS->quanta[rQ]);
	    scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQPrime], dotketS->quanta[dotQ]);
	    if (dotOp.get_fermion() && IsFermion(rketS->quanta[rQ])) scaleB *= -1;
	  }
	  
	  int parity = RIGHTOP.get_fermion() && IsFermion(lketS->quanta[lQ]) ? -1 : 1; //********
    if (!dmrginp.spinAdapted() && RIGHTOP.get_fermion()) parity *= -1;
	  char lc = RIGHTOP.conjugacy();
	  //const_cast<StackSparseMatrix&>(RIGHTOP).set_conjugacy(TransposeOf(lc));
	  if (lc == 'n') factor *=  getStandAlonescaling(-RIGHTOP.get_deltaQuantum(0), unCollectedrbraS->quanta[runcollectedQ], unCollectedrketS->quanta[runcollectedQPrime]);
	  //const_cast<StackSparseMatrix&>(RIGHTOP).set_conjugacy(lc);

	  if (fabs(factor*parity*scaleB) < TINY) continue; 

	  if (MatrixDotProduct(m,m) < TINY) {
	    if (leftOp.opName() == "OVERLAP")
	      copy(c.operator_element(lQ, runcollectedQPrime), m);
	    else
	      MatrixMultiply (lopm, TransposeOf(leftOp.conjugacy()), c.operator_element(lQ, runcollectedQPrime),
			      'n', m, 1.0, 0.0);
	  }

	  if (rightOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(lQPrime, runcollectedQ));
	  else
	    MatrixMultiply (m, 'n', rightOp.operator_element(rQ, rQPrime, TransposeOf(RIGHTOP.conjugacy())), 
			    rightOp.conjugacy() == 'n' ? RIGHTOP.conjugacy() : TransposeOf(RIGHTOP.conjugacy()), 
			    v[OMPRANK].operator_element(lQPrime, runcollectedQ), factor*parity*scaleB);	      

	}
      }

    }

    if (deallocate) Stackmem[omprank].deallocate(lopdata, lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]);

  }    
  Stackmem[OMPRANK].deallocate(dataArray1, maxlen);
}


void SpinAdapted::operatorfunctions::TensorMultiplyCDxCDsplitRightElementcopy(const StackSparseMatrix& leftOp, const StackSparseMatrix& rightOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, int lQPrime, double scale, bool doTranspose)
{
  long starttime = globaltimer.totalwalltime();

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();
  const std::vector<int>& colindsc = c.getActiveCols(lQPrime);
  if (colindsc.size() == 0) return;

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


  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray1;
  dataArray1 = Stackmem[OMPRANK].allocate(maxlen);


  for (int lQ = 0; lQ< leftBraOpSz; lQ++)
  if (allowed(leftOp.get_deltaQuantum(), lbraS->quanta[lQ], lketS->quanta[lQPrime])) {

    bool deallocate = leftOp.memoryUsed() == 0 ? true : false;
    StackMatrix lopm;
    //make the rightop element
    double* lopdata = deallocate ? Stackmem[omprank].allocate(lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]) :  0;
    if (leftOp.conjugacy() =='n') 
      lopm = deallocate ? StackMatrix(lopdata, lbraS->quantaStates[lQ], lketS->quantaStates[lQPrime]) : StackMatrix();    
    else
      lopm = deallocate ? StackMatrix(lopdata, lbraS->quantaStates[lQPrime], lketS->quantaStates[lQ]) : StackMatrix();    
    memset(lopdata, 0, lopm.Storage() * sizeof(double));
    const_cast<StackSparseMatrix&>(leftOp).build(lopm, lQ, lQPrime, *cblock->get_leftBlock());

    if (!doTranspose) {
      if (deallocate) Stackmem[omprank].deallocate(lopdata, lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]);
      continue;
    }
    const std::vector<int>& colindsc2 = c.getActiveCols(lQ);

    if (colindsc2.size() == 0) {
      if (deallocate) Stackmem[omprank].deallocate(lopdata, lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]);
      continue;
    }

    for (int rc=0; rc <colindsc2.size(); rc ++) {
      int runcollectedQPrime = colindsc2[rc];

      int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];
      StackMatrix m(dataArray1, lbraS->getquantastates(lQPrime), unCollectedrketS->getquantastates(runcollectedQPrime));
      ::Clear(m);
      
      const std::vector<int>& colinds2 = v[OMPRANK].getActiveCols(lQPrime);
      for (int r = 0; r < colinds2.size(); r++) {
	int runcollectedQ = colinds2[r];
	int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];

	//TRANSPOSE
	if (dotOp.allowed(dotQ, dotQPrime, TransposeOf(RIGHTOP.conjugacy())) && rightOp.allowed(rQ, rQPrime, TransposeOf(RIGHTOP.conjugacy())) &&
	    allowed(RIGHTOP.get_deltaQuantum(), unCollectedrbraS->quanta[runcollectedQPrime], unCollectedrketS->quanta[runcollectedQ])) {


	  //const_cast<StackSparseMatrix&>(leftOp).set_conjugacy('t');
	  double factor = scale*getStandAlonescaling(-leftOp.get_deltaQuantum(0), lbraS->quanta[lQPrime], lketS->quanta[lQ]);
	  //const_cast<StackSparseMatrix&>(leftOp).set_conjugacy('n');
	  factor *= dmrginp.get_ninej()(lketS->quanta[lQ].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					leftOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), opQ.get_s().getirrep(),
					lbraS->quanta[lQPrime].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  
	  
	  double scaleB = 1.0;
	  
	  if (RIGHTOP.conjugacy() == 't') {
	      scaleB = dmrginp.get_ninej()(rketS->quanta[rQPrime].get_s().getirrep() , dotketS->quanta[dotQPrime].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep(), 
						   rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
						   rbraS->quanta[rQ].get_s().getirrep() , dotbraS->quanta[dotQ].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep());

	      scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	      scaleB *= rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	      scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQ], dotketS->quanta[dotQPrime]);
	      if (dotOp.get_fermion() && IsFermion(rketS->quanta[rQPrime])) scaleB *= -1;
	    }
	  else {
	      scaleB = dmrginp.get_ninej()(rketS->quanta[rQ].get_s().getirrep() , dotketS->quanta[dotQ].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQ].get_s().getirrep(), 
						   rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
						   rbraS->quanta[rQPrime].get_s().getirrep() , dotbraS->quanta[dotQPrime].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQPrime].get_s().getirrep());
	      scaleB *= dotOp.operator_element(dotQ, dotQPrime)(1,1);
	      scaleB *= rightOp.get_scaling(rbraS->quanta[rQPrime], rketS->quanta[rQ]);
	      scaleB *= dotOp.get_scaling(dotbraS->quanta[dotQPrime], dotketS->quanta[dotQ]);
	      if (dotOp.get_fermion() && IsFermion(rketS->quanta[rQ])) scaleB *= -1;
	    }
	  
	  int parity = RIGHTOP.get_fermion() && IsFermion(lketS->quanta[lQ]) ? -1 : 1; //********
	  char lc = RIGHTOP.conjugacy();
	  //const_cast<StackSparseMatrix&>(RIGHTOP).set_conjugacy(TransposeOf(lc));
	  if (lc == 'n') factor *= getStandAlonescaling(-RIGHTOP.get_deltaQuantum(0), unCollectedrbraS->quanta[runcollectedQ], unCollectedrketS->quanta[runcollectedQPrime]);
	  //const_cast<StackSparseMatrix&>(RIGHTOP).set_conjugacy(lc);

	  if (fabs(factor*parity*scaleB) < TINY) continue; 

	  if (MatrixDotProduct(m,m) < TINY) {
	    if (leftOp.opName() == "OVERLAP")
	      copy(c.operator_element(lQ, runcollectedQPrime), m);
	    else
	      MatrixMultiply (lopm, TransposeOf(leftOp.conjugacy()), c.operator_element(lQ, runcollectedQPrime),
			      'n', m, 1.0, 0.0);
	  }

	  if (rightOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(lQPrime, runcollectedQ));
	  else
	    MatrixMultiply (m, 'n', rightOp.operator_element(rQ, rQPrime, TransposeOf(RIGHTOP.conjugacy())), 
			    rightOp.conjugacy() == 'n' ? RIGHTOP.conjugacy() : TransposeOf(RIGHTOP.conjugacy()), 
			    v[OMPRANK].operator_element(lQPrime, runcollectedQ), factor*parity*scaleB);	      

	}
	if (lQPrime == 19 && runcollectedQ == 1) {
	  pout << lQPrime<<"  "<<runcollectedQ<<endl;
	  pout <<"l "<< MatrixDotProduct(lopm, lopm)<<endl;
	  pout <<"r "<< MatrixDotProduct(rightOp.operator_element(rQ, rQPrime, TransposeOf(RIGHTOP.conjugacy())), rightOp.operator_element(rQ, rQPrime, TransposeOf(RIGHTOP.conjugacy())))<<endl;
	  pout <<"c "<< MatrixDotProduct(c.operator_element(lQ, runcollectedQPrime), c.operator_element(lQ, runcollectedQPrime))<<endl;
	  pout <<"v "<< MatrixDotProduct(v[omprank].operator_element(lQPrime, runcollectedQ), v[omprank].operator_element(lQPrime, runcollectedQ))<<endl;
	}
      }

    }

    if (deallocate) Stackmem[omprank].deallocate(lopdata, lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]);

  }    
  Stackmem[OMPRANK].deallocate(dataArray1, maxlen);
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


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();

  long maxlen = 0, maxrow=0, maxcol=0;
  for (int lQ=0; lQ <leftBraOpSz; lQ++)
    if (maxrow <lbraS->getquantastates(lQ)) maxrow = lbraS->getquantastates(lQ);
  for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime++)
    if (maxcol <rketS->getquantastates(rQPrime)) maxcol = rketS->getquantastates(rQPrime);

  maxlen = maxrow*maxcol;

  int OMPRANK = omprank;

  int quanta_thrds = dmrginp.quanta_thrds();

  double* dataArray[quanta_thrds];
  for (int q = 0; q < quanta_thrds; q++) {
    dataArray[q] = Stackmem[OMPRANK].allocate(maxlen);
  }

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int lQPrime = nonZeroBlocks[index].first.first, runcollectedQPrime = nonZeroBlocks[index].first.second;
    int rQPrime = unCollectedrbraS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrbraS->rightUnMapQuanta[runcollectedQPrime];

    const std::vector<int>& rowinds = leftOp.getActiveRows(lQPrime);
    for (int lrop=0; lrop <rowinds.size(); lrop ++) {
      int lQ = rowinds[lrop];

      StackMatrix m(dataArray[omprank], lketS->getquantastates(lQ), unCollectedrbraS->getquantastates(runcollectedQPrime));
      ::Clear(m);
      
      const std::vector<int>& colinds2 = v[OMPRANK].getActiveCols(lQ);
      for (int r = 0; r < colinds2.size(); r++) {
	int runcollectedQ = colinds2[r];
	int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];

	if (dotOp.allowed(dotQ, dotQPrime, RIGHTOP.conjugacy()) && rightOp.allowed(rQ, rQPrime, RIGHTOP.conjugacy())) {
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
	  if (fabs(factor*parity*scaleB) < TINY) continue; 

	  if (MatrixDotProduct(m,m) < TINY) {
	    if (leftOp.opName() == "OVERLAP")
	      copy(c.operator_element(lQPrime, runcollectedQPrime), m);
	    else
	      MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftOp.conjugacy(), c.operator_element(lQPrime, runcollectedQPrime),
			      'n', m, 1.0, 0.0);
	  }

	  if (rightOp.opName() == "OVERLAP")
	    MatrixScaleAdd(factor*parity*scaleB, m, v[OMPRANK].operator_element(lQ, runcollectedQ));
	  else
	    MatrixMultiply (m, 'n', rightOp.operator_element(rQ, rQPrime, RIGHTOP.conjugacy()), 
			    rightOp.conjugacy() == 'n' ? TransposeOf(RIGHTOP.conjugacy()) : RIGHTOP.conjugacy(), 
			    v[OMPRANK].operator_element(lQ, runcollectedQ), factor*parity*scaleB);	      

	}
      }
    }
  }
  
  
  for (int q = quanta_thrds-1; q > -1 ; q--) {
   Stackmem[OMPRANK].deallocate(dataArray[q], maxlen);
  }

  

}



//*****************************************************************
void SpinAdapted::operatorfunctions::TensorMultiplyLeftLeft(const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale)
{
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


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();

  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int luncollectedQPrime = nonZeroBlocks[index].first.first, rQPrime = nonZeroBlocks[index].first.second;
    int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime];

    const std::vector<int>& rowinds = v.getActiveRows(rQPrime);
    for (int rrq = 0; rrq < rowinds.size(); ++rrq)  {
      int luncollectedQ = rowinds[rrq];
      int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ];
      if (a.allowed(lQ, lQPrime))
	MatrixMultiply (a.operator_element(lQ, lQPrime), a.conjugacy(), c.operator_element(luncollectedQPrime, rQPrime), c.conjugacy(),
			v.operator_element(luncollectedQ, rQPrime), 1.0, 0.0);
      
    }
  }

}


void SpinAdapted::operatorfunctions::TensorMultiplyRight(const StackSparseMatrix& rightOp, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& LEFTOP, const StackSpinBlock *cblock, StackWavefunction& w, StackWavefunction* v, const SpinQuantum opQ, double scale)
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


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[omprank].get_nonZeroBlocks();

  int OMPRANK = omprank;
  int quanta_thrds = dmrginp.quanta_thrds();

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int luncollectedQ = nonZeroBlocks[index].first.first, rQ = nonZeroBlocks[index].first.second;
    int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];

    const std::vector<int>& colinds = rightOp.getActiveCols(rQ);
    for (int rrop=0; rrop <colinds.size(); rrop ++) {
      int rQPrime = colinds[rrop];
      
      const std::vector<int>& rowinds = w.getActiveRows(rQPrime);
      for (int l = 0; l < rowinds.size(); l++) {
	int luncollectedQPrime = rowinds[l];
	int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];
	if (dotOp.allowed(dotQ, dotQPrime, LEFTOP.conjugacy()) && leftOp.allowed(lQ, lQPrime, LEFTOP.conjugacy())) {

	  
	  double factor = scale*LEFTOP.get_scaling(unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]);	      
	  factor *= dmrginp.get_ninej()(unCollectedlketS->quanta[luncollectedQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep(), 
					LEFTOP.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					unCollectedlbraS->quanta[luncollectedQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
	  factor *= Symmetry::spatial_ninej(unCollectedlketS->quanta[luncollectedQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), v[OMPRANK].get_symm().getirrep(), 
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
	  if (fabs(factor*parity*scaleB) < TINY) continue; 
	  
	  MatrixMultiply (w.operator_element(luncollectedQPrime, rQPrime), 'n', rightOp.operator_element(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), 
			  v[OMPRANK].operator_element(luncollectedQ, rQ), factor*parity*scaleB);	      
	  
	}
      }
    }
  }
  

}



void SpinAdapted::operatorfunctions::TensorMultiplyRightLeft(const StackSparseMatrix& a, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, double scale)
{
    // can be used for situation with different bra and ket
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


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = c.get_nonZeroBlocks();

  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int lQPrime = nonZeroBlocks[index].first.first, runcollectedQPrime = nonZeroBlocks[index].first.second;
    int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime];

    const std::vector<int>& rowinds = v.getActiveCols(lQPrime);
    for (int rrq = 0; rrq < rowinds.size(); ++rrq)  {
      int runcollectedQ = rowinds[rrq];
      int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ];
      if (a.allowed(rQ, rQPrime))
	MatrixMultiply (c.operator_element(lQPrime, runcollectedQPrime), 'n', a.operator_element(rQ, rQPrime), TransposeOf(a.conjugacy()),
			v.operator_element(lQPrime, runcollectedQ), 1.0, 0.0);
      
    }
  }

}

void SpinAdapted::operatorfunctions::TensorMultiplyLeft(const StackSparseMatrix& leftOp, const StackSparseMatrix& rightOp, const StackSparseMatrix& dotOp, const StackSparseMatrix& RIGHTOP, const StackSpinBlock *cblock, StackWavefunction& w, StackWavefunction* v, const SpinQuantum opQ, double scale)
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


  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[omprank].get_nonZeroBlocks();

  int OMPRANK = omprank;
  int quanta_thrds = dmrginp.quanta_thrds();

#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int lQ = nonZeroBlocks[index].first.first, runcollectedQ = nonZeroBlocks[index].first.second;
    int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];

    const std::vector<int>& colinds = leftOp.getActiveCols(lQ);
    for (int lrop=0; lrop <colinds.size(); lrop ++) {
      int lQPrime = colinds[lrop];
      
	const std::vector<int>& colinds2 = w.getActiveCols(lQPrime);
	for (int r = 0; r < colinds2.size(); r++) {
	  int runcollectedQPrime = colinds2[r];
	  int rQrQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

	  if (dotOp.allowed(dotQ, dotQPrime, RIGHTOP.conjugacy())) {
	      const std::vector<int>& rowinds = RIGHTOP.conjugacy() == 'n' ? rightOp.getActiveRows(rQrQPrime) : rightOp.getActiveCols(rQrQPrime);
	      assert (rowinds.size() == 1);
	      int rQPrime = rowinds[0];
	      pout << lQ<<"  "<<runcollectedQ<<"  "<<lQPrime<<"  "<<runcollectedQPrime<<"      "<<dotQ<<"  "<<rQ<<"  "<<dotQPrime<<"  "<<rQPrime<<endl;
	      if (dotOp.allowed(dotQ, dotQPrime, RIGHTOP.conjugacy()) && rightOp.allowed(rQ, rQPrime, RIGHTOP.conjugacy())) {
		pout << "allowed "<<endl;
		
		double factor = scale*leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);	      
		factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep(), 
					      leftOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), opQ.get_s().getirrep(),
					      lbraS->quanta[lQ].get_s().getirrep(), unCollectedrbraS->quanta[runcollectedQ].get_s().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_s().getirrep());
		factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep(), unCollectedrketS->quanta[runcollectedQPrime].get_symm().getirrep() , v[OMPRANK].get_deltaQuantum(0).get_symm().getirrep(), 
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
		if (fabs(factor*parity*scaleB) < TINY) continue; 
		
		MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftOp.conjugacy(), w.operator_element(lQPrime, runcollectedQPrime), 'n',  v[OMPRANK].operator_element(lQ, runcollectedQ),
				factor*parity*scaleB);
	      }
	  }
	}
    }
  }
  

}



void SpinAdapted::operatorfunctions::TensorMultiplysplitLeftsplitRight(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& RIGHTOP, 
								       const StackSparseMatrix& leftOp, const StackSparseMatrix& ldotOp, 
								       const StackSparseMatrix& rightOp, const StackSparseMatrix& rdotOp, 
								       const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, 
								       double factor)
{
  long starttime = globaltimer.totalwalltime();

  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* ldotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* ldotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* rdotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* rdotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;

  const int leftBraOpSz = cblock->get_leftBlock()->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int dotBraOpSz = cblock->get_leftBlock()->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int dotKetOpSz = cblock->get_leftBlock()->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks2 = c.get_nonZeroBlocks();
  int OMPRANK = omprank;

  int first = 0, second = 0;

  for (int cc=0; cc<nonZeroBlocks2.size(); cc++) {
    int luncollectedQPrime = nonZeroBlocks2[cc].first.first;
    int runcollectedQPrime = nonZeroBlocks2[cc].first.second;

    int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], ldotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];
    int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], rdotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

    std::vector<int> rowinds;
    if (RIGHTOP.conjugacy() == 'n')
      rowinds = rightOp.getActiveRows(rQPrime);
    else
      rowinds = rightOp.getActiveCols(rQPrime);

    for (int r = 0; r < rowinds.size(); r++) {
      int rQ = rowinds[r];

      double* data = Stackmem[omprank].allocate(unCollectedlketS->getquantastates(luncollectedQPrime)* rbraS->getquantastates(rQ));
      StackMatrix blocks;
      blocks.allocate(data, unCollectedlketS->getquantastates(luncollectedQPrime), rbraS->getquantastates(rQ) );

      
      //L(l,l') c(l', dl', r, dr')     
      for (int lQ=0; lQ<leftBraOpSz; lQ++) 
      if (leftOp.allowed(lQ, lQPrime, LEFTOP.conjugacy()))  {
	double scale = factor;

	for (int ldotQ=0; ldotQ<dotBraOpSz; ldotQ++) 
	  if (ldotOp.allowed(ldotQ, ldotQPrime, LEFTOP.conjugacy())) {
	    std::vector<int>& luncollectedQvec = unCollectedlbraS->quantaMap(lQ, ldotQ);
	
	    for (int luncollectedQindex=0; luncollectedQindex<luncollectedQvec.size(); luncollectedQindex++) {
	      int luncollectedQ = luncollectedQvec[luncollectedQindex]; 
	      
	      //v(l, dl,  r, dr)
	      for (int runcollectedQ =0; runcollectedQ<unCollectedrbraS->quanta.size(); runcollectedQ++) { 
		if (v[omprank].allowed(luncollectedQ, runcollectedQ) && 
		    allowed(LEFTOP.get_deltaQuantum(), unCollectedlbraS->quanta[luncollectedQ], unCollectedlketS->quanta[luncollectedQPrime]) &&
		    rQ == unCollectedrbraS->leftUnMapQuanta[runcollectedQ]) {
		  int rdotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];
		  
		  if (rdotOp.allowed(rdotQ, rdotQPrime, RIGHTOP.conjugacy()) && rightOp.allowed(rQ, rQPrime, RIGHTOP.conjugacy()) && 
		      allowed(RIGHTOP.get_deltaQuantum(), unCollectedrbraS->quanta[runcollectedQ], unCollectedrketS->quanta[runcollectedQPrime])) {
		    
		    double scale2 = getScaling(LEFTOP,  leftOp,   ldotOp, 
					       RIGHTOP, rightOp, rdotOp, 
					       unCollectedlbraS->quanta[luncollectedQ],      lbraS->quanta[lQ],      ldotbraS->quanta[ldotQ],
					       unCollectedlketS->quanta[luncollectedQPrime], lketS->quanta[lQPrime], ldotketS->quanta[ldotQPrime],
					       unCollectedrbraS->quanta[runcollectedQ],      rbraS->quanta[rQ],      rdotbraS->quanta[rdotQ],
					       unCollectedrketS->quanta[runcollectedQPrime], rketS->quanta[rQPrime], rdotketS->quanta[rdotQPrime]);
		    scale2 *= ldotOp.operator_element(ldotQ, ldotQPrime, LEFTOP.conjugacy())(1,1);
		    scale2 *= rdotOp.operator_element(rdotQ, rdotQPrime, RIGHTOP.conjugacy())(1,1);

		    if (fabs(scale*scale2) > TINY) {
		      //R(r, r') c(l',dl', r', dr') ->  b(', dl', r, dr')
		      if (rightOp.opName() == "OVERLAP") 
			copy(c.operator_element(luncollectedQPrime, runcollectedQPrime), blocks);
		      else
			MatrixMultiply (c.operator_element(luncollectedQPrime, runcollectedQPrime), 'n', rightOp.operator_element(rQ, rQPrime, RIGHTOP.conjugacy()), rightOp.conjugacy() == 'n' ? TransposeOf(RIGHTOP.conjugacy()) : RIGHTOP.conjugacy(), blocks, 1.0, 0.);	      
		      if (leftOp.opName() == "OVERLAP")
			MatrixScaleAdd(scale*scale2, blocks, v[omprank].operator_element(luncollectedQ, runcollectedQ));
		      else
			MatrixMultiply (leftOp.operator_element(lQ, lQPrime, LEFTOP.conjugacy()), leftOp.conjugacy() == 'n' ? LEFTOP.conjugacy() : TransposeOf(LEFTOP.conjugacy()), 
					blocks, 'n', v[omprank].operator_element(luncollectedQ, runcollectedQ), scale*scale2);	      
		      second++;
		    
		    }
		  }
		}
	      }
	      
	      
	    }
	  }
      }
      
      Stackmem[omprank].deallocate(blocks.Store(), blocks.Storage());
    }
  }
  //pout << first<<"  "<<second<<endl;
}




double SpinAdapted::operatorfunctions::getScaling(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& leftOp, const StackSparseMatrix& ldotOp, 
						  const StackSparseMatrix& RIGHTOP, const StackSparseMatrix& rightOp, const StackSparseMatrix& rdotOp, 
						  const SpinQuantum& luncollectedQ, const SpinQuantum& lQ, const SpinQuantum& ldotQ,
						  const SpinQuantum& luncollectedQPrime, const SpinQuantum& lQPrime, const SpinQuantum& ldotQPrime,
						  const SpinQuantum& runcollectedQ, const SpinQuantum& rQ, const SpinQuantum& rdotQ,
						  const SpinQuantum& runcollectedQPrime, const SpinQuantum& rQPrime, const SpinQuantum& rdotQPrime)

{
  double factor = RIGHTOP.get_scaling(runcollectedQ, runcollectedQPrime);	      
  factor *= LEFTOP.get_scaling(luncollectedQ, luncollectedQPrime);

  factor *= dmrginp.get_ninej()(luncollectedQPrime.get_s().getirrep(), runcollectedQPrime.get_s().getirrep() , 0,  
				LEFTOP.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), 0,
				luncollectedQ.get_s().getirrep(), runcollectedQ.get_s().getirrep() , 0);
	  
  double scaleBr = 1.0;
  
  if (RIGHTOP.conjugacy() == 'n') {
    scaleBr = dmrginp.get_ninej()(rQPrime.get_s().getirrep() , rdotQPrime.get_s().getirrep(), runcollectedQPrime.get_s().getirrep(), 
				 rightOp.get_spin().getirrep(), rdotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
				 rQ.get_s().getirrep() , rdotQ.get_s().getirrep(), runcollectedQ.get_s().getirrep());

    scaleBr *= rightOp.get_scaling(rQ, rQPrime);
    scaleBr *= rdotOp.get_scaling(rdotQ, rdotQPrime);
    if (rdotOp.get_fermion() && IsFermion(rQPrime)) scaleBr *= -1;
  }
  else {
    scaleBr = dmrginp.get_ninej()(rQ.get_s().getirrep() , rdotQ.get_s().getirrep(), runcollectedQ.get_s().getirrep(), 
				 rightOp.get_spin().getirrep(), rdotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
				 rQPrime.get_s().getirrep() , rdotQPrime.get_s().getirrep(), runcollectedQPrime.get_s().getirrep());
    scaleBr *= rightOp.get_scaling(rQPrime, rQ);
    scaleBr *= rdotOp.get_scaling(rdotQPrime, rdotQ);
    if (rdotOp.get_fermion() && IsFermion(rQ)) scaleBr *= -1;
  }

  double scaleBl = 1.0;
  
  if (LEFTOP.conjugacy() == 'n') {
    scaleBl = dmrginp.get_ninej()(lQPrime.get_s().getirrep() , ldotQPrime.get_s().getirrep(), luncollectedQPrime.get_s().getirrep(), 
				 leftOp.get_spin().getirrep(), ldotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
				 lQ.get_s().getirrep() , ldotQ.get_s().getirrep(), luncollectedQ.get_s().getirrep());

    scaleBl *= leftOp.get_scaling(lQ, lQPrime);
    scaleBl *= ldotOp.get_scaling(ldotQ, ldotQPrime);
    if (ldotOp.get_fermion() && IsFermion(lQPrime)) scaleBl *= -1;
  }
  else {
    scaleBl = dmrginp.get_ninej()(lQ.get_s().getirrep() , ldotQ.get_s().getirrep(), luncollectedQ.get_s().getirrep(), 
				 leftOp.get_spin().getirrep(), ldotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
				 lQPrime.get_s().getirrep() , ldotQPrime.get_s().getirrep(), luncollectedQPrime.get_s().getirrep());
    scaleBl *= leftOp.get_scaling(lQPrime, lQ);
    scaleBl *= ldotOp.get_scaling(ldotQPrime, ldotQ);
    if (ldotOp.get_fermion() && IsFermion(lQ)) scaleBl *= -1;
  }

  
  int parity = RIGHTOP.get_fermion() && IsFermion(luncollectedQPrime) ? -1 : 1;
  //pout << factor<<"  "<<scaleBl<<"  "<<scaleBr<<"  "<<parity<<endl;
  return factor*scaleBr*scaleBl*parity;
  
}


void SpinAdapted::operatorfunctions::multiplyDotRightElement(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& leftOp, const StackSparseMatrix& dotOp, 
							     const StackSparseMatrix& rightOp, const StackMatrix& cMat, const StackMatrix& rightOpmat, StackMatrix& vMat,  
							     const SpinQuantum& luncollectedQ, const SpinQuantum& lQ, const SpinQuantum& dotQ, const SpinQuantum& rQ,
							     const SpinQuantum& luncollectedQPrime, const SpinQuantum& lQPrime, const SpinQuantum& dotQPrime, const SpinQuantum& rQPrime,
							     double scale)
{
  double factor = scale*LEFTOP.get_scaling(luncollectedQ, luncollectedQPrime);	      
  //c and v have spin 0, assume singlet embedding
  factor *= dmrginp.get_ninej()(luncollectedQPrime.get_s().getirrep(), rQPrime.get_s().getirrep() , 0,  
				LEFTOP.get_spin().getirrep(), rightOp.get_spin().getirrep(), 0,
				luncollectedQ.get_s().getirrep(), rQ.get_s().getirrep() , 0);
	  
  double scaleB = 1.0;
  
  if (LEFTOP.conjugacy() == 'n') {
    scaleB = dmrginp.get_ninej()(lQPrime.get_s().getirrep() , dotQPrime.get_s().getirrep(), luncollectedQPrime.get_s().getirrep(), 
				 leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
				 lQ.get_s().getirrep() , dotQ.get_s().getirrep(), luncollectedQ.get_s().getirrep());

    scaleB *= leftOp.get_scaling(lQ, lQPrime);
    scaleB *= dotOp.get_scaling(dotQ, dotQPrime);
    if (dotOp.get_fermion() && IsFermion(lQPrime)) scaleB *= -1;
  }
  else {
    scaleB = dmrginp.get_ninej()(lQ.get_s().getirrep() , dotQ.get_s().getirrep(), luncollectedQ.get_s().getirrep(), 
				 leftOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), LEFTOP.get_spin().getirrep(),
				 lQPrime.get_s().getirrep() , dotQPrime.get_s().getirrep(), luncollectedQPrime.get_s().getirrep());
    scaleB *= leftOp.get_scaling(lQPrime, lQ);
    scaleB *= dotOp.get_scaling(dotQPrime, dotQ);
    if (dotOp.get_fermion() && IsFermion(lQ)) scaleB *= -1;
  }
  
  int parity = rightOp.get_fermion() && IsFermion(luncollectedQPrime) ? -1 : 1;
  factor *=  rightOp.get_scaling(rQ, rQPrime);
  if (fabs(factor*parity*scaleB) < TINY) return; 

  MatrixMultiply (cMat, 'n', rightOpmat, TransposeOf(rightOp.conjugacy()), vMat, factor*parity*scaleB);	      
  
  
}



void SpinAdapted::operatorfunctions::multiplyDotRight(const StackSparseMatrix& LEFTOP, const StackSparseMatrix& leftop, const StackSparseMatrix& dotop, 
						      StackSparseMatrix& rightop, std::vector<StackMatrix>& lopCmat,
						      StackWavefunction* v,  const StackSpinBlock* cblock, int luncollectedQPrime, int rQPrime, double scale)
{
  const boost::shared_ptr<StateInfo> unCollectedlbraS = cblock->get_braStateInfo().leftStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedlketS = cblock->get_ketStateInfo().leftStateInfo->unCollectedStateInfo;
  const StateInfo* lbraS = cblock->get_leftBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* lketS = cblock->get_leftBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_leftBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_leftBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  StackSpinBlock *leftBlock = cblock->get_leftBlock()->get_leftBlock(), *dotBlock = cblock->get_leftBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock();

  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  int lQPrime = unCollectedlketS->leftUnMapQuanta[luncollectedQPrime], dotQPrime = unCollectedlketS->rightUnMapQuanta[luncollectedQPrime];

  std::vector<int> rowinds;
  if (LEFTOP.conjugacy() == 'n')
    rowinds = leftop.getActiveRows(lQPrime);
  else
    rowinds = leftop.getActiveCols(lQPrime);

  for (int r = 0; r < rowinds.size(); r++) {
    int lQ = rowinds[r];

    bool deallocate = rightop.memoryUsed() == 0 ? true : false;
    
    //R(r,r') c(l, d', r')  
    for (int rQ=0; rQ<rightBraOpSz; rQ++) {
      if (allowed(rightop.get_deltaQuantum(), rbraS->quanta[rQ], rketS->quanta[rQPrime])) {
	
	//v(l, d, r)
	for (int luncollectedQ =0; luncollectedQ<unCollectedlbraS->quanta.size(); luncollectedQ++) { 
	if (v[omprank].allowed(luncollectedQ, rQ) ) {

	  if (lQ == unCollectedlbraS->leftUnMapQuanta[luncollectedQ] && 
	      dotop.allowed(unCollectedlbraS->rightUnMapQuanta[luncollectedQ], dotQPrime, LEFTOP.conjugacy())) {
	    int lQ = unCollectedlbraS->leftUnMapQuanta[luncollectedQ], dotQ = unCollectedlbraS->rightUnMapQuanta[luncollectedQ];
	    
	    double* ropdata = deallocate ? Stackmem[omprank].allocate(rbraS->quantaStates[rQ]* rketS->quantaStates[rQPrime]) :  0;
	    StackMatrix ropm;
	    if (rightop.conjugacy() == 'n') 
	      ropm = deallocate ? StackMatrix(ropdata, rbraS->quantaStates[rQ], rketS->quantaStates[rQPrime]) : StackMatrix();
	    else
	      ropm = deallocate ? StackMatrix(ropdata, rketS->quantaStates[rQPrime], rbraS->quantaStates[rQ]) : StackMatrix();

	    memset(ropdata, 0, ropm.Storage() * sizeof(double));

	    rightop.build(ropm, rQ, rQPrime, *rightBlock);

	    operatorfunctions::multiplyDotRightElement(LEFTOP, leftop, dotop, rightop, 
						       lopCmat[r], ropm, v[omprank].operator_element(luncollectedQ, rQ),
						       unCollectedlbraS->quanta[luncollectedQ], lbraS->quanta[lQ], dotbraS->quanta[dotQ], rbraS->quanta[rQ],
						       unCollectedlketS->quanta[luncollectedQPrime], lketS->quanta[lQPrime], dotketS->quanta[dotQPrime], 
						       rketS->quanta[rQPrime], scale*dotop.operator_element(dotQ, dotQPrime, LEFTOP.conjugacy())(1,1));
	    if (deallocate) Stackmem[omprank].deallocate(ropdata, ropm.Storage());
	  }
	}
	}
      }
    }
  }
}


void SpinAdapted::operatorfunctions::multiplyDotLeftElement(const StackSparseMatrix& RIGHTOP, const StackSparseMatrix& rightOp, const StackSparseMatrix& dotOp, 
							     const StackSparseMatrix& leftOp, const StackMatrix& cMat, const StackMatrix& leftOpmat, StackMatrix& vMat,  
							     const SpinQuantum& lQ, const SpinQuantum& runcollectedQ, const SpinQuantum& rQ, const SpinQuantum& dotQ,
							     const SpinQuantum& lQPrime, const SpinQuantum& runcollectedQPrime, const SpinQuantum& rQPrime, const SpinQuantum& dotQPrime,
							     double scale)
{
  double factor = scale*RIGHTOP.get_scaling(runcollectedQ, runcollectedQPrime);	      
  //c and v have spin 0, assume singlet embedding
  factor *= dmrginp.get_ninej()(lQPrime.get_s().getirrep(), runcollectedQPrime.get_s().getirrep() , 0,  
				leftOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(), 0,
				lQ.get_s().getirrep(), runcollectedQ.get_s().getirrep() , 0);
	  
  double scaleB = 1.0;
  
  if (RIGHTOP.conjugacy() == 'n') {
    scaleB = dmrginp.get_ninej()(rQPrime.get_s().getirrep() , dotQPrime.get_s().getirrep(), runcollectedQPrime.get_s().getirrep(), 
				 rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
				 rQ.get_s().getirrep() , dotQ.get_s().getirrep(), runcollectedQ.get_s().getirrep());

    scaleB *= rightOp.get_scaling(rQ, rQPrime);
    scaleB *= dotOp.get_scaling(dotQ, dotQPrime);
    if (dotOp.get_fermion() && IsFermion(rQPrime)) scaleB *= -1;
  }
  else {
    scaleB = dmrginp.get_ninej()(rQ.get_s().getirrep() , dotQ.get_s().getirrep(), runcollectedQ.get_s().getirrep(), 
				 rightOp.get_spin().getirrep(), dotOp.get_spin().getirrep(), RIGHTOP.get_spin().getirrep(),
				 rQPrime.get_s().getirrep() , dotQPrime.get_s().getirrep(), runcollectedQPrime.get_s().getirrep());
    scaleB *= rightOp.get_scaling(rQPrime, rQ);
    scaleB *= dotOp.get_scaling(dotQPrime, dotQ);
    if (dotOp.get_fermion() && IsFermion(rQ)) scaleB *= -1;
  }
  
  int parity = RIGHTOP.get_fermion() && IsFermion(lQPrime) ? -1 : 1;
  factor *=  leftOp.get_scaling(lQ, lQPrime);
  if (fabs(factor*parity*scaleB) < TINY) return; 

  MatrixMultiply (leftOpmat, leftOp.conjugacy(), cMat, 'n', vMat, factor*parity*scaleB);	      
  
  
}

void SpinAdapted::operatorfunctions::multiplyDotLeft(const StackSparseMatrix& RIGHTOP, const StackSparseMatrix& rightop, const StackSparseMatrix& dotop, 
						      StackSparseMatrix& leftop, std::vector<StackMatrix>& ropCmat,
						      StackWavefunction* v,  const StackSpinBlock* cblock, int lQPrime, int runcollectedQPrime, double scale)
{
  const boost::shared_ptr<StateInfo> unCollectedrbraS = cblock->get_braStateInfo().rightStateInfo->unCollectedStateInfo;
  const boost::shared_ptr<StateInfo> unCollectedrketS = cblock->get_ketStateInfo().rightStateInfo->unCollectedStateInfo;
  const StateInfo* rbraS = cblock->get_rightBlock()->get_braStateInfo().leftStateInfo; 
  const StateInfo* rketS = cblock->get_rightBlock()->get_ketStateInfo().leftStateInfo;
  const StateInfo* dotbraS = cblock->get_rightBlock()->get_braStateInfo().rightStateInfo;
  const StateInfo* dotketS = cblock->get_rightBlock()->get_ketStateInfo().rightStateInfo;
  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;


  StackSpinBlock *leftBlock = cblock->get_leftBlock(), *dotBlock = cblock->get_rightBlock()->get_rightBlock();
  StackSpinBlock *rightBlock = cblock->get_rightBlock()->get_leftBlock();

  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();

  int rQPrime = unCollectedrketS->leftUnMapQuanta[runcollectedQPrime], dotQPrime = unCollectedrketS->rightUnMapQuanta[runcollectedQPrime];

  std::vector<int> rowinds;
  if (RIGHTOP.conjugacy() == 'n')
    rowinds = rightop.getActiveRows(rQPrime);
  else
    rowinds = rightop.getActiveCols(rQPrime);

  for (int r = 0; r < rowinds.size(); r++) {
    int rQ = rowinds[r];

    bool deallocate = leftop.memoryUsed() == 0 ? true : false;
    
    //L(l,l') c(l', dl'r, d')  
    for (int lQ=0; lQ<leftBraOpSz; lQ++) {
      if (allowed(leftop.get_deltaQuantum(), lbraS->quanta[lQ], lketS->quanta[lQPrime])) {
	
	//v(l, r, d)
	for (int runcollectedQ =0; runcollectedQ<unCollectedrbraS->quanta.size(); runcollectedQ++) { 
	if (v[omprank].allowed(lQ, runcollectedQ) ) {

	  if (rQ == unCollectedrbraS->leftUnMapQuanta[runcollectedQ] && 
	      dotop.allowed(unCollectedrbraS->rightUnMapQuanta[runcollectedQ], dotQPrime, RIGHTOP.conjugacy())) {
	    int rQ = unCollectedrbraS->leftUnMapQuanta[runcollectedQ], dotQ = unCollectedrbraS->rightUnMapQuanta[runcollectedQ];
	    
	    double* lopdata = deallocate ? Stackmem[omprank].allocate(lbraS->quantaStates[lQ]* lketS->quantaStates[lQPrime]) :  0;
	    StackMatrix lopm;
	    if (leftop.conjugacy() == 'n') 
	      lopm = deallocate ? StackMatrix(lopdata, lbraS->quantaStates[lQ], lketS->quantaStates[lQPrime]) : StackMatrix();
	    else
	      lopm = deallocate ? StackMatrix(lopdata, lketS->quantaStates[lQPrime], lbraS->quantaStates[lQ]) : StackMatrix();

	    memset(lopdata, 0, lopm.Storage() * sizeof(double));

	    leftop.build(lopm, lQ, lQPrime, *leftBlock);

	    operatorfunctions::multiplyDotLeftElement(RIGHTOP, rightop, dotop, leftop, 
						      ropCmat[r], lopm, v[omprank].operator_element(lQ, runcollectedQ),
						      lbraS->quanta[lQ], unCollectedrbraS->quanta[runcollectedQ], rbraS->quanta[rQ], dotbraS->quanta[dotQ],
						      lketS->quanta[lQPrime], unCollectedrketS->quanta[runcollectedQPrime], rketS->quanta[rQPrime],  
						      dotketS->quanta[dotQPrime], scale*dotop.operator_element(dotQ, dotQPrime, RIGHTOP.conjugacy())(1,1));
	    if (deallocate) Stackmem[omprank].deallocate(lopdata, lopm.Storage());
	  }
	}
	}
      }
    }
  }
}

void SpinAdapted::operatorfunctions::braTensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix &a, const StackSpinBlock *cblock, StackWavefunction &c, StackWavefunction &v, double scale, int num_thrds)
{
  //It get a result of <\Psi|O. 
  //It is similar to Transposeview(O)|\Psi>
  //However, spin coupling coefficients is different for transpose.
  //It is convenient for npdm with intermediate.
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *lketS = cblock->get_ketStateInfo().leftStateInfo;
  const StateInfo* rbraS = cblock->get_braStateInfo().rightStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  assert (cblock->get_leftBlock() == ablock || cblock->get_rightBlock() == ablock);
  if (cblock->get_leftBlock() == ablock)
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int lQ = 0; lQ < leftKetOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftBraOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQPrime, lQ))
              {
		const StackMatrix& aop = a.operator_element(lQPrime, lQ);
		  for (int rQ = 0; rQ < rightBraOpSz; ++rQ) 
		    if (c.allowed(lQPrime, rQ) && v.allowed(lQ, rQ))
		    {
                      double fac=scale;
		      fac *= dmrginp.get_ninej()(lbraS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						   (-a.get_spin()).getirrep(), 0, (-a.get_spin()).getirrep(),
						   lketS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		      fac *= Symmetry::spatial_ninej(lbraS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), c.get_symm().getirrep(), 
					   (-a.get_symm()).getirrep(), 0, (-a.get_symm()).getirrep(),
					   lketS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		      fac *= a.get_scaling(lbraS->quanta[lQPrime], lketS->quanta[lQ]);
		      MatrixMultiply (aop, TransposeOf(a.conjugacy()), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }

              }
	  }
      }
      }
    }
  else
    {
      //#pragma omp parallel default(shared)  num_threads(num_thrds)
      {
	//#pragma omp for schedule(dynamic)
      for (int rQ = 0; rQ < rightKetOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightBraOpSz; ++rQPrime)
	  if (a.allowed(rQPrime, rQ))
	    {
	      const StackMatrix& aop = a.operator_element(rQ, rQPrime);
	      for (int lQPrime = 0; lQPrime < leftBraOpSz; ++lQPrime) 
		if (v.allowed(lQPrime, rQ) && c.allowed(lQPrime, rQPrime)) {
                  double fac = scale;
		  fac *= dmrginp.get_ninej()(lbraS->quanta[lQPrime].get_s().getirrep(), rbraS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					       0, (-a.get_spin()).getirrep(), (-a.get_spin()).getirrep(),
					       lbraS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  fac *= Symmetry::spatial_ninej(lbraS->quanta[lQPrime].get_symm().getirrep() , rbraS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				      0, (-a.get_symm()).getirrep(), (-a.get_symm()).getirrep(),
				      lbraS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  fac *= a.get_scaling(rbraS->quanta[rQPrime], rketS->quanta[rQ]);
		  double parity = a.get_fermion() && IsFermion(lbraS->quanta[lQPrime]) ? -1 : 1;

		  MatrixMultiply (c.operator_element(lQPrime, rQPrime), c.conjugacy(),
				  aop, a.conjugacy(), v.operator_element(lQPrime, rQ), fac*parity);
		}

	    }
      }
      }
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
  if (op1.conjugacy() == 'n') {
    for (int lQ = 0; lQ< op2.nrows(); lQ++)
      for (int rQ = 0; rQ<op2.ncols(); rQ++)
        if (op2.allowed(lQ, rQ) && op1.allowed(lQ,rQ))
        {
	        double factor = op1.get_scaling(s.quanta[lQ], s.quanta[rQ]);
	        MatrixScaleAdd(scaleV*factor, op1.operator_element(lQ,rQ), op2.operator_element(lQ,rQ));
        }
  } else {
    for (int lQ = 0; lQ< op2.nrows(); lQ++)
      for (int rQ = 0; rQ<op2.ncols(); rQ++)
        if (op2.allowed(lQ, rQ) && op1.allowed(lQ,rQ)) {
          double scaling = getStandAlonescaling(op1.get_deltaQuantum(0), b.get_braStateInfo().quanta[lQ], b.get_ketStateInfo().quanta[rQ]);
          int nrows = op2.operator_element(lQ,rQ).Nrows();
          int ncols = op2.operator_element(lQ,rQ).Ncols();
          for (int row = 0; row < nrows; ++row)
            DAXPY(ncols, scaling*scaleV, op1.operator_element(lQ,rQ).Store() + row, nrows, &(op2.operator_element(lQ,rQ)(row+1, 1)), 1);
      }
  }
}

void SpinAdapted::operatorfunctions::OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const StackSparseMatrix& op1, StackSparseMatrix& op2, const std::vector<int>& rows, const std::vector<int>& cols) {
  const StateInfo& s = b.get_stateInfo();
  assert(rows.size() == cols.size());
  if (op1.conjugacy() == 'n') {
    for (int i = 0; i < rows.size(); ++i) {
      int lQ = rows.at(i), rQ = cols.at(i);
        if (op2.allowed(lQ, rQ) && op1.allowed(lQ,rQ))
        {
	        double factor = op1.get_scaling(s.quanta[lQ], s.quanta[rQ]);
	        MatrixScaleAdd(scaleV*factor, op1.operator_element(lQ,rQ), op2.operator_element(lQ,rQ));
        }
    }
  } else {
    for (int i = 0; i < rows.size(); ++i) {
      int lQ = rows.at(i), rQ = cols.at(i);
        if (op2.allowed(lQ, rQ) && op1.allowed(lQ,rQ)) {
          double scaling = getStandAlonescaling(op1.get_deltaQuantum(0), b.get_braStateInfo().quanta[lQ], b.get_ketStateInfo().quanta[rQ]);
          int nrows = op2.operator_element(lQ,rQ).Nrows();
          int ncols = op2.operator_element(lQ,rQ).Ncols();
          for (int row = 0; row < nrows; ++row)
            DAXPY(ncols, scaling*scaleV, op1.operator_element(lQ,rQ).Store() + row, nrows, &(op2.operator_element(lQ,rQ)(row+1, 1)), 1);
      }
    }
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

  
void SpinAdapted::operatorfunctions::TensorMultiply(const StackSparseMatrix& a, const StackSparseMatrix& b, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, bool aIsLeftOp, double scale)
{
  const int leftBraOpSz = brastateinfo->leftStateInfo->quanta.size ();
  const int leftKetOpSz = ketstateinfo->leftStateInfo->quanta.size ();
  const int rightBraOpSz = brastateinfo->rightStateInfo->quanta.size ();
  const int rightKetOpSz = ketstateinfo->rightStateInfo->quanta.size ();

  const StateInfo* lbraS = brastateinfo->leftStateInfo, *rbraS = brastateinfo->rightStateInfo;
  const StateInfo* lketS = ketstateinfo->leftStateInfo, *rketS = ketstateinfo->rightStateInfo;

  const char conjC = (aIsLeftOp) ? 'n' : 't';

  const StackSparseMatrix& leftOp = (conjC == 'n') ? a : b; // an ugly hack to support the release memory optimisation
  const StackSparseMatrix& rightOp = (conjC == 'n') ? b : a;
  const char leftConj = (conjC == 'n') ? a.conjugacy() : b.conjugacy();
  const char rightConj = (conjC == 'n') ? b.conjugacy() : a.conjugacy();
  const std::vector< std::pair<std::pair<int, int>, StackMatrix> >& nonZeroBlocks = v[omprank].get_nonZeroBlocks();

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


	      

void SpinAdapted::operatorfunctions::TensorMultiply(const StackSparseMatrix a, const StateInfo *brastateinfo, const StateInfo *ketstateinfo, const StackWavefunction& c, StackWavefunction& v, const SpinQuantum dQ, bool left, double scale)
{
  //Calculate O_{l or r} |\Psi> without building big block.
  const StateInfo* lbraS = brastateinfo->leftStateInfo, *lketS = ketstateinfo->leftStateInfo;
  const StateInfo* rbraS = brastateinfo->rightStateInfo, *rketS = ketstateinfo->rightStateInfo;
  const int leftBraOpSz = brastateinfo->leftStateInfo->quanta.size ();
  const int leftKetOpSz = ketstateinfo->leftStateInfo->quanta.size ();
  const int rightBraOpSz = brastateinfo->rightStateInfo->quanta.size ();
  const int rightKetOpSz = ketstateinfo->rightStateInfo->quanta.size ();

  if (left)
    {
      for (int lQ = 0; lQ < leftBraOpSz; ++lQ) {
	for (int lQPrime = 0; lQPrime < leftKetOpSz; ++lQPrime)
	  {
	    if (a.allowed(lQ, lQPrime))
	      {
		const StackMatrix& aop = a.operator_element(lQ, lQPrime);
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
		      MatrixMultiply (aop, a.conjugacy(), c.operator_element(lQPrime, rQ), c.conjugacy(),
				      v.operator_element(lQ, rQ), fac);
		    }
		
	      }
	  }
      }
    }
  else
    {
      for (int rQ = 0; rQ < rightBraOpSz; ++rQ) {
	for (int rQPrime = 0; rQPrime < rightKetOpSz; ++rQPrime)
	  if (a.allowed(rQ, rQPrime))
	    {
	      const StackMatrix& aop = a.operator_element(rQ, rQPrime);
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
				  aop, TransposeOf(a.conjugacy()), v.operator_element(lQPrime, rQ), fac*parity);
		}
	      
	    }
	
      }
    }
}


