/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include "operatorfunctions.h"
#include "Stackspinblock.h"
#include "Stackwavefunction.h"
#include "couplingCoeffs.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pario.h"


void SpinAdapted::operatorfunctions::TensorMultiply(const StackSpinBlock *ablock, const StackSparseMatrix& a, const StackSparseMatrix& b, const StackSpinBlock *cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum opQ, double scale)
{
  long starttime = globaltimer.totalwalltime();
  dmrginp.tensormultiply->start();
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

#pragma omp parallel for schedule(dynamic)
  for (int index = 0; index<nonZeroBlocks.size(); index++) {
    int lQ = nonZeroBlocks[index].first.first, rQ = nonZeroBlocks[index].first.second;

    for (int rQPrime=0; rQPrime <rightKetOpSz; rQPrime ++) {
      if (rightOp.allowed(rQ, rQPrime)) {

	const std::vector<int>& rowinds = c.getActiveRows(rQPrime);
	for (int l = 0; l < rowinds.size(); l++) {
	  int lQPrime = rowinds[l];
	  if (leftOp.allowed(lQ, lQPrime) ) {

	    long sum = 0;
	    double data[lbraS->getquantastates(lQ)* rketS->getquantastates(rQPrime)];
	    StackMatrix m(data, lbraS->getquantastates(lQ), rketS->getquantastates(rQPrime));
	    
	    double factor = leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftConj, c.operator_element(lQPrime, rQPrime), 'n',
			    m, factor, 0.);	      

	    {
	      double factor = scale;
	      
	      factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
					    leftOp.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
					    lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v[omprank].get_deltaQuantum(0).get_s().getirrep());
	      factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
						leftOp.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
						lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v[omprank].get_symm().getirrep());
	      int parity = rightOp.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
	      factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
	      MatrixMultiply (m, 'n', rightOp.operator()(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), v[omprank].operator_element(lQ, rQ), factor*parity);
	    }
	  }
	}
      }
    }
  }


  dmrginp.tensormultiply->stop();
  

}


void SpinAdapted::operatorfunctions::MultiplyWithOwnTranspose(const StackSparseMatrix& a, StackSparseMatrix& c, Real scale)
{
  if (fabs(scale) < TINY) return;
  const int aSz = a.nrows();
  const int aSzPrime = a.ncols();
  
    assert (c.nrows() == a.nrows() &&
	    c.ncols() == a.nrows());
    
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

void SpinAdapted::operatorfunctions::Product (const StackSpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, double scale)
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





/*
void SpinAdapted::operatorfunctions::braTensorMultiply(const SpinBlock *ablock, const Baseoperator<Matrix>& a, const SpinBlock *cblock, Wavefunction& c, Wavefunction& v, double scale, int num_thrds)
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
		const Matrix& aop = a.operator_element(lQPrime, lQ);
		  for (int rQ = 0; rQ < rightKetOpSz; ++rQ) 
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
	      const Matrix& aop = a.operator_element(rQ, rQPrime);
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
*/



/*
void SpinAdapted::operatorfunctions::TensorMultiply(const StackSpinBlock *ablock, const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, const StackSpinBlock *cblock, Wavefunction& c, Wavefunction& v, const SpinQuantum opQ, double scale)
{
  // can be used for situation with different bra and ket
  const int leftBraOpSz = cblock->get_leftBlock()->get_braStateInfo().quanta.size ();
  const int leftKetOpSz = cblock->get_leftBlock()->get_ketStateInfo().quanta.size ();
  const int rightBraOpSz = cblock->get_rightBlock()->get_braStateInfo().quanta.size ();
  const int rightKetOpSz = cblock->get_rightBlock()->get_ketStateInfo().quanta.size ();

  const StateInfo* lbraS = cblock->get_braStateInfo().leftStateInfo, *rbraS = cblock->get_braStateInfo().rightStateInfo;
  const StateInfo* lketS = cblock->get_ketStateInfo().leftStateInfo, *rketS = cblock->get_ketStateInfo().rightStateInfo;

  const char conjC = (cblock->get_leftBlock() == ablock) ? 'n' : 't';

  const Baseoperator<Matrix>& leftOp = (conjC == 'n') ? a : b; // an ugly hack to support the release memory optimisation
  const Baseoperator<Matrix>& rightOp = (conjC == 'n') ? b : a;
  const char leftConj = (conjC == 'n') ? a.conjugacy() : b.conjugacy();
  const char rightConj = (conjC == 'n') ? b.conjugacy() : a.conjugacy();

  Wavefunction u;
  u.resize(leftBraOpSz*leftKetOpSz, rightKetOpSz);

  int totalmem =0;

  {
    for (int lQrQPrime = 0; lQrQPrime<leftBraOpSz*rightKetOpSz; ++lQrQPrime)
    {
      int rQPrime = lQrQPrime%rightKetOpSz, lQ = lQrQPrime/rightKetOpSz;
	for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
	  if (leftOp.allowed(lQ, lQPrime) && c.allowed(lQPrime, rQPrime))
	  {
	    int lindex = lQ*leftKetOpSz+lQPrime;
	    u.allowed(lindex, rQPrime) = true;
            
	    u(lindex,rQPrime).ReSize(lbraS->getquantastates(lQ), rketS->getquantastates(rQPrime));
	    double factor = leftOp.get_scaling(lbraS->quanta[lQ], lketS->quanta[lQPrime]);
	    MatrixMultiply (leftOp.operator_element(lQ, lQPrime), leftConj, c.operator_element(lQPrime, rQPrime), 'n',
			    u.operator_element(lindex, rQPrime), factor, 0.);	      

	  }
    }
  }

  pout << "after first step in tensormultiply"<<endl;
      mcheck("before davidson but after all blocks are built");

  {
    for (int lQrQ = 0; lQrQ<leftBraOpSz*rightBraOpSz; ++lQrQ)
    {
      int rQ = lQrQ%rightBraOpSz, lQ=lQrQ/rightBraOpSz;
	if (v.allowed(lQ, rQ))
	  for (int rQPrime = 0; rQPrime < rightKetOpSz; rQPrime++)
	    if (rightOp.allowed(rQ, rQPrime))
	      for (int lQPrime = 0; lQPrime < leftKetOpSz; lQPrime++)
		if (leftOp.allowed(lQ, lQPrime) && u.allowed(lQ*leftKetOpSz+lQPrime, rQPrime))
		{
		  int lindex = lQ*leftKetOpSz+lQPrime;
		  double factor = scale;

		  factor *= dmrginp.get_ninej()(lketS->quanta[lQPrime].get_s().getirrep(), rketS->quanta[rQPrime].get_s().getirrep() , c.get_deltaQuantum(0).get_s().getirrep(), 
						leftOp.get_spin().getirrep(), rightOp.get_spin().getirrep(), opQ.get_s().getirrep(),
						lbraS->quanta[lQ].get_s().getirrep(), rbraS->quanta[rQ].get_s().getirrep() , v.get_deltaQuantum(0).get_s().getirrep());
		  factor *= Symmetry::spatial_ninej(lketS->quanta[lQPrime].get_symm().getirrep() , rketS->quanta[rQPrime].get_symm().getirrep(), c.get_symm().getirrep(), 
				       leftOp.get_symm().getirrep(), rightOp.get_symm().getirrep(), opQ.get_symm().getirrep(),
				       lbraS->quanta[lQ].get_symm().getirrep() , rbraS->quanta[rQ].get_symm().getirrep(), v.get_symm().getirrep());
		  int parity = rightOp.get_fermion() && IsFermion(lketS->quanta[lQPrime]) ? -1 : 1;
		  factor *=  rightOp.get_scaling(rbraS->quanta[rQ], rketS->quanta[rQPrime]);
		  MatrixMultiply (u.operator_element(lindex, rQPrime), 'n',
				  rightOp(rQ, rQPrime), TransposeOf(rightOp.conjugacy()), v.operator_element(lQ, rQ), factor*parity);
		}
    }
  }
	      
}
*/
void SpinAdapted::operatorfunctions::OperatorScaleAdd(double scaleV, const StackSpinBlock& b, const Baseoperator<Matrix>& op1, Baseoperator<Matrix>& op2)
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
void SpinAdapted::operatorfunctions::MultiplyProduct(const Baseoperator<Matrix>& a, const Baseoperator<Matrix>& b, Baseoperator<Matrix>& c, Real scale)
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

  

	      


