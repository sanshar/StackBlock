/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/
#include <iostream>
#include <newmatio.h>
#include <newmat.h>
#include <newmatap.h>
#include <cmath>
#include <vector>
#include "davidson.h"
#include "Stackwavefunction.h"
#include "linear.h"
#include "pario.h"
#include "global.h"
#include "MatrixBLAS.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#ifndef SERIAL
#include <boost/mpi.hpp>
#include "mpi.h"
#endif

using namespace boost;
using namespace std;

const double EPS=1.e-20;



void SpinAdapted::Linear::precondition(StackWavefunction& op, double e, DiagonalMatrix& diagonal, double levelshift)
{
  if (!mpigetrank())
  {
    int index = 1;
    for (int lQ = 0; lQ < op.nrows (); ++lQ)
      for (int rQ = 0; rQ < op.ncols (); ++rQ)
	if (op.allowed(lQ, rQ))
	  for (int lQState = 0; lQState < op.operator_element(lQ, rQ).Nrows (); ++lQState)
	    for (int rQState = 0; rQState < op.operator_element(lQ, rQ).Ncols (); ++rQState)
	      {
		if (fabs(e - diagonal(index)) > 1.e-12)		
		  op.operator_element(lQ, rQ)(lQState+1, rQState+1) /= (e - diagonal(index)+levelshift);
		++index;
	      }
  }
}



void SpinAdapted::Linear::olsenPrecondition(StackWavefunction& op, StackWavefunction& C0, double e, DiagonalMatrix& diagonal, double levelshift)
{

  StackWavefunction C0copy; 
  C0copy.initialise(C0);
  DCOPY(C0copy.memoryUsed(), C0.get_data(), 1, C0copy.get_data(), 1); 
  precondition(C0copy, e, diagonal, levelshift);
  double numerator = DotProduct(C0copy, op);
  double denominator = DotProduct(C0, C0copy);
  ScaleAdd(-numerator/denominator, C0, op);
  precondition(op, e, diagonal, levelshift);
  C0copy.deallocate();
}

void SpinAdapted::Linear::block_davidson(vector<StackWavefunction>& b, DiagonalMatrix& h_diag, double normtol, const bool &warmUp, Davidson_functor& h_multiply, bool& useprecond, int currentRoot, std::vector<StackWavefunction> &lowerStates)
{

#ifndef SERIAL
  mpi::communicator world;
#endif
  int iter = 0;
  double levelshift = 0.0;
  int nroots = dmrginp.setStateSpecific() ? 1 : dmrginp.nroots();
  //normalise all the guess roots

  double timer = globaltimer.totalwalltime();

  bool orthogonalSpace = true;
  if(mpigetrank() == 0) {
    for(int i=0;i<nroots;++i)
      {
	for(int j=0;j<i;++j)
	  {
	    double overlap = DotProduct(b[j], b[i]);
	    ScaleAdd(-overlap, b[j], b[i]);
	  }
	Normalise(b[i]);
      }
  
  
    //if we are doing state specific, lowerstates are lower energy states
    if (lowerStates.size() != 0) {
      for (int i=0; i<lowerStates.size(); i++) {
	double overlap = DotProduct(b[0], lowerStates[i]);
	if (DotProduct(lowerStates[i], lowerStates[i]) > NUMERICAL_ZERO) 
	  ScaleAdd(-overlap/DotProduct(lowerStates[i], lowerStates[i]), lowerStates[i], b[0]);
      }
      if (DotProduct(b[0], b[0]) > NUMERICAL_ZERO)
	Normalise(b[0]);
      else {
	b[0].Randomise();
	Normalise(b[0]);
	orthogonalSpace = false;
      }
    }

  }

#ifndef SERIAL
    mpi::broadcast(calc, orthogonalSpace, 0);
#endif
  if (!orthogonalSpace) {
    return;
  }

  vector<StackWavefunction> sigma(dmrginp.deflation_max_size());
  sigma[0].initialise(b[0]);
  for (int i=1; i<b.size(); i++) 
    if (mpigetrank() == 0)
      sigma[i].initialise(b[i]);

  StackWavefunction r; 
  if (mpigetrank() == 0)
    r.initialise(b[0]);

  if (mpigetrank() == 0) {
    printf("\t\t %15s  %5s  %15s  %9s  %10s %10s \n", "iter", "Root", "Energy", "Error", "Time", "FLOPS");
  }
  int sigmasize=0, bsize= currentRoot == -1 ? dmrginp.nroots() : 1;
  int converged_roots = 0;
  int maxiter = h_diag.Ncols() - lowerStates.size();
  while(true && iter < 100)
  {
    //p3out << "\t\t\t Davidson Iteration :: " << iter << endl;
    
    ++iter;
    dmrginp.hmultiply -> start();
    
#ifndef SERIAL
    mpi::broadcast(calc, sigmasize, 0);
    mpi::broadcast(calc, bsize, 0);
#endif
    //multiply all guess vectors with hamiltonian c = Hv

    for(int i=sigmasize;i<bsize;++i) {
      StackWavefunction* sigmaptr=&sigma[0], *bptr=&b[0];
	  
      if (mpigetrank() == 0) {
	sigmaptr = &sigma[i];
	bptr = &b[i];
      }

#ifndef SERIAL
      MPI_Bcast(bptr->get_data(), bptr->memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif
      sigmaptr->Clear();
      
      h_multiply(*bptr, *sigmaptr);
      sigmasize++;
    }
    dmrginp.hmultiply -> stop();


    DiagonalMatrix subspace_eigenvalues;

    double currentEnergy ;
    if (mpigetrank() == 0) {
      Matrix subspace_h(bsize, bsize);
      //#pragma omp parallel for schedule(dynamic)
      for (int i = 0; i < bsize; ++i) {
	for (int j = 0; j <= i; ++j) {
	  subspace_h.element(i, j) = DotProduct(b[i], sigma[j]);
	  subspace_h.element(j, i) = subspace_h.element(i, j);
	}
      }

      Matrix alpha;
      diagonalise(subspace_h, subspace_eigenvalues, alpha);

      currentEnergy = subspace_eigenvalues(converged_roots+1, converged_roots+1);
      //for (int i = 1; i <= nroots; ++i)
      //p3out << "\t\t\t " << iter << " ::  " << subspace_eigenvalues(i,i) << endl;

      //now calculate the ritz vectors which are approximate eigenvectors
      vector<StackWavefunction> tmp(bsize);
      for (int i=0; i<bsize; i++) {
	tmp[i].initialise(b[i]);
	tmp[i].copyData(b[i]);
	//copy(b[i].get_operatorMatrix(), tmp[i].get_operatorMatrix());
      }
      for (int i = 0; i < bsize; ++i)
	Scale(alpha.element(i, i), b[i]);

      //#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < bsize; ++j) {
      for (int i = 0; i < bsize; ++i) 
      if (i != j)
	ScaleAdd(alpha.element(i,j), tmp[i], b[j]);
      //DAXPY(b[0].memoryUsed(), alpha.element(i, j),  tmp[i], 1, b[j].get_data(), 1);
      }

      for (int i=0; i<bsize; i++) 
	tmp[i].copyData(sigma[i]);
      //copy(sigma[i].get_operatorMatrix(), tmp[i].get_operatorMatrix());


      for (int i = 0; i < bsize; ++i)
	Scale(alpha.element(i, i), sigma[i]);

      //#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < bsize; ++j) {
      for (int i = 0; i < bsize; ++i) 
      if (i != j)
	ScaleAdd(alpha.element(i,j), tmp[i], sigma[j]);
      }
      //DAXPY(sigma[j].memoryUsed(), alpha.element(i, j),  tmp[i], 1, sigma[j].get_data(), 1);

      for (int i=bsize-1; i>-1; i--)
	tmp[i].deallocate();
	
      // build residual 
      for (int i=0; i<converged_roots; i++) {
	r.copyData(sigma[i]);
	//copy(sigma[i].get_operatorMatrix(), r.get_operatorMatrix());
	ScaleAdd(-subspace_eigenvalues(i+1), b[i], r);
	double rnorm = DotProduct(r,r);
	if (rnorm > normtol) {
	  converged_roots = i;
	  p3out << "\t\t\t going back to converged root "<<i<<"  "<<rnorm<<" > "<<normtol<<endl;
	  continue;
	}
      }

      r.copyData(sigma[converged_roots]);
      //copy(sigma[converged_roots].get_operatorMatrix(), r.get_operatorMatrix());
      ScaleAdd(-subspace_eigenvalues(converged_roots+1), b[converged_roots], r);

      if (lowerStates.size() != 0) {
	for (int i=0; i<lowerStates.size(); i++) {
	  double overlap = DotProduct(r, lowerStates[i]);
	  ScaleAdd(-overlap/DotProduct(lowerStates[i], lowerStates[i]), lowerStates[i], r);
	}
      }
    }


    double rnorm;
    if (mpigetrank() == 0) {
      rnorm = DotProduct(r,r);  
      double totalFlops = 0.;
      for (int thrd=0; thrd<numthrds; thrd++) {
	totalFlops += dmrginp.matmultFlops[thrd];
	dmrginp.matmultFlops[thrd] = 0.0;
      }
      printf("\t\t %15i  %5i  %15.8f  %9.2e %10.2f (s)  %10.3e\n", iter, converged_roots, currentEnergy, rnorm, globaltimer.totalwalltime()-timer, totalFlops);

      timer = globaltimer.totalwalltime();
    }

#ifndef SERIAL
    mpi::broadcast(calc, converged_roots, 0);
    mpi::broadcast(calc, rnorm, 0);
#endif

    if (useprecond && mpigetrank() == 0)
      olsenPrecondition(r, b[converged_roots], subspace_eigenvalues(converged_roots+1), h_diag, levelshift);

    //p3out << "\t \t \t residual :: " << rnorm << endl;
    if (rnorm < normtol)
    {
      p3out << "\t\t\t Converged root " << converged_roots << endl;
      
      ++converged_roots;
      if (converged_roots == nroots)
      {
	if (mpigetrank() == 0) {
	  for (int i = 0; i < min((int)(bsize), h_diag.Ncols()); ++i)
	    h_diag.element(i) = subspace_eigenvalues.element(i);
	}
	break;
      }
    }
    else if (mpigetrank() == 0)
    {
      if(bsize >= dmrginp.deflation_max_size())
      {
	p3out << "\t\t\t Deflating block Davidson...\n";
	bsize = dmrginp.deflation_min_size();
	sigmasize = dmrginp.deflation_min_size();
      }
      for (int j = 0; j < bsize; ++j)
      {
	//Normalize
	double normalization = DotProduct(r, r);
	Scale(1./sqrt(normalization), r);
	
	double overlap = DotProduct(r, b[j]);
	ScaleAdd(-overlap, b[j], r);
      }

      //if we are doing state specific, lowerstates has lower energy states
      if (lowerStates.size() != 0) {
	for (int i=0; i<lowerStates.size(); i++) {
	  double overlap = DotProduct(r, lowerStates[i]);
	  ScaleAdd(-overlap/DotProduct(lowerStates[i], lowerStates[i]), lowerStates[i], r);
	  //ScaleAdd(-overlap, lowerStates[i], r);
	}
      }
      //double tau2 = DotProduct(r,r);

      Normalise(r);
      
      b[bsize].copyData(r);
      //copy(r.get_operatorMatrix(), b[bsize].get_operatorMatrix());
      bsize++;
    }
    
    
  }

  if (mpigetrank() == 0)
    r.deallocate();
  for (int i=b.size()-1; i>0; i--) 
    if (mpigetrank() == 0)
      sigma[i].deallocate();
  sigma[0].deallocate();
}


void makeOrthogonalToLowerStates(StackWavefunction& targetState, std::vector<StackWavefunction>& lowerStates) {
  for (int i=1; i<lowerStates.size(); i++) {
    double overlap2 = pow(DotProduct(lowerStates[i], lowerStates[i]), 0.5);
    if (fabs(overlap2) > NUMERICAL_ZERO) { 
      ScaleAdd(-DotProduct(targetState, lowerStates[i])/overlap2, 
	       lowerStates[i], targetState);
    }
  }
}


double SpinAdapted::Linear::ConjugateGradient(StackWavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<StackWavefunction>& lowerStates)
{
  setbuf(stdout, NULL);
  p3out.precision (12);
  int iter = 0, maxIter = 10;
  double levelshift = 0.0, overlap2 = 0.0, oldError=0.0, functional=0.0, Error=0.0;

  StackWavefunction& targetState = lowerStates[0];
  if (mpigetrank() == 0) {
    makeOrthogonalToLowerStates(targetState, lowerStates);
    makeOrthogonalToLowerStates(xi, lowerStates);
  }

#ifndef SERIAL
  mpi::communicator world;
  MPI_Bcast(xi.get_data(), xi.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

  StackWavefunction pi, ri; 
  ri.initialise(xi);
  ri.Clear();
  pi.initialise(ri);
  pi.Clear();
  h_multiply(xi, ri);  

  //Check if we should even perform CG or just exit with a zero vector.
  bool doCG = true;
  if (mpigetrank() == 0) {
    StackWavefunction ricopy; ricopy.initialise(ri); ricopy.Randomise();
    makeOrthogonalToLowerStates(ricopy, lowerStates);

    if (abs(DotProduct(ricopy, targetState)) < NUMERICAL_ZERO) {
      pout << "The problem is ill posed or the initial guess is very bad "<<DotProduct(ricopy, targetState)<<endl;
      doCG = false;
    }
    ricopy.deallocate();
  }
#ifndef SERIAL
  mpi::broadcast(calc, doCG, 0);
#endif
  if (!doCG) {
    xi.Clear();
    int success = 0;

    functional = 0.0;
    pi.deallocate();
    ri.deallocate();
    return functional;
  }

  if (mpigetrank() == 0) {
    ScaleAdd(-1.0, targetState, ri);
    Scale(-1.0, ri);
    
    makeOrthogonalToLowerStates(ri, lowerStates);    

    DCOPY(ri.memoryUsed(), ri.get_data(), 1, pi.get_data(), 1); 

    oldError = DotProduct(ri, ri);
    printf("\t\t\t %15s  %15s  %15s\n", "iter", "Functional", "Error");
  }


#ifndef SERIAL
  mpi::broadcast(calc, Error, 0);
  mpi::broadcast(calc, oldError, 0);
  mpi::broadcast(calc, functional, 0);
#endif

  if (oldError < normtol) {
    if (mpigetrank() == 0) {
      functional = DotProduct(xi, ri) - 2.*DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e\n", 0, functional, oldError);
    }
#ifndef SERIAL
    mpi::broadcast(calc, functional, 0);
#endif
    pi.deallocate();
    ri.deallocate();
    return functional;
  }

#ifndef SERIAL
  MPI_Bcast(pi.get_data(), pi.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

  StackWavefunction Hp; 
  Hp.initialise(pi); 

  while(true) {
    Hp.Clear();



    h_multiply(pi, Hp);

    if (mpigetrank() == 0) {
      makeOrthogonalToLowerStates(Hp, lowerStates);

      double alpha = oldError/DotProduct(pi, Hp);
      
      ScaleAdd(alpha, pi, xi);
      ScaleAdd(-alpha, Hp, ri);
      
      Error = DotProduct(ri, ri);
      functional = -DotProduct(xi, ri) - DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e\n", iter, functional, Error);
    }

    //exit(0);
#ifndef SERIAL
    mpi::broadcast(calc, Error, 0);
    mpi::broadcast(calc, functional, 0);
#endif

    if (Error < normtol || iter >maxIter) {
      Hp.deallocate();
      pi.deallocate();
      ri.deallocate();
      return functional;
    }
    else {      
      if (mpigetrank() == 0) {
	double beta = Error/oldError;
	oldError = Error;
	ScaleAdd(1.0/beta, ri, pi);
	Scale(beta, pi);	
	makeOrthogonalToLowerStates(pi, lowerStates);
      }
      iter ++;
    }
  }
}

double SpinAdapted::Linear::MinResMethod(StackWavefunction& xi, double normtol, Davidson_functor& h_multiply, std::vector<StackWavefunction>& lowerStates)
{
  setbuf(stdout, NULL);
  int iter = 0, maxIter = 20;
  double levelshift = 0.0, overlap2 = 0.0, oldError=0.0, functional=0.0, Error=0.0;

  StackWavefunction& targetState = lowerStates[0];
  if (mpigetrank() == 0) {
    makeOrthogonalToLowerStates(targetState, lowerStates);
    makeOrthogonalToLowerStates(xi, lowerStates);
  }

#ifndef SERIAL
  mpi::communicator world;
  MPI_Bcast(xi.get_data(), xi.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

  StackWavefunction pi, ri; 
  ri.initialise(xi);
  ri.Clear();
  
  h_multiply(xi, ri);  
  
  //Check if we should even perform CG or just exit with a zero vector.
  bool doCG = true;
  if (mpigetrank() == 0) {
    StackWavefunction ricopy; ricopy.initialise(ri);ricopy.Randomise();
    makeOrthogonalToLowerStates(ricopy, lowerStates);

    if (abs(DotProduct(ricopy, targetState)) < NUMERICAL_ZERO) {
      pout << "The problem is ill posed or the initial guess is very bad "<<DotProduct(ricopy, targetState)<<endl;
      doCG = false;
    }
    ricopy.deallocate();
  }
#ifndef SERIAL
    mpi::broadcast(calc, doCG, 0);
#endif
  if (!doCG) {
    xi.Clear();
    int success = 0;

    functional = 0.0;
    ri.deallocate();
    return functional;
  }

  if (mpigetrank() == 0) {
    ScaleAdd(-1.0, targetState, ri);
    Scale(-1.0, ri);
    
    makeOrthogonalToLowerStates(ri, lowerStates);
    pi.initialise(ri);

    DCOPY(ri.memoryUsed(), ri.get_data(), 1, pi.get_data(), 1); 
    //pi = ri;

    oldError = DotProduct(ri, ri);
    printf("\t\t\t %15s  %15s  %15s\n", "iter", "Functional", "Error");
  }

#ifndef SERIAL
  mpi::broadcast(calc, oldError, 0);
#endif
  
  if (oldError < normtol) {
    if (mpigetrank() == 0) {
      functional = -DotProduct(xi, ri) -  DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e\n", 0, functional, oldError);
    }
#ifndef SERIAL
    mpi::broadcast(calc, functional, 0);
#endif
    if (mpigetrank() == 0) pi.deallocate();
    ri.deallocate();
    return functional;
  }

#ifndef SERIAL
  MPI_Bcast(ri.get_data(), ri.memoryUsed(), MPI_DOUBLE, 0, Calc);
    //mpi::broadcast(world, ri, 0);
#endif
  
  double betaNumerator = 0, betaDenominator = 0;
  StackWavefunction Hr, Hp; Hr.initialise(ri); Hr.Clear();
  h_multiply(ri, Hr);
  betaDenominator = DotProduct(ri, Hr);

  if (mpigetrank() == 0) {
    makeOrthogonalToLowerStates(Hr, lowerStates);
    Hp.initialise(Hr);
    DCOPY(Hr.memoryUsed(), Hr.get_data(), 1, Hp.get_data(), 1); 
  }

  while(true) {

    if (mpigetrank() == 0) {
      double alpha = DotProduct(ri, Hr)/DotProduct(Hp, Hp);
      
      ScaleAdd(alpha, pi, xi);
      ScaleAdd(-alpha, Hp, ri);
      
      Error = DotProduct(ri, ri);

      functional = -DotProduct(xi, targetState);
      printf("\t\t\t %15i  %15.8e  %15.8e \n", iter, functional, Error);
    }

#ifndef SERIAL
    mpi::broadcast(calc, Error, 0);
    mpi::broadcast(calc, functional, 0);
    MPI_Bcast(ri.get_data(), ri.memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

    if (Error < normtol || iter >maxIter) {
      if (mpigetrank() == 0) Hp.deallocate();
      Hr.deallocate();
      if (mpigetrank() == 0) pi.deallocate();
      ri.deallocate();
      return functional;
    }
    else {      
      Hr.Clear();
      h_multiply(ri, Hr);
      if (mpigetrank() == 0) {
	makeOrthogonalToLowerStates(Hr, lowerStates);

	betaNumerator = DotProduct(ri, Hr);
	double beta = betaNumerator/betaDenominator;
	betaDenominator = betaNumerator;

	ScaleAdd(1./beta, ri, pi);
	Scale(beta, pi);

	ScaleAdd(1./beta, Hr, Hp);
	Scale(beta, Hp);
	
      }
      iter ++;
    }
  }
}

