/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "Stackspinblock.h"
#include "Stackdensity.h"
#include "Stackwavefunction.h"
#include "operatorloops.h"
#include "operatorfunctions.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "stackguess_wavefunction.h"
#include "distribute.h"
#include <boost/format.hpp>
#include "pario.h"


namespace SpinAdapted{
using namespace operatorfunctions;

void StackDensityMatrix::makedensitymatrix(std::vector<StackWavefunction>& wave_solutions, StackSpinBlock &big, 
				      const std::vector<double> &wave_weights, const double noise, const double additional_noise, bool warmup)
{

  //the density Matrix should already be allocated
  for(int i=0;i<wave_weights.size()&& mpigetrank() == 0;++i) {
    makedensitymatrix(wave_solutions[i], big, wave_weights[i]);
  }

#ifndef SERIAL
  boost::mpi::communicator world;
  //broadcast the data
  MPI_Bcast(this->get_data(), this->memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif

  if(noise > NUMERICAL_ZERO) {
    
    /* check normalisation */
    double norm = 0.0;
    for(int lQ=0;lQ<nrows();++lQ)
      if(allowed(lQ,lQ))
	for(int i=0;i<(*this)(lQ,lQ).Nrows();++i)
	  norm += (*this)(lQ,lQ)(i+1,i+1);
    p2out << "\t\t\t norm before modification " << norm << endl;
    
    
    int nroots = wave_solutions.size();
    
#ifndef SERIAL
    boost::mpi::communicator world;
    boost::mpi::broadcast(calc, nroots, 0);
#endif

    {

      double* backupData, *noiseMatrix;
      long requiredData;
      //make a backup of the actual density Matrix, only on the main root node
      if (mpigetrank() == 0) {
	requiredData = getRequiredMemory(*big.get_leftBlock(), get_deltaQuantum());
	backupData = Stackmem[omprank].allocate(requiredData);
	noiseMatrix = Stackmem[omprank].allocate(requiredData);
	memset(noiseMatrix, 0, requiredData * sizeof(double));
	DCOPY(requiredData, this->get_data(), 1, &backupData[0], 1);
      }
      mcheck("just before noise");
      StackWavefunction *wptr = &wave_solutions[0];

      for (int i=0; i<nroots; i++) {
	this->Clear();
	if (mpigetrank() == 0)
	  wptr = &wave_solutions[i];

#ifndef SERIAL
	MPI_Bcast(wptr->get_data(), wptr->memoryUsed(), MPI_DOUBLE, 0, Calc);
#endif
      
	this->add_onedot_noise(*wptr, big, (1.0*noise)/nroots);

	//add the noise from this wavefunction back to the noiseMatrix on the 0th proc
	if (mpigetrank() == 0)
	  DAXPY(requiredData, 1.0, this->get_data(), 1, noiseMatrix, 1);
      }
      mcheck("just after noise");
      //copy back the noiseMatrix to "this", and deallocate the noiseMatrix
      if (mpigetrank() == 0) {
	DCOPY(requiredData, noiseMatrix, 1, this->get_data(), 1);
	Stackmem[omprank].deallocate(noiseMatrix, requiredData);
      }
	
      if (mpigetrank() == 0) {
	norm = 0.0;
	for(int lQ=0;lQ<nrows();++lQ)
	  if(this->allowed(lQ,lQ))
	    for(int i=0;i<(*this)(lQ,lQ).Nrows();++i)
	      norm += (*this)(lQ,lQ)(i+1,i+1);
	//add the noise density matrix to the current matrix
	if (fabs(norm) > 1.0e-8)
	  DAXPY(requiredData, noise/norm, this->get_data(), 1, &backupData[0], 1);
	
	//copy it back to this
	DCOPY(requiredData, &backupData[0], 1, this->get_data(), 1);
	p2out << "\t\t\t norm after modification " << trace(*this) << endl;
      }

      if (mpigetrank() == 0)
	Stackmem[omprank].deallocate(backupData, requiredData);
    }
    
  }

}
  
void StackDensityMatrix::makedensitymatrix(StackWavefunction& wave_solution, StackSpinBlock &big, 
				      const double &wave_weight)
{
  MultiplyWithOwnTranspose (wave_solution, *this, wave_weight);  
}



void StackDensityMatrix::add_twodot_noise(const StackSpinBlock &big, const double noise)
{
  pout << "Twodot noise is not supported with StackDensityMatrix";
  exit(0);

}

StackDensityMatrix& StackDensityMatrix::operator+=(const StackDensityMatrix& other)
{
  DAXPY(totalMemory, 1.0, (const_cast<StackDensityMatrix&>(other)).get_data(), 1, get_data(), 1);

  return *this;
}



class onedot_noise_f 
{
private:
  const StackWavefunction& wavefunction;
  StackDensityMatrix*& dm;
  const StackSpinBlock& big; 
  const double scale;
  bool distributed;
  bool synced;
public:
  onedot_noise_f(StackDensityMatrix*& dm_, const StackWavefunction& wavefunction_, const StackSpinBlock& big_, const double scale_)
    : distributed(false), synced(true), wavefunction(wavefunction_), dm(dm_), big(big_), scale(scale_) { }
  
  void operator()(const boost::shared_ptr<StackSparseMatrix> op) const {
    vector<SpinQuantum> wQ = wavefunction.get_deltaQuantum();
    vector<SpinQuantum> oQ = op->get_deltaQuantum();
    vector<IrrepSpace> vec = wQ[0].get_symm() + oQ[0].get_symm();
    vector<SpinSpace> spinvec = wQ[0].get_s()+oQ[0].get_s();
    if (dmrginp.hamiltonian() == BCS) {
      for (int n = 0; n <= dmrginp.effective_molecule_quantum().get_n(); ++n) {
        bool valid_cre = false, valid_des = false;
        for (int k = 0; k < wQ.size(); ++k) {
          for (int l = 0; l < oQ.size(); ++l) {
            if (wQ[k].get_n() + oQ[l].get_n() == n) valid_cre = true;
            if (!big.get_leftBlock()->has(DES) && wQ[k].get_n() - oQ[l].get_n() == n) valid_des = true;
          }
        }
        if (!valid_cre && !valid_des) continue;
        if (!op->memoryUsed()) {
          op->allocate(big.get_leftBlock()->get_braStateInfo(), big.get_leftBlock()->get_ketStateInfo());
	        op->build(*big.get_leftBlock());
        }
        for (int j = 0; j < vec.size(); ++j) {
          for (int i = 0; i < spinvec.size(); ++i) {
            if (valid_cre) {
	            SpinQuantum q = SpinQuantum(n, spinvec[i], vec[j]);
	            StackWavefunction opxwave;
	            opxwave.initialise(std::vector<SpinQuantum>(1,q), *big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, wavefunction.get_onedot());
	            opxwave.set_onedot(wavefunction.get_onedot());
	            opxwave.Clear();
	            TensorMultiply(big.get_leftBlock(), *op, &big, const_cast<StackWavefunction&> (wavefunction), opxwave, dmrginp.molecule_quantum(), 1.0);
	            double norm = DotProduct(opxwave, opxwave);
	            if (abs(norm) > NUMERICAL_ZERO) {
		            Scale(1./sqrt(norm), opxwave);
		            MultiplyWithOwnTranspose (opxwave, dm[omprank], scale);  
	            }
	            opxwave.deallocate();
            }
            if (valid_des) {
	            SpinQuantum q = SpinQuantum(n, spinvec[i], vec[j]);
	            StackWavefunction opxwave2;
		          opxwave2.initialise(std::vector<SpinQuantum>(1,q), *big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, wavefunction.get_onedot());
		          opxwave2.set_onedot(wavefunction.get_onedot());
		          opxwave2.Clear();
		          TensorMultiply(big.get_leftBlock(), Transpose(*op), &big, const_cast<StackWavefunction&> (wavefunction), opxwave2, dmrginp.molecule_quantum(), 1.0);
		          double norm = DotProduct(opxwave2, opxwave2);
		          if (abs(norm) >NUMERICAL_ZERO) {
		            Scale(1./sqrt(norm), opxwave2);
		            MultiplyWithOwnTranspose (opxwave2, dm[omprank], scale);  
		            //MultiplyProduct(opxwave2, Transpose(opxwave2), dm[0], scale);
		          } 
		          opxwave2.deallocate();
            }
          }
        }
      }
	    op->deallocate();      
    } else {
      for (int k=0; k<wQ.size(); ++k)
      for (int l=0; l<oQ.size(); ++l)
	    for (int j=0; j<vec.size(); j++)
	    for (int i=0; i<spinvec.size(); i++) {
	      SpinQuantum q = SpinQuantum(wQ[k].get_n()+oQ[l].get_n(), spinvec[i], vec[j]);
	      op->allocate(big.get_leftBlock()->get_braStateInfo(), big.get_leftBlock()->get_ketStateInfo());
	      op->build(*big.get_leftBlock());
	      StackWavefunction opxwave;
	      opxwave.initialise(std::vector<SpinQuantum>(1,q), *big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, wavefunction.get_onedot());
	      opxwave.set_onedot(wavefunction.get_onedot());
	      opxwave.Clear();
	      TensorMultiply(big.get_leftBlock(), *op, &big, const_cast<StackWavefunction&> (wavefunction), opxwave, dmrginp.molecule_quantum(), 1.0);
	      double norm = DotProduct(opxwave, opxwave);
	      if (abs(norm) > NUMERICAL_ZERO) {
		      Scale(1./sqrt(norm), opxwave);
		      MultiplyWithOwnTranspose (opxwave, dm[omprank], scale);  
		      //MultiplyProduct(opxwave, Transpose(opxwave), dm[0], scale);
	      }
	      opxwave.deallocate();
	      
	      //this block has explicit transpose operators, so dont do this step
	      if (!big.get_leftBlock()->has(DES)) {
	        q = SpinQuantum(wQ[k].get_n()-oQ[l].get_n(), spinvec[i], vec[j]);
		      StackWavefunction opxwave2; //= Wavefunction(q, &big, wavefunction.get_onedot());
		      opxwave2.initialise(std::vector<SpinQuantum>(1,q), *big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, wavefunction.get_onedot());
		      opxwave2.set_onedot(wavefunction.get_onedot());
		      opxwave2.Clear();
		      TensorMultiply(big.get_leftBlock(), Transpose(*op), &big, const_cast<StackWavefunction&> (wavefunction), opxwave2, dmrginp.molecule_quantum(), 1.0);
		      double norm = DotProduct(opxwave2, opxwave2);
		      if (abs(norm) >NUMERICAL_ZERO) {
		        Scale(1./sqrt(norm), opxwave2);
		        MultiplyWithOwnTranspose (opxwave2, dm[omprank], scale);  
		        //MultiplyProduct(opxwave2, Transpose(opxwave2), dm[0], scale);
		      } 
		      opxwave2.deallocate();
	      }
	      
	      op->deallocate();
	    }
    }
  }
};



// accumulates into dm
void StackDensityMatrix::add_onedot_noise(StackWavefunction& wave_solution, StackSpinBlock& big, bool act2siteops)
{

  StackSpinBlock* leftBlock = big.get_leftBlock();
  //p1out << "\t\t\t Modifying density matrix " << endl;


  StackDensityMatrix* dm; 
  initiateMultiThread(this, dm, numthrds);

  onedot_noise_f onedot_noise(dm, wave_solution, big, 1.);

  std::vector<boost::shared_ptr<StackSparseMatrix> >  allops;


  if (leftBlock->has(CRE)) {
    for (int i=0; i<leftBlock->get_op_array(CRE).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(CRE).get_local_element(i).size(); j++) {
	allops.push_back(leftBlock->get_op_array(CRE).get_local_element(i)[j]);
      }
  }
  if (leftBlock->has(DES)) {
    for (int i=0; i<leftBlock->get_op_array(DES).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(DES).get_local_element(i).size(); j++) {
	allops.push_back(leftBlock->get_op_array(DES).get_local_element(i)[j]);
      }
  }

  //use overlap only when bra and ket are different i.e. when the block has des operator
  if (leftBlock->has(DES)&&leftBlock->has(OVERLAP) && mpigetrank() == 0) {
    for (int i=0; i<leftBlock->get_op_array(OVERLAP).get_size(); i++)
      for (int j=0; j<leftBlock->get_op_array(OVERLAP).get_local_element(i).size(); j++) {
	allops.push_back(leftBlock->get_op_array(OVERLAP).get_local_element(i)[j]);
      }
  }
  
  if (dmrginp.hamiltonian() != HUBBARD) {
    
    if (leftBlock->has(CRE_CRE)) {
      for (int i=0; i<leftBlock->get_op_array(CRE_CRE).get_size(); i++)
	for (int j=0; j<leftBlock->get_op_array(CRE_CRE).get_local_element(i).size(); j++) {
	  allops.push_back(leftBlock->get_op_array(CRE_CRE).get_local_element(i)[j]);
	}
      
      if (leftBlock->has(CRE_DES)) {
	for (int i=0; i<leftBlock->get_op_array(CRE_DES).get_size(); i++)
	  for (int j=0; j<leftBlock->get_op_array(CRE_DES).get_local_element(i).size(); j++) {
	    allops.push_back(leftBlock->get_op_array(CRE_DES).get_local_element(i)[j]);
	  }
      }
    } 
    else if (leftBlock->has(DES_DESCOMP)) {
      for (int i=0; i<leftBlock->get_op_array(DES_DESCOMP).get_size(); i++)
	for (int j=0; j<leftBlock->get_op_array(DES_DESCOMP).get_local_element(i).size(); j++) {
	  allops.push_back(leftBlock->get_op_array(DES_DESCOMP).get_local_element(i)[j]);
	}
      
      if (leftBlock->has(CRE_DESCOMP)) {
	for (int i=0; i<leftBlock->get_op_array(CRE_DESCOMP).get_size(); i++)
	  for (int j=0; j<leftBlock->get_op_array(CRE_DESCOMP).get_local_element(i).size(); j++) {
	    allops.push_back(leftBlock->get_op_array(CRE_DESCOMP).get_local_element(i)[j]);
	  }
      }
    }
    if (leftBlock->has(DES_DES)) {
      for (int i=0; i<leftBlock->get_op_array(DES_DES).get_size(); i++)
	for (int j=0; j<leftBlock->get_op_array(DES_DES).get_local_element(i).size(); j++) {
	  allops.push_back(leftBlock->get_op_array(DES_DES).get_local_element(i)[j]);
	}
      
      if (leftBlock->has(DES_CRE)) {
	for (int i=0; i<leftBlock->get_op_array(DES_CRE).get_size(); i++)
	  for (int j=0; j<leftBlock->get_op_array(DES_CRE).get_local_element(i).size(); j++) {
	    allops.push_back(leftBlock->get_op_array(DES_CRE).get_local_element(i)[j]);
	  }
      }
    }
    else if (leftBlock->has(CRE_CRECOMP)) {
      for (int i=0; i<leftBlock->get_op_array(CRE_CRECOMP).get_size(); i++)
	for (int j=0; j<leftBlock->get_op_array(CRE_CRECOMP).get_local_element(i).size(); j++) {
	  allops.push_back(leftBlock->get_op_array(CRE_CRECOMP).get_local_element(i)[j]);
	}
      
      if (leftBlock->has(DES_CRECOMP)) {
	for (int i=0; i<leftBlock->get_op_array(DES_CRECOMP).get_size(); i++)
	  for (int j=0; j<leftBlock->get_op_array(DES_CRECOMP).get_local_element(i).size(); j++) {
	    allops.push_back(leftBlock->get_op_array(DES_CRECOMP).get_local_element(i)[j]);
	  }   
      }   
    }
    
    
  }

  SplitStackmem();
  dmrginp.tensormultiply->start();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i<allops.size(); i++)  {
    onedot_noise(allops[i]);
  }
  dmrginp.tensormultiply->stop();  
  MergeStackmem();

  accumulateMultiThread(this, dm, numthrds);
  distributedaccumulate(*this);

  
}


}


  
