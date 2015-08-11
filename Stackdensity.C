/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

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
  MPI::COMM_WORLD.Bcast(this->get_data(), this->memoryUsed(), MPI_DOUBLE, 0);
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
    boost::mpi::broadcast(world, nroots, 0);
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
	MPI::COMM_WORLD.Bcast(wptr->get_data(), wptr->memoryUsed(), MPI_DOUBLE, 0);
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
  const int num_threads;
  opTypes optype, optype2;
  bool distributed;
  bool synced;
public:
  onedot_noise_f(StackDensityMatrix*& dm_, const StackWavefunction& wavefunction_, const StackSpinBlock& big_, const double scale_, const int num_threads_)
    : distributed(false), synced(true), wavefunction(wavefunction_), dm(dm_), big(big_), scale(scale_), num_threads(num_threads_) { }
  
  void set_opType(const opTypes &optype_)
  {
    optype = optype_;
    distributed = !big.get_leftBlock()->get_op_array(optype).is_local();
    if(distributed) synced = false;
  }
  void operator()(const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec) const
  {
    if ((mpigetrank() == 0 || distributed))// && op.get_deltaQuantum().particleNumber != 0)
    {
      for (int opind=0; opind<opvec.size(); opind++) {
	StackSparseMatrix& op = *opvec[opind];
#ifndef SERIAL
	boost::mpi::communicator world;
	if (op.get_orbs().size() == 1 && op.get_orbs()[0]%world.size() != mpigetrank())
	  continue;
#endif	  
	vector<SpinQuantum> wQ = wavefunction.get_deltaQuantum();
	vector<SpinQuantum> oQ = op.get_deltaQuantum();
	vector<IrrepSpace> vec = wQ[0].get_symm() + oQ[0].get_symm();
	vector<SpinSpace> spinvec = wQ[0].get_s()+oQ[0].get_s();
	for (int k=0; k<wQ.size(); ++k)
	  for (int l=0; l<oQ.size(); ++l)
	    for (int j=0; j<vec.size(); j++)
	      for (int i=0; i<spinvec.size(); i++)
		{
		  SpinQuantum q = SpinQuantum(wQ[k].get_n()+oQ[l].get_n(), spinvec[i], vec[j]);
		  op.allocate(big.get_leftBlock()->get_braStateInfo(), big.get_leftBlock()->get_ketStateInfo());
		  op.build(*big.get_leftBlock());
		  if (dmrginp.hamiltonian() != BCS || q.get_n() <= dmrginp.effective_molecule_quantum().get_n()) {
		    StackWavefunction opxwave;
		    opxwave.initialise(std::vector<SpinQuantum>(1,q), *big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, wavefunction.get_onedot());
		    opxwave.set_onedot(wavefunction.get_onedot());
		    opxwave.Clear();
		    TensorMultiply(big.get_leftBlock(), op, &big, const_cast<StackWavefunction&> (wavefunction), opxwave, dmrginp.molecule_quantum(), 1.0);
		    double norm = DotProduct(opxwave, opxwave);
		    if (abs(norm) > NUMERICAL_ZERO) {
		      Scale(1./sqrt(norm), opxwave);
		      MultiplyWithOwnTranspose (opxwave, dm[0], scale);  
		      //MultiplyProduct(opxwave, Transpose(opxwave), dm[0], scale);
		    }
		    opxwave.deallocate();
		  }

		  //this block has explicit transpose operators, so dont do this step
		  if (!big.get_leftBlock()->has(DES)) {
		    q = SpinQuantum(wQ[k].get_n()-oQ[l].get_n(), spinvec[i], vec[j]);
		    if (dmrginp.hamiltonian() != BCS || q.get_n() >= 0) {
		      StackWavefunction opxwave2; //= Wavefunction(q, &big, wavefunction.get_onedot());
		      opxwave2.initialise(std::vector<SpinQuantum>(1,q), *big.get_braStateInfo().leftStateInfo, *big.get_ketStateInfo().rightStateInfo, wavefunction.get_onedot());
		      opxwave2.set_onedot(wavefunction.get_onedot());
		      opxwave2.Clear();
		      op.set_conjugacy('t');
		      TensorMultiply(big.get_leftBlock(), op,&big, const_cast<StackWavefunction&> (wavefunction), opxwave2, dmrginp.molecule_quantum(), 1.0);
		      op.set_conjugacy('n');
		      double norm = DotProduct(opxwave2, opxwave2);
		      if (abs(norm) >NUMERICAL_ZERO) {
			Scale(1./sqrt(norm), opxwave2);
			MultiplyWithOwnTranspose (opxwave2, dm[0], scale);  
			//MultiplyProduct(opxwave2, Transpose(opxwave2), dm[0], scale);
		      } 
		      opxwave2.deallocate();
		    }   
		  }
		  
		  op.deallocate();
		}
      }
    }
  }
  
  void syncaccumulate(int toproc = 0)
  {
    for(int i=1;i<num_threads;++i)
      dm[0] += dm[i];

    distributedaccumulate<SpinAdapted::StackSparseMatrix>(dm[0]);
    synced = true;
  }
};



// accumulates into dm
void StackDensityMatrix::add_onedot_noise(StackWavefunction& wave_solution, StackSpinBlock& big, bool act2siteops)
{

  StackSpinBlock* leftBlock = big.get_leftBlock();
  //p1out << "\t\t\t Modifying density matrix " << endl;


  StackDensityMatrix *dm = 0;
  dm = this;

  onedot_noise_f onedot_noise(dm, wave_solution, big, 1., 1);

  if (leftBlock->has(CRE)) {
    onedot_noise.set_opType(CRE);
    for_all_singlethread(leftBlock->get_op_array(CRE), onedot_noise);
  }
  if (leftBlock->has(DES)) {
    onedot_noise.set_opType(DES);
    for_all_singlethread(leftBlock->get_op_array(DES), onedot_noise);
  }
  //use overlap only when bra and ket are different i.e. when the block has des operator
  if (leftBlock->has(DES)&&leftBlock->has(OVERLAP)) {
    onedot_noise.set_opType(OVERLAP);
    for_all_singlethread(leftBlock->get_op_array(OVERLAP), onedot_noise);
  }
  
  if (dmrginp.hamiltonian() != HUBBARD) {
    
    if (leftBlock->has(CRE_CRE)) {
      onedot_noise.set_opType(CRE_CRE);
      for_all_singlethread(leftBlock->get_op_array(CRE_CRE), onedot_noise);
      
      onedot_noise.set_opType(CRE_DES);
      for_all_singlethread(leftBlock->get_op_array(CRE_DES), onedot_noise);
    } 
    else if (leftBlock->has(DES_DESCOMP)) {
      onedot_noise.set_opType(DES_DESCOMP);
      for_all_singlethread(leftBlock->get_op_array(DES_DESCOMP), onedot_noise);
      
      onedot_noise.set_opType(CRE_DESCOMP);
      for_all_singlethread(leftBlock->get_op_array(CRE_DESCOMP), onedot_noise);
      
    }
    if (leftBlock->has(DES_DES)) {
      onedot_noise.set_opType(DES_DES);
      for_all_singlethread(leftBlock->get_op_array(DES_DES), onedot_noise);
      
      onedot_noise.set_opType(DES_CRE);
      for_all_singlethread(leftBlock->get_op_array(DES_CRE), onedot_noise);
    }
    else if (leftBlock->has(CRE_CRECOMP)) {
      onedot_noise.set_opType(CRE_CRECOMP);
      for_all_singlethread(leftBlock->get_op_array(CRE_CRECOMP), onedot_noise);
      
      onedot_noise.set_opType(DES_CRECOMP);
      for_all_singlethread(leftBlock->get_op_array(DES_CRECOMP), onedot_noise);
      
    }
    
    
  }
  onedot_noise.syncaccumulate();
  
}


}


  
