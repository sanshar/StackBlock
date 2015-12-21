#include "global.h"
#include "wrapper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include <iostream>
#include "fciqmchelper.h"

using namespace std;


int main(int argc, char* argv []) {

  int rank=0, size=1;
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Group orig_group, calc_group;
  MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
  MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &Calc);
  calc = boost::mpi::communicator(Calc, boost::mpi::comm_attach);

#endif
  initBoostMPI(argc, argv);
  ReadInputFromC(argv[1], 0);
  dmrginp.matmultFlops.resize(numthrds, 0.);

  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();

  if(!dmrginp.spinAdapted())
    MPS::sweepIters = dmrginp.last_site()/2-2;
  else
    MPS::sweepIters = dmrginp.last_site()-2;
  MPS::spinAdapted = false;


  int instate, outstate;
  std::vector<int> orbitals;
  int newN, newS, newIrrep;
  int nindices;
  std::map<int, string> ops;
  if (rank ==0) {
    printf("Reading file %s\n", argv[2]);
    ifstream file(argv[2]);
    file >> instate; //the state to read
    file >> outstate;

    file >> nindices;
    orbitals.resize(nindices);
    for (int i=0; i<nindices; i++) {
      file >> orbitals[i];
      if (orbitals[i] > 0) {
	if (ops.find(abs(orbitals[i])-1) == ops.end())
	  ops.insert(std::pair<int, string>(abs(orbitals[i])-1, "c"));
	else if (ops[abs(orbitals[i])-1] == "c")
	  ops[abs(orbitals[i])-1] = "cc";
	else 
	  ops[abs(orbitals[i])-1] = "cd";
      }
      else {
	if (ops.find(abs(orbitals[i])-1) == ops.end())
	  ops.insert(std::pair<int, string>(abs(orbitals[i])-1, "d"));
	else if (ops[abs(orbitals[i])-1] == "c")
	  ops[abs(orbitals[i])-1] = "cd";
	else 
	  ops[abs(orbitals[i])-1] = "dd";
      }

    }
    file.close();
  }

  MPS mps(instate);
  mps.writeToDiskForDMRG(outstate);mps.clearTensors();
  mps.initializeStateInfo(true, outstate);
  //mps.initializeStateInfo(false, outstate);
  mps.rightCanonicalize(outstate);


  SpinQuantum q = dmrginp.effective_molecule_quantum_vec()[0];
  q = (q-getSpinQuantum(3))[0];
  
  mps.Operator(ops, outstate, pow(2., 0.5));


#ifndef SERIAL
  MPI_Comm_free(&Calc);
  MPI_Finalize();
#endif
  return 0;
}
