#ifndef SERIAL
#include "mpi.h"
#endif
#include "wrapper.h"
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "fciqmchelper.h"
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "global.h"

using namespace std;

int main(int argc, char* argv []) {

  int rank=0, size=1;
#ifndef SERIAL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  initBoostMPI(argc, argv);

  ReadInputFromC(argv[1], -1);

  int mpsstate=0;
  
  //dmrginp.add_noninteracting_orbs() = false;
  readMPSFromDiskAndInitializeStaticVariables(false);
  //initializeGlobalMPS(mpsstate);


  double ovalue, hvalue;

  unsigned long temp=1, occ=0;

  std::vector<int> states;
  ifstream file("states");
  int stateindex ;
  while(file >> stateindex) {
    states.push_back(stateindex);
    if (rank == 0)
      printf("reading state %i\n", stateindex);
  }


  std::vector<MPS> mpsstates;
  for (int i=0; i<states.size(); i++)
    mpsstates.push_back(MPS(states[i]));


  std::vector<double> energies(2*states.size()+1, 0.0);

  //use the n+1 rule first
  for (int n=1; n<states.size(); n++) {
    hvalue = 0; ovalue = 0;
    calcHamiltonianAndOverlap(mpsstates[0], mpsstates[n], hvalue, ovalue);
    if (rank == 0)
      printf(" E(%3i) = %18.10e\n", n+1, hvalue);
    if (rank == 0)
      printf(" <0|H|%i> = %18.10e, <0|%i> = %18.10e\n", n+1, hvalue, n+1, ovalue);
  }
  if (rank == 0)
    printf("*** About to use the 2n+1 rule\n");

  //formula taken from JCP (112), Leininger
  for (int n=0; n<states.size(); n++) {

    //even value, which is to be done only if we dont need E0
    if (n != 0) {
      hvalue = 0; ovalue = 0;
      calcHamiltonianAndOverlap(mpsstates[n-1], mpsstates[n], hvalue, ovalue);
      energies[2*n] = hvalue;
      for (int i=1; i<n+1; i++)
	for (int j=1; j<n; j++) {
	  ovalue = calculateOverlap(mpsstates[i], mpsstates[j]);
	  energies[2*n] -= energies[2*n - i - j]*ovalue;
	}
      if (rank == 0)
	printf(" E(%3i) = %18.10e\n", 2*n, energies[2*n]);
    }

    //odd values
    hvalue = 0; ovalue = 0;
    calcHamiltonianAndOverlap(mpsstates[n], mpsstates[n], hvalue, ovalue);
    energies[2*n+1] = hvalue;
    for (int i=1; i<n+1; i++)
      for (int j=1; j<n+1; j++) {
	if (i==states.size() && j==states.size()) continue;
	ovalue = calculateOverlap(mpsstates[i], mpsstates[j]);
	energies[2*n+1] -= energies[2*n+1 - i - j]*ovalue;
      }

    if (rank == 0)
      printf(" E(%3i) = %18.10e\n", 2*n+1, energies[2*n+1]);
  }



#ifndef SERIAL
  MPI_Finalize();
#endif
  return 0;
}
