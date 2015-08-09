#include "wrapper.h"
#ifndef SERIAL
#include "mpi.h"
#endif
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
  ReadInputFromC(argv[1], 3);

  int mpsstate=0;
  
  //dmrginp.add_noninteracting_orbs() = false;
  readMPSFromDiskAndInitializeStaticVariables();
  //initializeGlobalMPS(mpsstate);


  double overlap, hvalue;

  unsigned long temp=1, occ=0;
  int msgsize=1000;
  char msgctr[msgsize];

  ifstream file("determinants");
  file.getline(msgctr, msgsize);
  string s( msgctr);
  vector<string> tok;
  boost::split(tok, s, is_any_of(", \t"), token_compress_on);
  int outMPS = atoi(tok[0].c_str());


  std::vector<MPS> dets;
  std::vector<bool> occupation;
  for (int i=0; i<1; i++) {
    file.getline(msgctr, msgsize);
    //string occstring = "1 1 0 0 1 1 0 0 1 1 0 0 0 0 1 1";

    string ss(msgctr);
    stringstream stream(ss);
    int n;
    while (stream >> n) {
      if (n==1) 
	occupation.push_back(1);
      else
	occupation.push_back(0);
    }
    dets.push_back(MPS(occupation));
  } 


  dets[0].writeToDiskForDMRG(outMPS, true);

#ifndef SERIAL
  MPI_Finalize();
#endif
  return 0;
}
