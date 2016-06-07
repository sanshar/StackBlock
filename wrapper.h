#ifndef WRAPPER_HEADER
#define WRAPPER_HEADER

#include <string>

#ifdef __cplusplus
extern "C" {
#endif

  void initBoostMPI(int argc, char* argv[]) ;
  void ReadInputFromC(char* conf, int outputlevel);
  void readMPSFromDiskAndInitializeStaticVariables(bool initializeDotBlocks=true);
  //void evaluateOverlapAndHamiltonian(unsigned long *occ, int length, double* o, double* h);
  void evaluateOverlapAndHamiltonian(int state1, int state2, double* o, double* h);
  void intFromString(unsigned long &occ, char* s);
  void writeToDisk(unsigned long &occ, int length, int stateIndex);
  void test(char* infile);
  void RDM(char* infile);
  void initializeGlobalMPS(int mpsstate);
  void CollapseToDeterminant(char* s, int stateIndex);
  void seedRandom(int seed);
  void ApplyCD(int i, int j);
  void AddMPSs(int* states, double* scale, int nstates, int outstate);
#ifdef __cplusplus
}
#endif


#endif
