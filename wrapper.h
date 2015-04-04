#ifndef WRAPPER_HEADER
#define WRAPPER_HEADER

#include <string>

#ifdef __cplusplus
extern "C" {
#endif

  void initBoostMPI(int argc, char* argv[]) ;
  void ReadInputFromC(char* conf, int outputlevel);
  void readMPSFromDiskAndInitializeStaticVariables(bool initializeDotBlocks=true);
  void evaluateOverlapAndHamiltonian(unsigned long *occ, int length, double* o, double* h);
  void intFromString(unsigned long &occ, const char* s);
  void test(char* infile);
  void initializeGlobalMPS(int mpsstate);
#ifdef __cplusplus
}
#endif


#endif
