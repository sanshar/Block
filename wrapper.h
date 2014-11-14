#ifndef WRAPPER_HEADER
#define WRAPPER_HEADER

#ifdef __cplusplus
extern "C" {
#endif

  void initBoostMPI(int argc, char* argv[]) ;
  void ReadInputFromC(char* conf, int outputlevel);
  void readMPSFromDiskAndInitializeStaticVariables();
  void evaluateOverlapAndHamiltonian(unsigned long *occ, int length, double* o, double* h);
  void test();
  void intFromString(unsigned long &occ, const char* s);
  void initializeGlobalMPS(int mpsstate);
#ifdef __cplusplus
}
#endif


#endif
