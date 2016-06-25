#ifndef MPS_NEVPT_HEADER
#define MPS_NEVPT_HEADER

namespace SpinAdapted{
  namespace mps_nevpt{
    extern std::vector<SpinBlock> siteBlocks_noDES;
    extern vector<double> ZeroEnergy;
    extern int sweepIters ;
    void readZeroEnergy();
    void mps_nevpt(double sweep_tol);
  }
}

#endif
