#ifndef SPIN_DIIS_HEADER
#define SPIN_DIIS_HEADER
#include "spinblock.h"
#include "sweep_params.h"

namespace SpinAdapted{
  class DIIS
  {
   public:
    int keepStates;
    int currentIndex;
    ColumnVector b;
    ColumnVector x;
    Matrix B;
    bool buildup;

    void BlockAndDecimate (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys);
    double do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize);
    void updateWavefunction(SweepParams &sweepParams);
    void updateHamiltonian(SweepParams &sweepParams, bool forward);
    double updateErrors(SweepParams &sweepParams, bool forward);
    double errorSweep (SweepParams &sweepParams, SpinBlock& system, SpinBlock& newSystem, const bool &useSlater, const bool& dot_with_sys);
    void updateB(SweepParams &sweepParams);
    void initialize(SweepParams &sweepParams);
    void LoadDIISWavefunction(Wavefunction& sigma, int currentIndex, int iter);
    void LoadDIISError(Wavefunction& sigma, int currentIndex, int iter);
    void SaveDIISWavefunction(const Wavefunction& sigma, int currentIndex, int iter);
    void SaveDIISError(const Wavefunction& sigma, int currentIndex, int iter);
    void maxOverlapPrevRotation(std::vector<Matrix>& rotation, const std::vector<int>& sites, std::vector<Matrix>& prevRotation);
    void LoadInterpolatedWavefunction(Wavefunction& sigma, int iter);
    void SaveInterpolatedWavefunction(const Wavefunction& sigma, int iter);
    void transformBlock(SpinBlock& big, const std::vector<Wavefunction>& wave, SweepParams &sweepParams);

    double assign_dm(std::vector<Matrix>& rotatematrix, std::vector<DiagonalMatrix>& eigenmatrix, 
		     SparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, 
		     vector<vector<int> >& wtsbyquanta, int totalstatesbydm, 
		     int totalstatesbyquanta, std::vector<int>& statesperquanta);

  };
}
#endif

