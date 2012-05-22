#include "guess_wavefunction.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif
namespace SpinAdapted{
void GuessWave::TransformLeftBlock(Wavefunction& oldwavefunction, const StateInfo& newstateinfo, const std::vector<Matrix>& RotationMatrix, Wavefunction& tempoldWave)
{  
  for (int a=0; a<tempoldWave.nrows(); a++)
    for (int b=0; b<tempoldWave.ncols(); b++)
    {      
      int olda = newstateinfo.leftStateInfo->leftStateInfo->newQuantaMap[a];
      Matrix& lM = oldwavefunction(olda, b);
      if(lM.Ncols() != 0)
      {
	Matrix tM = RotationMatrix[olda];
	Matrix& nM = tempoldWave.operator_element(a, b);
	MatrixMultiply(tM, 'n', lM, 'n', nM, 1.0);
      }
    }
}

void GuessWave::TransformRightBlock(const Wavefunction& tempnewWave, const StateInfo& oldStateInfo, const std::vector<Matrix>& RotationMatrix, Wavefunction& trial)
{
  for (int a=0; a<tempnewWave.nrows(); a++)
    for (int b=0; b<tempnewWave.ncols(); b++)
    {      
      int transB = oldStateInfo.rightStateInfo->leftStateInfo->newQuantaMap[b];
      Matrix& nM = trial.operator_element(a, transB);
      const Matrix& oM = tempnewWave.operator_element(a, b);
      if(oM.Ncols() != 0)
      {
	Matrix rM = RotationMatrix[transB];
	rM = rM.t();
	MatrixMultiply(oM, 'n', rM, 'n', nM, 1.0);
      }
    }

}




void GuessWave::transpose_previous_wavefunction(Wavefunction& trial, const SpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave)
{
  StateInfo oldStateInfo;
  Wavefunction oldWave;

  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t Transposing previous wavefunction" << endl;
  if(!onedot)
  {
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_rightBlock()->get_sites(), state);
    for (int i = 0; i < trial.nrows(); ++i)
      for (int j = 0; j < trial.ncols(); ++j)
        if (trial.allowed(i, j))
        {
          assert(oldWave.allowed(j, i));
          assert(trial(i, j).Nrows() == oldWave(j, i).Ncols());
          assert(trial(i, j).Ncols() == oldWave(j, i).Nrows());
          const StateInfo* s = &big.get_stateInfo();
	  Matrix tmp = oldWave(j, i);
	  tmp = tmp.t(); // this is really a transpose, not a hermitian conjugate...
	  trial(i, j) = tmp;
	  int parity = 1;
          if (IsFermion(s->leftStateInfo->quanta[i]) && IsFermion(s->rightStateInfo->quanta[j]))
	    parity = -1;
	  int A, B, J;
	  A = s->leftStateInfo->quanta[i].get_s();
	  B = s->rightStateInfo->quanta[j].get_s();
	  J = oldWave.get_deltaQuantum().get_s();
	  parity *= pow(-1.0, static_cast<int>( (3*A - B + J)/2));
	  if (parity == -1)
	    trial(i,j) *= -1.0;
        }
  }
  else
  {
    vector<int> wfsites = big.get_rightBlock()->get_sites();
    wfsites.insert(wfsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
    sort(wfsites.begin(), wfsites.end());
    oldWave.LoadWavefunctionInfo(oldStateInfo, wfsites, state);
    onedot_transpose_wavefunction(oldStateInfo, big.get_stateInfo(), oldWave, trial);
  }

}

/*!                                                                                                                                       
  @brief Transpose wavefunction from [s.][e] to [e.][s] form.                                                                             
*/

void GuessWave::onedot_transpose_wavefunction(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
					      const Wavefunction& guesswf, Wavefunction& transposewf)
{
  ObjectMatrix3D< vector<Matrix> > threewave;
  // first convert into three index wavefunction [s.][e] -> [s].[e]                                                                       
  onedot_twoindex_to_threeindex_wavefunction(guessstateinfo, guesswf, threewave);
  ObjectMatrix3D< vector<Matrix> > threewavetranspose(threewave.NDim2(), threewave.NDim1(), threewave.NDim0());

  // now perform transpose. [s].[e] -> [e].[s], so we need to move state |s>|.> past state |e>                                            
  // keeping track of the parity change                                                                                                   
  for (int a = 0; a < threewave.NDim0(); ++a)
    for (int b = 0; b < threewave.NDim1(); ++b)
      for (int c = 0; c < threewave.NDim2(); ++c)
	if (threewave(a, b, c).size () != 0)
	  {
	    threewavetranspose(c, b, a).resize( threewave(a, b, c).size());
	    int includedindex = 0;
	    for (int i=0; i< threewave(a,b,c).size(); i++) 
	    {
	      if (threewave(a,b,c)[i].Ncols() == 0)
		continue;
	      int ab = guessstateinfo.leftStateInfo->unCollectedStateInfo->quantaMap(a , b)[includedindex];
	      includedindex++;
	      copy(threewave(a, b, c)[i].t(), threewavetranspose(c, b, a)[i]);

	      // from |s.e> configuration, determine parity for |e.s> form                                                              

	      bool aodd = IsFermion(guessstateinfo.leftStateInfo->leftStateInfo->quanta[a]);
	      bool bodd = IsFermion(guessstateinfo.leftStateInfo->rightStateInfo->quanta[b]);
	      bool codd = IsFermion(guessstateinfo.rightStateInfo->quanta[c]);
	      int j1 = guessstateinfo.leftStateInfo->leftStateInfo->quanta[a].get_s();
	      int j2 = guessstateinfo.leftStateInfo->rightStateInfo->quanta[b].get_s();
	      int j3 = guessstateinfo.rightStateInfo->quanta[c].get_s();
	      int J = guessstateinfo.leftStateInfo->unCollectedStateInfo->quanta[ab].get_s();
	      
	      // first, from |s.e> -> |.se>                                                                                             
	      int parity = (aodd && bodd) ? -1 : 1;
	      parity = parity*pow(-1.0, static_cast<int>((3*j1 - j2 + J)/2));
	      // next, |.se> -> |e.s>                                                                                                   
	      bool xorabodd = ((aodd || bodd) && !(aodd && bodd));
	      parity *= (xorabodd && codd) ? -1: 1;
	      parity = parity*pow(-1.0, static_cast<int>((3*J - j3 + transposewf.get_deltaQuantum().get_s())/2));
	      if (parity == -1) 
		threewavetranspose(c, b, a)[i] *= -1.;
	    }
	  }

  onedot_threeindex_to_twoindex_wavefunction(transposestateinfo, threewavetranspose,
					     transposewf, *guessstateinfo.leftStateInfo->unCollectedStateInfo);
}



/*!                                                                                                                                         
@brief Convert wavefunction in three block form (with single dot) [s].[e] to two block form  [s.] [e].
  Note it always group the dot with [s], not with [e].
  Three block form means that the total wavefunction is written explicitly as
  |psi> = c_{s,dot,e} |s>|dot>|e>
  Two block form means that the total wavefunction is written as
  |psi> = c_{sdot,e} |sdot>|e>                                  
  @param twostateinfo: StateInfo for output form [s.][e]                                                                                   
  @param threewavefunction: input wavefunction [s].[e] notation
  @param twowavefunction: output wavefunction in [s.][e] notation                                                                         
*/
void GuessWave::onedot_threeindex_to_twoindex_wavefunction(const StateInfo& twostateinfo,
							   const ObjectMatrix3D< std::vector<Matrix> >& threewavefunction, 
							   Wavefunction& twowavefunction, const StateInfo& prevUnCollectedSI)
{
  int aSz = threewavefunction.NDim0();
  int bSz = threewavefunction.NDim1();
  int cSz = threewavefunction.NDim2();

  const StateInfo& uncollectedstateinfo = *(twostateinfo.leftStateInfo->unCollectedStateInfo);

  twowavefunction.AllowQuantaFor(uncollectedstateinfo, *twostateinfo.rightStateInfo,
                                 twowavefunction.get_deltaQuantum()); // twowavefunction should already know its own quantum number

  for (int a = 0; a < aSz; ++a)
    for (int b = 0; b < bSz; ++b)
      for (int c = 0; c < cSz; ++c)
	if (threewavefunction (a, b, c).size() != 0 && uncollectedstateinfo.allowedQuanta(a, b))
	  for (int i =0; i<uncollectedstateinfo.quantaMap(a,b).size(); i++)
	  {
	    int ab = uncollectedstateinfo.quantaMap (a, b)[i];

	    if (!twowavefunction.allowed(ab,c)) continue;
	    int A, B, AB, C, J, CB;
	    A = twostateinfo.leftStateInfo->leftStateInfo->quanta[a].get_s(); 
	    B = twostateinfo.leftStateInfo->rightStateInfo->quanta[b].get_s();
	    AB = uncollectedstateinfo.quanta[ab].get_s();
	    C =  twostateinfo.rightStateInfo->quanta[c].get_s(); 
	    J = twowavefunction.get_deltaQuantum().get_s(); 

	    int Al, Bl, ABl, Cl, Jl, CBl;
	    Al = twostateinfo.leftStateInfo->leftStateInfo->quanta[a].get_symm().getirrep(); 
	    Bl = twostateinfo.leftStateInfo->rightStateInfo->quanta[b].get_symm().getirrep();
	    ABl = uncollectedstateinfo.quanta[ab].get_symm().getirrep();
	    Cl =  twostateinfo.rightStateInfo->quanta[c].get_symm().getirrep(); 
	    Jl = twowavefunction.get_deltaQuantum().get_symm().getirrep(); 

	    for (int j=0; j<prevUnCollectedSI.quantaMap(c, b).size(); j++) {
	      int cb = prevUnCollectedSI.quantaMap(c,b)[j];
	      CB = prevUnCollectedSI.quanta[cb].get_s();
	      CBl = prevUnCollectedSI.quanta[cb].get_symm().getirrep();
	      int insertionNum = prevUnCollectedSI.quanta[cb].insertionNum(twostateinfo.leftStateInfo->rightStateInfo->quanta[b], twostateinfo.rightStateInfo->quanta[c]);
	      double scale = sixj(A, B, AB, C, J, CB);
	      scale *= pow((1.0*AB+1.0)*(1.0*CB+1.0), 0.5)* pow(-1.0, static_cast<int>((A +B +J +C)/2)); 
	      scale *= Symmetry::spatial_sixj(Al, Bl, ABl, Cl, Jl, CBl);

	      if (threewavefunction(a, b, c)[insertionNum].Nrows() != 0) 
		MatrixScaleAdd(scale, threewavefunction (a, b, c)[insertionNum], twowavefunction.operator_element (ab, c));
	    }
	  }


  // from tensor product form |a>|b> group together blocks with same quantum numbers                                                      
  twowavefunction.CollectQuantaAlongRows(uncollectedstateinfo, *twostateinfo.rightStateInfo);
}

void GuessWave::basic_guess_wavefunction(DiagonalMatrix& e, Wavefunction& trial, const StateInfo *stateinfo, const int state)
{
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t No trial vector" << endl;
  multimap<double, int> e_sort;
  for (int i = 0; i < e.Nrows (); ++i)
    e_sort.insert (pair<double, int> (e (i+1), i+1));

  multimap<double, int>::iterator e_iter = e_sort.begin ();

  int states = stateinfo->totalStates;
  RowVector trialvector(states);
  trialvector = 0.;
  for(int i=0;i<state;++i)
    ++e_iter;
  trialvector(e_iter->second) = 1.;
  trial.CollectFrom(trialvector);
}


void GuessWave::guess_wavefunctions(Wavefunction& solution, DiagonalMatrix& e, const SpinBlock &big,
				    const guessWaveTypes &guesswavetype, const bool &onedot, const int &state, const bool& transpose_guess_wave,
				    double additional_noise)
{
#ifndef SERIAL
  mpi::communicator world;
#endif
  solution.initialise(dmrginp.effective_molecule_quantum(), &big, onedot);
  
  if (!mpigetrank())
  {
    switch(guesswavetype)
    {
    case TRANSFORM:
      transform_previous_wavefunction(solution, big, state, onedot, transpose_guess_wave);
      break;
    case BASIC:
      basic_guess_wavefunction(e, solution, &big.get_stateInfo(), state);
      break;
    case TRANSPOSE:
      transpose_previous_wavefunction(solution, big, state, onedot, transpose_guess_wave);
      break;
    }
    
    double norm = DotProduct(solution, solution);
    if (additional_noise >1e-14) {
      Wavefunction noiseMatrix = solution;
      noiseMatrix.Randomise();
      double overlap = DotProduct(noiseMatrix, solution);
      ScaleAdd(-overlap, solution, noiseMatrix);
      double norm = DotProduct(noiseMatrix, noiseMatrix);
      if (abs(norm) >= 1e-14) {
	if (dmrginp.outputlevel() != 0) 
	  pout << "\t\t\t Norm is "<<norm<<". Adding noise of "<<additional_noise<<" to wavefunction "<<endl;
	ScaleAdd(additional_noise/sqrt(norm), noiseMatrix, solution);
      }
      Normalise(solution);
    }
    Normalise(solution);
    norm = DotProduct(solution, solution);
    if (dmrginp.outputlevel() != 0) 
      pout << "\t\t\t Norm of wavefunction :: "<<norm<<endl;
    
  }
#ifndef SERIAL
  broadcast(world, solution, 0);
#endif
}



void GuessWave::guess_wavefunctions(std::vector<Wavefunction>& solution, DiagonalMatrix& e, const SpinBlock &big, 
				    const guessWaveTypes &guesswavetype, const bool &onedot, const bool &transpose_guess_wave, double additional_noise)
{
  const int nroots = solution.size();

  for(int i=0;i<nroots;++i) 
    guess_wavefunctions(solution[i], e, big, guesswavetype, onedot, i, transpose_guess_wave, additional_noise);

}

/*!  
  @brief Convert wavefunction in twoblock form (with single dot) [s.] [e] to three block form [s].[e] 

  Three block form means that the total wavefunction is written explicitly as
  |psi> = c_{s,dot,e} |s>|dot>|e>

  Two block form means that the total wavefunction is written as
  |psi> = c_{sdot,e} |sdot>|e>
  
  @param stateinfo: StateInfo for [s.][e] superblock
  @param twowavefunction: input wavefunction in [s.][e] notation
  @param threewavefunction: output wavefunction [s].[e] notation
 */
void GuessWave::onedot_twoindex_to_threeindex_wavefunction(const StateInfo& stateinfo, const Wavefunction& twowavefunction, 
							   ObjectMatrix3D< vector<Matrix> >& threewavefunction)
{
  int NDimsys = stateinfo.leftStateInfo->leftStateInfo->quanta.size();
  int NDimdot = stateinfo.leftStateInfo->rightStateInfo->quanta.size();
  int NDimenv = stateinfo.rightStateInfo->quanta.size();

  threewavefunction.ReSize(NDimsys, NDimdot, NDimenv);

  Wavefunction uncollectedwf = twowavefunction;
  uncollectedwf.UnCollectQuantaAlongRows(*stateinfo.leftStateInfo, *stateinfo.rightStateInfo);

  const StateInfo& uncollectedstateinfo = *(stateinfo.leftStateInfo->unCollectedStateInfo);
  for (int ab = 0; ab < uncollectedwf.nrows(); ++ab)
    for (int c = 0; c < uncollectedwf.ncols(); ++c)
      if (uncollectedwf.allowed(ab, c))
	{
	  int a = uncollectedstateinfo.leftUnMapQuanta[ab];
	  int b = uncollectedstateinfo.rightUnMapQuanta[ab];
	  int AB = uncollectedstateinfo.quanta[ab].get_s();
	  int A = stateinfo.leftStateInfo->leftStateInfo->quanta[a].get_s();
	  int B = stateinfo.leftStateInfo->rightStateInfo->quanta[b].get_s();
	  int insertionNum = uncollectedstateinfo.quanta[ab].insertionNum(stateinfo.leftStateInfo->leftStateInfo->quanta[a],  stateinfo.leftStateInfo->rightStateInfo->quanta[b]);
	  vector<SpinQuantum> spq = stateinfo.leftStateInfo->leftStateInfo->quanta[a]+  stateinfo.leftStateInfo->rightStateInfo->quanta[b];
	  if (threewavefunction(a, b, c).size() != spq.size())
	    threewavefunction(a, b, c).resize(spq.size());

	  copy(uncollectedwf(ab, c), threewavefunction(a, b, c)[insertionNum]);
	}
}


void GuessWave::onedot_shufflesysdot(const StateInfo& guessstateinfo, const StateInfo& transposestateinfo,
				     const Wavefunction& guesswf, Wavefunction& transposewf)
{
  ObjectMatrix3D< vector<Matrix> > threewave;
  // first convert into three index wavefunction [s][e.] -> [s].[e]                                                         
  onedot_twoindex_to_threeindex_shufflesysdot(guessstateinfo, guesswf, threewave);
  onedot_threeindex_to_twoindex_wavefunction(transposestateinfo, threewave, transposewf, *guessstateinfo.rightStateInfo->unCollectedStateInfo);
}

void GuessWave::onedot_twoindex_to_threeindex_shufflesysdot(const StateInfo& stateinfo, const Wavefunction& twowavefunction, 
							    ObjectMatrix3D< vector<Matrix> >& threewavefunction)
{
  int NDimsys = stateinfo.leftStateInfo->quanta.size();
  int NDimdot = stateinfo.rightStateInfo->rightStateInfo->quanta.size();
  int NDimenv = stateinfo.rightStateInfo->leftStateInfo->quanta.size();

  threewavefunction.ReSize(NDimsys, NDimdot, NDimenv);

  Wavefunction uncollectedwf = twowavefunction;
  uncollectedwf.UnCollectQuantaAlongColumns(*stateinfo.leftStateInfo, *stateinfo.rightStateInfo);
  const StateInfo& uncollectedstateinfo = *(stateinfo.rightStateInfo->unCollectedStateInfo);

  for (int a = 0; a < uncollectedwf.nrows(); ++a)
    for (int bc = 0; bc < uncollectedwf.ncols(); ++bc)
      if (uncollectedwf.allowed(a, bc))
        {
          int b = uncollectedstateinfo.rightUnMapQuanta[bc];
          int c = uncollectedstateinfo.leftUnMapQuanta[bc];
          bool bodd = IsFermion(stateinfo.rightStateInfo->rightStateInfo->quanta[b]);
          bool codd = IsFermion(stateinfo.rightStateInfo->leftStateInfo->quanta[c]);
	  int parity = (bodd&codd) ? -1 : 1;
	  int j1 = stateinfo.rightStateInfo->leftStateInfo->quanta[c].get_s();
	  int j2 = stateinfo.rightStateInfo->rightStateInfo->quanta[b].get_s();
	  int J = uncollectedstateinfo.quanta[bc].get_s();

	  parity *= pow(-1.0, static_cast<int>( (3*j1 - j2 + J)/2));
	  int insertionNum = uncollectedstateinfo.quanta[bc].insertionNum(stateinfo.rightStateInfo->rightStateInfo->quanta[b], stateinfo.rightStateInfo->leftStateInfo->quanta[c]);
	  vector<SpinQuantum> spq = stateinfo.rightStateInfo->rightStateInfo->quanta[b]+ stateinfo.rightStateInfo->leftStateInfo->quanta[c];
	  if(threewavefunction(a, b, c).size() != spq.size())
	    threewavefunction(a, b, c).resize(spq.size());
          copy(uncollectedwf(a, bc), threewavefunction(a, b, c)[insertionNum]);
          if (parity == -1)
	    threewavefunction(a,b,c)[insertionNum] *= -1.0;
        }
}


void GuessWave::transform_previous_wavefunction(Wavefunction& trial, const SpinBlock &big, const int state, const bool &onedot, const bool& transpose_guess_wave)
{
  if (dmrginp.outputlevel() != 0) 
    pout << "\t\t\t Transforming previous wavefunction " << endl;
  ObjectMatrix3D< vector<Matrix> > oldTrialWavefunction;
  ObjectMatrix3D< vector<Matrix> > newTrialWavefunction;
  StateInfo oldStateInfo;
  Wavefunction oldWave;
  DiagonalMatrix D;
  Matrix U;
  Matrix V;
  std::vector<Matrix> inverseLeftRotationMatrix;
  if (transpose_guess_wave || !onedot){
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_leftBlock()->get_sites(), state);
    LoadRotationMatrix (big.get_leftBlock()->get_leftBlock()->get_sites(), inverseLeftRotationMatrix, state);
  }
  else{
    oldWave.LoadWavefunctionInfo (oldStateInfo, big.get_leftBlock()->get_sites(), state);
    LoadRotationMatrix (big.get_leftBlock()->get_sites(), inverseLeftRotationMatrix, state);
  }

  for (int q = 0; q < inverseLeftRotationMatrix.size (); ++q)
  {
    if (inverseLeftRotationMatrix [q].Nrows () > 0)
    {
      try
      {
        svd(inverseLeftRotationMatrix[q], D, U, V);
      }
      catch (Exception)
      {
        pout << Exception::what() << endl;
        pout << D << endl;
        pout << U << endl;
        pout << V << endl;
        abort();
      }
      Matrix vd = V;
      vd *= D.i();
      inverseLeftRotationMatrix[q].ReSize(V.Nrows(), U.Nrows());
      SpinAdapted::Clear(inverseLeftRotationMatrix[q]);
      MatrixMultiply(vd, 'n', U, 't', inverseLeftRotationMatrix[q], 1.);
    }
  }

  std::vector<Matrix> rightRotationMatrix;
  if(!onedot)
  {
    Wavefunction tempoldWave;
    tempoldWave.AllowQuantaFor(*big.get_stateInfo().leftStateInfo->leftStateInfo, *oldStateInfo.rightStateInfo, oldWave.get_deltaQuantum()); 
    TransformLeftBlock(oldWave, big.get_stateInfo(), inverseLeftRotationMatrix, tempoldWave);

    StateInfo tempoldStateInfo;
    TensorProduct (*(big.get_stateInfo().leftStateInfo->leftStateInfo), *oldStateInfo.rightStateInfo, tempoldStateInfo,
		   PARTICLE_SPIN_NUMBER_CONSTRAINT);
    tempoldStateInfo.CollectQuanta();

    Wavefunction tempnewWave;
    tempnewWave.AllowQuantaFor(*big.get_stateInfo().leftStateInfo, *oldStateInfo.rightStateInfo->leftStateInfo, oldWave.get_deltaQuantum()); 
    StateInfo tempnewStateInfo;
    TensorProduct (*(big.get_stateInfo().leftStateInfo), *oldStateInfo.rightStateInfo->leftStateInfo, tempnewStateInfo,
		   PARTICLE_SPIN_NUMBER_CONSTRAINT);
    tempnewStateInfo.CollectQuanta();
    onedot_shufflesysdot(tempoldStateInfo, tempnewStateInfo, tempoldWave, tempnewWave);
    
    LoadRotationMatrix (big.get_rightBlock()->get_sites(), rightRotationMatrix, state);
    trial.AllowQuantaFor(*big.get_stateInfo().leftStateInfo, *big.get_stateInfo().rightStateInfo, oldWave.get_deltaQuantum()); 
    TransformRightBlock(tempnewWave, oldStateInfo, rightRotationMatrix, trial);
    // from tensor product form |a>|b> group together blocks with same quantum numbers
    //trial.CollectQuantaAlongColumns(*big.get_stateInfo().leftStateInfo, *big.get_stateInfo().rightStateInfo->unCollectedStateInfo);
  }
  else
  {
    vector<int> rotsites;
    if (transpose_guess_wave) {
      rotsites = big.get_rightBlock()->get_sites();
      rotsites.insert(rotsites.end(), big.get_leftBlock()->get_rightBlock()->get_sites().begin(), big.get_leftBlock()->get_rightBlock()->get_sites().end());
      sort(rotsites.begin(), rotsites.end());
      LoadRotationMatrix(rotsites, rightRotationMatrix, state);
    }
    else {
      rotsites = big.get_rightBlock()->get_sites();
      LoadRotationMatrix(rotsites, rightRotationMatrix, state);
    }
    onedot_transform_wavefunction(oldStateInfo, big.get_stateInfo(), oldWave, inverseLeftRotationMatrix, rightRotationMatrix, trial, transpose_guess_wave);
  }

  double norm = DotProduct(trial, trial);
}


/*!  
  @brief Transform wavefunction from previous block configuration 
  [s.][e] -> [s'.][e'] (where s' = R[s.]; [R denotes renormalise] i.e. next block configuration)
  
  @param oldstateinfo: StateInfo for input configuration [s.][e] 
  @param newstateinfo: StateInfo for output configuration [s'.][e'] 
  @param oldWavefunction: input wavefunction
  @param inverseLeftRotationMatrix: inverse of transformation  from [s.]->[s']
  @param rightRotationMatrix: transformation from [.e'] -> [e]
  @param newwavefunction: output wavefunction. Should already have quantum number assigned.
 */

void GuessWave::onedot_transform_wavefunction(const StateInfo& oldstateinfo, const StateInfo& newstateinfo, 
					      const Wavefunction& oldwavefunction, const vector<Matrix>& inverseLeftRotationMatrix, 
					      const vector<Matrix>& rightRotationMatrix, Wavefunction& newwavefunction, 
					      const bool& transpose_guess_wave)
{
  assert (oldwavefunction.get_deltaQuantum() == newwavefunction.get_deltaQuantum());

  Timer transform;
  int oldASz = oldstateinfo.leftStateInfo->quanta.size ();
  int oldCSz = oldstateinfo.rightStateInfo->quanta.size ();

  // Old wavefunction is in [s.][e] configuration
  // Then [s.] -> [s'], and [.e'] -> [e]
  // Make wavefunction in [s'][.e'] configuration 
  // Then later transfer dot, [s'][.e'] -> [s'.][e']
  StateInfo newenvstateinfo;
  ObjectMatrix<char> preparedQuanta;
  Wavefunction tmpwavefunction;
  int aSz;
  int cSz;
  if (transpose_guess_wave){
    TensorProduct (*(newstateinfo.rightStateInfo), *(newstateinfo.leftStateInfo->rightStateInfo), newenvstateinfo,
		   NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    newenvstateinfo.CollectQuanta();
    tmpwavefunction.AllowQuantaFor (*(newstateinfo.leftStateInfo->leftStateInfo), newenvstateinfo, newwavefunction.get_deltaQuantum());
    aSz = newstateinfo.leftStateInfo->leftStateInfo->quanta.size (); 
    cSz = newenvstateinfo.quanta.size();
  }
  else {
    tmpwavefunction.AllowQuantaFor (*(newstateinfo.leftStateInfo), *(newstateinfo.rightStateInfo), newwavefunction.get_deltaQuantum());
    aSz = newstateinfo.leftStateInfo->quanta.size (); 
    cSz = newstateinfo.rightStateInfo->quanta.size();
  }

  ObjectMatrix<Matrix> tmp(oldASz, cSz); //cSz
  for (int a = 0; a < oldASz; ++a)
    for (int c = 0; c < oldCSz; ++c) // oldCSz <= cSz
    {
      const Matrix& oM = oldwavefunction.operator_element(a, c);
      if (oM.Ncols () != 0) // this quanta combination is not allowed
      {
	int transC = oldstateinfo.rightStateInfo->newQuantaMap [c];

	Matrix& tM = tmp (a, transC);
	Matrix rM = rightRotationMatrix[transC];
	tM.ReSize (oM.Nrows (), rM.Nrows ());
	
	double parity = 1.;
	
	SpinAdapted::Clear (tM);
	// we have the rotation matrices C_r'r, C_l'l, and we are transforming the wavefunction
	// coefficient d_lr -> consequently, we need to multiply by the transpose, but NOT hermitian
	// conjugate, of rM
	rM = rM.t();
	MatrixMultiply (oM, 'n', rM, 'n', tM, parity); // phi -> c psi : phi d -> psi c c^t d : consequently compute c^t d 
      }
    }

  for (int c = 0; c < cSz; ++c) //cSz
    for (int a = 0; a < aSz; ++a) //aSz // aSz <= oldASz
    {
      int oldA = 0;
      if (transpose_guess_wave)
	oldA = newstateinfo.leftStateInfo->leftStateInfo->newQuantaMap [a];
      else
	oldA = newstateinfo.leftStateInfo->newQuantaMap [a];

      Matrix& tM = tmp(oldA,c); //tmp (oldA, b, c);
      if (tM.Ncols () != 0) // this quanta combination is not allowed
      {
	//assert (newstateinfo.leftStateInfo->leftStateInfo->quanta [a] == oldstateinfo.leftStateInfo->quanta [oldA]);
	const Matrix& lM = inverseLeftRotationMatrix [oldA];
	//assert (newstateinfo.leftStateInfo->leftStateInfo->quantaStates [a] == lM.Nrows ());
	Matrix& nM = tmpwavefunction.operator_element(a, c);//tmpwavefunction (a, b, c);
	nM.ReSize (lM.Nrows (), tM.Ncols ());
	SpinAdapted::Clear (nM);
	MatrixMultiply (lM, 'n', tM, 'n', nM, 1.); 
      }
    }

  if (!transpose_guess_wave){
    newwavefunction = tmpwavefunction;
    newwavefunction.set_onedot(oldwavefunction.get_onedot());
    return;
  }

  // Now, wavefunction is in [s'][e'.] config, change to [s'.][e'] config.
  StateInfo tempoldStateInfo;
  TensorProduct (*(newstateinfo.leftStateInfo->leftStateInfo), newenvstateinfo, tempoldStateInfo,
		 PARTICLE_SPIN_NUMBER_CONSTRAINT);
  tempoldStateInfo.CollectQuanta();

  onedot_shufflesysdot(  tempoldStateInfo, newstateinfo, tmpwavefunction, newwavefunction);
}

}
