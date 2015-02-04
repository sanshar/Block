#include "fciqmchelper.h"
#include "rotationmat.h"
#include "operatorfunctions.h"
#include "spinblock.h"
#include "initblocks.h"
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

namespace SpinAdapted{

  //initializing the static variables
  int MPS::sweepIters ;
  bool MPS::spinAdapted ;
  std::vector<SpinBlock> MPS::siteBlocks;

void QSTensor::remove_empty() {
  for (int i = 0; i < nl; ++i) {
    for (int j = 0; j < nr; ++j) {
      if (allowed(i,j) && (data(i,j).Nrows() * data(i,j).Ncols() == 0)) {
        data(i,j).ReSize(0,0);
        allowedQuanta(i,j) = false;
      }
    }
  }
  int count = 0;
  std::vector<int> ltemp;
  for (int i = 0; i < nl; ++i) {
    bool empty = true;
    for (int j = 0; j < nr; ++j) {
      if (allowed(i,j)) empty = false;
    }
    if (empty) {
      ltemp.push_back(-1);
    } else {
      ltemp.push_back(count);
      ++count;
    }
  }
  int nl_new = count;
  count = 0;
  std::vector<int> rtemp;
  for (int i = 0; i < nr; ++i) {
    bool empty = true;
    for (int j = 0; j < nl; ++j) {
      if (allowed(j,i)) empty = false;
    }
    if (empty) {
      rtemp.push_back(-1);
    } else {
      rtemp.push_back(count);
      ++count;
    }
  }
  int nr_new = count;

  if (nl_new < nl || nr_new < nr) {
    auto allowedQuantaOld = allowedQuanta;
    auto dataOld = data;
    auto leftQuantaOld = leftQuanta;
    auto rightQuantaOld = rightQuanta;
    allowedQuanta.ReSize(nl_new, nr_new);
    data.ReSize(nl_new, nr_new);
    for (int i = 0; i < nl; ++i) {
      for (int j = 0; j < nr; ++j) {
        if (ltemp[i] >= 0 && rtemp[j] >= 0) {
          allowedQuanta(ltemp[i],rtemp[j]) = allowedQuantaOld(i,j);
          data(ltemp[i],rtemp[j]) = dataOld(i,j);
        }
      }
    }
    leftQuanta.clear();
    for (int i = 0; i < nl; ++i) {
     if (ltemp[i] >= 0) {
      leftQuanta.push_back(leftQuantaOld[i]);
     }
    }
    rightQuanta.clear();
    for (int i = 0; i < nr; ++i) {
     if (rtemp[i] >= 0) {
      rightQuanta.push_back(rightQuantaOld[i]);
     }
    }
    nl = nl_new; nr = nr_new;
  }
}

QSTensor TensorProduct(const QSTensor& A, const QSTensor& B) {
  QSTensor C(A.lQuanta(), B.rQuanta());
  for (int i = 0; i < C.lsize(); ++i) {
    for (int j = 0; j < C.rsize(); ++j) {
      for (int k = 0; k < A.rsize(); ++k) {
        if (A.allowed(i,k)) {
          for (int l = 0; l < B.lsize(); ++l) {
            if (B.allowed(l,j) && A.rQuanta()[k] == B.lQuanta()[l]) {
              if (C.allowed(i,j)) {
                C(i,j) += A(i,k) * B(l,j);
              } else {
                C.allowedMatrix()(i,j) = true;
                C(i,j) = A(i,k) * B(l,j);
              }
            }
          }
        }
      }
    }
  }
  C.remove_empty();
  return C;
}

ostream& operator<< (ostream& os, const QSTensor& q) {
  os << "Left SpinQuantums" << endl;
  for (int i = 0; i < q.nl; ++i) {
    os << q.leftQuanta[i] << "  ";
  }
  os << endl;
  os << "Right SpinQuantums" << endl;
  for (int i = 0; i < q.nr; ++i) {
    os << q.rightQuanta[i] << "  ";
  }
  os << endl;
  os << "Matrix Elements" << endl;
  for (int i = 0; i < q.nl; ++i) {
    for (int j = 0; j < q.nr; ++j) {
      if (q.allowed(i,j)) {
        os << i << " " << j << endl;
        os << q.data(i,j);
      } 
    }
  }
  return os;
}

void SaveQSTensor(const int& site, const QSTensor& m, int state) {
  Timer disktimer;
  int rank = mpigetrank();
  if (rank == 0) {
    char file [5000];
    if (state == -1)
      sprintf(file, "%s%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/MPS-", site, ".", mpigetrank(), ".state_average.tmp");
    else
      sprintf(file, "%s%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/MPS-", site, ".", mpigetrank(), ".state", state, ".tmp");
    std::ofstream ofs(file, std::ios::binary);
    boost::archive::binary_oarchive save_mat(ofs);
    save_mat << m;
    ofs.close();
  }
}

void LoadQSTensor(const int& site, QSTensor& m, int state) {
  Timer disktimer;
  int rank = mpigetrank();
  if (rank == 0) {
    char file [5000];
    if (state == -1)
      sprintf(file, "%s%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/MPS-", site, ".", mpigetrank(), ".state_average.tmp");
    else
      sprintf(file, "%s%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/MPS-", site, ".", mpigetrank(), ".state", state, ".tmp");
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> m;
    ifs.close();
  }
}
  void MPS::buildMPSrep() {
    // the first site
    SpinQuantum sq[] = {SpinQuantum(0, SpinSpace(0), IrrepSpace(0))};
    int qs[] = {1};
    StateInfo statel(1, sq, qs);

    for (int i = 0; i < MPS::sweepIters; ++i) {
      //cout << "----------------- Site " << i << " ---------------" << endl;      
      StateInfo statep = MPS::siteBlocks[i].get_stateInfo();
      StateInfo stater;
      TensorProduct(statel, statep, stater, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      stater.CollectQuanta();
      vector<QSTensor> A(statep.quanta.size(), QSTensor(statel.quanta, stater.quanta));
      for (int r_idx = 0; r_idx < stater.oldToNewState.size(); ++r_idx) { // index of qr
        auto OldtoNew = stater.oldToNewState[r_idx];
        int temp = 0;
        for (int k = 0; k < OldtoNew.size(); ++k) {
          int l_idx = stater.leftUnMapQuanta[OldtoNew[k]];
          int p_idx = stater.rightUnMapQuanta[OldtoNew[k]];
          A[p_idx].allowedMatrix()(l_idx, r_idx) = true;
          A[p_idx](l_idx, r_idx).ReSize(statel.quantaStates[l_idx], stater.quantaStates[r_idx]);
          A[p_idx](l_idx, r_idx) = 0.;
          IdentityMatrix I(statel.quantaStates[l_idx]);
          int temp_new = temp + statel.quantaStates[l_idx];
          A[p_idx](l_idx, r_idx).Columns(temp+1, temp_new) = I;
          temp = temp_new;
        }
      }
      for (int j = 0; j < A.size(); ++j) {
        A[j].remove_empty();
      }

      auto RotMat = getSiteTensors(i);
      
      QSTensor QSRotMat(stater.quanta, stater.quanta);
      for (int j = 0; j < stater.quanta.size(); ++j) {
        if (RotMat[j].Ncols() > 0) {
          QSRotMat.allowedMatrix()(j,j) = true;
          QSRotMat(j,j) = RotMat[j];
        }
      }
      QSRotMat.remove_empty();
      //cout << "Rotation Matrix" << endl;
      //cout << QSRotMat << endl;
      
      vector<QSTensor> B;
      for (int j = 0; j < A.size(); ++j) {
        B.push_back(TensorProduct(A[j], QSRotMat));
        //cout << "renormalised: Physical Index " << j << endl;
        //cout << B[j] << endl;
      }
      MPSrep.push_back(B);
      StateInfo renorm_stater;
      SpinAdapted::StateInfo::transform_state(RotMat, stater, renorm_stater);
      statel = renorm_stater;
    }
    // Second last site
    //cout << "----------------- Site " << MPS::sweepIters << " ---------------" << endl;
    StateInfo statep = MPS::siteBlocks[MPS::sweepIters].get_stateInfo();
    StateInfo stater;
    TensorProduct(statel, statep, stater, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    stater.CollectQuanta();
    vector<QSTensor> A(statep.quanta.size(), QSTensor(statel.quanta, stater.quanta));
    for (int r_idx = 0; r_idx < stater.oldToNewState.size(); ++r_idx) { // index of qr
      auto OldtoNew = stater.oldToNewState[r_idx];
      int temp = 0;
      for (int k = 0; k < OldtoNew.size(); ++k) {
        int l_idx = stater.leftUnMapQuanta[OldtoNew[k]];
        int p_idx = stater.rightUnMapQuanta[OldtoNew[k]];
        A[p_idx].allowedMatrix()(l_idx, r_idx) = true;
        A[p_idx](l_idx, r_idx).ReSize(statel.quantaStates[l_idx], stater.quantaStates[r_idx]);
        A[p_idx](l_idx, r_idx) = 0.;
        IdentityMatrix I(statel.quantaStates[l_idx]);
        int temp_new = temp + statel.quantaStates[l_idx];
        A[p_idx](l_idx, r_idx).Columns(temp+1, temp_new) = I;
        temp = temp_new;
      }
    }
    for (int j = 0; j < A.size(); ++j) {
      A[j].remove_empty();
    }
    statep = MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo();
    auto wfn = getw();
    QSTensor QSWfn(stater.quanta, statep.quanta);
    for (int i = 0; i < QSWfn.lsize(); ++i) {
      for (int j = 0; j < QSWfn.rsize(); ++j) {
        if (wfn.allowed(i,j)) {
          QSWfn.allowedMatrix()(i,j) = true;
          QSWfn(i,j) = wfn(i,j);
        }
      }
    }
    QSWfn.remove_empty();
    vector<QSTensor> B;
    for (int j = 0; j < A.size(); ++j) {
      B.push_back(TensorProduct(A[j], QSWfn));
      //cout << "renormalised: Physical Index " << j << endl;
      //cout << B[j] << endl;
    }
    MPSrep.push_back(B);
    // last site
    //cout << "----------------- Site " << MPS::sweepIters+1 << " ---------------" << endl;
    vector<SpinQuantum> qt(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
    Matrix m(1,1);
    m = 1;
    B.clear();
    B.resize(statep.quanta.size(), QSTensor(QSWfn.rQuanta(), qt));
    for (int j = 0; j < B.size(); ++j) {
      B[j].allowedMatrix()(j, 0) = true;
      B[j](j,0) = m;
      B[j].remove_empty();
      //cout << "renormalised: Physical Index " << j << endl;
      //cout << B[j] << endl;
    }
    MPSrep.push_back(B);    
  }

  double MPS::get_coefficient(const vector<bool>& occ_strings) {
    if (MPSrep.size() == 0) {
      buildMPSrep();
    }
    assert(occ_strings.size() == MPSrep.size()*2);
    int idx = occ_strings[0]*2 + occ_strings[1];
    QSTensor temp = MPSrep[0][idx];
    for (int i = 1; i < MPSrep.size(); ++i) {
      idx = occ_strings[i*2]*2 + occ_strings[i*2+1];
      temp = TensorProduct(temp, MPSrep[i][idx]);
    }
    if (temp.lsize() == 0 || temp.rsize() == 0 || !temp.allowed(0,0) || temp(0,0).Nrows() == 0 || temp(0,0).Ncols() == 0) {
      return 0.;
    } else {
      return temp(0,0)(1,1);
    }
  }



  //assumes that the state has been canonicalized in the left canonical form already and stored on disk
  MPS::MPS(int stateindex) {

    if (mpigetrank() == 0) {
      std::vector<int> rotSites(2,0);
      if (!dmrginp.spinAdapted()) {
	rotSites[0] = 0; rotSites[1] = 3; //this should be 1 for spinadapted
      }
      else {
	rotSites[0] = 0; rotSites[1] = 1; //this should be 1 for spinadapted
      }
      
      Matrix m(1,1); m=0;
      std::vector<Matrix> rotMat;
      if (!dmrginp.spinAdapted()) {
	rotMat = std::vector<Matrix>(4,m);
	rotMat[0](1,1) = 1.0;rotMat[1](1,1) = 1.0;rotMat[2](1,1) = 1.0;rotMat[3](1,1) = 1.0;
      }
      else {
	if (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0) {
	  rotMat = std::vector<Matrix>(4,m);
	  rotMat[0](1,1) = 1.0;rotMat[1](1,1) = 1.0;rotMat[2](1,1) = 1.0;rotMat[3](1,1) = 1.0;
	}
	else {
	  rotMat = std::vector<Matrix>(3,m);
	  rotMat[0](1,1) = 1.0;rotMat[1](1,1) = 1.0;rotMat[2](1,1) = 1.0;
	}
      }
      
      SiteTensors.push_back(rotMat);
      
      for (int i=0; i<MPS::sweepIters; i++) {
	LoadRotationMatrix(rotSites, rotMat, stateindex);
	
	SiteTensors.push_back(rotMat);
	if (!dmrginp.spinAdapted())
	  rotSites[1] += 2;
	else
	  rotSites[1] += 1;
      }
      
      //loading the final wavefunction
      if (!dmrginp.spinAdapted())
	rotSites[1] -=2;
      else
	rotSites[1] -=1;
      
      StateInfo s;
      w.LoadWavefunctionInfo(s, rotSites, stateindex);
    }

#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world, w, 0);
#endif
  }


  //this is a helper function to make a MPS from a occupation number representation of determinant
  void MPS::Init(std::vector<bool>& occnum)
  {
    //assert(occnum.size() == dmrginp.last_site());
    Matrix m(1,1); m=1.0;
    Matrix dummy;

    //first rotation matrix 
    std::vector<Matrix> rotMat; 
    int index = occnum[0]*2+occnum[1];

    if (!dmrginp.spinAdapted()) {
      rotMat.resize(4, dummy);
      rotMat[index] = m;
    }
    else {
      if (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0) {
	rotMat.resize(4, dummy);
	rotMat[index] = m;
      }
      else {
	if (index != 0) index--;
	rotMat.resize(3, dummy);
	rotMat[index] = m;
      }
    }

    SiteTensors.push_back(rotMat);
    //cout << MPS::siteBlocks[0]<<endl;
    SpinQuantum sTotal = MPS::siteBlocks[0].get_stateInfo().quanta.at(index);
    

    for (int i=0; i<MPS::sweepIters-1; i++) {
      //stateinfo of in incoming bond of dimension 1
      SpinQuantum sq[] = {sTotal}; int qs[] = {1}; int n = 1;
      StateInfo stateTotal(n, sq, qs), currentState;
      TensorProduct(stateTotal, const_cast<StateInfo&>(MPS::siteBlocks[i+1].get_stateInfo()), currentState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      std::vector<Matrix> rotMat; rotMat.resize(currentState.quanta.size(), dummy);

      int spinQuantumIndex = occnum[2*i+2]*2+occnum[2*i+3];

      SpinQuantum sqout;
      if (spinQuantumIndex == 0) sqout = sTotal;
      else if (spinQuantumIndex == 1) {
	assert(sTotal.totalSpin.getirrep() > 0);
	sqout = SpinQuantum(sTotal.particleNumber+1, SpinSpace(sTotal.totalSpin.getirrep()-1), (sTotal.orbitalSymmetry+SymmetryOfOrb(i+1))[0]);
      }
      else if (spinQuantumIndex == 2) {
	sqout = SpinQuantum(sTotal.particleNumber+1, SpinSpace(sTotal.totalSpin.getirrep()+1), (sTotal.orbitalSymmetry+SymmetryOfOrb(i+1))[0]);
      }
      else if (spinQuantumIndex == 3) {
	sqout = SpinQuantum(sTotal.particleNumber+2, sTotal.totalSpin, (((sTotal.orbitalSymmetry+SymmetryOfOrb(i+1))[0])+SymmetryOfOrb(i+1))[0]  );
      }
      pout << currentState<<endl;
      pout << sqout<<endl;
      std::vector<SpinQuantum>::iterator it = find(currentState.quanta.begin(), currentState.quanta.end(), sqout);
      if (it == currentState.quanta.end()) {
	pout << "Something is probably wrong with the determinant string. please check again."<<endl;
	exit(0);
      }
      int index = it - currentState.quanta.begin();
      rotMat[index]=m;

      sTotal = currentState.quanta[index];

      //cout << i+1<<" "<<sTotal<<endl;
      SiteTensors.push_back(rotMat);
    }
    
    //stateinfo of in incoming bond of dimension 1
    SpinQuantum sq[] = {sTotal}; int qs[] = {1}; int n = 1;
    StateInfo stateTotal(n, sq, qs), secondLastState;

    //the incoming bond x k-1 site stateinfo
    TensorProduct(stateTotal, const_cast<StateInfo&>(MPS::siteBlocks[MPS::sweepIters].get_stateInfo()), 
		  secondLastState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    int index1 =  occnum[2*MPS::sweepIters]*2+occnum[2*MPS::sweepIters+1];

    if (dmrginp.spinAdapted() && secondLastState.quanta.size() <= 3 && index1 > 1) index1--;


    //now make wavefunction with the big state A
    w.AllowQuantaFor(secondLastState, MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo(), 
		     dmrginp.effective_molecule_quantum_vec());

    int index2 = occnum[2*MPS::sweepIters+2]*2+occnum[2*MPS::sweepIters+3];
    if (dmrginp.spinAdapted() && index2 > 1) index2--;

    w(index1, index2) = m;
  }


  //this is a helper function to make a MPS from a occupation number representation of determinant
  //this representation is slightly different than the usual occupation, here each integer
  //element is a spatial orbital which can have a value 0, -1, 1, or 2.
  MPS::MPS(ulong *occnum, int length)
  {
    assert(length*64 >= dmrginp.last_site());
    
    //convert the int array into a vector<bool>
    std::vector<bool> occ(dmrginp.last_site(), 0);
    
    ulong temp = 1;
    int index = 0;
    for (int i=0; i <length ; i++) {
      long occtemp = occnum[i];
      for (int j=63; j>=0; j--) {
	if (dmrginp.spinAdapted() && index >=2*dmrginp.last_site()) break;
	if (!dmrginp.spinAdapted() && index >=dmrginp.last_site()) break;
	
	occ[index] = (occnum[i]>>j) & temp  ;
	
	index++;
      }
    }
    
    Init(occ);
  }
  
  
  MPS::MPS(std::vector<bool>& occ) {
    Init(occ);
  }

  //writes an MPS to the disk so that DMRG can use it as an initial guess
  void MPS::writeToDiskForDMRG(int stateindex, bool writeStateAverage)
  {
    
    StateInfo stateOfSites12, leftState, renormState, currentState;
    TensorProduct(const_cast<StateInfo&>(MPS::siteBlocks[0].get_stateInfo()), const_cast<StateInfo&>(MPS::siteBlocks[1].get_stateInfo()), stateOfSites12, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    stateOfSites12.CollectQuanta();

    bool isDeterminant = false;

    Matrix dummy, m(1,1);m=1.0;
    //first check if this MPS was produced by a determinant or read from DMRG output
    //if it is a determinant the first rotation matrix has some of the ncols 0
    for (int i=0; i<SiteTensors[0].size(); i++) {
      if (SiteTensors[0][i].Ncols() == 0) {
	isDeterminant = true;
	break;
      }
    }

    //if it was generated from a determinant we have to make a combined rotation matrix
    //for the first two sites
    if (mpigetrank() == 0 ) {
      Matrix dummy, m(1,1);m=1.0;
      std::vector<Matrix> Rotations12(stateOfSites12.quanta.size(), dummy);
      if (isDeterminant) {
	
	for (int I=0; I<stateOfSites12.quanta.size(); I++)
	{
	  const std::vector<int>& oldToNewI = stateOfSites12.oldToNewState.at(I); //all the uncollected state 
	  int stateIndex = 0;
	  for (int i=0; i<oldToNewI.size(); i++) {
	    int leftq = stateOfSites12.leftUnMapQuanta[ oldToNewI[i]];
	    int rightq = stateOfSites12.rightUnMapQuanta[ oldToNewI[i]]; 
	    
	    if (SiteTensors[1][rightq].Ncols() != 0) { //the right quanta is retained
	      
	      //right q is from a dot and only has single quantas, further we assume that the only options are empty or closed
	      assert (MPS::siteBlocks[1].get_stateInfo().quanta[rightq].get_s().getirrep() == 0);
	      
	      if (Rotations12[I].Ncols() == 0) {
		Rotations12[I] = Matrix(stateOfSites12.quantaStates[I], SiteTensors[0][leftq].Ncols());
		for (int c = 0; c<SiteTensors[0][leftq].Ncols(); c++)
		  Rotations12[I](c+1, c+1+stateIndex) = 1.0;
	      }
	      else {pout << "We cannot have multiple ways to arrive at the same quantum"<<endl;exit(0);}
	    }
	    stateIndex += MPS::siteBlocks[0].get_stateInfo().quantaStates[leftq]*MPS::siteBlocks[1].get_stateInfo().quantaStates[rightq];
	  }
	}
      }	      


      std::vector<int> rotSites(2,0);
      if (!dmrginp.spinAdapted()) {
	rotSites[0] = 0; rotSites[1] = 3; //this should be 1 for spinadapted
      }
      else {
	rotSites[0] = 0; rotSites[1] = 1; //this should be 1 for spinadapted
      }
      
      if (isDeterminant) {
	SaveRotationMatrix(rotSites, Rotations12, stateindex);
	if (writeStateAverage)
	  SaveRotationMatrix(rotSites, Rotations12, -1);
	SpinAdapted::StateInfo::transform_state(Rotations12, stateOfSites12, renormState);
      }
      else {
	SaveRotationMatrix(rotSites, SiteTensors[1], stateindex);
	if (writeStateAverage)
	  SaveRotationMatrix(rotSites, SiteTensors[1], -1);
	SpinAdapted::StateInfo::transform_state(SiteTensors[1], stateOfSites12, renormState);
      }

      if (!dmrginp.spinAdapted())
	rotSites[1] += 2;
      else
	rotSites[1] += 1;
      

      for (int i=1; i<MPS::sweepIters; i++) {
	leftState = renormState;
	TensorProduct(leftState, const_cast<StateInfo&>(MPS::siteBlocks[i+1].get_stateInfo()), currentState, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
	currentState.CollectQuanta();
	
	SaveRotationMatrix(rotSites, SiteTensors[i+1], stateindex);
	if (writeStateAverage)
	  SaveRotationMatrix(rotSites, SiteTensors[i+1], -1);
	
	if (!dmrginp.spinAdapted())
	  rotSites[1] += 2;
	else
	  rotSites[1] += 1;

	SpinAdapted::StateInfo::transform_state(SiteTensors[i+1], currentState, renormState);

      }
      StateInfo bigState;
      TensorProduct(currentState, const_cast<StateInfo&>(MPS::siteBlocks[MPS::sweepIters+1].get_stateInfo()), bigState, PARTICLE_SPIN_NUMBER_CONSTRAINT);

      if (!dmrginp.spinAdapted())
	rotSites[1] -= 2;
      else
	rotSites[1] -= 1;

      w.SaveWavefunctionInfo (bigState, rotSites, stateindex);

    }
    
  }

  double calculateOverlap(const MPS& statea, const MPS& stateb) {

    Overlap siteOverlap;
    SpinBlock system, siteblock;
    bool forward = true, restart=false, warmUp = false;
    int leftState=0, rightState=1, forward_starting_size=1, backward_starting_size=0, restartSize =0;
    InitBlocks::InitStartingBlock(system, forward, leftState, rightState, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0); 

    siteOverlap.makeIdentity(system.get_stateInfo());

    StateInfo stateA=system.get_stateInfo(), stateB=system.get_stateInfo();
    Overlap o = siteOverlap;

    std::vector<Matrix> Rotationa, Rotationb;
    for (int i=0; i<MPS::sweepIters; i++) {

      StateInfo renormA, renormB;
      if (mpigetrank() == 0) {
	Rotationa = statea.getSiteTensors(i);
	Rotationb = stateb.getSiteTensors(i);
      }

#ifndef SERIAL
      mpi::communicator world;
      mpi::broadcast(world, Rotationa, 0);
      mpi::broadcast(world, Rotationb, 0);
#endif

      SpinAdapted::StateInfo::transform_state(Rotationa, stateA, renormA);
      SpinAdapted::StateInfo::transform_state(Rotationb, stateB, renormB);
      o.renormalise_transform(Rotationa, &renormA, Rotationb, &renormB);

      SpinBlock dotsite(i+1, i+1, 0, false);
      //make overlap and state info for the current site
      siteOverlap.makeIdentity(dotsite.get_stateInfo());
      
      //take tensor product of stateinfo
      StateInfo A, B;
      TensorProduct(renormA, const_cast<StateInfo&>(dotsite.get_stateInfo()), 
		    A, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      TensorProduct(renormB, const_cast<StateInfo&>(dotsite.get_stateInfo()), 
		    B, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
      A.CollectQuanta(); B.CollectQuanta();
      stateA = A; stateB = B;

      //build the new Overlap matrix
      Overlap onew;
      onew.set_deltaQuantum(1, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)));
      onew.set_built()=true; onew.allocate(stateA, stateB); 
      onew.set_initialised() = true;
      SpinAdapted::operatorfunctions::TensorProduct(o, siteOverlap, &stateA, &stateB, onew, 1.0, true);
      o = onew;
    }
    Wavefunction temp = statea.getw();
    temp.Clear();

    SpinBlock dotsite(MPS::sweepIters+1, MPS::sweepIters+1, 0, false);
    siteOverlap.makeIdentity(dotsite.get_stateInfo());

    StateInfo braStateInfo, ketStateInfo;
    TensorProduct(stateA, const_cast<StateInfo&>(dotsite.get_stateInfo()), 
		  braStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);
    TensorProduct(stateB, const_cast<StateInfo&>(dotsite.get_stateInfo())
		  , ketStateInfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

    SpinAdapted::operatorfunctions::TensorMultiply(o, siteOverlap, &braStateInfo, &ketStateInfo, stateb.getw(), temp, SpinQuantum(0, SpinSpace(0), IrrepSpace(0)), true, 1.0);
    double overlap = DotProduct(statea.getw(), temp);
    return overlap;
  }


  void calcHamiltonianAndOverlap(const MPS& statea, const MPS& stateb, double& h, double& o, bool sameStates) {

    SpinBlock system, siteblock;
    bool forward = true, restart=false, warmUp = false;
    int leftState=0, rightState=1, forward_starting_size=1, backward_starting_size=0, restartSize =0;
    int statebindex = sameStates ? 0 : 1;
    InitBlocks::InitStartingBlock(system, forward, leftState, statebindex, forward_starting_size, backward_starting_size, restartSize, restart, warmUp, 0); 
    SpinQuantum hq(0, SpinSpace(0), IrrepSpace(0));

    p2out << system<<endl;

    std::vector<Matrix> Rotationa, Rotationb;
    if (mpigetrank() == 0) {
      Rotationa = statea.getSiteTensors(0);
      Rotationb = stateb.getSiteTensors(0);
    }

#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world, Rotationa, 0);
    mpi::broadcast(world, Rotationb, 0);
#endif

    system.transform_operators(const_cast<std::vector<Matrix>&>(Rotationa), 
			       const_cast<std::vector<Matrix>&>(Rotationb), false, false );

    int sys_add = true; bool direct = true; 

    for (int i=0; i<MPS::sweepIters-1; i++) {
      SpinBlock newSystem;
      system.addAdditionalCompOps();

      SpinBlock dotsite(i+1, i+1, 0, false);
      if (mpigetrank() == 0) {
	Rotationa = statea.getSiteTensors(i+1);
	Rotationb = stateb.getSiteTensors(i+1);
      }

#ifndef SERIAL
      mpi::broadcast(world, Rotationa, 0);
      mpi::broadcast(world, Rotationb, 0);
#endif
      InitBlocks::InitNewSystemBlock(system, dotsite, newSystem, 0, statebindex, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true);

      newSystem.transform_operators(const_cast<std::vector<Matrix>&>(Rotationa), 
				    const_cast<std::vector<Matrix>&>(Rotationb));

      system = newSystem;
    }

    SpinBlock newSystem, big;
    SpinBlock dotsite1(MPS::sweepIters, MPS::sweepIters, 0, false);
    SpinBlock dotsite2(MPS::sweepIters+1, MPS::sweepIters+1, 0, false);


    system.addAdditionalCompOps();
    InitBlocks::InitNewSystemBlock(system, dotsite1, newSystem, 0, statebindex, sys_add, direct, 0, DISTRIBUTED_STORAGE, false, true);
    
    newSystem.set_loopblock(false); system.set_loopblock(false);
    InitBlocks::InitBigBlock(newSystem, dotsite2, big); 
    
    Wavefunction temp = statea.getw();
    temp.Clear();

    big.multiplyH(const_cast<Wavefunction&>(stateb.getw()), &temp, 1);

    if (mpigetrank() == 0)
      h = DotProduct(statea.getw(), temp);

    temp.Clear();
    big.multiplyOverlap(const_cast<Wavefunction&>(stateb.getw()), &temp, 1);
    if (mpigetrank() == 0)
      o = DotProduct(statea.getw(), temp);

#ifndef SERIAL
    mpi::broadcast(world, h, 0);
    mpi::broadcast(world, o, 0);
#endif

    return;
  }


}

      
      


