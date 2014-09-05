
#include "nevpt2_renormalize.h"
#include "wavefunction.h"
#include "operatorfunctions.h"
#include "density.h"
#include "rotationmat.h"
#include "nevpt2_operators.h"
#include "nevpt2_info.h"
#include "nevpt2_util.h"

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

namespace SpinAdapted{
  
  //============================================================================
  // Sum up Wavefunctions that are generated from different processes
  //============================================================================
  void SumWavefunctionFromProcesses(Wavefunction &wave,SpinBlock &big){
    
#ifndef SERIAL
    mpi::communicator world;
    int rank = world.rank();
    int size = world.size();
    if (size==1){
      return;
    }
    else{
      //------
      //master
      //------
      if (rank==0){
        //generate a correctly size Wavefunction
        SpinQuantum WFQ(wave.get_deltaQuantum(0));
        Wavefunction tmprecv(WFQ,&big,true);
        for (int i=1;i<size;i++){
          //receive the contribtuions
          world.recv(i, i, tmprecv);
          //sum them up
          ScaleAdd(1.0,tmprecv,wave);
        }//slaves
      }//master
      //-----
      //slave
      //-----
      else{
        world.send(0, rank, wave);  
      }
      //eventually make the complete wavefunction available to everyone
      mpi::broadcast(world,wave,0);
    }//size>1
#endif
  }
  
  //============================================================================
  // Evaluate the rotation matrices for the states that are out of scope of the 
  // regular DMRG algorithm because of their quanta. These states are important
  // subsequent NEVPT2 calculations, in particular for the V(ij) class of perturber 
  // functions
  //============================================================================
  void NEVPT2_AddToRotationMat(SpinBlock &big, vector<Matrix> &RotMatrix, SpinBlock &system,
                               vector<Wavefunction> &WF, int sweepIter){
    
    char msg[512];
    int t,u;
    double twoS_,twoS,S,S_;
    double fac;
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();

    //--------------------------------------------------------------------------
    //generate the density that is used to determine the "important" states
    //--------------------------------------------------------------------------
    //the wavefunction quanta
    SpinQuantum WFQ = WF[0].get_deltaQuantum(0);
    twoS = (double) WFQ.get_s().getirrep();
    S = twoS/2.0;
    //the oeprator quanta of a creation operator
    SpinQuantum OpQT= SpinQuantum(2,SpinSpace(2),IrrepSpace(0));
    SpinQuantum OpQS= SpinQuantum(2,SpinSpace(0),IrrepSpace(0));
    //the possible quanta for the resulting wavefunction
    vector<SpinQuantum> VQT= OpQT+ WFQ;
    vector<SpinQuantum> VQS= OpQS+ WFQ;
    //generate the wavefunction that is used to generate the density
    vector<Wavefunction> wave;
    for (int iquanta=0;iquanta<VQT.size();iquanta++){
      Wavefunction wf(VQT[iquanta],&big,true);
      wave.push_back(wf);
    }
      for (int iquanta=0;iquanta<VQS.size();iquanta++){
        Wavefunction wf(VQS[iquanta],&big,true);
      wave.push_back(wf);
    }
    //the triplet prefactors
    vector<double> Fac;
    //evaluate the prefactors
    for (int iquanta=0;iquanta<VQT.size();iquanta++){
      EvalCompFactors(Fac,2,VQT[iquanta].get_s().getirrep(),WFQ.get_s().getirrep());
    }
    for (int iroot=0;iroot<WF.size();iroot++){
      //-------------
      //leftBlock
      //-------------
      for (int tucount=0;tucount<leftBlock->get_op_array(CRE_CRE).get_size();tucount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tucount)[0]->getworkingrepresentation(leftBlock);
        boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tucount)[1]->getworkingrepresentation(leftBlock);
        //get the orbital labels
        t = SingOp->get_orbs(0);
        u = SingOp->get_orbs(1);
        //the triplet case
        for (int iquanta=0;iquanta<VQT.size();iquanta++){
          //the wavefunction that holds the result of the multiplication 
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQT[iquanta],&big,true));
          //do the multiplication with the wavefunction
          operatorfunctions::TensorMultiply(leftBlock,*TripOp,&big,WF[iroot],*Vtu,OpQT,Fac[iquanta]);
          //add to the total wavefunction
          ScaleAdd(1.0,*Vtu,wave[iquanta]);
        }//iquanta
        //the singlet case
        for (int iquanta=0;iquanta<VQS.size();iquanta++){
          //the wavefunction that holds the result of the multiplication 
          boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQS[iquanta],&big,true));
          //do the multiplication with the wavefunction
          operatorfunctions::TensorMultiply(leftBlock,*SingOp,&big,WF[iroot],*Vtu,OpQS,1.0);
          //add to the total wavefunction
          int offset = VQT.size();
          ScaleAdd(1.0,*Vtu,wave[iquanta+offset]);
        }//iquanta
      }//t and u on the left side
    }//iroot
    //evaluate the norm of the density
    double Norm=0.0;
    for (int iquanta=0;iquanta<wave.size();iquanta++){
      //Sum contributions from all processes
      SumWavefunctionFromProcesses(wave[iquanta],big);
      //evaluate the total norm
      Norm += DotProduct(wave[iquanta],wave[iquanta]);
    }
    //--------------------
    //generate the density
    //--------------------
    vector<double> wt;
    wt.resize(wave.size(),1.0);//we do not put any weights on the different functions
    SpinBlock *oldSystem = big.get_leftBlock();
    DensityMatrix tracedMatrix;
    tracedMatrix.allocate(oldSystem->get_stateInfo());
    tracedMatrix.makedensitymatrix(wave,big,wt,0.0,0.0,0.0);//we do not add noise; this might be changed at some point
    Scale(1/Norm,tracedMatrix);
    //-----------------------
    //diagonalize the density
    //-----------------------
    if (mpi_rank()==0){
      DensityMatrix transformMatrix;
      transformMatrix.allocate(oldSystem->get_stateInfo());
      vector<DiagonalMatrix> eigenMatrix;
      diagonalise_dm(tracedMatrix, transformMatrix, eigenMatrix);
      //DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //print the eigenmatrix
      /*
      for (int iquanta=0;iquanta<eigenMatrix.size();iquanta++){
        int n = oldSystem->get_stateInfo().quanta[iquanta].get_n();
        int s = oldSystem->get_stateInfo().quanta[iquanta].get_s().getirrep();
        sprintf(msg,"\nN=%i   S=%i",n,s);pout << msg;
        for (int istate=0;istate<eigenMatrix[iquanta].Nrows();istate++){
          sprintf(msg," %4.12lf",eigenMatrix[iquanta].element(istate,istate));
          pout << msg;
        }//istate
      }//iquanta
      */
      //DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //--------------------------------------------------------------------------
      //determine the quanta that are out of scope of the regular DMRG algorithm
      //--------------------------------------------------------------------------
      int TotalN = dmrginp.effective_molecule_quantum().get_n();//particle number
      int TotalS = dmrginp.effective_molecule_quantum().get_s().getirrep();//total spin
      IrrepSpace TotalSymm = dmrginp.effective_molecule_quantum().get_symm();//symmetry

      vector<SpinQuantum> SystemQuanta = big.get_stateInfo().leftStateInfo->quanta;

      vector<bool> Include(SystemQuanta.size());
    
      //loop over all quanta in the system block
      for (int iquanta=0;iquanta<SystemQuanta.size();iquanta++){
        //get the particle number, spin and symmetry
        int n = SystemQuanta[iquanta].get_n();
        int s = SystemQuanta[iquanta].get_s().getirrep();
        IrrepSpace symm = SystemQuanta[iquanta].get_symm();

        //initialize the inclusion flag
        Include[iquanta]=false;

        //these conditions are reproduced from the Tensorproduct function in StateInfo.C
        if (n>TotalN) Include[iquanta] = true;
        if ((n==TotalN)&&(s!=TotalS)) Include[iquanta] = true;
        if (s+(TotalN-n)<TotalS) Include[iquanta] = true;
        if (s-(TotalN-n)>TotalS) Include[iquanta] = true;
        if ((n==TotalN)&&(symm!=TotalSymm)) Include[iquanta] = true;

        //if a set of quanta is included, make sure, the corresponding rotationmatrix 
        //is empty (we do not want to keep leftovers from the regular algorithm)
        if (Include[iquanta]){
          if (RotMatrix[iquanta].Ncols()>0){
            RotMatrix[iquanta].ReSize(0,0);
          }
        }

      }//iquanta

      //--------------------------------------------------------------------------
      //order the weights (this is largely reproduced and only slightly altered 
      //from the sort_weights functions
      //--------------------------------------------------------------------------
      // first sort weights with a multimap
      multimap<double, pair<int,int> > weightmap;
      int nquanta = eigenMatrix.size();
      vector<double> TotalWeightsOfQuanta(nquanta);
      for (int iquanta=0;iquanta<nquanta;iquanta++){
        //only include quanta if they are out of scope of the regular DMRG algorithm
        if (Include[iquanta]==false) continue;
        for (int istate=0;istate<eigenMatrix[iquanta].Nrows();istate++){
          //determine the weight
          double Weight=eigenMatrix[iquanta].element(istate,istate);
          //make an index pair of the quantum numbers and the state number
          pair<int,int> Indices(iquanta,istate);
          //make a pair of the weight and the index pair
          pair <double, pair<int,int> > Weight_Quanta_State(Weight,Indices);
          //insert this pair into a list (the insertion position is chosen according
          //to the size of weight)
          weightmap.insert(Weight_Quanta_State);
          TotalWeightsOfQuanta[iquanta] += Weight;
        }//istate
      }//iquanta
      // then create an ordered index list (ordered with respect to the weight)
      multimap<double, pair<int,int> >::reverse_iterator w = weightmap.rbegin();
      vector<pair<int, int> > orderedIndices;
      while (w!=weightmap.rend()){
        pair<int, int> IndexPair(w->second.first,w->second.second);
        orderedIndices.push_back(IndexPair);
        ++w;  
      }

      //--------------------------------------------------------------------------
      // generate the rotation matrices 
      //--------------------------------------------------------------------------
      int StatesLeft=0;
      StatesLeft = dmrginp.kept_nevpt2_states();
      //if the user has not specified the number of kept states, 10% of the number
      //of the "normal" states will be kept
      int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
      if (StatesLeft<0) StatesLeft = MaxM / 10;
      int iIndex=0;
      //loop over all kept states
      while (StatesLeft>0){
        //if we have included all states already, leave the loop
        if (iIndex==orderedIndices.size()) break;
        //get the quanta and state indices of the current state
        int iquanta = orderedIndices[iIndex].first;
        int istate  = orderedIndices[iIndex].second;
        //add the corresponding columns to the rotation matrices
        if (eigenMatrix[iquanta].element(istate,istate)>-0.01){//this restriction needs to be refined at some point
          if (RotMatrix[iquanta].Ncols()==0){
            RotMatrix[iquanta] = transformMatrix(iquanta,iquanta).Column(istate+1);
          }
          else{
            RotMatrix[iquanta] |= transformMatrix(iquanta,iquanta).Column(istate+1);
          }
        }
        StatesLeft--;
        iIndex++;
      }//kept states
    }//master
  }


  //============================================================================
  // Generate a (normalized) density from a set of perturber functions and the
  // reference function
  //============================================================================
  void NEVPT2MakeSpecialDensity(SpinBlock &big,vector <vector<WavefunctionArray> > &T,NEVPT2Info &Info, 
                                DensityMatrix &TotalDensity, vector<Wavefunction> &references){
    int dummy;
    int i;
    double TotalNorm=0.0;
    char msg[512];
    double RefDensWeight = Info.GetRefDensWeight();
    
    sprintf(msg,"\nMaking special Densities RefWeight=%4.12lf",RefDensWeight);pout << msg;
    //----------------------------------------------
    //build the density from the perturber functions
    //----------------------------------------------
    DensityMatrix PerturberDensity;
    PerturberDensity.allocate(big.get_leftBlock()->get_stateInfo());
    for (int iroot=0;iroot<T.size();iroot++){
      for (int iquanta=0;iquanta<T[iroot].size();iquanta++){
        //open the file
        T[iroot][iquanta].OpenFileRead();
        //get the quanta
        SpinQuantum Q = T[iroot][iquanta].GetQuanta();
        //generate the wavefunction that will generate the density
        boost::shared_ptr<Wavefunction> WF = T[iroot][iquanta].GetBasisOperator();
        //evaluate the norm
        double Norm = DotProduct(*WF,*WF);
        //build the partial density
        vector<Wavefunction> tmpwf;
        tmpwf.push_back(*WF);
        DensityMatrix D;
        D.allocate(big.get_leftBlock()->get_stateInfo());
        D.makedensitymatrix(tmpwf, big, std::vector<double>(1,1.0), 0.0, 0.0, false);
        //add it to the total density
        ScaleAdd(1.0/Norm,D,PerturberDensity);
        TotalNorm +=1.0;
        //close the file
        T[iroot][iquanta].CloseFileRead();
      }//iquanta
    }//iroot
    //if required, normalize the perturber density 
    if (TotalNorm != 1.0) Scale(1.0/TotalNorm,PerturberDensity);
    //----------------------------------------------
    //build the density from the reference functions
    //----------------------------------------------
    DensityMatrix RefDensity;
    RefDensity.allocate(big.get_leftBlock()->get_stateInfo());
    TotalNorm=0.0;
    for (int iroot=0;iroot<references.size();iroot++){
      double Norm = DotProduct(references[iroot],references[iroot]);
      //build the partial density
      vector<Wavefunction> tmpwf;
      tmpwf.push_back(references[iroot]);
      DensityMatrix D;
      D.allocate(big.get_leftBlock()->get_stateInfo());
      D.makedensitymatrix(tmpwf, big, std::vector<double>(1,1.0), 0.0, 0.0, false);
      //add it to the total density
      ScaleAdd(1.0/Norm,D,RefDensity);
      TotalNorm +=1.0;
    }//iroot
    if (TotalNorm != 1.0) Scale(1.0/TotalNorm,RefDensity);
    //-----------------------
    //Build the Total Density
    //-----------------------
    double fac = RefDensWeight;
    double fac_= 1-RefDensWeight;
    TotalDensity.allocate(big.get_leftBlock()->get_stateInfo());
    ScaleAdd(fac,RefDensity,TotalDensity);
    ScaleAdd(fac_,PerturberDensity,TotalDensity);
    
    
  }
  
  //============================================================================
  // Evaluate the Norm (trace) of a density matrix
  //============================================================================
  double DensityNorm(const DensityMatrix &D){
    double Norm = 0.0;
    for(int lQ=0;lQ<D.nrows();++lQ){
      for(int rQ=0;rQ<D.ncols();++rQ){
        if(D.allowed(lQ,rQ)){
          for(int i=0;i<D(lQ,rQ).Nrows();++i){
          Norm += D(lQ,rQ)(i+1,i+1);
          }//i
        }//allowed
      }//rq
    }//lq
    return Norm;
  }
  
  //============================================================================
  // Add the first order interacting space density to the density of the ground 
  // state (with a weight of 10%)
  //============================================================================
  void AddFOISDensity(SpinBlock &big, DensityMatrix &D, vector<Wavefunction> &WF,
                      NEVPT2Info &Info, int SweepIter, int BlockIter){
    double TotalNorm=0.0;
    char msg[512];
    double RefDensWeight = Info.GetRefDensWeight();
    int nroots = WF.size();
    
    if (1.0-RefDensWeight>1e-12){
      //Evaluate the Norm (trace) of the density
      double RefDensityNorm = DensityNorm(D);
      //sprintf(msg,"\nNorm of reference density N=%4.12lf",RefDensityNorm);pout << msg;
      //sprintf(msg,"\nAdding FOIS density with a weight of w=%4.12lf",1.0-RefDensWeight);pout << msg;
      //-------------------------
      //generate the FOIS density
      //-------------------------
      vector<Wavefunction> FOIS;
      for (int iroot=0;iroot<nroots;iroot++){
        SpinQuantum WFQ = WF[iroot].get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> FOISElement(new Wavefunction(WFQ,&big,true));
        //generate the action of the Hamiltonian on the reference function
        big.multiplyH(WF[iroot],FOISElement.get(),MAX_THRD);
        //store the resulting Wavefunction
        FOIS.push_back(*FOISElement);
      }
      DensityMatrix FOISDensity;
      FOISDensity.allocate(big.get_leftBlock()->get_stateInfo());
      FOISDensity.makedensitymatrix(FOIS,big,dmrginp.weights(SweepIter),0.0,0.0,false);
      double FOISDensityNorm = DensityNorm(FOISDensity);
      //sprintf(msg,"\nNorm of FOIS density N=%4.12lf",FOISDensityNorm);pout << msg;

      //-------------------------------------------
      //add the FOIS density to the regular density
      //-------------------------------------------
      //first normalize the FOISDensity
      Scale(RefDensityNorm/FOISDensityNorm,FOISDensity);
      //renormalize the reference density
      Scale(RefDensWeight,D);
      //then add them up
      ScaleAdd((1-RefDensWeight),FOISDensity,D);
      FOISDensity.Clear();

      double TotalDensityNorm = DensityNorm(D);
      //sprintf(msg,"\nNorm of total density N=%4.12lf",TotalDensityNorm);pout << msg;
    }//Reference weight < 1
  }
  
}





