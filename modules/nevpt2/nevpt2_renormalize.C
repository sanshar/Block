
#include "nevpt2_renormalize.h"
#include "wavefunction.h"
#include "operatorfunctions.h"
#include "density.h"
#include "rotationmat.h"
#include "nevpt2_operators.h"
#include "nevpt2_info.h"
#include "nevpt2_util.h"
#include "nevpt2_pal.h"
#include "distribute.h"

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
                      ThreeIndOpArray &CCD, ThreeIndOpArray &CDD, 
                      IntegralContainer &IKJL, IntegralContainer &IKJA, IntegralContainer &IAJB,
                      IntegralContainer &IJKA, NEVPT2Info &Info, int SweepIter, 
                      int BlockIter, bool WarmUp){
    double TotalNorm=0.0;
    char msg[512];
    double RefDensWeight = Info.GetRefDensWeight();
    int nroots = WF.size();
    SpinBlock* leftBlock = big.get_leftBlock();
    SpinBlock* rightBlock = big.get_rightBlock();
    int NOrbsLeft = leftBlock->get_sites().size();
    int NOrbsRight = rightBlock->get_sites().size();
    int OrbWin[6];
    Info.GetOrbWin(OrbWin);
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = big.get_sites().size();
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    double weight[3] ;weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
    double weight_[3];weight_[0]=0.0;weight_[1]=0.0;weight_[2]=0.0;
    //get the effective one-electron hamiltonian
    Matrix heff;
    Info.GetH(1,heff);
    
    if (1.0-RefDensWeight>1e-12){
      //Evaluate the Norm (trace) of the density
      double RefDensityNorm = DensityNorm(D);
      //sprintf(msg,"\nNorm of reference density N=%4.12lf",RefDensityNorm);pout << msg;
      //sprintf(msg,"\nAdding FOIS density with a weight of w=%4.12lf",1.0-RefDensWeight);pout << msg;
      //-------------------------
      //generate the FOIS density
      //-------------------------
      vector<Wavefunction> FOIS;
      DensityMatrix FOISDensity;
      FOISDensity.allocate(big.get_leftBlock()->get_stateInfo());
      for (int iroot=0;iroot<nroots;iroot++){
        SpinQuantum WFQ = WF[iroot].get_deltaQuantum(0);
        double twoS = (double) WFQ.get_s().getirrep();
        double S    = twoS/2.0;
        double twoS_,S_,fac;
        //generate all possible SpinQuanta
        vector<SpinQuantum> FOISQuanta;
        vector<SpinQuantum> OpQuanta;
        SpinQuantum opQ1 = SpinQuantum( 1,SpinSpace(1),IrrepSpace(0));OpQuanta.push_back(opQ1);//Cre
        SpinQuantum opQ2 = SpinQuantum(-1,SpinSpace(1),IrrepSpace(0));OpQuanta.push_back(opQ2);//Des
        SpinQuantum opQ3 = SpinQuantum( 2,SpinSpace(0),IrrepSpace(0));OpQuanta.push_back(opQ3);//CreCre
        SpinQuantum opQ4 = SpinQuantum( 2,SpinSpace(2),IrrepSpace(0));OpQuanta.push_back(opQ4);//CreCre
        SpinQuantum opQ5 = SpinQuantum( 0,SpinSpace(0),IrrepSpace(0));OpQuanta.push_back(opQ5);//CreDes
        SpinQuantum opQ6 = SpinQuantum( 0,SpinSpace(2),IrrepSpace(0));OpQuanta.push_back(opQ6);//CreDes
        SpinQuantum opQ7 = SpinQuantum(-2,SpinSpace(0),IrrepSpace(0));OpQuanta.push_back(opQ7);//DesDes
        SpinQuantum opQ8 = SpinQuantum(-2,SpinSpace(2),IrrepSpace(0));OpQuanta.push_back(opQ8);//DesDes
        for (int iq=0;iq<8;iq++){
          vector<SpinQuantum> tmp = WFQ + OpQuanta[iq];
          for (int itmp=0;itmp<tmp.size();itmp++){
            FOISQuanta.push_back(tmp[itmp]);
          }//itmp
        }//iq
        //-----------
        //Cre and DES
        //-----------
	//open the integral container
        IAJB.OpenFileRead();
        IKJL.OpenFileRead();
        IKJA.OpenFileRead();
        IJKA.OpenFileRead();
        //divide the loop over the parallel processes
        int tstart=0,tstop=0;
        PALDivideLoop(tstart,tstop,0,leftBlock->get_op_array(CRE).get_size());
        //for (int tcount=0;tcount<leftBlock->get_op_array(CRE).get_size();tcount++){
        for (int tcount=tstart;tcount<tstop;tcount++){
          //get the operator
          boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
          //get the orbital label
          int t = op->get_orbs(0);
          //get the operator quanta
          SpinQuantum opQ = (*op).get_deltaQuantum(0);
          SpinQuantum opQT = Transposeview(*op).get_deltaQuantum(0);
          //the possible quanta for the resulting wavefunction
          vector<SpinQuantum> VQ = opQ + WFQ;
          vector<SpinQuantum> VQT = opQT + WFQ;
          //the prefactors
          vector<double> Fac;
          //evaluate the prefactors for one-electron part
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            twoS_ = (double) VQ[iquanta].get_s().getirrep();
            S_ = twoS_/2.0;
            double fac = sqrt(2.0*S_+1.0)/sqrt(2.0*S+1.0);
            Fac.push_back(fac);
      }
          //-------------------
          //generate the weight
          //-------------------
          weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
          weight_[0]=0.0;weight_[1]=0.0;weight_[2]=0.0;
          int t_ = t+NInternal;
          //Va
          //accumulate the external one-electron contribution
          for (int a=a0;a<a1+1;a++){
            for (int iquanta=0;iquanta<VQ.size();iquanta++){
              weight[iquanta] += Fac[iquanta]*heff.element(t_,a);
            }//iquanta
          }//a
          //Vi
          //accumulate the internal one-electron contribution
          for (int i=i0;i<i1+1;i++){
            for (int iquanta=0;iquanta<VQ.size();iquanta++){
              weight_[iquanta]+= Fac[iquanta]*heff.element(t_,i);  
            }//iquanta
          }//i
          //--------------------------
          //generate the FOIS elements
          //--------------------------
          //the FOIS element and its density
          boost::shared_ptr<Wavefunction> FOISElement1(new Wavefunction(FOISQuanta,&big,true));
          boost::shared_ptr<DensityMatrix> tmp (new DensityMatrix);
          tmp->allocate(leftBlock->get_stateInfo());
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //the wavefunction that holds the result of the multiplication 
            boost::shared_ptr<Wavefunction> Vt(new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtT(new Wavefunction(VQT[iquanta],&big,true));
            //do the multiplication with the wavefunction
            operatorfunctions::TensorMultiply(leftBlock,*op,&big,WF[iroot],*Vt,opQ,1.0);
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*op),&big,WF[iroot],*VtT,opQT,1.0);
            //add the result to the FOIS
            ScaleAdd(weight_[iquanta],*Vt ,*FOISElement1);
            ScaleAdd(weight[iquanta],*VtT,*FOISElement1);
            //make the density from the FOIS function
          }//iquanta
          //make the density of the FOIS element
          tmp->makedensitymatrix(*FOISElement1,big,1.0);
          //add it to the total density
          ScaleAdd(1.0,*tmp,FOISDensity);
        }//t on the left side
        //------
        //CreCre
        //------
	for (int tu=0;tu<leftBlock->get_op_array(CRE_CRE).get_size();tu++){
          //get the triplet operator
          boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
          //get the orbital indices
          int t = TripOp->get_orbs(0);
          int u = TripOp->get_orbs(1);
          //the operator quanta
          SpinQuantum opQ=TripOp->get_deltaQuantum(0);
          //the vector that holds all possible combinations of wavefunction plus operator
          std::vector<SpinQuantum> VQ  = opQ  + WFQ;
          std::vector<SpinQuantum> VQ_ = opQ1 + WFQ;
          //-------------------
          //generate the weight
          //-------------------
          weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
          weight_[0]=0.0;weight_[1]=0.0;weight_[2]=0.0;
          int t_ = t+NInternal;
          int u_ = u+NInternal;
          //the prefactors (first part(Vij))
          vector<double> TripFac;
          for (int i=0;i<VQ.size();i++){
            EvalCompFactors(TripFac, 2, VQ[i].get_s().getirrep(), WFQ.get_s().getirrep());
          }
          //get the two-electron integrals
          boost::shared_ptr<Matrix> Ktu = IKJL.GetMatrix(t_,u_);
          //accumulate the two-electron contribution (first part (Vij))
          for (int i=0;i<NInternal;i++){
            for (int j=0;j<NInternal;j++){
              //the singlet part
              weight[0] += Ktu->element(j,i)*Ktu->element(j,i) + Ktu->element(i,j)*Ktu->element(i,j);
              //the triplet part
              for (int iquanta=0;iquanta<VQ.size();iquanta++){
                weight_[iquanta]+= TripFac[iquanta] * (Ktu->element(j,i)-Ktu->element(j,i));
              }//iquanta
            }//j
          }//i
          //--------------------------
          //generate the FOIS elements
          //--------------------------
          //the FOIS element and its density
          boost::shared_ptr<Wavefunction> FOISElement2(new Wavefunction(FOISQuanta,&big,true));
          boost::shared_ptr<DensityMatrix> tmp (new DensityMatrix);
          tmp->allocate(leftBlock->get_stateInfo());
          //Multiply the operator with the wavefunction
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunction that holds the result of the multiplication
            boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,*TripOp,&big,WF[iroot],*Vtu,opQ,1.0);
            //add the result to the FOIS
            ScaleAdd(weight_[iquanta],*Vtu,*FOISElement2);
          }//iquanta
          VQ.clear();
          //get the singlet operator
          boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
          //get the orbital indices
          t = SingOp->get_orbs(0);
          u = SingOp->get_orbs(1);
          //the operator quanta
          opQ = SingOp->get_deltaQuantum(0); 
          VQ = opQ + WFQ;
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
          //multiply the wavefunction with the singlet operator
          operatorfunctions::TensorMultiply(leftBlock,*SingOp,&big,WF[iroot],*VtuSing,opQ,1.0);
          //add the result to the FOIS
          ScaleAdd(weight[0],*VtuSing,*FOISElement2);
          //make the density of the FOIS element
          tmp->makedensitymatrix(*FOISElement2,big,1.0);
          //add it to the total density
          ScaleAdd(1.0,*tmp,FOISDensity);
        }//tu
        //------
        //CreDes
        //------
	for (int tu=0;tu<leftBlock->get_op_array(CRE_DES).get_size();tu++){
          //get the triplet operator
          boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_DES).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
          //get the orbital indices
          int t = TripOp->get_orbs(0);
          int u = TripOp->get_orbs(1);
          int t_ = t+NInternal;
          int u_ = u+NInternal;
          //the operator quanta
          SpinQuantum opQ=TripOp->get_deltaQuantum(0);
          //the vector that holds all possible combinations of wavefunction plus operator
          std::vector<SpinQuantum> VQ  = opQ  + WFQ;
          std::vector<SpinQuantum> VQ_ = opQ1 + WFQ;
          //-------------------
          //generate the weight
          //-------------------
          weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
          weight_[0]=0.0;weight_[1]=0.0;weight_[2]=0.0;
          double weightT[3],weight_T[3];
          weightT[0] =0.0;weightT[1] =0.0;weightT[2] =0.0;
          weight_T[0]=0.0;weight_T[1]=0.0;weight_T[2]=0.0;
          // the prefactors (Via)
          vector<double> TripFac;
          for (int i=0;i<VQ.size();i++){
            EvalCompFactors(TripFac, 2, VQ[i].get_s().getirrep(), WFQ.get_s().getirrep());
          }
          //get the two-electron integrals
          boost::shared_ptr<Matrix> Ktu = IKJA.GetMatrix(t_,u_);
          boost::shared_ptr<Matrix> Kut = IKJA.GetMatrix(u_,t_);
          boost::shared_ptr<Matrix> Jtu = IJKA.GetMatrix(t_,u_);
          //accumulate the two-electron contribution (first part, V(ia))
          for (int i=0;i<NInternal;i++){
            for (int a=0;a<NExternal;a++){
              weight[0] += 2.0 * Jtu->element(i,a)-Ktu->element(i,a);
              weightT[0]+= 2.0 * Jtu->element(i,a)-Ktu->element(i,a);
              for (int iquanta=0;iquanta<VQ.size();iquanta++){
                weight_[iquanta] += TripFac[iquanta] * Ktu->element(i,a);
                weight_T[iquanta]+= TripFac[iquanta] * Kut->element(i,a);
              }//iquanta
            }//b
          }//a
          //--------------------------
          //generate the FOIS elements
          //--------------------------
          //the FOIS element and its density
          boost::shared_ptr<Wavefunction> FOISElement3(new Wavefunction(FOISQuanta,&big,true));
          boost::shared_ptr<DensityMatrix> tmp (new DensityMatrix);
          tmp->allocate(leftBlock->get_stateInfo());
          //Multiply the operator with the wavefunction
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunction that holds the result of the multiplication
            boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtuT(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,*TripOp,&big,WF[iroot],*Vtu,opQ,1.0);
            if (t!=u){
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*TripOp),&big,WF[iroot],*VtuT,opQ,1.0);
            }//t!=u
            //add the result to the FOIS
            ScaleAdd(weight_[iquanta],*Vtu,*FOISElement3);
            ScaleAdd(weight_T[iquanta],*VtuT,*FOISElement3);
          }//iquanta
          VQ.clear();
          //get the singlet operator
          boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_DES).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
          //get the orbital indices
          t = SingOp->get_orbs(0);
          u = SingOp->get_orbs(1);
          //the operator quanta
          opQ = SingOp->get_deltaQuantum(0); 
          VQ = opQ + WFQ;
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
          boost::shared_ptr<Wavefunction> VtuSingT(new Wavefunction(VQ[0],&big,true));
          //multiply the wavefunction with the singlet operator
          operatorfunctions::TensorMultiply(leftBlock,*SingOp,&big,WF[iroot],*VtuSing,opQ,1.0);
          if (t!=u){
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*SingOp),&big,WF[iroot],*VtuSingT,opQ,1.0);
          }
          //add the result to the FOIS
          ScaleAdd(weight[0],*VtuSing,*FOISElement3);
          ScaleAdd(weightT[0],*VtuSingT,*FOISElement3);
          //make the density of the FOIS element
          tmp->makedensitymatrix(*FOISElement3,big,1.0);
          //add it to the total density
          ScaleAdd(1.0,*tmp,FOISDensity);
        }//tu
        boost::shared_ptr<DensityMatrix> tmp_ (new DensityMatrix);
        tmp_->allocate(leftBlock->get_stateInfo());
        //add the one-electron part
        weight[0] = 0.0;
        for (int i=0;i<NInternal;i++){
          for (int a=0;a<NExternal;a++){
            weight[0] += sqrt(2.0)*heff.element(i,a+a0);
          }//a
        }//i
        tmp_->makedensitymatrix(WF[iroot],big,weight[0]*weight[0]);
        //------
        //DesDes
        //------
	//open the integral container
        for (int tu=0;tu<leftBlock->get_op_array(CRE_CRE).get_size();tu++){
          //get the triplet operator
          boost::shared_ptr<SparseMatrix> TripOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
          //get the orbital indices: Note: since we will always deal with the conjugated
          //operator, the orbital indices are reversed
          int u = TripOp->get_orbs(0);
          int t = TripOp->get_orbs(1);
          //the operator quanta
          SpinQuantum opQ=Transposeview(*TripOp).get_deltaQuantum(0);
          //the vector that holds all possible combinations of wavefunction plus operator
          std::vector<SpinQuantum> VQ  = opQ  + WFQ;
          std::vector<SpinQuantum> VQ_ = opQ2 + WFQ;
          //-------------------
          //generate the weight
          //-------------------
          weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
          weight_[0]=0.0;weight_[1]=0.0;weight_[2]=0.0;
          int t_ = t+NInternal;
          int u_ = u+NInternal;
          //the prefactors (first part (Vab)
          vector<double> TripFac;
          for (int i=0;i<VQ.size();i++){
            EvalCompFactors(TripFac, 2, VQ[i].get_s().getirrep(), WFQ.get_s().getirrep());
          }
          //get the two-electron integrals
          boost::shared_ptr<Matrix> Ktu  = IAJB.GetMatrix(t_,u_);
          boost::shared_ptr<Matrix> Ktu_ = IKJA.GetMatrix(t_,u_);
          boost::shared_ptr<Matrix> Kut_ = IKJA.GetMatrix(u_,t_);
          //accumulate the two-electron contribution (first part,V(ab))
          for (int a=0;a<NExternal;a++){
            for (int b=a;b<NExternal;b++){
              //the singlet part
              weight[0] += Ktu->element(b,a) + Ktu->element(a,b);
              //the triplet part
              for (int iquanta=0;iquanta<VQ.size();iquanta++){
                weight_[iquanta]+= TripFac[iquanta]*(Ktu->element(b,a)+Ktu->element(a,b));
              }//iquanta
            }//b
          }//a
          //--------------------------
          //generate the FOIS elements
          //--------------------------
          //the FOIS element and its density
          boost::shared_ptr<Wavefunction> FOISElement4(new Wavefunction(FOISQuanta,&big,true));
          boost::shared_ptr<DensityMatrix> tmp (new DensityMatrix);
          tmp->allocate(leftBlock->get_stateInfo());
          //Multiply the operator with the wavefunction
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunction that holds the result of the multiplication
            boost::shared_ptr<Wavefunction> Vtu(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*TripOp),&big,WF[iroot],*Vtu,opQ,1.0);
            //add the result to the FOIS
            ScaleAdd(weight_[iquanta],*Vtu,*FOISElement4);
          }//iquanta
          VQ.clear();
          //get the singlet operator
          boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
          //get the orbital indices: Note: since we will always deal with the conjugated
          //operator, the orbital indices are reversed
          t = SingOp->get_orbs(1);
          u = SingOp->get_orbs(0);
          //the operator quanta
          opQ = Transposeview(*SingOp).get_deltaQuantum(0); 
          VQ = opQ + WFQ;
          //generate the wavefunction that holds the result of the multiplication
          boost::shared_ptr<Wavefunction> VtuSing(new Wavefunction(VQ[0],&big,true));
          //multiply the wavefunction with the singlet operator
          operatorfunctions::TensorMultiply(leftBlock,Transposeview(*SingOp),&big,WF[iroot],*VtuSing,opQ,-1.0);//the minus sign arises from the conjugation
          //add the result to the FOIS
          ScaleAdd(weight[0],*VtuSing,*FOISElement4);
          //make the density of the FOIS element
          tmp->makedensitymatrix(*FOISElement4,big,1.0);
          //add it to the total density
          ScaleAdd(1.0,*tmp,FOISDensity);
        }//tu
        //---------
        //CreCreDes
        //---------
	int iSweep = Info.GetNevSweep();
        int MaxIter= Info.GetMaxBlockIter(iSweep);
        if (!WarmUp&&BlockIter<=MaxIter+1){
          int t_,u_,v_;
          int t,u,v;
          boost::shared_ptr<Matrix> Kut;
          //Open the operator arrays and prepare the buffer
          CCD.OpenFileRead();
          CCD.ResetBuffer();
          //initialize the auxiliary indices
          t_ = -1;
          u_ = -1;
          v_ = -1;
          //start the loop
          bool EndOfArray = false;
          bool NeedTensorTrace = false;
          while (!EndOfArray){
            //get the operator
            boost::shared_ptr<Cre> Otuv = CCD.GetOpFromBuffer(t,u,v,NeedTensorTrace,EndOfArray);
            if (!EndOfArray){
              //if necessary, bring the operator in the current representation.
              //Note: This only needs to be done for the leftBlock. The operators on the right
              //Block are already in the correct representation
              if (NeedTensorTrace){
                boost::shared_ptr<Cre> NewOp(new Cre);
                NewOp->set_orbs() = Otuv->get_orbs();
                NewOp->set_initialised() = true;
                NewOp->resize_deltaQuantum(1);
                NewOp->set_deltaQuantum(0) = Otuv->get_deltaQuantum(0);
                NewOp->allocate(leftBlock->get_stateInfo());
                operatorfunctions::TensorTrace(leftBlock->get_leftBlock(),*Otuv,leftBlock,&(leftBlock->get_stateInfo()),*NewOp);
                //rename the new operator and free the old one
                Otuv.reset();
                Otuv = NewOp;
                NewOp.reset();
              }
              //if necessary, get the integral matrix
              if ((t!=t_)||(u!=u_)){
                 Kut = IKJL.GetMatrix(u+NInternal,t+NInternal);
              }
              //the possible quanta
              SpinQuantum opQ = Otuv->get_deltaQuantum(0);
              vector<SpinQuantum> VQ = WFQ + opQ;

              //generate the weight
              weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
              for (int i=0;i<NInternal;i++){
                for (int iquanta=0;iquanta<VQ.size();iquanta++){
                  //evaluate the prefactor
                  twoS_ = (double) VQ[iquanta].get_s().getirrep();
                  S_ = twoS_/2.0;
                  fac = sqrt((2.0*S_+1)/(2.0*S+1));
                  weight[iquanta] += fac * Kut->element(v+NInternal,i);
                }//iquanta
              }//a
              //the FOIS element and its density 
              boost::shared_ptr<Wavefunction> FOISElement5(new Wavefunction(FOISQuanta,&big,true));
              boost::shared_ptr<DensityMatrix> tmp (new DensityMatrix);
              tmp->allocate(leftBlock->get_stateInfo());
              for (int iquanta=0;iquanta<VQ.size();iquanta++){
                //generate V(tuv) = O(tuv)|psi>
                boost::shared_ptr<Wavefunction> Vtuv(new Wavefunction(VQ[iquanta],&big,true));
                operatorfunctions::TensorMultiply(leftBlock,*Otuv,&big,WF[iroot],*Vtuv,opQ,1.0);
                //add the result to the FOIS
                ScaleAdd(weight[iquanta],*Vtuv,*FOISElement5);
              }//iquanta
              //make the density of the FOIS element
              tmp->makedensitymatrix(*FOISElement5,big,1.0);
              //add it to the total density
              ScaleAdd(1.0,*tmp,FOISDensity);
            }//!EndOfArray
            t_ = t;
            u_ = u;
            v_ = v;
          }//tuv
          //close the Operator array
          CCD.CloseFileRead();
          //---------
          //CreDesDes
          //---------
	  boost::shared_ptr<Matrix> Ktv;
          //Open the operator arrays and prepare the buffer
          CDD.OpenFileRead();
          CDD.ResetBuffer();
          //initialize the auxiliary indices
          t_ = -1;
          u_ = -1;
          v_ = -1;
          //start the loop
          EndOfArray = false;
          NeedTensorTrace = false;
          while (!EndOfArray){
            boost::shared_ptr<Cre>  Otuv = CDD.GetOpFromBuffer(t,u,v,NeedTensorTrace,EndOfArray);
            if (!EndOfArray){
              //if necessary, bring the operator in the current representation.
              //Note: This only needs to be done for the leftBlock. The operators on the right
              //Block are already in the correct representation
              if (NeedTensorTrace){
                boost::shared_ptr<Cre> NewOp(new Cre);
                NewOp->set_orbs() = Otuv->get_orbs();
                NewOp->set_initialised() = true;
                NewOp->resize_deltaQuantum(1);
                NewOp->set_deltaQuantum(0) = Otuv->get_deltaQuantum(0);
                NewOp->allocate(leftBlock->get_stateInfo());
                operatorfunctions::TensorTrace(leftBlock->get_leftBlock(),*Otuv,leftBlock,&(leftBlock->get_stateInfo()),*NewOp);
                //rename the new operator and free the old one
                Otuv.reset();
                Otuv = NewOp;
                NewOp.reset();
              }
              //if necessary, get the integral matrix
              if ((t!=t_)||(v!=v_)){
                 Ktv = IKJA.GetMatrix(t+NInternal,v+NInternal);
              }
              //the possible quanta
              SpinQuantum opQ = Otuv->get_deltaQuantum(0);
              vector<SpinQuantum> VQ = WFQ + opQ;
              //generate the weight
              weight[0] =0.0;weight[1] =0.0;weight[2] =0.0;
              for (int a=0;a<NExternal;a++){
                for (int iquanta=0;iquanta<VQ.size();iquanta++){
                  //evaluate the prefactor
                  twoS_ = (double) VQ[iquanta].get_s().getirrep();
                  S_ = twoS_/2.0;
                  fac = sqrt((2.0*S_+1)/(2.0*S+1));
                  weight[iquanta] +=  fac * Ktv->element(u+NInternal,a);
                }//iquanta
              }//a
              //the FOIS element and its density
              boost::shared_ptr<Wavefunction> FOISElement6(new Wavefunction(FOISQuanta,&big,true));
              boost::shared_ptr<DensityMatrix> tmp (new DensityMatrix);
              tmp->allocate(leftBlock->get_stateInfo());
              for (int iquanta=0;iquanta<VQ.size();iquanta++){
                //generate V(tuv) = O(tuv)|psi>
                boost::shared_ptr<Wavefunction> Vtuv(new Wavefunction(VQ[iquanta],&big,true));
                operatorfunctions::TensorMultiply(leftBlock,*Otuv,&big,WF[iroot],*Vtuv,opQ,1.0);
                //add the result to the FOIS
                ScaleAdd(weight[iquanta],*Vtuv,*FOISElement6);
              }//iquanta
              //make the density of the FOIS element
              tmp->makedensitymatrix(*FOISElement6,big,1.0);
              //add it to the total density
              ScaleAdd(1.0,*tmp,FOISDensity);
              t_ = t;
              u_ = u;
              v_ = v;
            }
          }//tuv
          //close the Operator array
          CDD.CloseFileRead();
        }//take into account CCD and CDD
        //close the integral container
        IAJB.CloseFileRead();
        IKJA.CloseFileRead();
        IKJL.CloseFileRead();
        IJKA.CloseFileRead();
      }//root
      //gather density from different processes
      AddPalDensity(FOISDensity);
      double FOISDensityNorm = DensityNorm(FOISDensity);
      //sprintf(msg,"\nNorm of FOIS density N=%4.12lf",FOISDensityNorm);pout << msg;

      //-------------------------------------------
      //add the FOIS density to the regular density
      //-------------------------------------------
      if (FOISDensityNorm>1.0e-12){
        //first normalize the FOISDensity
        Scale(RefDensityNorm/FOISDensityNorm,FOISDensity);
        //renormalize the reference density
        Scale(RefDensWeight,D);
        //then add them up
        ScaleAdd((1-RefDensWeight),FOISDensity,D);
        FOISDensity.Clear();
      }
      
      double TotalDensityNorm = DensityNorm(D);
      //sprintf(msg,"\nNorm of total density N=%4.12lf",TotalDensityNorm);pout << msg;
    }//Reference weight < 1
  }
  
}





