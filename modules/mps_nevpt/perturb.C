#include "perturb.h"


//TODO
//Vtu in different processor should be stored together?

void SpinAdapted::mps_nevpt::DDsinglePerturb::store(const Wavefuntion& w, const SpinBlock& big)
{
  int NActive = dmrginp.act_size();

  SpinBlock *leftBlock  = big.get_leftBlock();
  SpinBlock *rightBlock = big.get_rightBlock();
  boost::shared_ptr<WavefunctionArray> VtuSing(new WavefunctionArray);

  Matrix Sign(NActive,NActive);
  for (int t=0;t<NActive;t++){
    for (int u=0;u<NActive;u++){
      Sign.element(t,u) = 1.0;
    }//u
  }//t

  int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
  VtuSing->CalcBufferSize(MaxM);

      //--------------------------------------------------------------------------
      // generate V(t,u) = a(t)a(u)|psi> 
      //--------------------------------------------------------------------------
      //open the container
      for (int i=0;i<VtuArray.size();i++){
        VtuArray[i]->OpenFileWrite();
      }

      //--------------------------------------------------------
      //case 1: both indices are on the left side of the lattice
      //--------------------------------------------------------
      for (tu=0;tu<leftBlock->get_op_array(DES_DES).get_size();tu++){

        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = leftBlock->get_op_array(DES_DES).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
        //get the orbital indices: Note: since we will always deal with the conjugated
        //operator, the orbital indices are reversed
        int t = SingOp->get_orbs(0);
        int u = SingOp->get_orbs(1);
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=w.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ = SingOp->get_deltaQuantum(0); 
        vector<SpinQuantum> VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> w1(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(leftBlock,*SingOp,&big,w,*VtuSing,opQ,1.0);//the minus sign arises from the conjugation
        //------------------
        //store the operator
        //------------------
        VtuSing->AppendOperator(w1,t,u);
        if (t!=u) VtuSing->DuplicateAddress(t,u,u,t);
      }//tu
      //--------------------------------------------------------
      //case 2: both indices are on the left side of the lattice
      //--------------------------------------------------------
      for (tu=0;tu<rightBlock->get_op_array(DES_DES).get_size();tu++){

        //------------------------
        boost::shared_ptr<SparseMatrix> SingOp = rightBlock->get_op_array(DES_DES).get_local_element(tu)[0]->getworkingrepresentation(rightBlock);
        int t = SingOp->get_orbs(0);
        int u = SingOp->get_orbs(1);
        //get the Quanta of the Wavefunction
        SpinQuantum WFQ=w.get_deltaQuantum(0);
        //the operator quanta
        SpinQuantum opQ = SingOp->get_deltaQuantum(0); 
        vector<SpinQuantum> VQ = opQ + WFQ;
        //generate the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> w1(new Wavefunction(VQ[0],&big,true));
        //multiply the wavefunction with the singlet operator
        operatorfunctions::TensorMultiply(rightBlock,*SingOp,&big,w,*VtuSing,opQ,1.0);//the minus sign arises from the conjugation
        //------------------
        //store the operator
        //------------------
        VtuSing->AppendOperator(w1,t,u);
        if (t!=u) VtuSing->DuplicateAddress(t,u,u,t);
      }//tu

      //--------------------------------------------------------
      //case 3: the indices are on the both sides of the lattice
      //--------------------------------------------------------
      int tstart,tstop;
      //divide the loop over the parallel processes
      DivideLoop(tstart,tstop,0,big.leftBlock->get_op_array(DES).get_size());
      //loop over local elements on the left side
      for (int tcount=tstart;tcount<tstop;tcount++){
        boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(DES).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
        //get the orbital index
        u = lop->get_orbs(0);
        for (int ucount=0;ucount<rightBlock->get_op_array(DES).size();ucount++){
          //get the operator
          boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(DES).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
          //get the orbital index
          t = rop->get_orbs(0);
          //get the Quanta of the Wavefunction
          std::vector<SpinQuantum> opQu = lop->get_deltaQuantum() + rop->get_deltaQuantum(); 
          std::vector<SpinQuantum> VQ = opQu[0] + ketquanta;
          for(int iquanta=0; iquanta<VQ.size(); iquanta++){
            boost::shared_ptr<Wavefunction> w1(new Wavefunction(VQ[iquanta],&big,true));
            //multiply the wavefunction with the triplet operator
            operatorfunctions::TensorMultiply(leftBlock,*lop,*rop,&big,WF,*Vtu,opQ,1.0);
            VtuSinge->AppendOperator(w1,t,u);
            VtuSinge->DuplicateAddress(t,u,u,t);
            Sign.element(u,t)=-1.0;
          }//iquanta
        }
      }
        VtuSing->CloseFileWrite();

}


  // Evaluate the action of H on a function V(tu), meaning:
  //            sigma(t,u) = H_0 * V(t,u)
  //                       = H_0 * O(t,u) |psi>
  //============================================================================
void SpinAdapted::mps_nevpt::DDsinglePerturb::run(const Wavefuntion& w, const SpinBlock& big, IntegralContainer &IAJB, IntegralContainer &IKJL)
{
  int NActive = dmrginp.act_size();

  SpinBlock *leftBlock  = big.get_leftBlock();
  SpinBlock *rightBlock = big.get_rightBlock();

  SigmaSing->Initialize(NActive,msg,mpi_world_size());
  SigmaSing->CalcBufferSize(MaxCore,MaxM);

}

  // Evaluate the action of H on a function V(tu), meaning:
  //            sigma(t,u) = H * V(t,u)
  //                       = H * O(t,u) |psi>
  //============================================================================
  void GenerateActionOfH(boost::shared_ptr<WavefunctionArray>  &Vtu,  boost::shared_ptr <WavefunctionArray>  &Sigmatu, perturb& pb, const SpinBlock& big){
    char msg[512];
    int t,u;
    //define the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    
    //-----------------------------------------------
    //get the quanta of the ground state Wavefunction
    //-----------------------------------------------
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    IrrepSpace Dummy(0);
    vector<SpinQuantum> VQTrip = WFQ + TripOpQ;
    vector<SpinQuantum> VQSing = WFQ + SingOpQ;
    
    //--------------------------------------------------------------
    // Generate the action of the Hamiltonian on the V(tu) functions
    //--------------------------------------------------------------
#ifndef SERIAL
    mpi::communicator world;
#endif
    //open the operator arrays
    for (int i=0;i<Vtu.size();i++){
      Vtu[i]->OpenFileRead();Vtu[i]->ResetBuffer();
      Sigmatu[i]->OpenFileWrite();
    }
    
    //get the number of processes
    int NumProc = mpi_world_size();
    //loop over all processes
    for (int iproc=0;iproc<NumProc;iproc++){
      bool Sender = (mpi_rank()==iproc);
      for (int iquanta=0;iquanta<Vtu.size();iquanta++){
        bool EndOfArray = false;
        while (!EndOfArray){
          boost::shared_ptr<Wavefunction> v;
          //----------------------
          //get the V(tu) operator
          //----------------------
          if (Sender){
            v = Vtu[iquanta]->GetOpFromBuffer(t,u,EndOfArray);
          }//sender
          else{
            if (iquanta==0){
              v = boost::make_shared<Wavefunction> (VQSing[0],&big,true);
            }//singlet
            else{
              int itmp = iquanta-1;
              v = boost::make_shared<Wavefunction> (VQTrip[itmp],&big,true);
            }//triplet
          }//receiver
#ifndef SERIAL
          mpi::broadcast(world,EndOfArray,iproc);
#endif
          //----------------------------
          //broadcast the V(tu) operator
          //----------------------------
          if (!EndOfArray){
#ifndef SERIAL
            mpi::broadcast(world,*v,iproc);
#endif
            //------------------
            //generate sigma(tu)
            //------------------
            boost::shared_ptr<Wavefunction> sigma;
            if (iquanta==0){
              //singlet functions
              sigma=boost::make_shared<Wavefunction> (VQSing[0],&big,true);
              //calculate the action of the Hamiltonian on V(tu)
              big.multiplyH_Q(*v,sigma.get(),MAX_THRD,VQSing[0]);
            }
            else{
              //triplet functions
              int itmp = iquanta-1;
              sigma=boost::make_shared<Wavefunction> (VQTrip[itmp],&big,true);
              //calculate the action of the Hamiltonian on V(tu)
              big.multiplyH_Q(*v,sigma.get(),MAX_THRD,VQTrip[itmp]);
            }
            //---------------
            //store sigma(tu)
            //---------------
            if (Sender) Sigmatu[iquanta]->AppendOperator(sigma,t,u);
          }
        }//tu (singlet)
      }//iquanta
    }//iproc
    //close the operator arrays
    for (int i=0;i<Vtu.size();i++){
      Vtu[i]->CloseFileRead();
      Sigmatu[i]->CloseFileWrite();
    }
  }
  
  //============================================================================
