
#include "spinblock.h"
#include "operatorfunctions.h"
#include "nevpt2_operators.h"
#include "nevpt2_util.h"
#include "nevpt2_opconstruct.h"
#include "nevpt2_info.h"
#include "nevpt2.h"
#include "nevpt2_pal.h"

#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif


namespace SpinAdapted{

  void CheckOperator_(SpinBlock &big, Wavefunction &WF, int t, int u, int v, SparseMatrix &Op, bool left){
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    char msg[512];
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    SpinQuantum opQ(-1,SpinSpace(1),IrrepSpace(0));
    vector<SpinQuantum> VQ = WFQ + opQ;
    SpinQuantum totOpQ(0,SpinSpace(0),IrrepSpace(0));
    //--------------
    //the left Block
    //--------------
    for (int i=0;i<leftBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
      if (left){
        //the first resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*tmpop),&big,WF,*tmpwf,opQ,sqrt(2.0));
        
        //the second resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Op,&big,*tmpwf,*tmpwf_,opQ,1.0);
        
        //build the element
        double val = DotProduct(WF,*tmpwf_);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e  left a",t,u,v,v_,val);
        pout << msg;
      }//Op is on left Block
      else{
        
        //the resulting wavefunction
        vector<SpinQuantum> totOpQ = opQ - tmpop->get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*tmpop),Op,&big,WF,*tmpwf,totOpQ[0],sqrt(2.0));
        //build the element
        double val = DotProduct(WF,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e  left a",t,u,v,v_,val);
        pout << msg;
      }
    }

    //---------------
    //the right Block
    //---------------
    for (int i=0;i<rightBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = rightBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(rightBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
      if (!left){
        //the first resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*tmpop),&big,WF,*tmpwf,opQ,sqrt(2.0));
        
        //the second resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Op,&big,*tmpwf,*tmpwf_,opQ,1.0);
        
        //build the element
        double val = DotProduct(WF,*tmpwf_);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e  left a",t,u,v,v_,val);
        pout << msg;
      }//Op is on left Block
      else{
        
        //the resulting wavefunction
        vector<SpinQuantum> totOpQ = opQ - tmpop->get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*tmpop),Op,&big,WF,*tmpwf,totOpQ[0],sqrt(2.0));
        //build the element
        double val = DotProduct(WF,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e  left a",t,u,v,v_,val);
        pout << msg;
      }
    }

    
  }
  
  void CheckOperator(SpinBlock &big, Wavefunction &WF, int t,SparseMatrix &Op, bool left){
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    char msg[512];
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    SpinQuantum opQ = Op.get_deltaQuantum(0);
    vector<SpinQuantum> VQ = WFQ + opQ;
    SpinQuantum totOpQ(0,SpinSpace(0),IrrepSpace(0));
    
    //--------------
    //the left Block
    //--------------
    for (int i=0;i<leftBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
      if (left){
        //the first resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Op,&big,WF,*tmpwf,opQ,sqrt(2.0));
        
        //the second resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(leftBlock,*tmpop,&big,WF,*tmpwf_,opQ,sqrt(2.0));
        
        //build the element
        double val = DotProduct(*tmpwf_,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i)=%4.8e  left a",v_,t,val);
        pout << msg;
      }//Op is on left Block
      
      else{
        
        //the resulting wavefunction
        vector<SpinQuantum> totOpQ = opQ - tmpop->get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*tmpop),Op,&big,WF,*tmpwf,totOpQ[0],sqrt(2.0));
        //build the element
        double val = DotProduct(WF,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i)=%4.8e  left b",v_,t,val);
        pout << msg;
      }
    }

    //---------------
    //the right Block
    //---------------
    for (int i=0;i<rightBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = rightBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(rightBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
      if (!left){
        //the first resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Op,&big,WF,*tmpwf,opQ,sqrt(2.0));
        
        //the second resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(rightBlock,*tmpop,&big,WF,*tmpwf_,opQ,sqrt(2.0));
        
        //build the element
        double val = DotProduct(*tmpwf_,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i)=%4.8e  right a",v_,t,val);
        pout << msg;
      }//Op is on left Block
      
      else{
        //the resulting wavefunction
        vector<SpinQuantum> totOpQ = opQ - tmpop->get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*tmpop),Op,&big,WF,*tmpwf,totOpQ[0],sqrt(2.0));
        //build the element
        double val = DotProduct(WF,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i)=%4.8e  right b",v_,t,val);
        pout << msg;
      }
    }
  }
  
  
  void CheckOperator(SpinBlock &big, Wavefunction &WF, int t, int u, int v, SparseMatrix &Op, bool left){
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    char msg[512];
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    SpinQuantum opQ = Op.get_deltaQuantum(0);
    vector<SpinQuantum> VQ = WFQ + opQ;
    SpinQuantum totOpQ(0,SpinSpace(0),IrrepSpace(0));
    
    //--------------
    //the left Block
    //--------------
    for (int i=0;i<leftBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
      if (left){
        //the first resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Op,&big,WF,*tmpwf,opQ,sqrt(2.0));
        
        //the second resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*tmpop),&big,WF,*tmpwf_,opQ,sqrt(2.0));
        
        //build the element
        double val = DotProduct(*tmpwf_,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e",v_,t,v,u,val);
        pout << msg;
      }//Op is on left Block
      
      else{
        
        //the resulting wavefunction
        vector<SpinQuantum> totOpQ = opQ + tmpop->get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(leftBlock,*tmpop,Op,&big,WF,*tmpwf,totOpQ[0],sqrt(2.0));
        //build the element
        double val = DotProduct(WF,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e",v_,t,v,u,val);
        pout << msg;
      }
    }

    //---------------
    //the right Block
    //---------------
    for (int i=0;i<rightBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = rightBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(rightBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
      if (!left){
        //the first resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Op,&big,WF,*tmpwf,opQ,sqrt(2.0));
        
        //the second resulting wavefunction
        boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(VQ[0],&big,true));
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*tmpop),&big,WF,*tmpwf_,opQ,sqrt(2.0));
        
        //build the element
        double val = DotProduct(*tmpwf_,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e",v_,t,v,u,val);
        pout << msg;
      }//Op is on left Block
      
      else{
        //the resulting wavefunction
        vector<SpinQuantum> totOpQ = opQ + tmpop->get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> tmpwf(new Wavefunction(WFQ,&big,true));
        operatorfunctions::TensorMultiply(leftBlock,Op,*tmpop,&big,WF,*tmpwf,totOpQ[0],sqrt(2.0));
        //build the element
        double val = DotProduct(WF,*tmpwf);
        sprintf(msg,"\nCheckOperator: D(%i,%i,%i,%i)=%4.8e",v_,t,v,u,val);
        pout << msg;
      }
    }
  }
  
  void CheckWF(SpinBlock &big, Wavefunction &WF, int t, int u, int v, Wavefunction &Op){
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    char msg[512];
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    SpinQuantum opQ = SpinQuantum(-1,SpinSpace(1),IrrepSpace(0));
    vector<SpinQuantum> VQ = WFQ + opQ;
    SpinQuantum totOpQ(0,SpinSpace(0),IrrepSpace(0));
    
    //--------------
    //the left Block
    //--------------
    for (int i=0;i<leftBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = leftBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(leftBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);
        
      //the second resulting wavefunction
      boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(VQ[0],&big,true));
      operatorfunctions::TensorMultiply(leftBlock,Transposeview(*tmpop),&big,WF,*tmpwf_,opQ,sqrt(2.0));
        
      //build the element
      double val = DotProduct(*tmpwf_,Op);
      sprintf(msg,"\nCheckWF: D(%i,%i,%i,%i)=%4.8e   Nleft=%i  NRight=%i",v_,t,v,u,val,tmpwf_->get_deltaQuantum(0).get_n(),Op.get_deltaQuantum(0).get_n());
      pout << msg;
        
    }

    //---------------
    //the right Block
    //---------------
    for (int i=0;i<rightBlock->get_op_array(CRE).size();i++){
      //get the operator
      boost::shared_ptr<SparseMatrix> tmpop = rightBlock->get_op_array(CRE).get_local_element(i)[0]->getworkingrepresentation(rightBlock);
      //get the orbital label
      int v_ = tmpop->get_orbs(0);

      //the second resulting wavefunction
      boost::shared_ptr<Wavefunction> tmpwf_(new Wavefunction(VQ[0],&big,true));
      operatorfunctions::TensorMultiply(rightBlock,Transposeview(*tmpop),&big,WF,*tmpwf_,opQ,sqrt(2.0));
        
      //build the element
      double val = DotProduct(*tmpwf_,Op);
      sprintf(msg,"\nCheckWF: D(%i,%i,%i,%i)=%4.8e",v_,t,v,u,val);
      pout << msg;
      
    }
  }
  
  
  //============================================================================
  // determine whether a given set of orbital indices is valid at this point in 
  // the sweep. This is to avoid double counting of certain index combinations
  //============================================================================
  bool CheckAllowedOperator(int t, int u, int v, int DotIndex, int iterCase, 
                            bool two_indexes_on_right){
    switch(iterCase){
      case _INITIAL_ITER_:
        return true;
        break;
      case _FINAL_ITER_:
        //make sure that one index is on the dot block
        if (t==DotIndex) return true;
        else if (u==DotIndex) return true;
        else if (v==DotIndex) return true;
        //this covers the 1 0 2 situation in the last iteration
        else if (two_indexes_on_right) return true;
        //if none of the other cases apply return a false
        else return false;
      case _REGULAR_ITER_:
        //make sure that one index is on the dot block
        if (t==DotIndex) return true;
        else if (u==DotIndex) return true;
        else if (v==DotIndex) return true;
        else return false;
    default:
      return false;
      
    } 
  }
  
  //============================================================================
  // Construct a+a+a where two indices are on the dot and the other on the
  // left block
  //============================================================================
  void ConstructCreCreDes_1_2_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CCD, int *OrbWin){
    char msg[512];
    //get the basename and stuff
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    int t,u,v,i,a,tu;
    int dummy=0;
    SpinBlock *leftBlock =  big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //open the operator array
    CCD.OpenFileWrite();
    //------------------------------------------------------
    // Case 1: a is on leftBlock and a+a+ is on the dotBlock
    //------------------------------------------------------
    for (tu=0;tu<dotBlock->get_op_array(CRE_CRE).get_size();tu++){
      //make sure we do not calculate the operator on multiple processors
      if (SkipOperator(big.get_sites(),dotBlock->get_sites()[0])) continue;
      //get the operator
      boost::shared_ptr<SparseMatrix> dotOpSing = dotBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(dotBlock);
      boost::shared_ptr<SparseMatrix> dotOpTrip = dotBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(dotBlock);
      //get the orbital labels
      t = dotOpSing->get_orbs(0);
      u = dotOpSing->get_orbs(1);
      //get the dot operator quanta
      SpinQuantum dOpQ = dotOpSing->get_deltaQuantum(0);
      
      for (int vcount=0;vcount<lBlock->get_op_array(CRE).get_size();vcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> leftOp = lBlock->get_op_array(CRE).get_local_element(vcount)[0]->getworkingrepresentation(lBlock);
        //get the orbital label
        v = leftOp->get_orbs(0);
        //get the left operator quanta
        SpinQuantum lopQ = Transposeview(*leftOp).get_deltaQuantum(0);
        //the possible quanta for the resulting operator
        SpinQuantum opQ = (dOpQ + lopQ)[0];
        //--------------
        //Construct Vtuv
        //--------------
        //Generate and construct the operator a+a+a (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        OpS.set_orbs() = dotOpSing->get_orbs();
        OpS.set_orbs().push_back(v);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = opQ;
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs() = dotOpSing->get_orbs();
        OpT.set_orbs().push_back(v);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = opQ;
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,Transposeview(*leftOp),*dotOpSing,leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,Transposeview(*leftOp),*dotOpTrip,leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operator
        boost::shared_ptr<Cre> Op (new Cre);
        Op->set_orbs().push_back(t);
        Op->set_orbs().push_back(u);
        Op->set_orbs().push_back(v);
        Op->set_initialised() = true;
        Op->set_fermion() = true;
        Op->resize_deltaQuantum(1);
        Op->set_deltaQuantum(0) = opQ;
        Op->allocate(leftBlock->get_stateInfo());

        //sum the two operators up
        ScaleAdd(-1.0/sqrt(2.0),OpS,*Op);
        ScaleAdd(-sqrt(3.0/2.0),OpT,*Op);//minus sign?
        //Store the operator
        CCD.AppendOperator(Op,t,u,v,false);
        
        //do the same thing for reversed operator indices
        if (t!=u){
          boost::shared_ptr<Cre> Op_ (new Cre);
          Op_->set_orbs().push_back(u);
          Op_->set_orbs().push_back(t);
          Op_->set_orbs().push_back(v);
          Op_->set_initialised() = true;
          Op_->set_fermion() = true;
          Op_->resize_deltaQuantum(1);
          Op_->set_deltaQuantum(0) = opQ;
          Op_->allocate(leftBlock->get_stateInfo());

          //sum the two operators up
          ScaleAdd(-1.0/sqrt(2.0),OpS,*Op);
          ScaleAdd(sqrt(3.0/2.0),OpT,*Op);//the additional minus sign arises from the change of orbital order
          
          //Store the operator
          CCD.AppendOperator(Op,u,t,v,false);
      
        }
        //clean up
        OpS.Clear();
        OpT.Clear();
      }//v
    }//tu
    
    //------------------------------------------------------
    // Case 2: a+ is on leftBlock and a+a is on the dotBlock
    //------------------------------------------------------
    for (int tv=0;tv<dotBlock->get_op_array(CRE_DES).get_size();tv++){
      //make sure we do not calculate the operator on multiple processors
      if (SkipOperator(big.get_sites(),dotBlock->get_sites()[0])) continue;
      //get the operator
      boost::shared_ptr<SparseMatrix> dotOpSing = dotBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(dotBlock);
      boost::shared_ptr<SparseMatrix> dotOpTrip = dotBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(dotBlock);
      //get the orbital labels
      t = dotOpSing->get_orbs(0);
      v = dotOpSing->get_orbs(1);
      //get the dot operator quanta
      SpinQuantum dOpQ = dotOpSing->get_deltaQuantum(0);
      
      for (int ucount=0;ucount<lBlock->get_op_array(CRE).get_size();ucount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> leftOp = lBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(lBlock);
        //get the orbital label
        u = leftOp->get_orbs(0);
        //get the left operator quanta
        SpinQuantum lopQ = leftOp->get_deltaQuantum(0);
        //the possible quanta for the resulting operator
        SpinQuantum opQ = (dOpQ + lopQ)[0];
        //--------------
        //Construct Vtuv
        //--------------
        //Generate and construct the operator a+a+a (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        Cre OpS_,OpT_;
        OpS.set_orbs() = dotOpSing->get_orbs();
        OpS.set_orbs().push_back(u);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = opQ;
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs() = dotOpSing->get_orbs();
        OpT.set_orbs().push_back(u);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = opQ;
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,*leftOp,*dotOpSing,leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,*leftOp,*dotOpTrip,leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operators
        boost::shared_ptr<Cre> Optuv (new Cre);
        Optuv->set_orbs().push_back(t);
        Optuv->set_orbs().push_back(u);
        Optuv->set_orbs().push_back(v);
        Optuv->set_initialised() = true;
        Optuv->set_fermion() = true;
        Optuv->resize_deltaQuantum(1);
        Optuv->set_deltaQuantum(0) = opQ;
        Optuv->allocate(leftBlock->get_stateInfo());

        boost::shared_ptr<Cre> Oputv (new Cre);
        Oputv->set_orbs().push_back(u);
        Oputv->set_orbs().push_back(t);
        Oputv->set_orbs().push_back(v);
        Oputv->set_initialised() = true;
        Oputv->set_fermion() = true;
        Oputv->resize_deltaQuantum(1);
        Oputv->set_deltaQuantum(0) = opQ;
        Oputv->allocate(leftBlock->get_stateInfo());

        //sum the two operators up
        ScaleAdd(-1.0/sqrt(2.0),OpS,*Optuv);
        ScaleAdd(-sqrt(3.0/2.0),OpT,*Optuv);
        ScaleAdd(sqrt(2.0),OpS,*Oputv);
        //Store the operator
        CCD.AppendOperator(Optuv,t,u,v,false);
        CCD.AppendOperator(Oputv,u,t,v,false);
        
        //clean up
        OpS.Clear();
        OpT.Clear();

        //do the same thing for the transposed operator
        if (t!=v){
          OpS_.set_orbs() = dotOpSing->get_orbs();
          OpS_.set_orbs().push_back(u);
          OpS_.set_initialised() = true;
          OpS_.set_fermion() = true;
          OpS_.resize_deltaQuantum(1);
          OpS_.set_deltaQuantum(0) = opQ;
          OpS_.allocate(leftBlock->get_stateInfo());
          OpT_.set_orbs() = dotOpSing->get_orbs();
          OpT_.set_orbs().push_back(u);
          OpT_.set_initialised() = true;
          OpT_.set_fermion() = true;
          OpT_.resize_deltaQuantum(1);
          OpT_.set_deltaQuantum(0) = opQ;
          OpT_.allocate(leftBlock->get_stateInfo());
          operatorfunctions::TensorProduct(lBlock,*leftOp,Transposeview(*dotOpSing),leftBlock,&(leftBlock->get_stateInfo()),OpS_,1.0);
          operatorfunctions::TensorProduct(lBlock,*leftOp,Transposeview(*dotOpTrip),leftBlock,&(leftBlock->get_stateInfo()),OpT_,1.0);
          
          //generate the product operators
          boost::shared_ptr<Cre> Opvut (new Cre);
          Opvut->set_orbs().push_back(v);
          Opvut->set_orbs().push_back(u);
          Opvut->set_orbs().push_back(t);
          Opvut->set_initialised() = true;
          Opvut->set_fermion() = true;
          Opvut->resize_deltaQuantum(1);
          Opvut->set_deltaQuantum(0) = opQ;
          Opvut->allocate(leftBlock->get_stateInfo());

          boost::shared_ptr<Cre> Opuvt (new Cre);
          Opuvt->set_orbs().push_back(u);
          Opuvt->set_orbs().push_back(v);
          Opuvt->set_orbs().push_back(t);
          Opuvt->set_initialised() = true;
          Opuvt->set_fermion() = true;
          Opuvt->resize_deltaQuantum(1);
          Opuvt->set_deltaQuantum(0) = opQ;
          Opuvt->allocate(leftBlock->get_stateInfo());
          
          //sum the two operators up
          ScaleAdd(-1.0/sqrt(2.0),OpS_,*Opvut);
          ScaleAdd(-sqrt(3.0/2.0),OpT_,*Opuvt);
          ScaleAdd(sqrt(2.0),OpS_,*Opuvt);
          //Store the operator
          CCD.AppendOperator(Opvut,v,u,t,false);
          CCD.AppendOperator(Opuvt,u,v,t,false);
    
          //clean up
          OpS_.Clear();
          OpT_.Clear();

        }
      }//u
    }//tv
    CCD.CloseFileWrite();
  }
  
  //============================================================================
  // Construct a+a+a where two indices are on the left and the other on the
  // dot block
  //============================================================================
  void ConstructCreCreDes_2_1_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CCD, int *OrbWin, Wavefunction &WF){
    char msg[512];
    //get the basename and stuff
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    int t,u,v,i,a,tu;
    int dummy=0;
    SpinBlock *leftBlock =  big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //open the operator array
    CCD.OpenFileWrite();
    //------------------------------------------------------
    // Case 1: a is on dotBlock and a+a+ is on the leftBlock
    //------------------------------------------------------
    for (tu=0;tu<lBlock->get_op_array(CRE_CRE).get_size();tu++){
      //get the operator
      boost::shared_ptr<SparseMatrix> lOpSing = lBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(lBlock);
      boost::shared_ptr<SparseMatrix> lOpTrip = lBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(lBlock);
      //get the orbital labels
      t = lOpSing->get_orbs(0);
      u = lOpSing->get_orbs(1);
      //get the left operator quanta
      SpinQuantum lOpQ = lOpSing->get_deltaQuantum(0);
      
      for (int vcount=0;vcount<dotBlock->get_op_array(CRE).get_size();vcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> dotOp = dotBlock->get_op_array(CRE).get_local_element(vcount)[0]->getworkingrepresentation(dotBlock);
        //get the orbital label
        v = dotOp->get_orbs(0);
        //get the dot operator quanta
        SpinQuantum dopQ = Transposeview(*dotOp).get_deltaQuantum(0);
        //the possible quanta for the resulting operator
        SpinQuantum opQ = (lOpQ + dopQ)[0];
        //--------------
        //Construct Vtuv
        //--------------
        //Generate and construct the operator a+a+a (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        OpS.set_orbs() = lOpSing->get_orbs();
        OpS.set_orbs().push_back(v);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = opQ;
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs() = lOpSing->get_orbs();
        OpT.set_orbs().push_back(v);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = opQ;
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,*lOpSing,Transposeview(*dotOp),leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,*lOpTrip,Transposeview(*dotOp),leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operator
        boost::shared_ptr<Cre> Op (new Cre);
        Op->set_orbs().push_back(t);
        Op->set_orbs().push_back(u);
        Op->set_orbs().push_back(v);
        Op->set_initialised() = true;
        Op->set_fermion() = true;
        Op->resize_deltaQuantum(1);
        Op->set_deltaQuantum(0) = opQ;
        Op->allocate(leftBlock->get_stateInfo());
        
        //sum the two operators up
        ScaleAdd(-1.0/sqrt(2.0),OpS,*Op);
        ScaleAdd(sqrt(3.0/2.0),OpT,*Op);//minus sign?

        //Store the operator
        CCD.AppendOperator(Op,t,u,v,false);
        //do the same thing for reversed operator indices
        if (t!=u){
          boost::shared_ptr<Cre> Op_ (new Cre);
          Op_->set_orbs().push_back(u);
          Op_->set_orbs().push_back(t);
          Op_->set_orbs().push_back(v);
          Op_->set_initialised() = true;
          Op_->set_fermion() = true;
          Op_->resize_deltaQuantum(1);
          Op_->set_deltaQuantum(0) = opQ;
          Op_->allocate(leftBlock->get_stateInfo());

          //sum the two operators up
          ScaleAdd(-1.0/sqrt(2.0),OpS,*Op_);
          ScaleAdd(-sqrt(3.0/2.0),OpT,*Op_);//the minus sign arises from the change of orbital order
          
          //Store the operator
          CCD.AppendOperator(Op_,u,t,v,false);
      
        }
        //clean up
        OpS.Clear();
        OpT.Clear();
      }//v
    }//tu
    
    //------------------------------------------------------
    // Case 2: a+ is on dotBlock and a+a is on the dotBlock
    //------------------------------------------------------
    for (int tv=0;tv<lBlock->get_op_array(CRE_DES).get_size();tv++){
      //get the operator
      boost::shared_ptr<SparseMatrix> lOpSing = lBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(lBlock);
      boost::shared_ptr<SparseMatrix> lOpTrip = lBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(lBlock);
      //get the orbital labels
      t = lOpSing->get_orbs(0);
      v = lOpSing->get_orbs(1);
      //get the dot operator quanta
      SpinQuantum lOpQ = lOpSing->get_deltaQuantum(0);
      
      for (int ucount=0;ucount<dotBlock->get_op_array(CRE).get_size();ucount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> dotOp = dotBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(dotBlock);
        //get the orbital label
        u = dotOp->get_orbs(0);
        //get the left operator quanta
        SpinQuantum dOpQ = dotOp->get_deltaQuantum(0);
        //the possible quanta for the resulting operator
        SpinQuantum opQ = (dOpQ + lOpQ)[0];
        //--------------
        //Construct Vtuv
        //--------------
        //Generate and construct the operator a+a+a (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        Cre OpS_,OpT_;
        OpS.set_orbs() = lOpSing->get_orbs();
        OpS.set_orbs().push_back(u);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = opQ;
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs() = lOpSing->get_orbs();
        OpT.set_orbs().push_back(u);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = opQ;
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,*lOpSing,*dotOp,leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,*lOpTrip,*dotOp,leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operators
        boost::shared_ptr<Cre> Optuv (new Cre);
        Optuv->set_orbs().push_back(t);
        Optuv->set_orbs().push_back(u);
        Optuv->set_orbs().push_back(v);
        Optuv->set_initialised() = true;
        Optuv->set_fermion() = true;
        Optuv->resize_deltaQuantum(1);
        Optuv->set_deltaQuantum(0) = opQ;
        Optuv->allocate(leftBlock->get_stateInfo());

        boost::shared_ptr<Cre> Oputv (new Cre);
        Oputv->set_orbs().push_back(u);
        Oputv->set_orbs().push_back(t);
        Oputv->set_orbs().push_back(v);
        Oputv->set_initialised() = true;
        Oputv->set_fermion() = true;
        Oputv->resize_deltaQuantum(1);
        Oputv->set_deltaQuantum(0) = opQ;
        Oputv->allocate(leftBlock->get_stateInfo());

        //sum the two operators up
        ScaleAdd(-1.0/sqrt(2.0),OpS,*Optuv);
        ScaleAdd(sqrt(3.0/2.0),OpT,*Optuv);
        ScaleAdd(sqrt(2.0),OpS,*Oputv);
        //Store the operator
        CCD.AppendOperator(Optuv,t,u,v,false);
        CCD.AppendOperator(Oputv,u,t,v,false);
        
        //clean up
        OpS.Clear();
        OpT.Clear();

        //do the same thing for the transposed operator
        if (t!=v){
          OpS_.set_orbs() = lOpSing->get_orbs();
          OpS_.set_orbs().push_back(u);
          OpS_.set_initialised() = true;
          OpS_.set_fermion() = true;
          OpS_.resize_deltaQuantum(1);
          OpS_.set_deltaQuantum(0) = opQ;
          OpS_.allocate(leftBlock->get_stateInfo());
          OpT_.set_orbs() = lOpSing->get_orbs();
          OpT_.set_orbs().push_back(u);
          OpT_.set_initialised() = true;
          OpT_.set_fermion() = true;
          OpT_.resize_deltaQuantum(1);
          OpT_.set_deltaQuantum(0) = opQ;
          OpT_.allocate(leftBlock->get_stateInfo());
          operatorfunctions::TensorProduct(lBlock,Transposeview(*lOpSing),*dotOp,leftBlock,&(leftBlock->get_stateInfo()),OpS_,1.0);
          operatorfunctions::TensorProduct(lBlock,Transposeview(*lOpTrip),*dotOp,leftBlock,&(leftBlock->get_stateInfo()),OpT_,1.0);
          
          //generate the product operators
          boost::shared_ptr<Cre> Opvut (new Cre);
          Opvut->set_orbs().push_back(v);
          Opvut->set_orbs().push_back(u);
          Opvut->set_orbs().push_back(t);
          Opvut->set_initialised() = true;
          Opvut->set_fermion() = true;
          Opvut->resize_deltaQuantum(1);
          Opvut->set_deltaQuantum(0) = opQ;
          Opvut->allocate(leftBlock->get_stateInfo());

          boost::shared_ptr<Cre> Opuvt (new Cre);
          Opuvt->set_orbs().push_back(u);
          Opuvt->set_orbs().push_back(v);
          Opuvt->set_orbs().push_back(t);
          Opuvt->set_initialised() = true;
          Opuvt->set_fermion() = true;
          Opuvt->resize_deltaQuantum(1);
          Opuvt->set_deltaQuantum(0) = opQ;
          Opuvt->allocate(leftBlock->get_stateInfo());
          
          //sum the two operators up
          ScaleAdd(-1.0/sqrt(2.0),OpS_,*Opvut);
          ScaleAdd(-sqrt(3.0/2.0),OpT_,*Opvut);//account for transposition and add a minus sign
          ScaleAdd(sqrt(2.0),OpS_,*Opuvt);
          
          //Store the operator
          CCD.AppendOperator(Opvut,v,u,t,false);
          CCD.AppendOperator(Opuvt,u,v,t,false);

          //clean up
          OpS_.Clear();
          OpT_.Clear();
          
        }
      }//u
    }//tv
    CCD.CloseFileWrite();
    
  }
  
  
  //============================================================================
  // Construct a+a+a on a single block
  //============================================================================
  void ConstructCreCreDesSingleBlock(SpinBlock &big, SpinBlock *block, ThreeIndOpArray &CCD,
                                     IntegralContainer &IKJL, int iterCase, Wavefunction &WF){
    char msg[512];
    char BaseName[512];
    int OrbWin[6];
    double ENuc=0.0;
    bool ConvOverlap;
    //get the basename and stuff
    ReadInput(BaseName,OrbWin,ENuc,ConvOverlap);
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = NInternal+NActive+NExternal;
    int t,u,v,i,a,tu;
    int dummy;
    SpinBlock *leftBlock =  big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();

    //open the integral container
    IKJL.OpenFileRead();
    //open the operator array
    CCD.OpenFileWrite();

    for (int uv=0;uv<block->get_op_array(CRE_DES).get_size();uv++){
      //if this is a regular or the last iteration make sure we do not calculate 
      //the operator on multiple processors
      if (iterCase!=_INITIAL_ITER_ && SkipOperator(big.get_sites(),block->get_sites()[0]))continue;
      //get the first operator
      boost::shared_ptr<SparseMatrix> Euv = block->get_op_array(CRE_DES).get_local_element(uv)[0]->getworkingrepresentation(block);
      //get the orbital labels
      u=Euv->get_orbs(0);
      v=Euv->get_orbs(1);
      //get the operator quanta
      SpinQuantum opQ1 = Euv->get_deltaQuantum(0);
      for (int tcount=0;tcount<block->get_op_array(CRE).get_size();tcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> at = block->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(block);
        //get the orbital label
        t = at->get_orbs(0);
        //get the operator quanta
        SpinQuantum opQ2 = at->get_deltaQuantum(0);
        SpinQuantum opQ = (opQ1 + opQ2)[0];
        //get the integrals
        boost::shared_ptr<Matrix> Kut = IKJL.GetMatrix(u+NInternal,t+NInternal);
        
        //---------------
        //Construct Vtuv
        //---------------
        //the resulting operator
        boost::shared_ptr<Cre> Op(new Cre);
        Op->set_orbs().push_back(t);
        Op->set_orbs().push_back(u);
        Op->set_orbs().push_back(v);
        Op->set_initialised()=true;
        Op->set_fermion()=true;
        Op->resize_deltaQuantum(1);
        Op->set_deltaQuantum(0)=opQ;
        Op->allocate(block->get_stateInfo());
        //multiply the two operators
        operatorfunctions::Product(block,*at,*Euv,*Op,sqrt(2.0));
        boost::shared_ptr<Cre> OpRep(new Cre);

        switch(iterCase){
          case _INITIAL_ITER_:
            //store the operator
            CCD.AppendOperator(Op,t,u,v);
            break;
          case _REGULAR_ITER_:
          case _FINAL_ITER_:
            //build the operator representation in the larger block
            OpRep->set_orbs() = Op->get_orbs();
            OpRep->set_initialised()=true;
            OpRep->resize_deltaQuantum(1);
            OpRep->set_deltaQuantum(0)=opQ;
            OpRep->allocate(leftBlock->get_stateInfo());
            operatorfunctions::TensorTrace(block,*Op,leftBlock,&(leftBlock->get_stateInfo()),*OpRep);
            //store the operator
            CCD.AppendOperator(OpRep,t,u,v);
            break;
        }//itercase
        //------------------------
        //do the same for the Vtvu
        //------------------------
        if (u!=v){
          //the resulting operator
          boost::shared_ptr<Cre> Op_(new Cre);
          Op_->set_orbs().push_back(t);
          Op_->set_orbs().push_back(v);
          Op_->set_orbs().push_back(u);
          Op_->set_initialised()=true;
          Op_->set_fermion()=true;
          Op_->resize_deltaQuantum(1);
          Op_->set_deltaQuantum(0)=opQ;
          Op_->allocate(block->get_stateInfo());
          //multiply the two operators
          operatorfunctions::Product(block,*at,Transposeview(*Euv),*Op_,sqrt(2.0));
          boost::shared_ptr<Cre> OpRep_(new Cre);

          switch (iterCase){
            case _INITIAL_ITER_:
              //store the operator
              CCD.AppendOperator(Op_,t,v,u);
              break;
            case _REGULAR_ITER_:
            case _FINAL_ITER_:
              //build the operator representation in the larger block
              OpRep_->set_orbs() = Op_->get_orbs();
              OpRep_->set_initialised()=true;
              OpRep_->resize_deltaQuantum(1);
              OpRep_->set_deltaQuantum(0)=opQ;
              OpRep_->allocate(leftBlock->get_stateInfo());
              operatorfunctions::TensorTrace(block,*Op_,leftBlock,&(leftBlock->get_stateInfo()),*OpRep_);
              //store the operator
              CCD.AppendOperator(OpRep_,t,v,u);
              break;
          }//itercase
        }//u!=v
      }//t
    }//uv
    //close the integral container
    IKJL.CloseFileRead();
    //close the array
    CCD.CloseFileWrite();
    
  }
  
  //============================================================================
  //Construct the contributions to the perturber functions V(i) where all indices
  //are on one block
  //============================================================================
  void ConstructViSingleBlock(SpinBlock &big, Wavefunction &WF, ThreeIndOpArray &CCD,
                              ThreeIndOpArray CCD_, vector<WavefunctionArray> &Ti, 
                              IntegralContainer &IKJL, NEVPT2Info &Info){
    
    char msg[512];
    //get the orbital spaces and stuff
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
    int t=1,u=-1,v=-1,i,a,tu;
    int t_,u_,v_;
    int dummy;
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    double twoS = (double) WFQ.get_s().getirrep();
    double S    = twoS/2.0;
    double twoS_,S_,fac;
    char BaseName[512];
    Info.getBaseName(BaseName);

    //---------------------------------------------------
    //first construct all contributions on the left Block
    //---------------------------------------------------
    //open the integral file
    IKJL.OpenFileRead();
    boost::shared_ptr<Matrix> Kut = IKJL.GetMatrix(0,0);
    //Open the operator arrays and prepare the buffer
    CCD.OpenFileRead();
    CCD.ResetBuffer();
    //open the V(i) buffer
    for (i=0;i<Ti.size();i++){
      Ti[i].OpenFileRead();
    }
    //initialize the auxiliary indices
    t_ = -1;
    u_ = -1;
    v_ = -1;
    //start the loop
    bool EndOfOuterArray = false;
    bool NeedTensorTrace = false;
    //make sure that all processes undergo the same number of cycles in the loop
    LoopSynchronizer LS(CCD.GetLength());//this might lead to trouble if we have stored the same operator twice;
    //get the dimension of the loop
    int Dim_tuv = LS.GetGlobalDim();
    int ituv=0;
    while (ituv<Dim_tuv){
      //determine if this is just a dummy iteration
      bool DummyIter = LS.DummyIter(ituv);
      //get the operator
      boost::shared_ptr<Cre> Otuv;
      if (!DummyIter) {
        Otuv = CCD.GetOpFromBuffer(t,u,v,NeedTensorTrace,EndOfOuterArray);
      }
      else{
        //if this is a dummy iteration, generate an empty dummy operator
        NeedTensorTrace=false;
        Otuv = boost::make_shared<Cre> ();
        Otuv->set_initialised() = true;
        Otuv->resize_deltaQuantum(1);
        Otuv->set_deltaQuantum(0) = SpinQuantum(1,SpinSpace(1),IrrepSpace(0));
      }
      //if we don't have an operator and this is not a dummy iteration, 
      //then something went wrong
      if (EndOfOuterArray&&!DummyIter){
        sprintf(msg,"\nERROR NEVPT2: in V(i): loop over (t,u,v) has not been properly synchronized!!!");
        pout << msg;
        DummyIter = true;
      }
      
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
      if (((t!=t_)||(u!=u_))  &&  !DummyIter){
         Kut = IKJL.GetMatrix(u+NInternal,t+NInternal);
      }
      //the possible quanta
      SpinQuantum opQ = Otuv->get_deltaQuantum(0);
      vector<SpinQuantum> VQ = WFQ + opQ;

      for (int iquanta=0;iquanta<VQ.size();iquanta++){

        //evaluate the prefactor
        twoS_ = (double) VQ[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt((2.0*S_+1)/(2.0*S+1));

        //generate V(tuv) = O(tuv)|psi>
        boost::shared_ptr<Wavefunction> Vtuv(new Wavefunction(VQ[iquanta],&big,true));
        if (!DummyIter) operatorfunctions::TensorMultiply(leftBlock,*Otuv,&big,WF,*Vtuv,opQ,fac);        

        //add it to the V(i) functions
        Ti[iquanta].ResetBuffer();
        bool EndOfInnerArray = false;
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Vi = Ti[iquanta].GetOpFromBuffer(dummy,i,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(Kut->element(v+NInternal,i),*Vtuv,*Vi);
            Ti[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
          }
        }//i
        //Send Vtuv around to other processes and receive theirs
        SendAroundVtuv_i(*Vtuv,*Kut,Ti[iquanta],v+NInternal,DummyIter);
      }//iquanta
      t_ = t;
      u_ = u;
      v_ = v;
      ituv++;
    }//tuv
    //close the Operator array
    CCD.CloseFileRead();
    
    //---------------------------------------------------
    //then construct all contributions on the right Block
    //---------------------------------------------------
    //Open the operator array and prepare the buffer
    CCD_.OpenFileRead();
    CCD_.ResetBuffer();
    //initialize the auxiliary indices
    t_ = -1;
    u_ = -1;
    v_ = -1;
    //start the loop
    EndOfOuterArray = false;
    //make sure that all processes undergo the same number of cycles in the loop
    LoopSynchronizer LS_(CCD_.GetLength());//this might lead to trouble if we have stored the same operator twice;
   
    //get the dimension of the loop
    Dim_tuv = LS_.GetGlobalDim();
    ituv=0;
    while (ituv<Dim_tuv){
      //determine if this is just a dummy iteration
      bool DummyIter = LS_.DummyIter(ituv);
    
      //get the operator
      boost::shared_ptr<Cre> Otuv;
      if (!DummyIter) {
        Otuv = CCD_.GetOpFromBuffer(t,u,v,NeedTensorTrace,EndOfOuterArray);
      }
      else{
        //if this is a dummy iteration, generate an empty dummy operator
        NeedTensorTrace=false;
        Otuv = boost::make_shared<Cre> ();
        Otuv->set_initialised() = true;
        Otuv->resize_deltaQuantum(1);
        Otuv->set_deltaQuantum(0) = SpinQuantum(1,SpinSpace(1),IrrepSpace(0));
      }
      //if we don't have an operator and this is not a dummy iteration, 
      //then something went wrong
      if (EndOfOuterArray&&!DummyIter){
        sprintf(msg,"\nERROR NEVPT2: in V(i): loop over (t,u,v) has not been properly synchronized!!!");
        mpi_message(msg);
        DummyIter = true;
      }

      //if necessary, get the integral matrix
      if (((t!=t_)||(u!=u_))   && !DummyIter){
         Kut = IKJL.GetMatrix(u+NInternal,t+NInternal);
      }
      //the possible quanta
      SpinQuantum opQ = Otuv->get_deltaQuantum(0);
      vector<SpinQuantum> VQ = WFQ + opQ;

      for (int iquanta=0;iquanta<VQ.size();iquanta++){
        //evaluate the prefactor
        twoS_ = (double) VQ[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt((2.0*S_+1)/(2.0*S+1));

        //generate V(tuv) = O(tuv)|psi>
        boost::shared_ptr<Wavefunction> Vtuv(new Wavefunction(VQ[iquanta],&big,true));
        SpinQuantum opQ = Vtuv->get_deltaQuantum(0);
        if (!DummyIter) operatorfunctions::TensorMultiply(rightBlock,*Otuv,&big,WF,*Vtuv,opQ,fac);        

        //add it to the V(i) functions
        Ti[iquanta].ResetBuffer();
        bool EndOfInnerArray = false;
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Vi = Ti[iquanta].GetOpFromBuffer(dummy,i,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(Kut->element(v+NInternal,i),*Vtuv,*Vi);
            Ti[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
          }
        }//i
        //Send Vtuv around to other processes and receive theirs
        SendAroundVtuv_i(*Vtuv,*Kut,Ti[iquanta],v+NInternal,DummyIter);
      }//iquanta
      t_ = t;
      u_ = u;
      v_ = v;
      ituv++;
    }//tuv
    //close the Operator array
    CCD_.CloseFileRead();
    //close the integral container
    IKJL.CloseFileRead();
    //close the V(i) array
    for (i=0;i<Ti.size();i++){
      Ti[i].CloseFileRead();
    }

  }

  
  
  //============================================================================
  // Construct a+a+a|psi> on two blocks block
  //============================================================================
  void ConstructCreCreDesTwoBlocks(SpinBlock &big, Wavefunction &WF, vector<WavefunctionArray> &T, 
                                 IntegralContainer &IKJL,int iterCase, int *OrbWin){
    char msg[512];
    char BaseName[512];
    //get the orbital spaces and stuff
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
    int t,u,v,i,a,tu;
    int dummy;
    
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *dotBlock = leftBlock->get_rightBlock();
    int dotIndex = dotBlock->get_sites()[0];//let's hope that this Block has only one member......
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    double twoS = (double) WFQ.get_s().getirrep();
    double S    = twoS/2.0;
    double twoS_,S_,fac;
    int NumProcs = mpi_world_size();

    //Note: we have to synchronize the loops over (t,u,v) and broadcast the
    //operators, because each processor needs every possible operator O(tuv).
    //The division between processes is made in index a
#ifndef SERIAL
    mpi::communicator world;
#endif


    //open the integral and operator container
    IKJL.OpenFileRead();
    for (int i=0;i<T.size();i++){
      T[i].OpenFileRead();
    }
    //--------------------------------------------------------------------------
    //Case 1: two creation operators on one block and the annihilation operator 
    //        on the other block
    //--------------------------------------------------------------------------
    //make sure that allprocesses undergo the same number of cycles in the loop
    LoopSynchronizer LS(leftBlock->get_op_array(CRE_CRE).get_size());
    int Dim_tu=LS.GetGlobalDim();
    //loop over a+a+ on the left block
    //for (tu=0;tu<leftBlock->get_op_array(CRE_CRE).get_size();tu++){
    for(int tu=0;tu<Dim_tu;tu++){
      bool DummyIter = LS.DummyIter(tu);
      //get the left operator
      boost::shared_ptr<SparseMatrix> LeftOpSing;
      boost::shared_ptr<SparseMatrix> LeftOpTrip;
      if (!DummyIter){
        LeftOpSing = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(leftBlock);
        LeftOpTrip = leftBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(leftBlock);
      }
      else{
        //if this is a dummy iteration generate dummy operators
        LeftOpSing = boost::make_shared<Cre>();
        LeftOpTrip = boost::make_shared<Cre>();
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpSing->resize_deltaQuantum(1);
        LeftOpTrip->resize_deltaQuantum(1);
        LeftOpSing->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(0),IrrepSpace(0));
        LeftOpTrip->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(2),IrrepSpace(0));
      }
      //get the orbital indices
      t = LeftOpSing->get_orbs(0);
      u = LeftOpSing->get_orbs(1);
      //the left operator quanta
      SpinQuantum lopQS = LeftOpSing->get_deltaQuantum(0);
      SpinQuantum lopQT = LeftOpTrip->get_deltaQuantum(0);
      //get the integrals
      boost::shared_ptr<Matrix> Ktu = IKJL.GetMatrix(t+NInternal,u+NInternal);
      //loop over a on the right block
      for (int vcount=0;vcount<rightBlock->get_op_array(CRE).get_size();vcount++){
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOp = rightBlock->get_op_array(CRE).get_local_element(vcount)[0]->getworkingrepresentation(rightBlock);
        //get the orbital index
        v = RightOp->get_orbs(0);
        //the right operator quanta
        SpinQuantum ropQ = Transposeview(*RightOp).get_deltaQuantum(0);
        //the total Operator Quanta (is the same for the left singlet and triplet operator)
        vector<SpinQuantum> VopQS = lopQS + ropQ;
        //the possible quanta for the resulting wavefunction
        vector <SpinQuantum> VQ = WFQ + VopQS[0];
        //evaluate the prefactors
        vector<double> TripFac;
        vector<double> SingFac;
        for (i=0;i<VQ.size();i++){
          twoS_ = (double)VQ[i].get_s().getirrep();
          S_    = twoS_/2.0;
          fac   = -sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
          SingFac.push_back(fac);
          TripFac.push_back(-fac*sqrt(3.0));
        }//resulting wavefunction quanta
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunctions that holds the result of the multiplications
          boost::shared_ptr<Wavefunction> VtuvS (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VtuvT (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> Vtuv (new Wavefunction(VQ[iquanta],&big,true));
          if (!DummyIter){
            //multiply the Wavefunction with the operators
            operatorfunctions::TensorMultiply(leftBlock,*LeftOpSing,Transposeview(*RightOp),&big,WF,*VtuvS,VopQS[0],SingFac[iquanta]);
            operatorfunctions::TensorMultiply(leftBlock,*LeftOpTrip,Transposeview(*RightOp),&big,WF,*VtuvT,VopQS[0],TripFac[iquanta]);
            //sum up the contributions
            ScaleAdd(1.0,*VtuvS,*Vtuv);
            ScaleAdd(1.0,*VtuvT,*Vtuv);
          }
          //do the same thing for Vutv
          boost::shared_ptr<Wavefunction> Vutv (new Wavefunction(VQ[iquanta],&big,true));
          if (t!=u&&!DummyIter){
            ScaleAdd(1.0,*VtuvS,*Vutv);
            ScaleAdd(-1.0,*VtuvT,*Vutv);//the minus sign arises from the change of orbital order
          }//t!=u
          if (!T[iquanta].GetBasisOnly()){
            //store the operator
            bool EndOfArray = false;
            T[iquanta].ResetBuffer();
            while (!EndOfArray){
              boost::shared_ptr<Wavefunction> Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
              if (!EndOfArray){
                //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                /*
                if (fabs(Ktu->element(i,v+NInternal))>1.e-16){
                  sprintf(msg,"\ncase 1: (%i%ii|%i%i) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock left a",t,i,u,v,t,u,v,iquanta,SingFac[iquanta],TripFac[iquanta],Ktu->element(i,v+NInternal),DotProduct(*Vtuv,*Vtuv));pout << msg;
                }
                if ((t!=u)&&(fabs(Ktu->element(v+NInternal,i))>1.e-16)){
                  sprintf(msg,"\ncase 1: (%i%i|%i%ii) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock left b",t,v,u,i,u,t,v,iquanta,SingFac[iquanta],TripFac[iquanta],Ktu->element(v+NInternal,i),DotProduct(*Vutv,*Vutv));pout << msg;
                }
                 */
                ScaleAdd(Ktu->element(i,v+NInternal),*Vtuv,*Vi);
                ScaleAdd(Ktu->element(v+NInternal,i),*Vutv,*Vi);
                T[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
              }
            }//i
            //Send Vtuv around to other processes and receive theirs
            SendAroundVtuv_i(*Vtuv,*Vutv,*Ktu,T[iquanta],v+NInternal,DummyIter);
          }//!BasisOnly
          T[iquanta].AddToBasisOperator(Vtuv);
          T[iquanta].AddToBasisOperator(Vutv);
        }//iquanta
      }//v
    }//tu
    LS.Reset();
    if (iterCase==_FINAL_ITER_){
      LS.SetLocalDim(rightBlock->get_op_array(CRE_CRE).get_size());
      Dim_tu = LS.GetGlobalDim();
      //loop over a+a+ on the right block
      //for (tu=0;tu<rightBlock->get_op_array(CRE_CRE).get_size();tu++){
      for (int tu=0;tu<Dim_tu;tu++){
        //determine whether this is a dummy iteration
        bool DummyIter = LS.DummyIter(tu);
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOpSing;
        boost::shared_ptr<SparseMatrix> RightOpTrip;
        if (!DummyIter){
          RightOpSing = rightBlock->get_op_array(CRE_CRE).get_local_element(tu)[0]->getworkingrepresentation(rightBlock);
          RightOpTrip = rightBlock->get_op_array(CRE_CRE).get_local_element(tu)[1]->getworkingrepresentation(rightBlock);
        }
        else{
          //if this is a dummy iteration generate dummy operators
          RightOpSing = boost::make_shared<Cre>();
          RightOpTrip = boost::make_shared<Cre>();
          RightOpSing->set_orbs().push_back(-1);
          RightOpSing->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpSing->resize_deltaQuantum(1);
          RightOpTrip->resize_deltaQuantum(1);
          RightOpSing->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(0),IrrepSpace(0));
          RightOpTrip->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(2),IrrepSpace(0));
        }
        //get the orbital indices
        t = RightOpSing->get_orbs(0);
        u = RightOpSing->get_orbs(1);
        //the right operator quanta
        SpinQuantum ropQS = RightOpSing->get_deltaQuantum(0);
        SpinQuantum ropQT = RightOpTrip->get_deltaQuantum(0);
        //get the integrals
        boost::shared_ptr<Matrix> Ktu = IKJL.GetMatrix(t+NInternal,u+NInternal);
        //loop over a on the left block
        for (int vcount=0;vcount<leftBlock->get_op_array(CRE).get_size();vcount++){
          //get the left operator
          boost::shared_ptr<SparseMatrix> LeftOp = leftBlock->get_op_array(CRE).get_local_element(vcount)[0]->getworkingrepresentation(leftBlock);
          //get the orbital index
          v = LeftOp->get_orbs(0);
          //check orbital indices for validity
          //if (!CheckAllowedOperator(t,u,v,dotIndex,iterCase,true)) continue;
          //the right operator quanta
          SpinQuantum lopQ = Transposeview(*LeftOp).get_deltaQuantum(0);
          //the total Operator Quanta (is the same for the left singlet and triplet operator)
          vector<SpinQuantum> VopQS = ropQS + lopQ;
          //the possible quanta for the resulting wavefunction
          vector <SpinQuantum> VQ = WFQ + VopQS[0];
          //evaluate the prefactors
          vector<double> TripFac;
          vector<double> SingFac;
          for (i=0;i<VQ.size();i++){
            twoS_ = (double)VQ[i].get_s().getirrep();
            S_    = twoS_/2.0;
            fac   = -sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
            SingFac.push_back(fac);
            ///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // NOTE:
            // for some strange reason I dont't understand we have to add a minus
            // sign if the triplet A operator is on the right Block. 
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            TripFac.push_back(fac*sqrt(3.0));
          }//resulting wavefunction quanta
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunctions that holds the result of the multiplications
            boost::shared_ptr<Wavefunction> VtuvS (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtuvT (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> Vtuv (new Wavefunction(VQ[iquanta],&big,true));
            //multiply the Wavefunction with the operators
            if (!DummyIter){
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOp),*RightOpSing,&big,WF,*VtuvS,VopQS[0],SingFac[iquanta]);
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOp),*RightOpTrip,&big,WF,*VtuvT,VopQS[0],TripFac[iquanta]);
              //sum up the contributions
              ScaleAdd(1.0,*VtuvS,*Vtuv);
              ScaleAdd(1.0,*VtuvT,*Vtuv);
            }
            //do the same thing for Vutv
            boost::shared_ptr<Wavefunction> Vutv (new Wavefunction(VQ[iquanta],&big,true));
            if (t!=u&&!DummyIter){
              ScaleAdd(1.0,*VtuvS,*Vutv);
              ScaleAdd(-1.0,*VtuvT,*Vutv);//the minus sign arises from the change of orbital order
            }//t!=u
            if (!T[iquanta].GetBasisOnly()){
              //store the operator
              T[iquanta].ResetBuffer();
              bool EndOfArray = false;
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
                if (!EndOfArray){
                  /*
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                  if (fabs(Ktu->element(i,v+NInternal))>1.e-18){
                    sprintf(msg,"\ncase 1: (%i%ii|%i%i) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock right a",t,i,u,v,t,u,v,iquanta,SingFac[iquanta],TripFac[iquanta],Ktu->element(i,v+NInternal),DotProduct(*Vtuv,*Vtuv));pout << msg;
                  }
                  if ((t!=u)&&(fabs(Ktu->element(v+NInternal,i))>1.e-18)){
                    sprintf(msg,"\ncase 1: (%i%i|%i%ii) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock right b",t,v,u,i,u,t,v,iquanta,SingFac[iquanta],TripFac[iquanta],Ktu->element(v+NInternal,i),DotProduct(*Vutv,*Vutv));pout << msg;
                  }
                  */
                  ScaleAdd(Ktu->element(i,v+NInternal),*Vtuv,*Vi);
                  ScaleAdd(Ktu->element(v+NInternal,i),*Vutv,*Vi);
                  T[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
                }
              }//i//store the operator
              //Send Vtuv around to other processes and receive theirs
              SendAroundVtuv_i(*Vtuv,*Vutv,*Ktu,T[iquanta],v+NInternal,DummyIter);
            }//!BasisOnly
            T[iquanta].AddToBasisOperator(Vtuv);
            T[iquanta].AddToBasisOperator(Vutv);
          }//iquanta
        }//v
      }//tu
      LS.Reset();
    }//final iteration
    
    //--------------------------------------------------------------------------
    //Case 2: a creation and an annihilation operator are on one block and the 
    //        the second creation operator is on the other block
    //--------------------------------------------------------------------------
    double fac_;
    //make sure that allprocesses undergo the same number of cycles in the loop
    LS.SetLocalDim(leftBlock->get_op_array(CRE_DES).get_size());
    int Dim_tv = LS.GetGlobalDim();
    //loop over a+a on left Block
    //for (int tv=0;tv<leftBlock->get_op_array(CRE_DES).get_size();tv++){
    for (int tv=0;tv<Dim_tv;tv++){
      //determine whether this is a dummy iteration
      bool DummyIter = LS.DummyIter(tv);
      //get the left operator
      boost::shared_ptr<SparseMatrix> LeftOpSing;
      boost::shared_ptr<SparseMatrix> LeftOpTrip;
      if (!DummyIter){
        LeftOpSing = leftBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(leftBlock);
        LeftOpTrip = leftBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(leftBlock);
      }
      else{
        //if this is a dummy iteration generate dummy operators
        LeftOpSing = boost::make_shared<Cre>();
        LeftOpTrip = boost::make_shared<Cre>();
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpSing->resize_deltaQuantum(1);
        LeftOpTrip->resize_deltaQuantum(1);
        LeftOpSing->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(0),IrrepSpace(0));
        LeftOpTrip->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(2),IrrepSpace(0));
      }
      //get the orbital indices
      t = LeftOpSing->get_orbs(0);
      v = LeftOpSing->get_orbs(1);
      //the left operator quanta
      SpinQuantum lopQS = LeftOpSing->get_deltaQuantum(0);
      SpinQuantum lopQT = LeftOpTrip->get_deltaQuantum(0);
      //get the integrals
      boost::shared_ptr<Matrix> Ktv = IKJL.GetMatrix(t+NInternal,v+NInternal);
      //loop over a+ on the right block
      for (int ucount=0;ucount<rightBlock->get_op_array(CRE).get_size();ucount++){
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOp = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
        //get the orbital index
        u = RightOp->get_orbs(0);
        //check orbital indices for validity
        //if (!CheckAllowedOperator(t,u,v,dotIndex,iterCase)) continue;
        //the right operator quanta
        SpinQuantum ropQ = RightOp->get_deltaQuantum(0);
        //the total Operator Quanta (is the same for the left singlet and triplet operator)
        vector<SpinQuantum> VopQS = lopQS + ropQ;
        //the possible quanta for the resulting wavefunction
        vector <SpinQuantum> VQ = WFQ + VopQS[0];
        //evaluate the prefactors
        vector<double> TripFac;
        vector<double> SingFac;
        vector<double> SingFac_;
        for (i=0;i<VQ.size();i++){
          twoS_ = (double)VQ[i].get_s().getirrep();
          S_    = twoS_/2.0;
          fac   = sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));\
          fac_ = sqrt(2.0)*sqrt(2.0*S_+1.0)/sqrt(2*S+1.0);//this second term arises from case 3
          SingFac.push_back(-fac);
          SingFac_.push_back(fac_);
          TripFac.push_back(fac*sqrt(3.0));
        }//resulting wavefunction quanta
        //get the integrals for case 3
        boost::shared_ptr<Matrix> Kut = IKJL.GetMatrix(u+NInternal,t+NInternal);
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //----------------------------
          //construct the operators Vtvu
          //----------------------------
          //generate the wavefunctions that holds the result of the multiplications
          boost::shared_ptr<Wavefunction> VtvuS (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VtvuT (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> Vtvu (new Wavefunction(VQ[iquanta],&big,true));
          //multiply the Wavefunction with the operators
          if (!DummyIter){
            operatorfunctions::TensorMultiply(leftBlock,*LeftOpSing,*RightOp,&big,WF,*VtvuS,VopQS[0],1.0);
            operatorfunctions::TensorMultiply(leftBlock,*LeftOpTrip,*RightOp,&big,WF,*VtvuT,VopQS[0],1.0);
            //sum up the contributions
            ScaleAdd(SingFac[iquanta],*VtvuS,*Vtvu);
            ScaleAdd(TripFac[iquanta],*VtvuT,*Vtvu);
          }
          //--------------------------
          //do the same thing for Vvtu
          //--------------------------
          boost::shared_ptr<Wavefunction> Vvtu (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VvtuS (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VvtuT (new Wavefunction(VQ[iquanta],&big,true));
          if (t!=v){
            //multiply the transposed operator
            if (!DummyIter){
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOpSing),*RightOp,&big,WF,*VvtuS,VopQS[0],1.0);
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOpTrip),*RightOp,&big,WF,*VvtuT,VopQS[0],1.0);
              //sum up the contributions
              ScaleAdd(SingFac[iquanta],*VvtuS,*Vvtu);
              ScaleAdd(-TripFac[iquanta],*VvtuT,*Vvtu);//account for transposition and add a minus sign
            }
          }//t!=u
          if (!T[iquanta].GetBasisOnly()){
            //store the operators
            T[iquanta].ResetBuffer();
            bool EndOfArray = false;
            while (!EndOfArray){
              boost::shared_ptr<Wavefunction> Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
              if (!EndOfArray){
                /*
                if (fabs(Ktv->element(i,u+NInternal))>1.e-18){
                  sprintf(msg,"\ncase 2: (%i%ii|%i%i) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock left a",t,i,v,u,t,u,v,iquanta,SingFac[iquanta],TripFac[iquanta],Ktv->element(i,u+NInternal),DotProduct(*Vtvu,*Vtvu));pout << msg;
                }
                if ((t!=v)&&(fabs(Ktv->element(u+NInternal,i))>1.e-18)){
                  sprintf(msg,"\ncase 2: (%i%i|%i%ii) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock left b",t,u,v,i,v,u,t,iquanta,SingFac[iquanta],TripFac[iquanta],Ktv->element(u+NInternal,i),DotProduct(*Vvtu,*Vvtu));pout << msg;
                }
                if (fabs(Kut->element(i,v+NInternal))>1.e-18){
                  sprintf(msg,"\ncase 3: (%i%ii|%i%i) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e int=%lf Norm=%lf twoBlock left a",u,i,t,v,u,t,v,iquanta,SingFac_[iquanta],Kut->element(i,v+NInternal),4.0*DotProduct(*VtvuS,*VtvuS));pout << msg;
                }
                if ((t!=v)&&(fabs(Kut->element(i,v+NInternal))>1.e-18)){
                  sprintf(msg,"\ncase 3: (%i%ii|%i%i) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e int=%lf Norm=%lf twoBlock left b",u,i,t,v,u,v,t,iquanta,SingFac_[iquanta],Kut->element(i,v+NInternal),4.0*DotProduct(*VvtuS,*VvtuS));pout << msg;
                }
                 */
                ScaleAdd(Ktv->element(i,u+NInternal),*Vtvu,*Vi);
                ScaleAdd(Kut->element(i,v+NInternal)*SingFac_[iquanta],*VtvuS,*Vi);//case 3
                ScaleAdd(Ktv->element(u+NInternal,i),*Vvtu,*Vi);
                ScaleAdd(Kut->element(i,v+NInternal)*SingFac_[iquanta],*VvtuS,*Vi);//case 3
                T[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
              }
            }//i
            //Send Vtuv around to other processes and receive theirs
            SendAroundVtuv_i(*Vtvu,*VtvuS,*Vvtu,*VvtuS,*Ktv,*Kut,T[iquanta],u+NInternal,v+NInternal,DummyIter,SingFac_[iquanta]);
          }//!BasisOnly
          T[iquanta].AddToBasisOperator(Vtvu);
          T[iquanta].AddToBasisOperator(VtvuS,SingFac_[iquanta]);
          T[iquanta].AddToBasisOperator(Vvtu);
          T[iquanta].AddToBasisOperator(VvtuS,SingFac_[iquanta]);
        }//iquanta
      }//u
    }//tv
    LS.Reset();
    if (iterCase==_FINAL_ITER_){
      //make sure that allprocesses undergo the same number of cycles in the loop
      LS.SetLocalDim(rightBlock->get_op_array(CRE_DES).get_size());
      Dim_tv = LS.GetGlobalDim();
      //loop over a+a on right Block
      //for (int tv=0;tv<rightBlock->get_op_array(CRE_DES).get_size();tv++){
      for (int tv=0;tv<Dim_tv;tv++){
        //determine whether this is a dummy iteration
        bool DummyIter = LS.DummyIter(tv);
        //get the left operator
        boost::shared_ptr<SparseMatrix> RightOpSing;
        boost::shared_ptr<SparseMatrix> RightOpTrip;
        if (!DummyIter){
          RightOpSing = rightBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(rightBlock);
          RightOpTrip = rightBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(rightBlock);
        }
        else{
          //if this is a dummy iteration generate dummy operators
          RightOpSing = boost::make_shared<Cre>();
          RightOpTrip = boost::make_shared<Cre>();
          RightOpSing->set_orbs().push_back(-1);
          RightOpSing->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpSing->resize_deltaQuantum(1);
          RightOpTrip->resize_deltaQuantum(1);
          RightOpSing->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(0),IrrepSpace(0));
          RightOpTrip->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(2),IrrepSpace(0));
       }
        //get the orbital indices
        t = RightOpSing->get_orbs(0);
        v = RightOpSing->get_orbs(1);
        //the left operator quanta
        SpinQuantum ropQS = RightOpSing->get_deltaQuantum(0);
        SpinQuantum ropQT = RightOpTrip->get_deltaQuantum(0);
        //get the integrals
        boost::shared_ptr<Matrix> Ktv = IKJL.GetMatrix(t+NInternal,v+NInternal);
        //loop over a+ on the right block
        for (int ucount=0;ucount<leftBlock->get_op_array(CRE).get_size();ucount++){
          //get the right operator
          boost::shared_ptr<SparseMatrix> LeftOp = leftBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(leftBlock);
          //get the orbital index
          u = LeftOp->get_orbs(0);
          //check orbital indices for validity
          //if (!CheckAllowedOperator(t,u,v,dotIndex,iterCase,true)) continue;
          //the right operator quanta
          SpinQuantum lopQ = LeftOp->get_deltaQuantum(0);
          //the total Operator Quanta (is the same for the left singlet and triplet operator)
          vector<SpinQuantum> VopQS = ropQS + lopQ;
          //the possible quanta for the resulting wavefunction
          vector <SpinQuantum> VQ = WFQ + VopQS[0];
          //evaluate the prefactors
          vector<double> TripFac;
          vector<double> SingFac;
          vector<double> SingFac_;
          for (i=0;i<VQ.size();i++){
            twoS_ = (double)VQ[i].get_s().getirrep();
            S_    = twoS_/2.0;
            fac   = sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
            fac_ = sqrt(2.0)*sqrt(2.0*S_+1.0)/sqrt(2*S+1.0);//this second term arises from case 3
            SingFac.push_back(-fac);
            SingFac_.push_back(fac_);
            ///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // NOTE:
            // for some strange reason I dont't understand we have to add a minus
            // sign if the triplet E operator is on the right Block. 
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            TripFac.push_back(-fac*sqrt(3.0));
          }//resulting wavefunction quanta
          //get the integrals for case 3
          boost::shared_ptr<Matrix> Kut = IKJL.GetMatrix(u+NInternal,t+NInternal);
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //----------------------------
            //construct the operators Vtvu
            //----------------------------
            //generate the wavefunctions that holds the result of the multiplications
            boost::shared_ptr<Wavefunction> VtvuS (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtvuT (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> Vtvu (new Wavefunction(VQ[iquanta],&big,true));
            //multiply the Wavefunction with the operators
            if (!DummyIter){
              operatorfunctions::TensorMultiply(leftBlock,*LeftOp,*RightOpSing,&big,WF,*VtvuS,VopQS[0],1.0);
              operatorfunctions::TensorMultiply(leftBlock,*LeftOp,*RightOpTrip,&big,WF,*VtvuT,VopQS[0],1.0);
              //sum up the contributions
              ScaleAdd(SingFac[iquanta],*VtvuS,*Vtvu);
              ScaleAdd(TripFac[iquanta],*VtvuT,*Vtvu);
            }
            //--------------------------
            //do the same thing for Vvtu
            //--------------------------
            boost::shared_ptr<Wavefunction> Vvtu (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VvtuS (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VvtuT (new Wavefunction(VQ[iquanta],&big,true));
            if (t!=v){
              if (!DummyIter){
                operatorfunctions::TensorMultiply(leftBlock,*LeftOp,Transposeview(*RightOpSing),&big,WF,*VvtuS,VopQS[0],1.0);
                operatorfunctions::TensorMultiply(leftBlock,*LeftOp,Transposeview(*RightOpTrip),&big,WF,*VvtuT,VopQS[0],1.0);
                ScaleAdd(SingFac[iquanta],*VvtuS,*Vvtu);
                ScaleAdd(-TripFac[iquanta],*VvtuT,*Vvtu);//account for transposition and add a minus sign
              }
            }//t!=u
            if (!T[iquanta].GetBasisOnly()){
              //store the operators
              T[iquanta].ResetBuffer();
              bool EndOfArray = false;
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
                if (!EndOfArray){
                  /*
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                  if (fabs(Ktv->element(i,u+NInternal))>1.e-18){
                    sprintf(msg,"\ncase 2: (%i%ii|%i%i) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock right a",t,i,v,u,t,u,v,iquanta,SingFac[iquanta],TripFac[iquanta],Ktv->element(i,u+NInternal),DotProduct(*Vtvu,*Vtvu));pout << msg;
                  }
                  if ((t!=v)&&(fabs(Ktv->element(u+NInternal,i))>1.e-18)){
                    sprintf(msg,"\ncase 2: (%i%i|%i%ii) V(%i,%i,%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%lf twoBlock right b",t,u,v,i,v,u,t,iquanta,SingFac[iquanta],TripFac[iquanta],Ktv->element(u+NInternal,i),DotProduct(*Vvtu,*Vvtu));pout << msg;
                  }

                  if (fabs(Kut->element(i,v+NInternal))>1.e-18){
                    sprintf(msg,"\ncase 3: (%i%ii|%i%i) iquanta=%i fac(singlet)=%4.6e int=%lf Norm=%lf twoBlock right a",u,i,t,v,u,t,v,iquanta,SingFac_[iquanta],Kut->element(i,v+NInternal),4.0*DotProduct(*VtvuS,*VtvuS));pout << msg;
                  }
                  if ((t!=v)&&(fabs(Kut->element(i,v+NInternal))>1.e-18)){
                    sprintf(msg,"\ncase 3: (%i%ii|%i%i) iquanta=%i fac(singlet)=%4.6e int=%lf Norm=%lf twoBlock right b",u,i,t,v,u,v,t,iquanta,SingFac_[iquanta],Kut->element(i,v+NInternal),4.0*DotProduct(*VvtuS,*VvtuS));pout << msg;
                  }
                   */
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                  ScaleAdd(Ktv->element(i,u+NInternal),*Vtvu,*Vi);
                  ScaleAdd(Kut->element(i,v+NInternal)*SingFac_[iquanta],*VtvuS,*Vi);//case 3
                  ScaleAdd(Ktv->element(u+NInternal,i),*Vvtu,*Vi);
                  ScaleAdd(Kut->element(i,v+NInternal)*SingFac_[iquanta],*VvtuS,*Vi);//case 3
                  T[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
                }
              }//i
              //Send Vtuv around to other processes and receive theirs
              SendAroundVtuv_i(*Vtvu,*VtvuS,*Vvtu,*VvtuS,*Ktv,*Kut,T[iquanta],u+NInternal,v+NInternal,DummyIter,SingFac_[iquanta]);
            }//!BasisOnly
            T[iquanta].AddToBasisOperator(Vtvu);
            T[iquanta].AddToBasisOperator(VtvuS,SingFac_[iquanta]);
            T[iquanta].AddToBasisOperator(Vvtu);
            T[iquanta].AddToBasisOperator(VvtuS,SingFac_[iquanta]);
          }//iquanta
        }//u
      }//tv
    }//final iteration
    
    //close the integral and operator container 
    IKJL.CloseFileRead();
    for (int i=0;i<T.size();i++){
      T[i].CloseFileRead();
    }
  }
  
  //============================================================================
  //Construct the one-electron part of the V(i) operator
  //============================================================================
  void ConstructCre(SpinBlock &big, Matrix &heff, Wavefunction &WF, vector<WavefunctionArray> &T){
    
    char msg[512];
    char BaseName[512];
    int OrbWin[6];
    double ENuc=0.0;
    bool ConvOverlap;
    //get the basename and stuff
    ReadInput(BaseName,OrbWin,ENuc,ConvOverlap);
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
    int t,u,v,i,a,tu;
    int dummy;

    Matrix D1_,D1;
    D1_.ReSize(NActive,NActive);
    D1.ReSize(NActive,NActive);
    
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    double twoS = (double) WFQ.get_s().getirrep();
    double S    = twoS/2.0;
    double twoS_,S_,fac;
    
    //open the operator container 
    for (int i=0;i<T.size();i++){
      T[i].OpenFileRead();
    }
    
    //-------------
    //leftBlock
    //-------------
    for (int tcount=0;tcount<leftBlock->get_op_array(CRE).get_size();tcount++){
      //get the operator
      boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
      //get the orbital label
      t = op->get_orbs(0);
      //get the operator quanta
      SpinQuantum opQ = (*op).get_deltaQuantum(0);
      SpinQuantum opQT = Transposeview(*op).get_deltaQuantum(0);
      //the possible quanta for the resulting wavefunction
      vector<SpinQuantum> VQ = opQ + WFQ;
      vector<SpinQuantum> VQT = opQT + WFQ;
      //the prefactors
      vector<double> Fac;
      //evaluate the prefactors
      for (int iquanta=0;iquanta<VQ.size();iquanta++){
        twoS_ = (double) VQ[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt(2.0*S_+1.0)/sqrt(2.0*S+1.0);
        Fac.push_back(fac);
      }
      for (int iquanta=0;iquanta<VQ.size();iquanta++){
        //the wavefunction that holds the result of the multiplication 
        boost::shared_ptr<Wavefunction> Vt(new Wavefunction(VQ[iquanta],&big,true));
        boost::shared_ptr<Wavefunction> VtT(new Wavefunction(VQT[iquanta],&big,true));
        //do the multiplication with the wavefunction
        operatorfunctions::TensorMultiply(leftBlock,*op,&big,WF,*Vt,opQ,Fac[iquanta]);
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*op),&big,WF,*VtT,opQT,Fac[iquanta]);
        if (!T[iquanta].GetBasisOnly()){
          //store the operator
          T[iquanta].ResetBuffer();
          bool EndOfArray = false;
          while (!EndOfArray){
            boost::shared_ptr<Wavefunction> Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
            if (!EndOfArray){
              ScaleAdd(heff.element(i,t+NInternal),*Vt,*Vi);
              T[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
            }
          }//i
        }//!BasisOnly
        T[iquanta].AddToBasisOperator(Vt);
      }//iquanta
    }//t on the left side
    //-------------
    //rightBlock
    //-------------
    for (int ucount=0;ucount<rightBlock->get_op_array(CRE).get_size();ucount++){
      //get the oeprator
      boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
      //get the orbital label
      u = rop->get_orbs(0);
      //the operator quanta
      SpinQuantum ropQ = (*rop).get_deltaQuantum(0);
      SpinQuantum ropQT = Transposeview(*rop).get_deltaQuantum(0);
      //the possible quanta for the resulting wavefunction
      vector<SpinQuantum> VQ = ropQ + WFQ;
      vector<SpinQuantum> VQT = ropQT + WFQ;
      //the prefactors
      vector<double> Fac;
      for (int iquanta=0;iquanta<VQ.size();iquanta++){
        twoS_ = (double) VQ[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt(2.0*S_+1.0)/sqrt(2.0*S+1.0);
        Fac.push_back(fac);
      }//iquanta
      for (int iquanta=0;iquanta<VQ.size();iquanta++){
        //the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> Vu (new Wavefunction(VQ[iquanta],&big,true));
        boost::shared_ptr<Wavefunction> VuT (new Wavefunction(VQT[iquanta],&big,true));
        //multiply the operator with the wavefunction
        operatorfunctions::TensorMultiply(rightBlock,*rop,&big,WF,*Vu,ropQ,Fac[iquanta]);
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*rop),&big,WF,*VuT,ropQ,Fac[iquanta]);
        if (!T[iquanta].GetBasisOnly()){
          //store the operator
          T[iquanta].ResetBuffer();
          bool EndOfArray = false;
          while (!EndOfArray){
            boost::shared_ptr<Wavefunction> Vi = T[iquanta].GetOpFromBuffer(dummy,i,EndOfArray);
            if (!EndOfArray){
              ScaleAdd(heff.element(i,u+NInternal),*Vu,*Vi);
              T[iquanta].ReplaceOpInBuffer(Vi,dummy,i);
            }
          }//i
        }//!BasisOnly
        T[iquanta].AddToBasisOperator(Vu);
      }//iquanta
    }//u on the right side
    //close operator container 
    for (int i=0;i<T.size();i++){
      T[i].CloseFileRead();
    }
    
    
  }
  //============================================================================
  //the driver for the construction of the a+a+a|psi> operators
  //============================================================================
  void ConstructCreCreDes(SpinBlock &big, vector<Wavefunction> &WF, ThreeIndOpArray &CCD,
                          IntegralContainer &IKJL, SweepParams &sweepParams,NEVPT2Info &Info){
    char msg[512];
    int OrbWin[6];
    //get the necessary infos
    char BaseName [512];
    Info.getBaseName(BaseName);
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
    int t,u,v,i,a,tu;
    int iSweep = Info.GetNevSweep();
    int MaxIter = Info.GetMaxBlockIter(iSweep);
    int iter = sweepParams.get_block_iter();
    
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *dotBlock = leftBlock->get_rightBlock();
    int MaxCore = Info.GetMaxCore();
    int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
    double StartTime=GetTime();

    //--------------------------------
    //build the two-electron operators
    //--------------------------------
    if (iter==0){
      //first iteration
      ConstructCreCreDesSingleBlock(big,big.get_leftBlock(),CCD,IKJL,_INITIAL_ITER_,WF[0]);
    }
    else if ((iSweep==1)&&(iter==MaxIter)){
      //last iteration
      ConstructCreCreDesSingleBlock(big,big.get_leftBlock()->get_rightBlock(),CCD,IKJL,_FINAL_ITER_,WF[0]);
      ConstructCreCreDes_1_2_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CCD,OrbWin);
      ConstructCreCreDes_2_1_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CCD,OrbWin,WF[0]);
    }
    else{
      //a regular iteration
      ConstructCreCreDesSingleBlock(big,big.get_leftBlock()->get_rightBlock(),CCD,IKJL,_REGULAR_ITER_,WF[0]);
      ConstructCreCreDes_1_2_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CCD,OrbWin);
      ConstructCreCreDes_2_1_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CCD,OrbWin,WF[0]);
    }
    
    //--------------------------------------------------------------------------
    //when in the middle of the sweep, actually build the V(i)|psi> functions
    //and evaluate the energy contribution
    //--------------------------------------------------------------------------
    if ((iSweep==1)&&(iter==MaxIter)){
      //get the one-electron effective Hamiltonian
      Matrix heff;
      Info.GetH(1,heff);
      //receive the OperatorArray from disk. Note: CDD_ here was CDD in function sweep_nevpt2.C:"do_one"
      ThreeIndOpArray CCD_;
      sprintf(msg,"%s.CCD.tmp",BaseName);
      CCD_.Retrieve(msg);
      for (int iroot=0;iroot<WF.size();iroot++){
        //------------------------------------
        //first fill the V(i) arrays with zero
        //------------------------------------
        SpinQuantum WFQ = WF[iroot].get_deltaQuantum(0);
        vector<SpinQuantum> VQ = WFQ + SpinQuantum(1,SpinSpace(1),IrrepSpace(0));
        //generate and initialize the array
        vector<WavefunctionArray> Vi;
        Vi.resize(VQ.size());
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          sprintf(msg,"%s.Ti.%i.%i.tmp",BaseName,iroot,iquanta);
          Vi[iquanta].Initialize(NInternal,msg,mpi_world_size());
          Vi[iquanta].CalcBufferSize(MaxCore,MaxM);
          //open the array
          Vi[iquanta].OpenFileWrite();
          Vi[iquanta].SetQuanta(VQ[iquanta]);
          //divide the internal orbital space between the processes
          int i_loc_start,i_loc_end;
          int loc_NInternal;
          PALDivideLoop(i_loc_start,i_loc_end,0,NInternal);
          loc_NInternal = i_loc_end-i_loc_start;
          for (i=i_loc_start;i<i_loc_end;i++){
            boost::shared_ptr<Wavefunction> Starter (new Wavefunction(VQ[iquanta],&big,true));
            Vi[iquanta].AppendOperator(Starter,0,i);
          }//a
          //close the array
          Vi[iquanta].CloseFileWrite();
          //create the empty basis state
          boost::shared_ptr<Wavefunction> EmptyBasis (new Wavefunction(VQ[iquanta],&big,true));
          Vi[iquanta].SetBasisOperator(EmptyBasis);
        }//iquanta
        //------------------------
        //then build the operators
        //------------------------
        ConstructViSingleBlock(big,WF[iroot],CCD,CCD_,Vi,IKJL,Info);
        ConstructCreCreDesTwoBlocks(big,WF[iroot],Vi,IKJL,_FINAL_ITER_,OrbWin);
        ConstructCre(big,heff,WF[iroot],Vi);
        //--------------------------------
        //evaluate the energy contribution
        //--------------------------------
        double Ei = V_i(Vi,big,WF[iroot],Info,iroot);
        Info.SetE(iroot,7,Ei);
        //--------
        //clean up
        //--------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          Vi[iquanta].Clear();
        }//iquanta
      }//iroot
      CCD.Clear();
      CCD_.Clear();
    }//final iteration
    //take the time and add it to the NEVPT2 time
    double FinishTime=GetTime();
    double time = FinishTime-StartTime;
    Info.AddTime(4,time);
    Info.AddTimePerClass(7,time);
  }
  
  
  //============================================================================
  // Construct a+aa|psi> where two indices are on the dot and the other on the
  // left block
  //============================================================================
  void ConstructCreDesDes_1_2_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CDD, IntegralContainer &IKJA, int *OrbWin){
    char msg[512];
    int t,u,v;
    SpinBlock *leftBlock =  big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //open the OperatorArray
    CDD.OpenFileWrite();
    //------------------------------------------------------
    // Case 3: a+ is on leftBlock and aa is on the dotBlock
    //------------------------------------------------------
    for (int uv=0;uv<dotBlock->get_op_array(CRE_CRE).get_size();uv++){
      //make sure we do not calculate the operator on multiple processors
      if (SkipOperator(big.get_sites(),dotBlock->get_sites()[0])) continue;
      //get the operator
      boost::shared_ptr<SparseMatrix> dotOpSing = dotBlock->get_op_array(CRE_CRE).get_local_element(uv)[0]->getworkingrepresentation(dotBlock);
      boost::shared_ptr<SparseMatrix> dotOpTrip = dotBlock->get_op_array(CRE_CRE).get_local_element(uv)[1]->getworkingrepresentation(dotBlock);
      //get the orbital labels, note: we assume that the dot has only one label
      u = dotOpSing->get_orbs(0);
      v = dotOpSing->get_orbs(1);
      //get the dot operator quanta
      SpinQuantum dOpQ = Transposeview(*dotOpSing).get_deltaQuantum(0);
      for (int tcount=0;tcount<lBlock->get_op_array(CRE).get_size();tcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> leftOp = lBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(lBlock);
        //get the orbital label
        t = leftOp->get_orbs(0);
        
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //sprintf(msg,"\nt=%i u=%i v=%i   1_2_0 case 3",t,u,v);pout << msg;
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        //get the left operator quanta
        SpinQuantum lopQ = leftOp->get_deltaQuantum(0);
        //the possible quanta for the resulting operator
        SpinQuantum opQ = (dOpQ + lopQ)[0];
        //Generate and construct the operator a+a+a (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        OpS.set_orbs().push_back(t);
        OpS.set_orbs().push_back(u);
        OpS.set_orbs().push_back(v);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = opQ;
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs().push_back(t);
        OpT.set_orbs().push_back(u);
        OpT.set_orbs().push_back(v);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = opQ;
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,*leftOp,Transposeview(*dotOpSing),leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,*leftOp,Transposeview(*dotOpTrip),leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operator
        boost::shared_ptr<Cre> Op (new Cre);
        Op->set_orbs().push_back(t);
        Op->set_orbs().push_back(u);
        Op->set_orbs().push_back(v);
        Op->set_initialised() = true;
        Op->set_fermion() = true;
        Op->resize_deltaQuantum(1);
        Op->set_deltaQuantum(0) = opQ;
        Op->allocate(leftBlock->get_stateInfo());
        //sum the two operators up
        ScaleAdd(1.0/sqrt(2.0),OpS,*Op);
        ScaleAdd(sqrt(3.0/2.0),OpT,*Op);//account for transposition and add a minus sign
        //Store the operator
        CDD.AppendOperator(Op,t,u,v,false);
        //clean up
        OpS.Clear();
        OpT.Clear();
      }//v
    }//tu
    
    //--------------------------------------------------------------
    // Cases 2 and 1: a is on leftBlock and a(+)a is on the dotBlock
    //--------------------------------------------------------------
    for (int tv=0;tv<dotBlock->get_op_array(CRE_DES).get_size();tv++){
      //make sure we do not calculate the operator on multiple processors
      if (SkipOperator(big.get_sites(),dotBlock->get_sites()[0])) continue;
      //get the left operator
      boost::shared_ptr<SparseMatrix> dotOpSing = dotBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(dotBlock);
      boost::shared_ptr<SparseMatrix> dotOpTrip = dotBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(dotBlock);
      //get the orbital indices
      t = dotOpSing->get_orbs(0);
      v = dotOpSing->get_orbs(1);
      //the left operator quanta
      SpinQuantum dopQS = dotOpSing->get_deltaQuantum(0);
      SpinQuantum dopQT = dotOpTrip->get_deltaQuantum(0);
      //loop over a+ on the right block
      for (int ucount=0;ucount<lBlock->get_op_array(CRE).get_size();ucount++){
        //get the right operator
        boost::shared_ptr<SparseMatrix> leftOp = lBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(lBlock);
        //get the orbital index
        u = leftOp->get_orbs(0);
        
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //sprintf(msg,"\nt=%i u=%i v=%i   1_2_0 case 2",t,u,v);pout << msg;
        //sprintf(msg,"\nt=%i u=%i v=%i   1_2_0 case 1",t,v,u);pout << msg;
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        //the right operator quanta
        SpinQuantum lopQ = Transposeview(*leftOp).get_deltaQuantum(0);
        //the total Operator Quanta (is the same for the left singlet and triplet operator)
        vector<SpinQuantum> VopQS = dopQS + lopQ;
        //Generate and construct the operator a+aa (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        Cre OpS_,OpT_;
        OpS.set_orbs().push_back(t);
        OpS.set_orbs().push_back(u);
        OpS.set_orbs().push_back(v);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = VopQS[0];
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs().push_back(t);
        OpT.set_orbs().push_back(u);
        OpT.set_orbs().push_back(v);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = VopQS[0];
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,Transposeview(*leftOp),*dotOpSing,leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,Transposeview(*leftOp),*dotOpTrip,leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operator (case 2)
        boost::shared_ptr<Cre> Op_tuv (new Cre);
        Op_tuv->set_orbs().push_back(t);
        Op_tuv->set_orbs().push_back(u);
        Op_tuv->set_orbs().push_back(v);
        Op_tuv->set_initialised() = true;
        Op_tuv->set_fermion() = true;
        Op_tuv->resize_deltaQuantum(1);
        Op_tuv->set_deltaQuantum(0) = VopQS[0];
        Op_tuv->allocate(leftBlock->get_stateInfo());
        //Sum the operators up
        ScaleAdd(-1.0/sqrt(2.0),OpS,*Op_tuv);
        ///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // NOTE:
        // for some strange reason I dont't understand we have to add a minus
        // sign if the triplet E operator is on the right Block. 
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ScaleAdd(sqrt(3.0/2.0),OpT,*Op_tuv);
        //Store the operator
        CDD.AppendOperator(Op_tuv,t,u,v,false);
        //generate the product operator (case 1)
        boost::shared_ptr<Cre> Op_tvu (new Cre);
        Op_tvu->set_orbs().push_back(t);
        Op_tvu->set_orbs().push_back(v);
        Op_tvu->set_orbs().push_back(u);
        Op_tvu->set_initialised() = true;
        Op_tvu->set_fermion() = true;
        Op_tvu->resize_deltaQuantum(1);
        Op_tvu->set_deltaQuantum(0) = VopQS[0];
        Op_tvu->allocate(leftBlock->get_stateInfo());
        //Sum the operators up
        ScaleAdd(sqrt(2.0),OpS,*Op_tvu);
        //Store the operator
        CDD.AppendOperator(Op_tvu,t,v,u);
        //clean up
        OpS.Clear();
        OpT.Clear();
        if (t!=v){
          //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //sprintf(msg,"\nt=%i u=%i v=%i   1_2_0 case 2b",v,u,t);pout << msg;
          //sprintf(msg,"\nt=%i u=%i v=%i   1_2_0 case 1b",v,t,u);pout << msg;
          //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //do the same thing for the transposed dot operator
          OpS_.set_orbs().push_back(v);
          OpS_.set_orbs().push_back(u);
          OpS_.set_orbs().push_back(t);
          OpS_.set_initialised() = true;
          OpS_.set_fermion() = true;
          OpS_.resize_deltaQuantum(1);
          OpS_.set_deltaQuantum(0) = VopQS[0];
          OpS_.allocate(leftBlock->get_stateInfo());
          OpT_.set_orbs().push_back(v);
          OpT_.set_orbs().push_back(u);
          OpT_.set_orbs().push_back(t);
          OpT_.set_initialised() = true;
          OpT_.set_fermion() = true;
          OpT_.resize_deltaQuantum(1);
          OpT_.set_deltaQuantum(0) = VopQS[0];
          OpT_.allocate(leftBlock->get_stateInfo());
          operatorfunctions::TensorProduct(dotBlock,Transposeview(*dotOpSing),Transposeview(*leftOp),leftBlock,&(leftBlock->get_stateInfo()),OpS_,1.0);
          operatorfunctions::TensorProduct(dotBlock,Transposeview(*dotOpTrip),Transposeview(*leftOp),leftBlock,&(leftBlock->get_stateInfo()),OpT_,1.0);
          //generate the product operator (case 2)
          boost::shared_ptr<Cre> Op_vut (new Cre);
          Op_vut->set_orbs().push_back(v);
          Op_vut->set_orbs().push_back(u);
          Op_vut->set_orbs().push_back(t);
          Op_vut->set_initialised() = true;
          Op_vut->set_fermion() = true;
          Op_vut->resize_deltaQuantum(1);
          Op_vut->set_deltaQuantum(0) = VopQS[0];
          Op_vut->allocate(leftBlock->get_stateInfo());
          //Sum the operators up
          ScaleAdd(-1.0/sqrt(2.0),OpS_,*Op_vut);
          ///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          // NOTE:
          // for some strange reason I dont't understand we have to add a minus
          // sign if the triplet E operator is on the right Block. 
          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ScaleAdd(-sqrt(3.0/2.0),OpT_,*Op_vut);//account for transposition and add a minus sign
          //Store the operator
          CDD.AppendOperator(Op_vut,v,u,t,false);
          //generate the product operator (case 1)
          boost::shared_ptr<Cre> Op_vtu (new Cre);
          Op_vtu->set_orbs().push_back(v);
          Op_vtu->set_orbs().push_back(t);
          Op_vtu->set_orbs().push_back(u);
          Op_vtu->set_initialised() = true;
          Op_vtu->set_fermion() = true;
          Op_vtu->resize_deltaQuantum(1);
          Op_vtu->set_deltaQuantum(0) = VopQS[0];
          Op_vtu->allocate(leftBlock->get_stateInfo());
          //Sum the operators up
          ScaleAdd(sqrt(2.0),OpS_,*Op_vtu);
          //Store the operator
          CDD.AppendOperator(Op_vtu,v,t,u);
          //clean up
          OpS.Clear();
          OpT.Clear();

        }
      }//u
    }//tv
    //close the operator array
    CDD.CloseFileWrite();
    
  }
  
  //============================================================================
  // Construct a+aa where two indices are on the leftBlock and the other on the
  // dot Block
  //============================================================================
  void ConstructCreDesDes_2_1_0(SpinBlock &big, SpinBlock *lBlock, SpinBlock *dotBlock,
                                ThreeIndOpArray &CDD, IntegralContainer &IKJA, int *OrbWin){
    char msg[512];
    int t,u,v;
    SpinBlock *leftBlock =  big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //open the OperatorArray
    CDD.OpenFileWrite();
    //------------------------------------------------------
    // Case 3: a+ is on dotBlock and aa is on the leftBlock
    //------------------------------------------------------
    for (int uv=0;uv<lBlock->get_op_array(CRE_CRE).get_size();uv++){
      //get the operator
      boost::shared_ptr<SparseMatrix> lOpSing = lBlock->get_op_array(CRE_CRE).get_local_element(uv)[0]->getworkingrepresentation(lBlock);
      boost::shared_ptr<SparseMatrix> lOpTrip = lBlock->get_op_array(CRE_CRE).get_local_element(uv)[1]->getworkingrepresentation(lBlock);
      //get the orbital labels, note: we assume that the dot has only one label
      u = lOpSing->get_orbs(0);
      v = lOpSing->get_orbs(1);
      //get the dot operator quanta
      SpinQuantum lOpQ = Transposeview(*lOpSing).get_deltaQuantum(0);
      for (int tcount=0;tcount<dotBlock->get_op_array(CRE).get_size();tcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> dotOp = dotBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(dotBlock);
        //get the orbital label
        t = dotOp->get_orbs(0);

        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //sprintf(msg,"\nt=%i u=%i v=%i   2_1_0 case 3",t,u,v);pout << msg;
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        //get the left operator quanta
        SpinQuantum dotopQ = dotOp->get_deltaQuantum(0);
        //the possible quanta for the resulting operator
        SpinQuantum opQ = (lOpQ + dotopQ)[0];
        //Generate and construct the operator a+a+a (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        OpS.set_orbs().push_back(t);
        OpS.set_orbs().push_back(u);
        OpS.set_orbs().push_back(v);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = opQ;
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs().push_back(t);
        OpT.set_orbs().push_back(u);
        OpT.set_orbs().push_back(v);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = opQ;
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,Transposeview(*lOpSing),*dotOp,leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,Transposeview(*lOpTrip),*dotOp,leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operator
        boost::shared_ptr<Cre> Op (new Cre);
        Op->set_orbs().push_back(t);
        Op->set_orbs().push_back(u);
        Op->set_orbs().push_back(v);
        Op->set_initialised() = true;
        Op->set_fermion() = true;
        Op->resize_deltaQuantum(1);
        Op->set_deltaQuantum(0) = opQ;
        Op->allocate(leftBlock->get_stateInfo());
        //sum the two operators up
        ScaleAdd(1.0/sqrt(2.0),OpS,*Op);
        ScaleAdd(sqrt(3.0/2.0),OpT,*Op);//account for the transposition and add a minus sign
        //Store the operator
        CDD.AppendOperator(Op,t,u,v,false);
        //do the same thing for the transposed operators
        if (u!=v){
          //generate the product operator
          boost::shared_ptr<Cre> Op_ (new Cre);
          Op_->set_orbs().push_back(t);
          Op_->set_orbs().push_back(v);
          Op_->set_orbs().push_back(u);
          Op_->set_initialised() = true;
          Op_->set_fermion() = true;
          Op_->resize_deltaQuantum(1);
          Op_->set_deltaQuantum(0) = opQ;
          Op_->allocate(leftBlock->get_stateInfo());
          //sum the two operators up
          ScaleAdd(1.0/sqrt(2.0),OpS,*Op_);
          ScaleAdd(-sqrt(3.0/2.0),OpT,*Op_);//account for the transposition and add a minus sign
          //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //sprintf(msg,"\nt=%i u=%i v=%i   2_1_0 case 3",t,v,u);pout << msg;
          //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //Store the operator
          CDD.AppendOperator(Op_,t,v,u,false);
        }//u!=v
        //clean up
        OpS.Clear();
        OpT.Clear();

      }//v
    }//tu
    
    //--------------------------------------------------------------
    // Cases 2 and 1: a is on dotBlock and a(+)a is on the leftBlock
    //--------------------------------------------------------------
    for (int tv=0;tv<lBlock->get_op_array(CRE_DES).get_size();tv++){
      //get the left operator
      boost::shared_ptr<SparseMatrix> lOpSing = lBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(lBlock);
      boost::shared_ptr<SparseMatrix> lOpTrip = lBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(lBlock);
      //get the orbital indices
      t = lOpSing->get_orbs(0);
      v = lOpSing->get_orbs(1);
      //the left operator quanta
      SpinQuantum lopQS = lOpSing->get_deltaQuantum(0);
      SpinQuantum lopQT = lOpTrip->get_deltaQuantum(0);
      //loop over a+ on the right block
      for (int ucount=0;ucount<dotBlock->get_op_array(CRE).get_size();ucount++){
        //get the right operator
        boost::shared_ptr<SparseMatrix> dotOp = dotBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(dotBlock);
        //get the orbital index
        u = dotOp->get_orbs(0);
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //sprintf(msg,"\nt=%i u=%i v=%i   2_1_0 case 2",t,u,v);pout << msg;
        //sprintf(msg,"\nt=%i u=%i v=%i   2_1_0 case 1",t,v,u);pout << msg;
        //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //the right operator quanta
        SpinQuantum dopQ = Transposeview(*dotOp).get_deltaQuantum(0);
        //the total Operator Quanta (is the same for the left singlet and triplet operator)
        vector<SpinQuantum> VopQS = lopQS + dopQ;
        //Generate and construct the operator a+aa (singlet parentage and triplet parentage)
        Cre OpS,OpT;
        Cre OpS_,OpT_;
        OpS.set_orbs().push_back(t);
        OpS.set_orbs().push_back(u);
        OpS.set_orbs().push_back(v);
        OpS.set_initialised() = true;
        OpS.set_fermion() = true;
        OpS.resize_deltaQuantum(1);
        OpS.set_deltaQuantum(0) = VopQS[0];
        OpS.allocate(leftBlock->get_stateInfo());
        OpT.set_orbs().push_back(t);
        OpT.set_orbs().push_back(u);
        OpT.set_orbs().push_back(v);
        OpT.set_initialised() = true;
        OpT.set_fermion() = true;
        OpT.resize_deltaQuantum(1);
        OpT.set_deltaQuantum(0) = VopQS[0];
        OpT.allocate(leftBlock->get_stateInfo());
        operatorfunctions::TensorProduct(lBlock,*lOpSing,Transposeview(*dotOp),leftBlock,&(leftBlock->get_stateInfo()),OpS,1.0);
        operatorfunctions::TensorProduct(lBlock,*lOpTrip,Transposeview(*dotOp),leftBlock,&(leftBlock->get_stateInfo()),OpT,1.0);
        //generate the product operator (case 2)
        boost::shared_ptr<Cre> Op_tuv (new Cre);
        Op_tuv->set_orbs().push_back(t);
        Op_tuv->set_orbs().push_back(u);
        Op_tuv->set_orbs().push_back(v);
        Op_tuv->set_initialised() = true;
        Op_tuv->set_fermion() = true;
        Op_tuv->resize_deltaQuantum(1);
        Op_tuv->set_deltaQuantum(0) = VopQS[0];
        Op_tuv->allocate(leftBlock->get_stateInfo());
        //Sum the operators up
        ScaleAdd(-1.0/sqrt(2.0),OpS,*Op_tuv);
        ScaleAdd(-sqrt(3.0/2.0),OpT,*Op_tuv);
        //Store the operator
        CDD.AppendOperator(Op_tuv,t,u,v,false);
        //generate the product operator (case 1)
        boost::shared_ptr<Cre> Op_tvu (new Cre);
        Op_tvu->set_orbs().push_back(t);
        Op_tvu->set_orbs().push_back(v);
        Op_tvu->set_orbs().push_back(u);
        Op_tvu->set_initialised() = true;
        Op_tvu->set_fermion() = true;
        Op_tvu->resize_deltaQuantum(1);
        Op_tvu->set_deltaQuantum(0) = VopQS[0];
        Op_tvu->allocate(leftBlock->get_stateInfo());
        //Sum the operators up
        ScaleAdd(sqrt(2.0),OpS,*Op_tvu);
        //Store the operator
        CDD.AppendOperator(Op_tvu,t,v,u);
        //clean up
        OpS.Clear();
        OpT.Clear();
        //do the same thing for the transposed dot operator
        if (t!=v){
          //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //sprintf(msg,"\nt=%i u=%i v=%i   2_1_0 case 2b",v,u,t);pout << msg;
          //sprintf(msg,"\nt=%i u=%i v=%i   2_1_0 case 1b",v,t,u);pout << msg;
          //!!!!!!!!!!!!!!!!DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          OpS_.set_orbs().push_back(v);
          OpS_.set_orbs().push_back(u);
          OpS_.set_orbs().push_back(t);
          OpS_.set_initialised() = true;
          OpS_.set_fermion() = true;
          OpS_.resize_deltaQuantum(1);
          OpS_.set_deltaQuantum(0) = VopQS[0];
          OpS_.allocate(leftBlock->get_stateInfo());
          OpT_.set_orbs().push_back(v);
          OpT_.set_orbs().push_back(u);
          OpT_.set_orbs().push_back(t);
          OpT_.set_initialised() = true;
          OpT_.set_fermion() = true;
          OpT_.resize_deltaQuantum(1);
          OpT_.set_deltaQuantum(0) = VopQS[0];
          OpT_.allocate(leftBlock->get_stateInfo());
          operatorfunctions::TensorProduct(lBlock,Transposeview(*lOpSing),Transposeview(*dotOp),leftBlock,&(leftBlock->get_stateInfo()),OpS_,1.0);
          operatorfunctions::TensorProduct(lBlock,Transposeview(*lOpTrip),Transposeview(*dotOp),leftBlock,&(leftBlock->get_stateInfo()),OpT_,1.0);
          //generate the product operator (case 2)
          boost::shared_ptr<Cre> Op_vut (new Cre);
          Op_vut->set_orbs().push_back(v);
          Op_vut->set_orbs().push_back(u);
          Op_vut->set_orbs().push_back(t);
          Op_vut->set_initialised() = true;
          Op_vut->set_fermion() = true;
          Op_vut->resize_deltaQuantum(1);
          Op_vut->set_deltaQuantum(0) = VopQS[0];
          Op_vut->allocate(leftBlock->get_stateInfo());
          //Sum the operators up
          ScaleAdd(-1.0/sqrt(2.0),OpS_,*Op_vut);
          ScaleAdd(sqrt(3.0/2.0),OpT_,*Op_vut);//account for the transposition and add a minus sign
          //Store the operator
          CDD.AppendOperator(Op_vut,v,u,t,false);
          //generate the product operator (case 1)
          boost::shared_ptr<Cre> Op_vtu (new Cre);
          Op_vtu->set_orbs().push_back(v);
          Op_vtu->set_orbs().push_back(t);
          Op_vtu->set_orbs().push_back(u);
          Op_vtu->set_initialised() = true;
          Op_vtu->set_fermion() = true;
          Op_vtu->resize_deltaQuantum(1);
          Op_vtu->set_deltaQuantum(0) = VopQS[0];
          Op_vtu->allocate(leftBlock->get_stateInfo());
          //Sum the operators up
          ScaleAdd(sqrt(2.0),OpS_,*Op_vtu);
          //Store the operator
          CDD.AppendOperator(Op_vtu,v,t,u);
          //clean up
          OpS.Clear();
          OpT.Clear();
        }
      }//u
    }//tv
    //close the operator array
    CDD.CloseFileWrite();
    
  }
  
  
  
  //============================================================================
  // Construct a+aa|psi> on two blocks block
  //============================================================================
  void ConstructCreDesDesTwoBlocks(SpinBlock &big, Wavefunction &WF, vector<WavefunctionArray> &T, 
                                 IntegralContainer &IKJA,int iterCase, int *OrbWin){
    char msg[512];
    //get the orbital spaces and stuff
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
    int t,u,v,i,a,tu;
    int dummy;
    
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *dotBlock = leftBlock->get_rightBlock();
    int dotIndex = dotBlock->get_sites()[0];//let's hope that this Block has only one member......
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    double twoS = (double) WFQ.get_s().getirrep();
    double S    = twoS/2.0;
    double twoS_,S_,fac;
    int NumProcs = mpi_world_size();

    //Note: we have to synchronize the loops over (t,u,v) and broadcast the
    //operators, because each processor needs every possible operator O(tuv).
    //The division between processes is made in index a
#ifndef SERIAL
    mpi::communicator world;
#endif

    //open the integral and operator container
    IKJA.OpenFileRead();
    for (int i=0;i<T.size();i++){
      T[i].OpenFileRead();
    }
    //--------------------------------------------------------------------------
    //Case 3: the creation operator is on one block while the two annihilation
    //        operators are on the other block
    //--------------------------------------------------------------------------
    //make sure that allprocesses undergo the same number of cycles in the loop
    LoopSynchronizer LS(leftBlock->get_op_array(CRE_CRE).get_size());
    int Dim_uv=LS.GetGlobalDim();
    //loop over a+a+ on the left block
    //for (int uv=0;uv<leftBlock->get_op_array(CRE_CRE).get_size();uv++){
    for (int uv=0;uv<Dim_uv;uv++){
      //determine whetehr this is a dummy iteration
      bool DummyIter = LS.DummyIter(uv);
      //get the left operator
      boost::shared_ptr<SparseMatrix> LeftOpSing;
      boost::shared_ptr<SparseMatrix> LeftOpTrip;
      if (!DummyIter){
        LeftOpSing = leftBlock->get_op_array(CRE_CRE).get_local_element(uv)[0]->getworkingrepresentation(leftBlock);
        LeftOpTrip = leftBlock->get_op_array(CRE_CRE).get_local_element(uv)[1]->getworkingrepresentation(leftBlock);
      }
      else{
        //if this is a dummy iteration generate dummy oeprators
        LeftOpSing = boost::make_shared<Cre>();
        LeftOpTrip = boost::make_shared<Cre>();
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpSing->resize_deltaQuantum(1);
        LeftOpTrip->resize_deltaQuantum(1);
        LeftOpSing->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(0),IrrepSpace(0));
        LeftOpTrip->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(2),IrrepSpace(0));
      }
      //get the orbital indices
      u = LeftOpSing->get_orbs(0);
      v = LeftOpSing->get_orbs(1);
      //the left operator quanta
      SpinQuantum lopQS = Transposeview(*LeftOpSing).get_deltaQuantum(0);
      SpinQuantum lopQT = Transposeview(*LeftOpTrip).get_deltaQuantum(0);
      //get the integrals
      boost::shared_ptr<Matrix> Kuv = IKJA.GetMatrix(u+NInternal,v+NInternal);
      boost::shared_ptr<Matrix> Kvu = IKJA.GetMatrix(v+NInternal,u+NInternal);
      //loop over a on the right block
      for (int tcount=0;tcount<rightBlock->get_op_array(CRE).get_size();tcount++){
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOp = rightBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(rightBlock);
        //get the orbital index
        t = RightOp->get_orbs(0);
        //the right operator quanta
        SpinQuantum ropQ = RightOp->get_deltaQuantum(0);
        //the total Operator Quanta (is the same for the left singlet and triplet operator)
        vector<SpinQuantum> VopQS = lopQS + ropQ;
        //the possible quanta for the resulting wavefunction
        vector <SpinQuantum> VQ = WFQ + VopQS[0];
        //evaluate the prefactors
        vector<double> TripFac;
        vector<double> SingFac;
        for (i=0;i<VQ.size();i++){
          twoS_ = (double)VQ[i].get_s().getirrep();
          S_    = twoS_/2.0;
          fac   = sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
          SingFac.push_back(fac);
          TripFac.push_back(fac*sqrt(3.0));//on account of the transposition of the triplet operator the minus sign is cancelled
        }//resulting wavefunction quanta

        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //generate the wavefunctions that holds the result of the multiplications
          boost::shared_ptr<Wavefunction> VtuvS (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VtuvT (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> Vtuv (new Wavefunction(VQ[iquanta],&big,true));
          if (!DummyIter){
            //multiply the Wavefunction with the operators
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOpSing),*RightOp,&big,WF,*VtuvS,VopQS[0],SingFac[iquanta]);
            operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOpTrip),*RightOp,&big,WF,*VtuvT,VopQS[0],TripFac[iquanta]);
            //sum up the contributions
            ScaleAdd(1.0,*VtuvS,*Vtuv);
            ScaleAdd(1.0,*VtuvT,*Vtuv);
          }
         //do the same thing for Vutv
          boost::shared_ptr<Wavefunction> Vtvu (new Wavefunction(VQ[iquanta],&big,true));
          if (u!=v&&!DummyIter){
            ScaleAdd(1.0,*VtuvS,*Vtvu);
            ScaleAdd(-1.0,*VtuvT,*Vtvu);//the minus sign arises from the change of orbital order
          }//t!=u
          if (!T[iquanta].GetBasisOnly()){
            //store the operator
            bool EndOfArray = false;
            T[iquanta].ResetBuffer();
            while (!EndOfArray){
              boost::shared_ptr<Wavefunction> Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
              if (!EndOfArray){
                ScaleAdd(Kuv->element(t+NInternal,a),*Vtuv,*Va);
                ScaleAdd(Kvu->element(t+NInternal,a),*Vtvu,*Va);
                T[iquanta].ReplaceOpInBuffer(Va,dummy,a);
                //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                /*
                int arel=a+NInternal+NActive;
                if ((fabs(Kuv->element(t+NInternal,a))>1e-12)){
                  sprintf(msg,"\ncase 3: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%4.12lf twoBlock left a",u,t,v,arel,iquanta,SingFac[iquanta],TripFac[iquanta],Kuv->element(t+NInternal,a),DotProduct(*Vtuv,*Vtuv));pout << msg;
                }
                if ((u!=v)&&(fabs(Kvu->element(t+NInternal,a))>1e-12)) {
                  sprintf(msg,"\ncase 3: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)= int=%lf Norm=%4.12lf twoBlock left b",v,t,u,arel,iquanta,SingFac[iquanta],Kvu->element(t+NInternal,a),DotProduct(*Vtvu,*Vtvu));pout << msg;
                }
                 */
                //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
              }
            }//a
            //Send Vtuv around to other processes and receive theirs
            SendAroundVtuv_a(*Vtuv,*Vtvu,*Kuv,*Kvu,T[iquanta],t+NInternal,DummyIter);
            T[iquanta].AddToBasisOperator(Vtuv);
            T[iquanta].AddToBasisOperator(Vtvu);
          }//!BasisOnly
        }//iquanta
      }//v
    }//tu
    LS.Reset();
    if (iterCase==_FINAL_ITER_){
      //make sure that allprocesses undergo the same number of cycles in the loop
      LS.SetLocalDim(rightBlock->get_op_array(CRE_CRE).get_size());
      Dim_uv=LS.GetGlobalDim();
      //loop over a+a+ on the right block
      //for (int uv=0;uv<rightBlock->get_op_array(CRE_CRE).get_size();uv++){
      for (int uv=0;uv<Dim_uv;uv++){
        //determine whetehr this is a dummy iteration
        bool DummyIter = LS.DummyIter(uv);
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOpSing;
        boost::shared_ptr<SparseMatrix> RightOpTrip;
        if (!DummyIter){
          RightOpSing = rightBlock->get_op_array(CRE_CRE).get_local_element(uv)[0]->getworkingrepresentation(rightBlock);
          RightOpTrip = rightBlock->get_op_array(CRE_CRE).get_local_element(uv)[1]->getworkingrepresentation(rightBlock);
        }//regular iteration
        else{
          //if this is a dummy iteration generate dummy oeprators
          RightOpSing = boost::make_shared<Cre>();
          RightOpTrip = boost::make_shared<Cre>();
          RightOpSing->set_orbs().push_back(-1);
          RightOpSing->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpSing->resize_deltaQuantum(1);
          RightOpTrip->resize_deltaQuantum(1);
          RightOpSing->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(0),IrrepSpace(0));
          RightOpTrip->set_deltaQuantum(0)=SpinQuantum(2,SpinSpace(2),IrrepSpace(0));
        }//dummy iteration
        //get the orbital indices
        u = RightOpSing->get_orbs(0);
        v = RightOpSing->get_orbs(1);
        //the right operator quanta
        SpinQuantum ropQS = Transposeview(*RightOpSing).get_deltaQuantum(0);
        SpinQuantum ropQT = Transposeview(*RightOpTrip).get_deltaQuantum(0);
        //get the integrals
        boost::shared_ptr<Matrix> Kuv = IKJA.GetMatrix(u+NInternal,v+NInternal);
        boost::shared_ptr<Matrix> Kvu = IKJA.GetMatrix(v+NInternal,u+NInternal);
        //loop over a on the left block
        for (int tcount=0;tcount<leftBlock->get_op_array(CRE).get_size();tcount++){
          //get the right operator
          boost::shared_ptr<SparseMatrix> LeftOp = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
          //get the orbital index
          t = LeftOp->get_orbs(0);
          //check orbital indices for validity
          //if (!CheckAllowedOperator(t,u,v,dotIndex,iterCase)) continue;
          //the right operator quanta
          SpinQuantum lopQ = LeftOp->get_deltaQuantum(0);
          //the total Operator Quanta (is the same for the left singlet and triplet operator)
          vector<SpinQuantum> VopQS = ropQS + lopQ;
          //the possible quanta for the resulting wavefunction
          vector <SpinQuantum> VQ = WFQ + VopQS[0];
          //evaluate the prefactors
          vector<double> TripFac;
          vector<double> SingFac;
          for (i=0;i<VQ.size();i++){
            twoS_ = (double)VQ[i].get_s().getirrep();
            S_    = twoS_/2.0;
            fac   = sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
            SingFac.push_back(fac);
            TripFac.push_back(fac*sqrt(3.0));//on account of the transposition of the triplet operator the minus sign is cancelled
          }//resulting wavefunction quanta
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //generate the wavefunctions that holds the result of the multiplications
            boost::shared_ptr<Wavefunction> VtuvS (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtuvT (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> Vtuv (new Wavefunction(VQ[iquanta],&big,true));
            //multiply the Wavefunction with the operators
            if (!DummyIter){
              operatorfunctions::TensorMultiply(rightBlock,Transposeview(*RightOpSing),*LeftOp,&big,WF,*VtuvS,VopQS[0],SingFac[iquanta]);
              operatorfunctions::TensorMultiply(rightBlock,Transposeview(*RightOpTrip),*LeftOp,&big,WF,*VtuvT,VopQS[0],TripFac[iquanta]);
              //sum up the contributions
              ScaleAdd(1.0,*VtuvS,*Vtuv);
              ///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              // NOTE:
              // for some strange reason I dont't understand we have to add a minus
              // sign if the triplet A operator is on the right Block.
              //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              ScaleAdd(-1.0,*VtuvT,*Vtuv);
            }
            //do the same thing for Vutv
            boost::shared_ptr<Wavefunction> Vtvu (new Wavefunction(VQ[iquanta],&big,true));
            if (u!=v&&!DummyIter){
              ScaleAdd(1.0,*VtuvS,*Vtvu);
              ScaleAdd(1.0,*VtuvT,*Vtvu);//the minus sign arises from the change of orbital order
            }//t!=u
            if (!T[iquanta].GetBasisOnly()){
              //store the operator
              bool EndOfArray = false;
              T[iquanta].ResetBuffer();
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
                if (!EndOfArray){
                  ScaleAdd(Kuv->element(t+NInternal,a),*Vtuv,*Va);
                  ScaleAdd(Kvu->element(t+NInternal,a),*Vtvu,*Va);
                  T[iquanta].ReplaceOpInBuffer(Va,dummy,a);
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                  /*
                  int arel=a+NInternal+NActive;
                  if (fabs(Kuv->element(t+NInternal,a))>1e-12){
                    sprintf(msg,"\ncase 3: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%4.12lf twoBlock right a",u,t,v,arel,iquanta,SingFac[iquanta],TripFac[iquanta],Kuv->element(t+NInternal,a),DotProduct(*Vtuv,*Vtuv));pout << msg;
                  }
                  if ((u!=v)&&(fabs(Kvu->element(t+NInternal,a))>1e-12)) {
                    sprintf(msg,"\ncase 3: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)= int=%lf Norm=%4.12lf twoBlock right b",v,t,u,arel,iquanta,SingFac[iquanta],Kvu->element(t+NInternal,a),DotProduct(*Vtvu,*Vtvu));pout << msg;
                  }
                  */
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                }
              }//a
              //Send Vtuv around to other processes and receive theirs
              SendAroundVtuv_a(*Vtuv,*Vtvu,*Kuv,*Kvu,T[iquanta],t+NInternal,DummyIter);
            }//!BasisOnly
            T[iquanta].AddToBasisOperator(Vtvu);
            T[iquanta].AddToBasisOperator(Vtvu);
          }//iquanta
        }//v
      }//tu
      LS.Reset();
    }//final iteration
    //--------------------------------------------------------------------------
    //Case 2 and 3: a creation and an annihilation operator are on one block and  
    //        the second annihilation operator is on the other block
    //--------------------------------------------------------------------------
    //make sure that allprocesses undergo the same number of cycles in the loop
    LS.SetLocalDim(leftBlock->get_op_array(CRE_DES).get_size());
    int Dim_tv = LS.GetGlobalDim();
    double fac_;
    //loop over a+a on left Block
    //for (int tv=0;tv<leftBlock->get_op_array(CRE_DES).get_size();tv++){
    for (int tv=0;tv<Dim_tv;tv++){
      //determine whetehr this is a dummy iteration
      bool DummyIter = LS.DummyIter(tv);
      //get the left operator
      boost::shared_ptr<SparseMatrix> LeftOpSing;
      boost::shared_ptr<SparseMatrix> LeftOpTrip;
      if (!DummyIter){
        LeftOpSing = leftBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(leftBlock);
        LeftOpTrip = leftBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(leftBlock);
      }//regular iteration
      else{
        //if this is a dummy iteration generate dummy oeprators
        LeftOpSing = boost::make_shared<Cre>();
        LeftOpTrip = boost::make_shared<Cre>();
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpSing->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpTrip->set_orbs().push_back(-1);
        LeftOpSing->resize_deltaQuantum(1);
        LeftOpTrip->resize_deltaQuantum(1);
        LeftOpSing->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(0),IrrepSpace(0));
        LeftOpTrip->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(2),IrrepSpace(0));
      }
      //get the orbital indices
      t = LeftOpSing->get_orbs(0);
      v = LeftOpSing->get_orbs(1);
      //the left operator quanta
      SpinQuantum lopQS = LeftOpSing->get_deltaQuantum(0);
      SpinQuantum lopQT = LeftOpTrip->get_deltaQuantum(0);
      //get the integrals
      boost::shared_ptr<Matrix> Ktv = IKJA.GetMatrix(t+NInternal,v+NInternal);
      boost::shared_ptr<Matrix> Kvt = IKJA.GetMatrix(v+NInternal,t+NInternal);
      //loop over a+ on the right block
      for (int ucount=0;ucount<rightBlock->get_op_array(CRE).get_size();ucount++){
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOp = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
        //get the orbital index
        u = RightOp->get_orbs(0);
        //check orbital indices for validity
        //if (!CheckAllowedOperator(t,u,v,dotIndex,iterCase)) continue;
        //the right operator quanta
        SpinQuantum ropQ = Transposeview(*RightOp).get_deltaQuantum(0);
        //the total Operator Quanta (is the same for the left singlet and triplet operator)
        vector<SpinQuantum> VopQS = lopQS + ropQ;
        //the possible quanta for the resulting wavefunction
        vector <SpinQuantum> VQ = WFQ + VopQS[0];
        //evaluate the prefactors
        vector<double> TripFac;
        vector<double> SingFac;
        vector<double> SingFac_;
        for (i=0;i<VQ.size();i++){
          twoS_ = (double)VQ[i].get_s().getirrep();
          S_    = twoS_/2.0;
          fac   = -sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
          fac_  = sqrt(2.0) * sqrt((2.0*S_+1)/(2.0*S+1));
          SingFac.push_back(fac);
          SingFac_.push_back(fac_);
          TripFac.push_back(fac*sqrt(3.0));
        }//resulting wavefunction quanta
        //get the integrals
        boost::shared_ptr<Matrix> Ktu = IKJA.GetMatrix(t+NInternal,u+NInternal);
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //----------------------------
          //construct the operators Vtvu
          //----------------------------
          //generate the wavefunctions that holds the result of the multiplications
          boost::shared_ptr<Wavefunction> VtvuS (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VtvuT (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> Vtvu (new Wavefunction(VQ[iquanta],&big,true));
          if (!DummyIter){
            //multiply the Wavefunction with the operators
            operatorfunctions::TensorMultiply(leftBlock,*LeftOpSing,Transposeview(*RightOp),&big,WF,*VtvuS,VopQS[0],1.0);
            operatorfunctions::TensorMultiply(leftBlock,*LeftOpTrip,Transposeview(*RightOp),&big,WF,*VtvuT,VopQS[0],1.0);
            //sum up the contributions
            ScaleAdd(SingFac[iquanta],*VtvuS,*Vtvu);
            ScaleAdd(TripFac[iquanta],*VtvuT,*Vtvu);
          }
          //--------------------------
          //do the same thing for Vvtu
          //--------------------------
          boost::shared_ptr<Wavefunction> Vvtu (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VvtuS (new Wavefunction(VQ[iquanta],&big,true));
          boost::shared_ptr<Wavefunction> VvtuT (new Wavefunction(VQ[iquanta],&big,true));
          if (t!=v){
            //multiply the transposed operator
            if (!DummyIter){
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOpSing),Transposeview(*RightOp),&big,WF,*VvtuS,VopQS[0],1.0);
              operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOpTrip),Transposeview(*RightOp),&big,WF,*VvtuT,VopQS[0],1.0);
              //sum up the contributions
              ScaleAdd(SingFac[iquanta],*VvtuS,*Vvtu);
              ScaleAdd(-TripFac[iquanta],*VvtuT,*Vvtu);//account for transposition and add a minus sign
            }
          }//t!=v
          //--------------------------
          //Evaluate Vtuv 
          //--------------------------
          boost::shared_ptr<Wavefunction> Vtuv (new Wavefunction(VQ[iquanta],&big,true));
          if (!DummyIter){
            ScaleAdd(SingFac_[iquanta],*VtvuS,*Vtuv);
            if (t!=v){
              ScaleAdd(SingFac_[iquanta],*VvtuS,*Vtuv);
            }
          }//!DummyIter
          if (!T[iquanta].GetBasisOnly()){
            //store the operators
            T[iquanta].ResetBuffer();
            bool EndOfArray = false;
            while (!EndOfArray){
              boost::shared_ptr<Wavefunction> Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
              if (!EndOfArray){
                //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                /*int arel = a + NInternal+NActive;
                if ((fabs(Ktv->element(u+NInternal,a))>1e-12)){
                  sprintf(msg,"\ncase 2: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%4.12llf twoBlock left a",t,u,v,arel,iquanta,SingFac[iquanta],TripFac[iquanta],Ktv->element(u+NInternal,a),DotProduct(*Vtvu,*Vtvu));pout << msg;
                }
                if ((t!=v)&&(fabs(Kvt->element(u+NInternal,a))>1e-12)) {
                  sprintf(msg,"\ncase 2: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)= int=%lf Norm=%4.12lf twoBlock left b",v,u,t,arel,iquanta,SingFac[iquanta],Kvt->element(u+NInternal,a),DotProduct(*Vvtu,*Vvtu));pout << msg;
                }
                if ((fabs(Ktu->element(v+NInternal,a))>1e-12)){
                  sprintf(msg,"\ncase 1: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  int=%lf Norm=%lf Norm2=%4.12lf twoBlock left a",t,v,u,arel,iquanta,SingFac_[iquanta],Ktv->element(u+NInternal,a),DotProduct(*VtvuS,*VtvuS),DotProduct(*Vtuv,*Vtuv));pout << msg;
                }
                if ((t!=v)&&(fabs(Ktu->element(v+NInternal,a))>1e-12)) {
                  sprintf(msg,"\ncase 1: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  int=%lf Norm=%lf Norm2=%4.12lf twoBlock left b",v,t,u,arel,iquanta,SingFac_[iquanta],Kvt->element(u+NInternal,a),DotProduct(*VvtuS,*VvtuS),DotProduct(*Vtuv,*Vtuv));pout << msg;
                }
                */
                //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                ScaleAdd(Ktv->element(u+NInternal,a),*Vtvu,*Va);
                ScaleAdd(Kvt->element(u+NInternal,a),*Vvtu,*Va);
                ScaleAdd(Ktu->element(v+NInternal,a),*Vtuv,*Va);
                T[iquanta].ReplaceOpInBuffer(Va,dummy,a);
              }
            }//a
            //Send Vtuv around to other processes and receive theirs
            SendAroundVtuv_a(*Vtvu,*Vvtu,*Vtuv,*Ktv,*Kvt,*Ktu,T[iquanta],u+NInternal,v+NInternal,DummyIter);
          }//!BasisOnly
          T[iquanta].AddToBasisOperator(Vtvu);
          T[iquanta].AddToBasisOperator(Vvtu);
          T[iquanta].AddToBasisOperator(Vtuv);
        }//iquanta
      }//u
    }//tv
    LS.Reset();
    if (iterCase==_FINAL_ITER_){
      //make sure that allprocesses undergo the same number of cycles in the loop
      LS.SetLocalDim(rightBlock->get_op_array(CRE_DES).get_size());
      int Dim_tv = LS.GetGlobalDim();
      //loop over a+a on right Block
      //for (int tv=0;tv<rightBlock->get_op_array(CRE_DES).get_size();tv++){
      for (int tv=0;tv<Dim_tv;tv++){
        //determine whetehr this is a dummy iteration
        bool DummyIter = LS.DummyIter(tv);
        //get the right operator
        boost::shared_ptr<SparseMatrix> RightOpSing;
        boost::shared_ptr<SparseMatrix> RightOpTrip;
        if (!DummyIter){
          RightOpSing = rightBlock->get_op_array(CRE_DES).get_local_element(tv)[0]->getworkingrepresentation(rightBlock);
          RightOpTrip = rightBlock->get_op_array(CRE_DES).get_local_element(tv)[1]->getworkingrepresentation(rightBlock);
        }//regular iteration
        else{
          //if this is a dummy iteration generate dummy oeprators
          RightOpSing = boost::make_shared<Cre>();
          RightOpTrip = boost::make_shared<Cre>();
          RightOpSing->set_orbs().push_back(-1);
          RightOpSing->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpTrip->set_orbs().push_back(-1);
          RightOpSing->resize_deltaQuantum(1);
          RightOpTrip->resize_deltaQuantum(1);
          RightOpSing->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(0),IrrepSpace(0));
          RightOpTrip->set_deltaQuantum(0)=SpinQuantum(0,SpinSpace(2),IrrepSpace(0));
        }//dummy iteration
        //get the orbital indices
        t = RightOpSing->get_orbs(0);
        v = RightOpSing->get_orbs(1);
        //the left operator quanta
        SpinQuantum ropQS = RightOpSing->get_deltaQuantum(0);
        SpinQuantum ropQT = RightOpTrip->get_deltaQuantum(0);
        //get the integrals
        boost::shared_ptr<Matrix> Ktv = IKJA.GetMatrix(t+NInternal,v+NInternal);
        boost::shared_ptr<Matrix> Kvt = IKJA.GetMatrix(v+NInternal,t+NInternal);
        //loop over a+ on the right block
        for (int ucount=0;ucount<leftBlock->get_op_array(CRE).get_size();ucount++){
          //get the right operator
          boost::shared_ptr<SparseMatrix> LeftOp = leftBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(leftBlock);
          //get the orbital index
          u = LeftOp->get_orbs(0);
          //check orbital indices for validity
          //if (!CheckAllowedOperator(t,u,v,dotIndex,iterCase,true)) continue;
          //the right operator quanta
          SpinQuantum lopQ = Transposeview(*LeftOp).get_deltaQuantum(0);
          //the total Operator Quanta (is the same for the left singlet and triplet operator)
          vector<SpinQuantum> VopQS = ropQS + lopQ;
          //the possible quanta for the resulting wavefunction
          vector <SpinQuantum> VQ = WFQ + VopQS[0];
          //evaluate the prefactors
          vector<double> TripFac;
          vector<double> SingFac;
          vector<double> SingFac_;
          for (i=0;i<VQ.size();i++){
            twoS_ = (double)VQ[i].get_s().getirrep();
            S_    = twoS_/2.0;
            fac   = -sqrt((2.0*S_+1)/(2.0*(2.0*S+1)));
            fac_  = sqrt(2.0) * sqrt((2.0*S_+1)/(2.0*S+1));
            SingFac.push_back(fac);
            SingFac_.push_back(fac_);
            ///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // NOTE:
            // for some strange reason I dont't understand we have to add a minus
            // sign if the triplet E operator is on the right Block.
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            TripFac.push_back(-fac*sqrt(3.0));
          }//resulting wavefunction quanta
          //get the integrals
          boost::shared_ptr<Matrix> Ktu = IKJA.GetMatrix(t+NInternal,u+NInternal);
          for (int iquanta=0;iquanta<VQ.size();iquanta++){
            //----------------------------
            //construct the operators Vtvu
            //----------------------------
            //generate the wavefunctions that holds the result of the multiplications
            boost::shared_ptr<Wavefunction> VtvuS (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VtvuT (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> Vtvu (new Wavefunction(VQ[iquanta],&big,true));
            if (!DummyIter){
              //multiply the Wavefunction with the operators
              operatorfunctions::TensorMultiply(rightBlock,*RightOpSing,Transposeview(*LeftOp),&big,WF,*VtvuS,VopQS[0],1.0);
              operatorfunctions::TensorMultiply(rightBlock,*RightOpTrip,Transposeview(*LeftOp),&big,WF,*VtvuT,VopQS[0],1.0);
              //sum up the contributions
              ScaleAdd(SingFac[iquanta],*VtvuS,*Vtvu);
              ScaleAdd(TripFac[iquanta],*VtvuT,*Vtvu);
            }
            //--------------------------
            //do the same thing for Vvtu
            //--------------------------
            boost::shared_ptr<Wavefunction> Vvtu (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VvtuS (new Wavefunction(VQ[iquanta],&big,true));
            boost::shared_ptr<Wavefunction> VvtuT (new Wavefunction(VQ[iquanta],&big,true));
            if (t!=v){
              if (!DummyIter){
                operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOp),Transposeview(*RightOpSing),&big,WF,*VvtuS,VopQS[0],1.0);
                operatorfunctions::TensorMultiply(leftBlock,Transposeview(*LeftOp),Transposeview(*RightOpTrip),&big,WF,*VvtuT,VopQS[0],1.0);
                ScaleAdd(SingFac[iquanta],*VvtuS,*Vvtu);
                ScaleAdd(-TripFac[iquanta],*VvtuT,*Vvtu);//account for transposition and add a minus sign
              }
            }//t!=u
            //--------------------------
            //Evaluate Vtuv 
            //--------------------------
            boost::shared_ptr<Wavefunction> Vtuv (new Wavefunction(VQ[iquanta],&big,true));
            if (!DummyIter){
              ScaleAdd(SingFac_[iquanta],*VtvuS,*Vtuv);
              if (t!=v){
                ScaleAdd(SingFac_[iquanta],*VvtuS,*Vtuv);
              }
            }//!DummyIter
            if (!T[iquanta].GetBasisOnly()){
              //store the operators
              T[iquanta].ResetBuffer();
              bool EndOfArray = false;
              while (!EndOfArray){
                boost::shared_ptr<Wavefunction> Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
                if (!EndOfArray){
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                  /*
                  int arel = a + NInternal+NActive;
                  if ((fabs(Ktv->element(u+NInternal,a))>1e-12)){
                    sprintf(msg,"\ncase 2: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)=%lf int=%lf Norm=%4.12lf twoBlock right a",t,u,v,arel,iquanta,SingFac[iquanta],TripFac[iquanta],Ktv->element(u+NInternal,a),DotProduct(*Vtvu,*Vtvu));pout << msg;
                  }
                  if ((t!=v)&&(fabs(Kvt->element(u+NInternal,a))>1e-12)) {
                    sprintf(msg,"\ncase 2: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  fac(triplet)= int=%lf Norm=%4.12lf twoBlock right b",v,u,t,arel,iquanta,SingFac[iquanta],Kvt->element(u+NInternal,a),DotProduct(*Vvtu,*Vvtu));pout << msg;
                  }
                  if ((fabs(Ktu->element(v+NInternal,a))>1e-12)){
                    sprintf(msg,"\ncase 1: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  int=%lf Norm=%lf Norm2=%4.12lf twoBlock right a",t,v,u,arel,iquanta,SingFac_[iquanta],Ktu->element(v+NInternal,a),DotProduct(*VtvuS,*VtvuS),DotProduct(*Vtuv,*Vtuv));pout << msg;
                  }
                  if ((t!=v)&&(fabs(Ktu->element(v+NInternal,a))>1e-12)) {
                    sprintf(msg,"\ncase 1: (%i%i|%i%i) iquanta=%i fac(singlet)=%4.6e  int=%lf Norm=%lf Norm2=%4.12lf twoBlock right b",v,t,u,arel,iquanta,SingFac_[iquanta],Ktu->element(v+NInternal,a),DotProduct(*VvtuS,*VvtuS),DotProduct(*Vtuv,*Vtuv));pout << msg;
                  }                  
                   */
                  //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBU
                  ScaleAdd(Ktv->element(u+NInternal,a),*Vtvu,*Va);
                  ScaleAdd(Kvt->element(u+NInternal,a),*Vvtu,*Va);
                  ScaleAdd(Ktu->element(v+NInternal,a),*Vtuv,*Va);
                  T[iquanta].ReplaceOpInBuffer(Va,dummy,a);
                }
              }//a
              //Send Vtuv around to other processes and receive theirs
              SendAroundVtuv_a(*Vtvu,*Vvtu,*Vtuv,*Ktv,*Kvt,*Ktu,T[iquanta],u+NInternal,v+NInternal,DummyIter);
            }//BasisOnly
            T[iquanta].AddToBasisOperator(Vtvu);
            T[iquanta].AddToBasisOperator(Vvtu);
            T[iquanta].AddToBasisOperator(Vtuv);
          }//iquanta
        }//u
      }//tv
    }//final iteration
    //close the integral and operator container 
    IKJA.CloseFileRead();
    for (int i=0;i<T.size();i++){
      T[i].CloseFileRead();
    }
  }
  
  
  //============================================================================
  // Construct a+aa on a single block
  //============================================================================
  void ConstructCreDesDesSingleBlock(SpinBlock &big, SpinBlock *block, 
                                     ThreeIndOpArray &CDD, int iterCase){
    char msg[512];
    int t,u,v;
    SpinBlock *leftBlock =  big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //open the OperatorArray
    CDD.OpenFileWrite();
    
    for (int tu=0;tu<block->get_op_array(CRE_DES).get_size();tu++){
      //if this is a regular or the last iteration make sure we do not calculate 
      //the operator on multiple processors
      if (iterCase!=_INITIAL_ITER_ && SkipOperator(big.get_sites(),block->get_sites()[0]))continue;
      //get the first operator
      boost::shared_ptr<SparseMatrix> Etu = block->get_op_array(CRE_DES).get_local_element(tu)[0]->getworkingrepresentation(block);
      //get the orbital labels
      t=Etu->get_orbs(0);
      u=Etu->get_orbs(1);
      //get the operator quanta
      SpinQuantum opQ1 = Etu->get_deltaQuantum(0);
      for (int vcount=0;vcount<block->get_op_array(CRE).get_size();vcount++){
        //get the operator
        boost::shared_ptr<SparseMatrix> av = block->get_op_array(CRE).get_local_element(vcount)[0]->getworkingrepresentation(block);
        //get the orbital label
        v = av->get_orbs(0);
        //get the operator quanta
        SpinQuantum opQ2 = Transposeview(*av).get_deltaQuantum(0);
        SpinQuantum opQ = (opQ1 + opQ2)[0];
        //---------------
        //Construct Vtuv
        //---------------
        //generate the product operator
        boost::shared_ptr<Cre> Op (new Cre);
        Op->set_orbs().push_back(t);
        Op->set_orbs().push_back(u);
        Op->set_orbs().push_back(v);
        Op->set_initialised() = true;
        Op->set_fermion() = true;
        Op->resize_deltaQuantum(1);
        Op->set_deltaQuantum(0) = opQ;
        Op->allocate(block->get_stateInfo());
        //multiply the two operators
        operatorfunctions::Product(block,*Etu,Transposeview(*av),*Op,sqrt(2.0));
        boost::shared_ptr<Cre> OpRep (new Cre);
        switch(iterCase){
          case _INITIAL_ITER_:
            //store the operator
            CDD.AppendOperator(Op,t,u,v);
            break;
          case _REGULAR_ITER_:
          case _FINAL_ITER_:
            //build the operator representation in the larger block
            OpRep->set_orbs() = Op->get_orbs();
            OpRep->set_initialised()=true;
            OpRep->resize_deltaQuantum(1);
            OpRep->set_deltaQuantum(0)=opQ;
            OpRep->allocate(leftBlock->get_stateInfo());
            operatorfunctions::TensorTrace(block,*Op,leftBlock,&(leftBlock->get_stateInfo()),*OpRep);
            //store the oeprator
            CDD.AppendOperator(OpRep,t,u,v);
            break;
        }//itercase
        //------------------------
        //do the same for the Vtvu
        //------------------------
        if (t!=u){
          //the resulting operator
          boost::shared_ptr<Cre> Op_ (new Cre);
          Op_->set_orbs().push_back(u);
          Op_->set_orbs().push_back(t);
          Op_->set_orbs().push_back(v);
          Op_->set_initialised()=true;
          Op_->set_fermion()=true;
          Op_->resize_deltaQuantum(1);
          Op_->set_deltaQuantum(0)=opQ;
          Op_->allocate(block->get_stateInfo());
          //multiply the two operators
          operatorfunctions::Product(block,Transposeview(*Etu),Transposeview(*av),*Op_,sqrt(2.0));
          boost::shared_ptr<Cre> OpRep_;
          switch (iterCase){
            case _INITIAL_ITER_:
              //store the operator
              CDD.AppendOperator(Op_,u,t,v);
              break;
            case _REGULAR_ITER_:
            case _FINAL_ITER_:
              //build the operator representation in the larger block
              OpRep_->set_orbs() = Op_->get_orbs();
              OpRep_->set_initialised()=true;
              OpRep_->resize_deltaQuantum(1);
              OpRep_->set_deltaQuantum(0)=opQ;
              OpRep_->allocate(leftBlock->get_stateInfo());
              operatorfunctions::TensorTrace(block,*Op_,leftBlock,&(leftBlock->get_stateInfo()),*OpRep_);
              //store the operator
              CDD.AppendOperator(OpRep_,u,t,v);
              break;
          }//itercase
        }//u!=v
      }//t
    }//uv
    //close the operator array
    CDD.CloseFileWrite();
  }
  

  
  
  //============================================================================
  //Construct the one-electron part of the V(a) operator
  //============================================================================
  void ConstructDes(SpinBlock &big, Matrix &heff_, Wavefunction &WF, vector<WavefunctionArray> &T,
                    IntegralContainer &IKJA, int *OrbWin){
    
    char msg[512];
    //get the orbital spaces and stuff
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    int NInternal = i1-i0+1;
    int NActive = big.get_sites().size();
    int NExternal = a1-a0+1;
    int NOcc = NInternal+NActive;
    int OrbDim = NInternal+NActive+NExternal;
    int t,u,v,i,a,tu;
    int dummy;

    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    double twoS = (double) WFQ.get_s().getirrep();
    double S    = twoS/2.0;
    double twoS_,S_,fac;
    
    //open the operator container 
    for (int i=0;i<T.size();i++){
      T[i].OpenFileRead();
    }
     
    //note: here we do not need to synchronize the loop and broadcast operators
    //because every process has the full set of CRE operators
    //-------------
    //leftBlock
    //-------------
    for (int tcount=0;tcount<leftBlock->get_op_array(CRE).get_size();tcount++){
      //get the operator
      boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE).get_local_element(tcount)[0]->getworkingrepresentation(leftBlock);
      //get the orbital label
      t = op->get_orbs(0);
      //get the operator quanta
      SpinQuantum opQT = Transposeview(*op).get_deltaQuantum(0);
      //the possible quanta for the resulting wavefunction
      vector<SpinQuantum> VQT = opQT + WFQ;
      //the prefactors
      vector<double> Fac;
      //evaluate the prefactors
      for (int iquanta=0;iquanta<VQT.size();iquanta++){
        twoS_ = (double) VQT[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt(2.0*S_+1.0)/sqrt(2.0*S+1.0);
        Fac.push_back(fac);
      }
      for (int iquanta=0;iquanta<VQT.size();iquanta++){
        //the wavefunction that holds the result of the multiplication 
        boost::shared_ptr<Wavefunction> VtT(new Wavefunction(VQT[iquanta],&big,true));
        //do the multiplication with the wavefunction
        operatorfunctions::TensorMultiply(leftBlock,Transposeview(*op),&big,WF,*VtT,opQT,Fac[iquanta]);
        if (!T[iquanta].GetBasisOnly()){
          //store the operator
          T[iquanta].ResetBuffer();
          bool EndOfArray = false;
          while (!EndOfArray){
            boost::shared_ptr<Wavefunction> Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
            if (!EndOfArray){
              ScaleAdd(heff_.element(a+NOcc,t+NInternal),*VtT,*Va);
              T[iquanta].ReplaceOpInBuffer(Va,dummy,a);
            }
          }//a
        }//!BasisOnly
        T[iquanta].AddToBasisOperator(VtT);
      }//iquanta
    }//t on the left side
    //-------------
    //rightBlock
    //-------------
    for (int ucount=0;ucount<rightBlock->get_op_array(CRE).get_size();ucount++){
      //get the oeprator
      boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(ucount)[0]->getworkingrepresentation(rightBlock);
      //get the orbital label
      u = rop->get_orbs(0);
      //the operator quanta
      SpinQuantum ropQT = Transposeview(*rop).get_deltaQuantum(0);
      //the possible quanta for the resulting wavefunction
      vector<SpinQuantum> VQT = ropQT + WFQ;
      //the prefactors
      vector<double> Fac;
      for (int iquanta=0;iquanta<VQT.size();iquanta++){
        twoS_ = (double) VQT[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt(2.0*S_+1.0)/sqrt(2.0*S+1.0);
        Fac.push_back(fac);
      }//iquanta
      for (int iquanta=0;iquanta<VQT.size();iquanta++){
        //the wavefunction that holds the result of the multiplication
        boost::shared_ptr<Wavefunction> VuT (new Wavefunction(VQT[iquanta],&big,true));
        //multiply the operator with the wavefunction
        operatorfunctions::TensorMultiply(rightBlock,Transposeview(*rop),&big,WF,*VuT,ropQT,Fac[iquanta]);
        if (!T[iquanta].GetBasisOnly()){
          //store the operator
          T[iquanta].ResetBuffer();
          bool EndOfArray = false;
          while (!EndOfArray){
            boost::shared_ptr<Wavefunction> Va = T[iquanta].GetOpFromBuffer(dummy,a,EndOfArray);
            if (!EndOfArray){
              ScaleAdd(heff_.element(a+NOcc,u+NInternal),*VuT,*Va);
              T[iquanta].ReplaceOpInBuffer(Va,dummy,a);
            }
          }//a
        }//!BasisOnly
        T[iquanta].AddToBasisOperator(VuT);
      }//iquanta
    }//u on the right side
    //close operator container 
    for (int i=0;i<T.size();i++){
      T[i].CloseFileRead();
    }
    
  }
  
  //============================================================================
  //Construct the contributions to the perturber functions V(a) where all indices
  //are on one block
  //============================================================================
  void ConstructVaSingleBlock(SpinBlock &big, Wavefunction &WF, ThreeIndOpArray &CDD,
                              ThreeIndOpArray &CDD_,vector<WavefunctionArray> &Ta, 
                              IntegralContainer &IKJA,NEVPT2Info &Info){
    
    char msg[512];
    //get the orbital spaces and stuff
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
    int NOcc = NInternal + NActive;
    int OrbDim = NInternal+NActive+NExternal;
    int t=-1,u=-1,v=-1,i,a,tu;
    int t_,u_,v_;
    int dummy;
    SpinQuantum WFQ = WF.get_deltaQuantum(0);
    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    double twoS = (double) WFQ.get_s().getirrep();
    double S    = twoS/2.0;
    double twoS_,S_,fac;
    char BaseName[512];
    Info.getBaseName(BaseName);
    //Note: we have to synchronize the loops over (t,u,v) and broadcast the
    //operators, because each processor needs every possible operator O(tuv).
    //The division between processes is made in index a
    //---------------------------------------------------
    //first construct all contributions on the left Block
    //---------------------------------------------------
    //open the integral file
    IKJA.OpenFileRead();
    boost::shared_ptr<Matrix> Ktv = IKJA.GetMatrix(0,0);
    //Open the operator arrays and prepare the buffer
    CDD.OpenFileRead();
    CDD.ResetBuffer();
    //open the V(a) buffer
    for (i=0;i<Ta.size();i++){
      Ta[i].OpenFileRead();
    }
    //initialize the auxiliary indices
    t_ = -1;
    u_ = -1;
    v_ = -1;
    //start the loop
    bool EndOfOuterArray = false;
    bool NeedTensorTrace = false;
    //make sure that all processes undergo the same number of cycles in the loop
    LoopSynchronizer LS(CDD.GetLength());//this might lead to trouble if we have stored the same operator twice;
    //get the dimension of the loop
    int Dim_tuv = LS.GetGlobalDim();
    int ituv=0;
    while (ituv<Dim_tuv){
      //determine if this is just a dummy iteration
      bool DummyIter = LS.DummyIter(ituv);
      //get the operator
      boost::shared_ptr<Cre> Otuv;
      if (!DummyIter) {
        Otuv = CDD.GetOpFromBuffer(t,u,v,NeedTensorTrace,EndOfOuterArray);
      }
      else{
        //if this is a dummy iteration, generate an empty dummy operator
        NeedTensorTrace=false;
        Otuv = boost::make_shared<Cre> ();
        Otuv->set_initialised() = true;
        Otuv->resize_deltaQuantum(1);
        Otuv->set_deltaQuantum(0) = SpinQuantum(-1,SpinSpace(1),IrrepSpace(0));
      }
      //if we don't have an operator and this is not a dummy iteration, 
      //then something went wrong
      if (EndOfOuterArray&&!DummyIter){
        sprintf(msg,"\nERROR NEVPT2: in V(a): loop over (t,u,v) has not been properly synchronized!!!");
        pout << msg;
        DummyIter = true;
      }
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
      if (((t!=t_)||(v!=v_))  && !DummyIter){
         Ktv = IKJA.GetMatrix(t+NInternal,v+NInternal);
      }
      //the possible quanta
      SpinQuantum opQ = Otuv->get_deltaQuantum(0);
      vector<SpinQuantum> VQ = WFQ + opQ;

      for (int iquanta=0;iquanta<VQ.size();iquanta++){

        //evaluate the prefactor
        twoS_ = (double) VQ[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt((2.0*S_+1)/(2.0*S+1));

        //generate V(tuv) = O(tuv)|psi>
        boost::shared_ptr<Wavefunction> Vtuv(new Wavefunction(VQ[iquanta],&big,true));
        if (!DummyIter) operatorfunctions::TensorMultiply(leftBlock,*Otuv,&big,WF,*Vtuv,opQ,fac);

        //add it to the V(a) functions
        Ta[iquanta].ResetBuffer();
        bool EndOfInnerArray = false;
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Va = Ta[iquanta].GetOpFromBuffer(dummy,a,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(Ktv->element(u+NInternal,a),*Vtuv,*Va);
            Ta[iquanta].ReplaceOpInBuffer(Va,dummy,a);
          }
        }//a
        //Send Vtuv around to other processes and receive theirs
        SendAroundVtuv_a(*Vtuv,*Ktv,Ta[iquanta],u+NInternal,DummyIter);
      }//iquanta
      t_ = t;
      u_ = u;
      v_ = v;
      ituv++;
    }//tuv
    //close the Operator array
    CDD.CloseFileRead();
    //---------------------------------------------------
    //then construct all contributions on the right Block
    //---------------------------------------------------
    //Open the operator array and prepare the buffer
    CDD_.OpenFileRead();
    CDD_.ResetBuffer();
    //initialize the auxiliary indices
    t_ = -1;
    u_ = -1;
    v_ = -1;
    //start the loop
    EndOfOuterArray = false;
    //make sure that all processes undergo the same number of cycles in the loop
    LoopSynchronizer LS_(CDD_.GetLength());//this might lead to trouble if we have stored the same operator twice;
    //get the dimension of the loop
    Dim_tuv = LS_.GetGlobalDim();
    ituv=0;
    while (ituv<Dim_tuv){
      //determine if this is just a dummy iteration
      bool DummyIter = LS_.DummyIter(ituv);
      //get the operator
      boost::shared_ptr<Cre> Otuv;
      if (!DummyIter) {
        Otuv = CDD_.GetOpFromBuffer(t,u,v,NeedTensorTrace,EndOfOuterArray);
      }
      else{
        //if this is a dummy iteration, generate an empty dummy operator
        NeedTensorTrace=false;
        Otuv = boost::make_shared<Cre> ();
        Otuv->set_initialised() = true;
        Otuv->resize_deltaQuantum(1);
        Otuv->set_deltaQuantum(0) = SpinQuantum(-1,SpinSpace(1),IrrepSpace(0));
      }
      //if we don't have an operator and this is not a dummy iteration, 
      //then something went wrong
      if (EndOfOuterArray&&!DummyIter){
        sprintf(msg,"\nERROR NEVPT2: in V(a): loop over (t,u,v) has not been properly synchronized!!!");
        mpi_message(msg);
        DummyIter = true;
      }
      //if necessary, get the integral matrix
      if (((t!=t_)||(v!=v_)) && !DummyIter){
         Ktv = IKJA.GetMatrix(t+NInternal,v+NInternal);
      }
      //the possible quanta
      SpinQuantum opQ = Otuv->get_deltaQuantum(0);
      vector<SpinQuantum> VQ = WFQ + opQ;

      for (int iquanta=0;iquanta<VQ.size();iquanta++){
        //evaluate the prefactor
        twoS_ = (double) VQ[iquanta].get_s().getirrep();
        S_ = twoS_/2.0;
        fac = sqrt((2.0*S_+1)/(2.0*S+1));

        //generate V(tuv) = O(tuv)|psi>
        boost::shared_ptr<Wavefunction> Vtuv(new Wavefunction(VQ[iquanta],&big,true));
        if (!DummyIter) operatorfunctions::TensorMultiply(rightBlock,*Otuv,&big,WF,*Vtuv,opQ,fac);        

        //add it to the V(a) functions
        Ta[iquanta].ResetBuffer();
        bool EndOfInnerArray = false;
        while (!EndOfInnerArray){
          boost::shared_ptr<Wavefunction> Va = Ta[iquanta].GetOpFromBuffer(dummy,a,EndOfInnerArray);
          if (!EndOfInnerArray){
            ScaleAdd(Ktv->element(u+NInternal,a),*Vtuv,*Va);
            Ta[iquanta].ReplaceOpInBuffer(Va,dummy,a);
          }
        }//a
        //Send Vtuv around to other processes and receive theirs
        SendAroundVtuv_a(*Vtuv,*Ktv,Ta[iquanta],u+NInternal,DummyIter);
      }//iquanta
      t_ = t;
      u_ = u;
      v_ = v;
      ituv++;
    }//tuv
    //close the Operator array
    CDD_.CloseFileRead();
    //close the integral container
    IKJA.CloseFileRead();
    //close the V(a) array
    for (i=0;i<Ta.size();i++){
      Ta[i].CloseFileRead();
    }
    
  }

  
  //============================================================================
  //the driver for the construction of the a+a+a|psi> operators
  //============================================================================
  void ConstructCreDesDes(SpinBlock &big, vector<Wavefunction> &WF,ThreeIndOpArray &CDD, 
                          IntegralContainer &IKJA, SweepParams &sweepParams, NEVPT2Info &Info){
    
    char msg[512];
    int OrbWin[6];
    //get the orbital spaces and stuff
    char BaseName[512];
    Info.getBaseName(BaseName);
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
    int iSweep = Info.GetNevSweep();
    int MaxIter = Info.GetMaxBlockIter(iSweep);
    int iter = sweepParams.get_block_iter();
    int t,u,v,i,a,tu;
    int dummy;

    SpinBlock *leftBlock  = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    SpinBlock *dotBlock = leftBlock->get_rightBlock();
    int MaxCore = Info.GetMaxCore();
    int MaxM = dmrginp.sweep_state_schedule()[dmrginp.sweep_qstate_schedule().size()];
    double StartTime = GetTime();
    
    
    //--------------------------------
    //build the two-electron operators
    //--------------------------------
    if (iter==0){
      //first iteration
      ConstructCreDesDesSingleBlock(big,big.get_leftBlock(),CDD,_INITIAL_ITER_);
    }
    else if ((iSweep==1)&&(iter==MaxIter)){
      //last iteration
      ConstructCreDesDesSingleBlock(big,big.get_leftBlock()->get_rightBlock(),CDD,_FINAL_ITER_);
      ConstructCreDesDes_1_2_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CDD,IKJA,OrbWin);
      ConstructCreDesDes_2_1_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CDD,IKJA,OrbWin);
    }
    else{
      //a regular iteration
      ConstructCreDesDesSingleBlock(big,big.get_leftBlock()->get_rightBlock(),CDD,_REGULAR_ITER_);
      ConstructCreDesDes_1_2_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CDD,IKJA,OrbWin);
      ConstructCreDesDes_2_1_0(big,big.get_leftBlock()->get_leftBlock(),big.get_leftBlock()->get_rightBlock(),CDD,IKJA,OrbWin);
    }
    
    //--------------------------------------------------------------------------
    //when in the middle of the second sweep, actually build the V(a)|psi> 
    //functions and evaluate the energy contribution
    //--------------------------------------------------------------------------
    if ((iSweep==1)&&(iter==MaxIter)){
      //get the effective one-electron hamiltonian
      Matrix heff_;
      Info.GetH(1,heff_);
      //receive the OperatorArray from disk. Note: CDD_ here was CDD in function sweep_nevpt2.C:"do_one"
      ThreeIndOpArray CDD_;
      sprintf(msg,"%s.CDD.tmp",BaseName);
      CDD_.Retrieve(msg);
      for (int iroot=0;iroot<WF.size();iroot++){
        //------------------------------------
        //first fill the V(a) arrays with zero
        //------------------------------------
        SpinQuantum WFQ = WF[iroot].get_deltaQuantum(0);
        vector<SpinQuantum> VQ = WFQ + SpinQuantum(-1,SpinSpace(1),IrrepSpace(0));
        //generate and initialize the array
        vector<WavefunctionArray> Va;
        Va.resize(VQ.size());
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          //divide the external orbital space between the processes
          int a_loc_start,a_loc_end;
          int loc_NExternal;
          PALDivideLoop(a_loc_start,a_loc_end,0,NExternal);
          loc_NExternal = a_loc_end-a_loc_start;
          
          sprintf(msg,"%s.Va.%i.%i.tmp",BaseName,iroot,iquanta);
          Va[iquanta].Initialize(NExternal,msg,mpi_world_size());
          Va[iquanta].CalcBufferSize(MaxCore,MaxM);
          //open the array
          Va[iquanta].OpenFileWrite();
          Va[iquanta].SetQuanta(VQ[iquanta]);
          for (a=a_loc_start;a<a_loc_end;a++){
            boost::shared_ptr<Wavefunction> Starter (new Wavefunction(VQ[iquanta],&big,true));
            Va[iquanta].AppendOperator(Starter,0,a);
          }//a
          //close the array
          Va[iquanta].CloseFileWrite();
          //create the empty basis state
          boost::shared_ptr<Wavefunction> EmptyBasis (new Wavefunction(VQ[iquanta],&big,true));
          Va[iquanta].SetBasisOperator(EmptyBasis);
        }//iquanta
        //------------------------
        //then build the operators
        //------------------------
        ConstructVaSingleBlock(big,WF[iroot],CDD,CDD_,Va,IKJA,Info);
        ConstructCreDesDesTwoBlocks(big,WF[iroot],Va,IKJA,_FINAL_ITER_,OrbWin);
        ConstructDes(big,heff_,WF[iroot],Va,IKJA,OrbWin);
        //--------------------------------
        //evaluate the energy contribution
        //--------------------------------
        double Ea = V_a(Va,big,WF[iroot],Info,iroot);
        Info.SetE(iroot,6,Ea);
        //--------
        //clean up
        //--------
        for (int iquanta=0;iquanta<VQ.size();iquanta++){
          Va[iquanta].Clear();
        }//iquanta
      }//iroot
      CDD.Clear();
      CDD_.Clear();
    }//final iteration
    //take the time and add it to the NEVPT2 time
    double FinishTime=GetTime();
    double time = FinishTime-StartTime;
    Info.AddTime(4,time);
    Info.AddTimePerClass(6,time);
  }
  
  
  
}





