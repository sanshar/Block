

#include "nevpt2_operators.h"
#include "operatorfunctions.h"
#include "ripdm.h"
#include "nevpt2_mpi.h"
#include "nevpt2_util.h"
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif


namespace SpinAdapted {
  
  //=================================
  // gather the local operator arrays
  //=================================
  void gather_terminal_operators(WavefunctionArray &T) {
#ifndef SERIAL
      WavefunctionArray tmp_recv;
      mpi::communicator world;
      if (world.size() == 1)
        return;
      if (world.rank()==0) {
        for (int i = 1; i < world.size(); ++i) {
          world.recv(i, i, tmp_recv);
          T.AppendArray(tmp_recv);
        }
      } else {
        world.send(0, world.rank(), T);
      }
#endif
  }
  
  //==========================================
  // synchronize the parallel 2pdm's or 4pdm's
  //==========================================
  void SynchronizePDM(array_4d<double> &pdm){
#ifndef SERIAL
    int i;
    int k,l,m,n;
    int norbs = pdm.dim1();
    char msg[512];
    mpi::communicator world;
    if (world.size()!=0){
      //----------------------
      //gather the information
      //----------------------
      if (world.rank()==0){
        array_4d<double> tmp;
        for (i=1;i<world.size();i++){
          world.recv(i,i,tmp);
          for (k=0;k<norbs;k++){
            for (l=0;l<norbs;l++){
              for (m=0;m<norbs;m++){
                for (n=0;n<norbs;n++){
                  if (fabs(tmp(k,l,m,n))>1e-12){
                    if (fabs(pdm(k,l,m,n))>1e-12){
                      //sprintf(msg,"\nERROR during synchronizing pdm's: trying to overwrite existing data!!!");
                      continue;
                    }
                    else{
                      pdm(k,l,m,n) = tmp(k,l,m,n);
                    }
                  }
                }//n
              }//k
            }//l
          }//k
        }//processes
      }//master
      else{
        world.send(0,world.rank(),pdm);
      }//slave
      //----------------------
      //spread the information
      //----------------------
      mpi::broadcast(world,pdm,0);
    }//size>1?
#endif
  }

  //================================
  // synchronize the parallel 3pdm's 
  //================================
  void SynchronizePDM(array_6d &pdm){
#ifndef SERIAL
    int i;
    int k,l,m,n,o,p;
    int norbs = pdm.get_dim1();
    char msg[512];
    double master, slave;
    mpi::communicator world;
    int foo=world.rank();
    if (world.size()!=0){
      //----------------------
      //gather the information
      //----------------------
      if (world.rank()==0){
        array_6d tmp(norbs);
        for (i=1;i<world.size();i++){
          world.recv(i,i,tmp);
          for (k=0;k<norbs;k++){
            for (l=0;l<norbs;l++){
              for (m=0;m<norbs;m++){
                for (n=0;n<norbs;n++){
                  for (o=0;o<norbs;o++){
                    for (p=0;p<norbs;p++){
                      if (fabs(tmp(k,l,m,n,o,p))>1e-12){
                        if (fabs(pdm(k,l,m,n,o,p))>1e-12){
                          //sprintf(msg,"\nERROR during synchronizing 3pdm's: trying to overwrite existing data!!!");pout << msg;
                          continue;
                        }
                        else{
                          pdm(k,l,m,n,o,p) = tmp(k,l,m,n,o,p);
                        }
                      }
                    }//p
                  }//o
                }//n
              }//k
            }//l
          }//k
        }//processes
      }//master
      else{
        world.send(0,world.rank(),pdm);
      }//slave        
      //----------------------
      //spread the information
      //----------------------
      mpi::broadcast(world,pdm,0);
    }//size>1?
#endif
  }

  //=====================================
  // Sum the parallel 2pdm's or 4pdm's up
  //=====================================
  void SumPDM(array_4d<double> &pdm){
#ifndef SERIAL
    int i;
    int k,l,m,n;
    int norbs = pdm.dim1();
    char msg[512];
    mpi::communicator world;
    if (world.size()!=0){
      //----------------------
      //gather the information
      //----------------------
      if (world.rank()==0){
        array_4d<double> tmp;
        for (i=1;i<world.size();i++){
          world.recv(i,i,tmp);
          for (k=0;k<norbs;k++){
            for (l=0;l<norbs;l++){
              for (m=0;m<norbs;m++){
                for (n=0;n<norbs;n++){
                  pdm(k,l,m,n) += tmp(k,l,m,n);
                }//n
              }//k
            }//l
          }//k
        }//processes
      }//master
      else{
        world.send(0,world.rank(),pdm);
      }//slave
      //----------------------
      //spread the information
      //----------------------
      mpi::broadcast(world,pdm,0);
    }//size>1?
#endif
  }

  //===========================
  // sum the parallel 3pdm's up
  //===========================
  void SumPDM(array_6d &pdm){
#ifndef SERIAL
    int i;
    int k,l,m,n,o,p;
    int norbs = pdm.get_dim1();
    char msg[512];
    mpi::communicator world;
    if (world.size()!=0){
      //----------------------
      //gather the information
      //----------------------
      if (world.rank()==0){
        array_6d tmp;
        for (i=1;i<world.size();i++){
          world.recv(i,i,tmp);
          for (k=0;k<norbs;k++){
            for (l=0;l<norbs;l++){
              for (m=0;m<norbs;m++){
                for (n=0;n<norbs;n++){
                  for (o=0;o<norbs;o++){
                    for (p=0;p<norbs;p++){
                      pdm(k,l,m,n,o,p) += tmp(k,l,m,n,o,p);
                    }//p
                  }//o
                }//n
              }//k
            }//l
          }//k
        }//processes
      }//master
      else{
        world.send(0,world.rank(),pdm);
      }//slave
      //----------------------
      //spread the information
      //----------------------
      mpi::broadcast(world,pdm,0);
    }//size>1?
#endif
  }

  void PAL_DivideLoop(int &start_loc, int &stop_loc, int start_global, int stop_global){
#ifndef SERIAL
    mpi::communicator world;
    if (world.size()==1){
      start_loc = start_global;
      stop_loc   = stop_global;
    }
    else{
      int loopsize = stop_global - start_global + 1;
      int numproc = world.size();
      int quot = loopsize / numproc;
      int iproc = world.rank();
      start_loc = iproc * quot;
      stop_loc = (iproc+1) * quot;
      if (world.rank()==numproc-1){
        stop_loc = stop_global;
      }
    }
#else
    start_loc = start_global;
    stop_loc   = stop_global;
#endif
  }
  
  void SumPalWF(Wavefunction &WF){
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
      Wavefunction tmp;
      for (int i=1;i<world.size();i++){
        world.recv(i,i,tmp);
        ScaleAdd(1.0,tmp,WF);
      }//i
    }//master
    else{
      world.send(0,world.rank(),WF);
    }//slave
#endif
  }
  
  
  //=========================================================
  // the function that generates the operators <psi|E(i,j)|M>
  //=========================================================
  void generate_terminal_replacement_operators(const SpinBlock& big, Wavefunction& WF,
          WavefunctionArray& T) {

    SpinBlock* leftBlock = big.get_leftBlock();
    SpinBlock* rightBlock = big.get_rightBlock();
    int ij = 0;
    //open the storage device
    T.OpenFileWrite();
    //----------------------------------------------------
    // case 1: i and j are on the left side of the lattice
    //----------------------------------------------------
    for (ij = 0; ij < leftBlock->get_op_array(CRE_DES).get_size(); ij++) {
      //get the singlet B operator
      boost::shared_ptr<SparseMatrix>sing_op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);
      //get the orbital labels
      int i = sing_op->get_orbs(0);
      int j = sing_op->get_orbs(1);
      // generate a wavefunction that will carry the tensor product of the B operator
      // and the solution wavefunction
      SpinQuantum dQ=WF.get_deltaQuantum(0);
      SpinQuantum opQ(sing_op->get_deltaQuantum(0).particleNumber,SpinSpace(0),IrrepSpace(0));
      boost::shared_ptr<Wavefunction> B_x_psi(new Wavefunction(dQ, &big, true));
      //tensor multiply the B operator with the solution wavefunction
      operatorfunctions::TensorMultiply(leftBlock, Transposeview(*sing_op), &big, WF, *B_x_psi, opQ, sqrt(2.0));
      //store the product
      T.AppendOperator(B_x_psi, i, j);
      if (i!=j){
        //do the same for the transposed operator
        boost::shared_ptr<Wavefunction> B_x_psi_T(new Wavefunction(dQ, &big, true));
        operatorfunctions::TensorMultiply(leftBlock, *sing_op, &big, WF, *B_x_psi_T, opQ, sqrt(2.0));
        T.AppendOperator(B_x_psi_T, j, i);
      }
    }//(i,j)
    //----------------------------------------------------
    // case 2: i and j are on the right side of the lattice
    //----------------------------------------------------
    for (ij = 0; ij < rightBlock->get_op_array(CRE_DES).get_size(); ij++) {
      //get the singlet B operator
      boost::shared_ptr<SparseMatrix>sing_op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(rightBlock);
      //get the orbital labels
      int i = sing_op->get_orbs(0);
      int j = sing_op->get_orbs(1);
      // generate a wavefunction that will carry the tensor product of the B operator
      // and the solution wavefunction
      SpinQuantum dQ=WF.get_deltaQuantum(0);
      SpinQuantum opQ(sing_op->get_deltaQuantum(0).particleNumber,SpinSpace(0),IrrepSpace(0));
      boost::shared_ptr<Wavefunction> B_x_psi(new Wavefunction(dQ, &big, true));
      //tensor multiply the B operator with the solution wavefunction
      operatorfunctions::TensorMultiply(rightBlock, Transposeview(*sing_op), &big, WF, *B_x_psi, opQ, sqrt(2.0));
      //store the product
      T.AppendOperator(B_x_psi, i, j);
      if (i!=j){
        //do the same for the transposed operator
        boost::shared_ptr<Wavefunction> B_x_psi_T(new Wavefunction(dQ, &big, true));
        operatorfunctions::TensorMultiply(rightBlock, *sing_op, &big, WF, *B_x_psi_T, opQ, sqrt(2.0));
        T.AppendOperator(B_x_psi_T, j, i);
      }
    }//(i,j)
    //-----------------------------------------------------------------
    // case 3: i and j are on the left and the right side, respectively
    //-----------------------------------------------------------------
    int li;
    int ri;
    int li_start,li_stop;
    //divide the loop over the parallel processes
    PAL_DivideLoop(li_start,li_stop,0,leftBlock->get_op_array(CRE).get_size());
    //loop over local elements on the left side
    for (li = li_start; li < li_stop; li++) {
      //get the left operator
      boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(li)[0]->getworkingrepresentation(leftBlock);
      //get the left orbital index
      int i = lop->get_orbs(0);
      //loop over local elements on the right side
      for (ri = 0; ri < rightBlock->get_op_array(CRE).get_size(); ri++) {
        //get the right operator
        boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(ri)[0]->getworkingrepresentation(rightBlock);
        //get the right orbital index
        int j = rop->get_orbs(0);
        //generate the possible spin quanta for the product operator
        vector<SpinQuantum> opQ = lop->get_deltaQuantum(0) + rop->get_deltaQuantum(0);
        // generate a wavefunction that will carry the tensor product of the two
        // operators and the solution wavefunction
        SpinQuantum dQ=WF.get_deltaQuantum(0);
        boost::shared_ptr<Wavefunction> B_x_psi(new Wavefunction(dQ, &big, true));
        //tensor multiply the two operators with the solution wavefunction
        operatorfunctions::TensorMultiply(leftBlock, *lop, Transposeview(*rop), &big, WF, *B_x_psi, opQ[0], sqrt(2.0));
        //store the product
        T.AppendOperator(B_x_psi, j, i);
        //do the same for the transposed matrix 
        boost::shared_ptr<Wavefunction> B_x_psi_T(new Wavefunction(dQ, &big, true));
        operatorfunctions::TensorMultiply(leftBlock, Transposeview(*lop), *rop, &big, WF, *B_x_psi_T, opQ[0], sqrt(2.0));
        T.AppendOperator(B_x_psi_T, i, j);
      }//right local index
    }//left local index
    //close the storage device
    T.CloseFileWrite();
    //gather the terminal vectors from all processors
    T.Synchronize();
#ifndef SERIAL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  }

  //============================================================================
  // the function that produces the two-particle reduced density matrices from 
  // the terminal operators and the one-particle reduced density matrices: 
  // Gamma(ij,kl) = <psi|E(i,k)|M><M|E(j,l)|psi> - delta(k,j) * Gamma(i,l)
  // the one-pdm is a by product of the procedure
  //============================================================================
  void generate_twopdm_RI_a(SpinBlock& big, WavefunctionArray& T, Wavefunction& WF,
                           array_4d<double>& twopdm, array_4d<double>& twopdm_C, Matrix& onepdm) {
    int ik, jl;
    int i, j, k, l;
    int num_orbs = big.size();
    int num_op = T.GetNumberOfOperators();
    char msg[512];
    std::pair <int, int> leftorbs, rightorbs;
    boost::shared_ptr<Wavefunction> leftop;
    boost::shared_ptr<Wavefunction> rightop;
    boost::shared_ptr<Wavefunction> one_op;
    double val = 0.0;
    
    //Initialize the density matrices
    for (i=0;i<num_orbs;i++){
      for (j=0;j<num_orbs;j++){
        for (k=0;k<num_orbs;k++){
          for (l=0;l<num_orbs;l++){
            twopdm(i,j,k,l) = 0.0;
            twopdm_C(i,j,k,l) = 0.0;
          }
        }
      }
    }
    onepdm.ReSize(num_orbs,num_orbs);
    for (i=0;i<num_orbs;i++){
      for (j=0;j<num_orbs;j++){
        onepdm.element(i,j) = 0.0;
      }
    }
        
    //open the operator storage device
    T.OpenFileRead();
    
    //now actually start to generate the one and two pdm's
    for (i = 0; i < num_orbs; i++) {
      for (k = 0; k < num_orbs; k++) {
        //check if the operator has been stored previously
        if (T.GetAdress(i, k).first == -1) {
          sprintf(msg, "\n\t\t\tERROR: Operator <psi|E(%i,%i)|M> was requested but is unavailable!!! 2PDM\n", i, k);
          pout << msg;
          continue;
        }
        //get the terminal operator <psi|E(i,k)|M>
        //leftop = T.GetOperator(i, k);
        leftop = T.GetOpPal(i, k);
        for (j = 0; j < num_orbs; j++) {
          for (l = 0; l < num_orbs; l++) {
            //check if the operator has been stored previously
            if (T.GetAdress(l, j).first == -1) {
              sprintf(msg, "\n\t\t\tERROR: Operator <psi|E(%i,%i)|M> was requested but is unavailable!!! 2PDM b\n", j, l);
              pout << msg;
              continue;
            }
            //get the right operator
            //rightop = T.GetOperator(l, j);
            rightop = T.GetOpPal(l, j);
            // gamma(ij,kl) = <psi|E(i,k)|M><M|E(j,l)|psi>
            val = DotProduct(*leftop, *rightop);
            //add the leading part to the two-pdm
            twopdm(i,j,k,l) += val;
            twopdm_C(i,j,k,l) += val;
          }//l
        }//j
        //obtain the one-particle reduced density matrix elements
        val = DotProduct(*leftop,WF);
        //add it to the appropriate places in the reduced twopdm
        for (j = 0; j < num_orbs; j++) {
          twopdm(i, j, j, k) -= val; // Gamma(ij,lk) -= delta(j,l)*gamma(i,k)
        }//j
        //add it to the one pdm
        onepdm.element(i,k) = val;
      }//k
    }//i
    //close the operator storage device
    T.CloseFileRead();
  };

  
  //============================================================================
  // the function that produces  three particle density matrices from terminal
  // and linking replacement operators. 
  //============================================================================
  void generate_threepdm_RI(SpinBlock &big, WavefunctionArray& T, array_6d& threepdm,
                            array_6d& threepdm_C, const array_4d<double>& twopdm, const Matrix& onepdm, Wavefunction& WF) {
    int ij,kl,mn;
    int i,j,k,l,m,n;
    int norbs = big.get_sites().size();
    double val = 0.0;
    int k0,k1;
    char msg[512];
    SpinBlock *leftBlock = big.get_leftBlock();
    SpinBlock *rightBlock = big.get_rightBlock();
    
    //Initialize the three pdm's
    threepdm.initialize();
    threepdm_C.initialize();
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //Timer tmain;
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    
    //open the operator storage device
    T.OpenFileRead();
    
    //start building the threepdm
    for (i=0;i<norbs;i++) {
      for (j=0;j<norbs;j++) {
        //check if the operator is available
        if (T.GetAdress(i, j).first == -1) {
          sprintf(msg, "\n\t\t\tERROR: Operator <psi|E(%i,%i)|M> was requested but is unavailable!!!\n", i, j);
          pout << msg;
          continue;
        }
        //get the left terminal operator
        //boost::shared_ptr<Wavefunction> left_op = T.GetOperator(i, j);
        boost::shared_ptr<Wavefunction> left_op = T.GetOpPal(i, j);
        for (m=0;m<norbs;m++) {
          for (n=0;n<norbs;n++) {
            //check if the operator is available
            if (T.GetAdress(n,m).first == -1) {
              sprintf(msg, "\n\t\t\tERROR: Operator <psi|E(%i,%i)|M> was requested but is unavailable!!!\n", n, m);
              pout << msg;
              continue;
            }
            //get the right terminal operator
            //boost::shared_ptr<Wavefunction> right_op = T.GetOperator(n, m);
            boost::shared_ptr<Wavefunction> right_op = T.GetOpPal(n, m);
            //----------------------------------------------------
            // case 1: k and l are on the left side of the lattice
            //----------------------------------------------------
            //loop over local array of operators
            for (kl=0;kl<leftBlock->get_op_array(CRE_DES).get_size();kl++) {
              //get the singlet replacement operator
              boost::shared_ptr<SparseMatrix> sing_op = leftBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(leftBlock);
              //get the orbital indices
              k = sing_op->get_orbs(0);
              l = sing_op->get_orbs(1);
              //generate the wavefunction that holds the result of the multiplication of 
              //the right terminal operator and the linking operator
              SpinQuantum dQ=WF.get_deltaQuantum(0);
              SpinQuantum opQ=sing_op->get_deltaQuantum(0);
              boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
              //multiply the linking operator with the right terminal operator
              operatorfunctions::TensorMultiply(leftBlock, *sing_op, &big, *right_op, *tmp, opQ, sqrt(2.0));
              //form the dot product of the resulting wavefunction and the left terminal operator
              val = DotProduct(*left_op, *tmp);
              //add the result to the three particle density matrix (it is the leading three-body term)
              threepdm(i,k,m,j,l,n) += val;
              threepdm_C(i,k,m,j,l,n) += val;
              //the same has to be done for the transposed linking operator
              if (k!=l){
                tmp->Clear();
                operatorfunctions::TensorMultiply(leftBlock, Transposeview(*sing_op), &big, *right_op, *tmp, opQ, sqrt(2.0));
                val = DotProduct(*left_op, *tmp);
                threepdm(i,l,m,j,k,n) += val;
                threepdm_C(i,l,m,j,k,n) += val;
              }
            }//kl
            //-----------------------------------------------------
            // case 2: k and l are on the right side of the lattice
            //-----------------------------------------------------
            //loop over local array of operators
            for (kl=0;kl<rightBlock->get_op_array(CRE_DES).get_size();kl++) {
              //get the singlet replacement operator
              boost::shared_ptr<SparseMatrix> sing_op = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);
              //get the orbital indices
              k = sing_op->get_orbs(0);
              l = sing_op->get_orbs(1);
              //generate the wavefunction that holds the result of the multiplication of 
              //the right terminal operator and the linking operator
              SpinQuantum dQ=WF.get_deltaQuantum(0);
              SpinQuantum opQ=sing_op->get_deltaQuantum(0);
              boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
              //multiply the linking operator with the right terminal operator
              operatorfunctions::TensorMultiply(rightBlock, *sing_op, &big, *right_op, *tmp, opQ, sqrt(2.0));
              //form the dot product of the resulting wavefunction and the left terminal operator
              val = DotProduct(*left_op, *tmp);
              //add the result to the three particle density matrix (it is the leading three-body term)
              threepdm(i,k,m,j,l,n) += val;
              threepdm_C(i,k,m,j,l,n) += val;
              //the same has to be done for the transposed linking operator
              if (k!=l){
                tmp->Clear();
                operatorfunctions::TensorMultiply(rightBlock, Transposeview(*sing_op), &big, *right_op, *tmp, opQ, sqrt(2.0));
                val = DotProduct(*left_op, *tmp);
                threepdm(i,l,m,j,k,n) += val;
                threepdm_C(i,l,m,j,k,n) += val;
              }
            }//kl
            //-----------------------------------------------------------------
            // case 3: k and l are on the left and the right side, respectively
            //-----------------------------------------------------------------
            PAL_DivideLoop(k0,k1,0,leftBlock->get_op_array(CRE).get_size());
            //loop over local array of operators
            for (int kcount=k0;kcount < k1;kcount++){
              //get the left operator
              boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(kcount)[0]->getworkingrepresentation(leftBlock);
              //get the left orbital index
              k = lop->get_orbs(0);
              //loop over local array of operators
              for (int lcount = 0; lcount < rightBlock->get_op_array(CRE).get_size(); lcount++) {
                //get the right operator
                boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(lcount)[0]->getworkingrepresentation(rightBlock);
                //get the right orbital index
                l = rop->get_orbs(0);
                //generate the wavefunction that holds the result of the multiplication of 
                //the right terminal operator and the linking operator
                SpinQuantum dQ=WF.get_deltaQuantum(0);
                SpinQuantum opQ(0,SpinSpace(0),IrrepSpace(0));
                boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
                //multiply the linking operator with the right terminal operator
                operatorfunctions::TensorMultiply(leftBlock, *lop, Transposeview(*rop), &big, *right_op, *tmp, opQ, sqrt(2.0));
                //form the dot product of the resulting wavefunction and the left terminal operator
                val = DotProduct(*left_op, *tmp);
                //add the result to the three particle density matrix (it is the leading three-body term)
                threepdm(i,k,m,j,l,n) += val;
                threepdm_C(i,k,m,j,l,n) += val;
                //the same has to be done for the transposed linking operator
                if (k!=l){
                  tmp->Clear();
                  operatorfunctions::TensorMultiply(leftBlock, Transposeview(*lop), *rop, &big, *right_op, *tmp, opQ, sqrt(2.0));
                  val = DotProduct(*left_op, *tmp);
                  threepdm(i,l,m,j,k,n) += val;
                  threepdm_C(i,l,m,j,k,n) += val;
                }//k!=l
              }//l
            }//k
          }//n
        }//m
      }//j
    }//i
    //synchronize the data
    SynchronizePDM(threepdm);
#ifndef SERIAL
    MPI_Barrier(MPI_COMM_WORLD);
#endif    
    SynchronizePDM(threepdm_C);
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //int main = tmain.elapsedwalltime();
    //Timer trest;
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    
    //--------------------------------------------------
    //add the one- and two particle density matrix parts
    //--------------------------------------------------
    for (i=0;i<norbs;i++) {
      for (j=0;j<norbs;j++) {
        for (m=0;m<norbs;m++) {
          for (n=0;n<norbs;n++) {
            //add the two-particle density matrix parts
            for (k=0;k<norbs;k++) {
              for (l=0;l<norbs;l++){
                if (m==l) threepdm(i,k,m,j,l,n) -= twopdm(i,k,j,n);
                if (m==j) threepdm(i,k,m,j,l,n) -= twopdm(i,k,n,l);
                if (k==j) threepdm(i,k,m,j,l,n) -= twopdm(i,m,l,n);
              }//l
            }//k
            //add the one-particle density matrix parts
            //add the one-pdm part only once for each pair of orbitals (m,n))
            threepdm(i,m,n,m,n,j) -= onepdm.element(i,j);
          }//n
        }//m
      }//j
    }//i
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //int rest = trest.elapsedwalltime();
    //sprintf(msg,"\nmain: %i s",main);pout<<msg;
    //sprintf(msg,"\nrest: %i s",rest);pout<<msg;
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    
    //close the operator storage device
    T.CloseFileRead();
  }
  
  
  //============================================================================
  // the function that generates the four pdm using the RI approximation from 
  // terminal and linking operators. 
  //============================================================================
  void generate_fourpdm_RI(SpinBlock &big, WavefunctionArray& T,array_4d<double>& fourpdm, Wavefunction &WF){
    int i,j,k,l,m,n,o,p;
    int aij,aji,akl,alk;
    int kl;
    int k0,k1;
    int norbs = big.get_sites().size();
    int norbs2 = norbs * norbs;
    double val = 0.0;
    char msg[512];
    SpinBlock *leftBlock=big.get_leftBlock();
    SpinBlock *rightBlock=big.get_rightBlock();
    
    //Initialize the fourpdm
    for (i=0;i<norbs2;i++){
      for (j=0;j<norbs2;j++){
        for (k=0;k<norbs2;k++){
          for (l=0;l<norbs2;l++){
            fourpdm(i,j,k,l) = 0.0;
          }
        }
      }
    }
    
    //open the operator storage device
    T.OpenFileRead();
    
    //---------------------------------------------------
    // generate intermediates: <psi|E(i,j)|M><M|E(k,l)|L>
    //---------------------------------------------------
    // create the storage device for the intermediate operators
    // note: since now we store four index quantities we address them by pairs of orbital index pairs
    WavefunctionArray T_IM;
    sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/IntermediateOperators.tmp.");
    T_IM.SetFileName(msg);
    T_IM.ResizeOrbSpace(norbs2,true);
    T_IM.OpenFileWrite();

      for (i=0;i<norbs;i++){
      for (j=0;j<norbs;j++){
        //define the address of pair (i,j)
        aij = i * norbs + j;
        aji = j * norbs + i;
        //check if the operator is available
        if (T.GetAdress(j, i).first == -1) {
          sprintf(msg, "\n\t\t\tERROR: Operator <psi|E(%i,%i)|M> was requested but is unavailable!!!\n", j, i);
          pout << msg;
          continue;
        }
        //get the left terminal operator
        boost::shared_ptr<Wavefunction> terminal = T.GetOpPal(j,i);
        //boost::shared_ptr<Wavefunction> terminal = T.GetOperator(j,i);
        //loop over locally stored operators
        //---------------------------------------------------------
        // case 1: both indices are on the left side of the lattice
        //---------------------------------------------------------
        for (kl=0;kl<leftBlock->get_op_array(CRE_DES).get_size();kl++){
          //get the local linking operator
          boost::shared_ptr<SparseMatrix> sing_op = leftBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(leftBlock);
          //get the orbital indices
          k = sing_op->get_orbs(0);
          l = sing_op->get_orbs(1);
          //define the addresses of pairs(k,l) and (l,k)
          akl = k * norbs + l;
          alk = l * norbs + k;          
          //generate the wavefunction that holds the result of the multiplication of 
          //the left terminal operator and the linking operator
          SpinQuantum dQ=WF.get_deltaQuantum(0);
          SpinQuantum opQ=sing_op->get_deltaQuantum(0);
          boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
          //multiply the linking operator with the right terminal operator
          operatorfunctions::TensorMultiply(leftBlock, *sing_op, &big, *terminal, *tmp, opQ, sqrt(2.0));
          //store the resulting wavefunction
          T_IM.AppendOperator(tmp,aji,alk);
          //do the same for the transposed operator
          if (k!=l){
            boost::shared_ptr<Wavefunction> tmp_T(new Wavefunction(dQ, &big, true));
            operatorfunctions::TensorMultiply(leftBlock, Transposeview(*sing_op), &big, *terminal, *tmp_T, opQ, sqrt(2.0));
            T_IM.AppendOperator(tmp_T,aji,akl);
          }//k!=l
        }//kl
        //----------------------------------------------------------
        // case 2: both indices are on the right side of the lattice
        //----------------------------------------------------------
        for (kl=0;kl<rightBlock->get_op_array(CRE_DES).get_size();kl++){
          //get the local linking operator
          boost::shared_ptr<SparseMatrix> sing_op = rightBlock->get_op_array(CRE_DES).get_local_element(kl)[0]->getworkingrepresentation(rightBlock);
          //get the orbital indices
          k = sing_op->get_orbs(0);
          l = sing_op->get_orbs(1);
          //define the addresses of pairs(k,l) and (l,k)
          akl = k * norbs + l;
          alk = l * norbs + k;
          //generate the wavefunction that holds the result of the multiplication of 
          //the left terminal operator and the linking operator
          SpinQuantum dQ=WF.get_deltaQuantum(0);
          SpinQuantum opQ=sing_op->get_deltaQuantum(0);
          boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
          //multiply the linking operator with the right terminal operator
          operatorfunctions::TensorMultiply(rightBlock, *sing_op, &big, *terminal, *tmp, opQ, sqrt(2.0));
          //store the resulting wavefunction
          T_IM.AppendOperator(tmp,aji,alk);
          //do the same for the transposed operator
          if (k!=l){
            boost::shared_ptr<Wavefunction> tmp_T(new Wavefunction(dQ, &big, true));
            operatorfunctions::TensorMultiply(rightBlock, Transposeview(*sing_op), &big, *terminal, *tmp_T, opQ, sqrt(2.0));
            T_IM.AppendOperator(tmp_T,aji,akl);
          }//k!=l
        }//kl
        //-------------------------------------------------------------------------
        // case 3: the two indices are on the left and the right side, respectively
        //-------------------------------------------------------------------------
        PAL_DivideLoop(k0,k1,0,leftBlock->get_op_array(CRE).get_size());
        for (int kcount=k0;kcount<k1;kcount++){
          //get the left operator
          boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(kcount)[0]->getworkingrepresentation(leftBlock);
          //get the orbital index
          k = lop->get_orbs(0);
          for (int lcount=0;lcount<rightBlock->get_op_array(CRE).get_size();lcount++){
            //get the right operator
            boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(lcount)[0]->getworkingrepresentation(rightBlock);
            //get the orbital index
            l = rop->get_orbs(0);
            //define the addresses of pairs(k,l) and (l,k)
            akl = k * norbs + l;
            alk = l * norbs + k;
            //generate the wavefunction that holds the result of the multiplication of 
            //the left terminal operator and the linking operator
            SpinQuantum dQ=WF.get_deltaQuantum(0);
            SpinQuantum opQ(0,SpinSpace(0),IrrepSpace(0));
            boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
            //multiply the linking operator with the right terminal operator
            operatorfunctions::TensorMultiply(leftBlock, *lop, Transposeview(*rop), &big, *terminal, *tmp, opQ, sqrt(2.0));
            //store the resulting wavefunction
            T_IM.AppendOperator(tmp,aji,alk);
            //do the same for the transposed operator
            if (k!=l){
              boost::shared_ptr<Wavefunction> tmp_T(new Wavefunction(dQ, &big, true));
              operatorfunctions::TensorMultiply(leftBlock, Transposeview(*lop), *rop, &big, *terminal, *tmp_T, opQ, sqrt(2.0));
              T_IM.AppendOperator(tmp_T,aji,akl);
            }//k!=l
          }//l
        }//k
      }//j
    }//i
    T_IM.CloseFileWrite();
    //gather all locally stored intermediate operators
    T_IM.Synchronize();
    T.CloseFileRead();
    T_IM.OpenFileRead();
    //--------------------------------
    // combine intermediates to 4pdm's
    //--------------------------------
    int amn,anm,aop,apo;
    for (i=0;i<norbs;i++){
      for (k=0;k<norbs;k++){
        for (j=0;j<norbs;j++){
          for (l=0;l<norbs;l++){
            //get the addresses
            aij = i * norbs + j;
            akl = k * norbs + l;
            //get the left operator
            //boost::shared_ptr<Wavefunction> Tikjl = T_IM.GetTerminalOperator(aij,akl);
            boost::shared_ptr<Wavefunction> Tikjl = T_IM.GetOpPal(aij,akl);
            for (m=0;m<norbs;m++){
              for (o=0;o<norbs;o++){
                for (n=0;n<norbs;n++){
                  for (p=0;p<norbs;p++){
                    //get the addresses
                    amn = m * norbs + n;
                    anm = n * norbs + m;
                    aop = o * norbs + p;
                    apo = p * norbs + o;
                    //get the right operator
                    //boost::shared_ptr<Wavefunction> Tpnom = T_IM.GetTerminalOperator(apo,anm);
                    boost::shared_ptr<Wavefunction> Tpnom = T_IM.GetOpPal(apo,anm);
                    //construct the 4pdms
                    val = DotProduct(*Tikjl,*Tpnom);
                    fourpdm(aij,akl,amn,aop) = val;
                  }//p
                }//n
              }//o
            }//m
          }//l
        }//j
      }//k
    }//i
    //Clean up
    T_IM.CloseFileRead();
    T_IM.Clear();
  }
  
  
  
  void generate_A_RI(SpinBlock &big, WavefunctionArray& T,array_6d& A, array_6d& A_, Wavefunction &WF,array_6d &threepdm_C){
    int i,j,k,l,m,n,o,p;
    int aij,aji,akl,alk;
    int kl;
    int norbs = big.get_sites().size();
    int norbs2 = norbs * norbs;
    double val = 0.0;
    char msg[512];
    SpinBlock *leftBlock=big.get_leftBlock();
    SpinBlock *rightBlock=big.get_rightBlock();
    
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //Timer ttot,tim;
    //int tot,im,main,rest;
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
#ifndef SERIAL
    mpi::communicator world;
#endif

    //Initialize the A matrices
    A.initialize();
    A_.initialize();
    
    //-----------------------------------------------------------------------------
    // generate intermediates: T'(l,m)  = (ij|kl) * <L|E(i,j)|M><M|E(k,m)|psi>
    //                         T''(l,m) = (ij|kl) * <L|E(i,j)|M><M|E(m,k)|psi>
    // where it is summed over i,j,k.
    //----------------------------------------------------------------------------
    int amk,akm,akj,ain,alm;
    int ij;
    double integral;
    int icount,jcount;
    double fac = 1.0;
    int i0,i1;
    SpinQuantum dQ=WF.get_deltaQuantum(0);
    T.OpenFileRead();
    //get the integral array
    TwoElectronArray IJKL = *(big.get_twoInt());
    // create, initialize and open the storage devices for the intermediate operators
    WavefunctionArray T_,T__,T_IM___;
    sprintf(msg,"T_.tmp");
    T_.Initialize(norbs,msg,mpi_world_size());
    T_.OpenFileWrite();

    sprintf(msg,"T__.tmp");
    T__.Initialize(norbs,msg,mpi_world_size());
    T__.OpenFileWrite();
    for (l=0;l<norbs;l++){
      for (m=0;m<norbs;m++){
        boost::shared_ptr<Wavefunction> intermediate1(new Wavefunction(dQ,&big,true));
        boost::shared_ptr<Wavefunction> intermediate2(new Wavefunction(dQ,&big,true));
        for (k=0;k<norbs;k++){
          //get the terminal operators
          //boost::shared_ptr<Wavefunction> terminal1 = T.GetOperator(m,k);
          //boost::shared_ptr<Wavefunction> terminal2 = T.GetOperator(k,m);
          boost::shared_ptr<Wavefunction> terminal1 = T.GetOpPal(m,k);
          boost::shared_ptr<Wavefunction> terminal2 = T.GetOpPal(k,m);
          //loop over locally stored operators
          //------------------------------------------------------------------
          // case 1: both indices of (i,j) are on the left side of the lattice
          //------------------------------------------------------------------
#ifndef SERIAL
          //if the leftblock is just a dotblock, this contribution will be added
          //by each process. Hence we have to divide it by the number of procs
          if (leftBlock->get_sites().size()==1){
            double numproc  = (double) world.size();
            fac = 1/numproc;
          }
          else fac = 1.0;
#endif  
          for (ij=0;ij<leftBlock->get_op_array(CRE_DES).get_size();ij++){
            //get the operator
            boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(leftBlock);
            //get the orbital indices
            i = op->get_orbs(0);
            j = op->get_orbs(1);
            //generate the wavefunction that holds the result of the multiplication of 
            //the terminal operator and the linking operator
            SpinQuantum opQ=op->get_deltaQuantum(0);
            boost::shared_ptr<Wavefunction> tmp1(new Wavefunction(dQ, &big, true));
            boost::shared_ptr<Wavefunction> tmp_T1(new Wavefunction(dQ, &big, true));
            boost::shared_ptr<Wavefunction> tmp2(new Wavefunction(dQ, &big, true));
            boost::shared_ptr<Wavefunction> tmp_T2(new Wavefunction(dQ, &big, true));
            //multiply the linking operator with the terminal operator
            operatorfunctions::TensorMultiply(leftBlock, *op, &big, *terminal1, *tmp1, opQ, sqrt(2.0));
            if (i!=j)operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op), &big, *terminal1, *tmp_T1, opQ, sqrt(2.0));
            operatorfunctions::TensorMultiply(leftBlock, *op, &big, *terminal2, *tmp2, opQ, sqrt(2.0));
            if (i!=j)operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op), &big, *terminal2, *tmp_T2, opQ, sqrt(2.0));
            //sum the product of operator strings and integrals
            ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp1,*intermediate1);
            if (i!=j)ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp_T1,*intermediate1);
            ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp2,*intermediate2);
            if (i!=j)ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp_T2,*intermediate2);
          }//ij
          //-------------------------------------------------------------------
          // case 2: both indices of (i,j) are on the right side of the lattice
          //-------------------------------------------------------------------
#ifndef SERIAL
          //if the leftblock is just a dotblock, this contribution will be added
          //by each process. Hence we have to divide it by the number of procs
          if (rightBlock->get_sites().size()==1){
            double numproc  = (double) world.size();
            fac = 1/numproc;
          }
          else fac = 1.0;
#endif  
          for (ij=0;ij<rightBlock->get_op_array(CRE_DES).get_size();ij++){
            //get the operator
            boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(CRE_DES).get_local_element(ij)[0]->getworkingrepresentation(rightBlock);
            //get the orbital indices
            i = op->get_orbs(0);
            j = op->get_orbs(1);
            //generate the wavefunction that holds the result of the multiplication of 
            //the terminal operator and the linking operator
            SpinQuantum opQ=op->get_deltaQuantum(0);
            boost::shared_ptr<Wavefunction> tmp1(new Wavefunction(dQ, &big, true));
            boost::shared_ptr<Wavefunction> tmp_T1(new Wavefunction(dQ, &big, true));
            boost::shared_ptr<Wavefunction> tmp2(new Wavefunction(dQ, &big, true));
            boost::shared_ptr<Wavefunction> tmp_T2(new Wavefunction(dQ, &big, true));
            //multiply the linking operator with the terminal operator
            operatorfunctions::TensorMultiply(rightBlock, *op, &big, *terminal1, *tmp1, opQ, sqrt(2.0));
            if (i!=j)operatorfunctions::TensorMultiply(rightBlock, Transposeview(*op), &big, *terminal1, *tmp_T1, opQ, sqrt(2.0));
            operatorfunctions::TensorMultiply(rightBlock, *op, &big, *terminal2, *tmp2, opQ, sqrt(2.0));
            if (i!=j)operatorfunctions::TensorMultiply(rightBlock, Transposeview(*op), &big, *terminal2, *tmp_T2, opQ, sqrt(2.0));
            //sum the product of operator strings and integrals
            ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp1,*intermediate1);
            if (i!=j)ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp_T1,*intermediate1);
            ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp2,*intermediate2);
            if (i!=j)ScaleAdd(fac*IJKL(2*i,2*k,2*j,2*l),*tmp_T2,*intermediate2);
          }//ij
          //---------------------------------------------------------------------------
          // case 3: indices i and j are on the both sides of the lattice, respectively
          //---------------------------------------------------------------------------
          PAL_DivideLoop(i0,i1,0,leftBlock->get_op_array(CRE).get_size());
          for (icount=i0;icount<i1;icount++){
            //get the left operator
            boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(icount)[0]->getworkingrepresentation(leftBlock);
            //get the orbital index i
            i = lop->get_orbs(0);
            for (jcount=0;jcount<rightBlock->get_op_array(CRE).get_size();jcount++){
              //get the operator
              boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(jcount)[0]->getworkingrepresentation(rightBlock);
              //get the orbital index j
              j = rop->get_orbs(0);
              //generate the wavefunction that holds the result of the multiplication of 
              //the terminal operator and the linking operator
              SpinQuantum opQ(0,SpinSpace(0),IrrepSpace(0));
              boost::shared_ptr<Wavefunction> tmp1(new Wavefunction(dQ, &big, true));
              boost::shared_ptr<Wavefunction> tmp_T1(new Wavefunction(dQ, &big, true));
              boost::shared_ptr<Wavefunction> tmp2(new Wavefunction(dQ, &big, true));
              boost::shared_ptr<Wavefunction> tmp_T2(new Wavefunction(dQ, &big, true));
              //multiply the linking operator with the terminal operator
              operatorfunctions::TensorMultiply(leftBlock,*lop,Transposeview(*rop),&big,*terminal1,*tmp1,opQ,sqrt(2.0));
              if (i!=j)operatorfunctions::TensorMultiply(leftBlock,Transposeview(*lop),*rop,&big,*terminal1,*tmp_T1,opQ,sqrt(2.0));
              operatorfunctions::TensorMultiply(leftBlock,*lop,Transposeview(*rop),&big,*terminal2,*tmp2,opQ,sqrt(2.0));
              if (i!=j)operatorfunctions::TensorMultiply(leftBlock,Transposeview(*lop),*rop,&big,*terminal2,*tmp_T2,opQ,sqrt(2.0));
              //sum the product of operator strings and integrals
              ScaleAdd(IJKL(2*i,2*k,2*j,2*l),*tmp1,*intermediate1);
              if (i!=j)ScaleAdd(IJKL(2*i,2*k,2*j,2*l),*tmp_T1,*intermediate1);
              ScaleAdd(IJKL(2*i,2*k,2*j,2*l),*tmp2,*intermediate2);
              if (i!=j)ScaleAdd(IJKL(2*i,2*k,2*j,2*l),*tmp_T2,*intermediate2);
            }//jcount
          }//icount
        }//k
        SumPalWF(*intermediate1);
        SumPalWF(*intermediate2);
#ifndef SERIAL
        if (world.rank()==0){
          T_.AppendOperator(intermediate1,l,m);
          T__.AppendOperator(intermediate2,l,m);
        }
#else
        T_.AppendOperator(intermediate1,l,m);
        T__.AppendOperator(intermediate2,l,m);
#endif
      }//m
    }//l
    //close the storage devices
    T_.CloseFileWrite();
    T__.CloseFileWrite();
    //gather all locally stored intermediate operators
    T_.Synchronize();
    T__.Synchronize();
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //im = tim.elapsedwalltime();
    //Timer tmain;
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //-------------------------------------------------
    // generate the auxiliary matrix A
    //-------------------------------------------------
    int i_,j_,k_,jj_;
    int ak_i_,aj_j,aik;
    int j_count;
    int j0,j1;
    double val1,val2,val3,val4;
    //open the storage devices for reading
    T_.OpenFileRead();
    T__.OpenFileRead();
    for (k_=0;k_<norbs;k_++){
      for (i_=0;i_<norbs;i_++){
        //get the address of pair(k_,i_)
        ak_i_ = k_ * norbs + i_;
        //get the terminal operator
        //boost::shared_ptr<Wavefunction> terminal = T.GetTerminalOperator(k_,i_);
        boost::shared_ptr<Wavefunction> terminal = T.GetOpPal(k_,i_);
        for (i=0;i<norbs;i++){
          for (k=0;k<norbs;k++){
            //boost::shared_ptr<Wavefunction> T_ik  = T_.GetOperator(i,k);
            //boost::shared_ptr<Wavefunction> T__ki = T__.GetOperator(k,i);
            boost::shared_ptr<Wavefunction> T_ik  = T_.GetOpPal(i,k);
            boost::shared_ptr<Wavefunction> T__ki = T__.GetOpPal(k,i);

            //loop over locally stored operators
            //-------------------------------------------------------------------
            // case 1: both indices of (j,j_) are on the left side of the lattice
            //-------------------------------------------------------------------
            for (jj_=0;jj_<leftBlock->get_op_array(CRE_DES).get_size();jj_++){
#ifndef SERIAL
              if (leftBlock->get_sites().size()==1){
                if (world.rank()!=0){
                  continue;
                }
              }
#endif  
              //get the operator
              boost::shared_ptr<SparseMatrix> op = leftBlock->get_op_array(CRE_DES).get_local_element(jj_)[0]->getworkingrepresentation(leftBlock);
              //get the orbital indices of i and j
              j  = op->get_orbs(0);
              j_ = op->get_orbs(1);
              //generate the wavefunction that holds the result of the multiplication of 
              //the terminal operator and the linking operator
              SpinQuantum dQ=WF.get_deltaQuantum(0);
              SpinQuantum opQ=op->get_deltaQuantum(0);
              boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
              boost::shared_ptr<Wavefunction> tmp_T(new Wavefunction(dQ, &big, true));
              //multiply the linking operator with the terminal operator
              operatorfunctions::TensorMultiply(leftBlock, *op, &big, *terminal, *tmp, opQ, sqrt(2.0));
              operatorfunctions::TensorMultiply(leftBlock, Transposeview(*op), &big, *terminal, *tmp_T, opQ, sqrt(2.0));

              val1 = DotProduct(*tmp,*T_ik);
              val2 = DotProduct(*tmp_T,*T_ik);
              val3 = DotProduct(*tmp,*T__ki);
              val4 = DotProduct(*tmp_T,*T__ki);
              
              A(i_,j_,k_,i,j,k) += val1;
              if (j_!=j) A(i_,j,k_,i,j_,k) += val2;
              A(i_,j_,k_,i,j,k) -= val3;
              if (j!=j_) A(i_,j,k_,i,j_,k) -= val4;
              A(i_,i,k_,j_,k,j) -= val3;
              if (j!=j_) A(i_,i,k_,j,k,j_) -= val4;

              A_(i_,j_,k_,i,j,k) -= val2;
              if (j!=j_) A_(i_,j,k_,i,j_,k) -= val1;
              A_(i_,j_,k_,i,j,k) += val4;
              if (j!=j_) A_(i_,j,k_,i,j_,k) += val3;
              A_(i_,k,k_,j_,i,j) -= val1;
              if (j!=j_) A_(i_,k,k_,j,i,j_) -= val2;
              
            }//jj_
            //--------------------------------------------------------------------
            // case 2: both indices of (j,j_) are on the right side of the lattice
            //--------------------------------------------------------------------
            for (jj_=0;jj_<rightBlock->get_op_array(CRE_DES).get_size();jj_++){
#ifndef SERIAL
              if (rightBlock->get_sites().size()==1){
                if (world.rank()!=0){
                  continue;
                }
              }
#endif  
              //get the operator
              boost::shared_ptr<SparseMatrix> op = rightBlock->get_op_array(CRE_DES).get_local_element(jj_)[0]->getworkingrepresentation(rightBlock);
              //get the orbital indices of i and j
              j  = op->get_orbs(0);
              j_ = op->get_orbs(1);
              //generate the wavefunction that holds the result of the multiplication of 
              //the terminal operator and the linking operator
              SpinQuantum dQ=WF.get_deltaQuantum(0);
              SpinQuantum opQ=op->get_deltaQuantum(0);
              boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
              boost::shared_ptr<Wavefunction> tmp_T(new Wavefunction(dQ, &big, true));
              //multiply the linking operator with the terminal operator
              operatorfunctions::TensorMultiply(rightBlock, *op, &big, *terminal, *tmp, opQ, sqrt(2.0));
              operatorfunctions::TensorMultiply(rightBlock, Transposeview(*op), &big, *terminal, *tmp_T, opQ, sqrt(2.0));
          
              val1 = DotProduct(*tmp,*T_ik);
              val2 = DotProduct(*tmp_T,*T_ik);
              val3 = DotProduct(*tmp,*T__ki);
              val4 = DotProduct(*tmp_T,*T__ki);
              
              A(i_,j_,k_,i,j,k) += val1;
              if (j_!=j) A(i_,j,k_,i,j_,k) += val2;
              A(i_,j_,k_,i,j,k) -= val3;
              if (j!=j_) A(i_,j,k_,i,j_,k) -= val4;
              A(i_,i,k_,j_,k,j) -= val3;
              if (j!=j_) A(i_,i,k_,j,k,j_) -= val4;

              A_(i_,j_,k_,i,j,k) -= val2;
              if (j!=j_) A_(i_,j,k_,i,j_,k) -= val1;
              A_(i_,j_,k_,i,j,k) += val4;
              if (j!=j_) A_(i_,j,k_,i,j_,k) += val3;
              A_(i_,k,k_,j_,i,j) -= val1;
              if (j!=j_) A_(i_,k,k_,j,i,j_) -= val2;
              
            }//jj_
            //-----------------------------------------------------------------------------
            // case 3: j and j_ are on the left and right side of the lattice, respectively
            //-----------------------------------------------------------------------------
            PAL_DivideLoop(j0,j1,0,leftBlock->get_op_array(CRE).get_size());
            for (jcount=j0;jcount<j1;jcount++){
              //get the operator
              boost::shared_ptr<SparseMatrix> lop = leftBlock->get_op_array(CRE).get_local_element(jcount)[0]->getworkingrepresentation(leftBlock);
              //get the orbital index
              j = lop->get_orbs(0);
              for (j_count=0;j_count<rightBlock->get_op_array(CRE).get_size();j_count++){
                boost::shared_ptr<SparseMatrix> rop = rightBlock->get_op_array(CRE).get_local_element(j_count)[0]->getworkingrepresentation(rightBlock);
                j_ = rop->get_orbs(0);
                //generate the wavefunction that holds the result of the multiplication of 
                //the terminal operator and the linking operator
                SpinQuantum dQ=WF.get_deltaQuantum(0);
                SpinQuantum opQ(0,SpinSpace(0),IrrepSpace(0));
                boost::shared_ptr<Wavefunction> tmp(new Wavefunction(dQ, &big, true));
                boost::shared_ptr<Wavefunction> tmp_T(new Wavefunction(dQ, &big, true));
                //multiply the linking operator with the terminal operator
                operatorfunctions::TensorMultiply(leftBlock, *lop, Transposeview(*rop), &big, *terminal, *tmp, opQ, sqrt(2.0));
                operatorfunctions::TensorMultiply(leftBlock, Transposeview(*lop), *rop, &big, *terminal, *tmp_T, opQ, sqrt(2.0));

                val1 = DotProduct(*tmp,*T_ik);
                val2 = DotProduct(*tmp_T,*T_ik);
                val3 = DotProduct(*tmp,*T__ki);
                val4 = DotProduct(*tmp_T,*T__ki);

                A(i_,j_,k_,i,j,k) += val1;
                if (j_!=j) A(i_,j,k_,i,j_,k) += val2;
                A(i_,j_,k_,i,j,k) -= val3;
                if (j!=j_) A(i_,j,k_,i,j_,k) -= val4;
                A(i_,i,k_,j_,k,j) -= val3;
                if (j!=j_) A(i_,i,k_,j,k,j_) -= val4;

                A_(i_,j_,k_,i,j,k) -= val2;
                if (j!=j_) A_(i_,j,k_,i,j_,k) -= val1;
                A_(i_,j_,k_,i,j,k) += val4;
                if (j!=j_) A_(i_,j,k_,i,j_,k) += val3;
                A_(i_,k,k_,j_,i,j) -= val1;
                if (j!=j_) A_(i_,k,k_,j,i,j_) -= val2;

              }//j_
            }//j
          }//k
        }//i
      }//k_
    }//i_
    //sum the contributions from different processes
    SumPDM(A);
    SumPDM(A_);
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //main = tmain.elapsedwalltime();
    //Timer trest;
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //now add the remaining three-body terms that arise from rearranging the third four-body term
    double val_;
    for (i_=0;i_<norbs;i_++){
      for (j_=0;j_<norbs;j_++){
        for (k_=0;k_<norbs;k_++){
          for (i=0;i<norbs;i++){
            for (j=0;j<norbs;j++){
              for (k=0;k<norbs;k++){
                val = 0.0;
                val_= 0.0;
                for (l=0;l<norbs;l++){
                  for (m=0;m<norbs;m++){
                    val -= IJKL(2*j,2*l,2*l,2*m)  * threepdm_C(k_,i,j_,i_,k,m);
                    val += IJKL(2*j_,2*j,2*l,2*m) * threepdm_C(k_,i,l,i_,k,m);
                    val -= IJKL(2*i,2*l,2*j,2*m)  * threepdm_C(k_,j_,l,i_,k,m);
                    A(i_,j_,k_,i,j,j_) += IJKL(2*j,2*k,2*l,2*m)  * threepdm_C(k_,i,k,i_,l,m);
                    val -= IJKL(2*i,2*j,2*l,2*m)  * threepdm_C(k_,j_,l,i_,m,k);
                    val += IJKL(2*j,2*k,2*l,2*m)  * threepdm_C(k_,j_,i,i_,l,m);
                    
                    val_ -= IJKL(2*j,2*j_,2*l,2*m) * threepdm_C(k_,i,l,i_,k,m);
                    val_ += IJKL(2*j,2*l,2*m,2*m)  * threepdm_C(k_,i,l,i_,k,j_);
                    A_(i_,j_,k_,j_,j,k) -= IJKL(2*j,2*i,2*l,2*m)  * threepdm_C(k_,l,i,i_,k,m);
                    val_ += IJKL(2*j,2*l,2*k,2*m)  * threepdm_C(k_,i,l,i_,j_,m);
                    val_ -= IJKL(2*j,2*l,2*m,2*i)  * threepdm_C(k_,m,l,i_,j_,k);
                    val_ += IJKL(2*j,2*k,2*l,2*m)  * threepdm_C(k_,l,i,i_,j_,m);
                  }//m
                }//l
                A(i_,j_,k_,i,j,k) += val;
                A_(i_,j_,k_,i,j,k) += val_;
              }//k
            }//j
          }//i
        }//k_
      }//j_
    }//i_
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    //rest = trest.elapsedwalltime();
    //tot = ttot.elapsedwalltime();
    //sprintf(msg,"\nintermediate: %i s",im);pout<<msg;
    //sprintf(msg,"\nmain:         %i s",main);pout<<msg;
    //sprintf(msg,"\nrest:         %i s",rest);pout<<msg;
    //sprintf(msg,"\nTOTAL:        %i s",tot);pout<<msg;
            
    //DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!DEBUG!!!!DEBUG!!!!DEBUG!!!!DEBUG!!!DEBUG!!
    
    //clean up
    T_.Clear();
    T__.Clear();
    T.CloseFileRead();
  }

  
  
  //============================================================================
  // the function that saves the one-particle reduced density matrix in readable
  // format
  //============================================================================
  void save_onepdm_RI_text(Matrix &onepdm, int root) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    char msg[512];
    int i, j;
    int norbs = onepdm.Nrows();
    sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_onepdm_RI.", root,".rsdm");
    FILE *f;
    f = fopen(FileName, "w");
    sprintf(msg, "Spatial one-particle reduced density matrix of job %s (root %i)\n", dmrginp.save_prefix().c_str(), root);
    fprintf(f,"%s",msg);
    for (i = 0; i < norbs; i++) {
      for (j = 0; j < norbs; j++) {
        sprintf(msg, "%4i %4i     %20.14e\n", i, j,  onepdm.element(i, j));
        fprintf(f,"%s",msg);
      }
    }
    fclose(f);
#ifndef SERIAL
    }
#endif
  }

  //============================================================================
  // the function that saves the two-particle reduced density matrix in readable
  // format
  //============================================================================
  void save_twopdm_RI_text(array_4d<double>& twopdm, int root, bool predensity) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    char msg[512];
    int i, j, k, l;
    int norbs = twopdm.dim1();
    if (!predensity){
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_twopdm_RI.", root,".rsdm");
    }
    else{
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_twopdm_RI_C.", root,".rsdm");
    }
    FILE *f;
    f = fopen(FileName, "w");
    sprintf(msg, "Spatial two-particle reduced density matrix of job %s (root %i)\n", dmrginp.save_prefix().c_str(), root);
    fprintf(f,"%s",msg);
    for (i = 0; i < norbs; i++) {
      for (j = 0; j < norbs; j++) {
        for (k = 0; k < norbs; k++) {
          for (l = 0; l < norbs; l++) {
            sprintf(msg, "%4i %4i %4i %4i    %20.14e\n", i, j, k , l, twopdm(i, j, k, l));
            fprintf(f,"%s",msg);
          }
        }
      }
    }
    fclose(f);
#ifndef SERIAL
    }
#endif
  }

  
  //============================================================================
  // the function that saves the three-particle reduced density matrix in readable
  // format
  //============================================================================
  void save_threepdm_RI_text(array_6d& threepdm, int root, bool predensity, int atype) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    char msg[512];
    int i,j,k,l,m,n;
    int norbs = threepdm.get_dim1();
    if (atype==_3PDM_NORMAL_){
      if (!predensity){
        sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm_RI.", root,".rsdm");
      }
      else{
        sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm_RI_C.", root,".rsdm");
      }
    }
    else if (atype==_3PDM_A_){
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/A.", root,".rsdm");
    }
    else if (atype==_3PDM_A_PRIME_){
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/A_prime.", root,".rsdm");
    }
    
    FILE *f;
    f = fopen(FileName, "w");
    sprintf(msg, "Spatial three-particle reduced density matrix of job %s (root %i)\n", dmrginp.save_prefix().c_str(), root);
    fprintf(f,"%s",msg);
    for (i = 0; i < norbs; i++) {
      for (j = 0; j < norbs; j++) {
        for (k = 0; k < norbs; k++) {
          for (l = 0; l < norbs; l++) {
            for (m = 0; m < norbs; m ++){
              for (n = 0; n < norbs; n ++){
                sprintf(msg, "%4i %4i %4i %4i %4i %4i    %20.14e\n", i, j, k, l, m, n, threepdm(i, j, k, l, m, n));
                fprintf(f,"%s",msg);
              }
            }
          }
        }
      }
    }
    fclose(f);
#ifndef SERIAL
    }
#endif

  }

  
  //============================================================================
  // the function that saves the four-particle reduced density matrix in readable
  // format
  //============================================================================
  void save_fourpdm_RI_text(SpinBlock& big, array_4d<double>& fourpdm, int root){
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    char msg[512];
    int i,j,k,l,m,n,o,p;
    int aim,ajn,ako,alp;
    int norbs = big.get_sites().size();
    sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_fourpdm_RI.", root,".rsdm");
    FILE *f;
    f = fopen(FileName, "w");
    sprintf(msg, "Spatial four-particle reduced density matrix of job %s (root %i)\n", dmrginp.save_prefix().c_str(), root);
    fprintf(f,"%s",msg);
    for (i = 0; i < norbs; i++) {
      for (j = 0; j < norbs; j++) {
        for (k = 0; k < norbs; k++) {
          for (l = 0; l < norbs; l++) {
            for (m = 0; m < norbs; m++){
              for (n = 0; n < norbs; n++){
                for (o = 0; o < norbs; o++){
                  for (p = 0;p < norbs; p++){
                    aim = i * norbs + m;
                    ajn = j * norbs + n;
                    ako = k * norbs + o;
                    alp = l * norbs + p;
                    sprintf(msg, "%4i %4i %4i %4i %4i %4i %4i %4i   %5.14e\n", i, j, k, l, m, n,o,p, fourpdm(aim,ajn,ako,alp));
                    fprintf(f,"%s",msg);
                  }
                }
              }
            }
          }
        }
      }
    }
    fclose(f);
#ifndef SERIAL
    }
#endif    
  }

  void save_onepdm_RI_binary(Matrix &onepdm, int root) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    int norbs = onepdm.Nrows();
    sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_onepdm_RI_bin.", root,".rsdm");
    FILE *f;
    f = fopen(FileName, "wb");
    int result = fwrite(&norbs, sizeof(int), 1, f);
    result = fwrite(onepdm.Store(), sizeof(double), norbs*norbs,f); 
    fclose(f);
#ifndef SERIAL
    }
#endif    
  }
  
  void save_twopdm_RI_binary(array_4d<double> &twopdm, int root, bool predensity) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    int norbs = twopdm.dim1();
    if (!predensity){
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_twopdm_RI_bin.", root,".rsdm");
    }
    else{
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_twopdm_RI_C_bin.", root,".rsdm");
    }
    FILE *f;
    f = fopen(FileName, "wb");
    int result = fwrite(&norbs, sizeof(int), 1, f);
    result = fwrite(&twopdm(0,0,0,0), sizeof(double), norbs*norbs*norbs*norbs, f);
    fclose(f);
#ifndef SERIAL
    }
#endif    
  }
  
  void save_threepdm_RI_binary(array_6d threepdm, int root, bool predensity, int atype) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    int norbs  = threepdm.get_dim1();
    int norbs2 = norbs * norbs;
    if (atype==_3PDM_NORMAL_){
      if (!predensity){
        sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm_RI_bin.", root,".rsdm");
      }
      else{
        sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_threepdm_RI_C_bin.", root,".rsdm");
      }
    }
    else if (atype==_3PDM_A_){
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/A_bin.", root,".rsdm");
    }
    else if (atype==_3PDM_A_PRIME_){
      sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/A_prime_bin.", root,".rsdm");
    }
    FILE *f;
    f = fopen(FileName, "wb");
    int result = fwrite(&norbs, sizeof(int), 1, f);
    result = fwrite(&threepdm(0,0,0,0,0,0), sizeof(double), norbs2*norbs2*norbs2, f);
    fclose(f);
#ifndef SERIAL
    }
#endif    
  }
  
  void save_fourpdm_RI_binary(array_4d<double> &fourpdm, int root) {
#ifndef SERIAL
    mpi::communicator world;
    if (world.rank()==0){
#endif    
    char FileName[512];
    int norbs2  = fourpdm.dim1();
    sprintf(FileName, "%s%s%i%s", dmrginp.save_prefix().c_str(), "/spatial_fourpdm_RI_bin.", root,".rsdm");
    FILE *f;
    f = fopen(FileName, "wb");
    int result = fwrite(&norbs2,  1, sizeof(int), f);
    result = fwrite(&fourpdm(0,0,0,0), sizeof(double), norbs2*norbs2*norbs2*norbs2, f);
    fclose(f);
#ifndef SERIAL
    }
#endif    
  }
  
  //============================================================================
  // The driver function for the generation of  all two-particle reduced density 
  // matrices using RI
  //============================================================================
  void generate_RI_density_matrices(std::vector<Wavefunction>& wavefunctions,
          SpinBlock& big) {
    int i;
    int nroots = wavefunctions.size();
    int norbs = big.get_sites().size();
    char OpFileName[512];
    
    //the array that holds all operators <psi|E(p,q)|M> (replacement operators)
    WavefunctionArray T_rep; 
    T_rep.ResizeOrbSpace(norbs, true);
    sprintf(OpFileName, "%s%s", dmrginp.save_prefix().c_str(), "/TerminalOperators.tmp");
    T_rep.SetFileName(OpFileName);
    T_rep.SetNumFiles(mpi_world_size());
    
    //the matrix that holds the one-particle density
    Matrix onepdm(norbs,norbs);
    
    // the arrays that hold the two-particle reduced density matrices 
    array_4d<double> twopdm(norbs, norbs, norbs, norbs);
    // and the array that holds only the leading terms of the two pdm (the predensity)
    array_4d<double> twopdm_C(norbs, norbs, norbs, norbs);

    // the array that holds the three-particle reduced density matrices
    array_6d threepdm(norbs);
    // and the array that holds only the leading terms of the three pdm (the predensity)
    array_6d threepdm_C(norbs);
    
    // the arrays that hold the nevpt2 intermediate matrices A and A_
    array_6d A(norbs);
    array_6d A_(norbs);
    
    // the array that holds the four-particle reduced density matrices
    int norbs2 = norbs * norbs;
    array_4d<double> fourpdm(norbs2,norbs2,norbs2,norbs2);
    //loop over all roots
    for (i = 0; i < nroots; i++) {
      //get the wavefunction
      Wavefunction WF = wavefunctions[i];
      //reinitialize the operator array
      T_rep.ResizeOrbSpace(norbs, true);
      
      //generate the terminal vectors <psi|E(p,q)|M>
      generate_terminal_replacement_operators(big,WF,T_rep);
#ifndef SERIAL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      //---------------------------------------------------------
      //create the two-particle reduced density matrices using RI
      //---------------------------------------------------------
      generate_twopdm_RI_a(big,T_rep, WF,twopdm,twopdm_C,onepdm);
#ifndef SERIAL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      //save the one- and two particle reduced density matrices
      save_onepdm_RI_binary(onepdm,i);
      save_twopdm_RI_binary(twopdm,i);
      save_twopdm_RI_binary(twopdm_C,i,true);
      //if requested, store the pdm's in readable format
      if (dmrginp.store_ripdm_readable()){
        save_onepdm_RI_text(onepdm,i);
        save_twopdm_RI_text(twopdm,i);
        save_twopdm_RI_text(twopdm_C,i,true);
      }
      
      //-----------------------------------------------------------
      //create the three-particle reduced density matrices using RI
      //-----------------------------------------------------------
      generate_threepdm_RI(big,T_rep,threepdm,threepdm_C,twopdm,onepdm, WF);
#ifndef SERIAL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
      //save the the three particle density matrices
      save_threepdm_RI_binary(threepdm,i);
      save_threepdm_RI_binary(threepdm_C,i,true);
      //if requested, store the pdm's in readable format
      if (dmrginp.store_ripdm_readable()){
        save_threepdm_RI_text(threepdm,i);
        save_threepdm_RI_text(threepdm_C,i,true);
      }
      
      //-------------------------------------------------------------------
      //if requested, create the four-particle predensity matrices using RI
      //-------------------------------------------------------------------
      if (dmrginp.calc_ri_4pdm()){
        generate_fourpdm_RI(big,T_rep,fourpdm, WF);
        //save the four particle predensity matrix
        if (dmrginp.store_ripdm_readable()) save_fourpdm_RI_text(big,fourpdm,i);
        save_fourpdm_RI_binary(fourpdm,i);
      }

      //---------------------------------------------
      //create the nevpt2 intermediates A and A_prime
      //---------------------------------------------
      generate_A_RI(big,T_rep,A,A_,WF,threepdm_C);
      //save the A matrices
      save_threepdm_RI_binary(A,i,false,_3PDM_A_);
      save_threepdm_RI_binary(A_,i,false,_3PDM_A_PRIME_);
      //if requested, store the A matrices in readable format
      if (dmrginp.store_ripdm_readable()){
        save_threepdm_RI_text(A,i,false,_3PDM_A_);
        save_threepdm_RI_text(A_,i,false,_3PDM_A_PRIME_);
      }
      
      //clear the operator arrays
      T_rep.Clear();
    }//roots
    //clean up
    onepdm.CleanUp();
    twopdm.Clear();
    twopdm_C.Clear();
    threepdm.clear();
    threepdm_C.clear();
    fourpdm.Clear();
    A.clear();
    A_.clear();
  };











}
