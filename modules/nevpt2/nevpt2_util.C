

#include "nevpt2_operators.h"
#include "nevpt2_mpi.h"
#include "MatrixBLAS.h"
#include <boost/serialization/serialization.hpp>
#include "nevpt2_util.h"
#include "tensor_operator.h"
#include "guess_wavefunction.h"
#include "nevpt2_info.h"
#include <sys/time.h>
#ifndef SERIAL
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>
#endif

namespace SpinAdapted{
  //===================
  //print the StateInfo
  //===================
  void PrintStateInfo(const SpinBlock &block){
    char msg[512];
    SpinBlock *leftBlock = block.get_leftBlock();
    SpinBlock *rightBlock = block.get_rightBlock();
    sprintf(msg,"\n");pout << msg;
    sprintf(msg,"\nPrinting the quanta of a SpinBlock");pout << msg;
    for (int i=0;i<block.get_stateInfo().quanta.size();i++){
      SpinQuantum q = block.get_stateInfo().quanta[i];
      sprintf(msg,"\nlN=%i   lS=%i",q.get_n(),q.get_s().getirrep());pout << msg;
    }//i
    
    for (int i=0;i<leftBlock->get_stateInfo().quanta.size();i++){
      SpinQuantum lQ =leftBlock->get_stateInfo().quanta[i];
      sprintf(msg,"\nlN=%i   lS=%i",lQ.get_n(),lQ.get_s().getirrep());pout << msg;
    }//i
    for (int j=0;j<rightBlock->get_stateInfo().quanta.size();j++){
      SpinQuantum rQ = rightBlock->get_stateInfo().quanta[j];
      sprintf(msg,"\nrN=%i   rS=%i",rQ.get_n(),rQ.get_s().getirrep());pout << msg;
    }//j
  }

  //===================
  //print the StateInfo
  //===================
  void PrintStateInfo(const StateInfo &Info){
    char msg[512];
    sprintf(msg,"\n");pout << msg;
    sprintf(msg,"\nPrinting the quanta of a SpinBlock");pout << msg;
    for (int i=0;i<Info.quanta.size();i++){
      SpinQuantum q = Info.quanta[i];
      sprintf(msg,"\nlN=%i   lS=%i  Dim=%i",q.get_n(),q.get_s().getirrep(),Info.quantaStates[i]);pout << msg;
    }//i
  }

  
  //====================================
  //get the Order of the active orbitals
  //====================================
  void GetOrder(vector<int> &reorder, int NActive){
    char s[512];
    char msg[512];
    vector<string> fields;
    bool ReOrder = false;
    //try to open the Fiedler vector file
    FILE *FiedlerFile;
    
    sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/fiedler_reorder.dat");
    FiedlerFile = fopen(msg,"r");
    if (FiedlerFile!=0){
      ReOrder = true;
      fscanf(FiedlerFile,"%s",&s[0]);
      boost::split(fields, s, is_any_of( " ," ) );
      pout << s;
      fclose(FiedlerFile);
    }
    else{
      FILE *GeneticFile;
      sprintf(msg, "%s%s", dmrginp.save_prefix().c_str(), "/genetic_reorder.dat");
      GeneticFile = fopen(msg,"r");
      if (GeneticFile !=0){
        ReOrder = true;
        fscanf(GeneticFile,"%s",&s[0]);
        boost::split(fields, s, is_any_of( " ," ) );
        pout << s;
        fclose(GeneticFile);
      }
    }
    if (ReOrder){
      if (fields.size()!=NActive){
        sprintf(msg,"\nERROR BLOCK-NEVPT2: Number of reordered orbitals does not equal number of active orbitals");
        pout << msg;
        exit(0);
      }
      else{
        for (int i=0;i<fields.size();i++){
          reorder[i] = atoi(fields[i].c_str());
        }
      }
    }//ReOrder
    else{
      for (int i=0;i<reorder.size();i++){
        reorder[i] = i+1;
      }
    }//use given order
  }
  
  //========================================
  //read the orbital spaces and the BaseName
  //========================================
  void ReadInput(char *BaseName, int *OrbWin, double& ENuc, bool& ConventionalOverlap){
    FILE *f;
    f = fopen("dmrg.nevpt2.inp","r");
    fscanf(f,"%s",BaseName);
    fscanf(f,"%i",&OrbWin[0]);
    fscanf(f,"%i",&OrbWin[1]);
    fscanf(f,"%i",&OrbWin[2]);
    fscanf(f,"%i",&OrbWin[3]);
    fscanf(f,"%i",&OrbWin[4]);
    fscanf(f,"%i",&OrbWin[5]);
    fscanf(f,"%lf",&ENuc);
    int ConvOverlap=0;
    fscanf(f,"%i",&ConvOverlap);
    ConventionalOverlap = (bool) ConvOverlap;
    fclose(f);
  }
  
  //===============================================
  //read the one-body and two-body density matrices
  //===============================================
  void ReadDensityMatrices(Matrix &D,array_4d<double> &D2, int root, const vector<int> &reorder){
    char FileName[512];
    char msg[512];
    int sites = D.Ncols();
    FILE *f;
    double val;
    int sites_from_file;
    if (mpi_rank()==0){
      //open the two-body density matrix file
      sprintf(FileName, "%s%s.%i.%i.bin", dmrginp.save_prefix().c_str(), "/spatial_binary_twopdm", root,root);
      //sprintf(FileName, "%s%s.%i.%i.bin", dmrginp.save_prefix().c_str(), "/spatial_twopdm", root,root);
      //sprintf (FileName, "%s%s%d%s",dmrginp.save_prefix().c_str(),"/spatial_twopdm_RI_bin.", root,".rsdm");
      f = fopen(FileName,"rb");
      assert(f);
      fread(&sites_from_file,sizeof(int),1,f);
      //check if the two active space sizes are equal
      if (sites!=sites_from_file){
        sprintf(msg,"\n\nERROR NEVPT2: The dimension of stored density matrices does not match the number of active orbitals!");pout << msg;
        sprintf(msg,"\nLeaving NEVPT2 section.......\n\n");pout << msg;
        exit(1);
      }

      //generate the reorder vector that starts at 0 instead of 1
      vector<int> ReOrder;
      ReOrder.resize(sites);
      for (int p=0;p<sites;p++){
        ReOrder[reorder[p]-1] = p;
      }//p
      //read the two-body density matrix
      for (int i=0;i<sites;i++){
        for (int j=0;j<sites;j++){
          for (int k=0;k<sites;k++){
            for (int l=0;l<sites;l++){
              fread(&val,sizeof(double),1,f);
              //D2(i,j,l,k) = 2.0 * val;//Note: D2(i,j,k,l) = E(i,k)E(j,l)-delta(k,j)E(i,l)
                                      //this induces a reversion of the last two indices 
                                      //and requires a factor of two
              D2(ReOrder[i],ReOrder[j],ReOrder[l],ReOrder[k]) = 2.0 * val;
            }//l
          }//k
        }//j
      }//i
      fclose(f);
      //read the one-body density matrix
      //sprintf(FileName, "%s%s.%i.%i.bin", dmrginp.save_prefix().c_str(), "/spatial_onepdm", root,root);
      sprintf(FileName, "%s%s.%i.%i.bin", dmrginp.save_prefix().c_str(), "/spatial_onepdm_bin", root,root);
      //sprintf (FileName, "%s%s%d%s",dmrginp.save_prefix().c_str(),"/spatial_onepdm_RI_bin.", root,".rsdm");
      f = fopen(FileName,"rb");
      assert(f);
      fread(&sites_from_file,sizeof(int),1,f);
      fread(&sites_from_file,sizeof(int),1,f);
      //check if the two active space sizes are equal
      if (sites!=sites_from_file){
        sprintf(msg,"\n\nERROR NEVPT2: The dimension of stored density matrices does not match the number of active orbitals!");pout << msg;
        sprintf(msg,"\nLeaving NEVPT2 section.......\n\n");pout << msg;
        exit(1);
      }
      //read the two-body density matrix
      for (int i=0;i<sites;i++){
        for (int j=0;j<sites;j++){
          fread(&val,sizeof(double),1,f);
          D.element(i,j) = val;
        }
      }
      //clean up
      ReOrder.clear();
    }
    //broadcast the density matrices
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,D2,0);
    mpi::broadcast(world,D,0);
#endif
  }
  
  //============================================================================
  // Read the prestored MO integrals
  //    IKJL => (i(1)k(1)|j(2)l(2)) stored as matrices of the form Kij(k,l)
  //    IKJA => (i(1)k(1)|j(2)a(2)) stored as matrices of the form Kij(k,a)
  //============================================================================
  void ReadK(IntegralContainer &IKJL,IntegralContainer &IKJA,int *OrbWin, const char *BaseName,
        const vector<int> &ReOrder){
    FILE *f;
    double val;
    char msg[512];
    int res;
    int i,j,k,l;
    int a,b;
    int dim1,dim2;
    //the orbital spaces
    int i0 = OrbWin[0];//core
    int i1 = OrbWin[1];//core
    int t0 = OrbWin[2];//core
    int t1 = OrbWin[3];//core
    int a0 = OrbWin[4];//core
    int a1 = OrbWin[5];//core
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    //create a reorder vector for the complete set of orbitals
    vector <int> reorder;
    reorder.resize(a1-i0+1,0.0);
    //the internal part
    for (i=0;i<NInternal;i++){
      reorder[i] = i;
    }
    //the active part
    for (i=NInternal;i<NInternal+NActive;i++){
      reorder[ReOrder[i-NInternal]+NInternal-1] = i;
    }
    //the external part
    for (i=NInternal+NActive;i<NInternal+NActive+NExternal;i++){
      reorder[i] = i;
    }
    if (mpi_rank()==0){
      //---------
      //read IJKL
      //---------
      /*IJKL.OpenFileWrite();
      sprintf(msg,"%s.IJKL.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      int dim1 = t1-i0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=i;j++){
          //read the matrix Mij (k,l) = (ij|kl)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim1));
          for (k=i0;k<=t1;k++){
            for (l=i0;l<=t1;l++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(k-i0,l-i0) = val;
            }//l
          }//k
          //store the matrix
          IJKL.SetMatrix(Mij,i,j);
          if ((i==3)&&(j==2)){
            PrintMatrix(*Mij,"J32int.tmp");
          }
        }//j<=i
      }//i
      fclose(f);
      IJKL.CloseFileWrite();
      */
      //---------
      //read IKJL
      //---------
      IKJL.OpenFileWrite();
      sprintf(msg,"%s.IKJL.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=i;j++){
          //read the matrix Mij (k,l) = (ij|kl)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim1));
          for (k=i0;k<=t1;k++){
            for (l=i0;l<=t1;l++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(reorder[k-i0],reorder[l-i0]) = val;
            }//l
          }//k
          //store the matrix
          IKJL.SetMatrix(Mij,reorder[i-i0],reorder[j-i0]);
        }//j
      }//i
      fclose(f);
      IKJL.CloseFileWrite();

      //---------
      //read IKJA
      //---------
      IKJA.OpenFileWrite();
      sprintf(msg,"%s.IKJA.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      dim2 = a1-a0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=t1;j++){
          //read the matrix Mij (k,a) = (ik|ja)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim2));
          for (k=i0;k<=t1;k++){
            for (a=a0;a<=a1;a++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(reorder[k-i0],a-a0) = val;
            }//l
          }//k
          //store the matrix
          IKJA.SetMatrix(Mij,reorder[i-i0],reorder[j-i0]);
        }//j<=i
      }//i
      fclose(f);
      IKJA.CloseFileWrite();
    }//rank == 0
    //broadcast the integrals
    IKJL.Broadcast(0);
    IKJA.Broadcast(0);
    
  }
  
  
  //============================================================================
  // Read the prestored MO integrals
  //    IJKL => (i(1)j(1)|k(2)l(2)) stored as matrices of the form Jij(k,l)
  //    IKJL => (i(1)k(1)|j(2)l(2)) stored as matrices of the form Kij(k,l)
  //    IAJB => (i(1)a(1)|j(2)b(2)) stored as matrices of the form Kij(a,b)
  //    IJKA => (i(1)j(1)|k(2)a(2)) stored as matrices of the form Jij(k,a)
  //    IKJA => (i(1)k(1)|j(2)a(2)) stored as matrices of the form Kij(k,a)
  //============================================================================
  void ReadIntegrals(IntegralContainer &IJKL, IntegralContainer &IKJL, 
                     IntegralContainer &IAJB, IntegralContainer &IJAB, 
                     IntegralContainer &IJKA, IntegralContainer &IKJA, 
                     Matrix &h, int *OrbWin, const char *BaseName,const vector<int> &ReOrder){
    FILE *f;
    double val;
    char msg[512];
    int res;
    int i,j,k,l;
    int a,b;
    int dim1,dim2;
    //the orbital spaces
    int i0 = OrbWin[0];//core
    int i1 = OrbWin[1];//core
    int t0 = OrbWin[2];//core
    int t1 = OrbWin[3];//core
    int a0 = OrbWin[4];//core
    int a1 = OrbWin[5];//core
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    //create a reorder vector for the complete set of orbitals
    vector <int> reorder;
    reorder.resize(a1-i0+1,0.0);
    //the internal part
    for (i=0;i<NInternal;i++){
      reorder[i] = i;
    }
    //the active part
    for (i=NInternal;i<NInternal+NActive;i++){
      reorder[ReOrder[i-NInternal]+NInternal-1] = i;
    }
    //the external part
    for (i=NInternal+NActive;i<NInternal+NActive+NExternal;i++){
      reorder[i] = i;
    }
    sprintf(msg,"ReOrder.tmp");
    PrintVector(reorder,msg);
    //---------
    //read IJKL
    //---------
    /*IJKL.OpenFileWrite();
    sprintf(msg,"%s.IJKL.tmp",BaseName);
    f = fopen(msg,"rb");
    assert(f);
    //res = fread(&i0, 1, sizeof(int),f);
    //res = fread(&i1, 1, sizeof(int),f);
    int dim1 = t1-i0+1;
    for (i=i0;i<=t1;i++){
      for (j=i0;j<=i;j++){
        //read the matrix Mij (k,l) = (ij|kl)
        boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim1));
        for (k=i0;k<=t1;k++){
          for (l=i0;l<=t1;l++){
            res = fread(&val,sizeof(double),1,f);
            (*Mij).element(k-i0,l-i0) = val;
          }//l
        }//k
        //store the matrix
        IJKL.SetMatrix(Mij,i,j);
        if ((i==3)&&(j==2)){
          PrintMatrix(*Mij,"J32int.tmp");
        }
      }//j<=i
    }//i
    fclose(f);
    IJKL.CloseFileWrite();
    */
    if (mpi_rank()==0){
      //---------
      //read IKJL
      //---------
      IKJL.OpenFileWrite();
      sprintf(msg,"%s.IKJL.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=i;j++){
          //read the matrix Mij (k,l) = (ij|kl)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim1));
          for (k=i0;k<=t1;k++){
            for (l=i0;l<=t1;l++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(reorder[k-i0],reorder[l-i0]) = val;
            }//l
          }//k
          //store the matrix
          IKJL.SetMatrix(Mij,reorder[i-i0],reorder[j-i0]);
        }//j
      }//i
      fclose(f);
      IKJL.CloseFileWrite();
      //remove(msg);

      //---------------------------------
      //bring IKJL integrals in IJKL form
      //---------------------------------
      IJKL.OpenFileWrite();
      IKJL.OpenFileRead();
      dim1 = t1-i0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=t1;j++){
          boost::shared_ptr<Matrix> Jij(new Matrix(dim1,dim1));
          for (k=i0;k<=t1;k++){
            boost::shared_ptr<Matrix> Kik = IKJL.GetMatrix(i-i0,k-i0);
            for (l=i0;l<=t1;l++){
              Jij->element(k-i0,l-i0) = Kik->element(j-i0,l-i0);
            }//l
          }//k
          IJKL.SetMatrix(Jij,i,j);
        }//j
      }//i
      IJKL.CloseFileWrite();
      IKJL.CloseFileRead();

      //---------
      //read IAJB
      //---------
      IAJB.OpenFileWrite();
      sprintf(msg,"%s.IAJB.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      dim2 = a1-a0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=i;j++){
          //read the matrix Mij (a,b) = (ia|jb)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim2,dim2));
          for (a=a0;a<=a1;a++){
            for (b=a0;b<=a1;b++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(a-a0,b-a0) = val;
            }//l
          }//k
          //store the matrix
          IAJB.SetMatrix(Mij,reorder[i],reorder[j]);
        }//j<=i
      }//i
      fclose(f);
      IAJB.CloseFileWrite();
      //remove(msg);

      //---------
      //read IJAB
      //---------
      /*IJAB.OpenFileWrite();
      sprintf(msg,"%s.IJAB.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      dim2 = a1-a0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=i;j++){
          //read the matrix Mij (a,b) = (ia|jb)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim2,dim2));
          for (a=a0;a<=a1;a++){
            for (b=a0;b<=a1;b++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(a-a0,b-a0) = val;
            }//l
          }//k
          //store the matrix
          IJAB.SetMatrix(Mij,i,j);
        }//j<=i
      }//i
      fclose(f);
      IJAB.CloseFileWrite();
      */
      //---------
      //read IJKA
      //---------
      IJKA.OpenFileWrite();
      sprintf(msg,"%s.IJKA.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      dim2 = a1-a0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=i;j++){
          //read the matrix Mij (k,a) = (ij|ka)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim2));
          for (k=i0;k<=t1;k++){
            for (a=a0;a<=a1;a++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(reorder[k-i0],a-a0) = val;
            }//l
          }//k
          //store the matrix
          IJKA.SetMatrix(Mij,reorder[i-i0],reorder[j-i0]);
          IJKA.SetMatrix(Mij,reorder[j-i0],reorder[i-i0]);
        }//j<=i
      }//i
      fclose(f);
      IJKA.CloseFileWrite();
      //remove(msg);

      //---------
      //read IKJA
      //---------
      IKJA.OpenFileWrite();
      sprintf(msg,"%s.IKJA.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      //res = fread(&i0, 1, sizeof(int),f);
      //res = fread(&i1, 1, sizeof(int),f);
      dim1 = t1-i0+1;
      dim2 = a1-a0+1;
      for (i=i0;i<=t1;i++){
        for (j=i0;j<=t1;j++){
          //read the matrix Mij (k,a) = (ik|ja)
          boost::shared_ptr<Matrix> Mij(new Matrix(dim1,dim2));
          for (k=i0;k<=t1;k++){
            for (a=a0;a<=a1;a++){
              res = fread(&val,sizeof(double),1,f);
              (*Mij).element(reorder[k-i0],a-a0) = val;
            }//l
          }//k
          //store the matrix
          IKJA.SetMatrix(Mij,reorder[i-i0],reorder[j-i0]);
        }//j<=i
      }//i
      fclose(f);
      IKJA.CloseFileWrite();
      //remove(msg);

      //------------------------
      //read one-electron matrix
      //------------------------
      sprintf(msg,"%s.H",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      for (i=i0;i<=a1;i++){
        for (j=i0;j<=a1;j++){
          res = fread(&val,sizeof(double),1,f);
          h.element(reorder[i-i0],reorder[j-i0]) = val;
        }
      }
      fclose(f);
      //remove(msg);
    }//rank == 0
    //broadcast all integrals
    mpi_barrier();
    IKJL.Broadcast(0);
    IJKL.Broadcast(0);
    IKJA.Broadcast(0);
    IJKA.Broadcast(0);
    IAJB.Broadcast(0);
    //IJAB.Broadcast(0);
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,h,0);
#endif
    mpi_barrier();
  }
  
  //============================================================================
  // Read the effective one-electron matrices from disk
  //    heff(p,q)  = h(p,q) + 2*(IJ|pq) - (Ip|Jq)
  //    heff_(p,q) = h(p,q) + 2*(IJ|pq) - (Ip|Jq) - (Tp|Tq)
  //============================================================================
  void ReadHeff(int *OrbWin, Matrix &Heff, Matrix &Heff_, const char *BaseName,const vector<int> &ReOrder){
    char msg[512];
    FILE *f;
    double val = 0.0;
    int res;
    int i,j;
    
    //the orbital spaces
    int i0 = OrbWin[0];//core
    int i1 = OrbWin[1];//core
    int t0 = OrbWin[2];//active
    int t1 = OrbWin[3];//active
    int a0 = OrbWin[4];//external
    int a1 = OrbWin[5];//external
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    
    //create a reorder vector for the complete set of orbitals
    vector <int> reorder;
    reorder.resize(a1-i0+1,0.0);
    //the internal part
    for (i=0;i<NInternal;i++){
      reorder[i] = i;
    }
    //the active part
    for (i=NInternal;i<NInternal+NActive;i++){
      reorder[ReOrder[i-NInternal]+NInternal-1] = i;
    }
    //the external part
    for (i=NInternal+NActive;i<NInternal+NActive+NExternal;i++){
      reorder[i] = i;
    }
    
    if (mpi_rank()==0){   
      sprintf(msg,"%s.Heff.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      for (i=i0;i<=a1;i++){
        for (j=i0;j<=a1;j++){
          res = fread(&val,sizeof(double),1,f);
          Heff.element(reorder[i],reorder[j]) = val;
        }
      }
      fclose(f);
      //remove(msg);

      sprintf(msg,"%s.Heff_.tmp",BaseName);
      f = fopen(msg,"rb");
      assert(f);
      for (i=i0;i<=a1;i++){
        for (j=i0;j<=a1;j++){
          res = fread(&val,sizeof(double),1,f);
          Heff_.element(reorder[i],reorder[j]) = val;
        }
      }
      fclose(f);
      //remove(msg);
    }//rank == 0
    //broadcast the matrices
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,Heff,0);
    mpi::broadcast(world,Heff_,0);
#endif
    sprintf(msg,"Heff_.tmp");
    PrintMatrix(Heff_,msg);

  }
  
  //============================================================================
  // Generate the two effective one-electron matrices
  //    heff(p,q)  = h(p,q) + 2*(IJ|pq) - (Ip|Jq)
  //    heff_(p,q) = h(p,q) + 2*(IJ|pq) - (Ip|Jq) - (Tp|Tq)
  //============================================================================
  void GenerateHeff(int *OrbWin,Matrix &h, Matrix &h_eff, Matrix &h_eff_, 
                    IntegralContainer &IJKL, IntegralContainer &IKJL,
                    IntegralContainer &IJAB, IntegralContainer &IAJB,
                    IntegralContainer &IJKA, IntegralContainer &IKJA){
    int p,q;
    int i,t;
    int p_,q_;
    int t_;
    //the orbital spaces
    int i0 = OrbWin[0];//core
    int i1 = OrbWin[1];//core
    int t0 = OrbWin[2];//active
    int t1 = OrbWin[3];//active
    int a0 = OrbWin[4];//external
    int a1 = OrbWin[5];//external
    
    boost::shared_ptr<Matrix> Jii;
    boost::shared_ptr<Matrix> Kii;
    boost::shared_ptr<Matrix> Ktt;
    boost::shared_ptr<Matrix> Jtu;
    boost::shared_ptr<Matrix> Ktu;
    
    if (mpi_rank()==0){
    
      //open the container
      IJKL.OpenFileRead();
      IKJL.OpenFileRead();
      IAJB.OpenFileRead();
      IJKA.OpenFileRead();
      IKJA.OpenFileRead();

      // copy the one-electron part
      h_eff  = h;
      h_eff_ = h;

      //do the common two-electron part
      for (i=i0;i<=i1;i++){
        //internal/internal
        Jii = IJKL.GetMatrix(i-i0,i-i0);
        Kii = IKJL.GetMatrix(i-i0,i-i0);
        for (p=i0;p<=t1;p++){
          for (q=i0;q<=t1;q++){
            p_ = p-i0;
            q_ = q-i0;
            h_eff.element(p,q)  += 2*(*Jii).element(p_,q_);
            h_eff_.element(p,q) += 2*(*Jii).element(p_,q_);
            h_eff.element(p,q)  -= (*Kii).element(p_,q_);
            h_eff_.element(p,q) -= (*Kii).element(p_,q_);
          }//q
        }//p

        //external/external: unfortunately we do not have the integrals for this part
        /*
        Jii = IJAB.GetMatrix(i,i);
        Kii = IAJB.GetMatrix(i,i);
        for (p=a0;p<=a1;p++){
          for (q=a0;q<=a1;q++){
            p_ = p-a0;
            q_ = q-a0;
            h_eff.element(p,q)  += 2*(*Jii).element(p_,q_);
            h_eff_.element(p,q) += 2*(*Jii).element(p_,q_);
            h_eff.element(p,q)  -= (*Kii).element(p_,q_);
            h_eff_.element(p,q) -= (*Kii).element(p_,q_);
          }//q
        }//p
        */
        //external/internal
        Jii = IJKA.GetMatrix(i-i0,i-i0);
        Kii = IKJA.GetMatrix(i-i0,i-i0);
        for (p=i0;p<=t1;p++){
          for (q=a0;q<=a1;q++){
            p_ = p-i0;
            q_ = q-a0;
            h_eff.element(p,q)  += 2*(*Jii).element(p_,q_);
            h_eff_.element(p,q) += 2*(*Jii).element(p_,q_);
            h_eff.element(p,q)  -= (*Kii).element(p_,q_);
            h_eff_.element(p,q) -= (*Kii).element(p_,q_);

            h_eff.element(q,p)  += 2*(*Jii).element(p_,q_);
            h_eff_.element(q,p) += 2*(*Jii).element(p_,q_);
            h_eff.element(q,p)  -= (*Kii).element(p_,q_);
            h_eff_.element(q,p) -= (*Kii).element(p_,q_);

          }//q
        }//p
      }//i

      //add the extra term to h_eff_
      for (t=t0;t<=t1;t++){
        //internal/internal
        Ktt = IKJL.GetMatrix(t-i0,t-i0);
        for (p=t0;p<=t1;p++){
          for (q=t0;q<=t1;q++){
            p_ = p-i0;
            q_ = q-i0;
            h_eff_.element(p,q) -= 0.5 * (*Ktt).element(p_,q_);
          }//q
        }//p

        //external/external
        Ktt = IAJB.GetMatrix(t-i0,t-i0);
        for (p=a0;p<=a1;p++){
          for (q=a0;q<=a1;q++){
            p_ = p-a0;
            q_ = q-a0;
            h_eff_.element(p,q) -= 0.5 * (*Ktt).element(p_,q_);
          }//q
        }//p

        //external/internal
        Ktt = IKJA.GetMatrix(t-i0,t-i0);
        for (p=i0;p<=t1;p++){
          for (q=a0;q<=a1;q++){
            p_ = p-i0;
            q_ = q-a0;
            //Note: here we add the full integral in order to create h bar and not h prime
            h_eff_.element(p,q) -= 1.0 * (*Ktt).element(p_,q_);
            h_eff_.element(q,p) -= 1.0 * (*Ktt).element(p_,q_);
          }//q
        }//p
      }//t

      //modify heff for debugging
      for (i=0;i<=t1;i++){
        for (int j=0;j<=t1;j++){
          //if ((i==0)&&(j==3)) h_eff.element(i,j) = 1.0;
          //else if ((i==3)&&(j==0)) h_eff.element(i,j) = 1.0;
          //else h_eff.element(i,j) = 0.0;
        }
      }

      //close the container
      IJKL.CloseFileRead();
      IKJL.CloseFileRead();
      IAJB.CloseFileRead();
      IJKA.CloseFileRead();
      IKJA.CloseFileRead();
    }//rank==0

    //broadcast the results
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,h_eff,0);
    mpi::broadcast(world,h_eff_,0);
#endif

  }
  
  //============================================================================
  // Build hole densities according to 
  //
  //    D1_(p,q)     = 2*delta(a,b) - D1(q,p)
  //    D2_(p,q,r,s) = D2(r,s,p,q) + delta(p,s)*D1(r,q) + delta(q,r)*D1(s,p)
  //                   -2*delta(p,r)*D1(s,q) -2*delta(q,s)D1(r,p)
  //============================================================================
  void BuildHoleDensities(int *OrbWin,const Matrix &D1, Matrix &D1_, 
                          const array_4d<double> &D2, array_4d<double> &D2_){
    int p,q,r,s;
    //--------------------------
    //build the one-hole density
    //--------------------------
    if (mpi_rank()==0){
      Initialize(D1_);
      //diagonal part
      for (p=0;p<D1_.Nrows();p++){
        D1_.element(p,p) += 2.0;
      }//p
      //general part
      for (p=0;p<D1_.Nrows();p++){
        for (q=0;q<D1_.Ncols();q++){
          D1_.element(p,q) -= D1.element(q,p);
        }//q
      }//p

      //--------------------------
      //build the two-hole density
      //--------------------------
      Initialize(D2_);
      int NActive = D2_.dim1();
      //diagonal parts
      for (p=0;p<NActive;p++){
        for (q=0;q<NActive;q++){
          D2_(p,q,p,q) += 4.0;
          D2_(p,q,q,p) -= 2.0;
        }//p
      }//q
      for (p=0;p<NActive;p++){
        for (q=0;q<NActive;q++){
          for (r=0;r<NActive;r++){
            for (s=0;s<NActive;s++){
              D2_(p,q,r,s) += D2(r,s,p,q);
              if (p==s) D2_(p,q,r,s) += D1.element(r,q);
              if (q==r) D2_(p,q,r,s) += D1.element(s,p);
              if (p==r) D2_(p,q,r,s) -= 2.0 * D1.element(s,q);
              if (s==q) D2_(p,q,r,s) -= 2.0 * D1.element(r,p);
            }//s
          }//r
        }//q
      }//p
    }//rank==0;
    //broadcast the resulting matrices
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,D1_,0);
    mpi::broadcast(world,D2_,0);
#endif
  }
  
  
  //============================================================================
  // Construct the leading term of the 3- and 2-PDM and the three hole density 
  //            DC3(i,j,k,l,m,n) = <psi|E(i,l)E(j,m)E(k,n)|psi>
  //            DC2(i,j,k,l)     = <psi|E(i,k)E(j,l)|psi>
  //============================================================================
  void ConstructAuxPDM(array_6d &DC3, array_6d &D3_, const array_6d &D3, 
                       const array_4d<double> &D2, array_4d<double> &DC2, 
                       const array_4d<double> &D2_, const Matrix &D1, const Matrix &D1_,
                       bool Conventional){
    char msg[512];
    double val=0.0;
    int i,j,k,l,m,n;
    int dim = D1.Ncols();
    if (mpi_rank()==0){
      if (Conventional){
        //--------------------------------------
        //construct the leading term of the 3pdm
        //--------------------------------------
        for (i=0;i<dim;i++){
          for (j=0;j<dim;j++){
            for (k=0;k<dim;k++){
              for (l=0;l<dim;l++){
                for (m=0;m<dim;m++){
                  for (n=0;n<dim;n++){
                    DC3(i,j,k,l,m,n) += D3(i,j,k,l,m,n);
                  }//n
                  DC3(i,j,k,l,k,m) += D2(i,j,l,m);
                  DC3(i,j,k,k,m,l) += D2(i,j,l,m);
                  DC3(i,j,k,j,m,l) += D2(i,k,m,l);
                }//m
                DC3(i,j,k,j,k,l) += D1.element(i,l);
              }//l
            }//k
          }//j
        }//i
      }//Conventional?
      //--------------------------------------
      //construct the leading term of the 2pdm
      //--------------------------------------
      Initialize(DC2);
      for (i=0;i<dim;i++){
        for (j=0;j<dim;j++){
          for (k=0;k<dim;k++){
            for (l=0;l<dim;l++){
              DC2(i,j,k,l) += D2(i,j,k,l);
            }//l
            DC2(i,j,j,k) += D1.element(i,k);
          }//k
        }//j
      }//i
      if (Conventional){
        //----------------------------
        //construct the 3-hole density
        //----------------------------
        for (i=0;i<dim;i++){
          for (j=0;j<dim;j++){
            for (k=0;k<dim;k++){
              for (l=0;l<dim;l++){
                for (m=0;m<dim;m++){
                  for (n=0;n<dim;n++){
                    D3_(i,j,k,l,m,n) -= DC3(l,m,n,i,j,k);
                  }//n
                  D3_(i,j,k,l,k,m) -= D2_(i,j,l,m);
                  D3_(i,j,k,k,m,l) -= D2_(i,j,l,m);
                  D3_(i,j,k,j,m,l) -= D2_(i,k,m,l);

                  D3_(i,j,k,i,m,l) += 2.0 * DC2(m,l,j,k);
                  D3_(i,j,k,l,j,m) += 2.0 * DC2(l,m,i,k);
                  D3_(i,j,k,l,m,k) += 2.0 * DC2(l,m,i,j);
                }//m
                D3_(i,j,k,j,k,l) -= D1_.element(i,l);

                D3_(i,j,k,i,j,l) -= 4.0 * D1.element(l,k);
                D3_(i,j,k,l,j,k) -= 4.0 * D1.element(l,i);
                D3_(i,j,k,i,l,k) -= 4.0 * D1.element(l,j);
              }//l
              D3_(i,j,k,i,j,k) += 8.0;
            }//k
          }//j
        }//i
      }//Conventional
    }//rank==0
    //broadcast the resulting matrices
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,DC2,0);
    if (Conventional) mpi::broadcast(world,DC3,0);
    if (Conventional) mpi::broadcast(world,D3_,0);
#endif
  }
  
  //============================================================================
  // Calculate the energy that arises from taking into account the core orbitals
  // E(core) = h(I,I) + 2*(II|JJ) - (IJ|IJ) + 2(II|TU) D(T,U) - (IT|JU) * D(T,U)
  //============================================================================
  void CalcCoreEnergy(double &CoreEnergy, Matrix &h, Matrix &D, int i0,
                      int i1, int t0, int t1, IntegralContainer &IJKL,
                      IntegralContainer &IKJL, array_4d<double> &D2,double &Eval){
    int i,j,t,u;
    int t_,u_,j_;
    int t__,u__;
    int v,w,v_,w_,v__,w__;
    if (mpi_rank()==0){
      boost::shared_ptr<Matrix> Jii;
      boost::shared_ptr<Matrix> Kii;
      boost::shared_ptr<Matrix> Jtu;
      IJKL.OpenFileRead();
      IKJL.OpenFileRead();
      double EActive = 0.0;
      double SingleEl = 0.0;
      double TwoEl = 0.0;
      CoreEnergy=0.0;
      Eval=0.0;
      char msg[512];
      //------------------
      // one-electron part
      //------------------
      //internal
      for (i=i0;i<=i1;i++){
        CoreEnergy += h.element(i,i) * 2.0;
      }//i
      //active
      for (t=t0;t<=t1;t++){
        for (u=t0;u<=t1;u++){
          t__ = t-t0;
          u__ = u-t0;
          CoreEnergy += h.element(t,u) * D.element(t__,u__);
          Eval       += h.element(t,u) * D.element(t__,u__);
        }//u
      }//t
      //------------------
      // two-electron part
      //------------------
      //internal/internal
      for (i=i0;i<=i1;i++){
        Jii = IJKL.GetMatrix(i-i0,i-i0);
        Kii = IKJL.GetMatrix(i-i0,i-i0);
        for (j=i0;j<=i1;j++){
          j_ = j-i0;
          CoreEnergy += 2.0 * (*Jii).element(j_,j_);//Coulomb
          CoreEnergy -=       (*Kii).element(j_,j_);//Exchange
        }//j<=i
      }//i
      //internal/active
      for (i=i0;i<=i1;i++){
        Jii = IJKL.GetMatrix(i-i0,i-i0);
        Kii = IKJL.GetMatrix(i-i0,i-i0);
        for (t=t0;t<=t1;t++){
          for (u=t0;u<=t1;u++){
            t_ = t-i0;
            u_ = u-i0;
            t__ = t-t0;
            u__ = u-t0;
            CoreEnergy += (*Jii).element(t_,u_) * 2.0 * D.element(t__,u__);//Coulomb
            CoreEnergy -= (*Kii).element(t_,u_) *       D.element(t__,u__);//Exchange
            Eval       += (*Jii).element(t_,u_) * 2.0 * D.element(t__,u__);//Coulomb
            Eval       -= (*Kii).element(t_,u_) *       D.element(t__,u__);//Exchange
          }//u
        }//t
      }//i
      //active/active
      for (t=t0;t<=t1;t++){
        for (u=t0;u<=t1;u++){
          Jtu = IJKL.GetMatrix(t-i0,u-i0);
          for (v=t0;v<=t1;v++){
            for (w=t0;w<=t1;w++){
              t_ = t-i0;
              u_ = u-i0;
              v_ = v-i0;
              w_ = w-i0;
              t__ = t-t0;
              u__ = u-t0;
              v__ = v-t0;
              w__ = w-t0;
              CoreEnergy += 0.5*(*Jtu).element(v_,w_) * D2(t__,v__,u__,w__);
              Eval       += 0.5*(*Jtu).element(v_,w_) * D2(t__,v__,u__,w__);
            }
          }
        }
      }
      IJKL.CloseFileRead();
      IKJL.CloseFileRead();
    }//rank==0;
    //broadcast the results
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,Eval,0);
    mpi::broadcast(world,CoreEnergy,0);
#endif
  }
  
  //============================================================================
  // Read the orbital energies from disk
  //============================================================================
  void ReadOrbEnergies(vector<double> &EOrb, const char* BaseName, int *OrbWin, const vector<int> &ReOrder){
    FILE *f;
    char msg[512];
    //the orbital spaces
    int i0 = OrbWin[0];//core
    int i1 = OrbWin[1];//core
    int t0 = OrbWin[2];//core
    int t1 = OrbWin[3];//core
    int a0 = OrbWin[4];//core
    int a1 = OrbWin[5];//core
    int NInternal = i1-i0+1;
    int NActive = t1-t0+1;
    int NExternal = a1-a0+1;
    int i;
    if (mpi_rank()==0){
      //create a reorder vector for the complete set of orbitals
      vector <int> reorder;
      reorder.resize(a1-i0+1,0.0);
      //the internal part
      for (i=0;i<NInternal;i++){
        reorder[i] = i;
      }
      //the active part
      for (i=NInternal;i<NInternal+NActive;i++){
        reorder[ReOrder[i-NInternal]+NInternal-1] = i;
      }
      //the external part
      for (i=NInternal+NActive;i<NInternal+NActive+NExternal;i++){
        reorder[i] = i;
      }
      sprintf(msg,"%s.EOrb.tmp",BaseName);
      f = fopen(msg,"r");
      double val = 0.0;
      for (i=0;i<EOrb.size();i++){
        fscanf(f,"%lf",&val);
        EOrb[reorder[i]] = val;
      }
      fclose(f);
    }//rank == 0
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,EOrb,0);
#endif

  }
  
  //============================================================================
  // Calculate the orbital energies from MO integrals according to 
  //    e(p) = h(p,p) + 2(pp|II) - (pI|pI) + D(T,T) {(pp|TT) - 0.5(pT|pT)}
  //============================================================================
  void  CalcOrbEnergies(int*OrbWin, Matrix &h, IntegralContainer &IJKL, IntegralContainer &IKJL, 
                        IntegralContainer &IJAB, IntegralContainer &IAJB, vector<double> &EOrb,
                        Matrix &D){
    
    int p,q;
    int i,j,t,u,a,b;
    char msg[512];
    
    //the orbital spaces
    int i0 = OrbWin[0];
    int i1 = OrbWin[1];
    int t0 = OrbWin[2];
    int t1 = OrbWin[3];
    int a0 = OrbWin[4];
    int a1 = OrbWin[5];
    
    int NInternal = i1-i0+1;
    int NActive   = t1-t0+1;
    int NExternal = a1-a0+1;
    int OrbDim = EOrb.size();
    if (mpi_rank()==0){
      //Open the integral container
      IJKL.OpenFileRead();
      IKJL.OpenFileRead();
      //-----------------
      //internal orbitals
      //-----------------
      for (i=i0;i<=i1;i++){
        //the one electron part
        EOrb[i] = h.element(i,i);

        //get the J and K integrals
        boost::shared_ptr<Matrix> Jii = IJKL.GetMatrix(i-i0,i-i0);
        boost::shared_ptr<Matrix> Kii = IKJL.GetMatrix(i-i0,i-i0);
        //internal two-electron part
        for (j=i0;j<=i1;j++){
          EOrb[i] += 2.0 * Jii->element(j-i0,j-i0);
          EOrb[i] -= 1.0 * Kii->element(j-i0,j-i0);
        }//j

        //the active two-electron part
        for (t=t0;t<=t1;t++){
          EOrb[i] += 1.0 * D.element(t-t0,t-t0) * Jii->element(t-i0,t-i0);
          EOrb[i] -= 0.5 * D.element(t-t0,t-t0) * Kii->element(t-i0,t-i0);
        }//t
      }//i

      //---------------
      //active orbitals
      //---------------
      for (t=t0;t<=t1;t++){
        //the one-electron part
        EOrb[t] = h.element(t,t);
        //get the J and K integrals
        boost::shared_ptr<Matrix> Jtt = IJKL.GetMatrix(t-i0,t-i0);
        boost::shared_ptr<Matrix> Ktt = IKJL.GetMatrix(t-i0,t-i0);
        //internal two-electron part
        for (j=i0;j<=i1;j++){
          EOrb[t] += 2.0 * Jtt->element(j-i0,j-i0);
          EOrb[t] -= 1.0 * Ktt->element(j-i0,j-i0);
        }//j
        //active two-electron part
        for (u=t0;u<=t1;u++){
          EOrb[t] += 1.0 * D.element(u-t0,u-t0) * Jtt->element(u-i0,u-i0);
          EOrb[t] -= 0.5 * D.element(u-t0,u-t0) * Ktt->element(u-i0,u-i0);        
        }//uu
      }//t

      //close internal integral container and open the external ones
      IJKL.CloseFileRead();
      IKJL.CloseFileRead();
      IJAB.OpenFileRead();
      IAJB.OpenFileRead();

      //------------------
      // external orbitals
      //------------------
      //the one-electron part
      for (a=a0;a<=a1;a++){
        EOrb[a] = h.element(a,a);
      }//a
    }//rank == 0;
    //broadcast the results
#ifndef SERIAL
    mpi::communicator world;
    mpi::broadcast(world,EOrb,0);
#endif

    /*
    //the internal two-electron part
    for (i=i0;i<=i1;i++){
      //get the J and K integrals
      boost::shared_ptr<Matrix> Jii = IJAB.GetMatrix(i-i0,i-i0);
      boost::shared_ptr<Matrix> Kii = IAJB.GetMatrix(i-i0,i-i0);
      for (a=a0;a<=a1;a++){
        EOrb[a] += 2.0 * Jii->element(a-a0,a-a0);
        EOrb[a] -= 1.0 * Kii->element(a-a0,a-a0);
      }//a
    }//i  
    
    //the active two-electron part
    for (t=t0;t<=t1;t++){
      //get the J and K integrals
      boost::shared_ptr<Matrix> Jtt = IJAB.GetMatrix(t-i0,t-i0);
      boost::shared_ptr<Matrix> Ktt = IAJB.GetMatrix(t-i0,t-i0);
      for (a=a0;a<=a1;a++){
        EOrb[a] += 1.0 * D.element(t-t0,t-t0) * Jtt->element(a-a0,a-a0);
        EOrb[a] -= 0.5 * D.element(t-t0,t-t0) * Ktt->element(a-a0,a-a0);
      }//a
    }//i  
    //close the remaining open integral container
    IJAB.CloseFileRead();
    IAJB.CloseFileRead();*/
  }
  
  
  //============================================================================
  // a function that evaluates factors for complementary operators in NEVPT2
  //============================================================================
  void EvalCompFactors(vector<double> &Fac, int Gamma, int S_, int S){
    char msg[512];
    double fac=0.0;
    int m,M_;
    double CGC=0.0;
    double prefac = 1.0;
    for (m=-Gamma;m<=Gamma;m+=2){
      //determine the prefactor
      if ((Gamma-m)%2==0) prefac = 1.0;
      else prefac = -1.0;
      for (M_=-S_;M_<=S_;M_+=2){
        //get the CGC
        CGC = clebsch(Gamma,m,S_,M_,S,S);
        fac += prefac * CGC * CGC;
        //sprintf(msg,"\nm=%i M_=%i fac=%lf",m,M_,fac);pout << msg;
      }//M_
    }//m
    double twoS_ = (double) S_;
    double twoS  = (double) S;
    double s_    = twoS_ / 2.0;
    double s     = twoS / 2.0;
    fac *= sqrt((2.0*s_+1)/(2.0*s+1));    
    Fac.push_back(fac);
  }
  
  //============================================================================
  //a little help function that facilitates the measurement of time
  //============================================================================
  double GetTime(){
    static long zsec = 0;
    static long zusec= 0;
    double ret;
    struct timeval tp;
    struct timezone tzp;

    gettimeofday(&tp, &tzp);
    if (zsec   == 0) zsec = tp.tv_sec;
    if ( zusec == 0) zusec= tp.tv_usec;

    ret= (tp.tv_sec-zsec) + (tp.tv_usec-zusec)*0.000001;
    
    return ret;
  }

  //===================================
  // Set all values in a matrix to zero
  //===================================
  void Initialize(Matrix& M){
    for (int row=0;row<M.Nrows();row++){
      for (int col=0;col<M.Ncols();col++){
        M.element(row,col) = 0.0;
      }//col
    }//row
  }
  void Initialize(array_4d<double> &M){
    for (int p=0;p<M.dim1();p++){
      for (int q=0;q<M.dim2();q++){
        for (int r=0;r<M.dim3();r++){
          for (int s=0;s<M.dim4();s++){
            M(p,q,r,s) = 0.0;
          }//s
        }//r
      }//q
    }//p
  }
  
  //========================
  // Read the 3PDM from disk
  //========================
  void Read3PDM(array_6d &D3, const char* FileName){
    char msg[512];
    double val=0.0;
    int res;
    int norbs=0;
    FILE *f;
    f = fopen(FileName,"rb");
    if (f==0){
      sprintf(msg,"\nERROR: #PDM File could not be opened for reading!!!");
      pout << msg;
      exit(EXIT_FAILURE);
    }
    res = fread(&norbs,1,sizeof(int),f);
    for (int i=0;i<norbs;i++){
      for (int j=0;j<norbs;j++){
        for (int k=0;k<norbs;k++){
          for (int l=0;l<norbs;l++){
            for (int m=0;m<norbs;m++){
              for (int n=0;n<norbs;n++){
                res = fread(&val,1,sizeof(double),f);
                D3(i,j,k,l,m,n) = val;
              }//n
            }//m
          }//l
        }//k
      }//j
    }//i
  }
  
  //===================
  // Transpose a Matrix
  //===================
  void Transpose(Matrix &M){
    double x;
    int i,j;
    for (i=0;i<M.Nrows();i++){
      for (j=0;j<=i;j++){
        x = M.element(i,j);
        M.element(i,j) = M.element(j,i);
        M.element(j,i) = x;
      }//j
    }//i
  }
  
  
  //===========================================
  // print a wavefunction or operator to a file
  //===========================================
  void PrintWavefunction(const Baseoperator<Matrix> &WF,const char *FileName){
    int iquanta,jquanta;
    int i,j;
    char msg[512];
    FILE *f;
    f = fopen(FileName,"w");
    //print the header
    sprintf(msg,"\n\n\nWavefunction\n");
    fprintf(f,"%s",msg);
    Matrix M;
    //print the Delta Quantum
    sprintf(msg,"Particle Number = %i      Spin = %i\n",WF.get_deltaQuantum(0).get_n(),WF.get_deltaQuantum(0).get_s().getirrep());
    fprintf(f,"%s",msg);
    for (iquanta=0;iquanta<WF.nrows();iquanta++){
      for (jquanta=0;jquanta<WF.ncols();jquanta++){
        if (WF.allowed(iquanta,jquanta)){
          sprintf(msg,"\n\nNRows=%i      NCols=%i",WF.operator_element(iquanta,jquanta).Nrows(),WF.operator_element(iquanta,jquanta).Ncols());
          fprintf(f,"%s",msg);
          //print the matrix
          PrintMatrix(WF.operator_element(iquanta,jquanta),f);
          /*for (i=0;i<WF.operator_element(iquanta,jquanta).Nrows();i++){
            for (j=0;j<WF.operator_element(iquanta,jquanta).Ncols();j++){
              double val = WF.operator_element(iquanta,jquanta).element(i,j);
              sprintf(msg,"\nM(%i,%i) = %5.7e",i,j,val);
              fprintf(f,msg);
            }
          }*/
        }
      }
    }
    fclose(f);
  }
  
  //===========================================
  // print a wavefunction or operator to a file
  //===========================================
  void PrintWavefunction(const Baseoperator<Matrix> &WF,const char *FileName, const SpinBlock &big){
    int iquanta,jquanta;
    int i,j;
    char msg[512];
    FILE *f;
    f = fopen(FileName,"w");
    StateInfo *lS = big.get_stateInfo().leftStateInfo;
    StateInfo *rS = big.get_stateInfo().rightStateInfo;
    //print the header
    sprintf(msg,"\n\n\nWavefunction\n");
    fprintf(f,"%s",msg);
    Matrix M;
    //print the Delta Quantum
    sprintf(msg,"Particle Number = %i      Spin = %i\n",WF.get_deltaQuantum(0).get_n(),WF.get_deltaQuantum(0).get_s().getirrep());
    fprintf(f,"%s",msg);
    for (iquanta=0;iquanta<WF.nrows();iquanta++){
      for (jquanta=0;jquanta<WF.ncols();jquanta++){
        if (WF.allowed(iquanta,jquanta)){
          sprintf(msg,"\n\n(lN=%i  lS=%i)   x  (rN=%i  rS=%i)",lS->quanta[iquanta].get_n(),lS->quanta[iquanta].get_s().getirrep(),
                                                               rS->quanta[jquanta].get_n(),rS->quanta[jquanta].get_s().getirrep());
          fprintf(f,"%s",msg);
          //print the matrix
          PrintMatrix(WF.operator_element(iquanta,jquanta),f);
          /*
          for (i=0;i<WF.operator_element(iquanta,jquanta).Nrows();i++){
            for (j=0;j<WF.operator_element(iquanta,jquanta).Ncols();j++){
              double val = WF.operator_element(iquanta,jquanta).element(i,j);
              sprintf(msg,"\nM(%i,%i) = %5.7e",i,j,val);
              fprintf(f,msg);
            }
          }*/
        }
      }
    }
    fclose(f);
  }
  
  //===========================================
  // print a wavefunction or operator to a file
  //===========================================
  void PrintWavefunction(const Baseoperator<Matrix> &WF,const char *FileName, const StateInfo &Info){
    int iquanta,jquanta;
    int i,j;
    char msg[512];
    FILE *f;
    f = fopen(FileName,"w");
    StateInfo *lS = Info.leftStateInfo;
    StateInfo *rS = Info.rightStateInfo;
    //print the header
    sprintf(msg,"\n\n\nWavefunction\n");
    fprintf(f,"%s",msg);
    Matrix M;
    //print the Delta Quantum
    sprintf(msg,"Particle Number = %i      Spin = %i\n",WF.get_deltaQuantum(0).get_n(),WF.get_deltaQuantum(0).get_s().getirrep());
    fprintf(f,"%s",msg);
    for (iquanta=0;iquanta<WF.nrows();iquanta++){
      for (jquanta=0;jquanta<WF.ncols();jquanta++){
        if (WF.allowed(iquanta,jquanta)){
          sprintf(msg,"\n\n(lN=%i  lS=%i)   x  (rN=%i  rS=%i)",lS->quanta[iquanta].get_n(),lS->quanta[iquanta].get_s().getirrep(),
                                                               rS->quanta[jquanta].get_n(),rS->quanta[jquanta].get_s().getirrep());
          fprintf(f,"%s",msg);
          //print the matrix
          PrintMatrix(WF.operator_element(iquanta,jquanta),f);
          /*
          for (i=0;i<WF.operator_element(iquanta,jquanta).Nrows();i++){
            for (j=0;j<WF.operator_element(iquanta,jquanta).Ncols();j++){
              double val = WF.operator_element(iquanta,jquanta).element(i,j);
              sprintf(msg,"\nM(%i,%i) = %5.7e",i,j,val);
              fprintf(f,msg);
            }
          }*/
        }
      }
    }
    fclose(f);
  }
  
  //===========================================
  // print a wavefunction or operator to a file
  //===========================================
  void PrintOperator(const Baseoperator<Matrix> &WF,const char *FileName, const StateInfo &Info){
    int iquanta,jquanta;
    int i,j;
    char msg[512];
    FILE *f;
    f = fopen(FileName,"w");
    //print the header
    sprintf(msg,"\n\n\nWavefunction\n");
    fprintf(f,"%s",msg);
    Matrix M;
    //print the Delta Quantum
    sprintf(msg,"Particle Number = %i      Spin = %i\n",WF.get_deltaQuantum(0).get_n(),WF.get_deltaQuantum(0).get_s().getirrep());
    fprintf(f,"%s",msg);
    for (iquanta=0;iquanta<WF.nrows();iquanta++){
      for (jquanta=0;jquanta<WF.ncols();jquanta++){
        if (WF.allowed(iquanta,jquanta)){
          sprintf(msg,"\n\n(lN=%i  lS=%i)   x  (rN=%i  rS=%i)",Info.quanta[iquanta].get_n(),Info.quanta[iquanta].get_s().getirrep(),
                                                               Info.quanta[jquanta].get_n(),Info.quanta[jquanta].get_s().getirrep());
          fprintf(f,"%s",msg);
          //print the matrix
          PrintMatrix(WF.operator_element(iquanta,jquanta),f);
          /*
          for (i=0;i<WF.operator_element(iquanta,jquanta).Nrows();i++){
            for (j=0;j<WF.operator_element(iquanta,jquanta).Ncols();j++){
              double val = WF.operator_element(iquanta,jquanta).element(i,j);
              sprintf(msg,"\nM(%i,%i) = %5.7e",i,j,val);
              fprintf(f,msg);
            }
          }*/
        }
      }
    }
    fclose(f);
  }
  
  void PrintWavefunctionProperties(const Baseoperator<Matrix> &WF,const StateInfo *lS,const StateInfo *rS){
    
    int iquanta,jquanta;
    int i,j;
    char msg[512];
    //print the header
    sprintf(msg,"\nWavefunction");
    pout << msg;
    Matrix M;
    //print the Delta Quantum
    sprintf(msg,"Particle Number = %i      Spin = %i\n",WF.get_deltaQuantum(0).get_n(),WF.get_deltaQuantum(0).get_s().getirrep());
    pout << msg;
    for (iquanta=0;iquanta<WF.nrows();iquanta++){
      for (jquanta=0;jquanta<WF.ncols();jquanta++){
        if (WF.allowed(iquanta,jquanta)){
          sprintf(msg,"\n(lN=%i  lS=%i  dim=%i)   x  (rN=%i  rS=%i  dim=%i)",lS->quanta[iquanta].get_n(),lS->quanta[iquanta].get_s().getirrep(),WF(iquanta,jquanta).Nrows(),
                                                                               rS->quanta[jquanta].get_n(),rS->quanta[jquanta].get_s().getirrep(),WF(iquanta,jquanta).Ncols());
          pout << msg;
        }
      }
    }
    
  }
  
  //=========================
  // print a vector to a file
  //=========================
  void PrintVector(const vector<double> &v,const char* name){
    int dim = v.size();
    char msg[512];
    FILE *f;
    f = fopen(name,"w");
    for (int i=0;i<dim;i++){
      sprintf(msg,"\n%s[%i] = %8.8lf",name,i,v[i]);
      fprintf(f,"%s",msg);
    }
    fclose(f);
  }
  void PrintVector(const vector<int> &v,const char* name){
    int dim = v.size();
    char msg[512];
    FILE *f;
    f = fopen(name,"w");
    for (int i=0;i<dim;i++){
      sprintf(msg,"\n%s[%i] = %i",name,i,v[i]);
      fprintf(f,"%s",msg);
    }
    fclose(f);
  }
  
  //=======================
  // Print a 2pdm to a file
  //=======================
  void Print2PDM(const array_4d<double> &D2,const char *name){
    FILE *f;
    char msg[512];
    int norbs = D2.dim1();
    f = fopen(name, "w");
    sprintf(msg, "Spatial two-particle reduced density matrix of job %s \n", dmrginp.save_prefix().c_str());
    fprintf(f,"%s",msg);
    int i,j,k,l;
    for (i = 0; i < norbs; i++) {
      for (j = 0; j < norbs; j++) {
        for (k = 0; k < norbs; k++) {
          for (l = 0; l < norbs; l++) {
            sprintf(msg, "%4i %4i %4i %4i    %20.14e\n", i, j, k , l, D2(i, j, k, l));
            fprintf(f,"%s",msg);
          }
        }
      }
    }
    fclose(f);
  }
  
  //=======================
  // Print a 3pdm to a file
  //=======================
  void Print3PDM(const array_6d &D3, const char *name){
    FILE *f;
    char msg[512];
    int norbs = D3.get_dim1();
    f = fopen(name, "w");
    sprintf(msg, "Spatial three-particle reduced density matrix of job %s \n", dmrginp.save_prefix().c_str());
    fprintf(f,"%s",msg);
    int i,j,k,l,m,n;
    for (i = 0; i < norbs; i++) {
      for (j = 0; j < norbs; j++) {
        for (k = 0; k < norbs; k++) {
          for (l = 0; l < norbs; l++) {
            for (m = 0; m < norbs; m++){
              for (n = 0; n < norbs; n++){
                sprintf(msg, "%4i %4i %4i %4i %4i %4i   %20.14e\n", i, j, k , l, m, n, D3(i, j, k, l, m, n));
                fprintf(f,"%s",msg);
              }//n
            }//m
          }//l
        }//k
      }//j
    }//i
    fclose(f);
  }
  
  //=========================
  // print a Matrix to a file
  //=========================
  void PrintMatrix(const Matrix &M, const char *Name){
    FILE *f;
    f = fopen(Name,"w");
    int i1 = M.Nrows();
    int j1 = M.Ncols();
    
    int a1 = j1 / 6;
    int a2 = j1 % 6;
    
    for (int block=0;block<a1;block++){
      //---------
      //head line
      //---------
      fprintf(f,"            ");
      for (int j=0;j<6;j++){
        //the column numbering
        fprintf(f,"    %3i    ",6*block+j);
      }//j
      fprintf(f,"\n");
      //------
      //values
      //------
      for (int i=0;i<i1;i++){
        //the row numbering
        fprintf(f,"    %3i    ",i);
        for (int j=0;j<6;j++){
          fprintf(f,"%11.6lf",M.element(i,6*block+j));
        }//j
        fprintf(f,"\n");
      }//i
    }//block
    
    //print the rest of the matrix
    if (a2>0){
      //---------
      //head line
      //---------
      fprintf(f,"            ");
      for (int j=0;j<a2;j++){
        //the column numbering
        fprintf(f,"    %3i    ",6*a1+j);
      }//j
      fprintf(f,"\n");
      //------
      //values
      //------
      for (int i=0;i<i1;i++){
        //the row numbering
        fprintf(f,"    %3i    ",i);
        for (int j=0;j<a2;j++){
          fprintf(f,"%11.6lf",M.element(i,6*a1+j));
        }//j
        fprintf(f,"\n");
      }//i
    }//have rest?
    
    fclose(f);
  }
  
  void PrintMatrix(const Matrix &M, FILE *f){
    int i1 = M.Nrows();
    int j1 = M.Ncols();
    
    int a1 = j1 / 6;
    int a2 = j1 % 6;
    
    for (int block=0;block<a1;block++){
      //---------
      //head line
      //---------
      fprintf(f,"            ");
      for (int j=0;j<6;j++){
        //the column numbering
        fprintf(f,"    %3i    ",6*block+j);
      }//j
      fprintf(f,"\n");
      //------
      //values
      //------
      for (int i=0;i<i1;i++){
        //the row numbering
        fprintf(f,"    %3i    ",i);
        for (int j=0;j<6;j++){
          fprintf(f,"%11.6lf",M.element(i,6*block+j));
        }//j
        fprintf(f,"\n");
      }//i
    }//block
    
    //print the rest of the matrix
    if (a2>0){
      //---------
      //head line
      //---------
      fprintf(f,"            ");
      for (int j=0;j<a2;j++){
        //the column numbering
        fprintf(f,"    %3i    ",6*a1+j);
      }//j
      fprintf(f,"\n");
      //------
      //values
      //------
      for (int i=0;i<i1;i++){
        //the row numbering
        fprintf(f,"    %3i    ",i);
        for (int j=0;j<a2;j++){
          fprintf(f,"%11.6lf",M.element(i,6*a1+j));
        }//j
        fprintf(f,"\n");
      }//i
    }//have rest?
  }
  
  
  //============================================================================
  //give the output for an unconventional NEVPT2 calculation
  //============================================================================
  void GiveNEVPT2Output(NEVPT2Info &Info){
    if (mpi_rank()==0){
      char msg[512];
      double TotEnergy=0.0;
      int nroots = Info.getNRoots();
      bool ConventionalOverlap = Info.ConventionalOverlap();
      int OrbWin[6];
      Info.GetOrbWin(OrbWin);
      int i0 = OrbWin[0];
      int i1 = OrbWin[1];
      int t0 = OrbWin[2];
      int t1 = OrbWin[3];
      int a0 = OrbWin[4];
      int a1 = OrbWin[5];
      //--------------
      //print a header
      //--------------
      sprintf(msg,"\n");pout << msg;
      sprintf(msg,"\n=============================================================");pout << msg;
      sprintf(msg,"\n                   BLOCK-NEVPT2 Calculation");pout << msg;
      sprintf(msg,"\n=============================================================");pout << msg;
      sprintf(msg,"\n");pout << msg;
      sprintf(msg,"\nInternal Orbitals    [%4i ; %4i]",i0,i1);pout << msg;
      sprintf(msg,"\nActive Orbitals      [%4i ; %4i]",t0,t1);pout << msg;
      sprintf(msg,"\nExternal Orbitals    [%4i ; %4i]",a0,a1);pout << msg;
      sprintf(msg,"\nOverlap is calculated in ");pout << msg;
      if (ConventionalOverlap) {sprintf(msg,"conventional ");pout << msg;}
      else {sprintf(msg,"non-conventional ");pout << msg;}
      sprintf(msg,"fashion!");pout << msg;
      //------------------
      //print the energies
      //------------------
      for (int iroot=0;iroot<nroots;iroot++){
        //read the energies
        //give the information
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"    **********\n");pout << msg;
        sprintf(msg,"      ROOT %i\n",iroot);pout << msg;
        sprintf(msg,"    **********\n");pout << msg;
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"-------------------------\n");pout << msg;
        sprintf(msg,"E0 = %6.8lf\n",Info.getE0(iroot));pout << msg;TotEnergy = Info.getE0(iroot);
        sprintf(msg,"-------------------------\n");pout << msg;
        sprintf(msg,"\nE(0,ijab) = %5.12lf",Info.getE(iroot,0));pout << msg;TotEnergy += Info.getE(iroot,0);
        sprintf(msg,"\nE(-1,iab) = %5.12lf",Info.getE(iroot,1));pout << msg;TotEnergy += Info.getE(iroot,1);
        sprintf(msg,"\nE(1,ija)  = %5.12lf",Info.getE(iroot,2));pout << msg;TotEnergy += Info.getE(iroot,2);
        sprintf(msg,"\nE(-2,ab)  = %5.12lf",Info.getE(iroot,3));pout << msg;TotEnergy += Info.getE(iroot,3);
        sprintf(msg,"\nE(2,ij)   = %5.12lf",Info.getE(iroot,4));pout << msg;TotEnergy += Info.getE(iroot,4);
        sprintf(msg,"\nE(0,ia)   = %5.12lf",Info.getE(iroot,5));pout << msg;TotEnergy += Info.getE(iroot,5);
        sprintf(msg,"\nE(-1,a)   = %5.12lf",Info.getE(iroot,6));pout << msg;TotEnergy += Info.getE(iroot,6);
        sprintf(msg,"\nE(1,i)    = %5.12lf",Info.getE(iroot,7));pout << msg;TotEnergy += Info.getE(iroot,7);
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"\nTotal Energy Correction: %6.12lf",TotEnergy-Info.getE0(iroot));pout << msg;
        sprintf(msg,"\n");pout << msg;
        sprintf(msg,"\n-------------------------------");pout << msg;
        sprintf(msg,"\nTOTAL ENERGY = %6.12lf",TotEnergy);pout << msg;
        sprintf(msg,"\n-------------------------------");pout << msg;
      }//iroot
      //-----------------
      //print the timings
      //-----------------
      //load the information
      sprintf(msg,"\n");pout << msg;
      sprintf(msg,"\nTimings within BLOCK-NEVPT2:");pout << msg;
      sprintf(msg,"\nPreparation Time          = %7.6lf sec",Info.getTime(0));pout << msg;
      sprintf(msg,"\nDensity Formation Time    = %7.6lf sec",Info.getTime(1));pout << msg;
      sprintf(msg,"\nCoreEnergy Formation Time = %7.6lf sec",Info.getTime(2));pout << msg;
      sprintf(msg,"\nOperator Transformation   = %7.6lf sec",Info.getTime(3));pout << msg;
      sprintf(msg,"\nNEVPT2 Time               = %7.6lf sec",Info.getTime(4));pout << msg;
      sprintf(msg,"\n----------------------------------------");pout << msg;
      double TotTime = Info.getTime(0)+Info.getTime(1)+Info.getTime(2)+Info.getTime(3)+Info.getTime(4);
      sprintf(msg,"\nTotal Time                = %lf sec",TotTime);pout << msg;
      sprintf(msg,"\n----------------------------------------");pout << msg;
      
      sprintf(msg,"\n");pout << msg;
      sprintf(msg,"\nTime required for perturbation classes:");pout << msg;
      sprintf(msg,"\nVijab                     = %7.6lf sec",Info.GetTimePerClass(0));pout << msg;
      sprintf(msg,"\nViab                      = %7.6lf sec",Info.GetTimePerClass(1));pout << msg;
      sprintf(msg,"\nVija                      = %7.6lf sec",Info.GetTimePerClass(2));pout << msg;
      sprintf(msg,"\nVab                       = %7.6lf sec",Info.GetTimePerClass(3));pout << msg;
      sprintf(msg,"\nVij                       = %7.6lf sec",Info.GetTimePerClass(4));pout << msg;
      sprintf(msg,"\nVia                       = %7.6lf sec",Info.GetTimePerClass(5));pout << msg;
      sprintf(msg,"\nVa                        = %7.6lf sec",Info.GetTimePerClass(6));pout << msg;
      sprintf(msg,"\nVi                        = %7.6lf sec",Info.GetTimePerClass(7));pout << msg;
    }//rank==0;
    
  }
  
  //============================================================================
  //the Kronecker delta
  //============================================================================
  double delta(const int&i, const int &j){
    if (i==j) return 1.0;
    else return 0.0;
  }
  
  
  //============================================================================
  // Print A set of rotation matrices
  //============================================================================
  void PrintRotMat(const vector<Matrix> &RotMat,const char *Name){
    for (int i=0;i<RotMat.size();i++){
      char msg[512];
      sprintf(msg,"%s.%i.tmp",Name,i);
      PrintMatrix(RotMat[i],msg);
    }
  }
  void PrintRotMatProperties(const vector<Matrix> &RotMat){
    char msg[512];
    sprintf(msg,"\nRotationMat: ");pout << msg;
    for (int imat=0;imat<RotMat.size();imat++){
      sprintf(msg,"  (%i x %i)",RotMat[imat].Nrows(),RotMat[imat].Ncols());
      pout << msg;
    }//imat
  }
  
  //============================================================================
  // Add a Matrix to a another matrix (taking into account transposition)
  //============================================================================
  void ScaleAdd_(double d, const SparseMatrix& a, SparseMatrix& b){
    int conj_case=0;
    if ((a.conjugacy()=='t')&&(b.conjugacy()=='n')) conj_case = 1;
    if ((a.conjugacy()=='n')&&(b.conjugacy()=='t')) conj_case = 2;
    if ((a.conjugacy()=='t')&&(b.conjugacy()=='t')) conj_case = 3;
    for (int lQ = 0; lQ < a.nrows(); ++lQ){
      for (int rQ = 0; rQ < a.ncols(); ++rQ){
        if (a.allowed(lQ, rQ)){
          if (!b.allowed(lQ, rQ)){
            pout <<"Not a valid addition"<<endl;
          }
          assert(b.allowed(lQ, rQ));
          switch(conj_case){
            case 0:
            case 3:
              MatrixScaleAdd(d, a.operator_element(lQ, rQ), b.operator_element(lQ, rQ));
              break;
            case 1:
              MatrixScaleAdd(d, a.operator_element(lQ, rQ).t(), b.operator_element(lQ, rQ));
              break;
            case 2:
              pout << "Not a valid addition"<<endl;
              break;
          }//conj_case
        }//a.allowed
      }//rQ
    }//lQ
  }
  
  //============================================================================
  // Transpose a wavefunction
  //============================================================================
  void Transpose(const Wavefunction &oldWave, Wavefunction &NewWave, SpinBlock &big){
    for (int i = 0; i < NewWave.nrows(); ++i){
      for (int j = 0; j < NewWave.ncols(); ++j){
        if (NewWave.allowed(i, j)){
          assert(oldWave.allowed(j, i));
          assert(NewWave(i, j).Nrows() == oldWave(j, i).Ncols());
          assert(NewWave(i, j).Ncols() == oldWave(j, i).Nrows());
          const StateInfo* s = &big.get_stateInfo();
	  Matrix tmp = oldWave(j, i);
	  tmp = tmp.t(); // this is really a transpose, not a hermitian conjugate...
	  NewWave(i, j) = tmp;
	  int parity = 1;
          if (IsFermion(s->leftStateInfo->quanta[i]) && IsFermion(s->rightStateInfo->quanta[j]))
	    parity = -1;
	  int A, B, J;
	  A = s->leftStateInfo->quanta[i].get_s().getirrep();
	  B = s->rightStateInfo->quanta[j].get_s().getirrep();
	  J = oldWave.get_deltaQuantum(0).get_s().getirrep();
	  parity *= pow(-1.0, static_cast<int>( (3*A - B + J)/2));
	  if (parity == -1)
	    NewWave(i,j) *= -1.0;
        }//allowed
      }//j
    }//i
  }
  
  
  
  //============================================================================
  // Evaluate the action of H on a function V(tu), meaning:
  //            sigma(t,u) = H * V(t,u)
  //                       = H * O(t,u) |psi>
  //============================================================================
  void GenerateActionOfH(Wavefunction &WF,vector<boost::shared_ptr<WavefunctionArray> > &Vtu, 
                         vector<boost::shared_ptr <WavefunctionArray> > &Sigmatu,
                         SpinQuantum &TripOpQ, SpinQuantum &SingOpQ, const SpinBlock &big,
                         int *OrbWin){
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
  //Determine whether this process should skip generating operators based on 
  //E(dotsite,dotsite)
  //============================================================================
bool SkipOperator(const vector<int> &sites, const int &dotsite){
  char msg[512];
  //determine the number of active orbitals
  int NActive = abs(sites[sites.size()-1]-sites[0])+1;
  int start,stop;
  PALDivideLoop(start,stop,0,NActive);
  if (dotsite>=start&&dotsite<stop) return false;
  else return true;
  
}

  //============================================================================
  // Determine the size of a file
  //============================================================================
int FileSize(char *FileName){
  ifstream ifs;
  ifs.open(FileName,std::ios::binary);
  if (!ifs.is_open()){
    pout << "Could not Open File: " << FileName;
    ifs.close();
    return -1;
  }
  ifs.seekg(0,ios::end);
  int  pos = ifs.tellg();
  ifs.close();
  return pos;
}

  //============================================================================
  // Check the size of an array
  //============================================================================
void CheckSize(ThreeIndOpArray &CCD){
  int t,u,v;
  char msg[512];
  //check the actual length of the array
  CCD.OpenFileRead();
  CCD.ResetBuffer();
  int size = 0;
  bool EOA=false;
  while (!EOA){
    boost::shared_ptr<Cre> Otuv = CCD.GetOpFromBuffer(t,u,v,EOA);
    if (!EOA){
      sprintf(msg,"\nReading O(%i,%i,%i)",t,u,v);
      //mpi_message(msg);
      size++;
    }
  }
  sprintf(msg,"\nSize of Array = %i",size);
  mpi_message(msg);
  CCD.CloseFileRead();
}

  //==========================================================================
  // Check if one of a given set of SpinQuanta equals another SpinQuantum
  // return the position of the first element in the vector that equals the
  // reference
  //==========================================================================
  int CheckEquality(const SpinQuantum &Reference, const vector<SpinQuantum> &Test){
    for (int i=0;i<Test.size();i++){
      if (Test[i]==Reference) return i;
    }//i
    //if none of the elements equal the reference, retur a dummy -1
    return -1;
  }
  
  //============================================================================
  //establish batching over second index b that runs from start to end-1
  //============================================================================
  void EstablishBatching(vector<pair<int,int> > &batches, int MaxCore, int M, int start,
                     int end, int NOperators){
        double ElementSize_ = (double) M * M * sizeof(double);
    //determine how many sets of N operators can be stored without exceeding MaxCore
    double tmp = (double) MaxCore;
    tmp *= 1024.0 * 1024.0;
    tmp /= (ElementSize_ * (double) NOperators);
    int BatchSize = (int) tmp;
    if (BatchSize<=0) BatchSize = 1;
    //DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    BatchSize = 1;
    //DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Evaluate the start and stop indices for each batch
    int NIters = end-start;
    int NBatches = NIters / BatchSize + 1;
    batches.resize(NBatches);
    for (int i=0;i<NBatches;i++){
      batches[i].first  = start + BatchSize * i;
      batches[i].second = start + BatchSize * (i+1);
      if (i==NBatches-1) batches[i].second = start + BatchSize * i + NIters % BatchSize;
    }//i
  }
  
  void NEVPT2SaveRotationMatrix (const std::vector<int>& sites, const std::vector<Matrix>& m1, int state, const char *name)
{
  Timer disktimer;
  int rank = mpi_rank();
  if (rank == 0)
    {

      char file [5000];
      sprintf (file, "%s%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/",name, sites [0], "-", *sites.rbegin (), ".", mpi_rank(),".state",state, ".tmp");
      p1out << "\t\t\t Saving Rotation Matrix :: " << file << endl;
      std::ofstream ofs(file, std::ios::binary);
      boost::archive::binary_oarchive save_mat(ofs);
      save_mat << m1;
      ofs.close();
    }
}

void NEVPT2LoadRotationMatrix (const std::vector<int>& sites, std::vector<Matrix>& m1, int state, const char *name)
{
  Timer disktimer;
  int rank = mpi_rank();
  if (rank == 0)
  {
    char file [5000];
    //sprintf (file, "%s%s%d%s%d%s%d%s", dmrginp.load_prefix().c_str(), "/Rotation-", sites [0], "-", *sites.rbegin (), ".", mpigetrank(), ".tmp");
    sprintf (file, "%s%s%s%d%s%d%s%d%s%d%s", dmrginp.save_prefix().c_str(), "/",name, sites [0], "-", *sites.rbegin (), ".", mpi_rank(),".state",state, ".tmp");
    p1out << "\t\t\t Loading Rotation Matrix :: " << file << endl;
    std::ifstream ifs(file, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> m1;
    ifs.close();
  }
}
}


